# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.agents.tsdenoise_agent
================================

:class:`TimeSeriesDenoisingAgent` — Raw field time-series denoising
(MMF-SVM-K-SVD, Gui et al., 2024).

Wraps :class:`~pycsamt.ai.processing.tsdenoise.TimeSeriesDenoiser`: one
stage *upstream* of everything else in :mod:`pycsamt.ai.processing` --
it works on a :class:`~pycsamt.ts.TSData` record (Ex, Ey, Hx, Hy, Hz
samples) before :func:`~pycsamt.ts.process.ts_to_spectra` turns it into
the impedance spectra :class:`~pycsamt.agents.denoising.DenoisingAgent`'s
``EMDenoiser`` operates on. Distinct from that agent: this one denoises
raw voltage/field samples, not frequency-domain impedance.

Combines mathematical-morphological filtering (MMF) to split a channel
into a smooth low-frequency component and a high-frequency residual, an
SVM classifier that flags windows of that residual as clean vs. noisy,
and K-SVD dictionary-learning denoising applied only to flagged windows.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in raw MT/AMT field time-series signal processing
(cultural noise, powerline harmonics, instrument transients).
Given a time-series denoising result summary, write 3-4 sentences that:
1. State how many windows were flagged noisy and denoised, out of the total.
2. Note the channel(s) processed and roughly how much of the record was
   affected.
3. Flag anything unusual (e.g. an unusually high flagged fraction suggesting
   a systematic acquisition problem, not isolated bursts).
4. Recommend whether the denoised record is ready for spectral processing
   (ts_to_spectra) or needs a closer look first.
Reply in plain English.
"""


class TimeSeriesDenoisingAgent(BaseAgent):
    """Denoise a raw MT field time-series record before spectral processing.

    Parameters
    ----------
    api_key, model, llm_provider : str
    mmf_size : int
        Morphological filter window (samples) separating the smooth
        low-frequency component from the high-frequency residual the
        classifier/K-SVD operate on (default 121).
    win_seconds : float
        Classification window length in seconds (default 60.0).
    channels : list of str or None
        Channels to denoise. ``None`` (default) processes every
        channel present in the loaded record.
    summary_channel : str or None
        Which channel the summary figure focuses on. ``None`` uses the
        first processed channel.

    Input keys
    ----------
    ``ts`` : TSData, optional
    ``ts_path`` : str, optional — read via :func:`pycsamt.ts.read_ts`
        (exactly one of ``ts`` / ``ts_path`` is required)
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``denoised_ts``    TSData with each processed channel denoised
    ``diagnostics``    pandas DataFrame — per-window label/score, from the
                        summary channel's classification pass
    ``n_windows``      int
    ``n_flagged``      int
    ``figures``        dict
    ``figure_paths``   dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        mmf_size: int = 121,
        win_seconds: float = 60.0,
        channels: list[str] | None = None,
        summary_channel: str | None = None,
    ) -> None:
        super().__init__(
            "TimeSeriesDenoisingAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.mmf_size = mmf_size
        self.win_seconds = win_seconds
        self.channels = channels
        self.summary_channel = summary_channel

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        try:
            from ..ai.processing import TimeSeriesDenoiser, mmf_split
        except ImportError as exc:
            return AgentResult.failed(
                f"TimeSeriesDenoisingAgent import failed: {exc}",
                elapsed=time.time() - t0,
            )

        ts = input_data.get("ts")
        ts_path = input_data.get("ts_path")
        if ts is None and ts_path is None:
            return AgentResult.failed(
                "No 'ts' or 'ts_path'.", elapsed=time.time() - t0
            )
        if ts is None:
            try:
                from ..ts import read_ts

                ts = read_ts(str(ts_path))
            except Exception as exc:
                return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        if ts.dt is None:
            return AgentResult.failed(
                "Loaded time series has no sampling interval (ts.dt is "
                "None); TimeSeriesDenoiser.apply() requires it.",
                elapsed=time.time() - t0,
            )

        output_dir = input_data.get("output_dir")
        chans = list(self.channels) if self.channels else list(ts.channels())
        if not chans:
            return AgentResult.failed(
                "Loaded time series has no channels.", elapsed=time.time() - t0
            )

        try:
            denoiser = TimeSeriesDenoiser(
                mmf_size=self.mmf_size,
                win_seconds=self.win_seconds,
                random_state=0,
            )
            denoised_ts = denoiser.apply(ts, channels=chans, inplace=False)
            diag = denoiser.diagnostics_
        except Exception as exc:
            return AgentResult.failed(
                f"TimeSeriesDenoiser.apply failed: {exc}",
                elapsed=time.time() - t0,
            )

        n_windows = int(len(diag))
        n_flagged = int((diag["label"] == "noisy").sum()) if n_windows else 0

        # ── figures ───────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}
        summary_chan = self.summary_channel or chans[0]
        try:
            from ..ai.processing import plot_ts_denoise_summary

            raw = ts.get(summary_chan)
            den = denoised_ts.get(summary_chan)
            t = np.arange(raw.size) * ts.dt
            low, high = mmf_split(raw, size=self.mmf_size)
            fig = plot_ts_denoise_summary(
                t,
                raw,
                den,
                diag,
                low=low,
                high=high,
                suptitle=f"TimeSeriesDenoiser -- {summary_chan}",
            )
            figures["ts_denoise_summary"] = fig
            p = self._save_figure(
                fig, output_dir, "ts_denoise_summary", warnings_list=warnings
            )
            if p:
                fig_paths["ts_denoise_summary"] = p
        except Exception as exc:
            warnings.append(f"Time-series denoise summary figure: {exc}")

        # ── LLM interpretation ──────────────────────────────────────────────
        interp: str | None = None
        if self.llm_available:
            prompt = (
                f"Time-series denoising summary:\n"
                f"  Channels processed: {chans}\n"
                f"  Summary channel: {summary_chan}\n"
                f"  Windows: {n_windows}, flagged noisy: {n_flagged} "
                f"({100 * n_flagged / max(n_windows, 1):.0f}%)\n"
                f"  mmf_size={self.mmf_size} samples, "
                f"win_seconds={self.win_seconds}\n\n"
                "Assess the record's noise conditions and recommend follow-up."
            )
            interp = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        return AgentResult(
            status="success",
            summary=(
                f"Denoised {len(chans)} channel(s); {n_flagged}/{n_windows} "
                f"window(s) flagged noisy on '{summary_chan}'. "
                f"{len(figures)} figure(s)."
            ),
            data={
                "denoised_ts": denoised_ts,
                "diagnostics": diag,
                "n_windows": n_windows,
                "n_flagged": n_flagged,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["TimeSeriesDenoisingAgent"]
