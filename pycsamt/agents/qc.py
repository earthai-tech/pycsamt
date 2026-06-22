"""
pycsamt.agents.qc
=================

:class:`DataQCAgent` — MT data quality control and frequency editing.

Wraps :mod:`pycsamt.emtools.qc`:

* :func:`~pycsamt.emtools.qc.build_qc_table`        — per-station metrics
* :func:`~pycsamt.emtools.qc.qc_flags`               — pass / fail flags
* :func:`~pycsamt.emtools.qc.frequency_confidence_table` — per-frequency scores
* :func:`~pycsamt.emtools.qc.plot_frequency_confidence_psection` — section figure
* :func:`~pycsamt.emtools.qc.plot_confidence_profile`            — profile figure
* :func:`~pycsamt.emtools.qc.station_confidence_table`           — per-station confidence

Output figures use :data:`~pycsamt.api.section.PYCSAMT_SECTION` so they are
consistent with all other pycsamt plots.
"""

from __future__ import annotations

import time
from typing import Any

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert MT/AMT/CSAMT data quality analyst for pycsamt v2.
Given a survey QC summary, write 3–4 sentences that:
1. State the overall data quality (good / moderate / poor).
2. Identify specific stations or frequency bands that need attention.
3. Explain the likely cause (instrument noise, EM interference, near-field).
4. Recommend the most important next processing step.
Reply in plain English. No bullet points, no markdown headings.
"""


class DataQCAgent(BaseAgent):
    """Run data quality control on a MT/AMT dataset.

    Parameters
    ----------
    api_key, model, llm_provider : str
        LLM configuration (optional).
    method : str
        Confidence scoring method: ``"composite"`` (default), ``"presence"``,
        ``"snr"``, or ``"spatial"``.
    min_frac_ok : float
        Minimum fraction of OK frequencies for a station to pass (0–1).
    min_snr_med : float
        Minimum median SNR for a station to pass.
    max_skew_med : float
        Maximum median |β| skewness for a station to pass.

    Input keys
    ----------
    ``sites`` : Sites or ``path`` : str
        EDI data to assess.
    ``output_dir`` : str, optional
        Where to save QC figures.
    ``period_range`` : [T_min, T_max], optional
        Restrict QC to this period window.

    Output data keys
    ----------------
    ``qc_table``            pandas DataFrame — per-station metrics
    ``flags``               pandas DataFrame — pass / fail per station
    ``confidence_table``    pandas DataFrame — per-station confidence scores
    ``freq_conf_table``     pandas DataFrame — per-frequency confidence
    ``n_flagged``           int
    ``flagged_stations``    list[str]
    ``figures``             dict — matplotlib Figure objects
    ``figure_paths``        dict — saved file paths (when output_dir set)

    Examples
    --------
    >>> agent  = DataQCAgent()
    >>> result = agent.execute({"path": "/data/L22PLT",
    ...                         "output_dir": "/out/qc"})
    >>> result["n_flagged"]
    2
    >>> result["figures"]["confidence_section"]
    <Figure …>
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        method: str = "composite",
        min_frac_ok: float = 0.6,
        min_snr_med: float = 2.0,
        max_skew_med: float = 6.0,
    ) -> None:
        super().__init__(
            "DataQCAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.method       = method
        self.min_frac_ok  = min_frac_ok
        self.min_snr_med  = min_snr_med
        self.max_skew_med = max_skew_med

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── resolve sites ─────────────────────────────────────────────────────
        from ..emtools._core import ensure_sites
        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed(
                "No 'sites' or 'path' key in input_data.",
                elapsed=time.time() - t0,
            )
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        output_dir = input_data.get("output_dir")

        # ── import emtools.qc functions ───────────────────────────────────────
        from ..emtools.qc import (
            build_qc_table,
            qc_flags,
            frequency_confidence_table,
            station_confidence_table,
            plot_frequency_confidence_psection,
            plot_confidence_profile,
        )

        # ── build tables ──────────────────────────────────────────────────────
        qc_table = conf_table = freq_conf = flags_df = None
        try:
            qc_table = build_qc_table(sites, verbose=0)
        except Exception as exc:
            warnings.append(f"build_qc_table: {exc}")

        try:
            flags_df = qc_flags(
                sites,
                min_frac_ok=self.min_frac_ok,
                min_snr_med=self.min_snr_med,
                max_skew_med=self.max_skew_med,
                verbose=0,
            )
        except Exception as exc:
            warnings.append(f"qc_flags: {exc}")

        try:
            conf_table = station_confidence_table(sites, method=self.method, verbose=0)
        except Exception as exc:
            warnings.append(f"station_confidence_table: {exc}")

        try:
            freq_conf = frequency_confidence_table(sites, method=self.method, verbose=0)
        except Exception as exc:
            warnings.append(f"frequency_confidence_table: {exc}")

        # ── flagged stations ──────────────────────────────────────────────────
        flagged: list[str] = []
        if flags_df is not None and hasattr(flags_df, "iterrows"):
            for _, row in flags_df.iterrows():
                if not bool(row.get("pass", True)):
                    flagged.append(str(row.get("station", "")))
        n_flagged = len(flagged)

        # ── figures ───────────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        # confidence pseudosection
        try:
            ax_conf = plot_frequency_confidence_psection(
                sites, method=self.method,
                section="pseudosection", verbose=0,
            )
            fig = ax_conf.get_figure() if hasattr(ax_conf, "get_figure") else ax_conf
            figures["confidence_section"] = fig
            p = self._save_figure(fig, output_dir, "qc_confidence_section",
                                  warnings_list=warnings)
            if p:
                fig_paths["confidence_section"] = p
        except Exception as exc:
            warnings.append(f"plot_frequency_confidence_psection: {exc}")

        # confidence profile (per-station score bar)
        try:
            ax_prof = plot_confidence_profile(sites, method=self.method, verbose=0)
            fig = ax_prof.get_figure() if hasattr(ax_prof, "get_figure") else ax_prof
            figures["confidence_profile"] = fig
            p = self._save_figure(fig, output_dir, "qc_confidence_profile",
                                  warnings_list=warnings)
            if p:
                fig_paths["confidence_profile"] = p
        except Exception as exc:
            warnings.append(f"plot_confidence_profile: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key:
            n_st = len(qc_table) if qc_table is not None else "?"
            prompt = (
                f"Survey QC summary:\n"
                f"  Stations evaluated: {n_st}\n"
                f"  Flagged (fail): {n_flagged}  — {flagged}\n"
                f"  Warnings: {warnings[:4] if warnings else 'none'}\n"
                f"  Method: {self.method}\n\n"
                "Assess data quality and recommend the next processing step."
            )
            interp = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        status = "success" if n_flagged == 0 else (
            "needs_review" if n_flagged < 5 else "needs_review"
        )
        return AgentResult(
            status=status,
            summary=(
                f"QC complete: {n_flagged} station(s) flagged out of "
                f"{len(qc_table) if qc_table is not None else '?'}. "
                f"{len(figures)} figure(s) produced."
            ),
            data={
                "qc_table":         qc_table,
                "flags":            flags_df,
                "confidence_table": conf_table,
                "freq_conf_table":  freq_conf,
                "n_flagged":        n_flagged,
                "flagged_stations": flagged,
                "figures":          figures,
                "figure_paths":     fig_paths,
                "sites":            sites,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["DataQCAgent"]
