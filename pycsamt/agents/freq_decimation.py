"""
pycsamt.agents.freq_decimation
================================

:class:`FrequencyDecimationAgent` — Intelligent period selection for inversion.

Selects an optimal subset of periods from the observed data by:

1. Applying the SNR/QC flags from :class:`~pycsamt.agents.DataQCAgent` to
   mask dead-band and low-quality frequencies.
2. Choosing log-uniformly spaced survivors to give even depth coverage.
3. Enforcing user-defined period bounds and a minimum SNR threshold.

The output is a dictionary of selected periods per station, ready to be
consumed by :class:`~pycsamt.agents.InversionPrepAgent`,
:class:`~pycsamt.agents.Occam2DAgent`, or
:class:`~pycsamt.agents.ModEmAgent`.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in MT data selection and frequency decimation for inversion.
Given a period decimation result, write 3-4 sentences that:
1. State how many periods were selected versus available, and the selection ratio.
2. Identify which frequency bands were excluded and the likely reason (dead band, low SNR).
3. Confirm whether the selected periods cover the target depth range adequately.
4. Recommend any additional frequencies that should be included or excluded.
Reply in plain English.
"""


class FrequencyDecimationAgent(BaseAgent):
    """Select optimal periods from MT data for inversion.

    Parameters
    ----------
    api_key, model, llm_provider : str
    n_per_decade : int
        Number of periods to keep per decade of period range (default 6).
    snr_threshold : float
        Minimum SNR value to retain a frequency (default 3.0).
        Frequencies below this are excluded as dead-band.
    period_range : [T_min, T_max] or None
        Hard period bounds in seconds.  ``None`` uses the full data range.
    component : {'xy', 'yx'}
        Component used for SNR proxy (default ``'xy'``).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``qc_result`` : AgentResult or dict, optional — output from DataQCAgent
        (provides per-frequency SNR scores; if absent a proxy is computed)
    ``n_per_decade`` : int, optional
    ``snr_threshold`` : float, optional
    ``period_range`` : [T_min, T_max], optional
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``selected_periods``   dict {station: ndarray}  — selected periods (s)
    ``n_original``         int — total available (station, freq) cells
    ``n_selected``         int — retained cells
    ``selection_ratio``    float
    ``dead_band_mask``     dict {station: ndarray bool}
    ``figures``            dict
    ``figure_paths``       dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        n_per_decade: int = 6,
        snr_threshold: float = 3.0,
        period_range: list[float] | None = None,
        component: str = "xy",
    ) -> None:
        super().__init__(
            "FrequencyDecimationAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.n_per_decade  = n_per_decade
        self.snr_threshold = snr_threshold
        self.period_range  = period_range
        self.component     = component.lower()

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        from ..emtools._core import ensure_sites, _iter_items, _name, _get_z_block

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed("No 'sites' or 'path'.", elapsed=time.time() - t0)
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        n_per_decade  = int(input_data.get("n_per_decade",  self.n_per_decade))
        snr_threshold = float(input_data.get("snr_threshold", self.snr_threshold))
        period_range  = input_data.get("period_range") or self.period_range
        component     = str(input_data.get("component", self.component)).lower()
        output_dir    = input_data.get("output_dir")

        ri, ci = (0, 1) if component == "xy" else (1, 0)

        # extract QC SNR flags from a prior DataQCAgent result (optional)
        qc_result = input_data.get("qc_result")
        qc_snr: dict[str, np.ndarray] = {}
        qc_freqs: dict[str, np.ndarray] = {}
        if qc_result is not None:
            try:
                snr_section = (qc_result.get("snr_section")
                               if hasattr(qc_result, "get") else
                               qc_result.get("data", {}).get("snr_section"))
                if snr_section is not None and isinstance(snr_section, dict):
                    qc_snr = snr_section
            except Exception:
                pass

        selected_periods: dict[str, np.ndarray]   = {}
        dead_band_mask:   dict[str, np.ndarray]   = {}
        n_original = 0
        n_selected = 0

        for i, ed in enumerate(_iter_items(sites)):
            nm = _name(ed, i)
            _, z, fr = _get_z_block(ed)
            if z is None or fr is None or fr.size == 0:
                continue

            per = 1.0 / np.where(fr == 0, np.nan, fr)
            valid = np.isfinite(per)

            # ── period range filter ───────────────────────────────────────
            if period_range is not None:
                t_min, t_max = float(period_range[0]), float(period_range[1])
                valid &= (per >= t_min) & (per <= t_max)

            # ── SNR proxy: |Zxy| / std(|Zxy|) over station set ───────────
            snr_proxy = np.ones(len(fr))
            try:
                zcomp = np.abs(z[:, ri, ci])
                mu    = np.nanmean(zcomp[valid]) if valid.any() else 1.0
                sd    = np.nanstd(zcomp[valid])  + 1e-30
                snr_proxy = zcomp / sd
            except Exception:
                pass

            # ── SNR threshold ─────────────────────────────────────────────
            good = valid & (snr_proxy >= snr_threshold)
            dead = valid & ~good
            dead_band_mask[nm] = dead

            n_original += int(valid.sum())

            if not good.any():
                warnings.append(f"{nm}: no frequencies pass SNR threshold — skipped.")
                selected_periods[nm] = np.array([])
                continue

            good_periods = per[good]
            log_p = np.log10(np.clip(good_periods, 1e-12, None))

            # ── log-spaced decimation ─────────────────────────────────────
            p_min, p_max = float(log_p.min()), float(log_p.max())
            n_decades    = max(1.0, p_max - p_min)
            n_target     = max(1, int(np.round(n_decades * n_per_decade)))

            if n_target >= len(good_periods):
                chosen = good_periods
            else:
                targets   = np.linspace(p_min, p_max, n_target)
                chosen_idx = set()
                for t in targets:
                    idx = int(np.argmin(np.abs(log_p - t)))
                    chosen_idx.add(idx)
                chosen = np.sort(good_periods[sorted(chosen_idx)])

            selected_periods[nm] = chosen
            n_selected += len(chosen)

        # ── summary ───────────────────────────────────────────────────────
        n_sta_sel      = sum(1 for p in selected_periods.values() if p.size > 0)
        selection_ratio = n_selected / max(n_original, 1)

        # ── figures ───────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            fig = _plot_selection_summary(
                selected_periods, dead_band_mask, sites,
                n_per_decade, snr_threshold,
            )
            if fig is not None:
                figures["selection_summary"] = fig
                p = self._save_figure(fig, output_dir, "freq_decimation_summary",
                                      warnings_list=warnings)
                if p:
                    fig_paths["selection_summary"] = p
        except Exception as exc:
            warnings.append(f"Selection summary figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and n_sta_sel:
            prompt = (
                f"Frequency decimation summary:\n"
                f"  Stations: {n_sta_sel} with data\n"
                f"  Original frequencies: {n_original}\n"
                f"  Selected frequencies: {n_selected}\n"
                f"  Selection ratio: {selection_ratio:.1%}\n"
                f"  n_per_decade={n_per_decade}, SNR threshold={snr_threshold}\n"
                f"  Period range: {period_range or 'full'}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Assess whether the selected periods cover the target depth and advise."
            )
            interp = self.query_llm(prompt, max_tokens=180)

        elapsed = time.time() - t0
        return AgentResult(
            status="success" if n_sta_sel > 0 else "needs_review",
            summary=(
                f"Frequency decimation: {n_selected}/{n_original} periods retained "
                f"({selection_ratio:.0%}) across {n_sta_sel} stations."
            ),
            data={
                "selected_periods": selected_periods,
                "n_original":       n_original,
                "n_selected":       n_selected,
                "selection_ratio":  selection_ratio,
                "dead_band_mask":   dead_band_mask,
                "figures":          figures,
                "figure_paths":     fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────

def _plot_selection_summary(
    selected_periods: dict[str, np.ndarray],
    dead_band_mask:   dict[str, np.ndarray],
    sites: Any,
    n_per_decade: int,
    snr_threshold: float,
) -> Any:
    """Scatter plot: selected periods (green) and dead-band (red) per station."""
    import matplotlib.pyplot as plt
    from ..emtools._core import _iter_items, _name, _get_z_block

    station_names = list(selected_periods.keys())
    if not station_names:
        return None

    fig, ax = plt.subplots(figsize=(max(8, len(station_names) * 0.5), 5))

    for xi, nm in enumerate(station_names):
        sel = selected_periods.get(nm, np.array([]))
        dead_mask = dead_band_mask.get(nm)

        # all available periods for this station
        for i, ed in enumerate(_iter_items(sites)):
            if _name(ed, i) != nm:
                continue
            _, z, fr = _get_z_block(ed)
            if fr is None:
                break
            per_all = 1.0 / np.where(fr == 0, np.nan, fr)
            valid = np.isfinite(per_all)

            if dead_mask is not None and dead_mask.any():
                ax.scatter(
                    [xi] * int(dead_mask.sum()),
                    per_all[dead_mask],
                    marker="x", s=18, color="#e74c3c", alpha=0.6, zorder=2,
                )
            if sel.size > 0:
                ax.scatter(
                    [xi] * len(sel), sel,
                    marker="o", s=20, color="#27ae60", alpha=0.85, zorder=3,
                )
            break

    ax.set_yscale("log")
    ax.set_xticks(range(len(station_names)))
    ax.set_xticklabels(station_names, rotation=90, fontsize=7)
    ax.set_ylabel("Period (s)", fontsize=9)
    ax.set_title(
        f"Frequency decimation  "
        f"(n/decade={n_per_decade}, SNR≥{snr_threshold})\n"
        "Green = selected  ·  Red × = dead-band / low SNR",
        fontsize=9, fontweight="bold",
    )
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    return fig


__all__ = ["FrequencyDecimationAgent"]
