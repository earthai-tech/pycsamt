"""
pycsamt.agents.static_shift
============================

:class:`StaticShiftAgent` — Detect and correct galvanic static shift.

Wraps :mod:`pycsamt.emtools.ss`:

* :func:`~pycsamt.emtools.ss.correct_ss_ama`             — AMA correction (default)
* :func:`~pycsamt.emtools.ss.estimate_ss_loess`          — LOESS smooth estimate
* :func:`~pycsamt.emtools.ss.estimate_ss_refmedian`      — reference-median estimate
* :func:`~pycsamt.emtools.ss.estimate_ss_bilateral`      — bilateral-filter estimate
* :func:`~pycsamt.emtools.ss.apply_ss_factors`           — applies loess/refmedian/bilateral factors
* :func:`~pycsamt.emtools.ss.ss_comparison_psection`     — before/after section
* :func:`~pycsamt.emtools.ss.plot_ss_summary`            — summary dashboard
* :func:`~pycsamt.emtools.ss.plot_ss_1d_curves`          — per-station 1-D curves

Supported methods: ``"ama"`` (default), ``"loess"``, ``"refmedian"``,
``"bilateral"``.
"""

from __future__ import annotations

import logging
import time
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in galvanic distortion and static-shift correction for \
magnetotelluric data.
Given a static-shift correction summary, write 3–4 sentences that:
1. State whether significant static shift was detected.
2. Identify stations with the largest corrections and their magnitude.
3. Assess whether the correction method was appropriate for this dataset.
4. Recommend any follow-up action (e.g., additional spatial filtering).
Reply in plain English. No bullet points or markdown.
"""


class StaticShiftAgent(BaseAgent):
    """Detect and correct galvanic static shift in MT/AMT data.

    Parameters
    ----------
    api_key, model, llm_provider : str
    method : {"ama", "loess", "refmedian", "bilateral"}
        Correction algorithm.  Default ``"ama"`` (adaptive moving average).
    half_window : int
        Spatial half-window for AMA / LOESS smoothing.
    pband : (T_min, T_max) or None
        Period band used to estimate shift factors.
    inplace : bool
        Modify the input Sites in-place.  Default ``False`` (returns a copy).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``method`` : str, optional  — overrides constructor default
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``corrected_sites``   Sites with static shift removed
    ``shift_factors``     dict {station: factor}
    ``rho_before``        ndarray (n_freq × n_sta) — log₁₀ ρa before
    ``rho_after``         ndarray — log₁₀ ρa after
    ``delta_stats``       dict — min/max/mean shift magnitude
    ``figures``           dict — matplotlib Figure objects
    ``figure_paths``      dict — saved file paths

    Examples
    --------
    >>> agent  = StaticShiftAgent(method="ama")
    >>> result = agent.execute({"path": "/data/L22PLT",
    ...                         "output_dir": "/out/ss"})
    >>> result["delta_stats"]
    {'mean': 0.18, 'max': 0.42, 'n_shifted': 7}
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        method: str = "ama",
        half_window: int = 3,
        pband: tuple[float, float] | None = None,
        inplace: bool = False,
    ) -> None:
        super().__init__(
            "StaticShiftAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.method = method
        self.half_window = half_window
        self.pband = pband
        self.inplace = inplace

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── resolve sites ─────────────────────────────────────────────────────
        from ..emtools._core import ensure_sites

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed(
                "No 'sites' or 'path'.", elapsed=time.time() - t0
            )
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        method = str(
            input_data.get("method")
            or input_data.get("ss_method")
            or self.method
        ).lower()
        output_dir = input_data.get("output_dir")

        # ── import ss functions ───────────────────────────────────────────────
        from ..emtools.ss import (
            apply_ss_factors,
            correct_ss_ama,
            estimate_ss_bilateral,
            estimate_ss_loess,
            estimate_ss_refmedian,
            plot_ss_1d_curves,
            plot_ss_summary,
            ss_comparison_psection,
        )

        # ── capture log₁₀ ρa before correction ───────────────────────────────
        rho_before, freqs_all, sta_labels = _collect_rho(sites)

        # ── apply correction ──────────────────────────────────────────────────
        corrected_sites = None
        shift_factors: dict[str, float] = {}

        try:
            if method in ("none", "skip"):
                # User explicitly chose no correction.
                corrected_sites = sites
                warnings.append(
                    "Static-shift correction skipped (method='none')."
                )
            elif method == "ama":
                corrected_sites = correct_ss_ama(
                    sites,
                    half_window=self.half_window,
                    pband=self.pband,
                    inplace=self.inplace,
                    verbose=0,
                )
            elif method == "loess":
                result_loess = estimate_ss_loess(
                    sites,
                    half_window=self.half_window,
                    pband=self.pband,
                    verbose=0,
                )
                corrected_sites = apply_ss_factors(
                    sites,
                    result_loess,
                    inplace=self.inplace,
                    verbose=0,
                )
            elif method == "bilateral":
                result_bilateral = estimate_ss_bilateral(
                    sites,
                    half_window=self.half_window,
                    pband=self.pband,
                    verbose=0,
                )
                corrected_sites = apply_ss_factors(
                    sites,
                    result_bilateral,
                    inplace=self.inplace,
                    verbose=0,
                )
            elif method == "refmedian":
                result_ss = estimate_ss_refmedian(
                    sites,
                    pband=self.pband,
                    verbose=0,
                )
                if hasattr(result_ss, "iterrows"):
                    for _, row in result_ss.iterrows():
                        st = str(row.get("station", ""))
                        fac = float(row.get("fac_z", row.get("factor", 1.0)))
                        shift_factors[st] = fac
                    corrected_sites = apply_ss_factors(
                        sites,
                        result_ss,
                        inplace=self.inplace,
                        verbose=0,
                    )
                else:
                    warnings.append(
                        "estimate_ss_refmedian returned"
                        " unexpected type; using AMA."
                    )
                    corrected_sites = correct_ss_ama(
                        sites,
                        half_window=self.half_window,
                        inplace=self.inplace,
                        verbose=0,
                    )
            else:
                warnings.append(
                    f"Unknown method {method!r}; falling back to AMA."
                )
                corrected_sites = correct_ss_ama(
                    sites,
                    half_window=self.half_window,
                    inplace=self.inplace,
                    verbose=0,
                )
        except Exception as exc:
            logger.warning(
                f"Static-shift correction ({method}) failed: {exc}",
                exc_info=True,
            )
            warnings.append(
                f"Static-shift correction"
                f" ({method}) failed: {exc}."
                " Raw (uncorrected) data used."
            )
            corrected_sites = sites

        # ── capture log₁₀ ρa after correction ────────────────────────────────
        rho_after = (
            _collect_rho(corrected_sites)[0]
            if corrected_sites is not None
            else rho_before
        )

        # ── delta statistics ──────────────────────────────────────────────────
        delta_stats: dict[str, float] = {}
        if rho_before is not None and rho_after is not None:
            try:
                delta = rho_after - rho_before
                col_delta = np.nanmedian(np.abs(delta), axis=0)
                delta_stats = {
                    "mean": float(np.nanmean(col_delta)),
                    "max": float(np.nanmax(col_delta)),
                    "min": float(np.nanmin(col_delta)),
                    "n_shifted": int(np.sum(col_delta > 0.05)),
                }
            except Exception:
                pass

        # ── figures ───────────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        if (
            rho_before is not None
            and rho_after is not None
            and freqs_all is not None
        ):
            # summary dashboard: before / after / delta
            try:
                # _collect_rho is (n_freq, n_sta); plotters want (n_st, n_f)
                fig_sum = plot_ss_summary(
                    rho_before.T,
                    rho_after.T,
                    freqs=freqs_all,
                    station_labels=sta_labels,
                )
                figures["ss_summary"] = fig_sum
                p = self._save_figure(
                    fig_sum, output_dir, "ss_summary", warnings_list=warnings
                )
                if p:
                    fig_paths["ss_summary"] = p
            except Exception as exc:
                warnings.append(f"plot_ss_summary: {exc}")

            # per-station 1-D ρa curves
            try:
                fig_1d = plot_ss_1d_curves(
                    rho_before.T,
                    rho_after.T,
                    freqs=freqs_all,
                    station_labels=sta_labels,
                )
                figures["ss_curves"] = fig_1d
                p = self._save_figure(
                    fig_1d, output_dir, "ss_1d_curves", warnings_list=warnings
                )
                if p:
                    fig_paths["ss_curves"] = p
            except Exception as exc:
                warnings.append(f"plot_ss_1d_curves: {exc}")

        # comparison pseudosection (before / after side-by-side)
        try:
            fig_cmp = ss_comparison_psection(
                sites,
                method=method,
                verbose=0,
            )
            if fig_cmp is not None:
                fig_c = (
                    fig_cmp
                    if hasattr(fig_cmp, "savefig")
                    else (
                        fig_cmp.get_figure()
                        if hasattr(fig_cmp, "get_figure")
                        else None
                    )
                )
                if fig_c is not None:
                    figures["ss_comparison"] = fig_c
                    p = self._save_figure(
                        fig_c,
                        output_dir,
                        "ss_comparison",
                        warnings_list=warnings,
                    )
                    if p:
                        fig_paths["ss_comparison"] = p
        except Exception as exc:
            warnings.append(f"ss_comparison_psection: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and delta_stats:
            prompt = (
                f"Static-shift correction summary ({method}):\n"
                f"  Mean correction magnitude: {delta_stats.get('mean', '?'):.3f} "
                f"  log₁₀(Ω·m)\n"
                f"  Max correction: {delta_stats.get('max', '?'):.3f}\n"
                f"  Stations significantly shifted (>0.05): "
                f"{delta_stats.get('n_shifted', '?')}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Assess the static-shift correction result."
            )
            interp = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        n_sig = delta_stats.get("n_shifted", 0)
        # When nothing was corrected, say so explicitly and explain why —
        # otherwise the identical before/after curves and empty delta
        # look like a bug rather than the expected "no shift detected".
        if n_sig == 0 and method not in ("none", "skip"):
            note = (
                "No significant static shift was detected, so the data is"
                " unchanged and the before/after curves are identical. If"
                " you expected corrections, the dataset may be strongly"
                " 3-D (high skew), too sparse in the chosen period band,"
                " or have too few stations for spatial comparison."
            )
            warnings.append(note)
            summary = (
                f"Static-shift ({method}): no significant shift detected"
                f" — data unchanged. {len(figures)} figure(s) produced."
            )
        else:
            summary = (
                f"Static-shift correction ({method}) applied. "
                f"{n_sig} station(s) shifted >0.05 log₁₀. "
                f"{len(figures)} figure(s) produced."
            )
        return AgentResult(
            status="success",
            summary=summary,
            data={
                "corrected_sites": corrected_sites,
                "shift_factors": shift_factors,
                "rho_before": rho_before,
                "rho_after": rho_after,
                "delta_stats": delta_stats,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── helpers ───────────────────────────────────────────────────────────────────


def _collect_rho(
    sites: Any,
) -> tuple[np.ndarray | None, np.ndarray | None, list[str]]:
    """Collect log₁₀ ρa_xy into a 2-D array (n_freq × n_sta).

    Returns (rho_mat, freqs, station_labels). Any of these may be None if
    data are unavailable.
    """
    from ..emtools._core import (
        _get_z_block,
        _iter_items,
        _name,
    )

    cols: list[np.ndarray] = []
    freqs_ref: np.ndarray | None = None
    labels: list[str] = []

    for i, ed in enumerate(_iter_items(sites)):
        nm = _name(ed, i)
        Z_obj, z, fr = _get_z_block(ed)
        if z is None or fr is None:
            continue
        # Always compute from Z — ed.rho is a
        # cached attribute that is stale after
        # impedance-tensor correction modifies Z.
        rho_xy = (0.2 / np.where(fr == 0, np.nan, fr)) * np.abs(
            z[:, 0, 1]
        ) ** 2

        log_rho = np.log10(np.clip(rho_xy, 1e-6, None))

        if freqs_ref is None:
            freqs_ref = fr
            cols.append(log_rho)
            labels.append(nm)
        else:
            # interpolate onto common freq grid
            if len(fr) == len(freqs_ref):
                cols.append(log_rho)
                labels.append(nm)
            else:
                # skip stations with different freq counts for now
                cols.append(
                    log_rho[: len(freqs_ref)]
                    if len(log_rho) >= len(freqs_ref)
                    else np.full(len(freqs_ref), np.nan)
                )
                labels.append(nm)

    if not cols or freqs_ref is None:
        return None, None, []

    n_freq = len(freqs_ref)
    mat = np.full((n_freq, len(cols)), np.nan)
    for j, col in enumerate(cols):
        n = min(len(col), n_freq)
        mat[:n, j] = col[:n]

    return mat, freqs_ref, labels


__all__ = ["StaticShiftAgent"]
