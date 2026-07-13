"""
pycsamt.agents.phase_analysis
==============================

:class:`PhaseAnalysisAgent` — Phase tensor, strike, and dimensionality analysis.

Wraps:

* :func:`~pycsamt.emtools.tensor.build_phase_tensor_table`
* :func:`~pycsamt.emtools.tensor.plot_phase_tensor_psection`
* :func:`~pycsamt.emtools.tensor.plot_phase_tensor_rose`
* :func:`~pycsamt.emtools.strike.estimate_strike_consensus`
* :func:`~pycsamt.emtools.strike.plot_strike_analysis`
* :func:`~pycsamt.emtools.strike.plot_strike_rose`
* :func:`~pycsamt.emtools.dimensionality.classify_dimensionality`
* :func:`~pycsamt.emtools.dimensionality.plot_dim_confidence_grid`
* :func:`~pycsamt.emtools.advanced.plot_impedance_mohr_circles`
* :func:`~pycsamt.emtools.advanced.plot_survey_fingerprint`

All figures are governed by :data:`~pycsamt.api.section.PYCSAMT_SECTION`
and :data:`~pycsamt.api.style.PYCSAMT_STYLE`.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in MT phase tensor analysis and geological interpretation.
Given a survey phase tensor summary, write 4–5 sentences that:
1. State the dominant dimensionality (1-D, 2-D, or 3-D) with evidence.
2. Report the consensus geoelectric strike direction and its reliability.
3. Identify periods / depth ranges where 3-D structure becomes significant.
4. Note any anomalous stations (high skew, low ellipticity).
5. Recommend whether to rotate data to strike before inversion.
Reply in plain English. No bullet points or markdown.
"""


class PhaseAnalysisAgent(BaseAgent):
    """Run a full phase tensor, strike, and dimensionality survey analysis.

    Parameters
    ----------
    api_key, model, llm_provider : str
    skew_th : float
        Skewness |β| threshold for 3-D classification (°).
    ellipt_th : float
        Ellipticity λ threshold for 2-D classification.
    band : (T_min, T_max) or None
        Period band for strike estimation.

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``period_range`` : [T_min, T_max], optional
    ``output_dir`` : str, optional
    ``run_mohr`` : bool, optional — also produce Mohr circles (default False)
    ``run_fingerprint`` : bool, optional — produce fingerprint grid (default True)

    Output data keys
    ----------------
    ``pt_table``          pandas DataFrame — full PT metrics per (station, period)
    ``strike_consensus``  float — consensus strike angle (°)
    ``strike_iqr``        float — IQR of strike across all stations
    ``dim_table``         pandas DataFrame — per-(station, period) classification
    ``n_1d``, ``n_2d``, ``n_3d``   int — count of observations per class
    ``figures``           dict — matplotlib Figure objects
    ``figure_paths``      dict — saved file paths

    Examples
    --------
    >>> agent  = PhaseAnalysisAgent()
    >>> result = agent.execute({"path": "/data/L22PLT",
    ...                         "output_dir": "/out/pt"})
    >>> result["strike_consensus"]
    42.5
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        skew_th: float = 5.0,
        ellipt_th: float = 0.1,
        band: tuple[float, float] | None = None,
    ) -> None:
        super().__init__(
            "PhaseAnalysisAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.skew_th = skew_th
        self.ellipt_th = ellipt_th
        self.band = band

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

        output_dir = input_data.get("output_dir")
        period_range = input_data.get("period_range")
        run_mohr = bool(input_data.get("run_mohr", False))
        run_fingerprint = bool(input_data.get("run_fingerprint", True))
        band = period_range or self.band

        # ── phase tensor table ────────────────────────────────────────────────
        from ..emtools.dimensionality import (
            classify_dimensionality,
            plot_dim_confidence_grid,
        )
        from ..emtools.strike import (
            estimate_strike_consensus,
            plot_strike_analysis,
        )
        from ..emtools.tensor import (
            build_phase_tensor_table,
            plot_phase_tensor_psection,
            plot_phase_tensor_rose,
        )

        pt_table = None
        dim_table = None
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            pt_table = build_phase_tensor_table(sites, verbose=0)
        except Exception as exc:
            warnings.append(f"build_phase_tensor_table: {exc}")

        # ── dimensionality classification ─────────────────────────────────────
        try:
            dim_table = classify_dimensionality(
                sites,
                skew_th=self.skew_th,
                ellipt_th=self.ellipt_th,
                verbose=0,
            )
        except Exception as exc:
            warnings.append(f"classify_dimensionality: {exc}")

        # compute class counts from pt_table
        n_1d = n_2d = n_3d = 0
        if pt_table is not None and not pt_table.empty:
            try:
                beta = np.abs(pt_table["beta"].to_numpy(float))
                ellipt = pt_table["ellipt"].to_numpy(float)
                u3d = np.clip(beta / self.skew_th, 0, 1)
                u1d = (1 - u3d) * np.clip(1 - ellipt / self.ellipt_th, 0, 1)
                u2d = 1 - u1d - u3d
                dom = np.argmax(np.column_stack([u1d, u2d, u3d]), axis=1)
                n_1d = int((dom == 0).sum())
                n_2d = int((dom == 1).sum())
                n_3d = int((dom == 2).sum())
            except Exception:
                pass

        # ── strike estimation ─────────────────────────────────────────────────
        strike_consensus = np.nan
        strike_iqr = np.nan
        try:
            st_result = estimate_strike_consensus(
                sites,
                band=band,
                verbose=0,
            )
            # returns a DataFrame with columns: station, ang, iqr, lo, hi, n
            if hasattr(st_result, "ang"):
                ang_vals = st_result["ang"].dropna()
                iqr_vals = st_result["iqr"].dropna()
                if not ang_vals.empty:
                    # circular median (double-angle trick → stay in [-90°, 90°])
                    rad2 = np.radians(2.0 * ang_vals.to_numpy())
                    med2 = 0.5 * np.degrees(
                        np.arctan2(np.sin(rad2).mean(), np.cos(rad2).mean())
                    )
                    strike_consensus = float(med2)
                    strike_iqr = float(np.nanmedian(iqr_vals.to_numpy()))
        except Exception as exc:
            warnings.append(f"estimate_strike_consensus: {exc}")

        # ── figures ───────────────────────────────────────────────────────────

        # 1. Phase tensor pseudosection
        try:
            ax_pt = plot_phase_tensor_psection(
                sites,
                period_range=band,
                figsize=(10.0, 5.5),
                verbose=0,
            )
            fig = (
                ax_pt.get_figure() if hasattr(ax_pt, "get_figure") else ax_pt
            )
            figures["pt_psection"] = fig
            p = self._save_figure(
                fig, output_dir, "pt_psection", warnings_list=warnings
            )
            if p:
                fig_paths["pt_psection"] = p
        except Exception as exc:
            warnings.append(f"plot_phase_tensor_psection: {exc}")

        # 2. Phase tensor rose
        try:
            fig_rose = plot_phase_tensor_rose(sites, band=band, verbose=0)
            if fig_rose is not None:
                f = (
                    fig_rose
                    if hasattr(fig_rose, "savefig")
                    else (
                        fig_rose.get_figure()
                        if hasattr(fig_rose, "get_figure")
                        else None
                    )
                )
                if f is not None:
                    figures["pt_rose"] = f
                    p = self._save_figure(
                        f, output_dir, "pt_rose", warnings_list=warnings
                    )
                    if p:
                        fig_paths["pt_rose"] = p
        except Exception as exc:
            warnings.append(f"plot_phase_tensor_rose: {exc}")

        # 3. Strike analysis (3-rose panel)
        try:
            fig_strike = plot_strike_analysis(sites, band=band, verbose=0)
            if fig_strike is not None:
                f = (
                    fig_strike
                    if hasattr(fig_strike, "savefig")
                    else (
                        fig_strike.get_figure()
                        if hasattr(fig_strike, "get_figure")
                        else None
                    )
                )
                if f is not None:
                    figures["strike_analysis"] = f
                    p = self._save_figure(
                        f,
                        output_dir,
                        "strike_analysis",
                        warnings_list=warnings,
                    )
                    if p:
                        fig_paths["strike_analysis"] = p
        except Exception as exc:
            warnings.append(f"plot_strike_analysis: {exc}")

        # 4. Dimensionality confidence grid
        try:
            ax_dim = plot_dim_confidence_grid(
                sites,
                skew_th=self.skew_th,
                ellipt_th=self.ellipt_th,
                verbose=0,
            )
            fig = (
                ax_dim.get_figure()
                if hasattr(ax_dim, "get_figure")
                else ax_dim
            )
            if fig is not None:
                figures["dim_confidence"] = fig
                p = self._save_figure(
                    fig,
                    output_dir,
                    "dim_confidence_grid",
                    warnings_list=warnings,
                )
                if p:
                    fig_paths["dim_confidence"] = p
        except Exception as exc:
            warnings.append(f"plot_dim_confidence_grid: {exc}")

        # 5. Optional: survey fingerprint (station × period, 4 metrics)
        if run_fingerprint:
            try:
                from ..emtools.advanced import (
                    plot_survey_fingerprint,
                )

                fig_fp = plot_survey_fingerprint(
                    sites,
                    quantities=["skew", "ellipt", "theta", "s1"],
                    period_range=band,
                )
                figures["survey_fingerprint"] = fig_fp
                p = self._save_figure(
                    fig_fp,
                    output_dir,
                    "survey_fingerprint",
                    warnings_list=warnings,
                )
                if p:
                    fig_paths["survey_fingerprint"] = p
            except Exception as exc:
                warnings.append(f"plot_survey_fingerprint: {exc}")

        # 6. Optional: Mohr circles for first station
        if run_mohr:
            try:
                from ..emtools.advanced import (
                    plot_impedance_mohr_circles,
                )

                fig_mohr = plot_impedance_mohr_circles(
                    sites,
                    n_periods=8,
                    verbose=0,
                )
                figures["mohr_circles"] = fig_mohr
                p = self._save_figure(
                    fig_mohr,
                    output_dir,
                    "mohr_circles",
                    warnings_list=warnings,
                )
                if p:
                    fig_paths["mohr_circles"] = p
            except Exception as exc:
                warnings.append(f"plot_impedance_mohr_circles: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        n_obs = (n_1d + n_2d + n_3d) or 1
        if self.api_key:
            prompt = (
                f"Phase tensor analysis summary:\n"
                f"  1-D observations: {n_1d} ({100 * n_1d / n_obs:.0f}%)\n"
                f"  2-D observations: {n_2d} ({100 * n_2d / n_obs:.0f}%)\n"
                f"  3-D observations: {n_3d} ({100 * n_3d / n_obs:.0f}%)\n"
                f"  Consensus strike: {strike_consensus:.1f}°  "
                f"  IQR: {strike_iqr:.1f}°\n"
                f"  Skew threshold: {self.skew_th}°, "
                f"  ellipticity threshold: {self.ellipt_th}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Interpret the dimensionality and strike, and recommend "
                "whether to rotate the data before inversion."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        elapsed = time.time() - t0
        frac_3d = 100 * n_3d / n_obs
        return AgentResult(
            status="success",
            summary=(
                f"Phase analysis complete: "
                f"{n_1d}/{n_2d}/{n_3d} (1-D/2-D/3-D). "
                f"Strike {strike_consensus:.1f}° ± {strike_iqr:.1f}°. "
                f"3-D fraction {frac_3d:.0f}%. "
                f"{len(figures)} figure(s) produced."
            ),
            data={
                "pt_table": pt_table,
                "dim_table": dim_table,
                "n_1d": n_1d,
                "n_2d": n_2d,
                "n_3d": n_3d,
                "strike_consensus": float(strike_consensus),
                "strike_iqr": float(strike_iqr),
                "figures": figures,
                "figure_paths": fig_paths,
                "sites": sites,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["PhaseAnalysisAgent"]
