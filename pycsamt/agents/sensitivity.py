"""
pycsamt.agents.sensitivity
============================

:class:`SensitivityAgent` — Compute vertical resolution and Bostick
sensitivity kernels for MT data.

Answers "which depth range does each datum actually constrain?" using the
Bostick penetration depth and the analytical vertical-resolution formula.
Produces a sensitivity pseudosection that makes the depth of investigation
(DOI) visually explicit — information that is absent from every standard
pseudosection plot.

Wraps
-----
* :func:`~pycsamt.emtools.csumt.vertical_resolution` — per-station ΔD table
* :func:`~pycsamt.emtools.advanced.plot_sensitivity_depth_section` — figure
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in MT data resolution and depth-of-investigation analysis.
Given a sensitivity analysis result, write 4-5 sentences that:
1. Describe the maximum Bostick depth reached by the data at the lowest frequency.
2. Identify depth ranges with poor resolution (large ΔD) and their cause.
3. State which stations have the best and worst overall depth coverage.
4. Advise on whether the target depth is within the data's sensitivity window.
5. Recommend frequency additions or deletions to improve depth resolution.
Reply in plain scientific English.
"""


class SensitivityAgent(BaseAgent):
    """Bostick sensitivity kernels and vertical resolution analysis.

    Parameters
    ----------
    api_key, model, llm_provider : str
    component : {'xy', 'yx'}
        Impedance component for ρa and Bostick depth (default ``'xy'``).
    rho_override : float or None
        Fixed background resistivity (Ω·m) for the analytical ΔD formula.
        ``None`` uses the measured ρa per frequency pair.
    depth_max : float or None
        Clip depth axis at this value in km (default: auto from data).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``component`` : str, optional
    ``rho_override`` : float, optional
    ``period_range`` : [T_min, T_max], optional
    ``depth_max`` : float, optional — km
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``resolution_table``  pandas.DataFrame — per-(station, freq-pair)
    ``doi_per_station``   dict {station: float}  — max Bostick depth (m)
    ``mean_doi_km``       float
    ``figures``           dict
    ``figure_paths``      dict

    Examples
    --------
    >>> agent = SensitivityAgent(component='xy')
    >>> r = agent.execute({"path": "/data/WILLY_EDIs"})
    >>> print(r["mean_doi_km"], "km mean DOI")
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        component: str = "xy",
        rho_override: float | None = None,
        depth_max: float | None = None,
    ) -> None:
        super().__init__(
            "SensitivityAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.component    = component.lower()
        self.rho_override = rho_override
        self.depth_max    = depth_max

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        try:
            from ..emtools.csumt   import vertical_resolution
            from ..emtools.advanced import plot_sensitivity_depth_section
        except ImportError as exc:
            return AgentResult.failed(
                f"emtools.csumt / emtools.advanced not available: {exc}",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import ensure_sites

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed("No 'sites' or 'path'.", elapsed=time.time() - t0)
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        component    = str(input_data.get("component",    self.component))
        rho_override = input_data.get("rho_override", self.rho_override)
        period_range = input_data.get("period_range")
        depth_max    = input_data.get("depth_max", self.depth_max)
        output_dir   = input_data.get("output_dir")

        # ── compute resolution table ──────────────────────────────────────
        try:
            res_table = vertical_resolution(
                sites,
                rho_override=rho_override,
                verbose=0,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"vertical_resolution failed: {exc}", elapsed=time.time() - t0,
            )

        # ── DOI per station (max Bostick depth_lo) ────────────────────────
        doi_per_station: dict[str, float] = {}
        try:
            if hasattr(res_table, "groupby"):
                doi_per_station = (
                    res_table.groupby("station")["depth_lo_m"].max().to_dict()
                )
        except Exception:
            pass

        mean_doi_km = float("nan")
        if doi_per_station:
            mean_doi_km = float(np.mean(list(doi_per_station.values()))) / 1000.0

        # ── figures ───────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            d_max_kw: dict = {}
            if depth_max is not None:
                d_max_kw = {"depth_max": depth_max}
            p_range_kw: dict = {}
            if period_range is not None:
                p_range_kw = {"period_range": tuple(period_range)}

            fig_sens = plot_sensitivity_depth_section(
                sites,
                component=component,
                depth_unit="km",
                **d_max_kw,
                **p_range_kw,
            )
            if fig_sens is not None:
                figures["sensitivity_section"] = fig_sens
                p = self._save_figure(fig_sens, output_dir,
                                      "sensitivity_depth_section",
                                      warnings_list=warnings)
                if p:
                    fig_paths["sensitivity_section"] = p
        except Exception as exc:
            warnings.append(f"plot_sensitivity_depth_section: {exc}")

        # resolution per station bar chart
        try:
            fig_doi = _plot_doi_bar(doi_per_station)
            if fig_doi is not None:
                figures["doi_bar"] = fig_doi
                p = self._save_figure(fig_doi, output_dir, "doi_per_station",
                                      warnings_list=warnings)
                if p:
                    fig_paths["doi_bar"] = p
        except Exception as exc:
            warnings.append(f"DOI bar figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and doi_per_station:
            max_doi_sta = max(doi_per_station, key=doi_per_station.get)
            min_doi_sta = min(doi_per_station, key=doi_per_station.get)
            prompt = (
                f"Sensitivity analysis summary:\n"
                f"  Component: Z{component}\n"
                f"  Stations analysed: {len(doi_per_station)}\n"
                f"  Mean DOI: {mean_doi_km:.2f} km\n"
                f"  Max DOI: {doi_per_station[max_doi_sta]/1000:.2f} km "
                f"at station {max_doi_sta}\n"
                f"  Min DOI: {doi_per_station[min_doi_sta]/1000:.2f} km "
                f"at station {min_doi_sta}\n"
                f"  ρ override: {rho_override or 'from data'}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Interpret the depth of investigation and resolution quality."
            )
            interp = self.query_llm(prompt, max_tokens=220)

        elapsed = time.time() - t0
        doi_str = f"{mean_doi_km:.2f} km" if not np.isnan(mean_doi_km) else "N/A"
        return AgentResult(
            status="success",
            summary=(
                f"Sensitivity analysis: {len(doi_per_station)} stations. "
                f"Mean DOI: {doi_str}. {len(figures)} figures."
            ),
            data={
                "resolution_table": res_table,
                "doi_per_station":  doi_per_station,
                "mean_doi_km":      mean_doi_km,
                "component":        component,
                "figures":          figures,
                "figure_paths":     fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────

def _plot_doi_bar(doi_per_station: dict[str, float]) -> Any:
    """Horizontal bar chart of DOI (km) per station."""
    if not doi_per_station:
        return None
    import matplotlib.pyplot as plt

    stations = list(doi_per_station.keys())
    depths   = [doi_per_station[s] / 1000.0 for s in stations]

    fig, ax = plt.subplots(figsize=(6, max(3, len(stations) * 0.35)))
    bars = ax.barh(range(len(stations)), depths, color="#2980b9", alpha=0.8,
                   edgecolor="white", linewidth=0.4)
    ax.set_yticks(range(len(stations)))
    ax.set_yticklabels(stations, fontsize=7)
    ax.set_xlabel("Depth of investigation (km)", fontsize=9)
    ax.set_title("Max Bostick DOI per station", fontsize=9, fontweight="bold")
    ax.tick_params(labelsize=8)
    ax.invert_yaxis()
    fig.tight_layout()
    return fig


__all__ = ["SensitivityAgent"]
