"""
pycsamt.agents.inversion_comparison
=====================================

:class:`InversionComparisonAgent` — Side-by-side comparison of two inversion
results.

Compares any two resistivity sections — e.g. AI-predicted vs Occam2D,
before/after regularisation change, or two different stations profiles — and
produces:

* A difference section (A − B) in log₁₀(Ω·m).
* Pearson correlation coefficient ρ between the two sections.
* RMSE between log₁₀ρ values.
* An LLM narrative describing which model fits better and where the two
  differ geologically.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in MT inversion model evaluation and comparison.
Given two inversion results, write 4-5 sentences that:
1. State the overall similarity (correlation, RMSE) between the two models.
2. Identify depth ranges or stations where the two models disagree most.
3. Discuss which model is more physically plausible and why.
4. Recommend which result to use for geological interpretation.
5. Suggest additional constraints (e.g. borehole, gravity) to discriminate them.
Reply in plain scientific English.
"""


class InversionComparisonAgent(BaseAgent):
    """Compare two resistivity inversion results.

    Each result is either an :class:`~pycsamt.agents._base.AgentResult`
    from an inversion agent, or a plain dict with keys ``pred_rho`` and
    ``depths_km``.

    Parameters
    ----------
    api_key, model, llm_provider : str

    Input keys
    ----------
    ``result_a`` : AgentResult or dict
        First model. Must contain ``"pred_rho"`` (n_stations × n_layers) or
        ``"predictions"`` dict {station: log₁₀ρ array}.
    ``result_b`` : AgentResult or dict
        Second model. Same structure as *result_a*.
    ``label_a`` : str — name for result_a (default ``"Model A"``)
    ``label_b`` : str — name for result_b (default ``"Model B"``)
    ``depths_km`` : ndarray, optional — shared depth axis (km)
    ``station_names`` : list[str], optional
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``correlation``     float — Pearson ρ between log₁₀ρ sections
    ``rmse``            float — RMSE in log₁₀(Ω·m)
    ``difference``      ndarray (n_layers, n_stations) — A − B
    ``label_a``         str
    ``label_b``         str
    ``figures``         dict
    ``figure_paths``    dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
    ) -> None:
        super().__init__(
            "InversionComparisonAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        result_a = input_data.get("result_a")
        result_b = input_data.get("result_b")
        if result_a is None or result_b is None:
            return AgentResult.failed(
                "Both 'result_a' and 'result_b' are required.",
                hint="Pass AgentResult objects from AIInversionAgent, "
                "Inv2DAgent, or any dict with 'pred_rho'.",
                elapsed=time.time() - t0,
            )

        label_a = str(input_data.get("label_a", "Model A"))
        label_b = str(input_data.get("label_b", "Model B"))
        output_dir = input_data.get("output_dir")

        # ── extract sections ──────────────────────────────────────────────
        mat_a, stations_a, depths_a = _extract_section(
            result_a, warnings, label_a
        )
        mat_b, stations_b, depths_b = _extract_section(
            result_b, warnings, label_b
        )

        if mat_a is None:
            return AgentResult.failed(
                f"Could not extract resistivity section from {label_a!r}.",
                elapsed=time.time() - t0,
            )
        if mat_b is None:
            return AgentResult.failed(
                f"Could not extract resistivity section from {label_b!r}.",
                elapsed=time.time() - t0,
            )

        # ── align shapes ──────────────────────────────────────────────────
        n_layers = min(mat_a.shape[0], mat_b.shape[0])
        n_sta = min(mat_a.shape[1], mat_b.shape[1])
        mat_a = mat_a[:n_layers, :n_sta]
        mat_b = mat_b[:n_layers, :n_sta]
        if mat_a.shape != mat_b.shape:
            warnings.append(
                f"Section shapes differ after alignment: "
                f"{mat_a.shape} vs {mat_b.shape}."
            )

        # ── depth axis ────────────────────────────────────────────────────
        depths_km = input_data.get("depths_km")
        if depths_km is None:
            depths_km = depths_a if depths_a is not None else depths_b
        if depths_km is None:
            depths_km = np.linspace(0, 2.0, n_layers + 1)

        # ── station names ─────────────────────────────────────────────────
        station_names = input_data.get("station_names")
        if station_names is None:
            station_names = (
                stations_a or stations_b or [f"S{i}" for i in range(n_sta)]
            )
        station_names = list(station_names)[:n_sta]

        # ── statistics ────────────────────────────────────────────────────
        flat_a = mat_a[np.isfinite(mat_a) & np.isfinite(mat_b)]
        flat_b = mat_b[np.isfinite(mat_a) & np.isfinite(mat_b)]

        correlation = float("nan")
        rmse = float("nan")
        if flat_a.size > 1:
            try:
                correlation = float(np.corrcoef(flat_a, flat_b)[0, 1])
                rmse = float(np.sqrt(np.mean((flat_a - flat_b) ** 2)))
            except Exception:
                pass

        difference = mat_a - mat_b  # (n_layers, n_sta)

        # ── figures ───────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            fig = _plot_comparison(
                mat_a,
                mat_b,
                difference,
                depths_km,
                station_names,
                label_a,
                label_b,
                correlation,
                rmse,
            )
            if fig is not None:
                figures["comparison"] = fig
                p = self._save_figure(
                    fig,
                    output_dir,
                    "inversion_comparison",
                    warnings_list=warnings,
                )
                if p:
                    fig_paths["comparison"] = p
        except Exception as exc:
            warnings.append(f"Comparison figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────
        interp: str | None = None
        if self.api_key:
            max_diff_idx = np.unravel_index(
                np.nanargmax(np.abs(difference)), difference.shape
            )
            prompt = (
                f"Inversion comparison:\n"
                f"  Model A: {label_a}  |  Model B: {label_b}\n"
                f"  Section size: {n_layers} layers × {n_sta} stations\n"
                f"  Pearson correlation: {correlation:.4f}\n"
                f"  RMSE: {rmse:.4f} log₁₀(Ω·m)\n"
                f"  Max difference at layer {max_diff_idx[0] + 1}, "
                f"station {station_names[max_diff_idx[1]] if max_diff_idx[1] < len(station_names) else '?'}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                f"Compare {label_a} and {label_b} and recommend which to use."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        elapsed = time.time() - t0
        corr_str = (
            f"{correlation:.4f}" if not np.isnan(correlation) else "N/A"
        )
        rmse_str = f"{rmse:.4f}" if not np.isnan(rmse) else "N/A"
        return AgentResult(
            status="success",
            summary=(
                f"Inversion comparison ({label_a} vs {label_b}): "
                f"Pearson ρ = {corr_str}, RMSE = {rmse_str} log₁₀(Ω·m). "
                f"{len(figures)} figures."
            ),
            data={
                "correlation": correlation,
                "rmse": rmse,
                "difference": difference,
                "label_a": label_a,
                "label_b": label_b,
                "section_a": mat_a,
                "section_b": mat_b,
                "depths_km": depths_km,
                "station_names": station_names,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────


def _extract_section(
    result: Any,
    warnings: list[str],
    label: str,
) -> tuple[np.ndarray | None, list[str] | None, np.ndarray | None]:
    """Return (mat (n_layers, n_sta), station_names, depths_km) or (None, ...)."""
    if result is None:
        return None, None, None

    # dict/AgentResult accessors
    get = (
        result.get
        if hasattr(result, "get")
        else lambda k, d=None: result.get(k, d)
    )

    # ── 2-D section: (n_layers, n_stations) ──────────────────────────────
    for key in (
        "pred_rho",
        "pred_section",
        "section_a",
        "section_b",
        "resistivity_section",
    ):
        val = get(key)
        if val is not None:
            mat = np.asarray(val, dtype=float)
            if mat.ndim == 2:
                depths_km = _get_depths(get)
                stations = _get_stations(get)
                return mat, stations, depths_km

    # ── 1-D predictions: dict {station: log_rho array} ───────────────────
    preds = get("predictions")
    if preds and isinstance(preds, dict):
        station_names = list(preds.keys())
        arrays = [np.asarray(preds[s], dtype=float) for s in station_names]
        n_layers = max(len(a) for a in arrays)
        mat = np.full((n_layers, len(station_names)), np.nan)
        for si, a in enumerate(arrays):
            n = min(len(a), n_layers)
            mat[:n, si] = a[:n]
        depths_km = _get_depths(get)
        return mat, station_names, depths_km

    warnings.append(f"{label}: no recognisable section key found.")
    return None, None, None


def _get_depths(get) -> np.ndarray | None:
    for key in ("depths_km", "depths", "depth_axis"):
        val = get(key)
        if val is not None:
            return np.asarray(val, dtype=float)
    return None


def _get_stations(get) -> list[str] | None:
    val = get("station_names") or get("stations")
    return list(val) if val is not None else None


def _plot_comparison(
    mat_a: np.ndarray,
    mat_b: np.ndarray,
    diff: np.ndarray,
    depths_km: np.ndarray,
    station_names: list[str],
    label_a: str,
    label_b: str,
    corr: float,
    rmse: float,
) -> Any:
    import matplotlib.pyplot as plt

    from ..api.section import PYCSAMT_SECTION
    from ..api.station import PYCSAMT_STATION_RENDERING

    n_layers, n_sta = mat_a.shape
    section = PYCSAMT_SECTION.style_for("inversion")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    d = np.asarray(depths_km)
    if d.size < n_layers + 1:
        d = np.linspace(0, float(d.max()) if d.size else 2.0, n_layers + 1)
    extent = (-0.5, n_sta - 0.5, float(d[-1]), float(d[0]))

    vv = np.concatenate(
        [mat_a[np.isfinite(mat_a)], mat_b[np.isfinite(mat_b)]]
    )
    vmin = float(np.percentile(vv, 5)) if vv.size else 0.0
    vmax = float(np.percentile(vv, 95)) if vv.size else 4.0

    for ax, mat, lbl in [
        (axes[0], mat_a, label_a),
        (axes[1], mat_b, label_b),
    ]:
        im = ax.imshow(
            mat,
            aspect="auto",
            origin="upper",
            extent=extent,
            cmap="jet_r",
            vmin=vmin,
            vmax=vmax,
            interpolation="bilinear",
        )
        PYCSAMT_STATION_RENDERING.apply(
            ax,
            np.arange(n_sta, dtype=float),
            station_names,
            preset="inversion",
            xlim=(-0.5, n_sta - 0.5),
        )
        ax.set_ylabel("Depth (km)", fontsize=9)
        ax.tick_params(labelsize=8)
        section.add_colorbar(im, ax, label="log₁₀ρ (Ω·m)")
        ax.set_title(lbl, fontsize=9, fontweight="bold")

    # difference panel
    dd = diff[np.isfinite(diff)]
    d_abs = float(np.percentile(np.abs(dd), 95)) if dd.size else 1.0
    im3 = axes[2].imshow(
        diff,
        aspect="auto",
        origin="upper",
        extent=extent,
        cmap="RdBu_r",
        vmin=-d_abs,
        vmax=d_abs,
        interpolation="bilinear",
    )
    PYCSAMT_STATION_RENDERING.apply(
        axes[2],
        np.arange(n_sta, dtype=float),
        station_names,
        preset="inversion",
        xlim=(-0.5, n_sta - 0.5),
    )
    axes[2].set_ylabel("Depth (km)", fontsize=9)
    axes[2].tick_params(labelsize=8)
    section.add_colorbar(im3, axes[2], label="Δlog₁₀ρ (A−B)")
    corr_str = f"{corr:.4f}" if not np.isnan(corr) else "N/A"
    rmse_str = f"{rmse:.4f}" if not np.isnan(rmse) else "N/A"
    axes[2].set_title(
        f"Difference  (ρ={corr_str}, RMSE={rmse_str})",
        fontsize=9,
        fontweight="bold",
    )

    fig.suptitle(
        f"Inversion comparison: {label_a}  vs  {label_b}",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout()
    return fig


__all__ = ["InversionComparisonAgent"]
