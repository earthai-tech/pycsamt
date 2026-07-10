"""
pycsamt.agents.resistivity_map
================================

:class:`ResistivityMapAgent` — Horizontal depth-slice resistivity maps.

Assembles pseudo-3D resistivity volumes from per-station 1-D (or 2-D profile)
inversion results.  At each requested depth, per-station log₁₀ρ values are
scattered onto a 2-D grid and interpolated with linear or IDW methods to
produce a plan-view map.

Typical use: feed the ``predictions`` output of
:class:`~pycsamt.agents.AIInversionAgent`, :class:`~pycsamt.agents.Inv3DAgent`,
or any compatible inversion result together with station coordinates.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in pseudo-3D resistivity interpretation from MT surveys.
Given a set of horizontal resistivity maps, write 4-5 sentences that:
1. Describe the dominant resistivity pattern at each depth level.
2. Identify lateral contrasts that suggest geological boundaries or structures.
3. Note any stations that appear anomalous at specific depths.
4. Discuss the reliability of the interpolation given station spacing.
5. Recommend drilling targets or geological follow-up based on the maps.
Reply in plain scientific English.
"""


class ResistivityMapAgent(BaseAgent):
    """Build horizontal resistivity depth-slice maps from 1-D inversion results.

    Parameters
    ----------
    api_key, model, llm_provider : str
    depth_indices : list[int] or None
        Layer indices to map (0-based).  ``None`` maps 3 evenly spaced layers.
    interp_method : {'linear', 'nearest', 'idw'}
        Interpolation method for gridding station values (default ``'linear'``).
    grid_n : int
        Number of grid cells per axis (default 50).

    Input keys
    ----------
    ``predictions`` : dict {station: ndarray} — log₁₀ρ per layer
        (from AIInversionAgent, Inv3DAgent, etc.)
    ``station_coords`` : dict {station: (x, y)} or ndarray (n_sta, 2) — metres
    ``depths_km`` : ndarray — depth axis (km)
    ``depth_indices`` : list[int], optional
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``depth_maps``       list[dict] — per-depth: {'depth_km', 'grid_x', 'grid_y', 'grid_rho'}
    ``depth_levels_km``  list[float]
    ``figures``          dict
    ``figure_paths``     dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        depth_indices: list[int] | None = None,
        interp_method: str = "linear",
        grid_n: int = 50,
    ) -> None:
        super().__init__(
            "ResistivityMapAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.depth_indices = depth_indices
        self.interp_method = interp_method.lower()
        self.grid_n        = grid_n

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        predictions = input_data.get("predictions")
        if not predictions:
            return AgentResult.failed(
                "No 'predictions' dict provided.",
                hint="Pass predictions from AIInversionAgent or Inv3DAgent.",
                elapsed=time.time() - t0,
            )

        coords_raw = input_data.get("station_coords")
        if coords_raw is None:
            return AgentResult.failed(
                "No 'station_coords' provided.",
                hint="Pass {'StationName': (x_m, y_m), ...} or ndarray (n_sta, 2).",
                elapsed=time.time() - t0,
            )

        depths_raw   = input_data.get("depths_km")
        depth_indices = input_data.get("depth_indices") or self.depth_indices
        output_dir    = input_data.get("output_dir")

        # ── parse coordinates ─────────────────────────────────────────────
        station_names = list(predictions.keys())
        n_sta = len(station_names)

        if isinstance(coords_raw, dict):
            xy = np.array([coords_raw[nm] for nm in station_names
                           if nm in coords_raw], dtype=float)
            present = [nm for nm in station_names if nm in coords_raw]
        else:
            xy = np.asarray(coords_raw, dtype=float)[:n_sta]
            present = station_names[:len(xy)]

        if len(xy) < 3:
            return AgentResult.failed(
                f"Need ≥ 3 stations with coordinates, got {len(xy)}.",
                elapsed=time.time() - t0,
            )

        # ── build section matrix ──────────────────────────────────────────
        n_layers = max(len(np.asarray(predictions[nm])) for nm in present)
        mat = np.full((n_sta, n_layers), np.nan)
        for si, nm in enumerate(present):
            v = np.asarray(predictions[nm], dtype=float)
            n = min(len(v), n_layers)
            mat[si, :n] = v[:n]

        # ── depth axis ────────────────────────────────────────────────────
        if depths_raw is not None:
            depths_km = np.asarray(depths_raw, dtype=float)
        else:
            depths_km = np.linspace(0.0, 2.0, n_layers + 1)

        # ── choose depth indices to map ───────────────────────────────────
        if depth_indices is None:
            depth_indices = list(
                np.linspace(0, n_layers - 1, min(3, n_layers), dtype=int)
            )
        depth_indices = [int(d) for d in depth_indices if 0 <= int(d) < n_layers]

        if not depth_indices:
            return AgentResult.failed(
                "No valid depth indices to map.", elapsed=time.time() - t0,
            )

        # ── build grid and interpolate ────────────────────────────────────
        cx, cy = xy[:len(present), 0], xy[:len(present), 1]
        xg = np.linspace(cx.min(), cx.max(), self.grid_n)
        yg = np.linspace(cy.min(), cy.max(), self.grid_n)
        Xg, Yg = np.meshgrid(xg, yg)

        depth_maps: list[dict] = []
        for li in depth_indices:
            vals = mat[:len(present), li]
            mask = np.isfinite(vals)
            if mask.sum() < 3:
                warnings.append(f"Layer {li}: insufficient valid data for interpolation.")
                continue
            try:
                rho_grid = _interpolate_map(
                    cx[mask], cy[mask], vals[mask],
                    Xg, Yg, method=self.interp_method,
                )
            except Exception as exc:
                warnings.append(f"Layer {li} interpolation failed: {exc}")
                continue

            d_km = float(depths_km[min(li, len(depths_km) - 1)])
            depth_maps.append({
                "layer_idx":  li,
                "depth_km":   d_km,
                "grid_x":     xg,
                "grid_y":     yg,
                "grid_rho":   rho_grid,
                "station_x":  cx,
                "station_y":  cy,
                "station_rho": vals,
            })

        if not depth_maps:
            return AgentResult.failed(
                "All depth layers failed to interpolate.", elapsed=time.time() - t0,
            )

        # ── figures ───────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            fig = _plot_depth_maps(depth_maps, present)
            if fig is not None:
                figures["depth_maps"] = fig
                p = self._save_figure(fig, output_dir, "resistivity_depth_maps",
                                      warnings_list=warnings)
                if p:
                    fig_paths["depth_maps"] = p
        except Exception as exc:
            warnings.append(f"Depth map figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────
        interp: str | None = None
        if self.api_key:
            depth_levels = [dm["depth_km"] for dm in depth_maps]
            mean_rhos = [float(np.nanmean(dm["grid_rho"])) for dm in depth_maps]
            prompt = (
                f"Resistivity depth maps:\n"
                f"  Stations: {len(present)}\n"
                f"  Depth levels mapped (km): {depth_levels}\n"
                f"  Mean log₁₀ρ per level: {[f'{v:.2f}' for v in mean_rhos]}\n"
                f"  Grid: {self.grid_n}×{self.grid_n}, method: {self.interp_method}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Interpret the horizontal resistivity maps and identify structures."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        elapsed = time.time() - t0
        return AgentResult(
            status="success",
            summary=(
                "Resistivity maps: {} depth levels ({} km). {} figures.".format(
                    len(depth_maps),
                    ", ".join(f"{dm['depth_km']:.2f}" for dm in depth_maps),
                    len(figures),
                )
            ),
            data={
                "depth_maps":       depth_maps,
                "depth_levels_km":  [dm["depth_km"] for dm in depth_maps],
                "station_names":    present,
                "station_coords":   xy,
                "figures":          figures,
                "figure_paths":     fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────

def _interpolate_map(
    cx: np.ndarray,
    cy: np.ndarray,
    vals: np.ndarray,
    Xg: np.ndarray,
    Yg: np.ndarray,
    method: str,
) -> np.ndarray:
    """Interpolate scattered station values onto a regular grid."""
    from scipy.interpolate import (
        NearestNDInterpolator,
        griddata,
    )

    points = np.column_stack([cx, cy])
    xi     = np.column_stack([Xg.ravel(), Yg.ravel()])

    if method == "idw":
        # Inverse-distance weighting (power = 2)
        dists = np.sqrt(((xi[:, None, :] - points[None, :, :]) ** 2).sum(-1))
        dists = np.maximum(dists, 1.0)  # avoid division by zero
        w     = 1.0 / dists ** 2
        rho_flat = (w * vals[None, :]).sum(-1) / w.sum(-1)
    else:
        rho_flat = griddata(points, vals, xi,
                            method="linear" if method == "linear" else "nearest")
        # fill NaN borders with nearest
        nan_mask = ~np.isfinite(rho_flat)
        if nan_mask.any():
            nn = NearestNDInterpolator(points, vals)
            rho_flat[nan_mask] = nn(xi[nan_mask])

    return rho_flat.reshape(Xg.shape)


def _plot_depth_maps(depth_maps: list[dict], station_names: list[str]) -> Any:
    """Panel of horizontal depth-slice maps."""
    import matplotlib.pyplot as plt

    n = len(depth_maps)
    if n == 0:
        return None

    fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 5.0))
    if n == 1:
        axes = [axes]

    all_vals = np.concatenate([dm["grid_rho"][np.isfinite(dm["grid_rho"])]
                                for dm in depth_maps])
    vmin = float(np.percentile(all_vals, 5))  if all_vals.size else 0.0
    vmax = float(np.percentile(all_vals, 95)) if all_vals.size else 4.0

    for ax, dm in zip(axes, depth_maps):
        xg  = dm["grid_x"] / 1000.0   # km
        yg  = dm["grid_y"] / 1000.0
        Xg, Yg = np.meshgrid(xg, yg)

        im = ax.pcolormesh(
            Xg, Yg, dm["grid_rho"],
            cmap="jet_r", vmin=vmin, vmax=vmax, shading="auto",
        )
        plt.colorbar(im, ax=ax, label="log₁₀ρ (Ω·m)", shrink=0.85)

        # station markers
        st_x = dm["station_x"] / 1000.0
        st_y = dm["station_y"] / 1000.0
        ax.scatter(st_x, st_y, c=dm["station_rho"], cmap="jet_r",
                   vmin=vmin, vmax=vmax, edgecolors="k", lw=0.5, s=50, zorder=5)
        for xi, yi, nm in zip(st_x, st_y, station_names):
            ax.text(xi, yi, nm, fontsize=5.5, ha="center", va="bottom", color="w",
                    fontweight="bold")

        ax.set_xlabel("E–W (km)", fontsize=8)
        ax.set_ylabel("N–S (km)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_title(f"Depth  {dm['depth_km']:.2f} km", fontsize=8, fontweight="bold")
        ax.set_aspect("equal")

    fig.suptitle("Pseudo-3D resistivity depth maps", fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig


__all__ = ["ResistivityMapAgent"]
