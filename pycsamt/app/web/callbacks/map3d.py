# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""callbacks/map3d.py — 3D resistivity map: fence, volume-block, depth-slices."""
from __future__ import annotations

import io
import base64

import numpy as np
import plotly.graph_objects as go

from dash import (Input, Output, State, ctx, no_update, dcc,
                  clientside_callback, html)
from dash.exceptions import PreventUpdate

from pycsamt.app.web.layout import IDs
from pycsamt.app.web.cache import cache_get, cache_get_inversion_result
from pycsamt.app.web.pages.map3d import _empty_3d_fig, _MAP3D_MODES, _DEFAULT_MODE
from pycsamt.map.geometry import normalize_offsets, resolve_offset, survey_uv


# ── Mode metadata ─────────────────────────────────────────────────────────────
_MODE_SLUGS = [slug for slug, _, _ in _MAP3D_MODES]
_MODE_ICON  = {slug: icon for slug, _, icon in _MAP3D_MODES}
_MODE_LABEL = {slug: lbl  for slug, lbl, _ in _MAP3D_MODES}

# ρ preset ranges (Ω·m)
_RHO_PRESETS = {
    "map3d-rho-preset-all":  (1,      100_000),
    "map3d-rho-preset-cond": (1,      100),
    "map3d-rho-preset-mid":  (100,    1_000),
    "map3d-rho-preset-res":  (1_000,  100_000),
}
# Depth preset ranges (m)
_DEPTH_PRESETS = {
    "map3d-depth-preset-full": (0, 5_000),
    "map3d-depth-preset-500":  (0,   500),
    "map3d-depth-preset-1k":   (0, 1_000),
    "map3d-depth-preset-2k":   (0, 2_000),
}


# ── Plotly scene themes ───────────────────────────────────────────────────────

def _scene_colors(dark: bool = True) -> dict:
    if dark:
        return {
            "scene": "#181825",
            "paper": "#1e1e2e",
            "grid": "#313244",
            "zero": "#45475a",
            "tick": "#a6adc8",
            "text": "#cdd6f4",
            "border": "#313244",
            "warn": "#f9e2af",
        }
    return {
        "scene": "#ffffff",
        "paper": "#eff1f5",
        "grid": "#dce0e8",
        "zero": "#bcc0cc",
        "tick": "#6c6f85",
        "text": "#4c4f69",
        "border": "#ccd0da",
        "warn": "#df8e1d",
    }


def _dark_scene(title: str = "", dark: bool = True) -> dict:
    colors = _scene_colors(dark)
    axis = dict(
        backgroundcolor=colors["scene"],
        gridcolor=colors["grid"],
        zerolinecolor=colors["zero"],
        tickfont=dict(color=colors["tick"], size=9),
        title=dict(font=dict(color=colors["text"], size=10)),
    )
    _t = lambda s: dict(text=s, font=dict(color=colors["text"], size=10))
    return dict(
        scene=dict(
            xaxis={**axis, "title": _t("Distance (m)")},
            yaxis={**axis, "title": _t("Profile offset (m)")},
            zaxis={**axis, "title": _t("Depth (m)")},
            bgcolor=colors["scene"],
            camera=dict(eye=dict(x=1.6, y=-1.6, z=0.9)),
            aspectmode="data",
        ),
        paper_bgcolor=colors["paper"],
        plot_bgcolor=colors["paper"],
        coloraxis_colorbar=dict(
            tickfont=dict(color=colors["tick"], size=9),
            title=dict(text="log₁₀ρ (Ω·m)", side="right",
                       font=dict(color=colors["text"], size=10)),
            bgcolor=colors["scene"],
            bordercolor=colors["border"],
        ),
        margin=dict(l=0, r=0, t=30 if title else 0, b=0),
        title=dict(text=title, font=dict(color=colors["text"], size=13),
                   x=0.5) if title else {},
    )


# ── ρ masking helper ──────────────────────────────────────────────────────────

def _apply_rho_mask(rho: np.ndarray, rho_lo: float, rho_hi: float) -> np.ndarray:
    """Return rho with values outside [rho_lo, rho_hi] set to NaN (transparent)."""
    if rho_lo <= 0 and rho_hi >= 1e8:
        return rho  # no filtering needed
    masked = rho.astype(float).copy()
    masked[(masked < rho_lo) | (masked > rho_hi)] = np.nan
    return masked


def _apply_depth_mask(z_arr: np.ndarray, rho_grid: np.ndarray,
                      depth_lo: float, depth_hi: float):
    """Clip z_arr to [depth_lo, depth_hi] and slice rho_grid accordingly."""
    z = np.abs(z_arr)
    mask = (z >= depth_lo) & (z <= depth_hi)
    if not mask.any():
        return z_arr, rho_grid
    z_clip   = z_arr[mask]
    # rho_grid is (n_z, n_x) — first axis is depth
    if rho_grid.ndim == 2 and rho_grid.shape[0] == len(z_arr):
        rho_clip = rho_grid[mask, :]
    elif rho_grid.ndim == 2 and rho_grid.shape[1] == len(z_arr):
        rho_clip = rho_grid[:, mask]
    else:
        rho_clip = rho_grid
    return z_clip, rho_clip


def _axis_fraction_mask(values: np.ndarray, fraction: float) -> np.ndarray:
    """Return a stable left-to-right clipping mask for any numeric axis."""
    arr = np.asarray(values, float)
    if arr.size == 0:
        return np.zeros_like(arr, dtype=bool)
    if fraction >= 0.999:
        return np.ones_like(arr, dtype=bool)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return np.zeros_like(arr, dtype=bool)
    lo = float(finite.min())
    hi = float(finite.max())
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.ones_like(arr, dtype=bool)
    stop = lo + (hi - lo) * max(0.0, min(1.0, float(fraction)))
    return arr <= stop


def _interp_finite_1d(x_new: np.ndarray, x: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Interpolate a 1-D vector without letting NaNs erase valid neighbours."""
    x = np.asarray(x, float)
    values = np.asarray(values, float)
    good = np.isfinite(x) & np.isfinite(values)
    if good.sum() >= 2:
        return np.interp(x_new, x[good], values[good])
    if good.sum() == 1:
        return np.full_like(x_new, values[good][0], dtype=float)
    return np.full_like(x_new, np.nan, dtype=float)


def _line_real_offsets(profiles: dict) -> dict | None:
    """Cross-strike offset (m) for each line, from real station
    lat/lon — see :mod:`pycsamt.map.geometry` (shared with the
    mapview 3-D builder). Requires every line to carry ``sta_lat``/
    ``sta_lon`` (only populated for the "pseudo" data source, not
    "inversion" — see ``_profiles_from_inversion_result``); returns
    ``None`` otherwise so callers fall back to the synthetic
    index-based line stack.
    """
    ids: list = []
    lats: list = []
    lons: list = []
    lines: list = []
    per_line_ids: dict[str, list] = {}
    for name, p in profiles.items():
        lat = p.get("sta_lat") or []
        lon = p.get("sta_lon") or []
        names = p.get("sta_names") or []
        if not lat or not lon or len(lat) != len(lon) or len(lat) != len(names):
            return None
        line_ids = [f"{name}::{n}" for n in names]
        per_line_ids[name] = line_ids
        ids.extend(line_ids)
        lats.extend(lat)
        lons.extend(lon)
        lines.extend([name] * len(line_ids))
    if len(ids) < 2:
        return None
    uv = survey_uv(ids, lats, lons, lines)
    if not uv:
        return None
    raw = {}
    for name, line_ids in per_line_ids.items():
        vs = [uv[i][1] for i in line_ids if i in uv]
        if not vs:
            return None
        raw[name] = float(np.median(vs))
    return normalize_offsets(raw)


# ── Figure builders ───────────────────────────────────────────────────────────

def _build_fence_fig(
    profiles: dict,
    line_spacing: float, azimuth_deg: float, panel_thick: float,
    cmap: str, vmin: float, vmax: float,
    opacity: float, contours: bool, n_contours: int,
    show_labels: bool, title: str,
    rho_lo: float = 1, rho_hi: float = 1e8,
    depth_lo: float = 0, depth_hi: float = 1e6,
    dark: bool = True,
    topo_src: str = "none",
    topo_opacity: float = 0.70,
    show_sta_markers: bool = True,
    apply_topo: bool = False,
    sta_symbol: str = "auto",
    sta_size: int = 8,
    sta_color: str = "black",
) -> go.Figure:
    traces = []
    annotations = []
    line_names = list(profiles.keys())
    n_lines = len(line_names)
    real_offsets = _line_real_offsets(profiles)
    az_rad = np.deg2rad(azimuth_deg)
    colors = _scene_colors(dark)
    use_topo = topo_src != "none"

    for i, name in enumerate(line_names):
        p = profiles[name]
        x_arr = np.asarray(p["x"])
        z_raw = np.asarray(p["z"])
        rho   = np.asarray(p["rho"])

        # Topo elevation for this line (Phase 3: source-aware resolution)
        sta_x          = np.asarray(p.get("sta_x", x_arr), float)
        sta_elev       = _resolve_profile_elevs(p, topo_src) if use_topo else np.zeros(len(sta_x))
        sta_names_line = p.get("sta_names", [f"S{j}" for j in range(len(sta_x))])
        elev_at_x      = _elev_along_line(x_arr, sta_x, sta_elev) if use_topo else None

        # Universal depth filter
        z_arr, rho = _apply_depth_mask(z_raw, rho, depth_lo, depth_hi)

        # Universal ρ mask (NaN → transparent gaps)
        rho = _apply_rho_mask(rho, rho_lo, rho_hi)

        # Handle shape: expect (n_z, n_x)
        if rho.ndim == 2 and rho.shape[0] == len(x_arr):
            rho = rho.T   # ensure (n_z, n_x)

        offset_m = resolve_offset(name, i, real_offsets, 1000.0, line_spacing)
        y_center = float(offset_m * np.cos(az_rad))

        # Phase 2: shift z by terrain elevation so depth is from surface
        if use_topo and apply_topo and elev_at_x is not None:
            # Z_mesh[iz, ix] = elev_at_x[ix] - z_arr[iz]
            X_mesh, Z_depth = np.meshgrid(x_arr, np.abs(z_arr))
            elev_broadcast  = np.tile(elev_at_x, (len(z_arr), 1))
            Z_mesh = elev_broadcast - Z_depth
            z_axis_title = "Elevation (m a.s.l.)"
        else:
            X_mesh, Z_mesh = np.meshgrid(x_arr, -np.abs(z_arr))
            z_axis_title = "Depth (m)"

        Y = np.full_like(X_mesh, y_center)

        surf_color = np.log10(np.clip(rho, 1e-3, None))
        c_lo = np.log10(max(vmin, 1e-3))
        c_hi = np.log10(max(vmax, vmin + 1e-3))

        # Plotly Surface does NOT render NaN surfacecolor as transparent —
        # it renders them at cmin colour (dark blue).  Transfer the ρ mask
        # to Z_mesh instead: NaN in z creates genuine transparent holes.
        outside = ~np.isfinite(surf_color)
        if outside.any():
            Z_mesh = Z_mesh.astype(float).copy()
            Z_mesh[outside] = np.nan

        cz = go.surface.Contours(
            z=go.surface.contours.Z(
                show=contours, usecolormap=True,
                highlightcolor=colors["text"], project_z=False,
                start=c_lo, end=c_hi,
                size=(c_hi - c_lo) / max(n_contours, 1),
            )
        )

        traces.append(go.Surface(
            x=X_mesh, y=Y, z=Z_mesh,
            surfacecolor=surf_color,
            customdata=surf_color,
            cmin=c_lo, cmax=c_hi,
            colorscale=cmap,
            opacity=opacity,
            contours=cz,
            showscale=(i == 0),
            colorbar=dict(
                title=dict(text="log₁₀ρ (Ω·m)", side="right",
                           font=dict(color=colors["text"], size=10)),
                tickfont=dict(color=colors["tick"], size=9),
            ) if i == 0 else None,
            name=name,
            hovertemplate=(f"<b>{name}</b><br>"
                           "x: %{x:.0f} m<br>z: %{z:.0f} m<br>"
                           "log₁₀ρ: %{customdata:.2f}<extra></extra>"),
        ))

        # Phase 1: terrain ribbon
        if use_topo and elev_at_x is not None:
            traces.append(_build_topo_strip_trace(
                sta_x=x_arr, sta_elev=elev_at_x,
                y_center=y_center, panel_thick=panel_thick,
                opacity=topo_opacity, dark=dark,
            ))
            if show_sta_markers:
                sta_elev_at_sta = _elev_along_line(sta_x, sta_x, sta_elev)
                traces.append(_build_station_markers_trace(
                    sta_x=sta_x, y_center=y_center,
                    sta_elev=sta_elev_at_sta,
                    sta_names=sta_names_line, dark=dark,
                    symbol=sta_symbol, size=sta_size, color=sta_color,
                ))

        if show_labels:
            label_z = (float(elev_at_x[len(elev_at_x) // 2])
                       if use_topo and elev_at_x is not None else 0.0)
            annotations.append(dict(
                x=float(x_arr[len(x_arr) // 2]),
                y=y_center,
                z=label_z,
                text=f"<b>{name}</b>",
                showarrow=False,
                font=dict(color=colors["text"], size=11),
            ))

    layout = _dark_scene(title, dark=dark)
    layout["scene"]["annotations"] = annotations
    if use_topo and apply_topo:
        layout["scene"]["zaxis"]["title"] = dict(
            text=z_axis_title, font=dict(color=colors["text"], size=10))
    return go.Figure(data=traces, layout=go.Layout(**layout))


def _build_block_fig(
    x_arr, y_arr, z_arr, rho_3d: np.ndarray,
    cmap: str, vmin: float, vmax: float,
    opacity: float,
    isovalue_lo: float, isovalue_hi: float,
    n_surfaces: int,
    clip_x: float, clip_y: float, clip_z: float,
    contours: bool, title: str,
    rho_lo: float = 1, rho_hi: float = 1e8,
    depth_lo: float = 0, depth_hi: float = 1e6,
    dark: bool = True,
    topo_src: str = "none",
    topo_opacity: float = 0.70,
    profiles: dict | None = None,
) -> go.Figure:
    colors = _scene_colors(dark)
    rho_vol = rho_3d.transpose(1, 0, 2)   # (n_x, n_lines, n_z)

    # Universal depth clipping on z_arr
    z = np.abs(z_arr)
    z_mask = (z >= depth_lo) & (z <= depth_hi)
    if z_mask.any():
        z_arr   = z_arr[z_mask]
        rho_vol = rho_vol[:, :, z_mask]

    # Do NOT call _apply_rho_mask here.  go.Volume natively makes values
    # outside [isomin, isomax] transparent, so the full dense scatter cloud
    # must be preserved.  Removing cells (NaN) collapses the cloud and leaves
    # go.Volume unable to reconstruct a volume — the block goes empty.

    # Upsample Y for smoother iso-surfaces
    n_lines = len(y_arr)
    n_y_dense = max(20, n_lines * 6)
    if n_lines >= 2:
        y_dense = np.linspace(y_arr[0], y_arr[-1], n_y_dense)
        rho_dense = np.zeros((rho_vol.shape[0], n_y_dense, rho_vol.shape[2]))
        for ix in range(rho_vol.shape[0]):
            for iz in range(rho_vol.shape[2]):
                rho_dense[ix, :, iz] = _interp_finite_1d(
                    y_dense, y_arr, rho_vol[ix, :, iz])
        rho_vol, y_arr = rho_dense, y_dense

    # Derive isovalue range from ρ filter if user wants specificity
    c_lo_data = np.log10(max(vmin, 1e-3))
    c_hi_data = np.log10(max(vmax, vmin + 1e-3))
    iso_lo = np.log10(max(rho_lo, 1e-3))
    iso_hi = np.log10(max(rho_hi, rho_lo + 1e-3))
    # Honour explicit isovalue overrides if they differ from defaults
    if isovalue_lo != 0.5 or isovalue_hi != 2.5:
        iso_lo = float(isovalue_lo)
        iso_hi = float(isovalue_hi)

    rho_log = np.where(
        np.isfinite(rho_vol) & (rho_vol > 0),
        np.log10(np.clip(rho_vol, 1e-3, None)),
        np.nan,
    )

    finite = rho_log[np.isfinite(rho_log)]
    if not finite.size:
        fig = _empty_3d_fig()
        fig.update_layout(**_dark_scene(title, dark=dark))
        fig.add_annotation(
            text="No finite resistivity values remain after the current filters.",
            x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False,
            font=dict(color=colors["warn"], size=12),
        )
        return fig

    data_lo = float(finite.min())
    data_hi = float(finite.max())
    iso_lo = max(float(iso_lo), data_lo)
    iso_hi = min(float(iso_hi), data_hi)
    if iso_hi <= iso_lo:
        if data_hi > data_lo:
            iso_lo, iso_hi = data_lo, data_hi
        else:
            delta = max(abs(data_lo) * 0.05, 0.05)
            iso_lo, iso_hi = data_lo - delta, data_hi + delta

    X, Y, Z_ = np.meshgrid(x_arr, y_arr, -np.abs(z_arr), indexing="ij")

    x_keep = _axis_fraction_mask(x_arr, clip_x)
    y_keep = _axis_fraction_mask(y_arr, clip_y)
    z_range = np.abs(z_arr).max() or 1.0
    mask = (
        x_keep[:, None, None] &
        y_keep[None, :, None] &
        (np.abs(Z_) <= z_range * max(0.0, min(1.0, float(clip_z)))) &
        np.isfinite(rho_log)
    )
    x_flat   = X.flatten()[mask.flatten()]
    y_flat   = Y.flatten()[mask.flatten()]
    z_flat   = Z_.flatten()[mask.flatten()]
    val_flat = rho_log.flatten()[mask.flatten()]

    if val_flat.size < 4:
        fig = _empty_3d_fig()
        fig.update_layout(**_dark_scene(title, dark=dark))
        fig.add_annotation(
            text="Not enough visible cells to build a 3D block. Relax the filters.",
            x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False,
            font=dict(color=colors["warn"], size=12),
        )
        return fig

    trace = go.Volume(
        x=x_flat, y=y_flat, z=z_flat,
        value=val_flat,
        isomin=iso_lo,
        isomax=iso_hi,
        surface_count=max(2, int(n_surfaces)),
        opacity=min(max(float(opacity), 0.05), 1.0),
        opacityscale=[
            [0.0, 0.0],
            [0.2, 0.3],
            [0.5, 0.7],
            [1.0, 1.0],
        ],
        colorscale=cmap,
        cmin=c_lo_data, cmax=c_hi_data,
        showscale=True,
        colorbar=dict(
            title=dict(text="log₁₀ρ (Ω·m)",
                       font=dict(color=colors["text"], size=10)),
            tickfont=dict(color=colors["tick"], size=9),
        ),
        caps=dict(
            x=dict(show=clip_x < 1, fill=0.8),
            y=dict(show=clip_y < 1, fill=0.8),
            z=dict(show=clip_z < 1, fill=0.8),
        ),
        hovertemplate="x: %{x:.0f} m<br>y: %{y:.0f} m<br>z: %{z:.0f} m<br>"
                      "log₁₀ρ: %{value:.2f}<extra></extra>",
    )
    data_traces = [trace]
    # Phase 1: terrain surface above block model
    if topo_src != "none" and profiles:
        topo_surf = _build_topo_2d_surface(
            profiles=profiles,
            x_arr=x_arr, y_arr=y_arr,
            topo_src=topo_src,
            opacity=topo_opacity, dark=dark,
        )
        if topo_surf is not None:
            data_traces.append(topo_surf)
    return go.Figure(data=data_traces, layout=go.Layout(**_dark_scene(title, dark=dark)))


def _build_depth_slices_fig(
    x_arr, y_arr,
    depth_lo: float, depth_hi: float, n_slices: int,
    rho_fn,
    cmap: str, vmin: float, vmax: float,
    opacity: float, contours: bool, n_contours: int,
    show_labels: bool, title: str,
    rho_lo: float = 1, rho_hi: float = 1e8,
    dark: bool = True,
    topo_src: str = "none",
    topo_opacity: float = 0.70,
    show_sta_markers: bool = True,
    apply_topo: bool = False,
    profiles: dict | None = None,
) -> go.Figure:
    colors  = _scene_colors(dark)
    depths  = np.linspace(depth_lo, depth_hi, n_slices)
    c_lo    = np.log10(max(vmin, 1e-3))
    c_hi    = np.log10(max(vmax, vmin + 1e-3))
    X, Y    = np.meshgrid(x_arr, y_arr)
    traces  = []
    use_topo = topo_src != "none"

    # Build 2-D elevation grid for depth-shift (Phase 2 + 3: source-aware)
    elev_2d: np.ndarray | None = None
    if use_topo and apply_topo and profiles:
        try:
            from scipy.interpolate import griddata
            pts_x, pts_y, pts_e = [], [], []
            for j, (name, p) in enumerate(profiles.items()):
                sta_x    = np.asarray(p.get("sta_x", []), float)
                sta_elev = _resolve_profile_elevs(p, topo_src)
                y_center = y_arr[j] if j < len(y_arr) else float(y_arr[-1])
                if sta_x.size and np.any(sta_elev != 0.0):
                    pts_x.extend(sta_x.tolist())
                    pts_y.extend([y_center] * len(sta_x))
                    pts_e.extend(sta_elev.tolist())
            if pts_x:
                pts  = np.column_stack([pts_x, pts_y])
                vals = np.asarray(pts_e, float)
                elev_2d = griddata(pts, vals, (X.T, Y.T), method="linear",
                                   fill_value=float(np.nanmean(vals)))
                elev_2d = elev_2d.T   # match (n_y, n_x) of X, Y
        except Exception:
            elev_2d = None

    for i, depth in enumerate(depths):
        rho_slice = rho_fn(depth)

        # Universal ρ mask
        rho_slice = _apply_rho_mask(rho_slice, rho_lo, rho_hi)

        surf_color = np.log10(np.clip(rho_slice, 1e-3, None))

        # Phase 2: drape slice onto terrain
        if use_topo and apply_topo and elev_2d is not None:
            Z = elev_2d - float(depth)
            z_label = f"z = elev − {depth:.0f} m"
        else:
            Z = np.full_like(X, -float(depth), dtype=float)
            z_label = f"z = {depth:.0f} m"

        # Plotly Surface does NOT render NaN surfacecolor as transparent —
        # it renders them at cmin colour (dark blue).  Transfer the ρ mask
        # to Z instead: NaN in z creates genuine transparent holes.
        outside = ~np.isfinite(surf_color)
        if outside.any():
            Z = Z.astype(float).copy()
            Z[outside] = np.nan

        cz = go.surface.Contours(
            z=go.surface.contours.Z(
                show=contours, usecolormap=True,
                start=c_lo, end=c_hi,
                size=(c_hi - c_lo) / max(n_contours, 1),
            )
        )
        traces.append(go.Surface(
            x=X, y=Y, z=Z,
            surfacecolor=surf_color,
            customdata=surf_color,
            cmin=c_lo, cmax=c_hi,
            colorscale=cmap,
            opacity=opacity,
            contours=cz,
            showscale=(i == n_slices - 1),
            colorbar=dict(
                title=dict(text="log₁₀ρ (Ω·m)", side="right",
                           font=dict(color=colors["text"], size=10)),
                tickfont=dict(color=colors["tick"], size=9),
            ) if i == n_slices - 1 else None,
            name=z_label,
            hovertemplate=f"<b>Depth {depth:.0f} m</b><br>"
                          "x: %{x:.0f}<br>y: %{y:.0f}<br>"
                          "log₁₀ρ: %{customdata:.2f}<extra></extra>",
        ))

    # Phase 1: add terrain surface on top
    if use_topo and profiles:
        topo_surf = _build_topo_2d_surface(
            profiles=profiles,
            x_arr=x_arr, y_arr=y_arr,
            topo_src=topo_src,
            opacity=topo_opacity, dark=dark,
        )
        if topo_surf is not None:
            traces.append(topo_surf)

    layout = _dark_scene(title, dark=dark)
    z_axis_lbl = ("Elevation (m a.s.l.)" if use_topo and apply_topo
                  else "Depth (m)")
    layout["scene"]["zaxis"]["title"] = dict(
        text=z_axis_lbl, font=dict(color=colors["text"], size=10))
    return go.Figure(data=traces, layout=go.Layout(**layout))


# ── Data extraction helpers ───────────────────────────────────────────────────

def _rho_xy_from_site(site) -> tuple:
    try:
        freq_raw = site.freq
        if freq_raw is None:
            return None, None
        freq   = np.asarray(freq_raw, float).ravel()
        if freq.size == 0:
            return None, None

        rho_raw = site.rho
        rho_xy  = None

        if rho_raw is not None:
            rho = np.asarray(rho_raw, float)
            if rho.ndim == 3 and rho.shape[-2:] == (2, 2):
                rho_xy = rho[:, 0, 1]
            elif rho.ndim == 2 and rho.shape[1] >= 2:
                rho_xy = rho[:, 1]
            elif rho.ndim == 1:
                rho_xy = rho
            else:
                rho_xy = rho.ravel()
            rho_xy = np.abs(rho_xy.ravel())
            if not np.any(np.isfinite(rho_xy) & (rho_xy > 0)):
                rho_xy = None

        if rho_xy is None:
            z_raw = site.z
            if z_raw is not None:
                z_arr = np.asarray(z_raw, complex)
                if z_arr.ndim == 3 and z_arr.shape[-2:] == (2, 2):
                    z_xy = z_arr[:, 0, 1]
                elif z_arr.ndim == 2 and z_arr.shape[1] >= 2:
                    z_xy = z_arr[:, 1]
                else:
                    z_xy = z_arr.ravel()
                n = min(len(freq), len(z_xy))
                rho_xy = (0.2 * np.abs(z_xy[:n]) ** 2
                          / np.clip(freq[:n], 1e-10, None))

        if rho_xy is None:
            return None, None
        rho_xy = rho_xy.ravel()
        if not np.any(np.isfinite(rho_xy)):
            return None, None
        n = min(len(freq), len(rho_xy))
        return freq[:n], rho_xy[:n]
    except Exception:
        return None, None


def _station_distance_m(coords_list: list) -> np.ndarray:
    n     = len(coords_list)
    dists = np.zeros(n, float)
    R     = 6_371_000.0
    for i in range(1, n):
        lat1 = np.radians(coords_list[i-1][0]); lon1 = np.radians(coords_list[i-1][1])
        lat2 = np.radians(coords_list[i][0]);   lon2 = np.radians(coords_list[i][1])
        dlat, dlon = lat2 - lat1, lon2 - lon1
        a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
        dists[i] = dists[i-1] + 2 * R * np.arcsin(min(np.sqrt(a), 1.0))
    if dists[-1] < 1.0:
        dists = np.arange(n, dtype=float) * 100.0
    return dists


# ── Topography helpers ────────────────────────────────────────────────────────

def _safe_elev(coords) -> float | None:
    """Extract scalar elevation (m) from a (lat, lon, elev) tuple or similar."""
    try:
        if coords is None:
            return None
        elev = float(coords[2])
        return elev if np.isfinite(elev) else None
    except Exception:
        return None


def _elev_along_line(x_arr: np.ndarray, sta_x: np.ndarray,
                     sta_elev: np.ndarray) -> np.ndarray:
    """Interpolate station elevations onto a uniform x_arr (1-D)."""
    if sta_x is None or sta_elev is None or len(sta_x) == 0:
        return np.zeros(len(x_arr), float)
    order = np.argsort(sta_x)
    sx, se = sta_x[order], sta_elev[order]
    finite = np.isfinite(se)
    if not finite.any():
        return np.zeros(len(x_arr), float)
    sx, se = sx[finite], se[finite]
    if len(sx) == 1:
        return np.full(len(x_arr), se[0], float)
    return np.interp(x_arr, sx, se)


def _build_topo_strip_trace(
    sta_x: np.ndarray, sta_elev: np.ndarray,
    y_center: float, panel_thick: float,
    opacity: float, dark: bool,
) -> go.Surface:
    """Return a terrain ribbon go.Surface for one fence line."""
    colors = _scene_colors(dark)
    half = max(float(panel_thick) / 2.0, 25.0)
    n_x  = len(sta_x)
    X    = np.tile(sta_x, (2, 1))
    Z    = np.tile(sta_elev, (2, 1))
    Y    = np.array([[y_center - half] * n_x,
                     [y_center + half] * n_x])
    topo_color = np.zeros_like(Z)
    return go.Surface(
        x=X, y=Y, z=Z,
        surfacecolor=topo_color,
        colorscale=[[0, "#7dc4e4"], [1, "#7dc4e4"]],   # sky-blue ribbon
        opacity=min(max(float(opacity), 0.05), 1.0),
        showscale=False,
        name="Terrain",
        hovertemplate="x: %{x:.0f} m<br>Elev: %{z:.1f} m<extra>Terrain</extra>",
    )


def _build_station_markers_trace(
    sta_x: np.ndarray, y_center: float,
    sta_elev: np.ndarray, sta_names: list[str],
    dark: bool,
    symbol: str = "diamond-open",
    size: int = 8,
    color: str = "black",
) -> go.Scatter3d:
    """Return station marker go.Scatter3d at terrain surface.

    API convention (api/station.py) mapped to valid Plotly 3D symbols
    (Plotly 3D supports: circle, circle-open, cross, diamond, diamond-open,
    square, square-open, x — triangles are 2D only):
      skin-depth / pseudo  → diamond-open  (hollow ◇, mirrors facecolor=white)
      inversion            → diamond       (filled ◆, mirrors facecolor=black)
    Callers pass the already-resolved symbol; this function just renders it.
    """
    colors = _scene_colors(dark)
    y_arr  = np.full(len(sta_x), y_center)
    text   = [f"<b>{n}</b><br>x: {x:.0f} m<br>z: {e:.1f} m a.s.l."
              for n, x, e in zip(sta_names, sta_x, sta_elev)]
    # Open (-open) variants are inherently hollow in Plotly 3D; line controls edge width
    return go.Scatter3d(
        x=sta_x, y=y_arr, z=sta_elev,
        mode="markers+text",
        marker=dict(
            symbol=symbol,
            size=int(size),
            color=color,
            line=dict(color=color, width=1.5),
        ),
        text=sta_names,
        textposition="top center",
        textfont=dict(size=8, color=colors["text"]),
        hovertext=text, hoverinfo="text",
        name="Stations",
    )


def _build_topo_2d_surface(
    profiles: dict,
    x_arr: np.ndarray, y_arr: np.ndarray,
    topo_src: str,
    opacity: float, dark: bool,
) -> go.Surface | None:
    """Build a full 2-D terrain go.Surface for block / depth-slice modes."""
    line_names = list(profiles.keys())
    elev_rows  = []
    for name in line_names:
        p        = profiles[name]
        sta_x    = np.asarray(p.get("sta_x", []), float)
        sta_elev = _resolve_profile_elevs(p, topo_src)   # Phase 3: source-aware
        if sta_x.size == 0 or not np.any(sta_elev != 0.0):
            return None
        elev_rows.append(_elev_along_line(x_arr, sta_x, sta_elev))

    if not elev_rows:
        return None

    elev_mat = np.vstack(elev_rows)   # (n_lines, n_x)
    if elev_mat.shape[0] == 1:
        elev_mat = np.tile(elev_mat, (2, 1))
        y_plot   = np.array([y_arr[0], y_arr[0]])
    else:
        y_plot = y_arr

    X, Y_ = np.meshgrid(x_arr, y_plot)
    surf_color = np.zeros_like(elev_mat)
    return go.Surface(
        x=X, y=Y_, z=elev_mat,
        surfacecolor=surf_color,
        colorscale=[[0, "#7dc4e4"], [1, "#7dc4e4"]],
        opacity=min(max(float(opacity), 0.05), 1.0),
        showscale=False,
        name="Terrain",
        hovertemplate="x: %{x:.0f} m<br>y: %{y:.0f} m<br>Elev: %{z:.1f} m<extra>Terrain</extra>",
    )


def _parse_topo_upload(contents: str, filename: str) -> list[dict]:
    """Parse a user-uploaded topography file into [{station, elev}, ...] records.

    Supported formats: CSV, H5/HDF5, NPZ.
    Expected columns/arrays: station (or ID/name) + elevation (or elev/Elevation).
    """
    import base64, io
    _, enc = contents.split(",", 1)
    raw = base64.b64decode(enc)
    fname = (filename or "").lower()
    records: list[dict] = []

    try:
        if fname.endswith(".csv"):
            import pandas as pd
            df = pd.read_csv(io.BytesIO(raw), dtype=str)
            df.columns = [c.strip() for c in df.columns]
            # Flexible column detection
            id_col   = next((c for c in df.columns
                             if c.lower() in ("station", "id", "name", "sta")), None)
            elev_col = next((c for c in df.columns
                             if c.lower() in ("elevation", "elev", "z", "alt",
                                              "altitude", "height")), None)
            if id_col is None or elev_col is None:
                return []
            for _, row in df.iterrows():
                name = str(row[id_col]).strip()
                try:
                    e = float(row[elev_col])
                except (TypeError, ValueError):
                    continue
                if name and np.isfinite(e):
                    records.append({"station": name, "elev": e})

        elif fname.endswith((".h5", ".hdf5")):
            import h5py
            with h5py.File(io.BytesIO(raw), "r") as f:
                # Try common dataset names for station IDs
                id_key  = next((k for k in f
                                if k.lower() in ("station", "station_names",
                                                 "id", "name", "sta")), None)
                elv_key = next((k for k in f
                                if k.lower() in ("elevation", "elev", "z",
                                                 "altitude", "height")), None)
                if id_key is None or elv_key is None:
                    return []
                names = [n.decode() if isinstance(n, bytes) else str(n)
                         for n in f[id_key][:]]
                elevs = np.asarray(f[elv_key][:], float)
                for name, e in zip(names, elevs):
                    if name.strip() and np.isfinite(e):
                        records.append({"station": name.strip(), "elev": float(e)})

        elif fname.endswith(".npz"):
            data = np.load(io.BytesIO(raw), allow_pickle=True)
            id_key  = next((k for k in data
                            if k.lower() in ("station", "station_names",
                                             "id", "name", "sta")), None)
            elv_key = next((k for k in data
                            if k.lower() in ("elevation", "elev", "z",
                                             "altitude", "height")), None)
            if id_key is None or elv_key is None:
                return []
            names = [n.decode() if isinstance(n, bytes) else str(n)
                     for n in data[id_key]]
            elevs = np.asarray(data[elv_key], float)
            for name, e in zip(names, elevs):
                if name.strip() and np.isfinite(e):
                    records.append({"station": name.strip(), "elev": float(e)})

    except Exception:
        return []

    return records


def _build_elevs_dict(
    elev_corr_store,
    elev_raw_store,
    station_records,
    topo_upload_store=None,
) -> dict:
    """Merge all elevation sources into {station_name: {elev_sta, elev_raw, elev_corr}}."""
    out: dict[str, dict] = {}

    def _try_float(v):
        try:
            f = float(v)
            return f if np.isfinite(f) else None
        except (TypeError, ValueError):
            return None

    for r in (station_records or []):
        name = str(r.get("ID") or r.get("station") or r.get("name", "")).strip()
        if not name:
            continue
        e = _try_float(r.get("Elevation") or r.get("elevation") or r.get("elev"))
        if e is not None:
            out.setdefault(name, {})["elev_sta"] = e

    for r in (elev_raw_store or []):
        name = str(r.get("station") or r.get("ID") or r.get("name", "")).strip()
        if not name:
            continue
        e = _try_float(r.get("elev"))
        if e is not None:
            out.setdefault(name, {})["elev_raw"] = e

    for r in (elev_corr_store or []):
        name = str(r.get("station") or r.get("ID") or r.get("name", "")).strip()
        if not name:
            continue
        e = _try_float(r.get("elev_corrected") or r.get("elev"))
        if e is not None:
            out.setdefault(name, {})["elev_corr"] = e

    for r in (topo_upload_store or []):
        name = str(r.get("station") or r.get("ID") or r.get("name", "")).strip()
        if not name:
            continue
        e = _try_float(r.get("elev") or r.get("elevation") or r.get("Elevation"))
        if e is not None:
            out.setdefault(name, {})["elev_upload"] = e

    return out


def _resolve_profile_elevs(profile: dict, topo_src: str) -> np.ndarray:
    """Return the elevation array for a profile based on the chosen topo source.

    Falls back gracefully: elev-corr → elev-raw → sta_elev → zeros.
    """
    sta_x = np.asarray(profile.get("sta_x", []), float)
    n = len(sta_x)
    if n == 0:
        return np.array([], float)

    if topo_src == "elev-corr":
        preferred = ("sta_elev_corr", "sta_elev_raw", "sta_elev")
    elif topo_src == "elev-raw":
        preferred = ("sta_elev_raw", "sta_elev")
    elif topo_src == "upload":
        preferred = ("sta_elev_upload", "sta_elev")
    else:  # "stations" or anything else
        preferred = ("sta_elev",)

    for key in preferred:
        arr = profile.get(key)
        if arr is not None:
            a = np.asarray(arr, float)
            if a.size == n and np.any(np.isfinite(a) & (a != 0.0)):
                # Replace NaN with linear interpolation / mean
                bad = ~np.isfinite(a)
                if bad.any():
                    good_idx = np.where(~bad)[0]
                    if good_idx.size:
                        a[bad] = np.interp(np.where(bad)[0], good_idx, a[good_idx])
                    else:
                        a = np.zeros(n, float)
                return a
    return np.zeros(n, float)


def _has_survey_line_profiles(store_data: dict | None) -> bool:
    """True when loaded survey metadata contains at least two named lines."""
    if not store_data:
        return False
    line_counts = store_data.get("line_counts") or {}
    named = [str(k).strip() for k, v in line_counts.items() if str(k).strip() and v]
    return len(named) >= 2


def _has_cached_inversion_result(session_id: str | None) -> bool:
    """True when a session inversion result exists and can become a 3D source."""
    result = cache_get_inversion_result(session_id or "")
    if result is None:
        return False
    try:
        _profiles_from_inversion_result(result)
        return True
    except Exception:
        return False


def _profiles_from_pseudo(
    sites,
    store_data: dict | None = None,
    *,
    require_line_metadata: bool = False,
) -> dict:
    try:
        site_list = list(sites)
    except Exception:
        return {}
    if not site_list:
        return {}

    sta_to_line: dict[str, str] = {}
    if store_data:
        for r in store_data.get("station_records", []):
            sid  = str(r.get("ID", "") or r.get("Station", "")).strip()
            line = str(r.get("Line", "") or "").strip() or "default"
            if sid:
                sta_to_line[sid] = line

    if require_line_metadata and not _has_survey_line_profiles(store_data):
        return {}

    raw: dict[str, list] = {}

    for site in site_list:
        try:
            freq, rho_xy = _rho_xy_from_site(site)
            if freq is None or rho_xy is None:
                continue
            depths = 503.0 * np.sqrt(
                np.clip(rho_xy, 1e-3, None) / np.clip(freq, 1e-10, None)
            )
            valid = np.isfinite(depths) & np.isfinite(rho_xy) & (depths > 0)
            if not valid.any():
                continue
            depths, rho_xy = depths[valid], rho_xy[valid]

            station = str(getattr(site, "name", None) or
                          getattr(site, "station", None) or f"S{len(raw)}")
            if require_line_metadata and station not in sta_to_line:
                continue
            line = (sta_to_line.get(station)
                    or str(getattr(site, "Line", None)
                           or getattr(site, "line", None) or "")
                    or "default")
            if require_line_metadata and line in {"", "default", "uploaded", "—"}:
                continue
            try:
                coords = site.coords
            except Exception:
                coords = (0.0, 0.0, 0.0)
            raw.setdefault(line, []).append((station, depths, rho_xy, coords))
        except Exception:
            continue

    if not raw:
        return {}

    result = {}
    for line, items in raw.items():
        if len(items) < 1:
            continue
        coords_list = [it[3] for it in items]
        x_arr       = _station_distance_m(coords_list)
        all_depths  = np.concatenate([it[1] for it in items])
        all_depths  = all_depths[np.isfinite(all_depths) & (all_depths > 0)]
        if all_depths.size == 0:
            continue
        z_arr    = np.linspace(all_depths.min(), np.percentile(all_depths, 95), 40)
        rho_grid = np.zeros((len(z_arr), len(items)), float)
        for j, (_, dep_j, rho_j, _) in enumerate(items):
            sort_idx          = np.argsort(dep_j)
            dep_s, rho_s      = dep_j[sort_idx], rho_j[sort_idx]
            rho_grid[:, j]    = np.interp(z_arr, dep_s, rho_s,
                                           left=rho_s[0], right=rho_s[-1])

        # Topo metadata
        sta_names = [it[0] for it in items]
        sta_lat   = [float(it[3][0]) if it[3] is not None else 0.0 for it in items]
        sta_lon   = [float(it[3][1]) if it[3] is not None else 0.0 for it in items]
        sta_elev  = [_safe_elev(it[3]) for it in items]
        # Fall back to 0.0 for missing elevations so arrays are always same length
        sta_elev  = [e if e is not None else 0.0 for e in sta_elev]

        result[line] = {
            "x": x_arr, "z": z_arr, "rho": rho_grid,
            "sta_x":     x_arr.tolist(),
            "sta_elev":  sta_elev,
            "sta_names": sta_names,
            "sta_lat":   sta_lat,
            "sta_lon":   sta_lon,
        }
    return result


def _rho_log_to_ohm_m(rho: np.ndarray) -> np.ndarray:
    """Convert pyCSAMT log10-rho sections to linear resistivity when needed."""
    arr = np.asarray(rho, float)
    finite = arr[np.isfinite(arr)]
    if finite.size and float(np.nanmax(finite)) <= 12.0:
        return 10.0 ** arr
    return arr


def _profiles_from_inversion_result(
    result,
    elevs_by_name: dict | None = None,
) -> dict:
    """Convert a cached InversionResult into map3D profile dictionaries."""
    model = result.to_resistivity_model()
    rho = _rho_log_to_ohm_m(np.asarray(model.rho_2d, float))
    if rho.ndim != 2 or rho.size == 0:
        return {}
    x = np.asarray(model.x_centers, float)
    z = np.asarray(model.z_centers, float)
    if x.size != rho.shape[1]:
        x = np.arange(rho.shape[1], dtype=float)
    if z.size != rho.shape[0]:
        z = np.arange(rho.shape[0], dtype=float)
    label = f"{result.method.upper()} {result.dimension.upper()} inversion"

    sta_x     = np.asarray(getattr(model, "station_x", x), float)
    sta_names = list(getattr(model, "station_names", None) or
                     [f"S{i:03d}" for i in range(len(sta_x))])

    profile: dict = {
        "x": x, "z": z, "rho": rho,
        "sta_x":     sta_x.tolist(),
        "sta_names": sta_names,
    }

    if elevs_by_name:
        def _pick(name, key):
            return float(elevs_by_name.get(name, {}).get(key, 0.0))

        profile["sta_elev"]        = [_pick(n, "elev_sta")    for n in sta_names]
        profile["sta_elev_raw"]    = [_pick(n, "elev_raw")    for n in sta_names]
        profile["sta_elev_corr"]   = [_pick(n, "elev_corr")   for n in sta_names]
        profile["sta_elev_upload"] = [_pick(n, "elev_upload")  for n in sta_names]
    else:
        profile["sta_elev"]        = [0.0] * len(sta_names)
        profile["sta_elev_raw"]    = [0.0] * len(sta_names)
        profile["sta_elev_corr"]   = [0.0] * len(sta_names)
        profile["sta_elev_upload"] = [0.0] * len(sta_names)

    return {label: profile}


def _profiles_from_inversion(session_id: str) -> dict:
    result = cache_get_inversion_result(session_id)
    if result is not None:
        return _profiles_from_inversion_result(result)

    # Fallback kept for legacy desktop-controller sessions.
    try:
        from pycsamt.app.desktop.controllers.inversion_controller import InversionController
        ctrl  = InversionController()
        model = getattr(getattr(ctrl, "state", None), "model", None)
        if model is None:
            return {}
        x   = np.asarray(getattr(model, "x_centers", np.arange(model.rho_2d.shape[1])))
        z   = np.asarray(getattr(model, "z_centers", np.arange(model.rho_2d.shape[0])))
        rho = np.asarray(model.rho_2d)
        return {"Line 1": {"x": x, "z": z, "rho": _rho_log_to_ohm_m(rho)}}
    except Exception:
        return {}


# ── Register callbacks ────────────────────────────────────────────────────────

def register_map3d(app) -> None:
    _register_mode_switch(app)
    _register_source_options(app)
    _register_presets(app)
    _register_generate(app)
    _register_display(app)
    _register_export_html(app)


def _register_mode_switch(app) -> None:
    """Mode pill buttons → store; clientside syncs classes + panel visibility."""

    @app.callback(
        Output(IDs.MAP3D_ACTIVE_MODE, "data"),
        [Input(f"map3d-mode-btn-{slug}", "n_clicks") for slug in _MODE_SLUGS],
        prevent_initial_call=True,
    )
    def switch_mode(*_clicks):
        if not any(_clicks):
            raise PreventUpdate
        triggered = ctx.triggered_id or ""
        for slug in _MODE_SLUGS:
            if triggered == f"map3d-mode-btn-{slug}":
                return slug
        raise PreventUpdate

    _slugs_js = '["' + '","'.join(_MODE_SLUGS) + '"]'
    _panel_ids = ["map3d-fence-panel", "map3d-block-panel", "map3d-depth-panel"]

    clientside_callback(
        f"""
        function(active) {{
            var slugs = {_slugs_js};
            var show  = {{display:"block"}};
            var none  = {{display:"none"}};
            var panelStyles = slugs.map(function(s) {{
                return s === active ? show : none;
            }});
            var btnClasses = slugs.map(function(s) {{
                return "fwd-dim-btn" + (s === active ? " active" : "");
            }});
            // iso-band card: only visible in block mode
            var isoStyle = active === "block" ? show : none;
            return panelStyles.concat(btnClasses).concat([isoStyle]);
        }}
        """,
        [Output(pid, "style")                            for pid in _panel_ids] +
        [Output(f"map3d-mode-btn-{slug}", "className")   for slug in _MODE_SLUGS] +
        [Output("map3d-iso-band-card", "style")],
        Input(IDs.MAP3D_ACTIVE_MODE, "data"),
    )

    # Context bar update when mode changes
    @app.callback(
        Output(IDs.MAP3D_INFO, "children", allow_duplicate=True),
        Input(IDs.MAP3D_ACTIVE_MODE, "data"),
        prevent_initial_call=True,
    )
    def update_ctx_on_mode(slug):
        if not slug:
            raise PreventUpdate


def _source_option(label: str, value: str, icon: str, color: str,
                   enabled: bool, reason: str):
    suffix = "" if enabled else f" - {reason}"
    return {
        "label": html.Span([
            html.I(className=f"bi {icon} me-2", style={"color": color}),
            html.Span(label),
            html.Span(suffix, className="map3d-source-reason ms-1"),
        ], className="d-flex align-items-center"),
        "value": value,
        "disabled": not enabled,
    }


def _register_source_options(app) -> None:
    """Enable only 3D data sources that have real session data behind them."""

    @app.callback(
        Output(IDs.MAP3D_DATA_SRC, "options"),
        Output(IDs.MAP3D_DATA_SRC, "value"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.SESSION_ID, "data"),
        State(IDs.MAP3D_DATA_SRC, "value"),
    )
    def update_source_options(store_data, session_id, current_value):
        has_sites = bool((store_data or {}).get("n_stations"))
        has_lines = _has_survey_line_profiles(store_data)
        has_inv = _has_cached_inversion_result(session_id)

        sources = [
            ("pseudo", "Skin-depth pseudo", "bi-database", "#89b4fa",
             has_sites, "load survey first"),
            ("inversion", "Session inversion result", "bi-layers", "#a6e3a1",
             has_inv, "run/save inversion first"),
            ("profiles", "Survey line profiles", "bi-diagram-3", "#fab387",
             has_lines, "need named lines"),
        ]
        options = [
            _source_option(label, value, icon, color, enabled, reason)
            for value, label, icon, color, enabled, reason in sources
        ]
        enabled_values = [value for value, *_rest, enabled, _reason in sources
                          if enabled]
        value = current_value if current_value in enabled_values else (
            enabled_values[0] if enabled_values else "pseudo"
        )
        return options, value
        icon  = _MODE_ICON.get(slug, "bi-box")
        label = _MODE_LABEL.get(slug, slug.title())
        return html.Div([
            html.I(className=f"bi {icon} me-2",
                   style={"fontSize": "12px", "color": "var(--blue)"}),
            html.Span(label,
                      style={"fontSize": "11px", "fontWeight": "600",
                             "color": "var(--txt)", "marginRight": "10px"}),
            html.Span("Load & Generate 3D to build the scene.",
                      style={"fontSize": "10px", "color": "var(--sub0)"}),
        ], className="d-flex align-items-center")


def _register_presets(app) -> None:
    """Quick-select preset buttons for ρ range and depth window."""

    @app.callback(
        Output(IDs.MAP3D_RHO_LO, "value"),
        Output(IDs.MAP3D_RHO_HI, "value"),
        [Input(pid, "n_clicks") for pid in _RHO_PRESETS],
        prevent_initial_call=True,
    )
    def set_rho_preset(*_):
        tid = ctx.triggered_id
        lo, hi = _RHO_PRESETS.get(tid, (no_update, no_update))
        return lo, hi

    @app.callback(
        Output(IDs.MAP3D_DEPTH_LO, "value"),
        Output(IDs.MAP3D_DEPTH_HI, "value"),
        [Input(pid, "n_clicks") for pid in _DEPTH_PRESETS],
        prevent_initial_call=True,
    )
    def set_depth_preset(*_):
        tid = ctx.triggered_id
        lo, hi = _DEPTH_PRESETS.get(tid, (no_update, no_update))
        return lo, hi

    # Auto-sync iso band lo/hi from ρ filter (convenience)
    @app.callback(
        Output(IDs.MAP3D_ISOVALUE_LO, "value"),
        Output(IDs.MAP3D_ISOVALUE_HI, "value"),
        Input(IDs.MAP3D_RHO_LO, "value"),
        Input(IDs.MAP3D_RHO_HI, "value"),
        prevent_initial_call=True,
    )
    def sync_iso_from_rho(rho_lo, rho_hi):
        if rho_lo is None or rho_hi is None:
            raise PreventUpdate
        rlo = float(rho_lo) if rho_lo else 1
        rhi = float(rho_hi) if rho_hi else 100_000
        return round(np.log10(max(rlo, 1e-3)), 2), round(np.log10(max(rhi, rlo + 1e-3)), 2)


def _register_generate(app) -> None:
    """
    'Load & Generate' button: extract Sites data → serialise to MAP3D_GRID_STORE.
    Only this callback touches the cache / Sites object.
    """
    @app.callback(
        Output(IDs.MAP3D_GRID_STORE, "data"),
        Output(IDs.MAP3D_RANGE_HINT, "children"),
        Output(IDs.MAP3D_SPINNER,    "children"),
        Output(IDs.TOAST_ERROR,      "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,       "children", allow_duplicate=True),
        Input(IDs.BTN_MAP3D_RUN,     "n_clicks"),
        State(IDs.MAP3D_DATA_SRC,    "value"),
        State(IDs.SESSION_ID,        "data"),
        State(IDs.STORE_DATA,        "data"),
        State("tool-elev-corrected-store",      "data"),
        State("tool-elev-raw-store",            "data"),
        State(IDs.MAP3D_TOPO_UPLOAD_STORE,      "data"),
        prevent_initial_call=True,
    )
    def generate_grid(n_clicks, data_src, session_id, store_data,
                      elev_corr_store, elev_raw_store, topo_upload_store):
        if not n_clicks:
            raise PreventUpdate

        data_src = data_src or "pseudo"
        sites    = cache_get(session_id)
        if sites is None:
            return no_update, "", "", True, "Load survey data first."

        # Pre-build merged elevation dict from all available sources (Phase 3)
        station_records = (store_data or {}).get("station_records", [])
        elevs_by_name   = _build_elevs_dict(elev_corr_store, elev_raw_store,
                                             station_records, topo_upload_store)

        try:
            if data_src == "inversion":
                if not _has_cached_inversion_result(session_id):
                    msg = ("Run an inversion first, then use "
                           "Session inversion result.")
                    return (
                        no_update,
                        msg,
                        "",
                        False,
                        "",
                    )
                result   = cache_get_inversion_result(session_id)
                profiles = _profiles_from_inversion_result(
                    result, elevs_by_name=elevs_by_name or None
                )
            elif data_src == "profiles":
                if not _has_survey_line_profiles(store_data):
                    msg = "Survey line profiles require at least two named lines."
                    return (
                        no_update,
                        msg,
                        "",
                        False,
                        "",
                    )
                profiles = _profiles_from_pseudo(
                    sites,
                    store_data=store_data,
                    require_line_metadata=True,
                )
            else:
                profiles = _profiles_from_pseudo(sites, store_data=store_data)
        except Exception as exc:
            return no_update, "", "", True, f"Data extraction failed: {exc}"

        # Phase 3: backfill sta_elev_raw / sta_elev_corr into pseudo profiles
        # using the tool stores so that source selection works at render time.
        if elevs_by_name:
            for p in profiles.values():
                sta_names_p = p.get("sta_names", [])
                if not sta_names_p:
                    continue
                if "sta_elev_raw" not in p:
                    p["sta_elev_raw"]  = [
                        float(elevs_by_name.get(n, {}).get("elev_raw",
                              elevs_by_name.get(n, {}).get("elev_sta", 0.0)))
                        for n in sta_names_p
                    ]
                if "sta_elev_corr" not in p:
                    p["sta_elev_corr"] = [
                        float(elevs_by_name.get(n, {}).get("elev_corr",
                              elevs_by_name.get(n, {}).get("elev_raw",
                              elevs_by_name.get(n, {}).get("elev_sta", 0.0))))
                        for n in sta_names_p
                    ]
                p["sta_elev_upload"] = [
                    float(elevs_by_name.get(n, {}).get("elev_upload", 0.0))
                    for n in sta_names_p
                ]

        if not profiles:
            return (no_update, "", "", True,
                    "No profile data found. Try 'Skin-depth pseudo' source.")

        all_rho = np.concatenate([
            np.asarray(p["rho"]).ravel() for p in profiles.values()
        ])
        all_rho = all_rho[np.isfinite(all_rho) & (all_rho > 0)]
        if all_rho.size:
            log_lo   = float(np.log10(all_rho.min()))
            log_hi   = float(np.log10(all_rho.max()))
            all_z    = np.concatenate([
                np.abs(np.asarray(p["z"])).ravel() for p in profiles.values()
            ])
            z_max    = float(all_z.max()) if all_z.size else 2000.0
            hint = (f"Data: ρ {all_rho.min():.0f}–{all_rho.max():.0f} Ω·m  "
                    f"({log_lo:.1f}–{log_hi:.1f} log₁₀)  ·  "
                    f"depth 0–{z_max:.0f} m  ·  {len(profiles)} profile(s)")
        else:
            log_lo = log_hi = 0.0
            hint   = ""

        grid_data = {
            name: {
                "x":            np.asarray(p["x"]).tolist(),
                "z":            np.asarray(p["z"]).tolist(),
                "rho":          np.asarray(p["rho"]).tolist(),
                "sta_x":        p.get("sta_x",        []),
                "sta_elev":       p.get("sta_elev",       []),
                "sta_elev_raw":   p.get("sta_elev_raw",   []),
                "sta_elev_corr":  p.get("sta_elev_corr",  []),
                "sta_elev_upload":p.get("sta_elev_upload", []),
                "sta_names":    p.get("sta_names",    []),
                "sta_lat":      p.get("sta_lat",      []),
                "sta_lon":      p.get("sta_lon",      []),
            }
            for name, p in profiles.items()
        }

        store = {
            "profiles":   grid_data,
            "n_profiles": len(profiles),
            "src":        data_src,
            "log_lo":     log_lo,
            "log_hi":     log_hi,
            "topo_corr":  elev_corr_store,
            "topo_raw":   elev_raw_store,
        }
        return store, hint, "", False, ""

    # ── Show / hide upload widget based on source selection ─────────────────
    @app.callback(
        Output("map3d-topo-upload-wrap", "style"),
        Input(IDs.MAP3D_TOPO_SRC, "value"),
        prevent_initial_call=False,
    )
    def _toggle_upload_widget(src):
        if src == "upload":
            return {"display": "block", "marginBottom": "4px"}
        return {"display": "none"}

    # ── Parse topography file and populate store ─────────────────────────────
    @app.callback(
        Output(IDs.MAP3D_TOPO_UPLOAD_STORE, "data"),
        Output(IDs.MAP3D_TOPO_UPLOAD_INFO,  "children"),
        Input(IDs.MAP3D_TOPO_UPLOAD,  "contents"),
        State(IDs.MAP3D_TOPO_UPLOAD,  "filename"),
        prevent_initial_call=True,
    )
    def _parse_topo_file(contents, filename):
        if not contents:
            raise PreventUpdate
        records = _parse_topo_upload(contents, filename)
        if not records:
            msg = html.Span(
                [html.I(className="bi bi-x-circle-fill me-1",
                        style={"color": "var(--red)"}),
                 f"Could not parse {filename}. "
                 "Expected columns: station + elevation."],
                className="small",
            )
            return None, msg
        n = len(records)
        msg = html.Span(
            [html.I(className="bi bi-check-circle-fill me-1",
                    style={"color": "var(--green)"}),
             f"{filename} — {n} station{'s' if n != 1 else ''} loaded. "
             "Re-generate 3D to apply."],
            className="small",
        )
        return records, msg


def _register_display(app) -> None:
    """Re-render whenever stored grid or ANY visual control changes."""
    @app.callback(
        Output(IDs.MAP3D_GRAPH, "figure"),
        Output(IDs.MAP3D_INFO,  "children",  allow_duplicate=True),
        Output(IDs.TOAST_ERROR, "is_open",   allow_duplicate=True),
        Output(IDs.TOAST_BODY,  "children",  allow_duplicate=True),
        Input(IDs.MAP3D_GRID_STORE,  "data"),
        Input(IDs.MAP3D_ACTIVE_MODE, "data"),
        # Universal filters
        Input(IDs.MAP3D_RHO_LO,      "value"),
        Input(IDs.MAP3D_RHO_HI,      "value"),
        Input(IDs.MAP3D_DEPTH_LO,    "value"),
        Input(IDs.MAP3D_DEPTH_HI,    "value"),
        # Color & rendering
        Input(IDs.MAP3D_CMAP,        "value"),
        Input(IDs.MAP3D_SCALE,       "value"),
        Input(IDs.MAP3D_VMIN,        "value"),
        Input(IDs.MAP3D_VMAX,        "value"),
        Input(IDs.MAP3D_OPACITY,     "value"),
        Input(IDs.MAP3D_CONTOURS,    "value"),
        Input(IDs.MAP3D_N_CONTOURS,  "value"),
        # Fence
        Input(IDs.MAP3D_LINE_SPACING,"value"),
        Input(IDs.MAP3D_AZIMUTH,     "value"),
        Input(IDs.MAP3D_PANEL_THICK, "value"),
        # Block
        Input(IDs.MAP3D_ISOVALUE_LO, "value"),
        Input(IDs.MAP3D_ISOVALUE_HI, "value"),
        Input(IDs.MAP3D_N_SURFACES,  "value"),
        Input(IDs.MAP3D_CLIP_X,      "value"),
        Input(IDs.MAP3D_CLIP_Y,      "value"),
        Input(IDs.MAP3D_CLIP_Z,      "value"),
        # Depth slices
        Input(IDs.MAP3D_N_SLICES,    "value"),
        # Annotations
        Input(IDs.MAP3D_SHOW_LABELS, "value"),
        Input(IDs.MAP3D_TITLE,       "value"),
        Input(IDs.STORE_THEME,       "data"),
        # Topography controls (Phase 1 + 2 + 3)
        Input(IDs.MAP3D_TOPO_SRC,      "value"),
        Input(IDs.MAP3D_TOPO_OPACITY,  "value"),
        Input(IDs.MAP3D_TOPO_STATIONS, "value"),
        Input(IDs.MAP3D_TOPO_APPLY,    "value"),
        # Station marker style
        Input(IDs.MAP3D_STA_SYMBOL,    "value"),
        Input(IDs.MAP3D_STA_SIZE,      "value"),
        Input(IDs.MAP3D_STA_COLOR,     "value"),
        prevent_initial_call=True,
    )
    def display_3d(
        grid_store, mode,
        rho_lo, rho_hi, depth_lo, depth_hi,
        cmap, scale, vmin, vmax,
        opacity, contours, n_contours,
        line_spacing, azimuth, panel_thick,
        iso_lo, iso_hi, n_surfaces, clip_x, clip_y, clip_z,
        n_slices,
        show_labels, title, theme,
        topo_src, topo_opacity, topo_stations, topo_apply,
        sta_symbol_ui, sta_size_ui, sta_color_ui,
    ):
        dark = (theme or "dark") == "dark"
        if not grid_store or not grid_store.get("profiles"):
            return (
                _empty_3d_fig(dark=dark),
                _info("Select a data source, then click Load & Generate 3D.",
                      "secondary"),
                False,
                "",
            )

        profiles = {
            name: {
                "x":             np.asarray(p["x"]),
                "z":             np.asarray(p["z"]),
                "rho":           np.asarray(p["rho"]),
                "sta_x":         np.asarray(p.get("sta_x",         []), float),
                "sta_elev":        np.asarray(p.get("sta_elev",        []), float),
                "sta_elev_raw":    np.asarray(p.get("sta_elev_raw",    []), float),
                "sta_elev_corr":   np.asarray(p.get("sta_elev_corr",   []), float),
                "sta_elev_upload": np.asarray(p.get("sta_elev_upload",  []), float),
                "sta_names":     p.get("sta_names", []),
                "sta_lat":       p.get("sta_lat",   []),
                "sta_lon":       p.get("sta_lon",   []),
            }
            for name, p in grid_store["profiles"].items()
        }

        mode      = mode       or _DEFAULT_MODE
        cmap      = cmap       or "RdYlBu_r"
        vmin      = float(vmin or 1)
        vmax      = float(vmax or 10_000)
        rho_lo_f  = float(rho_lo  or 1)
        rho_hi_f  = float(rho_hi  or 1e8)
        d_lo_f    = float(depth_lo or 0)
        d_hi_f    = float(depth_hi or 1e6)
        opacity   = float(opacity  if opacity   is not None else 0.85)
        contours  = bool(contours)
        n_contours= int(n_contours or 15)
        title     = title or ""
        n_prof    = len(profiles)
        src       = grid_store.get("src", "pseudo")

        # Topo kwargs
        topo_src_v   = topo_src     or "none"
        topo_opac_v  = float(topo_opacity  if topo_opacity  is not None else 0.70)
        show_sta_v   = bool(topo_stations) if topo_stations is not None else True
        apply_topo_v = bool(topo_apply)    if topo_apply    is not None else False

        # Station marker style — API convention as default (api/station.py)
        _src = grid_store.get("src", "pseudo")
        # Plotly 3D has no triangle symbols; use diamond as closest match
        # (mirrors api/station.py: facecolor=black→filled, facecolor=white→open)
        _default_sym = ("diamond" if _src == "inversion" else "diamond-open")
        sta_sym_v   = (_default_sym if (sta_symbol_ui or "auto") == "auto"
                       else str(sta_symbol_ui))
        sta_size_v  = int(sta_size_ui)  if sta_size_ui  is not None else 8
        sta_color_v = str(sta_color_ui) if sta_color_ui else "black"

        # Build filter description for context bar
        rho_desc = (f"ρ {rho_lo_f:.0f}–{rho_hi_f:.0f} Ω·m"
                    if rho_hi_f < 1e7 else "all ρ")
        dep_desc = (f"depth {d_lo_f:.0f}–{d_hi_f:.0f} m"
                    if d_hi_f < 1e5 else "full depth")

        try:
            if mode == "fence":
                fig = _build_fence_fig(
                    profiles=profiles,
                    line_spacing=float(line_spacing or 1.0),
                    azimuth_deg=float(azimuth or 0),
                    panel_thick=float(panel_thick or 50),
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    opacity=opacity, contours=contours, n_contours=n_contours,
                    show_labels=bool(show_labels), title=title,
                    rho_lo=rho_lo_f, rho_hi=rho_hi_f,
                    depth_lo=d_lo_f, depth_hi=d_hi_f,
                    dark=dark,
                    topo_src=topo_src_v, topo_opacity=topo_opac_v,
                    show_sta_markers=show_sta_v, apply_topo=apply_topo_v,
                    sta_symbol=sta_sym_v, sta_size=sta_size_v,
                    sta_color=sta_color_v,
                )
                topo_tag = f" · topo: {topo_src_v}" if topo_src_v != "none" else ""
                info = _info(
                    f"Fence · {n_prof} profiles · {rho_desc} · {dep_desc} · src: {src}{topo_tag}",
                    "success",
                )

            elif mode == "block":
                x_arr, y_arr, z_arr, rho_3d = _assemble_3d_grid(
                    profiles, line_spacing=float(line_spacing or 1.0)
                )
                iso_lo_f = float(iso_lo if iso_lo is not None else
                                 grid_store.get("log_lo", 0.5))
                iso_hi_f = float(iso_hi if iso_hi is not None else
                                 grid_store.get("log_hi", 3.0))
                fig = _build_block_fig(
                    x_arr=x_arr, y_arr=y_arr, z_arr=z_arr, rho_3d=rho_3d,
                    cmap=cmap, vmin=vmin, vmax=vmax, opacity=opacity,
                    isovalue_lo=iso_lo_f, isovalue_hi=iso_hi_f,
                    n_surfaces=int(n_surfaces or 10),
                    clip_x=float(clip_x or 1), clip_y=float(clip_y or 1),
                    clip_z=float(clip_z or 1),
                    contours=contours, title=title,
                    rho_lo=rho_lo_f, rho_hi=rho_hi_f,
                    depth_lo=d_lo_f, depth_hi=d_hi_f,
                    dark=dark,
                    topo_src=topo_src_v, topo_opacity=topo_opac_v,
                    profiles=profiles,
                )
                dlo = grid_store.get("log_lo", iso_lo_f)
                dhi = grid_store.get("log_hi", iso_hi_f)
                topo_tag = f" · topo: {topo_src_v}" if topo_src_v != "none" else ""
                info = _info(
                    f"Volume · {n_prof} profiles · {rho_desc} · {dep_desc} · "
                    f"data log: {dlo:.1f}–{dhi:.1f}{topo_tag}",
                    "success",
                )

            elif mode == "depth":
                x_arr, y_arr, z_arr, rho_3d = _assemble_3d_grid(
                    profiles, line_spacing=float(line_spacing or 1.0)
                )

                def rho_at_depth(d):
                    idx = int(np.argmin(np.abs(z_arr - d)))
                    return rho_3d[:, :, idx] if rho_3d.ndim == 3 else rho_3d

                fig = _build_depth_slices_fig(
                    x_arr=x_arr, y_arr=y_arr,
                    depth_lo=d_lo_f, depth_hi=d_hi_f,
                    n_slices=int(n_slices or 8),
                    rho_fn=rho_at_depth,
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    opacity=opacity, contours=contours, n_contours=n_contours,
                    show_labels=bool(show_labels), title=title,
                    rho_lo=rho_lo_f, rho_hi=rho_hi_f,
                    dark=dark,
                    topo_src=topo_src_v, topo_opacity=topo_opac_v,
                    show_sta_markers=show_sta_v, apply_topo=apply_topo_v,
                    profiles=profiles,
                )
                topo_tag = f" · topo: {topo_src_v}" if topo_src_v != "none" else ""
                info = _info(
                    f"Depth slices · {int(n_slices or 8)} levels · "
                    f"{rho_desc} · {dep_desc} · {n_prof} profiles{topo_tag}",
                    "success",
                )

            else:
                raise ValueError(f"Unknown mode: {mode}")

        except Exception as exc:
            return (_empty_3d_fig(),
                    _info(f"Render error: {exc}", "danger"),
                    True, str(exc))

        return fig, info, False, ""


def _register_export_html(app) -> None:
    @app.callback(
        Output(IDs.MAP3D_DOWNLOAD,       "data"),
        Input(IDs.BTN_MAP3D_EXPORT_HTML, "n_clicks"),
        State(IDs.MAP3D_GRAPH,           "figure"),
        State(IDs.MAP3D_TITLE,           "value"),
        prevent_initial_call=True,
    )
    def export_html(n_clicks, figure_dict, title):
        if not n_clicks or not figure_dict:
            raise PreventUpdate
        try:
            import plotly.io as pio
            fig      = go.Figure(figure_dict)
            html_str = pio.to_html(fig, full_html=True, include_plotlyjs="cdn")
            fname    = f"{(title or 'map3d').replace(' ', '_')}.html"
            return dcc.send_string(html_str, fname)
        except Exception:
            raise PreventUpdate


# ── Helpers ───────────────────────────────────────────────────────────────────

def _assemble_3d_grid(profiles: dict, line_spacing: float = 1.0):
    line_names = list(profiles.keys())
    n_lines    = len(line_names)
    real_offsets = _line_real_offsets(profiles)
    ref        = profiles[line_names[0]]
    x_arr      = np.asarray(ref["x"])
    z_arr      = np.asarray(ref["z"])
    n_x, n_z   = len(x_arr), len(z_arr)

    rho_3d = np.zeros((n_lines, n_x, n_z))
    for i, name in enumerate(line_names):
        p     = profiles[name]
        rho_i = np.asarray(p["rho"])
        if rho_i.shape == (n_z, n_x):
            rho_i = rho_i.T
        if rho_i.shape != (n_x, n_z):
            x_i = np.asarray(p["x"])
            z_i = np.asarray(p["z"])
            from scipy.interpolate import RegularGridInterpolator
            try:
                interp = RegularGridInterpolator(
                    (z_i, x_i), np.asarray(p["rho"]),
                    bounds_error=False, fill_value=np.nan,
                )
                xi    = np.array(np.meshgrid(z_arr, x_arr, indexing="ij")).T.reshape(-1, 2)
                rho_i = interp(xi).reshape(n_x, n_z)
            except Exception:
                rho_i = np.full((n_x, n_z), np.nanmean(rho_i))
        rho_3d[i] = rho_i

    y_arr = np.array(
        [
            resolve_offset(name, i, real_offsets, 1000.0, line_spacing)
            for i, name in enumerate(line_names)
        ],
        dtype=float,
    )
    return x_arr, y_arr, z_arr, rho_3d


def _info(msg: str, color: str = "secondary"):
    _bg = {
        "success":   "#1a7a45",
        "warning":   "#8a6200",
        "danger":    "#b02030",
        "secondary": "#4a5568",
    }
    return html.Span(
        msg,
        style={
            "fontSize":        "11px",
            "color":           "#ffffff",
            "backgroundColor": _bg.get(color, "#4a5568"),
            "padding":         "3px 8px",
            "borderRadius":    "6px",
            "display":         "inline-block",
            "lineHeight":      "1.6",
            "whiteSpace":      "normal",
        },
    )
