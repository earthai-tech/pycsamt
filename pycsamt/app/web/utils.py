# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
web/utils.py — shared helpers for the Dash web application.

No Qt imports.  Functions here are used by both layout.py and
callbacks.py.
"""

from __future__ import annotations

import base64
import io

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# ──────────────────────────────────────────────────────────────────────────────
# Plotly colorscale compatibility
# ──────────────────────────────────────────────────────────────────────────────

# Matplotlib colourmap names that are NOT valid Plotly named colorscales.
# Map them to the nearest Plotly equivalent so Scattermap never crashes.
_PLOTLY_CMAP_REMAP: dict[str, str] = {
    "coolwarm": "balance",
    "seismic": "rdbu",
    "terrain": "earth",
    "YlOrRd": "ylorrd",
    "RdBu_r": "rdbu_r",
    "copper": "thermal",
    "gnuplot2": "turbid",
    "tab10": "plotly3",
    "tab20": "plasma",
    "bwr": "balance",
    "PiYG": "piyg",
    "PRGn": "prgn",
    "RdYlGn": "rdylgn",
    "RdYlBu": "rdylbu",
    "Spectral": "spectral",
}


def _to_plotly_cmap(cmap: str | None, fallback: str = "plasma") -> str:
    """Return a Plotly-valid colorscale name, remapping matplotlib-only names."""
    if not cmap:
        return fallback
    return _PLOTLY_CMAP_REMAP.get(cmap, cmap)


# ──────────────────────────────────────────────────────────────────────────────
# Active-lines filtering
# ──────────────────────────────────────────────────────────────────────────────


def filter_sites_by_lines(sites, station_records: list, active_lines: list):
    """Return a Sites containing only stations from *active_lines*.

    Returns ``None`` when no stations remain after filtering (signals "no
    active lines" to the caller so it can show a warning figure).
    Returns *sites* unchanged when line information is unavailable.
    """
    if not active_lines or not station_records:
        return sites
    try:
        import pandas as pd

        df = pd.DataFrame(station_records)
        if "Line" not in df.columns or "ID" not in df.columns:
            return sites
        active_ids = set(df[df["Line"].isin(active_lines)]["ID"].tolist())
        if not active_ids:
            return None

        from pycsamt.emtools._core import (
            _iter_items,
            _name,
            ensure_sites,
        )

        items = list(_iter_items(sites))
        filtered = [
            getattr(ed, "edi", ed)
            for i, ed in enumerate(items)
            if _name(ed, i) in active_ids
        ]
        if not filtered:
            return None
        return ensure_sites(
            filtered,
            recursive=True,
            on_dup="replace",
            strict=False,
            verbose=0,
        )
    except Exception:
        return sites  # fail-safe: return unfiltered


def no_active_lines_src(dark: bool = True) -> str:
    """Data-URI PNG shown when all lines are muted."""
    import matplotlib.pyplot as plt

    bg = "#1e1e2e" if dark else "#eff1f5"
    fg = "#cdd6f4" if dark else "#4c4f69"
    acc = "#f9e2af" if dark else "#df8e1d"
    fig, ax = plt.subplots(figsize=(9, 4))
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)
    ax.set_axis_off()
    ax.text(
        0.5,
        0.60,
        "No active lines",
        ha="center",
        va="center",
        transform=ax.transAxes,
        color=acc,
        fontsize=18,
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.40,
        "All survey lines are muted.\n"
        "Enable at least one line in the Lines panel ( navbar → Lines ).",
        ha="center",
        va="center",
        transform=ax.transAxes,
        color=fg,
        fontsize=11,
        linespacing=1.7,
    )
    src = fig_to_src(fig)
    plt.close(fig)
    return src


# ──────────────────────────────────────────────────────────────────────────────
# Matplotlib ↔ web helpers
# ──────────────────────────────────────────────────────────────────────────────

#: Catppuccin Mocha palette — matches dark_theme.qss
_DARK_RC: dict = {
    "axes.facecolor": "#181825",
    "figure.facecolor": "#1e1e2e",
    "savefig.facecolor": "#1e1e2e",
    "axes.edgecolor": "#45475a",
    "axes.labelcolor": "#cdd6f4",
    "text.color": "#cdd6f4",
    "xtick.color": "#a6adc8",
    "ytick.color": "#a6adc8",
    "grid.color": "#313244",
    "grid.linestyle": "--",
    "grid.alpha": 0.35,
    "axes.grid": True,
    "lines.color": "#89b4fa",
    "legend.facecolor": "#1e1e2e",
    "legend.edgecolor": "#45475a",
    "legend.labelcolor": "#cdd6f4",
}

_LIGHT_RC: dict = {
    "axes.facecolor": "#eff1f5",
    "figure.facecolor": "#e6e9ef",
    "savefig.facecolor": "#e6e9ef",
    "axes.edgecolor": "#bcc0cc",
    "axes.labelcolor": "#4c4f69",
    "text.color": "#4c4f69",
    "xtick.color": "#6c6f85",
    "ytick.color": "#6c6f85",
    "grid.color": "#ccd0da",
    "grid.linestyle": "--",
    "grid.alpha": 0.5,
    "axes.grid": True,
    "lines.color": "#1e66f5",
}

#: Publication-quality override — white background, dark text, suitable for papers.
PUB_RC: dict = {
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.edgecolor": "#222222",
    "axes.labelcolor": "#111111",
    "text.color": "#111111",
    "xtick.color": "#222222",
    "ytick.color": "#222222",
    "xtick.labelcolor": "#222222",
    "ytick.labelcolor": "#222222",
    "grid.color": "#cccccc",
    "grid.linestyle": "--",
    "grid.alpha": 0.6,
    "axes.grid": True,
    "lines.color": "#1565c0",
    "legend.facecolor": "white",
    "legend.edgecolor": "#cccccc",
    "legend.labelcolor": "#111111",
    "axes.titlecolor": "#111111",
}


def apply_web_dark_theme() -> None:
    mpl.rcParams.update(_DARK_RC)


def apply_web_light_theme() -> None:
    mpl.rcParams.update(_LIGHT_RC)


def fig_to_src(fig, dpi: int = 130) -> str:
    """
    Convert a matplotlib Figure to a base64 ``data:image/png`` URI.

    The figure is closed after conversion.
    """
    buf = io.BytesIO()
    try:
        fig.savefig(
            buf,
            format="png",
            dpi=dpi,
            bbox_inches="tight",
            facecolor=fig.get_facecolor(),
        )
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("ascii")
    finally:
        plt.close(fig)
        buf.close()
    return f"data:image/png;base64,{encoded}"


def empty_src(dark: bool = True) -> str:
    """Return a 1×1 transparent PNG placeholder.

    The image is fully transparent so the CSS container background
    shows through, which is theme-aware via CSS variables (dark/light).
    The ``dark`` parameter is kept for API compatibility but ignored.
    """
    # Minimal valid 1×1 RGBA transparent PNG (hard-coded bytes, no Matplotlib)
    _TRANSPARENT_1X1 = (
        b"\x89PNG\r\n\x1a\n"  # PNG signature
        b"\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00"
        b"\x00\x01\x08\x06\x00\x00\x00\x1f\x15\xc4\x89"  # 1×1, RGBA
        b"\x00\x00\x00\nIDATx\x9cc\x00\x01\x00\x00\x05"
        b"\x00\x01\r\n-\xb4\x00\x00\x00\x00IEND\xaeB`\x82"
    )
    enc = base64.b64encode(_TRANSPARENT_1X1).decode("ascii")
    return f"data:image/png;base64,{enc}"


# ──────────────────────────────────────────────────────────────────────────────
# Plotly map builder
# ──────────────────────────────────────────────────────────────────────────────

_OVERLAY_LABELS = {
    # ── Survey station overlays ───────────────────────────────────────────────
    "Index": "Station index",
    "Apparent Resistivity": "ρₐ (Ω·m)",
    "Phase": "Phase (°)",
    "Quality Score": "Quality score",
    "Z-strike": "Strike (°)",
    "N_freq": "N frequencies",
    "Elevation": "Elevation (m)",
    # ── Computed AMT overlays (stored as temp columns) ────────────────────────
    "_rho_map": "log₁₀(ρₐ)  [Ω·m]",  # app. resistivity at freq
    "_depth_map": "Skin depth  (m)",  # Bostick skin depth
    # ── Inversion-derived overlays ────────────────────────────────────────────
    "_inv_rho": "log₁₀(ρ)  [Ω·m]",  # true resistivity from inversion
    "_inv_depth": "Depth  (m)",
}


def _colorbar_title(overlay: str) -> str:
    """Return a properly formatted colorbar/hover label for *overlay*."""
    return _OVERLAY_LABELS.get(overlay, overlay.lstrip("_").replace("_", " ").title())


# Per-line colour palette (Catppuccin Mocha accent colours)
_LINE_COLORS = [
    "#89b4fa",  # blue
    "#a6e3a1",  # green
    "#fab387",  # peach
    "#f38ba8",  # red
    "#cba6f7",  # mauve
    "#94e2d5",  # teal
    "#f9e2af",  # yellow
    "#89dceb",  # sky
]


def _auto_zoom(lats: np.ndarray, lons: np.ndarray) -> tuple[float, float, float, int]:
    """Return (center_lat, center_lon, bearing, zoom) for a Scattermap layout."""
    import math

    lats = lats[np.isfinite(lats)]
    lons = lons[np.isfinite(lons)]
    if len(lats) == 0 or len(lons) == 0:
        return 0.0, 0.0, 0, 3
    clat = float((lats.min() + lats.max()) / 2)
    clon = float((lons.min() + lons.max()) / 2)
    lat_span = float(lats.max() - lats.min())
    lon_span = float(lons.max() - lons.min())
    max_span = max(lat_span, lon_span, 1e-6)
    zoom = max(2, min(14, int(math.floor(8.0 - math.log2(max_span + 0.001)))))
    return clat, clon, 0, zoom


def build_station_map(
    df: pd.DataFrame,
    selected_id: str | None = None,
    overlay: str = "Index",
    dark: bool = True,
    line_filter: list | None = None,
    muted_lines: list | None = None,
    marker_size: int = 10,
    basemap_style: str | None = None,
    cmap: str = "plasma",
    opacity: float = 0.92,
    show_labels: bool = True,
    show_profiles: bool = True,
    uirevision: str = "map",
    bearing: float = 0.0,
) -> go.Figure:
    """
    Build a Plotly Scattermap figure for the station map.

    Parameters
    ----------
    df : DataFrame
        Columns: ID, Latitude, Longitude, Line (+ any overlay column).
    selected_id : str, optional
        Station to highlight with a larger gold marker.
    overlay : str
        Column to colour station markers by.
    dark : bool
        Use dark tile style (overridden by basemap_style if given).
    line_filter : list of str, optional
        If given, only stations on these survey lines are shown.
    marker_size : int
        Base marker size in pixels (default 10).
    basemap_style : str, optional
        Explicit Mapbox tile style (e.g. ``"carto-positron"``).
    cmap : str
        Plotly colorscale name (default ``"plasma"``).
    opacity : float
        Marker opacity 0–1 (default 0.92).
    show_labels : bool
        Render station name labels next to markers.
    show_profiles : bool
        Render survey-line polylines.
    """
    if df is None or df.empty:
        return _empty_map(dark)

    df = df.copy()
    aliases = {
        "ID": ("ID", "id", "station", "Station", "name", "Name"),
        "Latitude": ("Latitude", "latitude", "lat", "Lat"),
        "Longitude": ("Longitude", "longitude", "lon", "Lon", "long", "Long"),
        "Line": ("Line", "line", "profile", "Profile"),
        "Elevation": ("Elevation", "elevation", "elev", "Elev"),
    }
    for canonical, candidates in aliases.items():
        if canonical not in df.columns:
            match = next((c for c in candidates if c in df.columns), None)
            if match:
                df[canonical] = df[match]
    if not {"ID", "Latitude", "Longitude"}.issubset(df.columns):
        return _empty_map(dark)

    lats = pd.to_numeric(df["Latitude"], errors="coerce").values
    lons = pd.to_numeric(df["Longitude"], errors="coerce").values
    ids = df["ID"].astype(str).values

    valid = ~(np.isnan(lats) | np.isnan(lons))
    if not valid.any():
        return _empty_map(dark)

    lats, lons, ids = lats[valid], lons[valid], ids[valid]
    df_v = df[valid].reset_index(drop=True)

    # Apply line filter
    if line_filter and "Line" in df_v.columns:
        requested_lines = [str(line) for line in line_filter]
        mask = df_v["Line"].astype(str).isin(requested_lines)
        df_v = df_v[mask].reset_index(drop=True)
        lats = pd.to_numeric(df_v["Latitude"], errors="coerce").values
        lons = pd.to_numeric(df_v["Longitude"], errors="coerce").values
        ids = df_v["ID"].astype(str).values
        if len(ids) == 0:
            return _empty_map(dark)

    bg_color = "#1e1e2e" if dark else "#eff1f5"
    text_col = "#cdd6f4" if dark else "#4c4f69"
    sel_col = "#f9e2af"  # Catppuccin yellow — selected station

    # Overlay colour values
    if overlay in df_v.columns:
        colour_vals = pd.to_numeric(df_v[overlay], errors="coerce").fillna(0).values
    else:
        colour_vals = np.arange(len(df_v), dtype=float)

    sel_size = int(marker_size * 1.8)
    sizes = (
        np.where(ids == selected_id, sel_size, marker_size).tolist()
        if selected_id
        else [marker_size] * len(ids)
    )
    [sel_col if sid == selected_id else "#ffffff" for sid in ids]
    [2.5 if sid == selected_id else 0.8 for sid in ids]

    hover_label = _colorbar_title(overlay)

    # Guard: clamp opacity; remap any matplotlib-only colorscale name to Plotly
    opacity = float(max(0.1, min(1.0, opacity)))
    cmap = _to_plotly_cmap(cmap, "plasma")

    fig = go.Figure()

    _muted = set(muted_lines or [])

    # ── Survey line polylines (one trace per line) ──────────────────────────
    if show_profiles and "Line" in df_v.columns:
        all_lines_df = df_v.copy()  # includes muted — shown dimmed
        unique_lines = [ln for ln in all_lines_df["Line"].unique() if ln and ln != "—"]
        for idx, line_name in enumerate(unique_lines):
            is_muted = line_name in _muted
            grp = all_lines_df[all_lines_df["Line"] == line_name].copy()
            lon_span = grp["Longitude"].max() - grp["Longitude"].min()
            lat_span = grp["Latitude"].max() - grp["Latitude"].min()
            sort_col = "Longitude" if lon_span >= lat_span else "Latitude"
            grp = grp.sort_values(sort_col)
            lcol = "#45475a" if is_muted else _LINE_COLORS[idx % len(_LINE_COLORS)]
            lw = 1 if is_muted else 2
            lname = f"{line_name} (muted)" if is_muted else str(line_name)
            fig.add_trace(
                go.Scattermap(
                    lat=grp["Latitude"].tolist(),
                    lon=grp["Longitude"].tolist(),
                    mode="lines",
                    line=dict(color=lcol, width=lw),
                    opacity=0.3 if is_muted else 1.0,
                    name=lname,
                    hoverinfo="skip",
                    showlegend=len(unique_lines) > 1,
                )
            )

    # ── Muted station ghost markers ─────────────────────────────────────────
    if _muted and "Line" in df_v.columns:
        muted_mask = df_v["Line"].astype(str).isin({str(line) for line in _muted})
        df_muted = df_v[muted_mask]
        if not df_muted.empty:
            m_lats = pd.to_numeric(df_muted["Latitude"], errors="coerce").values
            m_lons = pd.to_numeric(df_muted["Longitude"], errors="coerce").values
            m_ids = df_muted["ID"].values
            fig.add_trace(
                go.Scattermap(
                    lat=m_lats.tolist(),
                    lon=m_lons.tolist(),
                    text=m_ids.tolist(),
                    mode="markers",
                    marker=dict(size=6, color="#6c7086", opacity=0.3),
                    name="Muted lines",
                    hovertemplate="<b>%{text}</b> (muted)<extra></extra>",
                    showlegend=False,
                )
            )
        df_v = df_v[~muted_mask].reset_index(drop=True)
        lats = pd.to_numeric(df_v["Latitude"], errors="coerce").values
        lons = pd.to_numeric(df_v["Longitude"], errors="coerce").values
        ids = df_v["ID"].astype(str).values
        if overlay in df_v.columns:
            colour_vals = pd.to_numeric(df_v[overlay], errors="coerce").fillna(0).values
        else:
            colour_vals = np.arange(len(df_v), dtype=float)
        sel_size = int(marker_size * 1.8)
        sizes = (
            np.where(ids == selected_id, sel_size, marker_size).tolist()
            if selected_id
            else [marker_size] * len(ids)
        )

    # ── Station markers ─────────────────────────────────────────────────────
    marker_mode = "markers+text" if show_labels else "markers"
    fig.add_trace(
        go.Scattermap(
            lat=lats.tolist(),
            lon=lons.tolist(),
            text=ids.tolist(),
            customdata=ids.tolist(),
            mode=marker_mode,
            textposition="top right",
            textfont=dict(size=9, color=text_col),
            marker=dict(
                size=sizes,
                color=colour_vals.tolist(),
                colorscale=cmap,
                showscale=True,
                colorbar=dict(
                    title=dict(
                        text=hover_label,
                        side="right",
                        font=dict(color=text_col, size=11),
                    ),
                    tickfont=dict(color=text_col, size=9),
                    bgcolor=bg_color,
                    bordercolor="#45475a" if dark else "#bcc0cc",
                    thickness=12,
                    len=0.7,
                    x=1.0,
                ),
                opacity=opacity,
            ),
            hovertemplate=(
                "<b>%{text}</b><br>"
                "Lat: %{lat:.4f}°<br>"
                "Lon: %{lon:.4f}°<br>"
                f"{hover_label}: %{{marker.color:.2f}}"
                "<extra></extra>"
            ),
            name="Stations",
            showlegend=False,
        )
    )

    # Selected station: gold ring marker
    if selected_id and selected_id in ids:
        sel_mask = ids == selected_id
        fig.add_trace(
            go.Scattermap(
                lat=lats[sel_mask].tolist(),
                lon=lons[sel_mask].tolist(),
                mode="markers",
                marker=dict(
                    size=int(marker_size * 2.4),
                    color=sel_col,
                    opacity=0.35,
                ),
                hoverinfo="skip",
                showlegend=False,
                name="_sel_ring",
            )
        )

    # ── Map layout ──────────────────────────────────────────────────────────
    # Auto-center on selected station when possible
    clat, clon, _, zoom = _auto_zoom(lats, lons)
    if selected_id and selected_id in ids:
        sel_mask = ids == selected_id
        sel_lat = lats[sel_mask]
        sel_lon = lons[sel_mask]
        if sel_lat.size > 0 and np.isfinite(sel_lat[0]) and np.isfinite(sel_lon[0]):
            clat = float(sel_lat[0])
            clon = float(sel_lon[0])
            zoom = max(zoom, 10)  # zoom in to the station

    _ESRI_URLS = {
        "esri-satellite": "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        "esri-topo": "https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}",
        "esri-natgeo": "https://server.arcgisonline.com/ArcGIS/rest/services/NatGeo_World_Map/MapServer/tile/{z}/{y}/{x}",
        "esri-ocean": "https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}",
        "esri-street": "https://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/{z}/{y}/{x}",
    }
    default_style = "carto-darkmatter" if dark else "carto-positron"
    tile_style = basemap_style or default_style

    bearing = float(bearing or 0.0) % 360.0

    if tile_style in _ESRI_URLS:
        map_cfg = dict(
            style="white-bg",
            layers=[
                {
                    "below": "traces",
                    "sourcetype": "raster",
                    "source": [_ESRI_URLS[tile_style]],
                    "opacity": 1.0,
                }
            ],
            center=dict(lat=clat, lon=clon),
            zoom=zoom,
            bearing=bearing,
        )
    else:
        map_cfg = dict(
            style=tile_style,
            center=dict(lat=clat, lon=clon),
            zoom=zoom,
            bearing=bearing,
        )

    fig.update_layout(
        map=map_cfg,
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor=bg_color,
        font=dict(color=text_col, size=11),
        showlegend=True,
        legend=dict(
            bgcolor="rgba(30,30,46,.75)" if dark else "rgba(239,241,245,.80)",
            bordercolor="#45475a" if dark else "#bcc0cc",
            borderwidth=1,
            font=dict(size=10, color=text_col),
            x=0.01,
            y=0.99,
            xanchor="left",
            yanchor="top",
        ),
        clickmode="event+select",
        uirevision=uirevision,
    )

    return fig


def build_contour_overlay(
    lats: np.ndarray,
    lons: np.ndarray,
    values: np.ndarray,
    *,
    cmap: str = "jet",
    n_levels: int = 12,
    interp_method: str = "cubic",
    extra_factor: float = 0.12,
    opacity: float = 0.60,
    mode: str = "filled+lines",
    show_labels: bool = False,
    line_width: float = 1.0,
    smooth_sigma: float = 1.0,
    log_scale: bool = False,
    grid_res: int = 150,
) -> dict | None:
    """
    Render a Surfer-style contour image over the station network.

    The function:
    1. Interpolates sparse station values onto a regular lat/lon grid.
    2. Optionally applies Gaussian smoothing.
    3. Renders a matplotlib filled-contour (and/or contour lines) to a
       transparent PNG.
    4. Applies an edge-fade alpha mask so the contour dissolves naturally
       near the survey boundary.

    Returns
    -------
    dict with keys:
        ``image_b64``   – ``data:image/png;base64,...`` PNG string
        ``coordinates`` – 4-corner [lon, lat] list for mapbox image layer
        ``vmin``        – minimum mapped value (original scale)
        ``vmax``        – maximum mapped value (original scale)
    or ``None`` if the overlay cannot be built (too few points, bad data).

    Parameters
    ----------
    lats, lons, values : ndarray
        Station coordinates and the scalar field to contour.
    cmap : str
        Matplotlib colormap name.
    n_levels : int
        Number of contour levels (4–30).
    interp_method : str
        ``"cubic"`` / ``"linear"`` / ``"nearest"`` (scipy griddata).
    extra_factor : float
        Fractional bbox expansion for mild extrapolation (0 = none).
    opacity : float
        Overall fill opacity 0–1.  Edge fade further reduces transparency
        near the boundary.
    mode : str
        ``"filled+lines"``, ``"filled"``, or ``"lines"``.
    show_labels : bool
        Annotate contour lines with numeric values.
    line_width : float
        Contour line width in pts.
    smooth_sigma : float
        Gaussian sigma applied to the interpolated grid (0 = off).
    log_scale : bool
        Work in log₁₀ space (useful for resistivity / depth).
    grid_res : int
        Number of grid cells per axis (e.g. 150 → 150×150 grid).
    """
    import base64
    import io

    import matplotlib.colors as mcolors
    import numpy as np
    from matplotlib.backends.backend_agg import (
        FigureCanvasAgg,
    )
    from matplotlib.figure import Figure
    from PIL import Image as _PILImage
    from scipy.interpolate import griddata
    from scipy.ndimage import gaussian_filter
    from scipy.spatial import cKDTree

    # ── Validate input ────────────────────────────────────────────────────
    lats = np.asarray(lats, dtype=float).ravel()
    lons = np.asarray(lons, dtype=float).ravel()
    values = np.asarray(values, dtype=float).ravel()

    mask = np.isfinite(lats) & np.isfinite(lons) & np.isfinite(values)
    lats, lons, values = lats[mask], lons[mask], values[mask]

    if len(values) < 3:
        return None

    if log_scale:
        pos = values > 0
        if pos.sum() < 3:
            return None
        lats, lons, values = lats[pos], lons[pos], values[pos]
        plot_values = np.log10(values)
    else:
        plot_values = values.copy()

    vmin_plot, vmax_plot = (
        float(np.nanmin(plot_values)),
        float(np.nanmax(plot_values)),
    )
    if vmin_plot >= vmax_plot:
        return None

    # ── Build interpolation grid ──────────────────────────────────────────
    lat_min, lat_max = lats.min(), lats.max()
    lon_min, lon_max = lons.min(), lons.max()

    lat_span = lat_max - lat_min or 0.001
    lon_span = lon_max - lon_min or 0.001

    lat_pad = lat_span * extra_factor
    lon_pad = lon_span * extra_factor

    g_lat = np.linspace(lat_min - lat_pad, lat_max + lat_pad, grid_res)
    g_lon = np.linspace(lon_min - lon_pad, lon_max + lon_pad, grid_res)
    glon_mg, glat_mg = np.meshgrid(g_lon, g_lat)  # shape (grid_res, grid_res)

    pts = np.column_stack([lons, lats])  # (N, 2) — x=lon, y=lat

    # Try requested method; fall back to linear for very sparse datasets
    method = interp_method
    if method == "cubic" and len(values) < 10:
        method = "linear"

    grid_z = griddata(pts, plot_values, (glon_mg, glat_mg), method=method)

    # Fill NaN (outside convex hull) with nearest-neighbour for extrapolation
    nan_mask = np.isnan(grid_z)
    if nan_mask.any():
        grid_z_nn = griddata(pts, plot_values, (glon_mg, glat_mg), method="nearest")
        grid_z[nan_mask] = grid_z_nn[nan_mask]

    # Gaussian smoothing
    if smooth_sigma > 0:
        grid_z = gaussian_filter(grid_z, sigma=smooth_sigma)

    # ── Alpha fade mask based on distance from nearest station ────────────
    tree = cKDTree(pts)
    dists, _ = tree.query(np.column_stack([glon_mg.ravel(), glat_mg.ravel()]))
    dists = dists.reshape(grid_res, grid_res)

    # Typical station spacing = median of each point's nearest-neighbour dist
    nn_d, _ = tree.query(pts, k=min(2, len(pts)))
    if nn_d.ndim > 1:
        nn_d = nn_d[:, -1]
    station_spacing = float(np.median(nn_d)) if len(nn_d) > 0 else 0.01
    # Alpha = 1 within 1× spacing, fades to 0 at 2× spacing + extra margin
    fade_start = station_spacing * 0.8
    fade_end = station_spacing * (2.5 + extra_factor * 6)
    alpha_mask = np.clip(
        1.0 - np.maximum(0.0, dists - fade_start) / max(fade_end - fade_start, 1e-12),
        0.0,
        1.0,
    )

    # ── Compute contour levels ────────────────────────────────────────────
    levels = np.linspace(vmin_plot, vmax_plot, n_levels)
    norm = (
        mcolors.LogNorm(vmin=10**vmin_plot, vmax=10**vmax_plot)
        if False
        else mcolors.Normalize(vmin=vmin_plot, vmax=vmax_plot)
    )

    # ── Render contour via matplotlib (Agg, no display) ──────────────────
    dpi = 100
    fig_inch = grid_res / dpi
    fig_mpl = Figure(figsize=(fig_inch, fig_inch), dpi=dpi, facecolor="none")
    canvas = FigureCanvasAgg(fig_mpl)
    ax = fig_mpl.add_axes([0, 0, 1, 1])
    ax.set_facecolor("none")
    ax.axis("off")
    ax.set_xlim(g_lon.min(), g_lon.max())
    ax.set_ylim(g_lat.min(), g_lat.max())

    do_fill = "filled" in mode
    do_lines = "lines" in mode

    if do_fill:
        ax.contourf(
            glon_mg,
            glat_mg,
            grid_z,
            levels=levels,
            cmap=cmap,
            norm=norm,
            alpha=1.0,
        )

    if do_lines:
        line_col = "white" if do_fill else None
        cs = ax.contour(
            glon_mg,
            glat_mg,
            grid_z,
            levels=levels,
            norm=norm,
            colors=line_col if line_col else None,
            cmap=None if line_col else cmap,
            linewidths=line_width,
        )
        if show_labels:
            fmt_str = "%.2g"
            if log_scale:

                def _fmt(v):
                    return f"10^{v:.1f}" if abs(v) < 4 else f"{10**v:.2g}"

                ax.clabel(cs, fmt=_fmt, fontsize=7, inline=True, inline_spacing=2)
            else:
                ax.clabel(cs, fmt=fmt_str, fontsize=7, inline=True, inline_spacing=2)

    canvas.draw()

    # Extract RGBA from canvas
    w, h = fig_inch * dpi, fig_inch * dpi
    raw = canvas.buffer_rgba()  # RGBA buffer
    img_arr = np.frombuffer(raw, dtype=np.uint8).reshape(int(h), int(w), 4).copy()

    # Apply opacity * alpha_mask to the A channel
    # Note: matplotlib sets A=0 where there's no artist (transparent bg)
    # and A=255 where contourf drew.  We multiply by (opacity * alpha_mask).
    alpha_combined = (opacity * alpha_mask * 255).astype(np.uint8)
    # Only reduce opacity where matplotlib painted (don't increase it)
    painted = img_arr[:, :, 3] > 0
    img_arr[:, :, 3] = np.where(
        painted,
        np.minimum(img_arr[:, :, 3], alpha_combined.astype(np.uint16)).astype(np.uint8),
        0,
    )

    pil_img = _PILImage.fromarray(img_arr, "RGBA")
    buf = io.BytesIO()
    pil_img.save(buf, format="PNG", optimize=False)
    image_b64 = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()

    # ── Return dict ───────────────────────────────────────────────────────
    # Mapbox image coordinates: [top-left, top-right, bottom-right, bottom-left]
    # each as [lon, lat]
    coordinates = [
        [float(g_lon.min()), float(g_lat.max())],  # top-left
        [float(g_lon.max()), float(g_lat.max())],  # top-right
        [float(g_lon.max()), float(g_lat.min())],  # bottom-right
        [float(g_lon.min()), float(g_lat.min())],  # bottom-left
    ]

    # Original-scale min/max for colorbar
    vmin_orig = float(values.min())
    vmax_orig = float(values.max())

    return {
        "image_b64": image_b64,
        "coordinates": coordinates,
        "vmin": vmin_orig,
        "vmax": vmax_orig,
        "log_scale": log_scale,
        "cmap": cmap,
    }


def _empty_map(dark: bool = True) -> go.Figure:
    bg = "#1e1e2e" if dark else "#eff1f5"
    text_col = "#585b70" if dark else "#8c8fa1"
    tile = "carto-darkmatter" if dark else "carto-positron"

    fig = go.Figure()
    fig.add_trace(
        go.Scattermap(
            lat=[],
            lon=[],
            mode="markers",
            showlegend=False,
        )
    )
    fig.add_annotation(
        text="Load survey data to display the station map",
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
        showarrow=False,
        font=dict(size=13, color=text_col),
        bgcolor="rgba(30,30,46,.65)" if dark else "rgba(239,241,245,.70)",
        borderpad=8,
    )
    fig.update_layout(
        map=dict(
            style=tile,
            center=dict(lat=20, lon=0),
            zoom=1,
        ),
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor=bg,
        showlegend=False,
        uirevision="empty",
    )
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# EDI file discovery
# ──────────────────────────────────────────────────────────────────────────────


def empty_state(
    icon: str = "bi-broadcast",
    title: str = "No data loaded",
    subtitle: str = "Load EDI files to get started",
    btn_label: str | None = None,
    btn_id: str | None = None,
):
    """
    Return a structured empty-state placeholder `html.Div`.

    Parameters
    ----------
    icon : str
        Bootstrap icon class without the ``bi-`` prefix (e.g. ``"radar"``).
    title : str
        Short headline shown under the icon.
    subtitle : str
        One-sentence description / call-to-action.
    btn_label : str, optional
        If given, a primary button is appended.
    btn_id : str, optional
        ``id`` for the optional button.
    """
    import dash_bootstrap_components as dbc
    from dash import html

    children = [
        html.I(className=f"bi bi-{icon} empty-state-icon"),
        html.P(title, className="empty-state-title"),
        html.P(subtitle, className="empty-state-sub"),
    ]
    if btn_label and btn_id:
        children.append(
            dbc.Button(
                [html.I(className="bi bi-cloud-upload me-2"), btn_label],
                id=btn_id,
                color="primary",
                size="sm",
            )
        )
    return html.Div(children, className="empty-state")


_SURVEY_EXTS = {".edi", ".avg", ".j"}


def find_edi_files(path: str) -> list[str]:
    """Recursively find survey files (.edi / .avg / .j) under *path*."""
    paths, _ = find_edi_files_by_line(path)
    return paths


def find_edi_files_by_line(
    path: str,
) -> tuple[list[str], dict[str, list[str]]]:
    """
    Recursively discover survey files (.edi / .avg / .j) under *path*,
    grouping by immediate parent sub-directory (= survey line).

    Matching is **case-insensitive** so ``.EDI`` and ``.edi`` are treated
    identically (important on Linux where ``rglob("*.edi")`` would silently
    miss upper-case variants).

    EDIs / AVGs that sit directly inside *path* (no sub-directory) are
    grouped under the root directory's own name.

    Returns
    -------
    paths : list[str]
        Flat list of every matched file, sorted by (line, filename).
    line_map : dict[str, list[str]]
        ``{line_name: [path, ...]}`` — insertion order follows directory
        sort order.  ``line_name`` is the immediate sub-directory name
        relative to *path* (or the directory's own name for root-level files).
    """
    from pathlib import Path

    p = Path(path)
    line_map: dict[str, list[str]] = {}

    def _matches(fp: Path) -> bool:
        return fp.suffix.lower() in _SURVEY_EXTS

    if p.is_file() and _matches(p):
        line_map[p.parent.name or "line"] = [str(p)]
    elif p.is_dir():
        # rglob("*") yields every file; filter by extension ourselves so the
        # match is case-insensitive on all platforms (Linux rglob is case-
        # sensitive, so "*.edi" would miss ".EDI" files).
        for survey_path in sorted(
            fp for fp in p.rglob("*") if fp.is_file() and _matches(fp)
        ):
            try:
                rel = survey_path.parent.relative_to(p)
            except ValueError:
                rel = survey_path.parent
            # Files directly under root → use the root dir's own name as line
            line_name = str(rel) if str(rel) != "." else p.name
            line_map.setdefault(line_name, []).append(str(survey_path))

    all_paths = [fp for plist in line_map.values() for fp in plist]
    return all_paths, line_map


# ──────────────────────────────────────────────────────────────────────────────
# dcc.Upload helpers
# ──────────────────────────────────────────────────────────────────────────────

_UPLOAD_EXTENSIONS = {".edi", ".avg", ".j"}


def decode_upload_to_tempdir(
    contents: list[str] | str,
    filenames: list[str] | str,
) -> tuple[list[str], str]:
    """
    Decode base64 payloads from ``dcc.Upload`` and write them to a temp dir.

    Parameters
    ----------
    contents :
        One or more ``"data:<mime>;base64,<data>"`` strings from
        ``dcc.Upload.contents``.
    filenames :
        Matching original file names from ``dcc.Upload.filename``.

    Returns
    -------
    paths : list[str]
        Absolute paths of the written temp files (recognised extensions only).
    tmpdir : str
        Temporary directory path.  Caller must remove it with
        ``shutil.rmtree(tmpdir, ignore_errors=True)`` after loading.
    """
    import base64
    import os
    import tempfile
    from pathlib import Path

    # dcc.Upload with multiple=True always returns lists, but guard anyway.
    if not isinstance(contents, list):
        contents = [contents]
        filenames = [filenames]

    tmpdir = tempfile.mkdtemp(prefix="pycsamt_upload_")
    paths: list[str] = []

    for content, name in zip(contents, filenames):
        if Path(name).suffix.lower() not in _UPLOAD_EXTENSIONS:
            continue  # skip unsupported files silently
        # Strip the data-URI header: "data:<mime>;base64,<payload>"
        try:
            _, b64 = content.split(",", 1)
            raw = base64.b64decode(b64)
        except Exception:
            continue
        fpath = os.path.join(tmpdir, name)
        with open(fpath, "wb") as fh:
            fh.write(raw)
        paths.append(fpath)

    return paths, tmpdir


def decode_folder_upload_to_tempdir(
    contents: list[str],
    filenames: list[str],
) -> tuple[list[str], str, dict[str, list[str]]]:
    """
    Decode base64 payloads that were collected by the JS folder-browse handler.

    Unlike ``decode_upload_to_tempdir``, the *filenames* list contains relative
    paths such as ``"WILLY_DATA/L18PLT/HBH01.edi"`` (or
    ``"L18PLT/HBH01.edi"`` if the root folder was stripped).  These paths are
    used to reconstruct the per-line grouping:

    * The deepest directory component that is a direct child of the root is
      treated as the **line name**.
    * Files whose relative path has no directory component (flat root) are
      grouped under the root folder's name derived from the first path found.

    Returns
    -------
    paths     : list[str]   — absolute temp-file paths, sorted by (line, name)
    tmpdir    : str         — caller must ``shutil.rmtree(tmpdir)``
    line_map  : dict        — ``{line_name: [path, ...]}``
    """
    import base64
    import os
    import tempfile
    from pathlib import Path

    if not isinstance(contents, list):
        contents = [contents]
        filenames = [filenames]

    tmpdir = tempfile.mkdtemp(prefix="pycsamt_folder_")
    line_map: dict[str, list[str]] = {}

    # Determine the common root prefix (first path component shared by all)
    # to strip it when it equals the selected folder name.
    # e.g. ["WILLY_DATA/L18PLT/X.edi", "WILLY_DATA/L22PLT/Y.edi"]
    # → strip "WILLY_DATA/" → use "L18PLT" / "L22PLT" as line names.
    def _first_part(rel: str) -> str:
        parts = Path(rel).parts
        return parts[0] if len(parts) > 1 else ""

    roots = {_first_part(fn) for fn in filenames if _first_part(fn)}
    common_root = roots.pop() if len(roots) == 1 else ""

    for content, rel_name in zip(contents, filenames):
        rel_path = Path(rel_name)
        if rel_path.suffix.lower() not in _UPLOAD_EXTENSIONS:
            continue
        try:
            _, b64 = content.split(",", 1)
            raw = base64.b64decode(b64)
        except Exception:
            continue

        # Strip the common root prefix so the line is the next component
        parts = rel_path.parts
        if common_root and parts[0] == common_root:
            parts = parts[1:]  # drop root folder name

        # Line name = first remaining directory component (or root label)
        if len(parts) > 1:
            line_name = parts[0]
        else:
            line_name = common_root or "line"

        # Write file preserving just the filename (no sub-dirs in tmpdir)
        fpath = os.path.join(tmpdir, rel_path.name)
        # Guard against filename collisions across lines
        if os.path.exists(fpath):
            stem, suffix = rel_path.stem, rel_path.suffix
            fpath = os.path.join(tmpdir, f"{line_name}_{stem}{suffix}")
        with open(fpath, "wb") as fh:
            fh.write(raw)

        line_map.setdefault(line_name, []).append(fpath)

    all_paths = [p for plist in line_map.values() for p in plist]
    return all_paths, tmpdir, line_map


# ──────────────────────────────────────────────────────────────────────────────
# Plotly pseudosection heatmap builders
# ──────────────────────────────────────────────────────────────────────────────


def _empty_plotly_heatmap(
    dark: bool = True, msg: str = "Load data to view pseudosection"
) -> go.Figure:
    """Empty Plotly figure shown before data is loaded — no grid lines."""
    bg = "#1e1e2e" if dark else "#e6e9ef"
    fig = go.Figure()
    fig.add_annotation(
        text=msg,
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
        showarrow=False,
        font=dict(size=13, color="#585b70"),
    )
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor=bg,
        plot_bgcolor=bg,
        height=340,
        xaxis=dict(visible=False, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, showgrid=False, zeroline=False),
    )
    return fig


def build_plotly_pseudosection(
    sites,
    quantity: str = "rho_xy",
    dark: bool = True,
) -> go.Figure:
    """
    Build an interactive Plotly heatmap for the apparent-resistivity or
    phase pseudosection.

    Reads data directly from Sites.as_list() -> EDIFile.Z to avoid the
    broken _df_resphase / _iter_items path in pycsamt v2.

    Parameters
    ----------
    sites : Sites
        Survey data object.
    quantity : str
        One of ``"rho_xy"``, ``"rho_yx"``, ``"phi_xy"``, ``"phi_yx"``.
    dark : bool
        Apply dark Catppuccin Mocha palette to the layout.

    Returns
    -------
    go.Figure
    """
    import numpy as np
    import pandas as pd

    from pycsamt.emtools._core import ensure_sites

    bg_color = "#1e1e2e" if dark else "#e6e9ef"
    paper_bg = "#181825" if dark else "#eff1f5"
    text_color = "#cdd6f4" if dark else "#4c4f69"
    grid_color = "#313244" if dark else "#ccd0da"
    cb_border = "#45475a" if dark else "#bcc0cc"

    # Map quantity name to Z tensor indices (row, col)
    _COMP = {
        "rho_xy": (0, 1),
        "phi_xy": (0, 1),
        "rho_yx": (1, 0),
        "phi_yx": (1, 0),
        "rho_xx": (0, 0),
        "phi_xx": (0, 0),
        "rho_yy": (1, 1),
        "phi_yy": (1, 1),
    }
    comp = _COMP.get(quantity, (0, 1))
    is_rho = quantity.startswith("rho")

    try:
        S = ensure_sites(sites, verbose=0)
        edis = S.as_list()
        if not edis:
            return _empty_plotly_heatmap(dark, "No stations found in survey")

        rows = []
        for edi in edis:
            Z = getattr(edi, "Z", None)
            if Z is None:
                continue
            freq = getattr(Z, "freq", None)
            if freq is None or len(freq) == 0:
                continue
            arr = Z.resistivity if is_rho else Z.phase
            if arr is None or arr.shape[0] != len(freq):
                continue
            vals = arr[:, comp[0], comp[1]]
            station = getattr(edi, "station", "??")
            period = np.where(freq > 0, 1.0 / freq, np.nan)
            for p, v in zip(period, vals):
                if np.isfinite(p) and np.isfinite(v):
                    if is_rho and v <= 0:
                        continue
                    rows.append(
                        {
                            "period": float(p),
                            "station": station,
                            "value": float(v),
                        }
                    )

        if not rows:
            return _empty_plotly_heatmap(dark, f"No valid data for '{quantity}'")

        df = pd.DataFrame(rows)
        piv = df.pivot_table(
            index="period",
            columns="station",
            values="value",
            aggfunc="median",
        ).sort_index()

        if piv.empty:
            return _empty_plotly_heatmap(dark, "Pivot table empty")

    except Exception as exc:
        return _empty_plotly_heatmap(dark, f"Data error: {exc}")

    Y_periods = piv.index.to_numpy(dtype=float)
    Z_mat = piv.to_numpy(dtype=float)
    station_names = list(piv.columns)

    if is_rho:
        with np.errstate(divide="ignore", invalid="ignore"):
            Z_plot = np.where(Z_mat > 0, np.log10(Z_mat), np.nan)
        colorscale = "Jet"
        colorbar_title = "log10(rho) Ohm.m"
        hover_label = "log10(rho)"
        zmid = None
    else:
        Z_plot = Z_mat
        colorscale = "RdBu_r"
        colorbar_title = "Phase (deg)"
        hover_label = "Phase"
        zmid = 45.0

    Y_log = np.log10(np.where(Y_periods > 0, Y_periods, np.nan))

    comp_label = quantity.split("_", 1)[1].upper() if "_" in quantity else "XY"
    ps_title = f"App. Resistivity ({comp_label})" if is_rho else f"Phase ({comp_label})"

    heatmap_kwargs: dict = dict(
        z=Z_plot,
        x=station_names,
        y=Y_log.tolist(),
        colorscale=colorscale,
        zsmooth="best",
        colorbar=dict(
            title=dict(
                text=colorbar_title,
                side="right",
                font=dict(color=text_color, size=11),
            ),
            tickfont=dict(color=text_color, size=9),
            bgcolor=bg_color,
            bordercolor=cb_border,
            thickness=14,
        ),
        hovertemplate=(
            "<b>%{x}</b><br>"
            "log10(T): %{y:.2f} s<br>" + hover_label + ": %{z:.2f}" + "<extra></extra>"
        ),
    )
    if zmid is not None:
        heatmap_kwargs["zmid"] = zmid

    fig = go.Figure(go.Heatmap(**heatmap_kwargs))

    fig.update_layout(
        title=dict(
            text=f"{ps_title} — Pseudosection",
            font=dict(size=12, color=text_color),
            x=0.0,
            pad=dict(l=6),
        ),
        margin=dict(l=10, r=10, t=36, b=10),
        paper_bgcolor=paper_bg,
        plot_bgcolor=bg_color,
        font=dict(color=text_color),
        xaxis=dict(
            title="Station",
            side="top",
            tickfont=dict(size=8, color=text_color),
            tickangle=-90,
            showgrid=False,
            linecolor=cb_border,
        ),
        yaxis=dict(
            title="log10(T) (s)",
            autorange="reversed",
            tickfont=dict(size=9, color=text_color),
            gridcolor=grid_color,
            gridwidth=0.5,
        ),
        height=340,
    )

    return fig


def build_multi_pseudosection(
    sites,
    components: list,
    plot_type: str = "rho",
    dark: bool = True,
) -> go.Figure:
    """
    Multi-component pseudosection: one stacked heatmap row per component.

    Parameters
    ----------
    sites : Sites
    components : list of str
        Subset of ["xy", "yx", "xx", "yy"].
    plot_type : "rho" or "phi"
    dark : bool
    """
    import numpy as np
    import pandas as pd
    from plotly.subplots import make_subplots

    if not components:
        components = ["xy"]

    bg_color = "#1e1e2e" if dark else "#e6e9ef"
    paper_bg = "#181825" if dark else "#eff1f5"
    text_color = "#cdd6f4" if dark else "#4c4f69"
    grid_color = "#313244" if dark else "#ccd0da"
    cb_border = "#45475a" if dark else "#bcc0cc"

    _COMP = {"xy": (0, 1), "yx": (1, 0), "xx": (0, 0), "yy": (1, 1)}
    is_rho = plot_type == "rho"

    _COMP_LABELS = {
        "xy": "Z_XY (TE)",
        "yx": "Z_YX (TM)",
        "xx": "Z_XX (Diag)",
        "yy": "Z_YY (Diag)",
    }
    prefix = "ρₐ " if is_rho else "φ "
    titles = [prefix + _COMP_LABELS.get(c, c.upper()) for c in components]

    n = len(components)
    v_space = max(0.04, 0.10 / max(n - 1, 1))

    fig = make_subplots(
        rows=n,
        cols=1,
        shared_xaxes=True,
        subplot_titles=titles,
        vertical_spacing=v_space,
    )

    try:
        from pycsamt.emtools._core import ensure_sites

        S = ensure_sites(sites, verbose=0)
        edis = S.as_list()
    except Exception:
        edis = []

    if not edis:
        return _empty_plotly_heatmap(dark, "No stations found in survey")

    any_data = False
    for row_i, comp in enumerate(components):
        cidx = _COMP.get(comp, (0, 1))
        rows = []
        for edi in edis:
            Z = getattr(edi, "Z", None)
            if Z is None:
                continue
            freq = getattr(Z, "freq", None)
            if freq is None or len(freq) == 0:
                continue
            arr = Z.resistivity if is_rho else Z.phase
            if arr is None or arr.shape[0] != len(freq):
                continue
            vals = arr[:, cidx[0], cidx[1]]
            station = getattr(edi, "station", "??")
            period = np.where(freq > 0, 1.0 / freq, np.nan)
            for p, v in zip(period, vals):
                if np.isfinite(p) and np.isfinite(v):
                    if is_rho and v <= 0:
                        continue
                    rows.append(
                        {
                            "period": float(p),
                            "station": station,
                            "value": float(v),
                        }
                    )

        if not rows:
            continue
        any_data = True

        df = pd.DataFrame(rows)
        piv = df.pivot_table(
            index="period",
            columns="station",
            values="value",
            aggfunc="median",
        ).sort_index()
        if piv.empty:
            continue

        Y_periods = piv.index.to_numpy(dtype=float)
        Z_mat = piv.to_numpy(dtype=float)
        stations = list(piv.columns)

        if is_rho:
            with np.errstate(divide="ignore", invalid="ignore"):
                Z_plot = np.where(Z_mat > 0, np.log10(Z_mat), np.nan)
            colorscale = "Jet"
            cb_title = "log₁₀(Ω·m)"
            hover_lbl = "log₁₀(ρₐ)"
            zmid = None
        else:
            Z_plot = Z_mat
            colorscale = "RdBu_r"
            cb_title = "Phase (°)"
            hover_lbl = "φ"
            zmid = 45.0

        Y_log = np.log10(np.where(Y_periods > 0, Y_periods, np.nan))
        cb_y = 1.0 - (row_i + 0.5) / n
        cb_len = max(0.15, 0.85 / n)
        trace_kw = dict(
            z=Z_plot,
            x=stations,
            y=Y_log.tolist(),
            colorscale=colorscale,
            zsmooth="best",
            showscale=True,
            colorbar=dict(
                title=dict(
                    text=cb_title,
                    side="right",
                    font=dict(color=text_color, size=10),
                ),
                tickfont=dict(color=text_color, size=8),
                bgcolor=bg_color,
                bordercolor=cb_border,
                thickness=12,
                len=cb_len,
                y=cb_y,
                yanchor="middle",
                x=1.01,
                xanchor="left",
            ),
            hovertemplate=(
                "<b>%{x}</b><br>"
                "log₁₀(T): %{y:.2f} s<br>" + hover_lbl + ": %{z:.2f}<extra></extra>"
            ),
        )
        if zmid is not None:
            trace_kw["zmid"] = zmid

        fig.add_trace(go.Heatmap(**trace_kw), row=row_i + 1, col=1)
        fig.update_yaxes(
            autorange="reversed",
            tickfont=dict(size=8, color=text_color),
            gridcolor=grid_color,
            gridwidth=0.5,
            row=row_i + 1,
            col=1,
        )
        fig.update_xaxes(
            showgrid=False,
            **(
                {
                    "side": "top",
                    "tickfont": dict(size=8, color=text_color),
                    "tickangle": -60,
                    "linecolor": cb_border,
                }
                if row_i == 0
                else {}
            ),
            row=row_i + 1,
            col=1,
        )

    if not any_data:
        return _empty_plotly_heatmap(
            dark, f"No valid data for: {', '.join(components)}"
        )

    row_h = 220 if n >= 3 else 270
    total_h = max(280, row_h * n)
    fig.update_layout(
        margin=dict(l=10, r=60, t=50, b=10),
        paper_bgcolor=paper_bg,
        plot_bgcolor=bg_color,
        font=dict(color=text_color, size=10),
        height=total_h,
        showlegend=False,
    )
    return fig
