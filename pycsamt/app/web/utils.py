# -*- coding: utf-8 -*-
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
from typing import Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go


# ──────────────────────────────────────────────────────────────────────────────
# Matplotlib ↔ web helpers
# ──────────────────────────────────────────────────────────────────────────────

#: Catppuccin Mocha palette — matches dark_theme.qss
_DARK_RC: dict = {
    "axes.facecolor":    "#181825",
    "figure.facecolor":  "#1e1e2e",
    "savefig.facecolor": "#1e1e2e",
    "axes.edgecolor":    "#45475a",
    "axes.labelcolor":   "#cdd6f4",
    "text.color":        "#cdd6f4",
    "xtick.color":       "#a6adc8",
    "ytick.color":       "#a6adc8",
    "grid.color":        "#313244",
    "grid.linestyle":    "--",
    "grid.alpha":        0.35,
    "axes.grid":         True,
    "lines.color":       "#89b4fa",
    "legend.facecolor":  "#1e1e2e",
    "legend.edgecolor":  "#45475a",
    "legend.labelcolor": "#cdd6f4",
}

_LIGHT_RC: dict = {
    "axes.facecolor":    "#eff1f5",
    "figure.facecolor":  "#e6e9ef",
    "savefig.facecolor": "#e6e9ef",
    "axes.edgecolor":    "#bcc0cc",
    "axes.labelcolor":   "#4c4f69",
    "text.color":        "#4c4f69",
    "xtick.color":       "#6c6f85",
    "ytick.color":       "#6c6f85",
    "grid.color":        "#ccd0da",
    "grid.linestyle":    "--",
    "grid.alpha":        0.5,
    "axes.grid":         True,
    "lines.color":       "#1e66f5",
}


def apply_web_dark_theme() -> None:
    mpl.rcParams.update(_DARK_RC)


def apply_web_light_theme() -> None:
    mpl.rcParams.update(_LIGHT_RC)


def fig_to_src(fig, dpi: int = 100) -> str:
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
    """Return a 1×1 transparent PNG as a data URI placeholder."""
    fig = plt.figure(figsize=(1, 0.5))
    fig.patch.set_facecolor("#1e1e2e" if dark else "#e6e9ef")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    buf.seek(0)
    enc = base64.b64encode(buf.read()).decode("ascii")
    plt.close(fig)
    return f"data:image/png;base64,{enc}"


# ──────────────────────────────────────────────────────────────────────────────
# Plotly map builder
# ──────────────────────────────────────────────────────────────────────────────

_OVERLAY_LABELS = {
    "Index":               "Station index",
    "Apparent Resistivity": "ρₐ (Ω·m)",
    "Phase":               "Phase (°)",
    "Quality Score":       "Quality",
    "Z-strike":            "Strike (°)",
    "N_freq":              "N frequencies",
}

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
    """Return (center_lat, center_lon, bearing, zoom) for a Scattermapbox layout."""
    import math
    clat = float((lats.min() + lats.max()) / 2)
    clon = float((lons.min() + lons.max()) / 2)
    lat_span = float(lats.max() - lats.min())
    lon_span = float(lons.max() - lons.min())
    max_span = max(lat_span, lon_span, 1e-6)
    zoom = max(2, min(14, int(math.floor(8.0 - math.log2(max_span + 0.001)))))
    return clat, clon, 0, zoom


def build_station_map(
    df: pd.DataFrame,
    selected_id: Optional[str] = None,
    overlay: str = "Index",
    dark: bool = True,
    line_filter: Optional[list] = None,
    marker_size: int = 10,
    basemap_style: Optional[str] = None,
) -> go.Figure:
    """
    Build a Plotly Scattermapbox figure for the station map.

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
    """
    if df is None or df.empty:
        return _empty_map(dark)

    lats = pd.to_numeric(df["Latitude"],  errors="coerce").values
    lons = pd.to_numeric(df["Longitude"], errors="coerce").values
    ids  = df["ID"].values

    valid = ~(np.isnan(lats) | np.isnan(lons))
    if not valid.any():
        return _empty_map(dark)

    lats, lons, ids = lats[valid], lons[valid], ids[valid]
    df_v = df[valid].reset_index(drop=True)

    # Apply line filter
    if line_filter and "Line" in df_v.columns:
        mask = df_v["Line"].isin(line_filter)
        df_v = df_v[mask].reset_index(drop=True)
        lats = pd.to_numeric(df_v["Latitude"],  errors="coerce").values
        lons = pd.to_numeric(df_v["Longitude"], errors="coerce").values
        ids  = df_v["ID"].values
        if len(ids) == 0:
            return _empty_map(dark)

    bg_color = "#1e1e2e" if dark else "#eff1f5"
    text_col = "#cdd6f4" if dark else "#4c4f69"
    sel_col  = "#f9e2af"   # Catppuccin yellow — selected station

    # Overlay colour values
    if overlay in df_v.columns:
        colour_vals = pd.to_numeric(df_v[overlay], errors="coerce").fillna(0).values
    else:
        colour_vals = np.arange(len(df_v), dtype=float)

    sel_size = int(marker_size * 1.8)
    sizes    = np.where(ids == selected_id, sel_size, marker_size).tolist() if selected_id else [marker_size] * len(ids)
    borders  = [sel_col if sid == selected_id else "#ffffff" for sid in ids]
    bwidths  = [2.5   if sid == selected_id else 0.8    for sid in ids]

    hover_label = _OVERLAY_LABELS.get(overlay, overlay)

    fig = go.Figure()

    # ── Survey line polylines (one trace per line) ──────────────────────────
    if "Line" in df_v.columns:
        unique_lines = [ln for ln in df_v["Line"].unique() if ln and ln != "—"]
        for idx, line_name in enumerate(unique_lines):
            grp = df_v[df_v["Line"] == line_name].copy()
            # Sort along the dominant axis so the polyline makes geographic sense
            lon_span = grp["Longitude"].max() - grp["Longitude"].min()
            lat_span = grp["Latitude"].max()  - grp["Latitude"].min()
            sort_col = "Longitude" if lon_span >= lat_span else "Latitude"
            grp = grp.sort_values(sort_col)
            lcol = _LINE_COLORS[idx % len(_LINE_COLORS)]
            fig.add_trace(go.Scattermapbox(
                lat=grp["Latitude"].tolist(),
                lon=grp["Longitude"].tolist(),
                mode="lines",
                line=dict(color=lcol, width=2),
                name=str(line_name),
                hoverinfo="skip",
                showlegend=len(unique_lines) > 1,
            ))

    # ── Station markers ─────────────────────────────────────────────────────
    fig.add_trace(go.Scattermapbox(
        lat=lats.tolist(),
        lon=lons.tolist(),
        text=ids.tolist(),
        customdata=ids.tolist(),
        mode="markers+text",
        textposition="top right",
        textfont=dict(size=9, color=text_col),
        marker=dict(
            size=sizes,
            color=colour_vals.tolist(),
            colorscale="Plasma",
            showscale=True,
            colorbar=dict(
                title=dict(text=hover_label, font=dict(color=text_col, size=11)),
                tickfont=dict(color=text_col, size=9),
                bgcolor=bg_color,
                bordercolor="#45475a" if dark else "#bcc0cc",
                thickness=12,
                len=0.7,
                x=1.0,
            ),
            opacity=0.92,
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
    ))

    # Selected station highlight ring
    if selected_id and selected_id in ids:
        sel_mask = ids == selected_id
        fig.add_trace(go.Scattermapbox(
            lat=lats[sel_mask].tolist(),
            lon=lons[sel_mask].tolist(),
            mode="markers",
            marker=dict(size=24, color="rgba(0,0,0,0)",
                        opacity=1,
                        # Draw a ring by overlaying a larger transparent marker
                        # with an outline colour set via the line field
                        ),
            hoverinfo="skip",
            showlegend=False,
        ))

    # ── Map layout ──────────────────────────────────────────────────────────
    clat, clon, _, zoom = _auto_zoom(lats, lons)
    tile_style = basemap_style or ("carto-darkmatter" if dark else "carto-positron")

    fig.update_layout(
        mapbox=dict(
            style=tile_style,
            center=dict(lat=clat, lon=clon),
            zoom=zoom,
        ),
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor=bg_color,
        font=dict(color=text_col, size=11),
        showlegend=True,
        legend=dict(
            bgcolor="rgba(30,30,46,.75)" if dark else "rgba(239,241,245,.80)",
            bordercolor="#45475a" if dark else "#bcc0cc",
            borderwidth=1,
            font=dict(size=10, color=text_col),
            x=0.01, y=0.99,
            xanchor="left", yanchor="top",
        ),
        clickmode="event+select",
        uirevision="map",
    )

    return fig


def _empty_map(dark: bool = True) -> go.Figure:
    bg       = "#1e1e2e" if dark else "#eff1f5"
    text_col = "#585b70" if dark else "#8c8fa1"
    tile     = "carto-darkmatter" if dark else "carto-positron"

    fig = go.Figure()
    fig.add_trace(go.Scattermapbox(
        lat=[], lon=[], mode="markers", showlegend=False,
    ))
    fig.add_annotation(
        text="Load survey data to display the station map",
        xref="paper", yref="paper",
        x=0.5, y=0.5, showarrow=False,
        font=dict(size=13, color=text_col),
        bgcolor="rgba(30,30,46,.65)" if dark else "rgba(239,241,245,.70)",
        borderpad=8,
    )
    fig.update_layout(
        mapbox=dict(
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
    btn_label: Optional[str] = None,
    btn_id: Optional[str] = None,
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
    from dash import html
    import dash_bootstrap_components as dbc

    children = [
        html.I(className=f"bi bi-{icon} empty-state-icon"),
        html.P(title,    className="empty-state-title"),
        html.P(subtitle, className="empty-state-sub"),
    ]
    if btn_label and btn_id:
        children.append(
            dbc.Button(
                [html.I(className="bi bi-cloud-upload me-2"), btn_label],
                id=btn_id, color="primary", size="sm",
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
        for survey_path in sorted(fp for fp in p.rglob("*") if fp.is_file() and _matches(fp)):
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
        contents  = [contents]
        filenames = [filenames]

    tmpdir = tempfile.mkdtemp(prefix="pycsamt_upload_")
    paths: list[str] = []

    for content, name in zip(contents, filenames):
        if Path(name).suffix.lower() not in _UPLOAD_EXTENSIONS:
            continue  # skip unsupported files silently
        # Strip the data-URI header: "data:<mime>;base64,<payload>"
        try:
            _, b64 = content.split(",", 1)
            raw    = base64.b64decode(b64)
        except Exception:
            continue
        fpath = os.path.join(tmpdir, name)
        with open(fpath, "wb") as fh:
            fh.write(raw)
        paths.append(fpath)

    return paths, tmpdir


# ──────────────────────────────────────────────────────────────────────────────
# Plotly pseudosection heatmap builders
# ──────────────────────────────────────────────────────────────────────────────

def _empty_plotly_heatmap(dark: bool = True, msg: str = "Load data to view pseudosection") -> go.Figure:
    """Empty Plotly figure shown before data is loaded — no grid lines."""
    bg = "#1e1e2e" if dark else "#e6e9ef"
    fig = go.Figure()
    fig.add_annotation(
        text=msg,
        xref="paper", yref="paper",
        x=0.5, y=0.5, showarrow=False,
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

    Reuses the same pivot-table pipeline as ``et.pseudosection()`` so the
    data is identical to the matplotlib version — only the renderer changes.

    Parameters
    ----------
    sites : Sites
        Survey data object.
    quantity : str
        ``"rho_xy"`` or ``"phi_xy"`` (any column produced by
        ``emtools.inspect._df_resphase``).
    dark : bool
        Apply dark Catppuccin Mocha palette to the layout.

    Returns
    -------
    go.Figure
    """
    import numpy as np
    from pycsamt.emtools._core import ensure_sites
    from pycsamt.emtools.inspect import _df_resphase

    bg_color   = "#1e1e2e" if dark else "#e6e9ef"
    paper_bg   = "#181825" if dark else "#eff1f5"
    text_color = "#cdd6f4" if dark else "#4c4f69"
    grid_color = "#313244" if dark else "#ccd0da"
    cb_border  = "#45475a" if dark else "#bcc0cc"

    try:
        S  = ensure_sites(sites, verbose=0)
        df = _df_resphase(S, kind="resphase")
    except Exception as exc:
        return _empty_plotly_heatmap(dark, f"Data error: {exc}")

    if df is None or df.empty or quantity not in df.columns:
        return _empty_plotly_heatmap(dark, f"No data for '{quantity}'")

    df = df.copy()
    if "period" not in df.columns:
        df["period"] = 1.0 / df["freq"].replace(0, np.nan)

    piv = (
        df.pivot_table(
            index="period",
            columns="station",
            values=quantity,
            aggfunc="median",
        )
        .sort_index()
    )
    if piv.empty:
        return _empty_plotly_heatmap(dark, "Pivot table empty — check quantity name")

    Y_periods     = piv.index.to_numpy(dtype=float)
    Z             = piv.to_numpy(dtype=float)           # [n_period, n_station]
    station_names = list(piv.columns)

    is_rho = quantity.startswith("rho")

    if is_rho:
        with np.errstate(divide="ignore", invalid="ignore"):
            Z_plot = np.where(Z > 0, np.log10(Z), np.nan)
        colorscale     = "Jet"
        colorbar_title = "log₁₀(ρₐ) (Ω·m)"
        hover_label    = "log₁₀(ρₐ)"
        zmid           = None
    else:
        Z_plot         = Z
        colorscale     = "RdBu_r"
        colorbar_title = "φ (°)"
        hover_label    = "φ"
        zmid           = 45.0          # centre diverging scale at 45 °

    Y_log = np.log10(np.where(Y_periods > 0, Y_periods, np.nan))

    heatmap_kwargs: dict = dict(
        z=Z_plot,
        x=station_names,
        y=Y_log,
        colorscale=colorscale,
        zsmooth="best",
        colorbar=dict(
            title=dict(text=colorbar_title, font=dict(color=text_color, size=11)),
            tickfont=dict(color=text_color, size=9),
            bgcolor=bg_color,
            bordercolor=cb_border,
            thickness=14,
        ),
        hovertemplate=(
            "<b>%{x}</b><br>"
            "log₁₀(T): %{y:.2f} s<br>"
            f"{hover_label}: %{{z:.2f}}"
            "<extra></extra>"
        ),
    )
    if zmid is not None:
        heatmap_kwargs["zmid"] = zmid

    fig = go.Figure(go.Heatmap(**heatmap_kwargs))

    ps_title = "ρₐ (XY)" if is_rho else "φ (XY)"
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
            side="top",                 # MT convention: stations at the top
            tickfont=dict(size=8, color=text_color),
            tickangle=-90,
            showgrid=False,
            linecolor=cb_border,
        ),
        yaxis=dict(
            title="log₁₀(T) (s)",
            autorange="reversed",       # short period (shallow) at the top
            tickfont=dict(size=9, color=text_color),
            gridcolor=grid_color,
            gridwidth=0.5,
        ),
        height=340,
    )

    return fig
