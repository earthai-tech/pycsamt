# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/tools.py — Tools off-canvas panel.

Responsibilities
----------------
- Open the Tools offcanvas when any item in the navbar Tools dropdown is clicked
  OR when a tool-picker rail button inside the offcanvas is clicked.
- Track the active tool in ACTIVE_TOOL store.
- Render the correct tool body panel inside the offcanvas.
- Execute each tool and stream results back into the panel.
"""
from __future__ import annotations

import json
import io

import dash_bootstrap_components as dbc
from dash import (
    ALL, Input, Output, State,
    callback_context as ctx,
    no_update, html, dcc, dash_table,
)

from pycsamt.app.web.layout import IDs, _TOOL_REGISTRY, _TOOL_BY_ID, _tool_panel, _icon


# ── Shared helpers ────────────────────────────────────────────────────────────

def _no_data_msg() -> html.Div:
    return html.Div(
        [html.I(className="bi bi-exclamation-triangle me-2 text-warning"),
         "Load survey data first."],
        className="tool-no-data",
    )


def _run_btn(btn_id: str, label: str, icon: str = "play-fill") -> html.Button:
    return html.Button(
        [html.I(className=f"bi bi-{icon} me-1"), label],
        id=btn_id,
        className="btn btn-primary btn-sm mt-2",
        n_clicks=0,
    )


def _out_area(div_id: str, min_h: int = 160, cls: str = "mt-3 fig-area") -> html.Div:
    return html.Div(id=div_id, className=cls, style={"minHeight": f"{min_h}px"})


# ── Tool body builders ────────────────────────────────────────────────────────

def _build_tool_body(tool_id: str, store_data: dict | None) -> html.Div:
    n  = (store_data or {}).get("n_stations", 0)
    nd = _no_data_msg()

    builders = {
        "strike":        lambda: _strike_body(n, nd),
        "validator":     lambda: _validator_body(n, nd),
        "converter":     lambda: _converter_body(n, nd),
        "coords":        lambda: _coords_body(),
        "dimensionality":lambda: _dimensionality_body(n, nd),
        "freq-editor":   lambda: _freq_editor_body(n, nd),
        "batch":         lambda: _batch_body(n, nd),
        "elevation":     lambda: _elevation_body(n, nd),
        "response":      lambda: _response_body(n, nd, store_data),
        "strike-profile":lambda: _strike_profile_body(n, nd, store_data),
        "phase-tensor":  lambda: _phase_tensor_body(n, nd),
        "layered-model": lambda: _layered_model_body(),
    }
    fn = builders.get(tool_id)
    return fn() if fn else _placeholder_body(tool_id)


# ── Implemented tools ─────────────────────────────────────────────────────────

def _strike_body(n: int, no_data: html.Div) -> html.Div:
    if not n:
        return no_data
    return html.Div([
        html.Label("Strike method", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-strike-method",
            options=[
                {"label": "Bahr (1991)",          "value": "bahr"},
                {"label": "Weaver (2000)",         "value": "weaver"},
                {"label": "Strike from Z",         "value": "z"},
            ],
            value="bahr", clearable=False, className="mb-2",
        ),
        html.Label("Frequency subset", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-strike-freq",
            options=[
                {"label": "All frequencies",          "value": "all"},
                {"label": "≥ 100 Hz (near-surface)",  "value": "high"},
                {"label": "≤ 1 Hz (deep crust)",      "value": "low"},
            ],
            value="all", clearable=False, className="mb-2",
        ),
        _run_btn("tool-strike-run", "Run Strike Analysis"),
        _out_area("tool-strike-out"),
    ])


def _validator_body(n: int, no_data: html.Div) -> html.Div:
    if not n:
        return no_data
    return html.Div([
        html.P("Run a per-station quality checklist on the loaded survey.",
               className="text-muted small"),
        _run_btn("tool-valid-run", "Run Validation", "clipboard-check"),
        _out_area("tool-valid-out", min_h=140, cls="mt-3 log-area"),
    ])


def _converter_body(n: int, no_data: html.Div) -> html.Div:
    if not n:
        return no_data
    return html.Div([
        html.Label("Export format", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-conv-fmt",
            options=[{"label": f, "value": f.lower()} for f in ["EDI", "CSV", "JSON", "HDF5"]],
            value="csv", clearable=False, className="mb-2",
        ),
        html.Label("Output path (server-side)", className="ctrl-label"),
        dbc.Input(id="tool-conv-path", placeholder="/data/exports/",
                  size="sm", className="mb-2"),
        _run_btn("tool-conv-run", "Convert & Save", "download"),
        _out_area("tool-conv-out", min_h=60, cls="mt-2 small"),
    ])


def _coords_body() -> html.Div:
    return html.Div([
        html.P("Convert coordinates between UTM and geographic Lat / Lon.",
               className="text-muted small"),
        html.Label("Input format", className="ctrl-label"),
        dbc.Input(id="tool-coord-easting",  type="number",
                  placeholder="Easting (m)", size="sm", className="mb-1"),
        dbc.Input(id="tool-coord-northing", type="number",
                  placeholder="Northing (m)", size="sm", className="mb-1"),
        dbc.Input(id="tool-coord-zone",     type="text",
                  placeholder="UTM zone  e.g. 49N", size="sm", className="mb-2"),
        _run_btn("tool-coord-run", "Convert", "arrow-left-right"),
        _out_area("tool-coord-out", min_h=60, cls="mt-2 small"),
    ])


def _dimensionality_body(n: int, no_data: html.Div) -> html.Div:
    if not n:
        return no_data
    return html.Div([
        html.Label("Method", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-dim-method",
            options=[
                {"label": "Bahr phase-sensitive skew", "value": "bahr"},
                {"label": "Swift skew",                "value": "swift"},
                {"label": "WAL (Weaver et al.)",       "value": "wal"},
            ],
            value="bahr", clearable=False, className="mb-2",
        ),
        _run_btn("tool-dim-run", "Classify", "grid"),
        _out_area("tool-dim-out", min_h=200),
    ])


def _freq_editor_body(n: int, no_data: html.Div) -> html.Div:
    if not n:
        return no_data
    return html.Div([
        html.P("Mask or drop frequency bands failing the confidence threshold.",
               className="text-muted small"),
        html.Label("Min confidence score (0 – 1)", className="ctrl-label"),
        dcc.Slider(id="tool-freq-thresh", min=0, max=1, step=0.05, value=0.50,
                   marks={0: "0", 0.5: "0.5", 1: "1"}, className="mb-3"),
        _run_btn("tool-freq-run", "Apply QC", "scissors"),
        _out_area("tool-freq-out", min_h=60, cls="mt-2 small"),
    ])


# ── New tool bodies ───────────────────────────────────────────────────────────

def _batch_body(n: int, no_data: html.Div) -> html.Div:
    if not n:
        return no_data
    return html.Div([
        html.P("Save all profile and pseudosection figures to disk in one operation.",
               className="text-muted small"),
        html.Label("Output folder (server-side)", className="ctrl-label"),
        dbc.Input(id="tool-batch-path", placeholder="/data/figures/", size="sm", className="mb-2"),
        dbc.Row([
            dbc.Col([
                html.Label("Format", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-batch-fmt",
                    options=[{"label": f, "value": f.lower()}
                             for f in ["PNG", "PDF", "SVG", "EPS"]],
                    value="png", clearable=False,
                ),
            ], width=6),
            dbc.Col([
                html.Label("DPI", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-batch-dpi",
                    options=[{"label": str(d), "value": d}
                             for d in [72, 96, 150, 300]],
                    value=150, clearable=False,
                ),
            ], width=6),
        ], className="mb-2"),
        html.Label("Include", className="ctrl-label"),
        dbc.Checklist(
            id="tool-batch-items",
            options=[
                {"label": "MT curves (ρ, φ)",     "value": "curves"},
                {"label": "Pseudosections",        "value": "pseudo"},
                {"label": "Station map",           "value": "map"},
                {"label": "Strike / dim. plots",   "value": "strike"},
            ],
            value=["curves", "pseudo"],
            inline=False,
            className="mb-2",
        ),
        _run_btn("tool-batch-run", "Export All Figures", "images"),
        _out_area("tool-batch-out", min_h=60, cls="mt-2 small"),
    ])


def _elevation_body(n: int, no_data: html.Div) -> html.Div:
    if not n:
        return no_data
    return html.Div([
        html.P(
            "Query SRTM elevation for every loaded station and update the "
            "Elevation column.  Requires an internet connection.",
            className="text-muted small",
        ),
        html.Label("DEM source", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-elev-source",
            options=[
                {"label": "Open-Elevation (free)",  "value": "open-elevation"},
                {"label": "SRTM 30 m (NASADEM)",    "value": "nasadem"},
                {"label": "SRTM 90 m (offline CSV)","value": "srtm90"},
            ],
            value="open-elevation", clearable=False, className="mb-2",
        ),
        _run_btn("tool-elev-run", "Fetch Elevations", "cloud-download"),
        _out_area("tool-elev-out", min_h=80, cls="mt-2 small"),
    ])


def _response_body(n: int, no_data: html.Div, store_data: dict | None) -> html.Div:
    if not n:
        return no_data
    records  = (store_data or {}).get("station_records", [])
    sta_opts = [{"label": r["ID"], "value": r["ID"]} for r in records]
    first    = sta_opts[0]["value"] if sta_opts else None
    return html.Div([
        html.Label("Station", className="ctrl-label"),
        dcc.Dropdown(id="tool-resp-station", options=sta_opts, value=first,
                     clearable=False, className="mb-2"),
        html.Label("Components", className="ctrl-label"),
        dbc.Checklist(
            id="tool-resp-comps",
            options=[
                {"label": "Zxy (TE mode)", "value": "xy"},
                {"label": "Zyx (TM mode)", "value": "yx"},
                {"label": "Zxx",           "value": "xx"},
                {"label": "Zyy",           "value": "yy"},
            ],
            value=["xy", "yx"],
            inline=True,
            className="mb-2",
        ),
        dbc.Row([
            dbc.Col([
                html.Label("Period axis", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-resp-xscale",
                    options=[{"label": "Log", "value": "log"},
                             {"label": "Linear", "value": "linear"}],
                    value="log", clearable=False,
                ),
            ], width=6),
            dbc.Col([
                html.Label("Plot", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-resp-ytype",
                    options=[{"label": "ρ & φ", "value": "rho_phi"},
                             {"label": "Re/Im Z", "value": "reim"}],
                    value="rho_phi", clearable=False,
                ),
            ], width=6),
        ], className="mb-2"),
        _run_btn("tool-resp-run", "Plot Response"),
        _out_area("tool-resp-out", min_h=320),
    ])


def _strike_profile_body(n: int, no_data: html.Div, store_data: dict | None) -> html.Div:
    if not n:
        return no_data
    line_counts = (store_data or {}).get("line_counts", {})
    line_opts   = [{"label": f"{ln} ({cnt} sta)", "value": ln}
                   for ln, cnt in line_counts.items()]
    if not line_opts:
        line_opts = [{"label": "All stations", "value": "__all__"}]
    return html.Div([
        html.Label("Survey line", className="ctrl-label"),
        dcc.Dropdown(id="tool-sprof-line", options=line_opts,
                     value=line_opts[0]["value"], clearable=False, className="mb-2"),
        html.Label("Strike method", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-sprof-method",
            options=[
                {"label": "Bahr (1991)",   "value": "bahr"},
                {"label": "Weaver (2000)", "value": "weaver"},
                {"label": "From Z",        "value": "z"},
            ],
            value="bahr", clearable=False, className="mb-2",
        ),
        html.Label("Period band", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-sprof-band",
            options=[
                {"label": "All periods",           "value": "all"},
                {"label": "< 0.01 s (shallow)",    "value": "shallow"},
                {"label": "0.01–1 s (mid-crust)",  "value": "mid"},
                {"label": "> 1 s (deep)",          "value": "deep"},
            ],
            value="all", clearable=False, className="mb-2",
        ),
        dbc.Switch(id="tool-sprof-iqr", label="Show IQR ribbon", value=True,
                   className="mb-2"),
        _run_btn("tool-sprof-run", "Plot Strike Profile"),
        _out_area("tool-sprof-out", min_h=280),
    ])


def _phase_tensor_body(n: int, no_data: html.Div) -> html.Div:
    if not n:
        return no_data
    return html.Div([
        dbc.Row([
            dbc.Col([
                html.Label("Period (s)", className="ctrl-label"),
                dbc.Input(id="tool-pt-period", type="number", value=1.0,
                          min=1e-5, step=0.1, size="sm"),
            ], width=6),
            dbc.Col([
                html.Label("Color by", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-pt-colorby",
                    options=[
                        {"label": "Φ_min",          "value": "phi_min"},
                        {"label": "Φ_max",          "value": "phi_max"},
                        {"label": "Skew β",         "value": "beta"},
                        {"label": "Azimuth α",      "value": "alpha"},
                        {"label": "Ellipticity λ",  "value": "ellipticity"},
                    ],
                    value="phi_min", clearable=False,
                ),
            ], width=6),
        ], className="mb-2"),
        html.Label("Colour scale", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-pt-cmap",
            options=[{"label": c, "value": c}
                     for c in ["RdYlBu_r", "viridis", "plasma", "coolwarm", "jet"]],
            value="RdYlBu_r", clearable=False, className="mb-2",
        ),
        dbc.Switch(id="tool-pt-tipper", label="Overlay tipper arrows", value=False,
                   className="mb-2"),
        _run_btn("tool-pt-run", "Plot Phase Tensor Map"),
        _out_area("tool-pt-out", min_h=360),
    ])


def _layered_model_body() -> html.Div:
    default_layers = [
        {"#": 1, "ρ (Ω·m)": 100.0,  "Thickness (m)": 50.0},
        {"#": 2, "ρ (Ω·m)": 500.0,  "Thickness (m)": 200.0},
        {"#": 3, "ρ (Ω·m)": 50.0,   "Thickness (m)": 500.0},
        {"#": 4, "ρ (Ω·m)": 1000.0, "Thickness (m)": None},  # halfspace
    ]
    return html.Div([
        html.Label("Preset", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-lm-preset",
            options=[
                {"label": "Custom",         "value": "custom"},
                {"label": "Simple 3-layer", "value": "simple"},
                {"label": "Marine",         "value": "marine"},
                {"label": "Continental",    "value": "continental"},
            ],
            value="simple", clearable=False, className="mb-2",
        ),
        html.Label("Layers  (last row = halfspace, leave Thickness blank)",
                   className="ctrl-label"),
        dash_table.DataTable(
            id="tool-lm-table",
            columns=[
                {"name": "#",              "id": "#",              "editable": False},
                {"name": "ρ (Ω·m)",       "id": "ρ (Ω·m)",       "editable": True,
                 "type": "numeric"},
                {"name": "Thickness (m)", "id": "Thickness (m)", "editable": True,
                 "type": "numeric"},
            ],
            data=default_layers,
            editable=True,
            row_deletable=True,
            style_table={"overflowX": "auto", "fontSize": "0.82rem"},
            style_cell={"padding": "4px 8px", "textAlign": "right",
                        "backgroundColor": "var(--mantle)",
                        "color": "var(--text)", "border": "1px solid var(--surface0)"},
            style_header={"backgroundColor": "var(--surface0)",
                          "color": "var(--subtext1)", "fontWeight": "600"},
            style_data_conditional=[
                {"if": {"state": "active"}, "backgroundColor": "var(--surface1)",
                 "border": "1px solid var(--blue)"},
            ],
            className="mb-1",
        ),
        html.Button(
            [html.I(className="bi bi-plus-lg me-1"), "Add layer"],
            id="tool-lm-add",
            className="btn btn-sm btn-outline-secondary me-2",
            n_clicks=0,
        ),
        html.Hr(className="tool-sep"),
        _run_btn("tool-lm-preview", "Preview Model", "eye"),
        _out_area("tool-lm-out", min_h=260),
    ])


def _placeholder_body(tool_id: str) -> html.Div:
    _, _, label, _ = _TOOL_BY_ID.get(tool_id, (tool_id, "tools", tool_id, ""))
    return html.Div([
        html.Div([
            html.I(className="bi bi-wrench-adjustable",
                   style={"fontSize": "2rem", "color": "var(--blue)", "marginBottom": "8px"}),
            html.P(f"{label} is available in the desktop application.",
                   className="text-muted small mb-1"),
            html.P("Web integration coming in a future release.",
                   className="text-muted small"),
        ], className="text-center py-4"),
    ], className="tool-placeholder")


# ── Callbacks ─────────────────────────────────────────────────────────────────

def register_tools(app) -> None:

    # 1. Open offcanvas when a navbar dropdown item or picker rail button is clicked
    @app.callback(
        Output(IDs.TOOLS_CANVAS, "is_open"),
        Output(IDs.ACTIVE_TOOL,  "data"),
        Input({"type": "tool-menu-item",  "index": ALL}, "n_clicks"),
        Input({"type": "tool-picker-btn", "index": ALL}, "n_clicks"),
        State(IDs.TOOLS_CANVAS, "is_open"),
        State(IDs.ACTIVE_TOOL,  "data"),
        prevent_initial_call=True,
    )
    def _open_tools(menu_clicks, picker_clicks, is_open, current_tool):
        triggered = ctx.triggered
        if not triggered:
            return no_update, no_update
        if not any(t["value"] for t in triggered if t["value"]):
            return no_update, no_update
        trigger_id = ctx.triggered_id
        if trigger_id is None:
            return no_update, no_update
        new_tool = (trigger_id.get("index", current_tool)
                    if isinstance(trigger_id, dict) else current_tool)
        return True, new_tool

    # 2. Render tool panel (header + body) when active tool or loaded data changes
    @app.callback(
        Output("tools-panel-content", "children"),
        Input(IDs.ACTIVE_TOOL, "data"),
        Input(IDs.STORE_DATA,  "data"),
        prevent_initial_call=False,
    )
    def _render_tool(tool_id, store_data):
        tool_id = tool_id or "strike"
        body = _build_tool_body(tool_id, store_data)
        return _tool_panel(tool_id, body=body)

    # 3. Highlight active picker button
    @app.callback(
        Output({"type": "tool-picker-btn", "index": ALL}, "className"),
        Input(IDs.ACTIVE_TOOL, "data"),
        prevent_initial_call=False,
    )
    def _highlight_picker(active):
        active = active or "strike"
        ids = [t[0] for t in _TOOL_REGISTRY]
        return [
            "tool-picker-btn active" if tid == active else "tool-picker-btn"
            for tid in ids
        ]

    # ── Tool run callbacks ─────────────────────────────────────────────────────

    # Coordinate transformer (no data needed)
    @app.callback(
        Output("tool-coord-out", "children"),
        Input("tool-coord-run", "n_clicks"),
        State("tool-coord-easting",  "value"),
        State("tool-coord-northing", "value"),
        State("tool-coord-zone",     "value"),
        prevent_initial_call=True,
    )
    def _run_coords(n, easting, northing, zone):
        if not n:
            return no_update
        if None in (easting, northing) or not zone:
            return html.Span("⚠ Fill in Easting, Northing, and UTM zone.",
                             className="text-warning")
        try:
            from pycsamt.gis.utils import utm_to_latlon
            lat, lon = utm_to_latlon(easting, northing, zone.strip())
            return html.Div([
                html.Span("Latitude: ", className="text-muted"),
                html.Span(f"{lat:.6f}°", className="fw-semibold me-3"),
                html.Span("Longitude: ", className="text-muted"),
                html.Span(f"{lon:.6f}°", className="fw-semibold"),
            ])
        except Exception as exc:
            return _err(exc)

    # Layered model: load preset into table
    @app.callback(
        Output("tool-lm-table", "data"),
        Input("tool-lm-preset", "value"),
        State("tool-lm-table", "data"),
        prevent_initial_call=True,
    )
    def _lm_preset(preset, current_data):
        presets = {
            "simple": [
                {"#": 1, "ρ (Ω·m)": 100.0,  "Thickness (m)": 50.0},
                {"#": 2, "ρ (Ω·m)": 500.0,  "Thickness (m)": 200.0},
                {"#": 3, "ρ (Ω·m)": 50.0,   "Thickness (m)": 500.0},
                {"#": 4, "ρ (Ω·m)": 1000.0, "Thickness (m)": None},
            ],
            "marine": [
                {"#": 1, "ρ (Ω·m)": 0.3,   "Thickness (m)": 1500.0},
                {"#": 2, "ρ (Ω·m)": 1.0,   "Thickness (m)": 500.0},
                {"#": 3, "ρ (Ω·m)": 10.0,  "Thickness (m)": 3000.0},
                {"#": 4, "ρ (Ω·m)": 100.0, "Thickness (m)": None},
            ],
            "continental": [
                {"#": 1, "ρ (Ω·m)": 200.0,  "Thickness (m)": 2000.0},
                {"#": 2, "ρ (Ω·m)": 50.0,   "Thickness (m)": 10000.0},
                {"#": 3, "ρ (Ω·m)": 500.0,  "Thickness (m)": 20000.0},
                {"#": 4, "ρ (Ω·m)": 1000.0, "Thickness (m)": None},
            ],
        }
        return presets.get(preset, current_data)

    # Layered model: add row
    @app.callback(
        Output("tool-lm-table", "data", allow_duplicate=True),
        Input("tool-lm-add", "n_clicks"),
        State("tool-lm-table", "data"),
        prevent_initial_call=True,
    )
    def _lm_add_row(n, rows):
        if not n:
            return no_update
        rows = list(rows or [])
        rows.append({"#": len(rows) + 1, "ρ (Ω·m)": 100.0, "Thickness (m)": 100.0})
        return rows

    # Layered model: preview resistivity-depth profile
    @app.callback(
        Output("tool-lm-out", "children"),
        Input("tool-lm-preview", "n_clicks"),
        State("tool-lm-table", "data"),
        prevent_initial_call=True,
    )
    def _lm_preview(n, rows):
        if not n or not rows:
            return no_update
        try:
            import plotly.graph_objects as go
            import numpy as np
            rhos, depths = [], [0.0]
            for row in rows:
                rho   = float(row.get("ρ (Ω·m)", 100) or 100)
                thick = row.get("Thickness (m)")
                rhos.append(rho)
                if thick is not None and str(thick).strip() not in ("", "None"):
                    depths.append(depths[-1] + float(thick))
                else:
                    depths.append(depths[-1] + rhos[-1] * 3)  # halfspace visual extent

            # Build step plot
            x_pts, y_pts = [], []
            for i, rho in enumerate(rhos):
                x_pts += [rho, rho]
                y_pts += [depths[i], depths[i + 1]]

            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=x_pts, y=y_pts,
                mode="lines",
                line={"color": "#89b4fa", "width": 2},
                name="ρ profile",
            ))
            fig.update_layout(
                xaxis={"title": "ρ (Ω·m)", "type": "log",
                       "gridcolor": "#313244", "zeroline": False},
                yaxis={"title": "Depth (m)", "autorange": "reversed",
                       "gridcolor": "#313244"},
                plot_bgcolor="#1e1e2e", paper_bgcolor="#1e1e2e",
                font={"color": "#cdd6f4"},
                margin={"l": 60, "r": 20, "t": 20, "b": 50},
                height=240,
            )
            return dcc.Graph(figure=fig, config={"displayModeBar": False})
        except Exception as exc:
            return _err(exc)

    # Strike analysis
    @app.callback(
        Output("tool-strike-out", "children"),
        Input("tool-strike-run", "n_clicks"),
        State("tool-strike-method", "value"),
        State("tool-strike-freq",   "value"),
        State(IDs.STORE_DATA,       "data"),
        State(IDs.SESSION_ID,       "data"),
        prevent_initial_call=True,
    )
    def _run_strike(n, method, freq_band, store_data, session_id):
        if not n:
            return no_update
        try:
            from pycsamt.app.web.cache import cache_get
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")
            from pycsamt.emtools.strike import compute_strike
            result = compute_strike(sites, method=method, freq_band=freq_band)
            import plotly.graph_objects as go
            fig = go.Figure()
            for line, angles in result.items():
                fig.add_trace(go.Box(y=angles, name=line, boxpoints="all",
                                     jitter=0.3, pointpos=-1.8))
            _style_fig(fig, "Strike angle (°)", height=240)
            return dcc.Graph(figure=fig, config={"displayModeBar": False})
        except ImportError:
            return _coming_soon("Strike analysis")
        except Exception as exc:
            return _err(exc)

    # EDI Validator
    @app.callback(
        Output("tool-valid-out", "children"),
        Input("tool-valid-run", "n_clicks"),
        State(IDs.STORE_DATA,  "data"),
        State(IDs.SESSION_ID,  "data"),
        prevent_initial_call=True,
    )
    def _run_validator(n, store_data, session_id):
        if not n:
            return no_update
        try:
            from pycsamt.app.web.cache import cache_get
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")
            from pycsamt.emtools.qc import validate_sites
            report = validate_sites(sites)
            lines = []
            for name, issues in report.items():
                tag = "✓" if not issues else f"✗ {len(issues)} issue(s)"
                lines.append(f"{name}: {tag}")
                for iss in issues:
                    lines.append(f"    · {iss}")
            return html.Pre("\n".join(lines), style={"fontSize": "0.78rem",
                                                      "whiteSpace": "pre-wrap"})
        except ImportError:
            return _coming_soon("EDI Validator")
        except Exception as exc:
            return _err(exc)

    # Format converter
    @app.callback(
        Output("tool-conv-out", "children"),
        Input("tool-conv-run", "n_clicks"),
        State("tool-conv-fmt",  "value"),
        State("tool-conv-path", "value"),
        State(IDs.STORE_DATA,   "data"),
        State(IDs.SESSION_ID,   "data"),
        prevent_initial_call=True,
    )
    def _run_converter(n, fmt, out_path, store_data, session_id):
        if not n:
            return no_update
        if not out_path or not out_path.strip():
            return _warn("Enter an output path first.")
        try:
            from pathlib import Path
            from pycsamt.app.web.cache import cache_get
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")
            dest = Path(out_path.strip())
            dest.mkdir(parents=True, exist_ok=True)
            from pycsamt.io.parsers import export_sites
            count = export_sites(sites, dest, fmt=fmt)
            return html.Span(f"✓ Exported {count} file(s) to {dest}",
                             className="text-success")
        except ImportError:
            return _coming_soon("Format Converter")
        except Exception as exc:
            return _err(exc)

    # Dimensionality classifier
    @app.callback(
        Output("tool-dim-out", "children"),
        Input("tool-dim-run",    "n_clicks"),
        State("tool-dim-method", "value"),
        State(IDs.STORE_DATA,    "data"),
        State(IDs.SESSION_ID,    "data"),
        prevent_initial_call=True,
    )
    def _run_dimensionality(n, method, store_data, session_id):
        if not n:
            return no_update
        try:
            from pycsamt.app.web.cache import cache_get
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")
            from pycsamt.emtools.dimensionality import classify_dimensionality
            df = classify_dimensionality(sites, method=method)
            import plotly.graph_objects as go
            counts = df["class"].value_counts()
            fig = go.Figure(go.Bar(
                x=counts.index.tolist(), y=counts.values.tolist(),
                marker_color=["#89b4fa", "#a6e3a1", "#fab387"],
            ))
            _style_fig(fig, "Count", height=200)
            return dcc.Graph(figure=fig, config={"displayModeBar": False})
        except ImportError:
            return _coming_soon("Dimensionality Classifier")
        except Exception as exc:
            return _err(exc)

    # Frequency editor / QC
    @app.callback(
        Output("tool-freq-out", "children"),
        Input("tool-freq-run",    "n_clicks"),
        State("tool-freq-thresh", "value"),
        State(IDs.STORE_DATA,     "data"),
        State(IDs.SESSION_ID,     "data"),
        prevent_initial_call=True,
    )
    def _run_freq_qc(n, thresh, store_data, session_id):
        if not n:
            return no_update
        try:
            from pycsamt.app.web.cache import cache_get
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")
            from pycsamt.emtools.qc import apply_confidence_qc
            report = apply_confidence_qc(sites, min_confidence=thresh)
            dropped = report.get("dropped_bands", 0)
            masked  = report.get("masked_bands", 0)
            return html.Span(
                f"✓ QC complete — {dropped} band(s) dropped, {masked} masked.",
                className="text-success",
            )
        except ImportError:
            return _coming_soon("Frequency Editor")
        except Exception as exc:
            return _err(exc)

    # Batch export
    @app.callback(
        Output("tool-batch-out", "children"),
        Input("tool-batch-run",   "n_clicks"),
        State("tool-batch-path",  "value"),
        State("tool-batch-fmt",   "value"),
        State("tool-batch-dpi",   "value"),
        State("tool-batch-items", "value"),
        State(IDs.STORE_DATA,     "data"),
        State(IDs.SESSION_ID,     "data"),
        prevent_initial_call=True,
    )
    def _run_batch(n, out_path, fmt, dpi, items, store_data, session_id):
        if not n:
            return no_update
        if not out_path or not out_path.strip():
            return _warn("Enter an output folder first.")
        if not items:
            return _warn("Select at least one figure type.")
        try:
            from pathlib import Path
            from pycsamt.app.web.cache import cache_get
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")
            from pycsamt.emtools.plot import batch_export
            dest  = Path(out_path.strip())
            dest.mkdir(parents=True, exist_ok=True)
            count = batch_export(sites, dest, fmt=fmt, dpi=int(dpi), items=items)
            return html.Span(f"✓ Saved {count} figure(s) to {dest}",
                             className="text-success")
        except ImportError:
            return _coming_soon("Batch Export")
        except Exception as exc:
            return _err(exc)

    # Elevation enrichment
    @app.callback(
        Output("tool-elev-out", "children"),
        Input("tool-elev-run",    "n_clicks"),
        State("tool-elev-source", "value"),
        State(IDs.STORE_DATA,     "data"),
        State(IDs.SESSION_ID,     "data"),
        prevent_initial_call=True,
    )
    def _run_elevation(n, source, store_data, session_id):
        if not n:
            return no_update
        try:
            from pycsamt.app.web.cache import cache_get
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")
            from pycsamt.gis.utils import enrich_elevations
            updated = enrich_elevations(sites, source=source)
            return html.Span(
                f"✓ Elevation updated for {updated} station(s).",
                className="text-success",
            )
        except ImportError:
            return _coming_soon("Elevation Enrichment")
        except Exception as exc:
            return _err(exc)

    # Station response inspector
    @app.callback(
        Output("tool-resp-out", "children"),
        Input("tool-resp-run",     "n_clicks"),
        State("tool-resp-station", "value"),
        State("tool-resp-comps",   "value"),
        State("tool-resp-xscale",  "value"),
        State("tool-resp-ytype",   "value"),
        State(IDs.SESSION_ID,      "data"),
        prevent_initial_call=True,
    )
    def _run_response(n, station, comps, xscale, ytype, session_id):
        if not n or not station:
            return no_update
        try:
            import numpy as np
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
            from pycsamt.app.web.cache import cache_get
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")
            site = sites.get(station)
            edi  = site.edi if hasattr(site, "edi") else site
            freqs = edi.Z.freq
            periods = 1.0 / freqs

            _COMP_IDX = {"xx": (0, 0), "xy": (0, 1), "yx": (1, 0), "yy": (1, 1)}
            _COLORS   = {"xx": "#f38ba8", "xy": "#89b4fa",
                         "yx": "#a6e3a1", "yy": "#fab387"}

            if ytype == "rho_phi":
                fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                                    vertical_spacing=0.06)
                for c in (comps or ["xy"]):
                    i, j = _COMP_IDX[c]
                    z    = edi.Z.z[:, i, j]
                    rho  = (np.abs(z) ** 2) / (freqs * 5 * np.pi * 4e-7)
                    phi  = np.degrees(np.angle(z))
                    col  = _COLORS[c]
                    fig.add_trace(go.Scatter(x=periods, y=rho, name=f"Z{c}",
                                             line={"color": col}), row=1, col=1)
                    fig.add_trace(go.Scatter(x=periods, y=phi, name=f"φ Z{c}",
                                             line={"color": col, "dash": "dot"},
                                             showlegend=False), row=2, col=1)
                fig.update_yaxes(title_text="ρ (Ω·m)", type="log", row=1, col=1)
                fig.update_yaxes(title_text="Phase (°)", row=2, col=1)
            else:
                fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                                    vertical_spacing=0.06)
                for c in (comps or ["xy"]):
                    i, j = _COMP_IDX[c]
                    z    = edi.Z.z[:, i, j]
                    col  = _COLORS[c]
                    fig.add_trace(go.Scatter(x=periods, y=z.real, name=f"Re Z{c}",
                                             line={"color": col}), row=1, col=1)
                    fig.add_trace(go.Scatter(x=periods, y=z.imag, name=f"Im Z{c}",
                                             line={"color": col, "dash": "dot"}),
                                  row=2, col=1)
                fig.update_yaxes(title_text="Re Z", row=1, col=1)
                fig.update_yaxes(title_text="Im Z", row=2, col=1)

            fig.update_xaxes(title_text="Period (s)", type=xscale,
                             gridcolor="#313244", row=2, col=1)
            fig.update_layout(
                plot_bgcolor="#1e1e2e", paper_bgcolor="#1e1e2e",
                font={"color": "#cdd6f4"},
                margin={"l": 60, "r": 20, "t": 10, "b": 50},
                height=300,
                legend={"orientation": "h", "y": 1.05},
            )
            return dcc.Graph(figure=fig, config={"displayModeBar": False})
        except Exception as exc:
            return _err(exc)

    # Strike profile viewer
    @app.callback(
        Output("tool-sprof-out", "children"),
        Input("tool-sprof-run",    "n_clicks"),
        State("tool-sprof-line",   "value"),
        State("tool-sprof-method", "value"),
        State("tool-sprof-band",   "value"),
        State("tool-sprof-iqr",    "value"),
        State(IDs.STORE_DATA,      "data"),
        State(IDs.SESSION_ID,      "data"),
        prevent_initial_call=True,
    )
    def _run_strike_profile(n, line, method, band, show_iqr, store_data, session_id):
        if not n:
            return no_update
        try:
            import numpy as np
            import plotly.graph_objects as go
            from pycsamt.app.web.cache import cache_get
            from pycsamt.emtools.strike import compute_strike_profile
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")
            result = compute_strike_profile(
                sites, line=line, method=method, band=band)
            stations = result["stations"]
            angles   = np.array(result["angles"])
            x = list(range(len(stations)))
            fig = go.Figure()
            if show_iqr and "q25" in result:
                q25 = np.array(result["q25"])
                q75 = np.array(result["q75"])
                fig.add_trace(go.Scatter(
                    x=x + x[::-1],
                    y=q75.tolist() + q25[::-1].tolist(),
                    fill="toself", fillcolor="rgba(137,180,250,0.15)",
                    line={"color": "rgba(0,0,0,0)"}, showlegend=False,
                    name="IQR",
                ))
            fig.add_trace(go.Scatter(
                x=x, y=angles.tolist(),
                mode="lines+markers",
                line={"color": "#89b4fa", "width": 2},
                marker={"size": 5},
                name="Strike (°)",
            ))
            fig.update_layout(
                xaxis={"tickvals": x, "ticktext": stations,
                       "tickangle": -45, "gridcolor": "#313244"},
                yaxis={"title": "Strike angle (°)", "gridcolor": "#313244"},
                plot_bgcolor="#1e1e2e", paper_bgcolor="#1e1e2e",
                font={"color": "#cdd6f4"},
                margin={"l": 60, "r": 20, "t": 20, "b": 80},
                height=260,
            )
            return dcc.Graph(figure=fig, config={"displayModeBar": False})
        except ImportError:
            return _coming_soon("Strike Profile Viewer")
        except Exception as exc:
            return _err(exc)

    # Phase tensor map
    @app.callback(
        Output("tool-pt-out", "children"),
        Input("tool-pt-run",     "n_clicks"),
        State("tool-pt-period",  "value"),
        State("tool-pt-colorby", "value"),
        State("tool-pt-cmap",    "value"),
        State("tool-pt-tipper",  "value"),
        State(IDs.SESSION_ID,    "data"),
        prevent_initial_call=True,
    )
    def _run_phase_tensor(n, period, colorby, cmap, show_tipper, session_id):
        if not n or period is None:
            return no_update
        try:
            from pycsamt.app.web.cache import cache_get
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")
            from pycsamt.emtools.tensor import plot_phase_tensor_map
            fig = plot_phase_tensor_map(
                sites,
                period=float(period),
                color_by=colorby,
                cmap=cmap,
                tipper=show_tipper,
                backend="plotly",
            )
            fig.update_layout(
                plot_bgcolor="#1e1e2e", paper_bgcolor="#1e1e2e",
                font={"color": "#cdd6f4"},
                margin={"l": 60, "r": 20, "t": 30, "b": 50},
                height=340,
            )
            return dcc.Graph(figure=fig, config={"displayModeBar": False})
        except ImportError:
            return _coming_soon("Phase Tensor Map")
        except Exception as exc:
            return _err(exc)


# ── Private helpers ───────────────────────────────────────────────────────────

def _style_fig(fig, ytitle: str, height: int = 220) -> None:
    fig.update_layout(
        plot_bgcolor="#1e1e2e", paper_bgcolor="#1e1e2e",
        font={"color": "#cdd6f4"},
        yaxis={"title": ytitle, "gridcolor": "#313244"},
        xaxis={"gridcolor": "#313244"},
        margin={"l": 50, "r": 20, "t": 20, "b": 50},
        height=height,
    )


def _err(exc: Exception) -> html.Span:
    return html.Span(f"✗ {exc}", className="text-danger small")


def _warn(msg: str) -> html.Span:
    return html.Span(f"⚠ {msg}", className="text-warning small")


def _coming_soon(name: str) -> html.Div:
    return html.Div([
        html.I(className="bi bi-clock-history me-2 text-muted"),
        html.Span(f"{name} backend not yet available in web mode.",
                  className="text-muted small"),
    ], className="mt-2")
