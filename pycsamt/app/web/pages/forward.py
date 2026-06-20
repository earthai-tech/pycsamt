# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Forward modelling page — 1D (MT/CSAMT/TEM), 2D MT, 3D MT.

Layout
------
Sidebar:
  • Sticky Run-Forward bar (always visible, top of sidebar)
  • Dimension selector (1-D | 2-D MT | 3-D MT) as pill buttons
  • Scrollable per-dimension controls + shared freq range + save

View area (right):
  • Fixed tab bar: 1-D · 2-D MT · 3-D MT  (persistent per-tab figures)
  • Context info bar
  • Per-tab panels — each with optional mini plot-selector bar + figure
"""
from __future__ import annotations

from dash import dash_table, dcc, html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.forward_controller import (
    ForwardController, GEOLOGY_PRESET_NAMES,
)
from pycsamt.app.web.layout import IDs, _command_bar, _icon
from pycsamt.app.web.utils import empty_src

PAGE_ID = "forward"

# ── Selector options ──────────────────────────────────────────────────────────

_METHOD_OPTS = [
    {"label": "MT 1-D",    "value": "MT1D"},
    {"label": "CSAMT 1-D", "value": "CSAMT1D"},
    {"label": "TEM 1-D",   "value": "TEM1D"},
]
_PRESET_OPTS = [{"label": p.replace("_", " ").title(), "value": p}
                for p in GEOLOGY_PRESET_NAMES]

_DEFAULT_LAYERS = [
    {"layer": 1, "resistivity": 50,   "thickness": 200},
    {"layer": 2, "resistivity": 500,  "thickness": 800},
    {"layer": 3, "resistivity": 20,   "thickness": 1000},
]

_2D_MODEL_OPTS = [
    {"label": "Halfspace",                  "value": "halfspace"},
    {"label": "Conductive/Resistive Anomaly", "value": "anomaly"},
    {"label": "Extend 1-D Layers",          "value": "from_layers"},
    {"label": "Random Layered",             "value": "random"},
]
_2D_PLOT_OPTS = [
    {"label": "2-D Resistivity Model",    "value": "model"},
    {"label": "ρₐ Pseudosection (TE)",    "value": "psd_te_rho"},
    {"label": "ρₐ Pseudosection (TM)",    "value": "psd_tm_rho"},
    {"label": "Phase Pseudosection (TE)", "value": "psd_te_phi"},
    {"label": "Phase Pseudosection (TM)", "value": "psd_tm_phi"},
    {"label": "Profiles (TE)",            "value": "profiles_te"},
    {"label": "Profiles (TM)",            "value": "profiles_tm"},
]

_3D_MODEL_OPTS = [
    {"label": "Halfspace",         "value": "halfspace"},
    {"label": "3-D Block Anomaly", "value": "block"},
    {"label": "Random Layered",    "value": "random"},
]
_3D_PLOT_OPTS = [
    {"label": "3-D Model Slices",  "value": "model"},
    {"label": "Response Map",      "value": "map"},
    {"label": "Response Section",  "value": "section"},
    {"label": "Tensor Components", "value": "tensor"},
]

_FWD_TABS = [
    ("fwd-tab-1d", "1-D",    "bi-reception-4"),
    ("fwd-tab-2d", "2-D MT", "bi-grid"),
    ("fwd-tab-3d", "3-D MT", "bi-cube"),
]
_FWD_DEFAULT = "fwd-tab-1d"

_IMG_STYLE = {"width": "100%", "height": "100%", "objectFit": "contain"}


# ── Small helpers ─────────────────────────────────────────────────────────────

def _num(id_, value, **kw):
    return dbc.Input(id=id_, type="number", value=value, size="sm", **kw)


def _lbl(text):
    return html.Small(text, className="text-muted")


def _clabel(text):
    return html.Div(text, className="ctrl-label")


# ── Layer table (inside sidebar 1D section) ───────────────────────────────────

def _layer_table() -> html.Div:
    tbl = dash_table.DataTable(
        id=IDs.FWD_LAYER_TABLE,
        columns=[
            {"name": "#",                 "id": "layer",       "editable": False},
            {"name": "Resistivity (Ω·m)", "id": "resistivity", "editable": True,
             "type": "numeric"},
            {"name": "Thickness (m)",     "id": "thickness",   "editable": True,
             "type": "numeric"},
        ],
        data=_DEFAULT_LAYERS,
        editable=True,
        row_deletable=True,
        style_table={"overflowY": "auto", "maxHeight": "180px"},
        style_header={
            "backgroundColor": "#11111b", "color": "#6c7086",
            "fontSize": "11px", "border": "none",
            "textTransform": "uppercase", "fontWeight": "600",
        },
        style_cell={
            "backgroundColor": "#181825", "color": "#cdd6f4",
            "border": "1px solid #313244", "fontSize": "11px",
            "padding": "4px 8px",
        },
        style_data_conditional=[
            {"if": {"state": "active"},
             "backgroundColor": "#313244", "border": "1px solid #89b4fa"},
        ],
    )
    return html.Div(tbl, className="fwd-layer-table")


# ── Sidebar sections (per dimension) ─────────────────────────────────────────

def _1d_controls():
    """Sidebar controls for 1-D forward (method + preset + layer model + method cards)."""
    return html.Div([
        # Method selector
        html.Div([
            _clabel("Forward Method"),
            dbc.Select(id=IDs.FWD_METHOD, options=_METHOD_OPTS,
                       value="MT1D", size="sm"),
            dcc.Store(id=IDs.FWD_DIM, data="MT1D"),
        ], className="ctrl-card"),

        # Geological preset
        html.Div([
            _clabel("Geological Preset"),
            dbc.Select(id=IDs.FWD_PRESET, options=_PRESET_OPTS,
                       value=_PRESET_OPTS[0]["value"], size="sm"),
            dbc.Button(
                [html.I(className="bi bi-magic me-1"), "Load Preset"],
                id=IDs.BTN_FWD_PRESET, color="secondary",
                size="sm", className="w-100 mt-2",
            ),
        ], className="ctrl-card"),

        # Layer model editor
        html.Div([
            _clabel("Layer Model"),
            _layer_table(),
            dbc.Button(
                [html.I(className="bi bi-plus-lg me-1"), "Add Layer"],
                id=IDs.BTN_FWD_ADD_LAYER, color="secondary", size="sm",
                className="mt-2",
            ),
            # Halfspace
            html.Div([
                html.Div([
                    html.Span(
                        [html.I(className="bi bi-infinity me-1"), "Halfspace"],
                        className="fwd-hs-label",
                    ),
                    html.Span(" — infinite depth",
                              style={"fontSize": "10px", "color": "var(--sub0)"}),
                ], className="mb-1 mt-2"),
                dbc.InputGroup([
                    dbc.InputGroupText("ρ", style={"fontSize": "11px"}),
                    _num(IDs.FWD_HALFSPACE_RHO, 1000, min=0.01, step=10),
                    dbc.InputGroupText("Ω·m", style={"fontSize": "11px"}),
                ], size="sm"),
            ], className="fwd-halfspace-row"),
        ], className="ctrl-card"),

        # CSAMT source offset (shown only for CSAMT1D)
        html.Div([
            _clabel("CSAMT Source Offset"),
            dbc.InputGroup([
                _num(IDs.FWD_OFFSET, 5000, min=100, step=100),
                dbc.InputGroupText("m", style={"fontSize": "11px"}),
            ], size="sm"),
        ], className="ctrl-card", id=IDs.FWD_CSAMT_CARD, style={"display": "none"}),

        # TEM loop radius (shown only for TEM1D)
        html.Div([
            _clabel("TEM Loop Radius"),
            dbc.InputGroup([
                _num(IDs.FWD_LOOP_R, 50, min=5, step=5),
                dbc.InputGroupText("m", style={"fontSize": "11px"}),
            ], size="sm"),
        ], className="ctrl-card", id=IDs.FWD_TEM_CARD, style={"display": "none"}),
    ], id="fwd-ctrl-1d")


def _2d_controls():
    """Sidebar controls for 2-D MT forward."""
    return html.Div([
        html.Div([
            _clabel("Model Type"),
            dbc.Select(id=IDs.FWD2_MODEL_TYPE, options=_2D_MODEL_OPTS,
                       value="halfspace", size="sm"),
        ], className="ctrl-card"),

        html.Div([
            _clabel("Background & Grid"),
            dbc.Row([
                dbc.Col([_lbl("ρ bg (Ω·m)"),
                         _num(IDs.FWD2_BG_RHO, 100, min=0.01)], width=6),
                dbc.Col([_lbl("Stations"),
                         _num(IDs.FWD2_N_STATIONS, 10, min=3, max=50, step=1)], width=6),
            ], className="g-1 mb-1"),
            dbc.Row([
                dbc.Col([_lbl("Nx"),  _num(IDs.FWD2_NX, 25, min=10, max=80, step=5)]),
                dbc.Col([_lbl("Nz"),  _num(IDs.FWD2_NZ, 20, min=10, max=60, step=5)]),
            ], className="g-1 mb-1"),
            dbc.Row([
                dbc.Col([_lbl("X max (m)"), _num(IDs.FWD2_X_MAX, 8000, min=1000, step=1000)]),
                dbc.Col([_lbl("Z max (m)"), _num(IDs.FWD2_Z_MAX, 5000, min=500,  step=500)]),
            ], className="g-1"),
        ], className="ctrl-card"),

        # Anomaly params (conditional)
        html.Div([
            _clabel("Anomaly"),
            dbc.Row([
                dbc.Col([_lbl("ρ (Ω·m)"), _num(IDs.FWD2_ANOMALY_RHO, 5, min=0.01)]),
            ], className="g-1 mb-1"),
            html.Div("X bounds (m):", className="fwd-bound-label"),
            dbc.Row([
                dbc.Col([_lbl("X lo"), _num(IDs.FWD2_AX_LO, 2000, step=100)]),
                dbc.Col([_lbl("X hi"), _num(IDs.FWD2_AX_HI, 4000, step=100)]),
            ], className="g-1 mb-1"),
            html.Div("Z bounds (m):", className="fwd-bound-label"),
            dbc.Row([
                dbc.Col([_lbl("Z top"), _num(IDs.FWD2_AZ_LO, 200, step=50)]),
                dbc.Col([_lbl("Z bot"), _num(IDs.FWD2_AZ_HI, 800, step=50)]),
            ], className="g-1"),
        ], className="ctrl-card", id=IDs.FWD2_ANOMALY_CARD, style={"display": "none"}),

        # Random params (conditional)
        html.Div([
            _clabel("Random Model"),
            dbc.Row([
                dbc.Col([_lbl("N layers"), _num(IDs.FWD2_N_LAYERS, 4, min=2, max=10, step=1)]),
                dbc.Col([_lbl("Seed"),     _num(IDs.FWD2_SEED, 42, min=0)]),
            ], className="g-1"),
        ], className="ctrl-card", id=IDs.FWD2_RANDOM_CARD, style={"display": "none"}),
    ], id="fwd-ctrl-2d", style={"display": "none"})


def _3d_controls():
    """Sidebar controls for 3-D MT forward."""
    return html.Div([
        html.Div([
            _clabel("Model Type"),
            dbc.Select(id=IDs.FWD3_MODEL_TYPE, options=_3D_MODEL_OPTS,
                       value="halfspace", size="sm"),
        ], className="ctrl-card"),

        html.Div([
            _clabel("Background & Grid"),
            dbc.Row([
                dbc.Col([_lbl("ρ bg (Ω·m)"), _num(IDs.FWD3_BG_RHO, 100, min=0.01)], width=6),
                dbc.Col([
                    _lbl("Stations (Nx×Ny)"),
                    dbc.Row([
                        dbc.Col(_num(IDs.FWD3_NX_ST, 4, min=2, max=15, step=1)),
                        dbc.Col(_num(IDs.FWD3_NY_ST, 4, min=2, max=15, step=1)),
                    ], className="g-1"),
                ], width=6),
            ], className="g-1 mb-1"),
            dbc.Row([
                dbc.Col([_lbl("Nx"), _num(IDs.FWD3_NX, 15, min=8, max=40, step=5)]),
                dbc.Col([_lbl("Ny"), _num(IDs.FWD3_NY, 15, min=8, max=40, step=5)]),
                dbc.Col([_lbl("Nz"), _num(IDs.FWD3_NZ, 12, min=6, max=30, step=3)]),
            ], className="g-1 mb-1"),
            dbc.Row([
                dbc.Col([_lbl("X max (m)"), _num(IDs.FWD3_X_MAX, 8000, min=1000, step=1000)]),
                dbc.Col([_lbl("Y max (m)"), _num(IDs.FWD3_Y_MAX, 8000, min=1000, step=1000)]),
                dbc.Col([_lbl("Z max (m)"), _num(IDs.FWD3_Z_MAX, 5000, min=500,  step=500)]),
            ], className="g-1"),
        ], className="ctrl-card"),

        # Block anomaly (conditional)
        html.Div([
            _clabel("Block Anomaly"),
            dbc.Row([
                dbc.Col([_lbl("ρ (Ω·m)"), _num(IDs.FWD3_ANOMALY_RHO, 5, min=0.01)]),
            ], className="g-1 mb-1"),
            html.Div("X bounds (m):", className="fwd-bound-label"),
            dbc.Row([
                dbc.Col([_lbl("X lo"), _num(IDs.FWD3_AX_LO, 2000, step=200)]),
                dbc.Col([_lbl("X hi"), _num(IDs.FWD3_AX_HI, 6000, step=200)]),
            ], className="g-1 mb-1"),
            html.Div("Y bounds (m):", className="fwd-bound-label"),
            dbc.Row([
                dbc.Col([_lbl("Y lo"), _num(IDs.FWD3_AY_LO, 2000, step=200)]),
                dbc.Col([_lbl("Y hi"), _num(IDs.FWD3_AY_HI, 6000, step=200)]),
            ], className="g-1 mb-1"),
            html.Div("Z bounds (m):", className="fwd-bound-label"),
            dbc.Row([
                dbc.Col([_lbl("Z top"), _num(IDs.FWD3_AZ_LO, 300, step=100)]),
                dbc.Col([_lbl("Z bot"), _num(IDs.FWD3_AZ_HI, 1500, step=100)]),
            ], className="g-1"),
        ], className="ctrl-card", id=IDs.FWD3_BLOCK_CARD, style={"display": "none"}),

        # Random params (conditional)
        html.Div([
            _clabel("Random Model"),
            dbc.Row([
                dbc.Col([_lbl("N layers"), _num(IDs.FWD3_N_LAYERS, 4, min=2, max=8, step=1)]),
                dbc.Col([_lbl("Seed"),     _num(IDs.FWD3_SEED, 42, min=0)]),
            ], className="g-1"),
        ], className="ctrl-card", id=IDs.FWD3_RANDOM_CARD, style={"display": "none"}),
    ], id="fwd-ctrl-3d", style={"display": "none"})


# ── View panel helpers ────────────────────────────────────────────────────────

def _view_panel_1d() -> html.Div:
    return html.Div(
        html.Div(
            html.Img(id=IDs.IMG_FWD_1D, src=empty_src(), style=_IMG_STYLE),
            className="fig-img-wrap profile-page-fig",
        ),
        id=IDs.FWD_1D_PANEL,
        className="prof-panel",
        style={"display": "flex"},
    )


def _view_panel_2d() -> html.Div:
    """2-D view: mini plot-type bar at top + persistent figure."""
    mini_bar = html.Div([
        html.Span("View:", className="fwd-plot-bar-label"),
        dbc.Select(
            id=IDs.FWD2_PLOT_TYPE,
            options=_2D_PLOT_OPTS,
            value="model",
            size="sm",
            style={"width": "260px"},
        ),
        html.Span(
            "Run the solver first, then switch views.",
            style={"fontSize": "10px", "color": "var(--sub0)", "marginLeft": "8px"},
        ),
    ], className="fwd-plot-bar")
    return html.Div(
        [
            mini_bar,
            html.Div(
                html.Img(id=IDs.IMG_FWD_2D, src=empty_src(), style=_IMG_STYLE),
                className="fig-img-wrap profile-page-fig",
            ),
        ],
        id=IDs.FWD_2D_PANEL,
        className="prof-panel",
        style={"display": "none"},
    )


def _view_panel_3d() -> html.Div:
    """3-D view: mini selectors bar + persistent figure."""
    mini_bar = html.Div([
        html.Span("Type:", className="fwd-plot-bar-label"),
        dbc.Select(id=IDs.FWD3_PLOT_TYPE, options=_3D_PLOT_OPTS,
                   value="model", size="sm", style={"width": "170px"}),
        html.Span("Comp:", className="fwd-plot-bar-label ms-3"),
        dbc.Select(id=IDs.FWD3_COMPONENT,
                   options=[{"label": "Zxy", "value": "xy"},
                             {"label": "Zyx", "value": "yx"}],
                   value="xy", size="sm", style={"width": "90px"}),
        html.Span("Freq #:", className="fwd-plot-bar-label ms-3"),
        dbc.Input(id=IDs.FWD3_FREQ_IDX, type="number", value=0, min=0, step=1,
                  size="sm", style={"width": "70px"}),
    ], className="fwd-plot-bar")
    return html.Div(
        [
            mini_bar,
            html.Div(
                html.Img(id=IDs.IMG_FWD_3D, src=empty_src(), style=_IMG_STYLE),
                className="fig-img-wrap profile-page-fig",
            ),
        ],
        id=IDs.FWD_3D_PANEL,
        className="prof-panel",
        style={"display": "none"},
    )


# ── Page layout ───────────────────────────────────────────────────────────────

def layout() -> html.Div:

    # ── Left sidebar ──────────────────────────────────────────────────────────

    # Active-tab store (inside sidebar so it's mounted with the page)
    store = dcc.Store(id=IDs.FWD_ACTIVE_TAB, data=_FWD_DEFAULT)

    # Sticky run bar — always visible regardless of scroll
    run_bar = html.Div([
        dbc.Button(
            [html.I(className="bi bi-play-fill me-1"), "Run Forward"],
            id=IDs.BTN_FWD_RUN,
            color="primary", size="sm",
            className="w-100 mb-1",
            n_clicks=0,
        ),
        dbc.Spinner(html.Div(id=IDs.FWD_SPINNER), size="sm", color="primary"),
        html.Div(id=IDs.FWD_FEEDBACK, className="fwd-feedback-mini"),
    ], className="fwd-run-bar")

    # Dimension selector pill-buttons
    dim_bar = html.Div([
        html.Button(
            [html.I(className=f"bi {icon} me-1"), label],
            id=f"fwd-dim-btn-{key.split('-')[2]}",
            className=f"fwd-dim-btn{' active' if key == _FWD_DEFAULT else ''}",
            n_clicks=0,
        )
        for key, label, icon in _FWD_TABS
    ], className="fwd-dim-bar")

    # Shared: frequency / time range
    freq_card = html.Div([
        html.Div(id=IDs.FWD_FREQ_LABEL,
                 children="Frequency range [log₁₀ Hz]",
                 className="ctrl-label"),
        dbc.Row([
            dbc.Col([_lbl("Min"), _num(IDs.FWD_FREQ_MIN, -3, step=0.5)]),
            dbc.Col([_lbl("Max"), _num(IDs.FWD_FREQ_MAX,  4, step=0.5)]),
            dbc.Col([_lbl("N pts"), _num(IDs.FWD_N_FREQ, 20, min=5, max=100, step=5)]),
        ], className="g-1"),
    ], className="ctrl-card")

    # Save card
    save_card = html.Div([
        _clabel("Save Model"),
        dbc.Row([
            dbc.Col(dbc.Input(id=IDs.FWD_MODEL_NAME,
                              placeholder="Model name…", size="sm")),
            dbc.Col(
                dbc.Button(
                    [_icon("save", size=12, cls=""), " Save"],
                    id=IDs.BTN_FWD_SAVE, color="secondary", size="sm",
                ),
                width="auto",
            ),
        ], className="g-1"),
    ], className="ctrl-card")

    # Scrollable controls area
    ctrl_scroll = html.Div([
        dim_bar,
        _1d_controls(),
        _2d_controls(),
        _3d_controls(),
        freq_card,
        save_card,
    ], className="fwd-ctrl-scroll")

    sidebar = html.Div(
        [store, run_bar, ctrl_scroll],
        className="analysis-controls fwd-sidebar",
    )

    # ── View area (right side) ────────────────────────────────────────────────

    tab_bar = html.Div([
        html.Button(
            [html.I(className=f"bi {icon} me-1"), label],
            id=f"fwd-view-tab-{key.split('-')[2]}",
            className=f"prof-tab-btn{' active' if key == _FWD_DEFAULT else ''}",
            n_clicks=0,
        )
        for key, label, icon in _FWD_TABS
    ], className="prof-tab-bar")

    ctx_bar = html.Div(
        [
            html.Span(
                [html.I(className="bi bi-reception-4 me-1"), "1-D"],
                className="adv-info-group",
            ),
            html.Span("·", className="adv-info-sep"),
            html.Span(
                "Configure the model then click Run Forward",
                style={"color": "var(--sub0)", "fontSize": "10.5px"},
            ),
        ],
        id=IDs.FWD_CTX_BAR,
        className="adv-info-bar",
    )

    panels = html.Div([
        _view_panel_1d(),
        _view_panel_2d(),
        _view_panel_3d(),
    ], className="prof-view-panel")

    view_area = html.Div(
        [tab_bar, ctx_bar, panels],
        className="adv-view-area",
    )

    # ── Final assembly ────────────────────────────────────────────────────────

    return html.Div([
        _command_bar(
            "forward", "Forward Modelling",
            "MT1D · CSAMT1D · TEM1D · MT2D · MT3D (quasi-3D)",
        ),
        html.Div(
            [sidebar, view_area],
            className="analysis-layout",
            style={"flex": "1"},
        ),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # handled by callbacks/forward.py
