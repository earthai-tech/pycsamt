# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Data correction / conditioning page — v2 redesign.

Layout anatomy
--------------
  ┌─ command bar ────────────────────────────────────────────────────────────┐
  │  (page title + subtitle — 36 px, flex-shrink:0)                         │
  └──────────────────────────────────────────────────────────────────────────┘
  ┌─ corr-main (flex-row, flex:1) ──────────────────────────────────────────┐
  │ ┌─ corr-sidebar ──────┐  ┌─ corr-view-area ──────────────────────────┐  │
  │ │  corr-run-bar       │  │  corr-ctx-bar  (chain info strip)         │  │
  │ │  (Preview/Apply/    │  │  corr-tab-bar  (4 pill buttons)           │  │
  │ │   Undo/Reset +      │  │  corr-panels-wrap (4 panels in DOM)       │  │
  │ │   spinner+feedback) │  └───────────────────────────────────────────┘  │
  │ │  fwd-ctrl-scroll    │                                                  │
  │ │  (method picker,    │                                                  │
  │ │   params, scope,    │                                                  │
  │ │   drop-freq card,   │                                                  │
  │ │   chain, export)    │                                                  │
  │ └─────────────────────┘                                                  │
  └──────────────────────────────────────────────────────────────────────────┘

All four view panels are always in the DOM; visibility is toggled by a
clientside callback driven by the `corr-active-tab` dcc.Store.
"""
from __future__ import annotations

from dash import dcc, html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.correction_controller import CATALOGUE, CATEGORIES
from pycsamt.app.web.layout import IDs, _command_bar
from pycsamt.app.web.utils import empty_src

PAGE_ID = "correction"

_CAT_OPTS    = [{"label": c, "value": c} for c in CATEGORIES]
_FIRST_METHS = [{"label": k, "value": k} for k in CATALOGUE[CATEGORIES[0]]]
_N_METHODS   = sum(len(v) for v in CATALOGUE.values())

_EXPORT_OPTS = [
    {"label": "CSV (tabular)",     "value": "csv"},
    {"label": "EDI (per station)", "value": "edi"},
    {"label": "XYZ (profile)",    "value": "xyz"},
]

_CORR_TABS = [
    ("ba",      "Before / After", "bi-layout-split"),
    ("rho-phi", "ρ/φ View",       "bi-graph-up-arrow"),
    ("overlay", "Overlay",        "bi-layers"),
    ("diff",    "Difference",     "bi-arrows-collapse"),
]


def _lbl(txt: str) -> html.Div:
    return html.Div(txt, style={"fontSize": "10px", "color": "#6c7086",
                                "marginBottom": "2px"})


# ── Sticky run bar ────────────────────────────────────────────────────────────

def _run_bar() -> html.Div:
    """Action bar fixed at the top of the sidebar — always visible."""
    return html.Div([
        # Row 1: Preview | Apply
        html.Div([
            dbc.Button(
                [html.I(className="bi bi-eye me-1"), "Preview"],
                id=IDs.BTN_CORR_PREVIEW, color="secondary", size="sm",
                className="w-100",
                title="Preview this correction — shows before vs after without committing",
            ),
            dbc.Button(
                [html.I(className="bi bi-check-lg me-1"), "Apply"],
                id=IDs.BTN_CORR_APPLY, color="success", size="sm",
                className="w-100",
                title="Add this correction step to the chain",
            ),
        ], className="d-flex gap-1 mb-1"),
        # Row 2: Undo | Reset All
        html.Div([
            dbc.Button(
                [html.I(className="bi bi-arrow-counterclockwise me-1"), "Undo"],
                id=IDs.BTN_CORR_UNDO, color="warning", size="sm",
                className="w-100",
                title="Remove the last correction step from the chain",
            ),
            dbc.Button(
                [html.I(className="bi bi-x-circle me-1"), "Reset All"],
                id=IDs.BTN_CORR_RESET, color="danger", size="sm",
                className="w-100",
                title="Clear all corrections and return to raw data",
            ),
        ], className="d-flex gap-1"),
        # Spinner + feedback
        html.Div(
            dbc.Spinner(html.Div(id=IDs.CORR_SPINNER), size="sm", color="success"),
            className="corr-run-spinner",
        ),
        html.Div(id=IDs.CORR_FEEDBACK, className="corr-feedback-text"),
    ], className="corr-run-bar")


# ── Scrollable controls ───────────────────────────────────────────────────────

def _ctrl_scroll() -> html.Div:
    """Scrollable controls panel below the sticky run bar."""
    return html.Div([

        # ── Category ──────────────────────────────────────────────────────
        html.Div([
            html.Div("Category", className="ctrl-label"),
            dbc.Select(id=IDs.CORR_CAT, options=_CAT_OPTS,
                       value=_CAT_OPTS[0]["value"], size="sm"),
        ], className="ctrl-card"),

        # ── Method + description ──────────────────────────────────────────
        html.Div([
            html.Div("Correction Method", className="ctrl-label"),
            dbc.Select(id=IDs.CORR_METHOD, options=_FIRST_METHS,
                       value=_FIRST_METHS[0]["value"], size="sm"),
            html.Div(id=IDs.CORR_DESC, className="corr-desc mt-1"),
        ], className="ctrl-card"),

        # ── Dynamic parameters ────────────────────────────────────────────
        html.Div([
            html.Div("Parameters", className="ctrl-label"),
            html.Div(id=IDs.CORR_PARAMS_PANEL),
        ], className="ctrl-card", id="corr-params-card"),

        # ── Correction scope (frequency band) ─────────────────────────────
        html.Div([
            html.Div("Frequency Scope", className="ctrl-label"),
            dbc.Row([
                dbc.Col([
                    _lbl("f low (Hz)"),
                    dbc.Input(
                        id=IDs.CORR_SCOPE_FREQ_LO, type="number",
                        placeholder="all", value=None,
                        min=1e-5, step=0.001, size="sm", debounce=True,
                    ),
                ]),
                dbc.Col([
                    _lbl("f high (Hz)"),
                    dbc.Input(
                        id=IDs.CORR_SCOPE_FREQ_HI, type="number",
                        placeholder="all", value=None,
                        min=1e-5, step=0.001, size="sm", debounce=True,
                    ),
                ]),
            ], className="g-1"),
            html.Div("Leave blank to apply across all frequencies.",
                     className="corr-hint-txt"),
        ], className="ctrl-card"),

        # ── Drop Bad Frequencies — quick card ─────────────────────────────
        html.Div([
            html.Div([
                html.I(className="bi bi-funnel me-1 corr-drop-icon"),
                html.Span("Drop Bad Frequencies", className="ctrl-label mb-0"),
                html.Span("Quick", className="corr-quick-badge"),
            ], className="d-flex align-items-center gap-1 mb-2"),
            dbc.Row([
                dbc.Col([
                    _lbl("SNR threshold"),
                    dbc.Input(
                        id=IDs.CORR_SNR_THRESH, type="number",
                        value=2.5, min=0.5, max=20.0, step=0.5,
                        size="sm", debounce=True,
                    ),
                ]),
                dbc.Col([
                    _lbl("Min stn frac."),
                    dbc.Input(
                        id=IDs.CORR_MIN_FRAC, type="number",
                        value=0.4, min=0.1, max=1.0, step=0.05,
                        size="sm", debounce=True,
                    ),
                ]),
            ], className="g-1 mb-2"),
            dbc.Button(
                [html.I(className="bi bi-funnel-fill me-1"), "Apply — Drop Bad Freq"],
                id=IDs.BTN_CORR_DROP_FREQ,
                color="info", size="sm", outline=True, className="w-100",
            ),
            html.Div(
                "Removes frequency rows with SNR < threshold in > (1 − min_frac) stations.",
                className="corr-hint-txt",
            ),
        ], className="ctrl-card corr-drop-freq-card"),

        # ── Smooth ρ/φ — quick card ────────────────────────────────────────
        html.Div([
            html.Div([
                html.I(className="bi bi-activity me-1 corr-smooth-icon"),
                html.Span("Smooth ρ/φ Trend", className="ctrl-label mb-0"),
                html.Span("Quick", className="corr-quick-badge"),
            ], className="d-flex align-items-center gap-1 mb-2"),
            dbc.Row([
                dbc.Col([
                    _lbl("Components"),
                    dbc.Select(
                        id=IDs.CORR_SMOOTH_COMP,
                        value="offdiag",
                        options=[
                            {"label": "Off-diag (XY+YX)", "value": "offdiag"},
                            {"label": "All",              "value": "all"},
                            {"label": "XY only",          "value": "xy"},
                            {"label": "YX only",          "value": "yx"},
                        ],
                        size="sm",
                    ),
                ], width=12, className="mb-1"),
            ], className="g-1"),
            dbc.Row([
                dbc.Col([
                    _lbl("Poly degree"),
                    dbc.Input(
                        id=IDs.CORR_SMOOTH_DEGREE,
                        type="number", value=3, min=0, max=8, step=1,
                        size="sm", debounce=True,
                    ),
                ]),
                dbc.Col([
                    _lbl("Blend"),
                    dbc.Input(
                        id=IDs.CORR_SMOOTH_BLEND,
                        type="number", value=1.0, min=0.0, max=1.0, step=0.05,
                        size="sm", debounce=True,
                    ),
                ]),
            ], className="g-1 mb-1"),
            dbc.Row([
                dbc.Col(
                    dbc.Checklist(
                        id=IDs.CORR_SMOOTH_ROBUST,
                        options=[{"label": "Robust fit (spike-resistant)", "value": "robust"}],
                        value=["robust"],
                        switch=True,
                        className="small",
                    ),
                    width=12,
                ),
            ], className="mb-2"),
            dbc.Row([
                dbc.Col(
                    dbc.Button(
                        [html.I(className="bi bi-eye me-1"), "Preview"],
                        id=IDs.BTN_CORR_SMOOTH_PREVIEW,
                        color="warning", size="sm", outline=True, className="w-100",
                    ),
                ),
                dbc.Col(
                    dbc.Button(
                        [html.I(className="bi bi-check2-circle me-1"), "Apply"],
                        id=IDs.BTN_CORR_SMOOTH_APPLY,
                        color="success", size="sm", outline=True, className="w-100",
                    ),
                ),
            ], className="g-1"),
            html.Div(
                "Fits a polynomial to log10(ρ_a) and phase; rebuilds Z so both remain consistent.",
                className="corr-hint-txt mt-1",
            ),
        ], className="ctrl-card corr-smooth-card"),

        # ── Manual Frequency Selection — quick card ────────────────────────
        html.Div([
            html.Div([
                html.I(className="bi bi-sliders me-1 corr-mfreq-icon"),
                html.Span("Manual Frequency Drop", className="ctrl-label mb-0"),
                html.Span("Quick", className="corr-quick-badge"),
            ], className="d-flex align-items-center gap-1 mb-2"),
            # Load + select-all / none row
            dbc.Row([
                dbc.Col(
                    dbc.Button(
                        [html.I(className="bi bi-arrow-repeat me-1"), "Load Frequencies"],
                        id=IDs.BTN_CORR_MFREQ_LOAD,
                        color="secondary", size="sm", outline=True, className="w-100",
                    ),
                    width=12,
                ),
            ], className="g-1 mb-1"),
            # Status line + select-all / none
            html.Div([
                html.Span(
                    id=IDs.CORR_MFREQ_STATUS,
                    children="Click 'Load Frequencies' to populate the list.",
                    className="small text-muted flex-grow-1",
                ),
                html.A("All",  id=IDs.BTN_CORR_MFREQ_ALL,  href="#",
                       className="small corr-mfreq-sel-link me-1"),
                html.Span("·", className="small text-muted"),
                html.A("None", id=IDs.BTN_CORR_MFREQ_NONE, href="#",
                       className="small corr-mfreq-sel-link ms-1"),
            ], className="d-flex align-items-center gap-1 mb-1"),
            # Scrollable frequency checklist
            html.Div(
                dcc.Checklist(
                    id=IDs.CORR_MFREQ_LIST,
                    options=[],
                    value=[],
                    inline=False,
                    labelStyle={"display": "flex", "alignItems": "center",
                                "gap": "6px", "padding": "1px 0"},
                    inputStyle={"cursor": "pointer", "accentColor": "var(--red)"},
                    className="small corr-mfreq-checklist",
                ),
                className="corr-mfreq-scroll",
            ),
            # Action buttons
            dbc.Row([
                dbc.Col(
                    dbc.Button(
                        [html.I(className="bi bi-eye me-1"), "Preview"],
                        id=IDs.BTN_CORR_MFREQ_PREVIEW,
                        color="warning", size="sm", outline=True, className="w-100",
                    ),
                ),
                dbc.Col(
                    dbc.Button(
                        [html.I(className="bi bi-trash3 me-1"), "Apply Drop"],
                        id=IDs.BTN_CORR_MFREQ_APPLY,
                        color="danger", size="sm", outline=True, className="w-100",
                    ),
                ),
            ], className="g-1 mt-2"),
            html.Div(
                "Checked frequencies will be set to NaN in all stations when applied.",
                className="corr-hint-txt mt-1",
            ),
        ], className="ctrl-card corr-mfreq-card"),

        # ── Applied corrections chain ──────────────────────────────────────
        html.Div([
            html.Div([
                html.Div("Correction Chain", className="ctrl-label mb-0"),
                html.Span(id="corr-chain-badge", className="corr-chain-badge"),
            ], className="d-flex justify-content-between align-items-center mb-1"),
            html.Div(
                id=IDs.CORR_STACK,
                children=[html.Span("No corrections applied yet.",
                                    className="text-muted small")],
            ),
        ], className="ctrl-card"),

        # ── Export ────────────────────────────────────────────────────────
        html.Div([
            html.Div("Export Corrected Data", className="ctrl-label"),
            dbc.Row([
                dbc.Col([
                    _lbl("Format"),
                    dbc.Select(id=IDs.CORR_EXPORT_FMT, options=_EXPORT_OPTS,
                               value="csv", size="sm"),
                ]),
                dbc.Col([
                    _lbl(" "),
                    dbc.Button(
                        [html.I(className="bi bi-download me-1"), "Export"],
                        id=IDs.BTN_CORR_EXPORT, color="info", size="sm",
                        outline=True, className="w-100",
                        style={"marginTop": "14px"},
                    ),
                ]),
            ], className="g-1"),
        ], className="ctrl-card"),

    ], className="fwd-ctrl-scroll")


# ── Context info bar ──────────────────────────────────────────────────────────

def _ctx_bar() -> html.Div:
    return html.Div(
        id="corr-ctx-bar",
        children=[
            html.I(className="bi bi-info-circle me-1"),
            "No corrections applied — chain is empty.",
        ],
        className="corr-ctx-bar",
    )


# ── Tab pill bar ──────────────────────────────────────────────────────────────

def _tab_bar() -> html.Div:
    """Four pill buttons; clientside callback shows the matching panel."""
    return html.Div([
        html.Button(
            [html.I(className=f"bi {icon} me-1"), label],
            id=f"corr-tab-btn-{slug}",
            className="corr-tab-btn" + (" active" if slug == "ba" else ""),
            n_clicks=0,
        )
        for slug, label, icon in _CORR_TABS
    ], className="corr-tab-bar")


# ── View panels ───────────────────────────────────────────────────────────────

def _view_panels() -> html.Div:
    """All four panels always rendered; clientside toggling via display."""
    def _img(img_id):
        return html.Img(
            id=img_id, src=empty_src(dark=True),
            style={"width": "100%", "height": "100%", "objectFit": "contain"},
        )

    return html.Div([

        # ── Before / After ─────────────────────────────────────────────────
        html.Div([
            html.Div([
                html.Div([
                    html.Div(id="corr-label-left",
                             children=[html.I(className="bi bi-camera me-1"),
                                       "Raw / Before"],
                             className="corr-panel-label"),
                    html.Div(_img(IDs.IMG_CORR_LEFT),
                             className="fig-area", style={"flex": "1"}),
                ], className="corr-half-panel"),
                html.Div([
                    html.Div(id="corr-label-right",
                             children=[html.I(className="bi bi-stars me-1"),
                                       "Corrected / After"],
                             className="corr-panel-label"),
                    html.Div(_img(IDs.IMG_CORR_RIGHT),
                             className="fig-area", style={"flex": "1"}),
                ], className="corr-half-panel"),
            ], className="dual-fig corr-ba-split"),
        ], id="corr-panel-ba", className="corr-panel corr-panel-active"),

        # ── ρ/φ View ───────────────────────────────────────────────────────
        html.Div([
            # Row 1: view style + mode controls
            html.Div([
                html.Div("Style", className="corr-rho-phi-label"),
                dbc.Select(
                    id=IDs.CORR_RHO_PHI_STYLE,
                    options=[
                        {"label": "ρ/φ Grid",          "value": "grid"},
                        {"label": "Station Pairs",      "value": "pairs"},
                        {"label": "Freq Bands",         "value": "freq-bands"},
                    ],
                    value="grid",
                    size="sm",
                    style={"width": "148px"},
                ),
                html.Div(className="corr-rho-phi-sep"),
                html.Div("Show", className="corr-rho-phi-label"),
                dbc.RadioItems(
                    id=IDs.CORR_RHO_PHI_MODE,
                    options=[
                        {"label": "Both",        "value": "both"},
                        {"label": "Corrected",   "value": "corrected"},
                        {"label": "Raw",         "value": "raw"},
                    ],
                    value="both",
                    inline=True,
                    className="corr-rho-phi-radio",
                ),
                html.Div(className="corr-rho-phi-sep"),
                html.Div("Comp.", className="corr-rho-phi-label"),
                dbc.RadioItems(
                    id=IDs.CORR_COMP_SEL,
                    options=[
                        {"label": "XY",   "value": "xy"},
                        {"label": "YX",   "value": "yx"},
                        {"label": "Both", "value": "both"},
                    ],
                    value="both",
                    inline=True,
                    className="corr-rho-phi-radio corr-comp-radio",
                ),
                html.Div(style={"flex": "1"}),
                html.Div([
                    _lbl("Max stn."),
                    dbc.Input(
                        id=IDs.CORR_RHO_PHI_STN,
                        type="number", value=6,
                        min=1, max=24, step=1, size="sm",
                        debounce=True,
                        style={"width": "58px"},
                    ),
                ]),
            ], className="corr-rho-phi-bar"),
            # Row 2: line + station selectors
            html.Div([
                html.Div([
                    html.I(className="bi bi-diagram-3 me-1 corr-rho-phi-sel-icon"),
                    html.Span("Lines", className="corr-rho-phi-label"),
                ], className="d-flex align-items-center flex-shrink-0"),
                dcc.Dropdown(
                    id=IDs.CORR_LINE_SEL,
                    options=[],
                    value=None,
                    placeholder="All lines",
                    multi=True,
                    clearable=True,
                    className="corr-rho-phi-drop",
                    style={"minWidth": "140px", "flex": "1", "maxWidth": "260px"},
                ),
                html.Div(className="corr-rho-phi-sep"),
                html.Div([
                    html.I(className="bi bi-pin-map me-1 corr-rho-phi-sel-icon"),
                    html.Span("Stations", className="corr-rho-phi-label"),
                ], className="d-flex align-items-center flex-shrink-0"),
                dcc.Dropdown(
                    id=IDs.CORR_STN_SEL,
                    options=[],
                    value=None,
                    placeholder="All stations in selected lines",
                    multi=True,
                    clearable=True,
                    className="corr-rho-phi-drop",
                    style={"flex": "2", "minWidth": "180px"},
                ),
            ], className="corr-rho-phi-bar corr-rho-phi-bar-sel"),
            # Style hint row
            html.Div(
                id="corr-rho-phi-hint",
                className="corr-rho-phi-hint",
                children="Grid: ρ_a/φ per station. Switch to Station Pairs for static-shift inspection, "
                         "or Freq Bands to see removed frequency ranges.",
            ),
            # Figure area
            html.Div(
                _img(IDs.IMG_CORR_RHO_PHI),
                className="fig-area corr-rho-phi-fig", style={"flex": "1"},
            ),
        ], id="corr-panel-rho-phi", className="corr-panel",
           style={"display": "none", "flexDirection": "column"}),

        # ── Overlay ─────────────────────────────────────────────────────────
        html.Div([
            # Overlay control bar
            html.Div([
                html.Span([
                    html.I(className="bi bi-layers me-1"),
                    html.Span("Style", className="me-1"),
                ], className="corr-rho-phi-sep small text-muted"),
                dbc.Select(
                    id=IDs.CORR_OVERLAY_STYLE,
                    value="median-band",
                    options=[
                        {"label": "Median ± IQR band",  "value": "median-band"},
                        {"label": "Spaghetti + median",  "value": "spaghetti"},
                        {"label": "Per-line panels",     "value": "per-line"},
                    ],
                    size="sm",
                    style={"width": "180px", "display": "inline-block"},
                ),
                html.Span(className="flex-grow-1"),
                html.Small([
                    html.I(className="bi bi-info-circle me-1"),
                    html.Span(id="corr-overlay-hint",
                              children="Median ± IQR: statistical summary across all stations — "
                                       "gray=raw, color=corrected"),
                ], className="text-muted"),
            ], className="corr-rho-phi-bar d-flex align-items-center gap-2 px-2 py-1"),
            # Figure area
            html.Div(_img(IDs.IMG_CORR_OVERLAY),
                     className="fig-area", style={"flex": "1"}),
        ], id="corr-panel-overlay", className="corr-panel",
           style={"display": "none", "flexDirection": "column"}),

        # ── Difference ──────────────────────────────────────────────────────
        html.Div([
            html.Div(_img(IDs.IMG_CORR_DIFF),
                     className="fig-area", style={"flex": "1"}),
            html.Div(id=IDs.CORR_DIFF_STATS, className="corr-diff-stats"),
        ], id="corr-panel-diff", className="corr-panel",
           style={"display": "none", "flexDirection": "column"}),

    ], className="corr-panels-wrap")


# ── Full page layout ──────────────────────────────────────────────────────────

def layout() -> html.Div:

    sidebar = html.Div([
        _run_bar(),
        _ctrl_scroll(),
    ], className="corr-sidebar fwd-sidebar")

    view_area = html.Div([
        _ctx_bar(),
        _tab_bar(),
        _view_panels(),
    ], className="corr-view-area")

    return html.Div([
        dcc.Store(id=IDs.CORR_ACTIVE_TAB, data="ba"),
        _command_bar(
            "correction",
            "Data Correction",
            f"{_N_METHODS} corrections · {len(CATEGORIES)} categories · "
            "Non-destructive chain · per-step undo · export",
        ),
        html.Div(
            [sidebar, view_area],
            className="corr-main",
            style={"flex": "1", "display": "flex", "minHeight": "0",
                   "overflow": "hidden"},
        ),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # handled by callbacks/correction.py
