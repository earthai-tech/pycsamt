# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Geological interpretation page — 9 categories, 40+ workflows."""

from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc, html

from pycsamt.app.desktop.controllers.interp_controller import (
    CATEGORIES,
    WORKFLOW_CATALOGUE,
)
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import empty_src

PAGE_ID = "interpretation"

# ── Category metadata ─────────────────────────────────────────────────────────
# (slug, full-name, short-label, bootstrap-icon)
_INTERP_CATS = [
    ("setup", "Setup & Model", "Setup", "bi-diagram-3"),
    ("geo", "Geological", "Geology", "bi-layers"),
    ("hydro", "Hydrology", "Hydrology", "bi-droplet-half"),
    ("field", "Field Constraints", "Field", "bi-pin-map"),
    ("emdiag", "EM Diagnostics", "EM Diag", "bi-activity"),
    ("unc", "Uncertainty (MC)", "Uncert.", "bi-bar-chart-steps"),
    ("adv", "Advanced Plots", "Advanced", "bi-graph-up-arrow"),
    ("fus", "Fusion & Time-Lapse", "Fusion", "bi-clock-history"),
    ("exp", "Export", "Export", "bi-box-arrow-up"),
]

_DEFAULT_CAT_SLUG = "setup"
_CAT_BY_SLUG = {slug: full for slug, full, _, _ in _INTERP_CATS}
_SLUG_BY_CAT = {full: slug for slug, full, _, _ in _INTERP_CATS}

# ── Catalogue helpers ─────────────────────────────────────────────────────────


def _plot_opts(cat: str) -> list:
    return [
        {"label": lbl, "value": fn}
        for lbl, fn, _ in WORKFLOW_CATALOGUE.get(cat, [])
    ]


_FIRST_CAT = _CAT_BY_SLUG[_DEFAULT_CAT_SLUG]
_FIRST_OPTS = _plot_opts(_FIRST_CAT)

_CMAPS = [
    "viridis",
    "plasma",
    "inferno",
    "magma",
    "RdBu_r",
    "seismic",
    "coolwarm",
    "jet",
    "hot",
    "gray",
]
_FIGSIZE_OPTS = [
    {"label": "Auto", "value": "auto"},
    {"label": "Wide", "value": "wide"},
    {"label": "Square", "value": "square"},
    {"label": "Tall", "value": "tall"},
    {"label": "Large", "value": "large"},
]

# ── Helpers ───────────────────────────────────────────────────────────────────


def _ctrl_label(txt: str) -> html.Div:
    return html.Div(txt, className="ctrl-label")


def _hint(txt: str) -> html.Div:
    return html.Div(txt, className="fwd-feedback-mini")


# ── Sub-components ────────────────────────────────────────────────────────────


def _run_bar() -> html.Div:
    return html.Div(
        [
            dbc.Button(
                [html.I(className="bi bi-play-fill me-1"), "Generate"],
                id=IDs.BTN_INTERP_RUN,
                color="primary",
                size="sm",
                className="w-100 mb-2",
            ),
            dbc.Spinner(
                html.Div(id=IDs.INTERP_SPINNER), size="sm", color="primary"
            ),
        ],
        className="fwd-run-bar",
    )


def _cat_bar() -> html.Div:
    btns = []
    for slug, _full, short, icon in _INTERP_CATS:
        active = slug == _DEFAULT_CAT_SLUG
        btns.append(
            html.Button(
                [html.I(className=f"bi {icon}"), html.Span(short)],
                id=f"interp-cat-btn-{slug}",
                className="interp-cat-btn" + (" active" if active else ""),
                n_clicks=0,
            )
        )
    return html.Div(btns, className="interp-cat-bar")


def _controls_scroll() -> html.Div:
    _first_desc = ""
    if WORKFLOW_CATALOGUE.get(_FIRST_CAT):
        _first_desc = WORKFLOW_CATALOGUE[_FIRST_CAT][0][2]

    return html.Div(
        [
            # ── Plot / Analysis selector ──────────────────────────────────────
            html.Div(
                [
                    _ctrl_label("Plot / Analysis"),
                    dbc.Select(
                        id=IDs.INTERP_PLOT,
                        options=_FIRST_OPTS,
                        value=_FIRST_OPTS[0]["value"]
                        if _FIRST_OPTS
                        else None,
                        size="sm",
                    ),
                    html.Div(
                        id=IDs.INTERP_DESC,
                        children=_first_desc,
                        className="interp-desc-text mt-1",
                    ),
                ],
                className="ctrl-card",
            ),
            # ── Data source ───────────────────────────────────────────────────
            html.Div(
                [
                    _ctrl_label("Data Source"),
                    dbc.RadioItems(
                        id=IDs.INTERP_DATA_SRC,
                        options=[
                            {
                                "label": html.Span(
                                    [
                                        html.I(
                                            className="bi bi-database me-2",
                                            style={"color": "#89b4fa"},
                                        ),
                                        "Raw (loaded)",
                                    ],
                                    className="d-flex align-items-center",
                                ),
                                "value": "raw",
                            },
                            {
                                "label": html.Span(
                                    [
                                        html.I(
                                            className="bi bi-funnel me-2",
                                            style={"color": "#a6e3a1"},
                                        ),
                                        "Corrected",
                                    ],
                                    className="d-flex align-items-center",
                                ),
                                "value": "corrected",
                            },
                        ],
                        value="raw",
                        className="interp-src-radio",
                        inputStyle={"cursor": "pointer"},
                        labelStyle={
                            "fontSize": "12px",
                            "cursor": "pointer",
                            "padding": "3px 0",
                        },
                    ),
                    html.Div(
                        id="interp-src-hint",
                        children="Using raw loaded data.",
                        className="fwd-feedback-mini mt-1",
                    ),
                ],
                className="ctrl-card",
            ),
            # ── Display options ───────────────────────────────────────────────
            html.Div(
                [
                    _ctrl_label("Display Options"),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    _hint("Figure size"),
                                    dbc.Select(
                                        id=IDs.INTERP_FIGSIZE,
                                        options=_FIGSIZE_OPTS,
                                        value="auto",
                                        size="sm",
                                    ),
                                ],
                                width=6,
                            ),
                            dbc.Col(
                                [
                                    _hint("Colormap"),
                                    dbc.Select(
                                        id=IDs.INTERP_CMAP,
                                        options=[
                                            {"label": c, "value": c}
                                            for c in _CMAPS
                                        ],
                                        value="viridis",
                                        size="sm",
                                    ),
                                ],
                                width=6,
                            ),
                        ],
                        className="g-1 mb-2",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    _hint("Depth max (m)"),
                                    dbc.Input(
                                        id=IDs.INTERP_DEPTH_MAX,
                                        type="number",
                                        value=2000,
                                        min=10,
                                        max=50000,
                                        step=50,
                                        size="sm",
                                        debounce=True,
                                    ),
                                ],
                                width=6,
                            ),
                        ],
                        className="g-1",
                    ),
                ],
                className="ctrl-card",
            ),
            # ── Frequency band ────────────────────────────────────────────────
            html.Div(
                [
                    _ctrl_label("Frequency Band (Hz)"),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    _hint("Low (min)"),
                                    dbc.Input(
                                        id=IDs.INTERP_FREQ_LO,
                                        type="number",
                                        value=None,
                                        placeholder="all",
                                        min=1e-4,
                                        step=0.001,
                                        size="sm",
                                        debounce=True,
                                    ),
                                ]
                            ),
                            dbc.Col(
                                [
                                    _hint("High (max)"),
                                    dbc.Input(
                                        id=IDs.INTERP_FREQ_HI,
                                        type="number",
                                        value=None,
                                        placeholder="all",
                                        min=1e-4,
                                        step=0.001,
                                        size="sm",
                                        debounce=True,
                                    ),
                                ]
                            ),
                        ],
                        className="g-1",
                    ),
                ],
                className="ctrl-card",
            ),
            # ── Prerequisites (shown per category) ────────────────────────────
            html.Div(
                id="interp-prereq-card",
                style={"display": "none"},
                className="ctrl-card",
                children=[
                    html.Div("Prerequisites", className="ctrl-label"),
                    html.Div(
                        id="interp-prereq-hint",
                        className="fwd-feedback-mini mb-2",
                        style={"color": "#f9e2af", "lineHeight": "1.4"},
                    ),
                    dbc.Button(
                        [
                            html.I(className="bi bi-lightning-fill me-1"),
                            html.Span(
                                id="interp-prereq-btn-label",
                                children="Run Estimation",
                            ),
                        ],
                        id=IDs.BTN_INTERP_PREREQ,
                        color="warning",
                        size="sm",
                        className="w-100 mb-2",
                        outline=True,
                    ),
                    html.Div(
                        id=IDs.INTERP_PREREQ_OUT,
                        className="fwd-feedback-mini",
                        style={"color": "#a6e3a1", "minHeight": "14px"},
                    ),
                ],
            ),
            # ── Export ────────────────────────────────────────────────────────
            html.Div(
                [
                    _ctrl_label("Export Figure"),
                    dbc.RadioItems(
                        id=IDs.INTERP_EXPORT_FMT,
                        options=[
                            {"label": "PNG", "value": "png"},
                            {"label": "SVG", "value": "svg"},
                            {"label": "PDF", "value": "pdf"},
                            {"label": "EPS", "value": "eps"},
                        ],
                        value="png",
                        inline=True,
                        inputStyle={"cursor": "pointer"},
                        labelStyle={
                            "fontSize": "11px",
                            "cursor": "pointer",
                            "marginRight": "10px",
                        },
                        className="mb-2",
                    ),
                    dbc.Button(
                        [html.I(className="bi bi-download me-1"), "Export"],
                        id=IDs.BTN_INTERP_EXPORT,
                        color="secondary",
                        size="sm",
                        outline=True,
                        className="w-100",
                    ),
                ],
                className="ctrl-card",
            ),
        ],
        className="fwd-ctrl-scroll",
    )


def _view_panels() -> list:
    panels = []
    for slug, _full, _short, _ in _INTERP_CATS:
        visible = slug == _DEFAULT_CAT_SLUG
        panels.append(
            html.Div(
                html.Img(
                    id=f"img-interp-{slug}",
                    src=empty_src(dark=True),
                    style={
                        "width": "100%",
                        "maxHeight": "100%",
                        "objectFit": "contain",
                    },
                ),
                id=f"interp-panel-{slug}",
                style={
                    "display": "flex" if visible else "none",
                    "flexDirection": "column",
                    "flex": "1",
                    "height": "100%",
                    "minHeight": "0",
                    "padding": "8px",
                },
            )
        )
    return panels


def layout() -> html.Div:
    from pycsamt.app.web.layout import _command_bar

    total_plots = sum(len(v) for v in WORKFLOW_CATALOGUE.values())

    store = dcc.Store(id=IDs.INTERP_ACTIVE_CAT, data=_DEFAULT_CAT_SLUG)
    last_src = dcc.Store(id="interp-last-src", data="")

    sidebar = html.Div(
        [store, last_src, _run_bar(), _cat_bar(), _controls_scroll()],
        className="analysis-controls fwd-sidebar",
    )

    ctx_bar = html.Div(
        id=IDs.INTERP_INFO_STRIP,
        className="interp-info-strip",
        children=html.Span(
            "Select a workflow and click Generate.",
            style={"fontSize": "11px", "color": "#6c7086"},
        ),
    )

    view_area = html.Div(
        [
            ctx_bar,
            html.Div(
                _view_panels(),
                style={
                    "flex": "1",
                    "minHeight": "0",
                    "display": "flex",
                    "overflow": "hidden",
                },
            ),
        ],
        className="analysis-output interp-output",
        style={"display": "flex", "flexDirection": "column"},
    )

    return html.Div(
        [
            _command_bar(
                "interpret",
                "Geological Interpretation",
                f"{total_plots} workflows · {len(CATEGORIES)} categories "
                "· Hydrology · EM Diagnostics · Advanced · Fusion",
            ),
            html.Div(
                [sidebar, view_area],
                className="analysis-layout",
                style={"flex": "1"},
            ),
        ],
        style={
            "display": "flex",
            "flexDirection": "column",
            "height": "100%",
        },
    )


def register_callbacks(app) -> None:
    pass  # handled by callbacks/interp.py
