# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""QC page v2 — quality-control dashboard.

Layout mirrors the Correction page design:
  · Left sidebar  — run bar (top), group pills, plot list, params, figsize
  · Right view    — ctx-bar (Lines/Stations/Comp), tab bar, panels
"""
from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc, html

from pycsamt.app.desktop.controllers.qc_controller import (
    ALL_GROUPS,
)
from pycsamt.app.web.layout import IDs, _command_bar
from pycsamt.app.web.utils import empty_src

PAGE_ID = "qc"

_GROUPS = [(gname, plots) for gname, plots in ALL_GROUPS]

_FIGSIZE_OPTS = [
    {"label": "Compact", "value": "compact"},
    {"label": "Wide",    "value": "wide"},
    {"label": "Tall",    "value": "tall"},
    {"label": "Pub",     "value": "pub"},
]

_METRIC_OPTS = [
    {"label": "coherence",   "value": "coherence"},
    {"label": "coverage",    "value": "coverage"},
    {"label": "snr",         "value": "snr"},
    {"label": "reliability", "value": "reliability"},
]

_METHOD_OPTS = [
    {"label": "iqr",      "value": "iqr"},
    {"label": "variance", "value": "variance"},
    {"label": "zscore",   "value": "zscore"},
]


def _plot_opts(group_name: str) -> list[dict]:
    for g, plots in _GROUPS:
        if g == group_name:
            return [{"label": lbl, "value": fn} for lbl, fn, _ in plots]
    return []


def layout() -> html.Div:
    first_group = _GROUPS[0][0]
    first_plots = _plot_opts(first_group)

    # ── Group pill rows (2 pills per row) ─────────────────────────────────────
    pill_rows = []
    for i in range(0, len(_GROUPS), 2):
        row_btns = []
        for gname, _plots in _GROUPS[i : i + 2]:
            is_active = gname == first_group
            row_btns.append(
                html.Button(
                    gname,
                    id={"type": "qc-grp-btn", "index": gname},
                    className="qc-grp-btn" + (" active" if is_active else ""),
                    n_clicks=0,
                )
            )
        pill_rows.append(html.Div(row_btns, className="qc-grp-row"))

    # ── Left sidebar ──────────────────────────────────────────────────────────
    sidebar = html.Div(
        [
            # ── Run bar at sidebar top (mirrors fwd-run-bar / corr-run-bar) ──
            html.Div(
                [
                    # Row: Refresh button + Auto toggle side by side
                    html.Div(
                        [
                            dbc.Button(
                                [html.I(className="bi bi-arrow-repeat me-1"), "Refresh"],
                                id=IDs.BTN_QC_RUN,
                                color="info",
                                size="sm",
                                className="flex-grow-1",
                            ),
                            dbc.Checklist(
                                id=IDs.QC_AUTO_TOGGLE,
                                options=[{"label": "Auto", "value": "auto"}],
                                value=["auto"],
                                switch=True,
                                inline=True,
                                className="qc-auto-toggle mb-0",
                            ),
                        ],
                        className="d-flex align-items-center gap-2 mb-1",
                    ),
                    # Spinner + feedback below
                    dbc.Spinner(
                        html.Div(id=IDs.QC_SPINNER),
                        size="sm",
                        color="info",
                    ),
                    html.Div(id=IDs.QC_FEEDBACK, className="qc-feedback-sidebar"),
                ],
                className="qc-run-bar",
            ),

            # ── Scrollable area (everything below run bar scrolls) ───────────
            html.Div(
                [
                    # Group pills
                    html.Div(
                        [
                            html.Div("Group", className="ctrl-label"),
                            html.Div(pill_rows, className="qc-grp-pills"),
                        ],
                        className="ctrl-card",
                    ),

                    # Plot RadioItems list
                    html.Div(
                        [
                            html.Div("Plot", className="ctrl-label"),
                            html.Div(
                                dbc.RadioItems(
                                    id=IDs.QC_PLOT_SELECT,
                                    options=first_plots,
                                    value=first_plots[0]["value"] if first_plots else None,
                                    className="qc-plot-list",
                                ),
                                className="qc-plot-list-wrap",
                            ),
                        ],
                        className="ctrl-card qc-plot-card",
                    ),

                    # Contextual parameters panel
                    html.Div(
                        [
                            html.Div("Parameters", className="ctrl-label"),
                            html.Div(
                                [
                                    html.Span("Metric", className="qc-param-label"),
                                    dbc.Select(
                                        id=IDs.QC_PARAM_METRIC,
                                        options=_METRIC_OPTS,
                                        value="coherence",
                                        size="sm",
                                    ),
                                ],
                                className="qc-param-row mb-1",
                            ),
                            html.Div(
                                [
                                    html.Span("Method", className="qc-param-label"),
                                    dbc.Select(
                                        id=IDs.QC_PARAM_METHOD,
                                        options=_METHOD_OPTS,
                                        value="iqr",
                                        size="sm",
                                    ),
                                ],
                                className="qc-param-row mb-1",
                            ),
                            html.Div(
                                [
                                    html.Span("Threshold", className="qc-param-label"),
                                    dbc.Input(
                                        id=IDs.QC_PARAM_THRESH,
                                        type="number",
                                        value=0.3,
                                        min=0,
                                        max=1,
                                        step=0.05,
                                        size="sm",
                                    ),
                                ],
                                className="qc-param-row mb-1",
                            ),
                            html.Div(
                                [
                                    html.Span("N-max stns", className="qc-param-label"),
                                    dbc.Input(
                                        id=IDs.QC_PARAM_NMAX,
                                        type="number",
                                        value=20,
                                        min=1,
                                        max=500,
                                        step=1,
                                        size="sm",
                                    ),
                                ],
                                className="qc-param-row mb-1",
                            ),
                            html.Div(
                                [
                                    html.Span("Mains Hz", className="qc-param-label"),
                                    dbc.RadioItems(
                                        id=IDs.QC_PARAM_MAINS,
                                        options=[
                                            {"label": "50 Hz", "value": 50},
                                            {"label": "60 Hz", "value": 60},
                                        ],
                                        value=50,
                                        inline=True,
                                        className="qc-comp-radio",
                                    ),
                                ],
                                className="qc-param-row",
                            ),
                        ],
                        id=IDs.QC_PARAMS_PANEL,
                        className="ctrl-card qc-params-card",
                    ),

                    # Figure size
                    html.Div(
                        [
                            html.Div("Figure size", className="ctrl-label"),
                            dbc.RadioItems(
                                id=IDs.QC_FIGSIZE,
                                options=_FIGSIZE_OPTS,
                                value="compact",
                                inline=True,
                                className="qc-figsize-radio",
                            ),
                        ],
                        className="ctrl-card",
                    ),

                    # Hidden stores live here so they scroll with content (non-visible)
                    dcc.Store(id=IDs.QC_GROUP_SEL,  data=first_group),
                    dcc.Store(id=IDs.QC_ACTIVE_TAB, data="plot"),
                    html.Div(
                        dbc.Select(
                            id=IDs.QC_GROUP_SELECT,
                            options=[{"label": g, "value": g} for g, _ in _GROUPS],
                            value=first_group,
                        ),
                        style={"display": "none"},
                    ),
                ],
                className="qc-ctrl-scroll",
            ),
        ],
        className="qc-sidebar",
    )

    # ── Right view area ───────────────────────────────────────────────────────
    view_area = html.Div(
        [
            # ── Context bar: Lines / Stations / Component ──────────────────────
            html.Div(
                [
                    html.Span("Lines", className="qc-ctx-label"),
                    dcc.Dropdown(
                        id=IDs.QC_LINE_SEL,
                        placeholder="All lines…",
                        multi=True,
                        clearable=True,
                        className="qc-ctx-drop",
                        style={"minWidth": "130px", "maxWidth": "220px"},
                    ),
                    html.Span("Stations", className="qc-ctx-label ms-2"),
                    dcc.Dropdown(
                        id=IDs.QC_STN_SEL,
                        placeholder="All stations…",
                        multi=True,
                        clearable=True,
                        className="qc-ctx-drop",
                        style={"minWidth": "130px", "maxWidth": "220px"},
                    ),
                    html.Div(
                        [
                            html.Span("Comp", className="qc-ctx-label ms-2"),
                            dbc.RadioItems(
                                id=IDs.QC_COMP_SEL,
                                options=[
                                    {"label": "XY",   "value": "xy"},
                                    {"label": "YX",   "value": "yx"},
                                    {"label": "Both", "value": "both"},
                                ],
                                value="both",
                                inline=True,
                                className="qc-comp-radio",
                            ),
                        ],
                        className="d-flex align-items-center gap-1 flex-shrink-0",
                    ),
                ],
                className="qc-ctx-bar",
            ),

            # ── Tab bar ────────────────────────────────────────────────────────
            html.Div(
                [
                    html.Button(
                        [html.I(className="bi bi-bar-chart me-1"), "Plot"],
                        id="qc-tab-btn-plot",
                        className="qc-tab-btn active",
                        n_clicks=0,
                    ),
                    html.Button(
                        [html.I(className="bi bi-grid-3x3-gap me-1"), "Overview"],
                        id="qc-tab-btn-overview",
                        className="qc-tab-btn",
                        n_clicks=0,
                    ),
                ],
                className="qc-tab-bar",
            ),

            # ── Panels ────────────────────────────────────────────────────────
            html.Div(
                [
                    html.Div(
                        html.Img(
                            id=IDs.IMG_QC_PLOT,
                            src=empty_src(dark=True),
                            style={
                                "width": "100%",
                                "height": "100%",
                                "objectFit": "contain",
                            },
                        ),
                        id="qc-panel-plot",
                        className="qc-panel qc-panel-active",
                    ),
                    html.Div(
                        html.Img(
                            id=IDs.IMG_QC_OVERVIEW,
                            src=empty_src(dark=True),
                            style={
                                "width": "100%",
                                "height": "100%",
                                "objectFit": "contain",
                            },
                        ),
                        id="qc-panel-overview",
                        className="qc-panel",
                        style={"display": "none"},
                    ),
                ],
                className="qc-panels-wrap",
            ),
        ],
        className="qc-view-area",
    )

    return html.Div(
        [
            _command_bar(
                "qc",
                "Quality Control",
                f"{sum(len(p) for _, p in _GROUPS)} plots · {len(_GROUPS)} groups",
            ),
            html.Div(
                [sidebar, view_area],
                className="qc-body",
            ),
        ],
        style={"display": "flex", "flexDirection": "column", "height": "100%"},
    )


def register_callbacks(app) -> None:
    pass  # all callbacks in pycsamt.app.web.callbacks.qc
