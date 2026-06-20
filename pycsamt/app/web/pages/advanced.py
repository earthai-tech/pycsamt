# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Advanced emtools plots page — fixed tab bar + per-group persistent figures."""
from __future__ import annotations

from dash import dcc, html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.advanced_controller import (
    ADVANCED_GROUPS,
    ADVANCED_PLOT_DESCRIPTIONS,
    STRIKE_PLOTS,
    PHASE_TENSOR_PLOTS,
    INDUCTION_PLOTS,
    IMPEDANCE_PLOTS,
    DEPTH_PLOTS,
    SURVEY_PLOTS,
)
from pycsamt.app.web.layout import IDs, _icon
from pycsamt.app.web.utils import empty_src

PAGE_ID = "advanced"

# ── Tab definitions (one per non-empty plot group) ────────────────────────────
_ADV_TABS = [
    ("strike",    "Strike",        "bi-compass",           STRIKE_PLOTS),
    ("pt",        "Phase Tensor",  "bi-hexagon",           PHASE_TENSOR_PLOTS),
    ("induction", "Induction",     "bi-arrows-expand",     INDUCTION_PLOTS),
    ("impedance", "Impedance/Z",   "bi-lightning-charge",  IMPEDANCE_PLOTS),
    ("depth",     "Depth",         "bi-layers",            DEPTH_PLOTS),
    ("survey",    "Survey Tools",  "bi-clipboard-data",    SURVEY_PLOTS),
]
_ADV_DEFAULT_TAB = "strike"

# Map tab-id → group label (for info bar)
_ADV_TAB_TO_GROUP = {
    "strike":    "Strike Analysis",
    "pt":        "Phase Tensor",
    "induction": "Induction / Tipper",
    "impedance": "Impedance / Z",
    "depth":     "Depth Imaging",
    "survey":    "Survey Tools",
}

# Map tab-id → per-group image component ID
_ADV_TAB_IMG = {
    "strike":    IDs.IMG_ADV_STRIKE,
    "pt":        IDs.IMG_ADV_PT,
    "induction": IDs.IMG_ADV_INDUCTION,
    "impedance": IDs.IMG_ADV_IMPEDANCE,
    "depth":     IDs.IMG_ADV_DEPTH,
    "survey":    IDs.IMG_ADV_SURVEY,
}

# Plots that expose a period-index parameter (ip=)
_PERIOD_FNS = {
    "plot_phase_tensor_map", "plot_induction_arrows", "plot_induction_map",
    "plot_xyyx_crossover_map", "plot_dim_map", "plot_dim_occupancy_area",
    "plot_anisotropy", "plot_ellipticity_psection", "plot_depth_section",
    "plot_apparent_depth_psection", "plot_gradient_section",
    "plot_induction_section", "plot_theta_vs_period",
    "plot_theta_stability_stripe", "plot_strike_ribbon", "plot_strike_mapsticks",
}

_FIGSIZE_OPTIONS = [
    {"label": "Compact  8×4",  "value": "8x4"},
    {"label": "Standard 10×6", "value": "10x6"},
    {"label": "Wide     14×6", "value": "14x6"},
    {"label": "Tall     10×8", "value": "10x8"},
]

# First plot of the default group
_FIRST_PLOTS = [{"label": l, "value": f} for l, f, _ in STRIKE_PLOTS]
_FIRST_FN    = _FIRST_PLOTS[0]["value"]


def _ctrl_label(text: str) -> html.Div:
    return html.Div(text, className="ctrl-label")


def layout() -> html.Div:
    _img_style = {"width": "100%", "height": "100%", "objectFit": "contain"}

    # ── Left sidebar ──────────────────────────────────────────────────────────
    controls = html.Div(
        [
            # ── Run bar — FIXED at top ─────────────────────────────────────
            html.Div(
                [
                    dbc.Button(
                        [html.I(className="bi bi-play-fill me-1"), "Generate"],
                        id=IDs.BTN_ADV_RUN,
                        color="primary", size="sm",
                        className="w-100 mb-1",
                        n_clicks=0,
                    ),
                    dbc.Spinner(
                        html.Div(id=IDs.ADV_SPINNER), size="sm", color="primary",
                    ),
                    html.Div(id=IDs.ADV_DATA_BAR, className="adv-data-bar mt-1"),
                ],
                className="fwd-run-bar",
            ),

            # ── Scrollable controls ────────────────────────────────────────
            html.Div(
                [
                    # Plot selector
                    html.Div(
                        [
                            _ctrl_label("Plot"),
                            dbc.Select(
                                id=IDs.ADV_PLOT,
                                options=_FIRST_PLOTS,
                                value=_FIRST_FN,
                                size="sm",
                            ),
                        ],
                        className="ctrl-card",
                    ),

                    # Active Lines filter
                    html.Div(
                        [
                            _ctrl_label("Active Lines"),
                            dcc.Dropdown(
                                id=IDs.ADV_LINE_FILTER,
                                options=[],
                                value=None,
                                multi=True,
                                placeholder="All active lines",
                                className="adv-filter-drop",
                                optionHeight=26,
                            ),
                        ],
                        className="ctrl-card",
                    ),

                    # Station filter
                    html.Div(
                        [
                            html.Div(
                                [
                                    _ctrl_label("Stations"),
                                    html.Div(
                                        [
                                            html.Button(
                                                "All", id=IDs.ADV_STN_ALL,
                                                className="agents-mini-btn", n_clicks=0,
                                            ),
                                            html.Button(
                                                "None", id=IDs.ADV_STN_NONE,
                                                className="agents-mini-btn", n_clicks=0,
                                            ),
                                        ],
                                        className="agents-filter-btn-row",
                                    ),
                                ],
                                className="agents-filter-header-row",
                            ),
                            dcc.Dropdown(
                                id=IDs.ADV_STN_FILTER,
                                options=[],
                                value=None,
                                multi=True,
                                placeholder="All stations from active lines",
                                className="adv-filter-drop",
                                optionHeight=26,
                            ),
                        ],
                        className="ctrl-card",
                    ),

                    # Parameters
                    html.Div(
                        [
                            _ctrl_label("Parameters"),
                            html.Div(
                                [
                                    html.Div("Period index (ip)", className="adv-param-label"),
                                    dbc.Input(
                                        id="adv-period-idx",
                                        type="number", value=0, min=0, step=1,
                                        size="sm", className="adv-param-input",
                                    ),
                                    html.Div("0 = highest frequency",
                                             className="adv-param-hint"),
                                ],
                                id="adv-period-section",
                                style={"display": "none"},
                            ),
                            html.Div(
                                [
                                    html.Div("Figure size", className="adv-param-label"),
                                    dbc.Select(
                                        id="adv-figsize",
                                        options=_FIGSIZE_OPTIONS,
                                        value="10x6",
                                        size="sm",
                                    ),
                                ],
                                className="mb-2",
                            ),
                            html.Div(
                                [
                                    html.Hr(className="adv-divider"),
                                    html.Div("Dictionary model (ATOM)",
                                             className="adv-param-label"),
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Div("n_atoms", className="adv-atom-label"),
                                                    dbc.Input(
                                                        id="adv-atom-atoms",
                                                        type="number",
                                                        value=6, min=2, max=64, step=1,
                                                        size="sm", className="adv-param-input",
                                                    ),
                                                ],
                                                className="adv-atom-row",
                                            ),
                                            html.Div(
                                                [
                                                    html.Div("n_iter", className="adv-atom-label"),
                                                    dbc.Input(
                                                        id="adv-atom-iter",
                                                        type="number",
                                                        value=40, min=5, max=500, step=5,
                                                        size="sm", className="adv-param-input",
                                                    ),
                                                ],
                                                className="adv-atom-row",
                                            ),
                                        ],
                                        className="mb-2",
                                    ),
                                    dbc.Button(
                                        [html.I(className="bi bi-cpu me-1"), "Train Model"],
                                        id=IDs.ADV_BTN_TRAIN,
                                        color="info", size="sm", outline=True,
                                        className="w-100 mb-1",
                                        n_clicks=0,
                                    ),
                                    dbc.Spinner(
                                        html.Div(id=IDs.ADV_TRAIN_SPINNER), size="sm",
                                        color="info",
                                    ),
                                    html.Div(id=IDs.ADV_TRAIN_STATUS,
                                             className="adv-train-status"),
                                ],
                                id="adv-atom-section",
                                style={"display": "none"},
                            ),
                        ],
                        className="ctrl-card",
                    ),

                    # Description
                    html.Div(
                        [
                            _ctrl_label("Description"),
                            html.Div(
                                ADVANCED_PLOT_DESCRIPTIONS.get(_FIRST_FN, ""),
                                id=IDs.ADV_DESCRIPTION,
                                className="adv-description",
                            ),
                        ],
                        className="ctrl-card",
                    ),
                ],
                className="fwd-ctrl-scroll",
            ),
        ],
        className="analysis-controls fwd-sidebar",
        style={"width": "260px", "minWidth": "260px"},
    )

    # ── Fixed tab bar ─────────────────────────────────────────────────────────
    tab_bar = html.Div(
        [
            html.Button(
                [html.I(className=f"bi {icon} me-1"), label],
                id=f"adv-tab-{tid}",
                className=(
                    f"prof-tab-btn{' active' if tid == _ADV_DEFAULT_TAB else ''}"
                ),
                n_clicks=0,
                title=_ADV_TAB_TO_GROUP.get(tid, label),
            )
            for tid, label, icon, _ in _ADV_TABS
        ],
        className="prof-tab-bar",
    )

    # ── Context info bar ──────────────────────────────────────────────────────
    info_bar = html.Div(
        [
            html.Span(
                [html.I(className="bi bi-compass me-1"), "Strike Analysis"],
                className="adv-info-group",
            ),
            html.Span("·", className="adv-info-sep"),
            html.Span("Select a plot and click Generate",
                      style={"color": "var(--sub0)", "fontSize": "10.5px"}),
        ],
        id=IDs.ADV_INFO_BAR,
        className="adv-info-bar",
    )

    # ── Per-group figure panels ───────────────────────────────────────────────
    def _fig_panel(tab_id: str, img_id: str, visible: bool = False) -> html.Div:
        return html.Div(
            html.Div(
                html.Img(
                    id=img_id,
                    src=empty_src(dark=True),
                    style=_img_style,
                ),
                className="fig-img-wrap profile-page-fig",
            ),
            id=f"adv-panel-{tab_id}",
            className="prof-panel",
            style={"display": "flex"} if visible else {"display": "none"},
        )

    panels = html.Div(
        [
            _fig_panel("strike",    IDs.IMG_ADV_STRIKE,    visible=True),
            _fig_panel("pt",        IDs.IMG_ADV_PT),
            _fig_panel("induction", IDs.IMG_ADV_INDUCTION),
            _fig_panel("impedance", IDs.IMG_ADV_IMPEDANCE),
            _fig_panel("depth",     IDs.IMG_ADV_DEPTH),
            _fig_panel("survey",    IDs.IMG_ADV_SURVEY),
        ],
        className="prof-view-panel",
    )

    # ── View area (tab bar + info bar + panels) ───────────────────────────────
    view_area = html.Div(
        [
            dcc.Store(id=IDs.ADV_ACTIVE_TAB, data=_ADV_DEFAULT_TAB),
            tab_bar,
            info_bar,
            panels,
        ],
        className="adv-view-area",
    )

    from pycsamt.app.web.layout import _command_bar
    return html.Div(
        [
            _command_bar(
                "advanced-tools",
                "Advanced Analysis",
                "Strike · Phase Tensor · Induction · Impedance · Depth",
            ),
            html.Div(
                [controls, view_area],
                className="analysis-layout",
                style={"flex": "1"},
            ),
        ],
        style={"display": "flex", "flexDirection": "column", "height": "100%"},
    )


def register_callbacks(app) -> None:
    pass  # handled by callbacks/advanced.py
