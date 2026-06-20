# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""TDEM analysis page — folder browser, fixed tab bar, per-category figures,
rich parameter controls, and high-resolution export."""
from __future__ import annotations

from pathlib import Path

from dash import dcc, html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.tdem_controller import (
    DECAY_PLOTS, SECTION_PLOTS, MAP_PLOTS, DASHBOARD_PLOTS,
)
from pycsamt.app.web.layout import IDs, _icon
from pycsamt.app.web.utils import empty_src

PAGE_ID = "tdem"

# ── Locate bundled JIANGSU sample data ────────────────────────────────────────
def _sample_path() -> str:
    pkg_root = Path(__file__).resolve().parents[4]
    jiangsu  = pkg_root / "data" / "TEMAVG" / "JIANGSU"
    return str(jiangsu) if jiangsu.is_dir() else ""

_SAMPLE = _sample_path()

# ── Tab definitions ────────────────────────────────────────────────────────────
_TDEM_TABS = [
    ("decay",     "Decay / Rho",    "bi-activity",  DECAY_PLOTS),
    ("section",   "Survey Section", "bi-layers",    SECTION_PLOTS),
    ("map",       "Map & Overview", "bi-map",       MAP_PLOTS),
    ("dashboard", "Dashboard",      "bi-grid-3x3",  DASHBOARD_PLOTS),
]
_TDEM_DEFAULT_TAB = "decay"

_TDEM_TAB_GROUP = {
    "decay":     "Decay / Rho",
    "section":   "Survey Section",
    "map":       "Map & Overview",
    "dashboard": "Dashboard",
}

# First plot options for the default tab
_FIRST_PLOTS = [{"label": lbl, "value": cls} for lbl, cls, *_ in DECAY_PLOTS]
_FIRST_CLS   = _FIRST_PLOTS[0]["value"]

_FIGSIZE_OPTIONS = [
    {"label": "Compact  8×4",  "value": "8x4"},
    {"label": "Standard 10×5", "value": "10x5"},
    {"label": "Wide     13×5", "value": "13x5"},
    {"label": "Tall     10×7", "value": "10x7"},
]

_CMAP_OPTIONS = [
    {"label": "RdYlBu_r (default)", "value": "RdYlBu_r"},
    {"label": "viridis",            "value": "viridis"},
    {"label": "terrain",            "value": "terrain"},
    {"label": "coolwarm",           "value": "coolwarm"},
    {"label": "jet",                "value": "jet"},
]

def _ctrl_label(text: str) -> html.Div:
    return html.Div(text, className="ctrl-label")


def layout() -> html.Div:
    from pycsamt.app.web.layout import _command_bar

    _img_style = {"width": "100%", "height": "100%", "objectFit": "contain"}

    # ── Left sidebar ──────────────────────────────────────────────────────────
    controls = html.Div(
        [
            # ── Run bar — FIXED at top ─────────────────────────────────────
            html.Div(
                [
                    dbc.Button(
                        [html.I(className="bi bi-play-fill me-1"), "Generate"],
                        id=IDs.BTN_TDEM_RUN,
                        color="primary", size="sm",
                        className="w-100 mb-1",
                        n_clicks=0,
                    ),
                    dbc.Spinner(
                        html.Div(id=IDs.TDEM_SPINNER), size="sm", color="primary",
                    ),
                ],
                className="fwd-run-bar",
            ),

            # ── Scrollable controls ────────────────────────────────────────
            html.Div(
                [
                    # 1. Folder loader
                    html.Div(
                        [
                            _ctrl_label("TEM Data Folder"),
                            dbc.Button(
                                [html.I(className="bi bi-folder2-open me-2"), "Browse folder…"],
                                id=IDs.BTN_TDEM_BROWSE,
                                color="secondary", size="sm",
                                className="w-100 mb-2 tdem-browse-btn",
                                n_clicks=0,
                            ),
                            dbc.Input(
                                id=IDs.TDEM_FOLDER,
                                placeholder="Folder path will appear here",
                                type="text", size="sm",
                                className="mb-2 tdem-path-input",
                            ),
                            dbc.Button(
                                [html.I(className="bi bi-box-arrow-in-down me-1"), "Load"],
                                id=IDs.BTN_TDEM_LOAD,
                                color="outline-teal", size="sm",
                                className="w-100 btn-tdem-load",
                                n_clicks=0,
                            ),
                            html.Div(id=IDs.TDEM_INFO, className="mt-2 tdem-info-text"),
                        ],
                        className="ctrl-card",
                    ),

                    # 2. Demo data card
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Span(
                                        [html.I(className="bi bi-box-seam me-1"), "Demo dataset"],
                                        className="tdem-demo-badge",
                                    ),
                                ],
                                className="tdem-demo-header",
                            ),
                            html.Div(
                                "JIANGSU · 9 stations · .AVG / .Z files",
                                className="tdem-demo-desc",
                            ),
                            dbc.Button(
                                [html.I(className="bi bi-play-circle me-1"), "Load demo"],
                                id=IDs.BTN_TDEM_SAMPLE,
                                color="outline-secondary", size="sm",
                                className="w-100 mt-2",
                                disabled=not bool(_SAMPLE),
                                title=(
                                    _SAMPLE if _SAMPLE
                                    else "Sample data not bundled in this installation"
                                ),
                                n_clicks=0,
                            ),
                            html.Div(
                                "For testing only — not your survey data",
                                className="tdem-demo-note",
                            ),
                        ],
                        className="ctrl-card tdem-demo-card",
                        style={"display": "block"} if _SAMPLE else {"display": "none"},
                    ),

                    # 3. Plot selector
                    html.Div(
                        [
                            _ctrl_label("Plot"),
                            dbc.Select(
                                id=IDs.TDEM_PLOT,
                                options=_FIRST_PLOTS,
                                value=_FIRST_CLS,
                                size="sm",
                            ),
                        ],
                        className="ctrl-card",
                    ),

                    # 4. Parameters
                    html.Div(
                        [
                            _ctrl_label("Parameters"),
                            html.Div(
                                [
                                    html.Div("Figure size", className="adv-param-label"),
                                    dbc.Select(
                                        id=IDs.TDEM_FIGSIZE,
                                        options=_FIGSIZE_OPTIONS,
                                        value="10x5",
                                        size="sm",
                                    ),
                                ],
                                className="mb-2",
                            ),
                            html.Div(
                                [
                                    html.Hr(className="adv-divider"),
                                    html.Div("Sounding index", className="adv-param-label"),
                                    dbc.Input(
                                        id=IDs.TDEM_SOUNDING_IDX,
                                        type="number", value=0, min=0, step=1,
                                        size="sm", className="adv-param-input",
                                    ),
                                    html.Div(
                                        "0 = first sounding · blank = all",
                                        className="adv-param-hint",
                                    ),
                                ],
                                id="tdem-sounding-section",
                                style={"display": "none"},
                            ),
                            html.Div(
                                [
                                    html.Hr(className="adv-divider"),
                                    html.Div("Colormap", className="adv-param-label"),
                                    dbc.Select(
                                        id=IDs.TDEM_CMAP,
                                        options=_CMAP_OPTIONS,
                                        value="RdYlBu_r",
                                        size="sm",
                                    ),
                                ],
                                id="tdem-cmap-section",
                                style={"display": "none"},
                            ),
                            html.Div(
                                [
                                    html.Hr(className="adv-divider"),
                                    html.Div("Gate index", className="adv-param-label"),
                                    dbc.Input(
                                        id=IDs.TDEM_GATE_IDX,
                                        type="number", value=0, min=0, step=1,
                                        size="sm", className="adv-param-input",
                                    ),
                                    html.Div(
                                        "0 = first gate · blank = all",
                                        className="adv-param-hint",
                                    ),
                                ],
                                id="tdem-gate-section",
                                style={"display": "none"},
                            ),
                        ],
                        className="ctrl-card",
                    ),
                ],
                className="fwd-ctrl-scroll",
            ),
        ],
        className="analysis-controls fwd-sidebar",
        style={"width": "240px", "minWidth": "240px"},
    )

    # ── Fixed tab bar ─────────────────────────────────────────────────────────
    tab_bar = html.Div(
        [
            html.Button(
                [html.I(className=f"bi {icon} me-1"), label],
                id=f"tdem-tab-{tid}",
                className=(
                    f"prof-tab-btn{' active' if tid == _TDEM_DEFAULT_TAB else ''}"
                ),
                n_clicks=0,
                title=_TDEM_TAB_GROUP.get(tid, label),
            )
            for tid, label, icon, _ in _TDEM_TABS
        ],
        className="prof-tab-bar",
    )

    # ── Context info bar ─────────────────────────────────────────────────────
    ctx_bar = html.Div(
        [
            html.Span(
                [html.I(className="bi bi-activity me-1"), "Decay / Rho"],
                className="adv-info-group",
            ),
            html.Span("·", className="adv-info-sep"),
            html.Span(
                "Load a TDEM folder then click Generate",
                style={"color": "var(--sub0)", "fontSize": "10.5px"},
            ),
        ],
        id=IDs.TDEM_CTX_BAR,
        className="adv-info-bar",
    )

    # ── Per-tab figure panels ─────────────────────────────────────────────────
    def _fig_panel(tab_id: str, img_id: str, visible: bool = False) -> html.Div:
        return html.Div(
            html.Div(
                html.Img(
                    id=img_id,
                    src=empty_src(),
                    style=_img_style,
                ),
                className="fig-img-wrap profile-page-fig",
            ),
            id=f"tdem-panel-{tab_id}",
            className="prof-panel",
            style={"display": "flex"} if visible else {"display": "none"},
        )

    panels = html.Div(
        [
            _fig_panel("decay",     IDs.IMG_TDEM_DECAY,   visible=True),
            _fig_panel("section",   IDs.IMG_TDEM_SECTION),
            _fig_panel("map",       IDs.IMG_TDEM_MAP),
            _fig_panel("dashboard", IDs.IMG_TDEM_DASH),
        ],
        className="prof-view-panel",
    )

    # ── View area ─────────────────────────────────────────────────────────────
    view_area = html.Div(
        [
            dcc.Store(id=IDs.TDEM_ACTIVE_TAB, data=_TDEM_DEFAULT_TAB),
            tab_bar,
            ctx_bar,
            panels,
        ],
        className="adv-view-area",
    )

    return html.Div(
        [
            _command_bar(
                "tdem", "TDEM Analysis",
                "Time-domain EM · Decay, Section, Map & Dashboard",
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
    pass  # handled by callbacks/tdem.py
