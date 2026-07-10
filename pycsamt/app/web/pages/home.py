# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Home page — station map + profile view (the original main layout)."""
from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html

from pycsamt.app.desktop.agent_registry import agent_names
from pycsamt.app.web.layout import IDs, _icon
from pycsamt.app.web.utils import (
    _empty_map,
    _empty_plotly_heatmap,
    empty_src,
)

PAGE_ID = "home"

# ── overlay options ────────────────────────────────────────────────────────────
_OVERLAY_OPTIONS = [
    {"label": lbl, "value": lbl}
    for lbl in ["Index", "Apparent Resistivity", "Phase",
                "Quality Score", "Z-strike", "N_freq"]
]

# ── helpers ────────────────────────────────────────────────────────────────────
def _station_panel():
    _hidden_cols = [{"name": c, "id": c}
                    for c in ["ID", "Line", "Latitude", "Longitude",
                               "Elevation", "N_freq"]]
    return dbc.Card(
        [
            dbc.CardHeader(
                html.Div([
                    html.Span([_icon("station-response"), "Stations"],
                              className="panel-title"),
                    html.Span(id=IDs.STATION_SUMMARY, className="stat-chip",
                              children="0 sites"),
                ], style={"display": "flex", "alignItems": "center",
                          "justifyContent": "space-between"}),
                className="py-2 px-3",
            ),
            # line filter pills
            html.Div(id=IDs.STATION_LINE_PILLS,
                     className="sta-line-pills-row", children=[]),
            # scrollable station list
            html.Div(
                html.Div(id=IDs.STATION_LIST_BODY, className="sta-list-body"),
                style={"overflowY": "auto", "flex": "1", "minHeight": "0",
                       "padding": "0"},
            ),
            # hidden DataTable keeps existing callbacks wired
            dash_table.DataTable(
                id=IDs.STATION_TABLE,
                columns=_hidden_cols,
                data=[],
                row_selectable="single",
                style_table={"display": "none", "height": "0"},
            ),
        ],
        style={"height": "100%", "display": "flex", "flexDirection": "column",
               "overflow": "hidden"},
    )


def _map_panel():
    return dbc.Card([
        dbc.CardHeader(
            dbc.Row([
                dbc.Col(html.Span([_icon("map-view"), "Station Map"],
                                  className="panel-title"), width="auto"),
                dbc.Col(
                    dbc.Select(id=IDs.HOME_MAP_OVERLAY, options=_OVERLAY_OPTIONS,
                               value="Index", size="sm",
                               style={"maxWidth": "200px"}),
                    width="auto", className="ms-auto",
                ),
            ], align="center"),
            className="py-1 px-3",
        ),
        dbc.CardBody(
            dcc.Graph(
                id=IDs.HOME_MAP_GRAPH, figure=_empty_map(dark=True),
                config={"displayModeBar": True,
                        "modeBarButtonsToRemove": ["lasso2d"],
                        "scrollZoom": True, "displaylogo": False},
                style={"height": "320px"},
                clear_on_unhover=True,
            ),
            className="p-0",
        ),
    ], className="mb-2")


def _profile_panel():
    _img_style = {"width": "100%", "maxHeight": "340px", "objectFit": "contain"}
    _bg        = {"backgroundColor": "#1e1e2e"}

    def _img_tab(img_id):
        return html.Div(html.Img(id=img_id, src=empty_src(dark=True),
                                 style=_img_style), style=_bg)

    def _graph_tab(graph_id):
        return html.Div(dcc.Graph(
            id=graph_id, figure=_empty_plotly_heatmap(dark=True),
            config={"displayModeBar": True,
                    "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                    "scrollZoom": True, "displaylogo": False,
                    "toImageButtonOptions": {"format": "png",
                                            "filename": "pycsamt_pseudosection",
                                            "scale": 2}},
            style={"height": "340px"},
        ), style=_bg)

    _freq_slider = dcc.RangeSlider(
        id=IDs.FREQ_SLIDER, min=-4, max=4, step=0.25, value=[-4, 4],
        marks=None, tooltip={"placement": "bottom", "always_visible": False},
        disabled=True, className="px-1",
    )

    from pycsamt.app.web.pages.qc_page import _qc_tab_content
    return dbc.Card([
        dbc.CardHeader([
            dbc.Row([
                dbc.Col(html.Span([_icon("profile-view"), "Profile View"],
                                  className="panel-title"),
                        width="auto", className="align-self-center"),
                dbc.Col(html.Span("log₁₀(T) range:",
                                  className="text-muted small",
                                  style={"whiteSpace": "nowrap"}),
                        width="auto", className="align-self-center pe-1"),
                dbc.Col(_freq_slider),
            ], align="center", className="g-0"),
        ], className="py-1 px-3"),
        dbc.CardBody(
            dbc.Tabs([
                dbc.Tab(_img_tab(IDs.IMG_RHO_PHI),   label="ρₐ / φ",           tab_id="tab-rho-phi"),
                dbc.Tab(_graph_tab(IDs.IMG_RHO_PS),   label="Pseudosection ρₐ", tab_id="tab-rho-ps"),
                dbc.Tab(_graph_tab(IDs.IMG_PHASE_PS), label="Pseudosection φ",  tab_id="tab-phase-ps"),
                dbc.Tab(_img_tab(IDs.IMG_TIPPER),     label="Tipper",           tab_id="tab-tipper"),
                dbc.Tab(_img_tab(IDs.IMG_PT),         label="Phase Tensor",     tab_id="tab-pt"),
                dbc.Tab(_qc_tab_content(),            label="QC",               tab_id="tab-qc"),
            ], id=IDs.PROFILE_TABS, active_tab="tab-rho-ps"),
            className="p-1",
        ),
    ])


def layout() -> html.Div:
    return html.Div([
        dbc.Row([
            dbc.Col(_station_panel(), width=2, className="pe-1"),
            dbc.Col([_map_panel(), _profile_panel()], width=7, className="px-1"),
            dbc.Col(_agents_sidebar(), width=3, className="ps-1"),
        ], className="g-0 mt-1", style={"minHeight": "calc(100vh - 120px)"}),
    ], style={"width": "100%"})


def _agents_sidebar():
    """Compact agent panel — stays in home view as a side panel."""
    _names = agent_names()
    return dbc.Card([
        dbc.CardHeader(
            html.Span([_icon("agents"), "AI Agents"], className="panel-title"),
            className="py-2 px-3",
        ),
        dbc.CardBody([
            dbc.Select(id=IDs.AGENT_SELECT,
                       options=[{"label": n, "value": n} for n in _names],
                       value=_names[0] if _names else None,
                       size="sm", className="mb-2"),
            dbc.Button([html.I(className="bi bi-play-fill me-1"), "Run Agent"],
                       id=IDs.BTN_RUN_AGENT, color="success",
                       size="sm", className="mb-2 me-2"),
            dbc.Spinner(html.Div(id=IDs.AGENT_SPINNER), size="sm", color="success"),
            html.Hr(style={"borderColor": "#313244", "margin": "6px 0"}),
            html.P([html.I(className="bi bi-terminal-fill me-1 text-muted"), "Log"],
                   className="mb-1 text-muted small"),
            html.Pre("", id=IDs.AGENT_LOG, style={
                "backgroundColor": "#181825", "color": "#cdd6f4",
                "fontSize": "11px", "height": "140px",
                "overflowY": "auto", "padding": "6px", "borderRadius": "4px",
                "whiteSpace": "pre-wrap",
            }),
            html.P([html.I(className="bi bi-image me-1 text-muted"), "Result"],
                   className="mb-1 mt-2 text-muted small"),
            html.Img(id=IDs.AGENT_RESULT, src=empty_src(dark=True), style={
                "width": "100%", "maxHeight": "200px", "objectFit": "contain",
                "backgroundColor": "#181825",
            }),
            html.Div("", id=IDs.AGENT_SUMMARY,
                     className="mt-1 text-muted small", style={"fontSize": "11px"}),
        ], className="p-2", style={"overflowY": "auto", "maxHeight": "680px"}),
    ], style={"height": "100%"})


def register_callbacks(app) -> None:
    pass  # home callbacks live in the main callbacks/ package
