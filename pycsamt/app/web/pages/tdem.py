# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""TDEM (time-domain EM) analysis page."""
from __future__ import annotations

from dash import dcc, html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.tdem_controller import TDEM_GROUPS
from pycsamt.app.web.layout import IDs, _icon
from pycsamt.app.web.utils import empty_src

PAGE_ID = "tdem"


def layout() -> html.Div:
    group_opts  = [{"label": g, "value": g} for g, _ in TDEM_GROUPS]
    first_plots = [{"label": lbl, "value": cls} for lbl, cls, *_ in TDEM_GROUPS[0][1]]

    controls = html.Div([
        html.Div([
            html.Div("TEM Folder", className="ctrl-label"),
            dbc.InputGroup([
                dbc.Input(id=IDs.TDEM_FOLDER, placeholder="/path/to/TEM/folder",
                          type="text", size="sm"),
                dbc.Button([_icon("open", size=12, cls=""), " Load"],
                           id=IDs.BTN_TDEM_LOAD, color="secondary",
                           size="sm", className="btn-nav-icon"),
            ], size="sm"),
            html.Div(id=IDs.TDEM_INFO,
                     className="mt-1 text-muted small",
                     style={"fontSize": "10px"}),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Category", className="ctrl-label"),
            dbc.Select(id=IDs.TDEM_CAT, options=group_opts,
                       value=group_opts[0]["value"], size="sm"),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Plot", className="ctrl-label"),
            dbc.Select(id=IDs.TDEM_PLOT, options=first_plots,
                       value=first_plots[0]["value"], size="sm"),
        ], className="ctrl-card"),
        html.Div([
            dbc.Button([html.I(className="bi bi-play-fill me-1"), "Generate"],
                       id=IDs.BTN_TDEM_RUN, color="info", size="sm",
                       className="w-100 mb-2"),
            dbc.Spinner(html.Div(id=IDs.TDEM_SPINNER), size="sm", color="info"),
        ], className="ctrl-card"),
    ], className="analysis-controls")

    output = html.Div([
        html.Div([
            html.Img(id=IDs.IMG_TDEM, src=empty_src(dark=True),
                     style={"width": "100%", "maxHeight": "100%",
                            "objectFit": "contain"}),
        ], className="fig-area", style={"minHeight": "500px"}),
    ], className="analysis-output")

    from pycsamt.app.web.layout import _command_bar
    return html.Div([
        _command_bar("tdem", "TDEM Analysis",
                     "Time-domain EM sounding · 16 plots · 4 categories"),
        html.Div([controls, output], className="analysis-layout", style={"flex": "1"}),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # handled by callbacks/tdem.py
