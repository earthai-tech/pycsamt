# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""QC page — full-page quality-control dashboard."""
from __future__ import annotations

from dash import html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.qc_controller import ALL_GROUPS
from pycsamt.app.web.layout import IDs, _icon
from pycsamt.app.web.utils import empty_src

PAGE_ID = "qc"

_GROUP_MAP = {gname: plots for gname, plots in ALL_GROUPS}


def layout() -> html.Div:
    group_opts = [{"label": g, "value": g} for g, _ in ALL_GROUPS]
    first_plots = [{"label": l, "value": f} for l, f, _ in ALL_GROUPS[0][1]]

    controls = html.Div([
        html.Div([
            html.Div("Group", className="ctrl-label"),
            dbc.Select(id=IDs.QC_GROUP_SELECT, options=group_opts,
                       value=group_opts[0]["value"], size="sm"),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Plot", className="ctrl-label"),
            dbc.Select(id=IDs.QC_PLOT_SELECT, options=first_plots,
                       value=first_plots[0]["value"], size="sm"),
        ], className="ctrl-card"),
        html.Div([
            dbc.Button([html.I(className="bi bi-play-fill me-1"), "Run Plot"],
                       id=IDs.BTN_QC_RUN, color="info", size="sm",
                       className="w-100 mb-2"),
            dbc.Spinner(html.Div(id=IDs.QC_SPINNER), size="sm", color="info"),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Categories", className="ctrl-label"),
            html.Ul([
                html.Li(f"{g}  ({len(plots)})",
                        style={"fontSize": "11px", "color": "#a6adc8",
                               "lineHeight": "1.8"})
                for g, plots in ALL_GROUPS
            ], style={"paddingLeft": "14px", "margin": 0}),
        ], className="ctrl-card"),
    ], className="analysis-controls")

    output = html.Div([
        html.Div([
            html.Img(id=IDs.IMG_QC, src=empty_src(dark=True),
                     style={"width": "100%", "maxHeight": "100%",
                            "objectFit": "contain"}),
        ], className="fig-area", style={"minHeight": "500px"}),
    ], className="analysis-output")

    from pycsamt.app.web.layout import _command_bar
    return html.Div([
        _command_bar("qc", "Quality Control",
                     "25 diagnostic plots · 6 categories"),
        html.Div([controls, output], className="analysis-layout", style={"flex": "1"}),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # handled by callbacks/qc.py
