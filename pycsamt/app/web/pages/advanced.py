# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Advanced emtools plots page."""
from __future__ import annotations

from dash import html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.advanced_controller import ADVANCED_GROUPS
from pycsamt.app.web.layout import IDs, _icon
from pycsamt.app.web.utils import empty_src

PAGE_ID = "advanced"

_GROUPS = [(g, plots) for g, plots in ADVANCED_GROUPS if plots]  # skip topo/conv stubs


def layout() -> html.Div:
    group_opts  = [{"label": g, "value": g} for g, _ in _GROUPS]
    first_plots = [{"label": l, "value": f} for l, f, _ in _GROUPS[0][1]]

    controls = html.Div([
        html.Div([
            html.Div("Analysis Group", className="ctrl-label"),
            dbc.Select(id=IDs.ADV_GROUP, options=group_opts,
                       value=group_opts[0]["value"], size="sm"),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Plot", className="ctrl-label"),
            dbc.Select(id=IDs.ADV_PLOT, options=first_plots,
                       value=first_plots[0]["value"], size="sm"),
        ], className="ctrl-card"),
        html.Div([
            dbc.Button([html.I(className="bi bi-play-fill me-1"), "Generate"],
                       id=IDs.BTN_ADV_RUN, color="primary", size="sm",
                       className="w-100 mb-2"),
            dbc.Spinner(html.Div(id=IDs.ADV_SPINNER), size="sm", color="primary"),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Available Groups", className="ctrl-label"),
            html.Ul([
                html.Li(f"{g}  ({len(p)})",
                        style={"fontSize": "11px", "color": "#a6adc8",
                               "lineHeight": "1.8"})
                for g, p in _GROUPS
            ], style={"paddingLeft": "14px", "margin": 0}),
        ], className="ctrl-card"),
    ], className="analysis-controls")

    output = html.Div([
        html.Div([
            html.Img(id=IDs.IMG_ADV, src=empty_src(dark=True),
                     style={"width": "100%", "maxHeight": "100%",
                            "objectFit": "contain"}),
        ], className="fig-area", style={"minHeight": "500px"}),
    ], className="analysis-output")

    from pycsamt.app.web.layout import _command_bar
    return html.Div([
        _command_bar("advanced-tools", "Advanced Analysis",
                     "Strike · Phase Tensor · Induction · Impedance · Depth"),
        html.Div([controls, output], className="analysis-layout", style={"flex": "1"}),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # handled by callbacks/advanced.py
