# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Geological interpretation page (InterpController — 10 workflow categories)."""
from __future__ import annotations

from dash import html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.interp_controller import (
    WORKFLOW_CATALOGUE, CATEGORIES,
)
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import empty_src

PAGE_ID = "interpretation"

_CAT_OPTS = [{"label": c, "value": c} for c in CATEGORIES]

def _first_plots(cat: str):
    return [{"label": label, "value": func}
            for label, func, _ in WORKFLOW_CATALOGUE[cat]]

_FIRST_PLOT_OPTS = _first_plots(CATEGORIES[0])


def layout() -> html.Div:
    controls = html.Div([
        html.Div([
            html.Div("Workflow Category", className="ctrl-label"),
            dbc.Select(
                id=IDs.INTERP_CAT,
                options=_CAT_OPTS,
                value=_CAT_OPTS[0]["value"],
                size="sm",
            ),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Plot", className="ctrl-label"),
            dbc.Select(
                id=IDs.INTERP_PLOT,
                options=_FIRST_PLOT_OPTS,
                value=_FIRST_PLOT_OPTS[0]["value"],
                size="sm",
            ),
            html.Div(id=IDs.INTERP_DESC,
                     children=WORKFLOW_CATALOGUE[CATEGORIES[0]][0][2],
                     className="mt-1 text-muted",
                     style={"fontSize": "10px", "lineHeight": "1.4"}),
        ], className="ctrl-card"),
        html.Div([
            dbc.Button(
                [html.I(className="bi bi-play-fill me-1"), "Generate"],
                id=IDs.BTN_INTERP_RUN, color="primary", size="sm",
                className="w-100 mb-2",
            ),
            dbc.Spinner(html.Div(id=IDs.INTERP_SPINNER),
                        size="sm", color="primary"),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Category summary", className="ctrl-label"),
            html.Ul([
                html.Li(
                    f"{cat}: {len(plots)} plot{'s' if len(plots)!=1 else ''}",
                    className="text-muted small",
                    style={"fontSize": "10px"},
                )
                for cat, plots in WORKFLOW_CATALOGUE.items()
            ], style={"paddingLeft": "14px", "margin": "0"}),
        ], className="ctrl-card", style={"flex": "1", "overflowY": "auto"}),
    ], className="analysis-controls")

    output = html.Div([
        html.Div([
            html.Img(
                id=IDs.IMG_INTERP,
                src=empty_src(dark=True),
                style={"width": "100%", "maxHeight": "100%",
                       "objectFit": "contain"},
            ),
        ], className="fig-area", style={"minHeight": "500px"}),
    ], className="analysis-output")

    total_plots = sum(len(v) for v in WORKFLOW_CATALOGUE.values())

    from pycsamt.app.web.layout import _command_bar
    return html.Div([
        _command_bar("interpret", "Geological Interpretation",
                     f"{total_plots} workflows · {len(CATEGORIES)} categories"
                     " · Hydrology · EM Diagnostics · Advanced · Fusion"),
        html.Div([controls, output], className="analysis-layout",
                 style={"flex": "1"}),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # handled by callbacks/interp.py
