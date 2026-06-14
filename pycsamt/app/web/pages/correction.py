# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Data correction / conditioning page (CorrectionController)."""
from __future__ import annotations

from dash import html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.correction_controller import CATALOGUE, CATEGORIES
from pycsamt.app.web.layout import IDs, _icon
from pycsamt.app.web.utils import empty_src

PAGE_ID = "correction"

_CAT_OPTS    = [{"label": c, "value": c} for c in CATEGORIES]
_FIRST_METHS = [{"label": k, "value": k} for k in CATALOGUE[CATEGORIES[0]]]


def layout() -> html.Div:
    controls = html.Div([
        # Category + method
        html.Div([
            html.Div("Category", className="ctrl-label"),
            dbc.Select(id=IDs.CORR_CAT, options=_CAT_OPTS,
                       value=_CAT_OPTS[0]["value"], size="sm"),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Correction", className="ctrl-label"),
            dbc.Select(id=IDs.CORR_METHOD, options=_FIRST_METHS,
                       value=_FIRST_METHS[0]["value"], size="sm"),
        ], className="ctrl-card"),
        # Action buttons
        html.Div([
            dbc.Button([html.I(className="bi bi-eye me-1"), "Preview"],
                       id=IDs.BTN_CORR_PREVIEW, color="secondary",
                       size="sm", className="w-100 mb-2"),
            dbc.Button([html.I(className="bi bi-check-lg me-1"), "Apply"],
                       id=IDs.BTN_CORR_APPLY, color="success",
                       size="sm", className="w-100 mb-2"),
            dbc.Row([
                dbc.Col(dbc.Button([html.I(className="bi bi-arrow-counterclockwise me-1"),
                                    "Undo"],
                                   id=IDs.BTN_CORR_UNDO, color="warning",
                                   size="sm", className="w-100")),
                dbc.Col(dbc.Button([html.I(className="bi bi-x-circle me-1"), "Reset"],
                                   id=IDs.BTN_CORR_RESET, color="danger",
                                   size="sm", className="w-100")),
            ], className="g-1"),
            html.Div(
                dbc.Spinner(html.Div(id=IDs.CORR_SPINNER), size="sm", color="success"),
                className="mt-2",
            ),
            html.Div(id=IDs.CORR_FEEDBACK,
                     className="mt-2 text-muted small"),
        ], className="ctrl-card"),
        # Correction stack
        html.Div([
            html.Div("Applied corrections", className="ctrl-label"),
            html.Div(id=IDs.CORR_STACK,
                     children=[html.Span("None yet",
                                         className="text-muted small")]),
        ], className="ctrl-card", style={"flex": "1", "overflowY": "auto"}),
    ], className="analysis-controls")

    output = html.Div([
        dbc.Tabs([
            dbc.Tab(
                html.Div([
                    html.Div([
                        html.Div("Raw", className="fig-label"),
                        html.Div(
                            html.Img(id=IDs.IMG_CORR_LEFT,
                                     src=empty_src(dark=True),
                                     style={"width": "100%", "height": "100%",
                                            "objectFit": "contain"}),
                            className="fig-area",
                        ),
                    ], style={"flex": "1"}),
                    html.Div([
                        html.Div("Corrected", className="fig-label"),
                        html.Div(
                            html.Img(id=IDs.IMG_CORR_RIGHT,
                                     src=empty_src(dark=True),
                                     style={"width": "100%", "height": "100%",
                                            "objectFit": "contain"}),
                            className="fig-area",
                        ),
                    ], style={"flex": "1"}),
                ], className="dual-fig", style={"minHeight": "480px"}),
                label="Before / After", tab_id="corr-tab-ba",
            ),
            dbc.Tab(
                html.Div(
                    html.Img(id=IDs.IMG_CORR_OVERLAY,
                             src=empty_src(dark=True),
                             style={"width": "100%", "maxHeight": "500px",
                                    "objectFit": "contain"}),
                    className="fig-area", style={"minHeight": "480px"},
                ),
                label="Overlay", tab_id="corr-tab-overlay",
            ),
            dbc.Tab(
                html.Div(
                    html.Img(id=IDs.IMG_CORR_DIFF,
                             src=empty_src(dark=True),
                             style={"width": "100%", "maxHeight": "500px",
                                    "objectFit": "contain"}),
                    className="fig-area", style={"minHeight": "480px"},
                ),
                label="Difference", tab_id="corr-tab-diff",
            ),
        ], id=IDs.CORR_VIEW, active_tab="corr-tab-ba"),
    ], className="analysis-output")

    from pycsamt.app.web.layout import _command_bar
    return html.Div([
        _command_bar("sites-correction", "Data Correction",
                     "13 corrections · Static Shift · Noise · Tensor Rotation · Coordinates"),
        html.Div([controls, output], className="analysis-layout", style={"flex": "1"}),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # handled by callbacks/correction.py
