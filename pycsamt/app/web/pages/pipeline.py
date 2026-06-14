# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Processing pipeline page (PipelineController — 8 sequential steps)."""
from __future__ import annotations

from dash import dcc, html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.pipeline_controller import _build_steps
from pycsamt.app.web.layout import IDs, _command_bar, _icon
from pycsamt.app.web.utils import empty_src

PAGE_ID = "pipeline"

_STEPS = _build_steps()
_STEP_OPTS = [{"label": f"{s.id + 1}. {s.name}", "value": s.id} for s in _STEPS]


def _horizontal_stepper() -> html.Div:
    """Visual progress stepper rendered above the output area."""
    nodes = []
    for s in _STEPS:
        nodes.append(
            html.Div(
                [
                    html.Div(
                        str(s.id + 1),
                        className="pipe-step-circle",
                        id=f"pipe-circle-{s.id}",
                    ),
                    html.Div(
                        s.name,
                        className="pipe-step-label",
                        id=f"pipe-label-{s.id}",
                    ),
                ],
                className="pipe-step-node" + (" step-active" if s.id == 0 else " step-pending"),
                id=f"pipe-node-{s.id}",
                n_clicks=0,
                **{"data-step": str(s.id)},
            )
        )
    return html.Div(nodes, className="pipe-stepper")


def layout() -> html.Div:
    # Sidebar step list (compact, shows status dots)
    step_list = html.Div([
        html.Div(
            [
                html.Div(className="pipe-step-dot dot-pending",
                         id=f"pipe-dot-{s.id}"),
                html.Span(s.name, style={"fontSize": "11.5px"}),
            ],
            className="pipe-step-item" + (" active" if s.id == 0 else ""),
            id=f"pipe-step-row-{s.id}",
            **{"data-step": s.id},
        )
        for s in _STEPS
    ])

    controls = html.Div([
        html.Div([
            html.Div("Pipeline Step", className="ctrl-label"),
            dbc.Select(id=IDs.PIPE_STEP, options=_STEP_OPTS,
                       value=0, size="sm"),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Step Description", className="ctrl-label"),
            html.Div(id=IDs.PIPE_STEP_INFO,
                     children=_STEPS[0].description,
                     style={"fontSize": "11px", "color": "#a6adc8",
                            "lineHeight": "1.5"}),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Method", className="ctrl-label"),
            dbc.Select(id=IDs.PIPE_METHOD,
                       options=[{"label": m.label, "value": m.name}
                                for m in _STEPS[0].methods],
                       value=_STEPS[0].default_method, size="sm"),
        ], className="ctrl-card"),
        html.Div([
            dbc.Button([html.I(className="bi bi-play-fill me-1"), "Run Step"],
                       id=IDs.BTN_PIPE_RUN, color="success", size="sm",
                       className="w-100 mb-2"),
            dbc.Button([html.I(className="bi bi-fast-forward-fill me-1"), "Run All"],
                       id=IDs.BTN_PIPE_ALL, color="primary", size="sm",
                       className="w-100 mb-2"),
            dbc.Row([
                dbc.Col(dbc.Button([html.I(className="bi bi-skip-forward me-1"), "Skip"],
                                   id=IDs.BTN_PIPE_SKIP,
                                   color="secondary", size="sm", className="w-100")),
                dbc.Col(dbc.Button([html.I(className="bi bi-arrow-counterclockwise me-1"),
                                    "Reset"],
                                   id=IDs.BTN_PIPE_RESET,
                                   color="danger", size="sm", className="w-100")),
            ], className="g-1"),
            html.Div(
                dbc.Spinner(html.Div(id=IDs.PIPE_SPINNER), size="sm", color="success"),
                className="mt-2 text-center",
            ),
        ], className="ctrl-card"),
        html.Div(
            [
                html.Div([html.I(className="bi bi-list-check me-1"), "Steps"],
                         className="ctrl-label"),
                step_list,
            ],
            className="ctrl-card",
            style={"flex": "1", "overflowY": "auto"},
        ),
        dcc.Store(id=IDs.PIPE_STORE, storage_type="memory"),
    ], className="analysis-controls")

    output = html.Div([
        # Horizontal stepper at the top of the output area
        _horizontal_stepper(),
        dbc.Tabs([
            dbc.Tab(
                html.Pre("", id=IDs.PIPE_LOG, className="log-area",
                         style={"height": "460px", "margin": "8px"}),
                label="Execution Log", tab_id="pipe-tab-log",
            ),
            dbc.Tab(
                html.Div(
                    html.Img(id=IDs.IMG_PIPE, src=empty_src(dark=True),
                             style={"width": "100%", "maxHeight": "100%",
                                    "objectFit": "contain"}),
                    className="fig-area",
                    style={"minHeight": "440px", "margin": "8px"},
                ),
                label="Step Preview", tab_id="pipe-tab-preview",
            ),
            dbc.Tab(
                html.Div(id=IDs.PIPE_STATUS,
                         style={"padding": "10px", "fontSize": "12px"}),
                label="Step Status", tab_id="pipe-tab-status",
            ),
        ], active_tab="pipe-tab-log"),
    ], className="analysis-output", style={"padding": "0"})

    return html.Div([
        _command_bar(
            "pipeline", "Processing Pipeline",
            f"{len(_STEPS)}-step workflow · QC → Freq Edit → Static Shift → Noise → Strike → Export",
        ),
        html.Div([controls, output], className="analysis-layout", style={"flex": "1"}),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # handled by callbacks/pipeline.py
