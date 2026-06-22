# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Processing pipeline page v2 — redesigned to match Correction / QC layout.

Layout anatomy
--------------
  ┌─ command bar ──────────────────────────────────────────────────────────┐
  ├─ pipe-body (flex-row, flex:1) ─────────────────────────────────────────┤
  │ ┌─ pipe-sidebar ──────┐  ┌─ pipe-view-area ───────────────────────── ┐ │
  │ │  pipe-run-bar (top, │  │  pipe-stepper (fixed — progress nodes)    │ │
  │ │   FIXED):           │  │  pipe-tab-bar (Log | Preview | Status)    │ │
  │ │  [Run Step][Run All]│  │  pipe-panels-wrap                         │ │
  │ │  [Skip]    [Reset]  │  │    pipe-panel-log  <pre PIPE_LOG>         │ │
  │ │  ⟳ spinner  text    │  │    pipe-panel-preview step-bar + IMG_PIPE │ │
  │ ├─────────────────────┤  │    pipe-panel-status PIPE_STATUS rows     │ │
  │ │  pipe-ctrl-scroll   │  └───────────────────────────────────────────┘ │
  │ │   (SCROLLABLE)      │                                                 │
  │ │   Step select       │                                                 │
  │ │   Method select     │                                                 │
  │ │   Step description  │                                                 │
  │ │   Output folder     │                                                 │
  │ │   Step status list  │                                                 │
  │ └─────────────────────┘                                                 │
  └─────────────────────────────────────────────────────────────────────────┘
"""
from __future__ import annotations

import os

from dash import dcc, html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.pipeline_controller import _build_steps
from pycsamt.app.web.layout import IDs, _command_bar
from pycsamt.app.web.utils import empty_src

PAGE_ID = "pipeline"

_STEPS = _build_steps()
_STEP_OPTS = [{"label": f"{s.id + 1}. {s.name}", "value": str(s.id)} for s in _STEPS]


# ── Folder-browser modal ──────────────────────────────────────────────────────

def _browse_modal() -> dbc.Modal:
    return dbc.Modal(
        [
            dbc.ModalHeader(
                dbc.ModalTitle(
                    [html.I(className="bi bi-folder2-open me-2"), "Select Export Folder"]
                ),
                close_button=True,
            ),
            dbc.ModalBody([
                html.Div([
                    html.Span(
                        [html.I(className="bi bi-hdd me-1"), "Current folder:"],
                        style={"fontSize": "11px", "color": "#a6adc8",
                               "marginRight": "6px"},
                    ),
                    html.Code(
                        id=IDs.PIPE_BROWSE_CWD,
                        children=os.path.expanduser("~"),
                        style={"fontSize": "11px", "wordBreak": "break-all"},
                    ),
                ], className="mb-2"),
                dbc.Button(
                    [html.I(className="bi bi-arrow-up me-1"), "Parent folder"],
                    id=IDs.BTN_PIPE_BROWSE_UP,
                    color="secondary", outline=True, size="sm", className="mb-2",
                ),
                html.Div(
                    id=IDs.PIPE_BROWSE_LIST,
                    style={
                        "maxHeight": "320px", "overflowY": "auto",
                        "border": "1px solid var(--surface1)",
                        "borderRadius": "6px", "padding": "4px",
                    },
                ),
                html.Hr(style={"margin": "12px 0"}),
                html.Div([
                    html.Div("Create new folder inside current location:",
                             style={"fontSize": "11px", "color": "#a6adc8",
                                    "marginBottom": "4px"}),
                    dbc.InputGroup([
                        dbc.Input(
                            id=IDs.PIPE_BROWSE_NEW, placeholder="new_folder_name",
                            type="text", size="sm",
                        ),
                        dbc.Button(
                            [html.I(className="bi bi-folder-plus me-1"), "Create"],
                            id=IDs.BTN_PIPE_BROWSE_MK, color="secondary", size="sm",
                        ),
                    ], size="sm"),
                ]),
                dcc.Store(id=IDs.PIPE_BROWSE_PATH, storage_type="memory",
                          data=os.path.expanduser("~")),
            ]),
            dbc.ModalFooter([
                dbc.Button("Cancel", id="pipe-browse-cancel",
                           color="secondary", size="sm", className="me-auto"),
                dbc.Button(
                    [html.I(className="bi bi-check2-circle me-1"), "Select this folder"],
                    id=IDs.BTN_PIPE_BROWSE_SEL, color="success", size="sm",
                ),
            ]),
        ],
        id=IDs.PIPE_BROWSE_MODAL, is_open=False, size="lg", backdrop="static",
    )


# ── Horizontal stepper (fixed in view area) ───────────────────────────────────

def _stepper() -> html.Div:
    nodes = []
    for s in _STEPS:
        nodes.append(
            html.Div(
                [
                    html.Div(str(s.id + 1), className="pipe-step-circle",
                             id=f"pipe-circle-{s.id}"),
                    html.Div(s.name, className="pipe-step-label",
                             id=f"pipe-label-{s.id}"),
                ],
                className="pipe-step-node" + (" step-active" if s.id == 0 else " step-pending"),
                id=f"pipe-node-{s.id}",
                n_clicks=0,
                **{"data-step": str(s.id)},
            )
        )
    return html.Div(nodes, className="pipe-stepper")


# ── Sidebar step status list ──────────────────────────────────────────────────

def _step_list() -> html.Div:
    return html.Div(
        [
            html.Div(
                [
                    html.Div(className="pipe-step-dot dot-pending",
                             id=f"pipe-dot-{s.id}"),
                    html.Span(s.name, style={"fontSize": "11px"}),
                ],
                className="pipe-step-item" + (" active" if s.id == 0 else ""),
                id=f"pipe-step-row-{s.id}",
                **{"data-step": s.id},
            )
            for s in _STEPS
        ]
    )


# ── Page layout ───────────────────────────────────────────────────────────────

def layout() -> html.Div:

    # ── Sidebar ───────────────────────────────────────────────────────────────
    sidebar = html.Div(
        [
            # Run bar — FIXED at top (mirrors fwd-run-bar / corr-run-bar)
            html.Div(
                [
                    # Row 1: Run Step | Run All
                    html.Div(
                        [
                            dbc.Button(
                                [html.I(className="bi bi-play-fill me-1"), "Run Step"],
                                id=IDs.BTN_PIPE_RUN, color="success", size="sm",
                                className="w-100",
                            ),
                            dbc.Button(
                                [html.I(className="bi bi-fast-forward-fill me-1"), "Run All"],
                                id=IDs.BTN_PIPE_ALL, color="primary", size="sm",
                                className="w-100",
                            ),
                        ],
                        className="d-flex gap-1 mb-1",
                    ),
                    # Row 2: Skip | Reset
                    html.Div(
                        [
                            dbc.Button(
                                [html.I(className="bi bi-skip-forward me-1"), "Skip"],
                                id=IDs.BTN_PIPE_SKIP, color="secondary", size="sm",
                                className="w-100",
                            ),
                            dbc.Button(
                                [html.I(className="bi bi-arrow-counterclockwise me-1"), "Reset"],
                                id=IDs.BTN_PIPE_RESET, color="danger", size="sm",
                                className="w-100",
                            ),
                        ],
                        className="d-flex gap-1 mb-1",
                    ),
                    # Spinner + feedback
                    dbc.Spinner(
                        html.Div(id=IDs.PIPE_SPINNER),
                        size="sm", color="success",
                    ),
                    html.Div(id=IDs.PIPE_FEEDBACK, className="pipe-feedback-txt"),
                ],
                className="pipe-run-bar",
            ),

            # Scrollable controls
            html.Div(
                [
                    # Step selector
                    html.Div(
                        [
                            html.Div("Pipeline Step", className="ctrl-label"),
                            dbc.Select(
                                id=IDs.PIPE_STEP,
                                options=_STEP_OPTS,
                                value="0",
                                size="sm",
                            ),
                        ],
                        className="ctrl-card",
                    ),

                    # Method selector
                    html.Div(
                        [
                            html.Div("Method", className="ctrl-label"),
                            dbc.Select(
                                id=IDs.PIPE_METHOD,
                                options=[{"label": m.label, "value": m.name}
                                         for m in _STEPS[0].methods],
                                value=_STEPS[0].default_method,
                                size="sm",
                            ),
                        ],
                        className="ctrl-card",
                    ),

                    # Step description
                    html.Div(
                        [
                            html.Div("Description", className="ctrl-label"),
                            html.Div(
                                id=IDs.PIPE_STEP_INFO,
                                children=_STEPS[0].description,
                                className="pipe-desc-txt",
                            ),
                        ],
                        className="ctrl-card",
                    ),

                    # Output folder (step 7 only — shown/hidden via callback)
                    html.Div(
                        [
                            html.Div("Output Folder", className="ctrl-label"),
                            dbc.InputGroup(
                                [
                                    dbc.Input(
                                        id=IDs.PIPE_EXPORT_FOLDER,
                                        placeholder="/path/to/export/",
                                        type="text",
                                        size="sm",
                                    ),
                                    dbc.Button(
                                        html.I(className="bi bi-folder2-open"),
                                        id=IDs.BTN_PIPE_BROWSE,
                                        color="secondary",
                                        size="sm",
                                    ),
                                ],
                                size="sm",
                            ),
                            html.Div(
                                "Folder is created if it does not exist.",
                                className="pipe-hint-txt",
                            ),
                        ],
                        className="ctrl-card",
                        id=IDs.PIPE_EXPORT_CARD,
                        style={"display": "none"},
                    ),

                    # Step progress list
                    html.Div(
                        [
                            html.Div(
                                [html.I(className="bi bi-list-check me-1"), "Steps"],
                                className="ctrl-label",
                            ),
                            _step_list(),
                        ],
                        className="ctrl-card",
                    ),

                    # Stores live here (non-visible)
                    dcc.Store(id=IDs.PIPE_STORE, storage_type="memory", data={}),
                    dcc.Store(id=IDs.PIPE_ACTIVE_TAB, data="log"),
                ],
                className="pipe-ctrl-scroll",
            ),
        ],
        className="pipe-sidebar",
    )

    # ── View area ─────────────────────────────────────────────────────────────
    view_area = html.Div(
        [
            # Horizontal stepper — fixed at top of view area
            _stepper(),

            # Tab bar
            html.Div(
                [
                    html.Button(
                        [html.I(className="bi bi-terminal me-1"), "Log"],
                        id="pipe-tab-btn-log",
                        className="pipe-tab-btn active",
                        n_clicks=0,
                    ),
                    html.Button(
                        [html.I(className="bi bi-image me-1"), "Preview"],
                        id="pipe-tab-btn-preview",
                        className="pipe-tab-btn",
                        n_clicks=0,
                    ),
                    html.Button(
                        [html.I(className="bi bi-clipboard-check me-1"), "Status"],
                        id="pipe-tab-btn-status",
                        className="pipe-tab-btn",
                        n_clicks=0,
                    ),
                ],
                className="pipe-tab-bar",
            ),

            # Panels
            html.Div(
                [
                    # Log panel
                    html.Div(
                        html.Pre(
                            "",
                            id=IDs.PIPE_LOG,
                            className="log-area pipe-log-pre",
                        ),
                        id="pipe-panel-log",
                        className="pipe-panel pipe-panel-active pipe-log-wrap",
                    ),

                    # Preview panel
                    html.Div(
                        [
                            # Step picker bar
                            html.Div(
                                [
                                    html.Span(
                                        "View step:",
                                        className="pipe-preview-label",
                                    ),
                                    dbc.Select(
                                        id=IDs.PIPE_VIEW_STEP,
                                        options=_STEP_OPTS,
                                        value="",
                                        size="sm",
                                        style={"width": "200px"},
                                    ),
                                ],
                                className="pipe-preview-bar",
                            ),
                            html.Div(
                                html.Img(
                                    id=IDs.IMG_PIPE,
                                    src=empty_src(dark=True),
                                    style={
                                        "width": "100%",
                                        "height": "100%",
                                        "objectFit": "contain",
                                    },
                                ),
                                className="pipe-preview-img-wrap",
                            ),
                        ],
                        id="pipe-panel-preview",
                        className="pipe-panel pipe-preview-panel",
                        style={"display": "none"},
                    ),

                    # Status panel
                    html.Div(
                        html.Div(
                            id=IDs.PIPE_STATUS,
                            className="pipe-status-scroll",
                        ),
                        id="pipe-panel-status",
                        className="pipe-panel",
                        style={"display": "none"},
                    ),
                ],
                className="pipe-panels-wrap",
            ),
        ],
        className="pipe-view-area",
    )

    return html.Div(
        [
            _command_bar(
                "pipeline",
                "Processing Pipeline",
                f"{len(_STEPS)}-step workflow · QC → Freq Edit → Static Shift"
                " → Noise → Strike → Export",
            ),
            html.Div(
                [sidebar, view_area],
                className="pipe-body",
            ),
            _browse_modal(),
        ],
        style={"display": "flex", "flexDirection": "column", "height": "100%"},
    )


def register_callbacks(app) -> None:
    pass  # handled by callbacks/pipeline.py
