# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the Processing Pipeline page."""
from __future__ import annotations

import traceback

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.pipeline_controller import (
    PipelineController, _build_steps, StepStatus,
)
from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import apply_web_dark_theme, empty_src, fig_to_src

_CTRL  = PipelineController()
_STEPS = _build_steps()

_STATUS_CLASS: dict = {
    StepStatus.PENDING: "step-pending",
    StepStatus.RUNNING: "step-running",
    StepStatus.DONE:    "step-done",
    StepStatus.ERROR:   "step-failed",
    StepStatus.SKIPPED: "step-skipped",
}


def _stepper_classes() -> list[str]:
    """Return pipe-step-node className for each step based on controller state."""
    classes = []
    active_assigned = False
    for step in _CTRL.steps:
        sc = _STATUS_CLASS.get(step.status, "step-pending")
        if step.status == StepStatus.PENDING and not active_assigned:
            prev_done = all(
                s.status in (StepStatus.DONE, StepStatus.SKIPPED)
                for s in _CTRL.steps[: step.id]
            )
            if prev_done or step.id == 0:
                sc = "step-active"
                active_assigned = True
        classes.append(f"pipe-step-node {sc}")
    return classes


def register_pipeline(app) -> None:

    @app.callback(
        Output(IDs.PIPE_STEP_INFO, "children"),
        Output(IDs.PIPE_METHOD,    "options"),
        Output(IDs.PIPE_METHOD,    "value"),
        Input(IDs.PIPE_STEP,       "value"),
    )
    def sync_step_info(step_id):
        if step_id is None:
            raise PreventUpdate
        step_id = int(step_id)
        step = _STEPS[step_id]
        desc = step.description
        opts = [{"label": m.label, "value": m.name} for m in step.methods]
        return desc, opts, step.default_method

    @app.callback(
        Output(IDs.PIPE_LOG,      "children"),
        Output(IDs.IMG_PIPE,      "src"),
        Output(IDs.PIPE_SPINNER,  "children"),
        Output(IDs.TOAST_ERROR,   "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,    "children", allow_duplicate=True),
        Input(IDs.BTN_PIPE_RUN,   "n_clicks"),
        State(IDs.PIPE_STEP,      "value"),
        State(IDs.PIPE_METHOD,    "value"),
        State(IDs.SESSION_ID,     "data"),
        prevent_initial_call=True,
    )
    def run_step(n_clicks, step_id, method, session_id):
        if not n_clicks:
            raise PreventUpdate
        log_lines = []
        try:
            sites = cache_get(session_id)
            _CTRL.set_input_sites(sites)
            step_id = int(step_id)
            step = _CTRL.steps[step_id]
            step.set_active_method(method)

            _CTRL.execute_step(step_id, log_cb=lambda m: log_lines.append(m))
            log_text = "\n".join(log_lines)

            # Try to render a plot from the output sites
            result_sites = _CTRL._sites_chain[step_id]
            src = empty_src(dark=True)
            if result_sites is not None:
                try:
                    apply_web_dark_theme()
                    fig = plt.figure(figsize=(10, 4))
                    ax = fig.add_subplot(111)
                    result_sites.plot_pseudosection(ax=ax, component="res")
                    src = fig_to_src(fig)
                except Exception:
                    pass

            return log_text, src, "", False, ""
        except Exception as exc:
            log_lines.append(f"ERROR: {exc}")
            return "\n".join(log_lines), empty_src(dark=True), "", True, str(exc)

    @app.callback(
        Output(IDs.PIPE_LOG,      "children",  allow_duplicate=True),
        Output(IDs.PIPE_SPINNER,  "children",  allow_duplicate=True),
        Output(IDs.TOAST_ERROR,   "is_open",   allow_duplicate=True),
        Output(IDs.TOAST_BODY,    "children",  allow_duplicate=True),
        Input(IDs.BTN_PIPE_ALL,   "n_clicks"),
        State(IDs.SESSION_ID,     "data"),
        prevent_initial_call=True,
    )
    def run_all(n_clicks, session_id):
        if not n_clicks:
            raise PreventUpdate
        log_lines = []
        try:
            sites = cache_get(session_id)
            _CTRL.set_input_sites(sites)
            for i in range(len(_CTRL.steps)):
                try:
                    _CTRL.execute_step(i, log_cb=lambda m: log_lines.append(m))
                    log_lines.append(f"── Step {i+1} done ──")
                except Exception as exc:
                    log_lines.append(f"Step {i+1} ERROR: {exc}")
                    break
            return "\n".join(log_lines), "", False, ""
        except Exception as exc:
            return "\n".join(log_lines), "", True, str(exc)

    @app.callback(
        Output(IDs.PIPE_LOG,     "children",  allow_duplicate=True),
        Input(IDs.BTN_PIPE_SKIP, "n_clicks"),
        State(IDs.PIPE_STEP,     "value"),
        prevent_initial_call=True,
    )
    def skip_step(n_clicks, step_id):
        if not n_clicks:
            raise PreventUpdate
        step_id = int(step_id)
        _CTRL.steps[step_id].status = StepStatus.SKIPPED
        return f"Step {step_id + 1} skipped."

    @app.callback(
        Output(IDs.PIPE_LOG,      "children",  allow_duplicate=True),
        Output(IDs.PIPE_STATUS,   "children"),
        *[Output(f"pipe-node-{s.id}", "className", allow_duplicate=True) for s in _STEPS],
        Input(IDs.BTN_PIPE_RESET, "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_pipeline(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        _CTRL.reset()
        # First step active, rest pending after reset
        node_classes = ["pipe-step-node step-active"] + [
            "pipe-step-node step-pending" for _ in _STEPS[1:]
        ]
        return ("Pipeline reset.", "", *node_classes)

    _STEP_COLOUR = {
        StepStatus.PENDING: "#a6adc8",
        StepStatus.RUNNING: "#f9e2af",
        StepStatus.DONE:    "#a6e3a1",
        StepStatus.ERROR:   "#f38ba8",
        StepStatus.SKIPPED: "#cba6f7",
    }

    @app.callback(
        Output(IDs.PIPE_STATUS,  "children",  allow_duplicate=True),
        *[Output(f"pipe-node-{s.id}", "className", allow_duplicate=True) for s in _STEPS],
        Input(IDs.PIPE_LOG,      "children"),
        prevent_initial_call=True,
    )
    def update_status(_log):
        rows = []
        for step in _CTRL.steps:
            colour = _STEP_COLOUR.get(step.status, "#a6adc8")
            rows.append(
                html.Div(
                    [
                        html.Span(f"{step.id + 1}. {step.name}",
                                  style={"color": colour, "fontWeight": "500"}),
                        html.Span(f"  [{step.status.value}]",
                                  style={"color": colour, "fontSize": "10px"}),
                        html.Span(
                            f"  {step.elapsed_s:.2f}s" if step.elapsed_s else "",
                            style={"color": "#585b70", "fontSize": "10px"},
                        ),
                    ],
                    style={"marginBottom": "4px"},
                )
            )
        return (rows, *_stepper_classes())
