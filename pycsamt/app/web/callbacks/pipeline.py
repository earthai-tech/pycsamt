# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the Processing Pipeline page."""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import dash_bootstrap_components as dbc
from dash import ALL, Input, Output, State, ctx, html, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.pipeline_controller import (
    PipelineController, _build_steps, StepStatus,
)
from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    apply_web_light_theme,
    empty_src,
    fig_to_src,
    filter_sites_by_lines,
    no_active_lines_src,
)

_CTRL  = PipelineController()
_STEPS = _build_steps()

_STATUS_CLASS: dict = {
    StepStatus.PENDING: "step-pending",
    StepStatus.RUNNING: "step-running",
    StepStatus.DONE:    "step-done",
    StepStatus.ERROR:   "step-failed",
    StepStatus.SKIPPED: "step-skipped",
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def _stepper_classes() -> list[str]:
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


def _result_to_figure(obj):
    """Normalise an Axes or Figure return value to a Figure."""
    try:
        import matplotlib.figure as _mf
        if isinstance(obj, _mf.Figure):
            return obj
        get_fig = getattr(obj, "get_figure", None)
        if callable(get_fig):
            f = get_fig()
            if isinstance(f, _mf.Figure):
                return f
    except Exception:
        pass
    return None


def _generate_preview(step, sites, dark: bool) -> str:
    """Generate a diagnostic preview for a completed step.

    Tries the step's diag_fn first; falls back to a ρ pseudo-section.
    Always returns a valid base64 src (empty placeholder on total failure).
    """
    import pycsamt.emtools as et

    if dark:
        apply_web_dark_theme()
    else:
        apply_web_light_theme()

    if step.diag_fn and sites is not None:
        fn = getattr(et, step.diag_fn, None)
        if fn is not None:
            try:
                result = fn(sites)
                fig = _result_to_figure(result)
                if fig is not None:
                    src = fig_to_src(fig)
                    plt.close(fig)
                    return src
            except Exception:
                plt.close("all")

    # Fallback: simple pseudosection
    if sites is not None:
        try:
            fig = plt.figure(figsize=(10, 4))
            ax  = fig.add_subplot(111)
            sites.plot_pseudosection(ax=ax, component="res")
            src = fig_to_src(fig)
            plt.close(fig)
            return src
        except Exception:
            plt.close("all")

    return empty_src(dark=dark)


def _filter_by_active_lines(sites, store_data, active_lines_store):
    """Return (sites, all_muted).  Unchanged when nothing is muted."""
    _als = active_lines_store or {}
    _all = _als.get("all", [])
    _active = set(_als.get("active", _all))
    if not _all or _active == set(_all):
        return sites, False
    muted = [ln for ln in _all if ln not in _active]
    if not muted:
        return sites, False
    records = (store_data or {}).get("station_records", [])
    filtered = filter_sites_by_lines(sites, records, list(_active))
    if filtered is None:
        return None, True
    return filtered, False


# ── Callbacks ─────────────────────────────────────────────────────────────────

def register_pipeline(app) -> None:

    # ── 0. Clientside tab switcher ────────────────────────────────────────
    app.clientside_callback(
        """
        function(n_log, n_prev, n_stat, cur_tab) {
            const tid = dash_clientside.callback_context.triggered_id;
            let tab = cur_tab || 'log';
            if (tid === 'pipe-tab-btn-log')     tab = 'log';
            if (tid === 'pipe-tab-btn-preview') tab = 'preview';
            if (tid === 'pipe-tab-btn-status')  tab = 'status';
            const show  = {display: 'flex'};
            const hide  = {display: 'none'};
            const isLog  = tab === 'log';
            const isPrev = tab === 'preview';
            const isStat = tab === 'status';
            return [
                tab,
                isLog  ? show : hide,
                isPrev ? show : hide,
                isStat ? show : hide,
                isLog  ? 'pipe-tab-btn active' : 'pipe-tab-btn',
                isPrev ? 'pipe-tab-btn active' : 'pipe-tab-btn',
                isStat ? 'pipe-tab-btn active' : 'pipe-tab-btn',
            ];
        }
        """,
        Output(IDs.PIPE_ACTIVE_TAB,      "data"),
        Output("pipe-panel-log",         "style"),
        Output("pipe-panel-preview",     "style"),
        Output("pipe-panel-status",      "style"),
        Output("pipe-tab-btn-log",       "className"),
        Output("pipe-tab-btn-preview",   "className"),
        Output("pipe-tab-btn-status",    "className"),
        Input("pipe-tab-btn-log",        "n_clicks"),
        Input("pipe-tab-btn-preview",    "n_clicks"),
        Input("pipe-tab-btn-status",     "n_clicks"),
        State(IDs.PIPE_ACTIVE_TAB,       "data"),
        prevent_initial_call=True,
    )

    # ── 1. Sync step info + show/hide export folder ───────────────────────
    @app.callback(
        Output(IDs.PIPE_STEP_INFO,    "children"),
        Output(IDs.PIPE_METHOD,       "options"),
        Output(IDs.PIPE_METHOD,       "value"),
        Output(IDs.PIPE_EXPORT_CARD,  "style"),
        Input(IDs.PIPE_STEP,          "value"),
    )
    def sync_step_info(step_id):
        if step_id is None or step_id == "":
            raise PreventUpdate
        step_id = int(step_id)
        step = _STEPS[step_id]
        opts = [{"label": m.label, "value": m.name} for m in step.methods]
        export_style = (
            {"display": "block"} if step_id == 7
            else {"display": "none"}
        )
        return step.description, opts, step.default_method, export_style

    # ── 2. Run a single step ──────────────────────────────────────────────
    @app.callback(
        Output(IDs.PIPE_LOG,           "children"),
        Output(IDs.IMG_PIPE,           "src"),
        Output(IDs.PIPE_STORE,         "data"),
        Output(IDs.PIPE_VIEW_STEP,     "value"),
        Output(IDs.PIPE_STEP,          "value"),       # auto-advance
        Output(IDs.PIPE_SPINNER,       "children"),
        Output(IDs.TOAST_ERROR,        "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,         "children", allow_duplicate=True),
        Input(IDs.BTN_PIPE_RUN,        "n_clicks"),
        State(IDs.PIPE_STEP,           "value"),
        State(IDs.PIPE_METHOD,         "value"),
        State(IDs.PIPE_EXPORT_FOLDER,  "value"),
        State(IDs.SESSION_ID,          "data"),
        State(IDs.PIPE_STORE,          "data"),
        State(IDs.STORE_ACTIVE_LINES,  "data"),
        State(IDs.STORE_DATA,          "data"),
        State(IDs.STORE_THEME,         "data"),
        prevent_initial_call=True,
    )
    def run_step(n_clicks, step_id, method, export_folder, session_id,
                 preview_store, active_lines_store, store_data, theme):
        if not n_clicks:
            raise PreventUpdate

        dark = (theme or "dark") == "dark"
        log_lines = []
        store = dict(preview_store or {})

        try:
            sites = cache_get(session_id)

            sites, all_muted = _filter_by_active_lines(
                sites, store_data, active_lines_store
            )
            if all_muted:
                msg = "No active lines — unmute at least one line to run the pipeline."
                return msg, no_active_lines_src(dark), store, no_update, no_update, "", False, ""

            step_id = int(step_id)
            step = _CTRL.steps[step_id]
            step.active_method = method
            step.reset_params()

            # Inject export folder into params for step 7
            if step_id == 7 and export_folder:
                step.params["folder"] = export_folder.strip()

            _CTRL.set_input_sites(sites)
            _CTRL.execute_step(step_id, log_cb=lambda m: log_lines.append(m))

            result_sites = _CTRL._sites_chain[step_id]

            # Generate and cache the diagnostic preview for this step
            if step_id != 7:
                src = _generate_preview(step, result_sites, dark)
            else:
                src = empty_src(dark=dark)

            store[str(step_id)] = src

            # Auto-advance to next step (stay on last step)
            next_step = str(step_id + 1) if step_id < len(_STEPS) - 1 else str(step_id)

            return (
                "\n".join(log_lines),
                src,
                store,
                str(step_id),   # preview picker shows just-completed step
                next_step,      # PIPE_STEP advances to next
                "",
                False, "",
            )

        except Exception as exc:
            log_lines.append(f"ERROR: {exc}")
            plt.close("all")
            return (
                "\n".join(log_lines),
                empty_src(dark=dark),
                store,
                no_update,
                no_update,      # stay on current step on error
                "",
                True, str(exc),
            )

    # ── 3. Run all steps ──────────────────────────────────────────────────
    @app.callback(
        Output(IDs.PIPE_LOG,           "children",  allow_duplicate=True),
        Output(IDs.PIPE_STORE,         "data",      allow_duplicate=True),
        Output(IDs.PIPE_SPINNER,       "children",  allow_duplicate=True),
        Output(IDs.TOAST_ERROR,        "is_open",   allow_duplicate=True),
        Output(IDs.TOAST_BODY,         "children",  allow_duplicate=True),
        Input(IDs.BTN_PIPE_ALL,        "n_clicks"),
        State(IDs.PIPE_EXPORT_FOLDER,  "value"),
        State(IDs.SESSION_ID,          "data"),
        State(IDs.PIPE_STORE,          "data"),
        State(IDs.STORE_ACTIVE_LINES,  "data"),
        State(IDs.STORE_DATA,          "data"),
        State(IDs.STORE_THEME,         "data"),
        prevent_initial_call=True,
    )
    def run_all(n_clicks, export_folder, session_id,
                preview_store, active_lines_store, store_data, theme):
        if not n_clicks:
            raise PreventUpdate

        dark = (theme or "dark") == "dark"
        log_lines = []
        store = dict(preview_store or {})

        try:
            sites = cache_get(session_id)

            sites, all_muted = _filter_by_active_lines(
                sites, store_data, active_lines_store
            )
            if all_muted:
                return (
                    "No active lines — unmute at least one line.",
                    store, "", False, "",
                )

            _CTRL.set_input_sites(sites)

            for i in range(len(_CTRL.steps)):
                step = _CTRL.steps[i]
                # Inject export folder for step 7
                if i == 7 and export_folder:
                    step.active_method = step.default_method
                    step.reset_params()
                    step.params["folder"] = export_folder.strip()
                try:
                    _CTRL.execute_step(i, log_cb=lambda m: log_lines.append(m))
                    log_lines.append(f"── Step {i + 1} done ──")
                    # Cache preview for each step
                    result_sites = _CTRL._sites_chain[i]
                    if i != 7 and result_sites is not None:
                        try:
                            src = _generate_preview(step, result_sites, dark)
                            store[str(i)] = src
                        except Exception:
                            plt.close("all")
                except Exception as exc:
                    log_lines.append(f"Step {i + 1} ERROR: {exc}")
                    break

            return "\n".join(log_lines), store, "", False, ""

        except Exception as exc:
            return "\n".join(log_lines), store, "", True, str(exc)

    # ── 4. Browse per-step previews ───────────────────────────────────────
    @app.callback(
        Output(IDs.IMG_PIPE,       "src",       allow_duplicate=True),
        Input(IDs.PIPE_VIEW_STEP,  "value"),
        State(IDs.PIPE_STORE,      "data"),
        State(IDs.STORE_THEME,     "data"),
        prevent_initial_call=True,
    )
    def view_step_preview(step_id, preview_store, theme):
        if step_id is None or step_id == "":
            raise PreventUpdate
        store = preview_store or {}
        dark  = (theme or "dark") == "dark"
        src   = store.get(str(step_id))
        return src if src else empty_src(dark=dark)

    # ── 5. Skip a step ────────────────────────────────────────────────────
    @app.callback(
        Output(IDs.PIPE_LOG, "children", allow_duplicate=True),
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

    # ── 6. Reset pipeline ─────────────────────────────────────────────────
    @app.callback(
        Output(IDs.PIPE_LOG,      "children",  allow_duplicate=True),
        Output(IDs.PIPE_STATUS,   "children"),
        Output(IDs.PIPE_STORE,    "data",      allow_duplicate=True),
        *[Output(f"pipe-node-{s.id}", "className", allow_duplicate=True) for s in _STEPS],
        Input(IDs.BTN_PIPE_RESET, "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_pipeline(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        _CTRL.reset()
        node_classes = ["pipe-step-node step-active"] + [
            "pipe-step-node step-pending" for _ in _STEPS[1:]
        ]
        return ("Pipeline reset.", "", {}, *node_classes)

    _STEP_COLOUR = {
        StepStatus.PENDING: "#a6adc8",
        StepStatus.RUNNING: "#f9e2af",
        StepStatus.DONE:    "#a6e3a1",
        StepStatus.ERROR:   "#f38ba8",
        StepStatus.SKIPPED: "#cba6f7",
    }

    # ── 7. Refresh stepper nodes after each log update ────────────────────
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

    # ─────────────────────────────────────────────────────────────────────
    # Folder-browser modal callbacks
    # ─────────────────────────────────────────────────────────────────────

    def _list_dirs(path: str) -> list:
        """Return sorted list of subdirectory names under *path*."""
        try:
            entries = sorted(
                e for e in os.listdir(path)
                if not e.startswith(".")
                and os.path.isdir(os.path.join(path, e))
            )
        except PermissionError:
            return []
        return entries

    def _dir_listing(path: str) -> list:
        """Build dbc.ListGroup items for all subdirectories at *path*."""
        dirs = _list_dirs(path)
        if not dirs:
            return [
                html.Div(
                    "No sub-folders here.",
                    style={"color": "#6c7086", "fontSize": "11px",
                           "padding": "8px"},
                )
            ]
        items = []
        for name in dirs:
            full = os.path.join(path, name)
            items.append(
                dbc.ListGroupItem(
                    [
                        html.I(className="bi bi-folder-fill me-2",
                               style={"color": "#f9e2af"}),
                        name,
                    ],
                    id={"type": "pipe-dir-item", "index": full},
                    action=True,
                    style={"cursor": "pointer", "fontSize": "12px",
                           "padding": "6px 10px"},
                    n_clicks=0,
                )
            )
        return [dbc.ListGroup(items, flush=True)]

    # ── 8a. Open browse modal ─────────────────────────────────────────────
    @app.callback(
        Output(IDs.PIPE_BROWSE_MODAL, "is_open"),
        Output(IDs.PIPE_BROWSE_PATH,  "data"),
        Input(IDs.BTN_PIPE_BROWSE,    "n_clicks"),
        Input("pipe-browse-cancel",   "n_clicks"),
        State(IDs.PIPE_EXPORT_FOLDER, "value"),
        State(IDs.PIPE_BROWSE_MODAL,  "is_open"),
        prevent_initial_call=True,
    )
    def toggle_browse_modal(n_open, n_cancel, current_folder, is_open):
        trigger = ctx.triggered_id
        if trigger == "pipe-browse-cancel":
            return False, no_update
        # Open: start from current folder's parent if it exists, else home
        start = os.path.expanduser("~")
        if current_folder:
            cf = current_folder.strip()
            candidate = cf if os.path.isdir(cf) else os.path.dirname(cf)
            if os.path.isdir(candidate):
                start = candidate
        return True, start

    # ── 8b. Navigate: click a folder item ────────────────────────────────
    @app.callback(
        Output(IDs.PIPE_BROWSE_PATH, "data", allow_duplicate=True),
        Input({"type": "pipe-dir-item", "index": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def navigate_into(n_clicks_list):
        if not any(n_clicks_list):
            raise PreventUpdate
        triggered = ctx.triggered_id
        if triggered and isinstance(triggered, dict):
            new_path = triggered.get("index", "")
            if new_path and os.path.isdir(new_path):
                return new_path
        raise PreventUpdate

    # ── 8c. Navigate up to parent ─────────────────────────────────────────
    @app.callback(
        Output(IDs.PIPE_BROWSE_PATH, "data", allow_duplicate=True),
        Input(IDs.BTN_PIPE_BROWSE_UP, "n_clicks"),
        State(IDs.PIPE_BROWSE_PATH,   "data"),
        prevent_initial_call=True,
    )
    def navigate_up(n_clicks, current_path):
        if not n_clicks or not current_path:
            raise PreventUpdate
        parent = os.path.dirname(current_path)
        if parent and parent != current_path:
            return parent
        raise PreventUpdate

    # ── 8d. Refresh listing when path changes ─────────────────────────────
    @app.callback(
        Output(IDs.PIPE_BROWSE_LIST, "children"),
        Output(IDs.PIPE_BROWSE_CWD,  "children"),
        Input(IDs.PIPE_BROWSE_PATH,  "data"),
    )
    def refresh_listing(path):
        if not path:
            path = os.path.expanduser("~")
        return _dir_listing(path), path

    # ── 8e. Create new folder ─────────────────────────────────────────────
    @app.callback(
        Output(IDs.PIPE_BROWSE_PATH, "data",     allow_duplicate=True),
        Output(IDs.PIPE_BROWSE_NEW,  "value"),
        Input(IDs.BTN_PIPE_BROWSE_MK, "n_clicks"),
        State(IDs.PIPE_BROWSE_NEW,    "value"),
        State(IDs.PIPE_BROWSE_PATH,   "data"),
        prevent_initial_call=True,
    )
    def make_dir(n_clicks, name, current_path):
        if not n_clicks or not name or not current_path:
            raise PreventUpdate
        name = name.strip().replace("/", "_").replace("\\", "_")
        new_dir = os.path.join(current_path, name)
        try:
            os.makedirs(new_dir, exist_ok=True)
        except Exception:
            raise PreventUpdate
        return new_dir, ""     # navigate into the new folder and clear input

    # ── 8f. Confirm selection → populate export folder + close modal ──────
    @app.callback(
        Output(IDs.PIPE_EXPORT_FOLDER, "value",   allow_duplicate=True),
        Output(IDs.PIPE_BROWSE_MODAL,  "is_open", allow_duplicate=True),
        Input(IDs.BTN_PIPE_BROWSE_SEL, "n_clicks"),
        State(IDs.PIPE_BROWSE_PATH,    "data"),
        prevent_initial_call=True,
    )
    def select_folder(n_clicks, path):
        if not n_clicks or not path:
            raise PreventUpdate
        return path, False
