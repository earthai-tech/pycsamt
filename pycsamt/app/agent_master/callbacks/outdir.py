# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Output directory browser modal callbacks."""

from __future__ import annotations

import os
from pathlib import Path

from dash import (
    ALL,
    Input,
    Output,
    State,
    ctx,
    html,
    no_update,
)
from dash.exceptions import PreventUpdate

from .._ids import IDs


def _list_dirs(path: str) -> list[str]:
    """Return sorted subdirectory names in path."""
    try:
        entries = sorted(
            e
            for e in os.listdir(path)
            if os.path.isdir(os.path.join(path, e)) and not e.startswith(".")
        )
        return entries
    except PermissionError:
        return []


def _dir_listing(path: str) -> list:
    """Build directory listing children."""
    dirs = _list_dirs(path)
    if not dirs:
        return [
            html.Div(
                html.Em(
                    "No subdirectories",
                    style={
                        "fontSize": "12px",
                        "color": "var(--am-muted)",
                    },
                )
            )
        ]
    return [
        html.Button(
            [
                html.I(
                    className=("bi bi-folder me-2"),
                    style={
                        "color": "#e8a020",
                    },
                ),
                html.Span(
                    name,
                    style={"fontSize": "13px"},
                ),
            ],
            id={
                "type": "am-outdir-entry",
                "index": i,
                "name": name,
            },
            className="am-outdir-entry",
            n_clicks=0,
            style={
                "display": "block",
                "width": "100%",
                "textAlign": "left",
                "background": "none",
                "border": "none",
                "padding": "4px 8px",
                "cursor": "pointer",
                "borderRadius": "4px",
            },
        )
        for i, name in enumerate(dirs)
    ]


def register_outdir(app) -> None:

    # 1. Open modal on browse button click
    @app.callback(
        Output(IDs.MODAL_OUTPUT_BROWSE, "is_open"),
        Output(IDs.OUTPUT_BROWSE_STORE, "data"),
        Output(IDs.OUTPUT_BROWSE_PATH, "children"),
        Output(IDs.OUTPUT_BROWSE_LIST, "children"),
        Input(IDs.BTN_OUTPUT_BROWSE, "n_clicks"),
        State(IDs.OUTPUT_DIR, "value"),
        prevent_initial_call=True,
    )
    def open_output_modal(n, current_val):
        if not n:
            raise PreventUpdate
        # Start from current value, EDI path,
        # or home directory
        start = current_val or os.path.expanduser("~")
        if not os.path.isdir(start):
            start = os.path.expanduser("~")
        return (
            True,
            {"path": start},
            start,
            _dir_listing(start),
        )

    # 2. Navigate into a subdirectory
    @app.callback(
        Output(
            IDs.OUTPUT_BROWSE_STORE,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.OUTPUT_BROWSE_PATH,
            "children",
            allow_duplicate=True,
        ),
        Output(
            IDs.OUTPUT_BROWSE_LIST,
            "children",
            allow_duplicate=True,
        ),
        Output(
            IDs.OUTPUT_BROWSE_STATUS,
            "children",
            allow_duplicate=True,
        ),
        Input(
            {"type": "am-outdir-entry", "index": ALL, "name": ALL},
            "n_clicks",
        ),
        State(IDs.OUTPUT_BROWSE_STORE, "data"),
        prevent_initial_call=True,
    )
    def enter_subdir(n_clicks, store):
        if not any(n_clicks or []):
            raise PreventUpdate
        triggered = ctx.triggered_id
        if not isinstance(triggered, dict):
            raise PreventUpdate
        name = triggered.get("name", "")
        cur = (store or {}).get("path", os.path.expanduser("~"))
        new_path = os.path.join(cur, name)
        if not os.path.isdir(new_path):
            raise PreventUpdate
        return (
            {"path": new_path},
            new_path,
            _dir_listing(new_path),
            "",
        )

    # 3. Navigate up one level
    @app.callback(
        Output(
            IDs.OUTPUT_BROWSE_STORE,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.OUTPUT_BROWSE_PATH,
            "children",
            allow_duplicate=True,
        ),
        Output(
            IDs.OUTPUT_BROWSE_LIST,
            "children",
            allow_duplicate=True,
        ),
        Input(IDs.BTN_OUTPUT_UP, "n_clicks"),
        State(IDs.OUTPUT_BROWSE_STORE, "data"),
        prevent_initial_call=True,
    )
    def go_up(n, store):
        if not n:
            raise PreventUpdate
        cur = (store or {}).get("path", os.path.expanduser("~"))
        parent = str(Path(cur).parent)
        if parent == cur:
            raise PreventUpdate
        return (
            {"path": parent},
            parent,
            _dir_listing(parent),
        )

    # 4. Create new subfolder
    @app.callback(
        Output(
            IDs.OUTPUT_BROWSE_STORE,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.OUTPUT_BROWSE_PATH,
            "children",
            allow_duplicate=True,
        ),
        Output(
            IDs.OUTPUT_BROWSE_LIST,
            "children",
            allow_duplicate=True,
        ),
        Output(
            IDs.OUTPUT_BROWSE_STATUS,
            "children",
            allow_duplicate=True,
        ),
        Output(
            IDs.OUTPUT_MKDIR_INPUT,
            "value",
            allow_duplicate=True,
        ),
        Input(IDs.BTN_OUTPUT_MKDIR, "n_clicks"),
        State(IDs.OUTPUT_MKDIR_INPUT, "value"),
        State(IDs.OUTPUT_BROWSE_STORE, "data"),
        prevent_initial_call=True,
    )
    def mkdir(n, name, store):
        if not n or not (name or "").strip():
            raise PreventUpdate
        cur = (store or {}).get("path", os.path.expanduser("~"))
        name = name.strip()
        new_dir = os.path.join(cur, name)
        try:
            Path(new_dir).mkdir(parents=True, exist_ok=True)
            status = html.Span(
                f"Created: {name}",
                style={"color": "#2f9e44"},
            )
        except Exception as exc:
            return (
                no_update,
                no_update,
                no_update,
                html.Span(
                    str(exc),
                    style={"color": "red"},
                ),
                name,
            )
        return (
            {"path": new_dir},
            new_dir,
            _dir_listing(new_dir),
            status,
            "",
        )

    # 5. Confirm: set OUTPUT_DIR and close
    @app.callback(
        Output(IDs.OUTPUT_DIR, "value"),
        Output(
            IDs.MODAL_OUTPUT_BROWSE,
            "is_open",
            allow_duplicate=True,
        ),
        Input(IDs.BTN_OUTPUT_CONFIRM, "n_clicks"),
        State(IDs.OUTPUT_BROWSE_STORE, "data"),
        prevent_initial_call=True,
    )
    def confirm_output_dir(n, store):
        if not n:
            raise PreventUpdate
        path = (store or {}).get("path", "")
        return path, False
