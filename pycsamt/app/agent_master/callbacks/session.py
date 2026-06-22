# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Session save / load callbacks."""

from __future__ import annotations

import json
import time
from datetime import datetime
from pathlib import Path

from dash import Input, Output, State
from dash import dcc, html, no_update
from dash.exceptions import PreventUpdate

from .._ids import IDs

_SESSION_DIR = (
    Path.home() / ".local" / "share"
    / "pycsamt" / "sessions"
)


def _session_path(name: str = "") -> Path:
    _SESSION_DIR.mkdir(
        parents=True, exist_ok=True
    )
    if not name:
        ts = datetime.now().strftime(
            "%Y%m%d_%H%M%S"
        )
        name = f"session_{ts}.json"
    return _SESSION_DIR / name


def register_session(app) -> None:

    @app.callback(
        Output(IDs.SAVE_STATUS, "children"),
        Input(IDs.BTN_SAVE_SESSION, "n_clicks"),
        State(IDs.STORE_MESSAGES, "data"),
        State(IDs.STORE_EDI, "data"),
        State(IDs.STORE_SETTINGS, "data"),
        prevent_initial_call=True,
    )
    def save_session(
        n, messages, edi_store, settings
    ):
        if not n:
            raise PreventUpdate
        payload = {
            "saved_at": datetime.now().isoformat(),
            "edi": edi_store or {},
            "settings": {
                k: v
                for k, v in (
                    settings or {}
                ).items()
                if not k.startswith("key_")
            },
            "message_count": len(
                messages or []
            ),
        }
        path = _session_path()
        path.write_text(
            json.dumps(payload, indent=2)
        )
        return html.Span(
            [
                html.I(
                    className=(
                        "bi bi-check-circle-fill"
                        " me-1"
                    )
                ),
                f"Saved to {path.name}",
            ],
            style={
                "fontSize": "11px",
                "color": "var(--tag-ok)",
            },
        )
