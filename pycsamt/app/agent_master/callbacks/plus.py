# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the "+" resource/action menu."""

from __future__ import annotations

from dash import Input, Output
from dash.exceptions import PreventUpdate

from .._ids import IDs


def register_plus(app) -> None:

    @app.callback(
        Output(
            IDs.CANVAS_EDI,
            "is_open",
            allow_duplicate=True,
        ),
        Input(IDs.PLUS_LOAD, "n_clicks"),
        prevent_initial_call=True,
    )
    def _open_edi(n):
        if not n:
            raise PreventUpdate
        return True

    @app.callback(
        Output(
            IDs.CANVAS_SETTINGS,
            "is_open",
            allow_duplicate=True,
        ),
        Input(IDs.PLUS_SETTINGS, "n_clicks"),
        prevent_initial_call=True,
    )
    def _open_settings(n):
        if not n:
            raise PreventUpdate
        return True

    @app.callback(
        Output(
            IDs.INPUT,
            "value",
            allow_duplicate=True,
        ),
        Input(IDs.PLUS_PATH, "n_clicks"),
        prevent_initial_call=True,
    )
    def _paste_path(n):
        if not n:
            raise PreventUpdate
        return "Load /path/to/your/EDIs, then "

    @app.callback(
        Output(
            IDs.INPUT,
            "value",
            allow_duplicate=True,
        ),
        Input(IDs.PLUS_WEB, "n_clicks"),
        prevent_initial_call=True,
    )
    def _inject_web_launch(n):
        if not n:
            raise PreventUpdate
        return "open web app"
