# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Help / About modal toggle."""

from __future__ import annotations

from dash import Input, Output, State
from dash.exceptions import PreventUpdate

from .._ids import IDs


def register_help(app) -> None:

    @app.callback(
        Output(IDs.MODAL_HELP, "is_open"),
        Input(IDs.BTN_HELP, "n_clicks"),
        Input(IDs.BTN_HELP_CLOSE, "n_clicks"),
        State(IDs.MODAL_HELP, "is_open"),
        prevent_initial_call=True,
    )
    def toggle_help(_open, _close, is_open):
        if not (_open or _close):
            raise PreventUpdate
        return not is_open
