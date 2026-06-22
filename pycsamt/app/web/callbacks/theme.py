# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""callbacks/theme.py — dark / light theme toggle."""

from __future__ import annotations

from dash import Input, Output, State, no_update

from pycsamt.app.web.layout import IDs


def register_theme(app) -> None:
    @app.callback(
        Output(IDs.STORE_THEME, "data"),
        Output(IDs.BTN_THEME,   "children"),
        Input(IDs.BTN_THEME,    "n_clicks"),
        State(IDs.STORE_THEME,  "data"),
        prevent_initial_call=True,
    )
    def toggle_theme(_n_clicks, current_theme):
        new_theme = "light" if (current_theme or "dark") == "dark" else "dark"
        label     = "☀ Theme" if new_theme == "light" else "☾ Theme"
        return new_theme, label
