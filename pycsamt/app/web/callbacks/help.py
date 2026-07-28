# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/help.py — Help menu: About modal toggle.
"""

from __future__ import annotations

from dash import Input, Output, State, no_update

from pycsamt.app.web.layout import IDs


def register_help(app) -> None:
    @app.callback(
        Output(IDs.ABOUT_MODAL, "is_open"),
        Input(IDs.BTN_ABOUT, "n_clicks"),
        Input("about-close-btn", "n_clicks"),
        State(IDs.ABOUT_MODAL, "is_open"),
        prevent_initial_call=True,
    )
    def _toggle_about(n_open, n_close, is_open):
        if not (n_open or n_close):
            return no_update
        triggered = [
            t["prop_id"] for t in __import__("dash").callback_context.triggered
        ]
        if any("about-close-btn" in t for t in triggered):
            return False
        return True
