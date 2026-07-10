# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Help / About modal toggle."""

from __future__ import annotations

from dash import ALL, Input, Output, State, ctx
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

    # Click an example prompt → drop it into the chat input and close
    # the Help modal so the user can review / edit / send it.
    @app.callback(
        Output(IDs.INPUT, "value", allow_duplicate=True),
        Output(
            IDs.MODAL_HELP, "is_open",
            allow_duplicate=True,
        ),
        Input(
            {"type": IDs.HELP_CHIP, "index": ALL},
            "n_clicks",
        ),
        prevent_initial_call=True,
    )
    def use_example(_clicks):
        if not ctx.triggered or not ctx.triggered[0].get(
            "value"
        ):
            raise PreventUpdate
        trig = ctx.triggered_id
        if not isinstance(trig, dict):
            raise PreventUpdate
        from ..layout import _HELP_EXAMPLES
        try:
            text = _HELP_EXAMPLES[trig["index"]]
        except (IndexError, KeyError, TypeError):
            raise PreventUpdate
        return text, False
