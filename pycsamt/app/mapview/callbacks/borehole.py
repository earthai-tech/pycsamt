# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Standalone PCBH upload callback for Map View."""

from __future__ import annotations

from dash import Input, Output, State, html
from dash.exceptions import PreventUpdate

from pycsamt.app._borehole import decode_pcbh_upload

from .._ids import IDs


def register_borehole(app) -> None:
    @app.callback(
        Output(IDs.PCBH_STORE, "data"),
        Output(IDs.PCBH_UPLOAD_INFO, "children"),
        Input(IDs.PCBH_UPLOAD, "contents"),
        State(IDs.PCBH_UPLOAD, "filename"),
        prevent_initial_call=True,
    )
    def load_pcbh(contents, filename):
        if not contents:
            raise PreventUpdate
        try:
            store = decode_pcbh_upload(contents, filename)
        except (TypeError, ValueError) as error:
            return None, html.Span(str(error), className="text-danger")
        count = store["n_boreholes"]
        return store, html.Span(
            f"{store['filename']} — {count} borehole(s)",
            className="text-success",
        )
