# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Export the current view as a standalone interactive HTML file."""

from __future__ import annotations

from dash import Input, Output, State, no_update

from .._ids import IDs
from ..cache import get_view
from .._render import figure_for


def register_export(app) -> None:
    @app.callback(
        Output(IDs.EXPORT_DL, "data"),
        Input(IDs.BTN_EXPORT, "n_clicks"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_VIEW, "data"),
        State(IDs.STORE_CONTROLS, "data"),
        State(IDs.STORE_THEME, "data"),
        State(IDs.STORE_LINES, "data"),
        prevent_initial_call=True,
    )
    def export(_n, session_id, view_name, controls, theme, lines):
        view = get_view(session_id)
        if view is None:
            return no_update
        fig = figure_for(
            view_name or "map",
            view,
            controls,
            theme=theme or "light",
            active_lines=(lines or {}).get("active"),
        )
        if not hasattr(fig, "to_html"):
            return no_update
        html_str = fig.to_html(include_plotlyjs="cdn", full_html=True)
        return {
            "content": html_str,
            "filename": f"mapview_{view_name or 'station'}.html",
        }
