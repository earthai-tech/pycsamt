# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""View-rail switching and the main canvas figure render."""

from __future__ import annotations

from dash import Input, Output, State, ctx

from .._ids import IDs
from .._render import VIEW_TITLES, empty_figure, figure_for
from ..cache import get_view

_RAIL = {
    IDs.RAIL_MAP: "map",
    IDs.RAIL_3D: "map3d",
}


def register_view(app) -> None:
    _register_rail(app)
    _register_render(app)


def _register_rail(app) -> None:
    @app.callback(
        Output(IDs.STORE_VIEW, "data"),
        Output(IDs.RAIL_MAP, "className"),
        Output(IDs.RAIL_3D, "className"),
        Output(IDs.CANVAS_TITLE, "children"),
        Input(IDs.RAIL_MAP, "n_clicks"),
        Input(IDs.RAIL_3D, "n_clicks"),
        prevent_initial_call=True,
    )
    def switch(*_clicks):
        active_id = ctx.triggered_id
        view = _RAIL.get(active_id, "map")

        def cls(bid):
            return "mv-rail-btn mv-rail-active" if bid == active_id else "mv-rail-btn"

        return (
            view,
            cls(IDs.RAIL_MAP),
            cls(IDs.RAIL_3D),
            VIEW_TITLES.get(view, "Map"),
        )


def _register_render(app) -> None:
    @app.callback(
        Output(IDs.CANVAS_GRAPH, "figure"),
        Output(IDs.WELCOME, "style"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_VIEW, "data"),
        Input(IDs.STORE_CONTROLS, "data"),
        Input(IDs.STORE_THEME, "data"),
        Input(IDs.STORE_LINES, "data"),
        Input(IDs.STORE_FIT, "data"),
        Input(IDs.STORE_MASKED, "data"),
        State(IDs.SESSION_ID, "data"),
        prevent_initial_call=False,
    )
    def render(store, view_name, controls, theme, lines, fit, masked, session_id):
        theme = theme or "light"
        if not store or not store.get("n_stations"):
            return empty_figure(theme), {"display": "flex"}
        view = get_view(session_id)
        if view is None:
            return empty_figure(theme, "Session data unavailable — reload lines."), {
                "display": "flex"
            }
        active = (lines or {}).get("active")
        return (
            figure_for(
                view_name or "map",
                view,
                controls,
                theme=theme,
                active_lines=active,
                masked=masked,
                fit=int(fit or 0),
            ),
            {"display": "none"},
        )
