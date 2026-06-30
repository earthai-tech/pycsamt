# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Chrome callbacks: seed adoption, theme, help, sidebar, dock."""

from __future__ import annotations

from dash import Input, Output, State, no_update

from .._ids import IDs
from ..cache import set_view, take_seed
from .._render import store_from_view


def register_chrome(app) -> None:
    _register_seed(app)
    _register_theme(app)
    _register_help(app)
    _register_sidebar(app)
    _register_inspector(app)
    _register_dock(app)


# Plotly's responsive mode only re-lays out on a window resize event, so we
# fire one after toggling a panel to make the canvas reclaim/yield space.
_RESIZE = (
    "setTimeout(function(){"
    "window.dispatchEvent(new Event('resize'));}, 80);"
)


def _register_seed(app) -> None:
    """Adopt a view handed in by ``MapView.launch`` on first load."""
    @app.callback(
        Output(IDs.STORE_DATA, "data"),
        Output(IDs.DATA_BADGE_TEXT, "children"),
        Output(IDs.DATA_BADGE, "className"),
        Input(IDs.SESSION_ID, "data"),
        prevent_initial_call=False,
    )
    def adopt_seed(session_id):
        view = take_seed()
        if view is None or not session_id:
            return no_update, no_update, no_update
        set_view(session_id, view)
        store = store_from_view(view)
        badge = f"{store['n_stations']} stations · {store['n_lines']} line(s)"
        return store, badge, "mv-data-badge visible"


def _register_theme(app) -> None:
    @app.callback(
        Output(IDs.STORE_THEME, "data"),
        Output("mv-theme-icon", "className"),
        Input(IDs.BTN_THEME, "n_clicks"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def toggle_theme(_n, theme):
        new = "dark" if (theme or "light") == "light" else "light"
        icon = "bi bi-sun" if new == "dark" else "bi bi-moon-stars"
        return new, icon


def _register_help(app) -> None:
    @app.callback(
        Output(IDs.MODAL_HELP, "is_open"),
        Input(IDs.BTN_HELP, "n_clicks"),
        Input(IDs.BTN_HELP_CLOSE, "n_clicks"),
        State(IDs.MODAL_HELP, "is_open"),
        prevent_initial_call=True,
    )
    def toggle_help(_o, _c, is_open):
        return not is_open


def _register_sidebar(app) -> None:
    app.clientside_callback(
        """
        function(n) {
            var hidden = (n || 0) % 2 === 1;
            """ + _RESIZE + """
            return hidden ? 'mv-datapanel mv-datapanel--hidden'
                          : 'mv-datapanel';
        }
        """,
        Output("mv-datapanel", "className"),
        Input(IDs.BTN_SIDEBAR, "n_clicks"),
        prevent_initial_call=True,
    )


def _register_inspector(app) -> None:
    app.clientside_callback(
        """
        function(n) {
            var collapsed = (n || 0) % 2 === 1;
            """ + _RESIZE + """
            return collapsed ? 'mv-inspector mv-inspector--collapsed'
                             : 'mv-inspector';
        }
        """,
        Output(IDs.INSPECTOR, "className"),
        Input(IDs.BTN_INSPECTOR, "n_clicks"),
        prevent_initial_call=True,
    )


def _register_dock(app) -> None:
    app.clientside_callback(
        """
        function(n) {
            var open = (n || 0) % 2 === 1;
            """ + _RESIZE + """
            return [
                open ? {display:'block'} : {display:'none'},
                open ? 'bi bi-chevron-down ms-2' : 'bi bi-chevron-up ms-2'
            ];
        }
        """,
        Output(IDs.DOCK_BODY, "style"),
        Output("mv-dock-chevron", "className"),
        Input(IDs.DOCK_TOGGLE, "n_clicks"),
        prevent_initial_call=True,
    )
