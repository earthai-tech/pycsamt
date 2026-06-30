# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Canvas toolbar (Fit / layer toggles / basemap / markers) and the
coordinate-system panel.

The toolbar mirrors the inspector inputs — clicking a toolbar button
flips the matching inspector control value, so there is a single source
of truth and the two stay in sync.
"""

from __future__ import annotations

from dash import ctx, Input, Output, State, no_update, html

from .._ids import IDs

_BM_BTN_STYLE = {
    IDs.TB_BM_DARK: "carto-darkmatter",
    IDs.TB_BM_LIGHT: "carto-positron",
    IDs.TB_BM_SAT: "esri-satellite",
    IDs.TB_BM_STREET: "esri-street",
    IDs.TB_BM_TOPO: "esri-topo",
}


def register_toolbar(app) -> None:
    _register_info(app)
    _register_fit(app)
    _register_toggles(app)
    _register_active_state(app)
    _register_basemap_quickswitch(app)
    _register_markers(app)
    _register_crs(app)


def _register_info(app) -> None:
    @app.callback(
        Output(IDs.TB_INFO, "children"),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=False,
    )
    def info(store):
        if not store or not store.get("n_stations"):
            return ""
        n, m = store.get("n_stations", 0), store.get("n_lines", 0)
        return f"{n} stations · {m} line(s)"


def _register_fit(app) -> None:
    app.clientside_callback(
        "function(n){ return n || 0; }",
        Output(IDs.STORE_FIT, "data"),
        Input(IDs.TB_FIT, "n_clicks"),
        prevent_initial_call=True,
    )


def _register_toggles(app) -> None:
    """Toolbar Labels/Profiles/Contour flip the matching inspector switch."""
    for btn, target in (
        (IDs.TB_LABELS, IDs.CTL_LABELS),
        (IDs.TB_PROFILES, IDs.CTL_PROFILES),
        (IDs.TB_CONTOUR, IDs.CTL_CONTOUR_ENABLE),
    ):
        app.clientside_callback(
            "function(n, cur){ if(!n) return window.dash_clientside.no_update;"
            " return !cur; }",
            Output(target, "value", allow_duplicate=True),
            Input(btn, "n_clicks"),
            State(target, "value"),
            prevent_initial_call=True,
        )


def _register_active_state(app) -> None:
    """Reflect the switch value as the toolbar button's active class."""
    for btn, target in (
        (IDs.TB_LABELS, IDs.CTL_LABELS),
        (IDs.TB_PROFILES, IDs.CTL_PROFILES),
        (IDs.TB_CONTOUR, IDs.CTL_CONTOUR_ENABLE),
    ):
        app.clientside_callback(
            "function(v){ return v ? 'mv-tb-btn active' : 'mv-tb-btn'; }",
            Output(btn, "className"),
            Input(target, "value"),
            prevent_initial_call=False,
        )


def _register_basemap_quickswitch(app) -> None:
    @app.callback(
        Output(IDs.CTL_BASEMAP, "value"),
        Input(IDs.TB_BM_DARK, "n_clicks"),
        Input(IDs.TB_BM_LIGHT, "n_clicks"),
        Input(IDs.TB_BM_SAT, "n_clicks"),
        Input(IDs.TB_BM_STREET, "n_clicks"),
        Input(IDs.TB_BM_TOPO, "n_clicks"),
        prevent_initial_call=True,
    )
    def quickswitch(*_clicks):
        return _BM_BTN_STYLE.get(ctx.triggered_id, no_update)


def _register_markers(app) -> None:
    app.clientside_callback(
        """
        function(dec, inc, cur) {
            var v = cur || 10;
            var t = window.dash_clientside.callback_context.triggered;
            if (!t.length) return window.dash_clientside.no_update;
            var id = t[0].prop_id.split('.')[0];
            if (id.indexOf('mark-inc') >= 0) v = Math.min(24, v + 2);
            else if (id.indexOf('mark-dec') >= 0) v = Math.max(4, v - 2);
            return v;
        }
        """,
        Output(IDs.CTL_MARKER_SIZE, "value", allow_duplicate=True),
        Input(IDs.TB_MARK_DEC, "n_clicks"),
        Input(IDs.TB_MARK_INC, "n_clicks"),
        State(IDs.CTL_MARKER_SIZE, "value"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        "function(v){ return String(v || 10); }",
        Output(IDs.TB_MARK_VAL, "children"),
        Input(IDs.CTL_MARKER_SIZE, "value"),
        prevent_initial_call=False,
    )


def _register_crs(app) -> None:
    # show/hide UTM vs EPSG inputs by mode
    app.clientside_callback(
        """
        function(mode) {
            var m = mode || 'geo';
            return [
                m === 'utm'    ? {display:'block'} : {display:'none'},
                m === 'custom' ? {display:'block'} : {display:'none'}
            ];
        }
        """,
        Output(IDs.GRP_UTM, "style"),
        Output(IDs.GRP_EPSG, "style"),
        Input(IDs.CTL_CRS_MODE, "value"),
        prevent_initial_call=False,
    )

    @app.callback(
        Output(IDs.CRS_INFO, "children"),
        Input(IDs.CTL_CRS_MODE, "value"),
        Input(IDs.CTL_UTM_ZONE, "value"),
        Input(IDs.CTL_UTM_HEM, "value"),
        Input(IDs.CTL_EPSG, "value"),
        prevent_initial_call=False,
    )
    def crs_info(mode, zone, hem, epsg):
        from pycsamt.map import resolve_crs_info

        mode = mode or "geo"
        try:
            text = resolve_crs_info(
                mode,
                zone=int(zone or 50),
                hemisphere=hem or "N",
                epsg=epsg or 4326,
            )
        except Exception:
            text = "Invalid CRS settings"
        note = (
            " — coordinates are reprojected to lon/lat for the basemap."
            if mode != "geo" else ""
        )
        return html.Span([html.I(className="bi bi-info-circle me-1"),
                          text + note])
