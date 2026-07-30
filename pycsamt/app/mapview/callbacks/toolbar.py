# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Canvas toolbar (Fit / layer toggles / basemap / markers) and the
coordinate-system panel.

The toolbar mirrors the inspector inputs — clicking a toolbar button
flips the matching inspector control value, so there is a single source
of truth and the two stay in sync.
"""

from __future__ import annotations

from dash import Input, Output, State, ctx, html, no_update

from .._ids import IDs

_BM_BTN_STYLE = {
    IDs.TB_BM_DARK: "carto-darkmatter",
    IDs.TB_BM_LIGHT: "carto-positron",
    IDs.TB_BM_SAT: "esri-satellite",
    IDs.TB_BM_STREET: "esri-street",
    IDs.TB_BM_TOPO: "esri-topo",
}


_MODE3D_BTN = {
    IDs.TB3D_MODE_FENCE: "fence",
    IDs.TB3D_MODE_BLOCK: "block",
    IDs.TB3D_MODE_DEPTH: "depth",
    IDs.TB3D_MODE_SURFACE: "surface",
}


def register_toolbar(app) -> None:
    _register_info(app)
    _register_fit(app)
    _register_toggles(app)
    _register_active_state(app)
    _register_basemap_quickswitch(app)
    _register_markers(app)
    _register_crs(app)
    _register_view_visibility(app)
    _register_mode3d_quick(app)
    _register_topo_quick(app)


def _register_view_visibility(app) -> None:
    """Show the 2-D toolbar in map mode, the 3-D toolbar in map3d mode —
    the two are mutually exclusive (mirrors
    ``controls._register_group_visibility`` for the inspector panels)."""
    app.clientside_callback(
        """
        function(view) {
            var v = view || 'map';
            var show = {display:'flex'}, hide = {display:'none'};
            return [
                v === 'map'   ? show : hide,
                v === 'map3d' ? show : hide
            ];
        }
        """,
        Output(IDs.TOOLBAR_2D, "style"),
        Output(IDs.TOOLBAR_3D, "style"),
        Input(IDs.STORE_VIEW, "data"),
        prevent_initial_call=False,
    )


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
        "function(n1, n2){ return (n1 || 0) + (n2 || 0); }",
        Output(IDs.STORE_FIT, "data"),
        Input(IDs.TB_FIT, "n_clicks"),
        Input(IDs.TB3D_RESET, "n_clicks"),
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


def _register_mode3d_quick(app) -> None:
    """3-D toolbar mode buttons write CTL_MODE3D (same select the Inspector's
    'Mode & quantity' accordion section uses) and reflect the active mode
    via className."""

    @app.callback(
        Output(IDs.CTL_MODE3D, "value"),
        Input(IDs.TB3D_MODE_FENCE, "n_clicks"),
        Input(IDs.TB3D_MODE_BLOCK, "n_clicks"),
        Input(IDs.TB3D_MODE_DEPTH, "n_clicks"),
        Input(IDs.TB3D_MODE_SURFACE, "n_clicks"),
        prevent_initial_call=True,
    )
    def quickswitch(*_clicks):
        return _MODE3D_BTN.get(ctx.triggered_id, no_update)

    app.clientside_callback(
        """
        function(mode, fence_id, block_id, depth_id, surface_id) {
            var m = mode || 'fence';
            var map = {};
            map[fence_id] = 'fence';
            map[block_id] = 'block';
            map[depth_id] = 'depth';
            map[surface_id] = 'surface';
            var ids = [fence_id, block_id, depth_id, surface_id];
            return ids.map(function(id){
                return map[id] === m ? 'mv-tb-btn active' : 'mv-tb-btn';
            });
        }
        """,
        Output(IDs.TB3D_MODE_FENCE, "className"),
        Output(IDs.TB3D_MODE_BLOCK, "className"),
        Output(IDs.TB3D_MODE_DEPTH, "className"),
        Output(IDs.TB3D_MODE_SURFACE, "className"),
        Input(IDs.CTL_MODE3D, "value"),
        State(IDs.TB3D_MODE_FENCE, "id"),
        State(IDs.TB3D_MODE_BLOCK, "id"),
        State(IDs.TB3D_MODE_DEPTH, "id"),
        State(IDs.TB3D_MODE_SURFACE, "id"),
        prevent_initial_call=False,
    )


def _register_topo_quick(app) -> None:
    """3-D toolbar Topo button flips CTL_TOPO (same switch the Inspector's
    'Topography' accordion section uses)."""
    app.clientside_callback(
        "function(n, cur){ if(!n) return window.dash_clientside.no_update;"
        " return !cur; }",
        Output(IDs.CTL_TOPO, "value", allow_duplicate=True),
        Input(IDs.TB3D_TOPO, "n_clicks"),
        State(IDs.CTL_TOPO, "value"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        "function(v){ return v ? 'mv-tb-btn active' : 'mv-tb-btn'; }",
        Output(IDs.TB3D_TOPO, "className"),
        Input(IDs.CTL_TOPO, "value"),
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
            " — station coordinates are shown in this system "
            "in the inspector and table (the basemap stays lon/lat)."
            if mode != "geo"
            else ""
        )
        return html.Span(
            [html.I(className="bi bi-info-circle me-1"), text + note]
        )
