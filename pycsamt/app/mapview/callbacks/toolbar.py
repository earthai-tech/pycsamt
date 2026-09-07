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
    _register_geology_legend_quick(app)
    _register_viewport(app)
    _register_spin(app)


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
                v === 'map3d' ? show : hide,
                v === 'bh'    ? show : hide
            ];
        }
        """,
        Output(IDs.TOOLBAR_2D, "style"),
        Output(IDs.TOOLBAR_3D, "style"),
        Output(IDs.TOOLBAR_BH, "style"),
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


def _register_viewport(app) -> None:
    """Persist the user's live camera / pan-zoom so a control change never
    moves the scene.

    The main render callback rebuilds the whole figure on every control
    change. ``uirevision`` is supposed to keep the user's view across that,
    but it leaks here: the 2-D ``map`` (MapLibre) subplot re-applies its
    ``center``/``zoom`` on each rebuild, and ``dcc.Loading`` can blank the
    3-D camera mid-update. So we capture every ``relayoutData`` into a
    per-view store (client-side, no server round-trip) and
    ``_render._apply_viewport`` replays it onto the fresh figure. Fit /
    Reset view clear the store, which is what lets those buttons — and
    only those buttons — snap back to the data.
    """
    app.clientside_callback(
        """
        function(relayout, view, store) {
            if (!relayout) return window.dash_clientside.no_update;
            var s = Object.assign({}, store || {});
            var v = view || 'map';
            var cur = Object.assign({}, s[v] || {});
            var touched = false;
            if (v === 'map3d') {
                var cam = relayout['scene.camera'];
                if (!cam && relayout['scene.camera.eye']) {
                    cam = Object.assign({}, cur.camera || {});
                    cam.eye = relayout['scene.camera.eye'];
                    if (relayout['scene.camera.center'])
                        cam.center = relayout['scene.camera.center'];
                    if (relayout['scene.camera.up'])
                        cam.up = relayout['scene.camera.up'];
                }
                if (cam) { cur.camera = cam; touched = true; }
            } else {
                ['map', 'mapbox'].forEach(function(p) {
                    if (relayout[p + '.center'] !== undefined) {
                        cur.center = relayout[p + '.center']; touched = true;
                    }
                    if (relayout[p + '.center.lon'] !== undefined) {
                        cur.center = {
                            lon: relayout[p + '.center.lon'],
                            lat: relayout[p + '.center.lat']
                        };
                        touched = true;
                    }
                    if (relayout[p + '.zoom'] !== undefined) {
                        cur.zoom = relayout[p + '.zoom']; touched = true;
                    }
                    if (relayout[p + '.bearing'] !== undefined) {
                        cur.bearing = relayout[p + '.bearing']; touched = true;
                    }
                    if (relayout[p + '.pitch'] !== undefined) {
                        cur.pitch = relayout[p + '.pitch']; touched = true;
                    }
                });
                if (relayout['xaxis.range[0]'] !== undefined) {
                    cur.xrange = [
                        relayout['xaxis.range[0]'], relayout['xaxis.range[1]']
                    ];
                    touched = true;
                }
                if (relayout['yaxis.range[0]'] !== undefined) {
                    cur.yrange = [
                        relayout['yaxis.range[0]'], relayout['yaxis.range[1]']
                    ];
                    touched = true;
                }
            }
            if (!touched) return window.dash_clientside.no_update;
            s[v] = cur;
            return s;
        }
        """,
        Output(IDs.STORE_VIEWPORT, "data"),
        Input(IDs.CANVAS_GRAPH, "relayoutData"),
        State(IDs.STORE_VIEW, "data"),
        State(IDs.STORE_VIEWPORT, "data"),
        prevent_initial_call=True,
    )
    # Fit / Reset view (they bump STORE_FIT) forget the saved viewport so
    # the render falls back to the data-fitted default.
    app.clientside_callback(
        "function(_fit){ return {}; }",
        Output(IDs.STORE_VIEWPORT, "data", allow_duplicate=True),
        Input(IDs.STORE_FIT, "data"),
        prevent_initial_call=True,
    )


def _register_spin(app) -> None:
    """3-D turntable: a toolbar toggle that orbits the camera eye about
    the vertical axis, one small step per :class:`dcc.Interval` tick.

    Pure client-side (``Plotly.relayout`` on the live graph div) so it
    never hits the server or fights the render callback; each frame's
    ``relayoutData`` flows through ``_register_viewport`` like a manual
    drag, so stopping the spin just leaves the view where it landed.
    """
    app.clientside_callback(
        "function(n, cur){"
        " if (!n) return window.dash_clientside.no_update;"
        " return !cur; }",
        Output(IDs.STORE_SPIN, "data"),
        Input(IDs.TB3D_SPIN, "n_clicks"),
        State(IDs.STORE_SPIN, "data"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        """
        function(spin, view) {
            var on = !!spin && (view === 'map3d');
            return [on ? 'mv-tb-btn active' : 'mv-tb-btn', !on];
        }
        """,
        Output(IDs.TB3D_SPIN, "className"),
        Output(IDs.SPIN_INTERVAL, "disabled"),
        Input(IDs.STORE_SPIN, "data"),
        Input(IDs.STORE_VIEW, "data"),
        prevent_initial_call=False,
    )
    app.clientside_callback(
        """
        function(_n) {
            try {
                var gd = document.getElementById('mv-canvas-graph');
                if (gd && !gd._fullLayout) {
                    gd = gd.querySelector('.js-plotly-plot');
                }
                if (!gd || !gd._fullLayout || !gd._fullLayout.scene
                    || !window.Plotly) {
                    return window.dash_clientside.no_update;
                }
                var e = gd._fullLayout.scene.camera.eye;
                var a = 0.03;
                var nx = e.x * Math.cos(a) - e.y * Math.sin(a);
                var ny = e.x * Math.sin(a) + e.y * Math.cos(a);
                window.Plotly.relayout(
                    gd, {'scene.camera.eye': {x: nx, y: ny, z: e.z}}
                );
            } catch (err) { /* graph not ready yet */ }
            return window.dash_clientside.no_update;
        }
        """,
        Output(IDs.SPIN_TICK, "data"),
        Input(IDs.SPIN_INTERVAL, "n_intervals"),
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


def _register_geology_legend_quick(app) -> None:
    """3-D toolbar Legend button flips GEO_LEGEND_VISIBLE (same switch
    the Inspector's Interpretation section uses) -- hides the on-canvas
    geology legend to give the plot back its full width without
    discarding the applied legend itself."""
    app.clientside_callback(
        "function(n, cur){ if(!n) return window.dash_clientside.no_update;"
        " return !cur; }",
        Output(IDs.GEO_LEGEND_VISIBLE, "value", allow_duplicate=True),
        Input(IDs.TB3D_LEGEND, "n_clicks"),
        State(IDs.GEO_LEGEND_VISIBLE, "value"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        "function(v){ return v ? 'mv-tb-btn active' : 'mv-tb-btn'; }",
        Output(IDs.TB3D_LEGEND, "className"),
        Input(IDs.GEO_LEGEND_VISIBLE, "value"),
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
