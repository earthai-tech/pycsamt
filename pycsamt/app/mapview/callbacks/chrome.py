# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Chrome callbacks: seed adoption, theme, help, sidebar, dock."""

from __future__ import annotations

from dash import Input, Output, State, no_update

from .._ids import IDs
from .._render import store_from_view
from ..cache import set_view, take_seed, take_seed_state


def register_chrome(app) -> None:
    _register_seed(app)
    _register_seed_state(app)
    _register_theme(app)
    _register_help(app)
    _register_sidebar(app)
    _register_inspector(app)
    _register_dock(app)


# Plotly's responsive mode only re-lays out on a window resize event, so we
# fire one after toggling a panel to make the canvas reclaim/yield space.
_RESIZE = (
    "setTimeout(function(){window.dispatchEvent(new Event('resize'));}, 80);"
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


# widget -> the seed-state control feeding it (see pcsf_scene.mapview_state)
_SEED_WIDGETS = (
    (IDs.CTL_MODE3D, "mode3d"),
    (IDs.CTL_CMAP, "cmap"),
    (IDs.CTL_DEPTH_LO, "depth_lo"),
    (IDs.CTL_DEPTH_HI, "depth_hi"),
    (IDs.CTL_TOPO, "topography"),
    (IDs.CTL_TERRAIN, "terrain"),
    (IDs.CTL_OPACITY, "opacity"),
    (IDs.CTL_SHOW_STA, "show_stations"),
    (IDs.CTL_STA_LABELS, "station_labels"),
    (IDs.GEO_LEGEND_VISIBLE, "geology_legend"),
    (IDs.TB3D_GEO_FILL, "geology_fill"),
    (IDs.CTL_ASPECT, "aspect"),
    (IDs.CTL_VE, "vertical_exaggeration"),
    (IDs.CTL_STA_SYMBOL, "station_symbol"),
    (IDs.CTL_STA_SIZE, "station_size"),
    (IDs.CTL_STA_COLOR, "station_color"),
    (IDs.CTL_STA_MAX, "station_max"),
    (IDs.CTL_STA_LABEL_ANGLE, "station_label_angle"),
    (IDs.CTL_STA_LABEL_DENSITY, "station_label_density"),
    (IDs.CTL_STA_LABEL_NAMES, "station_label_names"),
    (IDs.CTL_LABELS, "labels"),
    (IDs.CTL_SCALE, "scale"),
    (IDs.CTL_VMIN, "vmin"),
    (IDs.CTL_VMAX, "vmax"),
    (IDs.CTL_CRANGE_PLO, "crange_plo"),
    (IDs.CTL_CRANGE_PHI, "crange_phi"),
    (IDs.CTL_RHO_LO, "rho_lo"),
    (IDs.CTL_RHO_HI, "rho_hi"),
    (IDs.CTL_RHO_CUTOFF, "rho_cutoff"),
    (IDs.CTL_CONTOURS, "contours"),
    (IDs.CTL_NSLICES, "n_slices"),
    (IDs.CTL_SURFACES, "surface_count"),
    (IDs.CTL_SPACING, "line_spacing"),
    (IDs.CTL_AZIMUTH, "azimuth"),
    (IDs.CTL_X_UNIT, "x_unit"),
    (IDs.CTL_DEPTH_UNIT, "depth_unit"),
    (IDs.CTL_SMOOTH, "smooth_sections"),
    (IDs.CTL_SECTION_RES, "section_res"),
    (IDs.CTL_VOL_SMOOTH, "volume_smoothing"),
)


def _register_seed_state(app) -> None:
    """Open on the scene the launcher had (desktop "Open in Map View").

    Runs once, when the seeded survey lands in STORE_DATA: sets the 3-D
    widgets (so Map View's own control gathering keeps them), the
    overlays, spin, theme and camera, and switches to the 3-D view by
    "clicking" its rail button.
    """

    @app.callback(
        *[Output(wid, "value", allow_duplicate=True)
          for wid, _k in _SEED_WIDGETS],
        Output(IDs.GEO_STORE, "data", allow_duplicate=True),
        Output(IDs.PCBH_STORE, "data", allow_duplicate=True),
        Output(IDs.STRUCT_STORE, "data", allow_duplicate=True),
        Output(IDs.STORE_SPIN, "data", allow_duplicate=True),
        Output(IDs.STORE_THEME, "data", allow_duplicate=True),
        Output(IDs.STORE_VIEWPORT, "data", allow_duplicate=True),
        Output(IDs.RAIL_3D, "n_clicks", allow_duplicate=True),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def apply_seed_state(store):
        if not store:
            return seed_state_outputs(None)
        return seed_state_outputs(take_seed_state())


def seed_state_outputs(state: dict | None) -> tuple:
    """Callback outputs for a seed *state* (``no_update`` where unset):
    the 3-D widgets, GEO/PCBH/STRUCT stores, spin, theme, viewport and
    the 3-D rail click."""
    n = len(_SEED_WIDGETS) + 7
    if not state:
        return (no_update,) * n
    c = state.get("controls") or {}
    widgets = [c[k] if k in c else no_update for _w, k in _SEED_WIDGETS]
    camera = state.get("camera")
    return (
        *widgets,
        state.get("geo") or no_update,
        state.get("pcbh") or no_update,
        state.get("struct") or no_update,
        bool(state.get("spin", False)),
        state.get("theme") or no_update,
        {"map3d": {"camera": camera}} if camera else no_update,
        1 if state.get("view") == "map3d" else no_update,
    )


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
            """
        + _RESIZE
        + """
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
            """
        + _RESIZE
        + """
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
        function(toggleClicks, closeClicks, currentStyle) {
            var trig = window.dash_clientside.callback_context.triggered;
            if (!trig.length) return window.dash_clientside.no_update;
            var id = trig[0].prop_id.split('.')[0];
            var isOpen = currentStyle && currentStyle.display === 'block';
            var open = id === 'mv-dock-close' ? false : !isOpen;
            """
        + _RESIZE
        + """
            return [
                open ? {display:'block'} : {display:'none'},
                open ? 'bi bi-chevron-down ms-2' : 'bi bi-chevron-up ms-2'
            ];
        }
        """,
        Output(IDs.DOCK_BODY, "style"),
        Output("mv-dock-chevron", "className"),
        Input(IDs.DOCK_TOGGLE, "n_clicks"),
        Input(IDs.DOCK_CLOSE, "n_clicks"),
        State(IDs.DOCK_BODY, "style"),
        prevent_initial_call=True,
    )
