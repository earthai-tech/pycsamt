# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Inspector controls → STORE_CONTROLS, and frequency-slider population."""

from __future__ import annotations

from dash import Input, Output, State, ctx

from .._ids import IDs

_DEPTH_PRESETS = {
    IDs.BTN_DEPTH_FULL: (None, None),
    IDs.BTN_DEPTH_500: (0, 500),
    IDs.BTN_DEPTH_1K: (0, 1000),
    IDs.BTN_DEPTH_2K: (0, 2000),
    # 3-D toolbar duplicates of the same presets (kept in sync — one dict,
    # one callback — see _register_depth_presets).
    IDs.TB3D_DEPTH_FULL: (None, None),
    IDs.TB3D_DEPTH_500: (0, 500),
    IDs.TB3D_DEPTH_1K: (0, 1000),
    IDs.TB3D_DEPTH_2K: (0, 2000),
}

_RHO_PRESETS = {
    IDs.BTN_RHO_ALL: (1, 100_000),
    IDs.BTN_RHO_COND: (1, 100),
    IDs.BTN_RHO_MID: (100, 1_000),
    IDs.BTN_RHO_RES: (1_000, 100_000),
}


def register_controls(app) -> None:
    _register_freq_slider(app)
    _register_gather(app)
    _register_group_visibility(app)
    _register_depth_presets(app)
    _register_rho_presets(app)


def _register_depth_presets(app) -> None:
    @app.callback(
        Output(IDs.CTL_DEPTH_LO, "value"),
        Output(IDs.CTL_DEPTH_HI, "value"),
        Input(IDs.BTN_DEPTH_FULL, "n_clicks"),
        Input(IDs.BTN_DEPTH_500, "n_clicks"),
        Input(IDs.BTN_DEPTH_1K, "n_clicks"),
        Input(IDs.BTN_DEPTH_2K, "n_clicks"),
        Input(IDs.TB3D_DEPTH_FULL, "n_clicks"),
        Input(IDs.TB3D_DEPTH_500, "n_clicks"),
        Input(IDs.TB3D_DEPTH_1K, "n_clicks"),
        Input(IDs.TB3D_DEPTH_2K, "n_clicks"),
        prevent_initial_call=True,
    )
    def apply_preset(*_clicks):
        lo, hi = _DEPTH_PRESETS.get(ctx.triggered_id, (None, None))
        return lo, hi


def _register_rho_presets(app) -> None:
    @app.callback(
        Output(IDs.CTL_RHO_LO, "value"),
        Output(IDs.CTL_RHO_HI, "value"),
        Input(IDs.BTN_RHO_ALL, "n_clicks"),
        Input(IDs.BTN_RHO_COND, "n_clicks"),
        Input(IDs.BTN_RHO_MID, "n_clicks"),
        Input(IDs.BTN_RHO_RES, "n_clicks"),
        prevent_initial_call=True,
    )
    def apply_rho_preset(*_clicks):
        lo, hi = _RHO_PRESETS.get(ctx.triggered_id, (None, None))
        return lo, hi


def _register_group_visibility(app) -> None:
    """Show the relevant inspector control groups for the active view."""
    app.clientside_callback(
        """
        function(view) {
            var v = view || 'map';
            var show = {display:'block'}, hide = {display:'none'};
            return [
                v === 'map'   ? show : hide,   // Map controls
                v === 'map3d' ? show : hide    // 3-D controls
            ];
        }
        """,
        Output(IDs.GRP_MAP, "style"),
        Output(IDs.GRP_3D, "style"),
        Input(IDs.STORE_VIEW, "data"),
        prevent_initial_call=False,
    )


def _fmt_freq(freq) -> str:
    if freq is None:
        return "—"
    if freq >= 1000:
        return f"{freq / 1000:.3g} kHz"
    if freq >= 1:
        return f"{freq:.4g} Hz"
    return f"{freq:.4g} Hz ({1.0 / freq:.4g} s)"


def _register_freq_slider(app) -> None:
    @app.callback(
        Output(IDs.CTL_FREQUENCY, "min"),
        Output(IDs.CTL_FREQUENCY, "max"),
        Output(IDs.CTL_FREQUENCY, "marks"),
        Output(IDs.CTL_FREQUENCY, "value"),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def populate(store):
        freqs = (store or {}).get("frequencies", [])
        if not freqs:
            return 0, 0, None, 0
        n = len(freqs)
        # sparse marks: ends + middle
        idxs = sorted({0, n // 2, n - 1})
        marks = {
            i: {"label": _fmt_freq(freqs[i]), "style": {"fontSize": "9px"}}
            for i in idxs
        }
        return 0, n - 1, marks, 0


def _register_gather(app) -> None:
    @app.callback(
        Output(IDs.STORE_CONTROLS, "data"),
        Output(IDs.CTL_FREQ_LABEL, "children"),
        Input(IDs.CTL_OVERLAY, "value"),
        Input(IDs.CTL_QUANTITY, "value"),
        Input(IDs.CTL_COMPONENT, "value"),
        Input(IDs.CTL_FREQUENCY, "value"),
        Input(IDs.CTL_CMAP, "value"),
        Input(IDs.CTL_LOG, "value"),
        Input(IDs.CTL_MODE3D, "value"),
        Input(IDs.CTL_LABELS, "value"),
        Input(IDs.CTL_OPACITY, "value"),
        Input(IDs.CTL_AZIMUTH, "value"),
        Input(IDs.CTL_SPACING, "value"),
        Input(IDs.CTL_DEPTH_LO, "value"),
        Input(IDs.CTL_DEPTH_HI, "value"),
        Input(IDs.CTL_NSLICES, "value"),
        Input(IDs.CTL_SURFACES, "value"),
        Input(IDs.CTL_CONTOURS, "value"),
        Input(IDs.CTL_SCALE, "value"),
        Input(IDs.CTL_VMIN, "value"),
        Input(IDs.CTL_VMAX, "value"),
        Input(IDs.CTL_RHO_LO, "value"),
        Input(IDs.CTL_RHO_HI, "value"),
        Input(IDs.CTL_TOPO, "value"),
        Input(IDs.CTL_TERRAIN, "value"),
        Input(IDs.CTL_BASEMAP, "value"),
        Input(IDs.CTL_CONTOUR_LEVELS, "value"),
        Input(IDs.CTL_CONTOUR_MODE, "value"),
        Input(IDs.CTL_SHOW_STA, "value"),
        Input(IDs.CTL_STA_LABELS, "value"),
        Input(IDs.CTL_STA_SYMBOL, "value"),
        Input(IDs.CTL_STA_SIZE, "value"),
        Input(IDs.CTL_STA_COLOR, "value"),
        Input(IDs.CTL_CONTOUR_ENABLE, "value"),
        Input(IDs.CTL_MARKER_SIZE, "value"),
        Input(IDs.CTL_OPACITY_MAP, "value"),
        Input(IDs.CTL_PROFILES, "value"),
        Input(IDs.CTL_CRS_MODE, "value"),
        Input(IDs.CTL_UTM_ZONE, "value"),
        Input(IDs.CTL_UTM_HEM, "value"),
        Input(IDs.CTL_EPSG, "value"),
        Input(IDs.CTL_CONTOUR_INTERP, "value"),
        Input(IDs.CTL_CONTOUR_SMOOTH, "value"),
        Input(IDs.CTL_CONTOUR_RES, "value"),
        Input(IDs.CTL_ASPECT, "value"),
        Input(IDs.CTL_X_UNIT, "value"),
        Input(IDs.CTL_DEPTH_UNIT, "value"),
        Input(IDs.CTL_SMOOTH, "value"),
        Input(IDs.CTL_SECTION_RES, "value"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def gather(
        overlay,
        quantity,
        component,
        freq_idx,
        cmap,
        log,
        mode3d,
        labels,
        opacity,
        azimuth,
        spacing,
        depth_lo,
        depth_hi,
        n_slices,
        surfaces,
        contours,
        scale,
        vmin,
        vmax,
        rho_lo,
        rho_hi,
        topo,
        terrain,
        basemap,
        contour_levels,
        contour_mode,
        show_sta,
        sta_labels,
        sta_symbol,
        sta_size,
        sta_color,
        contour_enable,
        marker_size,
        map_opacity,
        profiles,
        crs_mode,
        utm_zone,
        utm_hem,
        epsg,
        contour_interp,
        contour_smooth,
        contour_res,
        aspect,
        x_unit,
        depth_unit,
        smooth_sections,
        section_res,
        store,
    ):
        freqs = (store or {}).get("frequencies", [])
        freq = None
        if freqs:
            i = int(freq_idx or 0)
            i = max(0, min(i, len(freqs) - 1))
            freq = float(freqs[i])
        controls = {
            "overlay": overlay or "index",
            "quantity": quantity or "rho",
            "component": component or "xy",
            "frequency": freq,
            "cmap": cmap or "plasma",
            "log": bool(log),
            "mode3d": mode3d or "fence",
            "labels": bool(labels),
            "opacity": opacity if opacity is not None else 0.85,
            "azimuth": azimuth if azimuth is not None else 0.0,
            "line_spacing": spacing if spacing is not None else 1.0,
            "depth_lo": depth_lo,
            "depth_hi": depth_hi,
            "n_slices": n_slices if n_slices is not None else 8,
            "surface_count": surfaces if surfaces is not None else 12,
            "contours": bool(contours),
            "scale": scale or "log",
            "vmin": vmin,
            "vmax": vmax,
            "rho_lo": rho_lo,
            "rho_hi": rho_hi,
            "topography": bool(topo),
            "terrain": bool(terrain),
            "basemap": basemap or "esri-satellite",
            "contour_levels": contour_levels if contour_levels else 12,
            "contour_mode": contour_mode or "filled+lines",
            "show_stations": bool(show_sta),
            "station_labels": bool(sta_labels),
            "station_symbol": sta_symbol or "diamond",
            "station_size": sta_size if sta_size else 4,
            "station_color": sta_color or "#1f2937",
            "contour_enable": bool(contour_enable),
            "marker_size": marker_size if marker_size else 10,
            "map_opacity": map_opacity if map_opacity else 92,
            "profiles": bool(profiles),
            "crs_mode": crs_mode or "geo",
            "utm_zone": utm_zone,
            "utm_hem": utm_hem or "N",
            "epsg": epsg,
            "contour_interp": contour_interp or "cubic",
            "contour_smooth": contour_smooth if contour_smooth is not None else 1.0,
            "contour_res": int(contour_res) if contour_res else 150,
            "aspect": aspect or "data",
            "x_unit": x_unit or "m",
            "depth_unit": depth_unit or "m",
            "smooth_sections": bool(smooth_sections),
            "section_res": int(section_res) if section_res else 100,
        }
        return controls, _fmt_freq(freq)
