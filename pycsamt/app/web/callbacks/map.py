# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""callbacks/map.py — station map figure, map-type controls, selection sync."""

from __future__ import annotations

import numpy as np
import pandas as pd
from dash import ALL, ctx, clientside_callback, Input, Output, State, no_update

from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import build_station_map, build_contour_overlay


def _station_dataframe(store_data) -> pd.DataFrame:
    """Return normalized station records ready for map rendering."""
    records = (store_data or {}).get("station_records", [])
    df = pd.DataFrame(records)
    if df.empty:
        return df

    aliases = {
        "ID": ("ID", "id", "station", "Station", "name", "Name"),
        "Latitude": ("Latitude", "latitude", "lat", "Lat"),
        "Longitude": ("Longitude", "longitude", "lon", "Lon", "long", "Long"),
        "Line": ("Line", "line", "profile", "Profile"),
        "Elevation": ("Elevation", "elevation", "elev", "Elev"),
    }
    for canonical, candidates in aliases.items():
        if canonical in df.columns:
            continue
        match = next((c for c in candidates if c in df.columns), None)
        if match:
            df[canonical] = df[match]

    if "ID" in df.columns:
        df["ID"] = df["ID"].astype(str)
    if "Line" in df.columns:
        df["Line"] = df["Line"].map(
            lambda v: str(v).strip() if pd.notna(v) and str(v).strip() else "—"
        )
    for col in ("Latitude", "Longitude", "Elevation"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _usable_line_filter(line_filter, df: pd.DataFrame) -> list[str] | None:
    """Return only line-filter values that actually match the station records."""
    if not line_filter or "Line" not in df.columns:
        return None
    requested = [str(v).strip() for v in line_filter if str(v).strip()]
    if not requested:
        return None
    available = set(df["Line"].dropna().astype(str))
    keep = [line for line in requested if line in available]
    return keep or None


def _layout_has_id(layout, target_id: str) -> bool:
    """Return True when *target_id* exists in the current Dash layout tree."""
    if layout is None:
        return False
    if isinstance(layout, (list, tuple)):
        return any(_layout_has_id(child, target_id) for child in layout)
    if str(getattr(layout, "id", "")) == target_id:
        return True
    return _layout_has_id(getattr(layout, "children", None), target_id)


# ── EPSG info helper ──────────────────────────────────────────────────────────

def _resolve_crs_info(mode: str, zone: int = 50, hem: str = "N",
                      epsg: str = "4326") -> str:
    """Return a short human-readable CRS description."""
    if mode == "geo":
        return "EPSG:4326  Geographic lat/lon (WGS 84)"
    if mode == "utm":
        code = 32600 + int(zone or 50) if (hem or "N") == "N" else 32700 + int(zone or 50)
        return f"EPSG:{code}  UTM Zone {zone}{hem} (WGS 84)"
    # custom
    try:
        from pyproj import CRS
        c = CRS.from_user_input(f"EPSG:{epsg}")
        return f"EPSG:{epsg}  {c.name}"
    except Exception:
        return f"EPSG:{epsg}  (cannot validate — pyproj needed)"


# ── Resistivity / depth helpers ───────────────────────────────────────────────

def _extract_rho_at_freq(sites, freq_hz: float, comp: str) -> dict[str, float]:
    """
    Return {station_id: rho_at_freq} for the closest frequency to *freq_hz*.
    comp must be one of "xy", "yx", "det".
    """
    _IDX = {"xy": (0, 1), "yx": (1, 0), "xx": (0, 0), "yy": (1, 1)}
    result: dict[str, float] = {}
    for edi in sites.as_list():
        try:
            Z    = edi.Z
            freq = np.asarray(Z.freq, dtype=float)
            if freq.size == 0:
                continue
            i    = int(np.argmin(np.abs(freq - freq_hz)))
            if comp == "det":
                rho = float(Z.resistivity[i, 0, 1] * Z.resistivity[i, 1, 0]) ** 0.5
            else:
                r, c = _IDX.get(comp, (0, 1))
                rho  = float(Z.resistivity[i, r, c])
            if np.isfinite(rho) and rho > 0:
                result[edi.station] = rho
        except Exception:
            pass
    return result


def _extract_depth_at_freq(sites, freq_hz: float,
                            comp: str, target_depth_m: float = 0.0
                            ) -> dict[str, float]:
    """
    Return {station_id: skin_depth_m} where δ = 503 * sqrt(ρ / f).
    If target_depth_m > 0, return period T at which δ ≈ target_depth instead.
    """
    rho_map = _extract_rho_at_freq(sites, freq_hz, comp)
    if target_depth_m > 0:
        result = {}
        for sid, rho in rho_map.items():
            # T = (depth / 503)^2 / rho
            T = (target_depth_m / 503.0) ** 2 / rho if rho > 0 else float("nan")
            result[sid] = T if np.isfinite(T) else float("nan")
        return result
    return {sid: 503.0 * np.sqrt(rho / max(freq_hz, 1e-9))
            for sid, rho in rho_map.items()}


# ──────────────────────────────────────────────────────────────────────────────


def register_map(app) -> None:
    has_home_map = (
        _layout_has_id(app.layout, IDs.HOME_MAP_GRAPH)
        and _layout_has_id(app.layout, IDs.HOME_MAP_OVERLAY)
    )
    if has_home_map:
        _register_home_map(app)
    _register_map_update(app)
    _register_freq_options(app)
    _register_depth_info(app)
    _register_inv_controls(app)
    _register_crs_info(app)
    _register_selection_sync(app, include_home_map=has_home_map)
    _register_status_station(app)
    _register_map_type_visibility(app)
    _register_crs_visibility(app)
    _register_contour_visibility(app)
    _register_toolbar(app)


# ── Home-page mini-map (separate ID to avoid duplicate-ID conflicts) ──────────

def _register_home_map(app) -> None:
    @app.callback(
        Output(IDs.HOME_MAP_GRAPH,  "figure"),
        Input(IDs.STORE_DATA,       "data"),
        Input(IDs.HOME_MAP_OVERLAY, "value"),
        Input(IDs.STORE_SELECTION,  "data"),
        Input(IDs.STORE_THEME,      "data"),
    )
    def update_home_map(store_data, overlay, selection, theme):
        dark = (theme or "dark") == "dark"
        if not store_data:
            return build_station_map(pd.DataFrame(), dark=dark)
        df = _station_dataframe(store_data)
        if df.empty:
            return build_station_map(df, dark=dark)
        selected_id = (selection or {}).get("station_id")
        return build_station_map(df, selected_id=selected_id,
                                 overlay=overlay or "Index", dark=dark)


# ── Map View full page ────────────────────────────────────────────────────────

def _register_map_update(app) -> None:
    @app.callback(
        Output(IDs.MAP_GRAPH,                 "figure"),
        # Page-visibility trigger — fires when user navigates to Map View
        Input(IDs.NAV_SECTION,                "data"),
        Input(IDs.STORE_DATA,                 "data"),
        Input(IDs.MAP_OVERLAY,                "value"),
        Input(IDs.MAP_PAGE_MAP_TYPE,          "value"),
        Input(IDs.MAP_PAGE_FREQ_SEL,          "value"),
        Input(IDs.MAP_PAGE_COMP,              "value"),
        Input(IDs.MAP_PAGE_DEPTH_TARGET,      "value"),
        Input(IDs.STORE_SELECTION,            "data"),
        Input(IDs.STORE_THEME,                "data"),
        Input(IDs.MAP_PAGE_LINE_FILTER,       "value"),
        Input(IDs.MAP_PAGE_MARKER_SIZE,       "value"),
        Input(IDs.MAP_PAGE_BASEMAP,           "value"),
        Input(IDs.MAP_PAGE_CMAP,              "value"),
        Input(IDs.MAP_PAGE_OPACITY,           "value"),
        Input(IDs.MAP_PAGE_SHOW_LABELS,       "value"),
        # CRS — triggers reprojection if user switches coordinate system
        Input(IDs.MAP_PAGE_CRS_MODE,          "value"),
        Input(IDs.MAP_PAGE_UTM_ZONE,          "value"),
        Input(IDs.MAP_PAGE_UTM_HEM,           "value"),
        Input(IDs.MAP_PAGE_EPSG,              "value"),
        # Contour inputs
        Input(IDs.MAP_PAGE_CONTOUR_ENABLE,    "value"),
        Input(IDs.MAP_PAGE_CONTOUR_MODE,      "value"),
        Input(IDs.MAP_PAGE_CONTOUR_CMAP,      "value"),
        Input(IDs.MAP_PAGE_CONTOUR_LEVELS,    "value"),
        Input(IDs.MAP_PAGE_CONTOUR_OPACITY,   "value"),
        Input(IDs.MAP_PAGE_CONTOUR_INTERP,    "value"),
        Input(IDs.MAP_PAGE_CONTOUR_EXTRA,     "value"),
        Input(IDs.MAP_PAGE_CONTOUR_SMOOTH,    "value"),
        Input(IDs.MAP_PAGE_CONTOUR_LWIDTH,    "value"),
        Input(IDs.MAP_PAGE_CONTOUR_OPTS,      "value"),
        Input(IDs.MAP_PAGE_CONTOUR_RES,       "value"),
        # Inversion overlay
        Input(IDs.MAP_PAGE_INV_DEPTH_M,       "value"),
        Input(IDs.MAP_PAGE_INV_LINE_SEL,      "value"),
        # Bearing / rotation
        Input(IDs.MAP_PAGE_BEARING,           "value"),
        # Toolbar buttons that trigger map rebuild
        Input(IDs.MAP_PAGE_REFRESH,           "n_clicks"),
        Input("map-tb-fit",                   "n_clicks"),
        State(IDs.SESSION_ID,                 "data"),
        State(IDs.STORE_ACTIVE_LINES,         "data"),
    )
    def update_map(nav_section,
                   store_data, overlay, map_type, freq_sel, comp, depth_target,
                   selection, theme, line_filter, marker_size, basemap_style,
                   cmap, opacity_pct, show_opts,
                   crs_mode, utm_zone, utm_hem, epsg,
                   contour_enable, contour_mode, contour_cmap,
                   contour_levels, contour_opacity_pct, contour_interp,
                   contour_extra, contour_smooth, contour_lwidth,
                   contour_opts, contour_res,
                   inv_depth_m, inv_line_sel,
                   bearing,
                   _refresh, _fit,
                   session_id, active_lines_store):

        # Skip re-render when user navigates away from the map page.
        # Other inputs (data change, param change, refresh click) always proceed.
        if ctx.triggered_id == IDs.NAV_SECTION and nav_section != "map":
            return no_update

        dark          = (theme or "dark") == "dark"
        selected_id   = (selection or {}).get("station_id")
        show_labels   = "labels"   in (show_opts or [])
        show_profiles = "profiles" in (show_opts or [])
        opacity       = float(opacity_pct or 90) / 100.0
        mtype         = map_type or "station"

        if not store_data:
            return build_station_map(pd.DataFrame(), dark=dark)

        df = _station_dataframe(store_data)
        if df.empty:
            return build_station_map(df, dark=dark)

        # ── CRS projection ───────────────────────────────────────────────────
        # If the user specified a non-geographic CRS, reproject X/Y → lat/lon.
        crs_mode = crs_mode or "geo"
        if crs_mode != "geo" and "Latitude" in df.columns and "Longitude" in df.columns:
            xs = pd.to_numeric(df["Longitude"], errors="coerce").values
            ys = pd.to_numeric(df["Latitude"],  errors="coerce").values
            valid_crs = np.isfinite(xs) & np.isfinite(ys)
            if valid_crs.any():
                try:
                    from pyproj import Transformer
                    if crs_mode == "utm":
                        zone = int(utm_zone or 50)
                        hem  = (utm_hem or "N").upper()
                        epsg_code = (32600 + zone) if hem == "N" else (32700 + zone)
                    else:
                        epsg_code = int(epsg or 4326)
                    transformer = Transformer.from_crs(
                        f"EPSG:{epsg_code}", "EPSG:4326", always_xy=True
                    )
                    lon_proj, lat_proj = transformer.transform(xs[valid_crs], ys[valid_crs])
                    df = df.copy()
                    df.loc[valid_crs, "Longitude"] = lon_proj
                    df.loc[valid_crs, "Latitude"]  = lat_proj
                except Exception:
                    pass  # keep original coordinates on pyproj error

        # ── uirevision — "fit" click resets zoom; all other triggers preserve it
        _fit_rev = f"fit-{_fit}" if ctx.triggered_id == "map-tb-fit" else "map"

        # ── Map type processing ──────────────────────────────────────────────
        sites     = None
        value_map = None   # {station_id: float} used later by contour

        if mtype == "elevation":
            overlay = "Elevation"

        elif mtype in ("resistivity", "depth") and freq_sel and session_id:
            sites = cache_get(session_id)
            if sites:
                freq_hz = float(freq_sel)
                if mtype == "resistivity":
                    rho_map = _extract_rho_at_freq(sites, freq_hz, comp or "xy")
                    if rho_map:
                        # Store as log10 for sensible colour stretch
                        log_rho = {k: float(np.log10(max(v, 1e-9))) for k, v in rho_map.items()}
                        df["_rho_map"] = df["ID"].map(log_rho)
                        overlay   = "_rho_map"
                        value_map = log_rho
                else:
                    dep_map = _extract_depth_at_freq(
                        sites, freq_hz, comp or "xy",
                        float(depth_target or 0),
                    )
                    if dep_map:
                        df["_depth_map"] = df["ID"].map(dep_map)
                        overlay   = "_depth_map"
                        value_map = dep_map

        elif mtype in ("inv_profile", "inv_depth_slice") and session_id:
            from pycsamt.app.web.cache import cache_get_inversion_result
            result = cache_get_inversion_result(session_id)
            if result is not None:
                try:
                    model  = result.to_resistivity_model()
                    rho_2d = np.asarray(model.rho_2d, float)        # (n_z, n_x) log10 or Ω·m
                    z_arr  = np.asarray(model.z_centers, float)      # depth centres (m)
                    # Normalise to Ω·m if stored in log10 scale (max < 12)
                    finite = rho_2d[np.isfinite(rho_2d)]
                    if finite.size and float(np.nanmax(finite)) <= 12.0:
                        rho_2d = 10.0 ** rho_2d

                    sta_names = list(getattr(model, "station_names", None) or
                                     [f"S{i:03d}" for i in range(rho_2d.shape[1])])

                    if mtype == "inv_profile":
                        # Colour by mean log10(ρ) across all depth layers
                        mean_rho = np.nanmean(rho_2d, axis=0)       # (n_x,)
                        rho_col  = {n: float(np.log10(max(r, 1e-9)))
                                    for n, r in zip(sta_names, mean_rho)}
                    else:
                        # Nearest depth layer to the requested depth
                        target_d = float(inv_depth_m or 500)
                        iz       = int(np.argmin(np.abs(z_arr - target_d)))
                        rho_col  = {n: float(np.log10(max(r, 1e-9)))
                                    for n, r in zip(sta_names, rho_2d[iz, :])}

                    if rho_col:
                        df["_inv_rho"] = df["ID"].map(rho_col)
                        overlay        = "_inv_rho"
                        value_map      = rho_col
                        # Depth-slice: force contour overlay on (user can still disable)
                        if mtype == "inv_depth_slice" and "on" not in (contour_enable or []):
                            contour_enable = ["on"]
                except Exception:
                    pass  # fall back to default station view if inversion result is malformed

        # Derive muted lines: lines in "all" that are not in "active"
        _als = active_lines_store or {}
        _all_lines = _als.get("all", [])
        _active     = set(_als.get("active", _all_lines))
        muted_lines = [ln for ln in _all_lines if ln not in _active] or None

        fig = build_station_map(
            df,
            selected_id=selected_id,
            overlay=overlay or "Index",
            dark=dark,
            line_filter=_usable_line_filter(line_filter, df),
            marker_size=int(marker_size) if marker_size else 10,
            basemap_style=basemap_style or None,
            cmap=cmap or "plasma",
            opacity=opacity,
            show_labels=show_labels,
            show_profiles=show_profiles,
            muted_lines=muted_lines,
            uirevision=_fit_rev,
            bearing=float(bearing or 0),
        )

        # ── Contour overlay ──────────────────────────────────────────────────
        if "on" in (contour_enable or []):
            fig = _inject_contour(
                fig, df, overlay, value_map,
                contour_mode    = contour_mode or "filled+lines",
                contour_cmap    = contour_cmap or "jet",
                n_levels        = int(contour_levels or 12),
                c_opacity       = float(contour_opacity_pct or 60) / 100.0,
                interp_method   = contour_interp or "cubic",
                extra_factor    = float(contour_extra or 0.12),
                smooth_sigma    = float(contour_smooth or 1.0),
                line_width      = float(contour_lwidth or 1.0),
                contour_opts    = contour_opts or [],
                grid_res        = int(contour_res or 150),
                dark            = dark,
            )

        return fig


# ── Frequency dropdown population ─────────────────────────────────────────────

def _register_freq_options(app) -> None:
    @app.callback(
        Output(IDs.MAP_PAGE_FREQ_SEL, "options"),
        Output(IDs.MAP_PAGE_FREQ_SEL, "value"),
        Input(IDs.STORE_DATA,         "data"),
        Input(IDs.MAP_PAGE_FREQ_UNIT, "value"),
        State(IDs.SESSION_ID,         "data"),
    )
    def populate_freq_options(store_data, freq_unit, session_id):
        if not store_data or not session_id:
            return [], None
        sites = cache_get(session_id)
        if not sites:
            return [], None
        # Collect all unique frequencies across all stations
        freqs: set[float] = set()
        for edi in sites.as_list():
            try:
                fr = np.asarray(edi.Z.freq, dtype=float)
                freqs.update(float(f) for f in fr if np.isfinite(f) and f > 0)
            except Exception:
                pass
        if not freqs:
            return [], None
        sorted_freqs = sorted(freqs, reverse=True)
        use_period = (freq_unit or "hz") == "period"
        if use_period:
            opts = [{"label": f"{1/f:.5g} s", "value": str(f)} for f in sorted_freqs]
        else:
            opts = [{"label": f"{f:.5g} Hz", "value": str(f)} for f in sorted_freqs]
        return opts, opts[0]["value"] if opts else None


# ── Depth info label ──────────────────────────────────────────────────────────

def _register_depth_info(app) -> None:
    @app.callback(
        Output(IDs.MAP_PAGE_DEPTH_INFO, "children"),
        Input(IDs.MAP_PAGE_FREQ_SEL,    "value"),
        Input(IDs.MAP_PAGE_DEPTH_TARGET, "value"),
        State(IDs.SESSION_ID,           "data"),
    )
    def update_depth_info(freq_val, depth_m, session_id):
        if not freq_val or not depth_m or not session_id:
            return ""
        try:
            f    = float(freq_val)
            d    = float(depth_m)
            # Rough estimate: δ = 503 * sqrt(ρ/f) → ρ ≈ (d/503)² * f
            rho  = (d / 503.0) ** 2 * f
            T    = 1.0 / f
            return f"≈ {f:.4g} Hz / T={T:.4g} s  (ρ̄≈{rho:.3g} Ω·m → δ≈{d:.0f} m)"
        except Exception:
            return ""


# ── Inversion result: populate line selector + depth info ────────────────────

def _register_inv_controls(app) -> None:
    @app.callback(
        Output(IDs.MAP_PAGE_INV_LINE_SEL, "options"),
        Output(IDs.MAP_PAGE_INV_LINE_SEL, "value"),
        Input(IDs.STORE_DATA,             "data"),
        State(IDs.SESSION_ID,             "data"),
    )
    def populate_inv_lines(store_data, session_id):
        if not session_id:
            return [], None
        from pycsamt.app.web.cache import cache_get_inversion_result
        result = cache_get_inversion_result(session_id)
        if result is None:
            return [{"label": "No inversion result cached", "value": "none", "disabled": True}], None
        try:
            label = f"{result.method.upper()} {result.dimension.upper()} inversion"
        except Exception:
            label = "Inversion result"
        opts = [{"label": label, "value": "default"}]
        return opts, "default"

    @app.callback(
        Output("map-inv-depth-info", "children"),
        Input(IDs.MAP_PAGE_INV_DEPTH_M, "value"),
        State(IDs.SESSION_ID,           "data"),
    )
    def update_inv_depth_info(depth_m, session_id):
        if not depth_m or not session_id:
            return ""
        from pycsamt.app.web.cache import cache_get_inversion_result
        result = cache_get_inversion_result(session_id)
        if result is None:
            return "No inversion result cached."
        try:
            model = result.to_resistivity_model()
            z_arr = np.asarray(model.z_centers, float)
            iz    = int(np.argmin(np.abs(z_arr - float(depth_m))))
            actual = float(z_arr[iz])
            n_z, n_x = np.asarray(model.rho_2d).shape[:2] if hasattr(model, "rho_2d") else (0, 0)
            return f"Layer {iz + 1}/{n_z}  ·  actual depth {actual:.0f} m  ·  {n_x} stations"
        except Exception as exc:
            return f"⚠ {exc}"


# ── CRS info display ──────────────────────────────────────────────────────────

def _register_crs_info(app) -> None:
    @app.callback(
        Output(IDs.MAP_PAGE_CRS_INFO, "children"),
        Input(IDs.MAP_PAGE_CRS_MODE,  "value"),
        Input(IDs.MAP_PAGE_UTM_ZONE,  "value"),
        Input(IDs.MAP_PAGE_UTM_HEM,   "value"),
        Input(IDs.MAP_PAGE_EPSG,      "value"),
    )
    def update_crs_info(mode, zone, hem, epsg):
        return _resolve_crs_info(mode or "geo",
                                  int(zone or 50), hem or "N", epsg or "4326")


# ── Clientside: show/hide frequency & depth sections ─────────────────────────

def _register_map_type_visibility(app) -> None:
    clientside_callback(
        """
        function(map_type) {
            var freq_style  = {display: 'none'};
            var depth_style = {display: 'none'};
            var inv_style   = {display: 'none'};
            var inv_d_style = {display: 'none'};
            if (map_type === 'resistivity' || map_type === 'depth') {
                freq_style = {display: 'block'};
            }
            if (map_type === 'depth') {
                depth_style = {display: 'block'};
            }
            if (map_type === 'inv_profile' || map_type === 'inv_depth_slice') {
                inv_style = {display: 'block'};
            }
            if (map_type === 'inv_depth_slice') {
                inv_d_style = {display: 'block'};
            }
            return [freq_style, depth_style, inv_style, inv_d_style];
        }
        """,
        Output(IDs.MAP_PAGE_FREQ_SEC,   "style"),
        Output(IDs.MAP_PAGE_DEPTH_SEC,  "style"),
        Output(IDs.MAP_PAGE_INV_SEC,    "style"),
        Output(IDs.MAP_PAGE_INV_DEPTHS, "style"),
        Input(IDs.MAP_PAGE_MAP_TYPE,    "value"),
    )

    # Bearing reset: clicking "N↑" sets the bearing slider back to 0
    clientside_callback(
        """
        function(n) {
            if (!n) return window.dash_clientside.no_update;
            return 0;
        }
        """,
        Output(IDs.MAP_PAGE_BEARING, "value"),
        Input("map-bearing-reset",   "n_clicks"),
        prevent_initial_call=True,
    )


# ── Clientside: show/hide UTM/EPSG sub-sections ──────────────────────────────

def _register_crs_visibility(app) -> None:
    clientside_callback(
        """
        function(mode) {
            var utm   = {display: mode === 'utm'    ? 'block' : 'none'};
            var epsg  = {display: mode === 'custom' ? 'block' : 'none'};
            return [utm, epsg];
        }
        """,
        Output(IDs.MAP_PAGE_UTM_SEC,  "style"),
        Output(IDs.MAP_PAGE_EPSG_SEC, "style"),
        Input(IDs.MAP_PAGE_CRS_MODE,  "value"),
    )


# ── Selection sync (both maps + table) ───────────────────────────────────────

def _register_selection_sync(app, *, include_home_map: bool) -> None:
    callback_inputs = [Input(IDs.MAP_GRAPH, "clickData")]
    if include_home_map:
        callback_inputs.append(Input(IDs.HOME_MAP_GRAPH, "clickData"))
    callback_inputs.extend([
        Input(IDs.STATION_TABLE, "selected_rows"),
        Input({"type": "sta-row", "index": ALL}, "n_clicks"),
    ])

    @app.callback(
        Output(IDs.STORE_SELECTION, "data", allow_duplicate=True),
        *callback_inputs,
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def sync_selection(*args):
        if include_home_map:
            click_data, home_click_data, selected_rows, _row_clicks, store_data = args
        else:
            click_data, selected_rows, _row_clicks, store_data = args
            home_click_data = None

        trigger = ctx.triggered_id
        if trigger == IDs.MAP_GRAPH or (include_home_map and trigger == IDs.HOME_MAP_GRAPH):
            cd = click_data if trigger == IDs.MAP_GRAPH else home_click_data
            if cd:
                pts = cd.get("points", [])
                if pts:
                    sid = pts[0].get("customdata") or pts[0].get("text")
                    return {"station_id": str(sid)}
        if trigger == IDs.STATION_TABLE and selected_rows and store_data:
            records = store_data.get("station_records", [])
            if selected_rows[0] < len(records):
                return {"station_id": records[selected_rows[0]]["ID"]}
        if isinstance(trigger, dict) and trigger.get("type") == "sta-row":
            sid = trigger.get("index")
            if sid:
                return {"station_id": str(sid)}
        return no_update


# ── Contour overlay injector ─────────────────────────────────────────────────

def _inject_contour(
    fig,
    df: "pd.DataFrame",
    overlay: str,
    value_map: "dict | None",
    *,
    contour_mode: str,
    contour_cmap: str,
    n_levels: int,
    c_opacity: float,
    interp_method: str,
    extra_factor: float,
    smooth_sigma: float,
    line_width: float,
    contour_opts: list,
    grid_res: int,
    dark: bool,
):
    """Build a contour overlay and inject it into an existing mapbox figure."""
    import numpy as np

    # Resolve the value array for each station
    lats  = pd.to_numeric(df.get("Latitude",  pd.Series()), errors="coerce").values
    lons  = pd.to_numeric(df.get("Longitude", pd.Series()), errors="coerce").values
    ids   = df["ID"].values if "ID" in df.columns else np.array([])

    # Prefer the pre-computed value_map (resistivity/depth), else use the
    # overlay column already embedded in df.
    if value_map and len(value_map) > 0:
        vals = np.array([float(value_map.get(sid, np.nan)) for sid in ids])
    elif overlay in df.columns:
        vals = pd.to_numeric(df[overlay], errors="coerce").values
    elif overlay == "Elevation" and "Elevation" in df.columns:
        vals = pd.to_numeric(df["Elevation"], errors="coerce").values
    else:
        # Fall back to station index
        vals = np.arange(len(df), dtype=float)

    # Drop rows with missing coordinates or values
    good = (np.isfinite(lats) & np.isfinite(lons) & np.isfinite(vals))
    if good.sum() < 3:
        return fig

    result = build_contour_overlay(
        lats[good], lons[good], vals[good],
        cmap          = contour_cmap,
        n_levels      = max(4, n_levels),
        interp_method = interp_method,
        extra_factor  = extra_factor,
        opacity       = c_opacity,
        mode          = contour_mode,
        show_labels   = "labels" in contour_opts,
        line_width    = line_width,
        smooth_sigma  = smooth_sigma,
        log_scale     = "log" in contour_opts,
        grid_res      = grid_res,
    )

    if result is None:
        return fig

    # ── Inject image layer into map config ──────────────────────────────
    contour_layer = {
        "sourcetype": "image",
        "source":      result["image_b64"],
        "coordinates": result["coordinates"],
        "opacity":     1.0,    # transparency baked into the PNG alpha channel
    }

    # Preserve existing layers (e.g. ESRI tile layer).
    existing = []
    if fig.layout.map and fig.layout.map.layers:
        existing = [ly.to_plotly_json() for ly in fig.layout.map.layers]

    fig.update_layout(map={"layers": existing + [contour_layer]})

    # ── Value-range annotation (replaces second colorbar) ────────────────
    # The station markers already carry the full colorbar for the same
    # quantity.  A compact text stamp in the bottom-left corner gives the
    # contour's min/max without cluttering the figure with a duplicate bar.
    vmin, vmax = result["vmin"], result["vmax"]
    log_scale  = result["log_scale"]
    text_col   = "#cdd6f4" if dark else "#4c4f69"
    bg_rgba    = "rgba(30,30,46,.72)" if dark else "rgba(239,241,245,.80)"

    label = str(overlay).lstrip("_")
    if log_scale:
        stamp = f"Contour (log₁₀) {vmin:.3g} – {vmax:.3g}"
    else:
        stamp = f"Contour  {vmin:.3g} – {vmax:.3g}"

    fig.add_annotation(
        text=stamp,
        xref="paper", yref="paper",
        x=0.01, y=0.01,
        xanchor="left", yanchor="bottom",
        showarrow=False,
        font=dict(size=10, color=text_col),
        bgcolor=bg_rgba,
        borderpad=4,
        opacity=0.9,
    )

    return fig


# ── Clientside: show/hide contour controls ─────────────────────────────────

def _register_contour_visibility(app) -> None:
    clientside_callback(
        """
        function(enable_val) {
            var on = enable_val && enable_val.indexOf('on') >= 0;
            return {display: on ? 'block' : 'none'};
        }
        """,
        Output(IDs.MAP_PAGE_CONTOUR_SEC, "style"),
        Input(IDs.MAP_PAGE_CONTOUR_ENABLE, "value"),
    )


# ── Status bar station label ──────────────────────────────────────────────────

def _register_status_station(app) -> None:
    @app.callback(
        Output(IDs.STATUS_STATION, "children"),
        Input(IDs.STORE_SELECTION, "data"),
        State(IDs.STORE_DATA,      "data"),
    )
    def update_status_station(selection, store_data):
        if not selection:
            return "No station selected"
        station_id = selection.get("station_id", "")
        if not station_id:
            return "No station selected"
        label = f"Station: {station_id}"
        if store_data:
            records = store_data.get("station_records", [])
            match = next((r for r in records if str(r.get("ID")) == str(station_id)), None)
            if match:
                lat = match.get("Latitude", "")
                lon = match.get("Longitude", "")
                if lat and lon:
                    label = f"Station: {station_id}  ({lat:.4f}°N, {lon:.4f}°E)"
        return label


# ── Map toolbar — info chip, toggles, basemap quick-switch, export ────────────

def _register_toolbar(app) -> None:
    from dash import html as _html, clientside_callback

    # ── Station / line count chip ─────────────────────────────────────────────
    @app.callback(
        Output("map-tb-info", "children"),
        Input(IDs.STORE_DATA,  "data"),
    )
    def _tb_info(store_data):
        if not store_data or not store_data.get("n_stations"):
            return _html.Span("No data loaded",
                               className="text-muted", style={"fontSize": "11px"})
        n  = store_data.get("n_stations", 0)
        nl = store_data.get("n_lines",    0)
        return [
            _html.I(className="bi bi-geo-alt-fill me-1",
                    style={"color": "var(--green)", "fontSize": "12px"}),
            _html.Span(f"{n}", className="map-tb-count"),
            _html.Span(" stations", className="map-tb-dim"),
            _html.Span("·", className="map-tb-dot"),
            _html.I(className="bi bi-slash-lg me-1",
                    style={"color": "var(--blue)", "fontSize": "12px"}),
            _html.Span(f"{nl}", className="map-tb-count"),
            _html.Span(" lines", className="map-tb-dim"),
        ]

    # ── Toolbar toggle → Labels checklist (clientside) ─────────────────────────
    clientside_callback(
        """
        function(n, cur) {
            if (!n) return window.dash_clientside.no_update;
            var v = Array.isArray(cur) ? cur.slice() : [];
            var i = v.indexOf("labels");
            if (i > -1) v.splice(i, 1); else v.push("labels");
            return v;
        }
        """,
        Output(IDs.MAP_PAGE_SHOW_LABELS, "value"),
        Input("map-tb-labels", "n_clicks"),
        State(IDs.MAP_PAGE_SHOW_LABELS, "value"),
        prevent_initial_call=True,
    )

    # ── Toolbar toggle → Profiles checklist (clientside) ─────────────────────
    clientside_callback(
        """
        function(n, cur) {
            if (!n) return window.dash_clientside.no_update;
            var v = Array.isArray(cur) ? cur.slice() : [];
            var i = v.indexOf("profiles");
            if (i > -1) v.splice(i, 1); else v.push("profiles");
            return v;
        }
        """,
        Output(IDs.MAP_PAGE_SHOW_LABELS, "value", allow_duplicate=True),
        Input("map-tb-profiles", "n_clicks"),
        State(IDs.MAP_PAGE_SHOW_LABELS, "value"),
        prevent_initial_call=True,
    )

    # ── Toolbar toggle → Contour enable switch (clientside) ───────────────────
    clientside_callback(
        """
        function(n, cur) {
            if (!n) return window.dash_clientside.no_update;
            var v = Array.isArray(cur) ? cur.slice() : [];
            var i = v.indexOf("on");
            if (i > -1) v.splice(i, 1); else v.push("on");
            return v;
        }
        """,
        Output(IDs.MAP_PAGE_CONTOUR_ENABLE, "value"),
        Input("map-tb-contour", "n_clicks"),
        State(IDs.MAP_PAGE_CONTOUR_ENABLE, "value"),
        prevent_initial_call=True,
    )

    # ── Basemap quick-switch (clientside) ──────────────────────────────────────
    clientside_callback(
        """
        function(n1, n2, n3, n4, n5) {
            var ctx = window.dash_clientside.callback_context;
            if (!ctx.triggered || !ctx.triggered.length)
                return window.dash_clientside.no_update;
            var id = ctx.triggered[0].prop_id.split('.')[0];
            var bm = {
                'map-tb-bm-dark':   'carto-darkmatter',
                'map-tb-bm-light':  'carto-positron',
                'map-tb-bm-sat':    'esri-satellite',
                'map-tb-bm-street': 'open-street-map',
                'map-tb-bm-topo':   'esri-topo'
            };
            return bm[id] !== undefined ? bm[id] : window.dash_clientside.no_update;
        }
        """,
        Output(IDs.MAP_PAGE_BASEMAP, "value"),
        Input("map-tb-bm-dark",   "n_clicks"),
        Input("map-tb-bm-light",  "n_clicks"),
        Input("map-tb-bm-sat",    "n_clicks"),
        Input("map-tb-bm-street", "n_clicks"),
        Input("map-tb-bm-topo",   "n_clicks"),
        prevent_initial_call=True,
    )

    # ── Marker-size stepper (clientside) ──────────────────────────────────────
    clientside_callback(
        """
        function(nUp, nDown, curSize) {
            var ctx = window.dash_clientside.callback_context;
            if (!ctx.triggered || !ctx.triggered.length)
                return [window.dash_clientside.no_update,
                        window.dash_clientside.no_update];
            var id  = ctx.triggered[0].prop_id.split('.')[0];
            var cur = typeof curSize === 'number' ? curSize : 10;
            var next = id === 'map-tb-sz-up'
                       ? Math.min(cur + 2, 24)
                       : Math.max(cur - 2, 4);
            return [next, String(next)];
        }
        """,
        Output(IDs.MAP_PAGE_MARKER_SIZE, "value"),
        Output("map-tb-sz-val",          "children"),
        Input("map-tb-sz-up",   "n_clicks"),
        Input("map-tb-sz-down", "n_clicks"),
        State(IDs.MAP_PAGE_MARKER_SIZE, "value"),
        prevent_initial_call=True,
    )

    # ── Active-state sync for layer-toggle buttons (clientside) ───────────────
    clientside_callback(
        """
        function(show_opts, contour_enable) {
            var so = show_opts || [];
            var ce = contour_enable || [];
            var on  = 'map-tb-btn active';
            var off = 'map-tb-btn';
            return [
                so.indexOf('labels')   > -1 ? on : off,
                so.indexOf('profiles') > -1 ? on : off,
                ce.indexOf('on')       > -1 ? on : off
            ];
        }
        """,
        Output("map-tb-labels",  "className"),
        Output("map-tb-profiles","className"),
        Output("map-tb-contour", "className"),
        Input(IDs.MAP_PAGE_SHOW_LABELS,     "value"),
        Input(IDs.MAP_PAGE_CONTOUR_ENABLE,  "value"),
    )

    # ── Active-state sync for basemap quick-switch buttons (clientside) ────────
    clientside_callback(
        """
        function(bm) {
            var map = {
                'carto-darkmatter': 0,
                'carto-positron':   1,
                'esri-satellite':   2,
                'open-street-map':  3,
                'esri-topo':        4
            };
            var clss = ['map-tb-btn','map-tb-btn','map-tb-btn',
                         'map-tb-btn','map-tb-btn'];
            var idx = map[bm || 'carto-darkmatter'];
            if (idx !== undefined) clss[idx] = 'map-tb-btn active';
            return clss;
        }
        """,
        Output("map-tb-bm-dark",   "className"),
        Output("map-tb-bm-light",  "className"),
        Output("map-tb-bm-sat",    "className"),
        Output("map-tb-bm-street", "className"),
        Output("map-tb-bm-topo",   "className"),
        Input(IDs.MAP_PAGE_BASEMAP, "value"),
    )

    # ── Export PNG via Plotly.downloadImage (clientside) ─────────────────────
    clientside_callback(
        """
        function(n) {
            if (!n) return window.dash_clientside.no_update;
            var gd = document.getElementById('map-graph');
            if (gd) {
                Plotly.downloadImage(gd, {
                    format: 'png', filename: 'pycsamt_map',
                    height: 900, width: 1400, scale: 2
                });
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output("map-tb-export", "n_clicks"),
        Input("map-tb-export",  "n_clicks"),
        prevent_initial_call=True,
    )

    # ── Re-centre without zoom-reset: just trigger a soft refresh (clientside)
    clientside_callback(
        """
        function(n, curClicks) {
            if (!n) return window.dash_clientside.no_update;
            return (curClicks || 0) + 1;
        }
        """,
        Output(IDs.MAP_PAGE_REFRESH, "n_clicks"),
        Input("map-tb-center", "n_clicks"),
        State(IDs.MAP_PAGE_REFRESH, "n_clicks"),
        prevent_initial_call=True,
    )

    # ── Force MapLibre resize when map page becomes visible ──────────────────
    # MapLibre initialises at 0×0 when the container is display:none.
    # Dispatching 'resize' + calling Plotly.relayout({autosize:true}) after
    # the container becomes visible forces MapLibre to re-measure and load tiles.
    clientside_callback(
        """
        function(nav_section) {
            if (nav_section !== 'map') return window.dash_clientside.no_update;
            requestAnimationFrame(function() {
                requestAnimationFrame(function() {
                    // Let the display:flex layout settle, then force a full resize
                    window.dispatchEvent(new Event('resize'));
                    var gd = document.getElementById('map-graph');
                    if (gd && window.Plotly) {
                        try {
                            Plotly.relayout(gd, {autosize: true});
                        } catch(e) {}
                    }
                });
            });
            return window.dash_clientside.no_update;
        }
        """,
        Output("map-resize-anchor", "className"),
        Input(IDs.NAV_SECTION, "data"),
    )
