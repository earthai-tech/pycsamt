# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the Inversion Results Viewer page.

Three solver branches:
  • modem    — ModEM 3-D MT   (InversionResult from pycsamt.models.modem)
  • occam2d  — Occam2D 2-D MT (InversionResult from pycsamt.models.occam2d)
  • mare2dem — MARE2DEM 2.5D  (InversionResult from pycsamt.models.mare2dem)

All heavy I/O runs inside the Load callback; the result object is stored in
``_RESULT_CACHE`` keyed by session_id. A lightweight metadata dict is pushed
into the Dash store so the Generate callback can reconstruct params without
re-loading files on every click.
"""
from __future__ import annotations

import base64
import io
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import dash_bootstrap_components as dbc
import matplotlib.pyplot as plt
from dash import (
    ALL,
    Input,
    Output,
    State,
    clientside_callback,
    ctx,
    html,
    no_update,
)
from dash.exceptions import PreventUpdate

from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import PUB_RC, fig_to_src

# ── Tab slugs and solver→tab mapping ─────────────────────────────────────────

_TAB_SLUGS = [
    "conv", "section", "depthmap", "profiles", "covariance",
    "model", "response", "pseudo", "sounding1d", "sitemisfit",
    "respgrid", "survey",
]

_SOLVER_TABS = {
    "modem":    ["conv", "section", "depthmap", "profiles",
                 "covariance", "response", "pseudo"],
    "occam2d":  ["conv", "model", "response", "pseudo",
                 "sounding1d", "sitemisfit", "respgrid"],
    "mare2dem": ["conv", "model", "response", "survey"],
}

# ── Folder browser helpers ────────────────────────────────────────────────────

def _list_dirs(path: str) -> list[str]:
    """Sorted subdirectory names (no hidden) under path."""
    try:
        return sorted(
            e for e in os.listdir(path)
            if not e.startswith(".") and os.path.isdir(os.path.join(path, e))
        )
    except (PermissionError, OSError):
        return []


def _dir_listing(path: str) -> list:
    """Build dbc.ListGroup items for subdirectories at path."""
    dirs = _list_dirs(path)
    if not dirs:
        return [html.Div("No sub-folders here.",
                         style={"color": "var(--sub0)", "fontSize": "11px",
                                "padding": "8px"})]
    items = [
        dbc.ListGroupItem(
            [html.I(className="bi bi-folder-fill me-2",
                    style={"color": "var(--yellow)"}),
             name],
            id={"type": "invr-dir-item", "index": os.path.join(path, name)},
            action=True,
            style={"cursor": "pointer", "fontSize": "12px", "padding": "6px 10px"},
            n_clicks=0,
        )
        for name in dirs
    ]
    return [dbc.ListGroup(items, flush=True)]

# ── Solver detection ──────────────────────────────────────────────────────────

def _detect_solver(path: Path) -> str:
    if not path.is_dir():
        return "unknown"
    names = [f.name.lower() for f in path.iterdir() if f.is_file()]
    if any(n.endswith(".rho") for n in names):
        return "modem"
    if "startup" in names or any(n.endswith(".iter") for n in names):
        return "occam2d"
    if any(n.endswith(".emdata") for n in names):
        return "mare2dem"
    return "unknown"


# ── Loader functions ──────────────────────────────────────────────────────────

def _load_modem(workdir: Path):
    from pycsamt.models.modem.results import InversionResult
    return InversionResult(workdir)


def _load_occam2d(workdir: Path):
    from pycsamt.models.occam2d.results import InversionResult
    return InversionResult(workdir=workdir)


def _load_mare2dem(workdir: Path):
    from pycsamt.models.mare2dem.results import (
        InversionResult,
    )
    return InversionResult(workdir=workdir)


_LOADERS = {
    "modem":   _load_modem,
    "occam2d": _load_occam2d,
    "mare2dem":_load_mare2dem,
}

_DEFAULT_TABS = {
    "modem": "conv", "occam2d": "conv", "mare2dem": "conv",
}

# ── Plot dispatchers ──────────────────────────────────────────────────────────

def _make_fig_modem(result, tab: str, controls: dict) -> plt.Figure:
    from pycsamt.models.modem.plot import (
        PlotAllProfiles,
        PlotCovariance,
        PlotDepthMap,
        PlotMisfit,
        PlotPseudo,
        PlotResponse,
        PlotSection,
    )
    rho_min    = float(controls.get("rho_min",    1.0)   or 1.0)
    rho_max    = float(controls.get("rho_max",    1000.) or 1000.)
    depth_max  = float(controls.get("depth_max",  2000.) or 2000.)
    cmap       = controls.get("cmap",      "jet_r") or "jet_r"
    show_sta   = bool(controls.get("show_sta",   True))
    show_names = bool(controls.get("show_names", True))
    sta_tol    = float(controls.get("sta_tol",   300.)  or 300.)
    which      = controls.get("which",     "final")
    iteration  = controls.get("iteration", None)

    model = None
    if which == "final":
        model = result.model_final
    elif which == "initial":
        model = result.model_initial
    elif which == "iteration" and iteration is not None:
        iter_keys = sorted(
            [k for k in result.models if k.startswith("iter_")],
            key=lambda k: int(k.split("_")[-1]) if k.split("_")[-1].isdigit() else 0,
        )
        if iter_keys:
            idx   = max(0, min(int(iteration) - 1, len(iter_keys) - 1))
            model = result.models[iter_keys[idx]]
    if model is None:
        model = result.model_final or result.model_initial

    if tab == "conv":
        return PlotMisfit(result=result).plot()

    if tab == "section":
        mode = controls.get("sect_mode", "axial")
        if mode == "arbitrary":
            s_lat = controls.get("sect_start_lat")
            s_lon = controls.get("sect_start_lon")
            e_lat = controls.get("sect_end_lat")
            e_lon = controls.get("sect_end_lon")
            o_lat = controls.get("sect_orig_lat")
            o_lon = controls.get("sect_orig_lon")
            n_samp = int(controls.get("sect_n_samp", 200) or 200)
            if None in (s_lat, s_lon, e_lat, e_lon):
                raise ValueError(
                    "Provide all four lat/lon coordinates for arbitrary cut."
                )
            return PlotSection(
                result=result, model=model,
                start_point=(float(s_lat), float(s_lon)),
                end_point  =(float(e_lat), float(e_lon)),
                use_latlon=True,
                origin_lat=float(o_lat) if o_lat is not None else None,
                origin_lon=float(o_lon) if o_lon is not None else None,
                n_samples=n_samp,
                rho_min=rho_min, rho_max=rho_max, depth_max=depth_max,
                cmap=cmap, show_stations=show_sta, show_names=show_names,
                station_tol=sta_tol,
            ).plot()
        direction = controls.get("sect_dir",    "NS")
        offset    = float(controls.get("sect_offset", 0) or 0)
        return PlotSection(
            result=result, model=model,
            direction=direction, offset=offset,
            rho_min=rho_min, rho_max=rho_max, depth_max=depth_max,
            cmap=cmap, show_stations=show_sta, show_names=show_names,
            station_tol=sta_tol,
        ).plot()

    if tab == "depthmap":
        raw    = controls.get("dm_depths", "50,200,500,1000") or "50,200,500,1000"
        depths = [float(v.strip()) for v in raw.split(",") if v.strip()]
        ncols  = int(controls.get("dm_ncols", 2) or 2)
        dm_lat = controls.get("dm_lat")
        dm_lon = controls.get("dm_lon")
        return PlotDepthMap(
            result=result, model=model,
            depths=depths, ncols=ncols,
            origin_lat=float(dm_lat) if dm_lat is not None else None,
            origin_lon=float(dm_lon) if dm_lon is not None else None,
            rho_min=rho_min, rho_max=rho_max,
            cmap=cmap, show_stations=show_sta, show_names=show_names,
        ).plot()

    if tab == "profiles":
        direction  = controls.get("ap_dir", "NS")
        n_prof     = int(controls.get("ap_n", 5) or 5)
        raw_off    = controls.get("ap_offsets", "") or ""
        raw_names  = controls.get("ap_names",   "") or ""
        offsets    = ([float(v.strip()) for v in raw_off.split(",") if v.strip()]
                      if raw_off.strip() else None)
        names      = ([v.strip() for v in raw_names.split(",") if v.strip()]
                      if raw_names.strip() else None)
        ncols      = int(controls.get("ap_ncols", 3) or 3)
        return PlotAllProfiles(
            result=result, model=model,
            direction=direction, n_profiles=n_prof,
            offsets=offsets, profile_names=names, ncols=ncols,
            rho_min=rho_min, rho_max=rho_max, depth_max=depth_max,
            cmap=cmap, show_stations=show_sta, show_names=show_names,
            station_tol=sta_tol,
        ).plot()

    if tab == "covariance":
        return PlotCovariance(
            result=result,
            mask_cmap=controls.get("cov_cmap", "Blues") or "Blues",
        ).plot()

    if tab == "response":
        return PlotResponse(
            result=result,
            station  =controls.get("resp_station") or None,
            component=controls.get("resp_comp", "ZXY") or "ZXY",
        ).plot()

    if tab == "pseudo":
        return PlotPseudo(
            result=result, model=model,
            component=controls.get("ps_comp", "ZXY") or "ZXY",
            rho_min=rho_min, rho_max=rho_max,
        ).plot()

    raise ValueError(f"Unknown ModEM tab: {tab!r}")


def _make_fig_occam2d(result, tab: str, controls: dict) -> plt.Figure:
    from pycsamt.models.occam2d.plot import (
        PlotMisfit,
        PlotModel,
        PlotPseudo,
        PlotResponse,
        PlotResponseGrid,
        PlotSiteMisfit,
        PlotSounding1D,
    )
    rho_min   = float(controls.get("rho_min",   1.0)   or 1.0)
    rho_max   = float(controls.get("rho_max",   1000.) or 1000.)
    depth_max = float(controls.get("depth_max", 2000.) or 2000.)
    cmap      = controls.get("cmap", "jet_r") or "jet_r"

    if tab == "conv":
        return PlotMisfit(result=result).plot()
    if tab == "model":
        return PlotModel(result=result, rho_min=rho_min, rho_max=rho_max,
                         depth_max=depth_max, cmap=cmap).plot()
    if tab == "response":
        return PlotResponse(result=result,
                            station  =controls.get("resp_station") or None,
                            component=controls.get("resp_comp", "ZXY") or "ZXY").plot()
    if tab == "pseudo":
        return PlotPseudo(result=result,
                          component=controls.get("ps_comp", "ZXY") or "ZXY",
                          rho_min=rho_min, rho_max=rho_max).plot()
    if tab == "sounding1d":
        return PlotSounding1D(result=result,
                              station=controls.get("s1d_station") or None).plot()
    if tab == "sitemisfit":
        return PlotSiteMisfit(result=result).plot()
    if tab == "respgrid":
        return PlotResponseGrid(result=result,
                                component=controls.get("rg_comp", "ZXY") or "ZXY").plot()

    raise ValueError(f"Unknown Occam2D tab: {tab!r}")


def _make_fig_mare2dem(result, tab: str, controls: dict) -> plt.Figure:
    from pycsamt.models.mare2dem.plot import (
        PlotConvergence,
        PlotModel,
        PlotResponse,
        PlotSurveyLayout,
    )

    if tab == "conv":
        return PlotConvergence(result).plot()

    if tab == "model":
        fig, ax = plt.subplots(figsize=(8, 5))
        PlotModel(result).plot(ax=ax)
        fig.tight_layout()
        return fig

    if tab == "survey":
        fig, ax = plt.subplots(figsize=(8, 5))
        em_arg = result.em if hasattr(result, "em") else result
        PlotSurveyLayout(em_arg).plot(ax=ax)
        fig.tight_layout()
        return fig

    if tab == "response":
        fig, ax = plt.subplots(figsize=(8, 5))
        PlotResponse(result).plot(ax=ax)
        fig.tight_layout()
        return fig

    raise ValueError(f"Unknown MARE2DEM tab: {tab!r}")


_DISPATCHERS = {
    "modem":   _make_fig_modem,
    "occam2d": _make_fig_occam2d,
    "mare2dem":_make_fig_mare2dem,
}

# ── In-process result cache ───────────────────────────────────────────────────

_RESULT_CACHE: dict[str, object] = {}


def _store_result(session_id: str, result) -> None:
    _RESULT_CACHE[f"invr::{session_id}"] = result


def _fetch_result(session_id: str):
    return _RESULT_CACHE.get(f"invr::{session_id}")


# ── Info strip builder ────────────────────────────────────────────────────────

def _build_info_strip(result, solver: str, path: str) -> list:
    from dash import html

    def _badge(text, color="secondary"):
        return html.Span(text, className=f"badge bg-{color} me-1")

    items = [
        _badge(solver.upper(), "primary"),
        html.Span(Path(path).name, className="adv-info-group me-2 text-muted"),
    ]
    try:
        if solver == "modem":
            mdl  = result.model_final or result.model_initial
            grid = f"{mdl.nx}×{mdl.ny}×{mdl.nz}" if mdl else "?"
            rms  = result.final_rms
            items += [
                _badge(f"mode={result.mode}", "dark"),
                _badge(f"grid={grid}", "dark"),
                _badge(f"iter={result.n_iter}", "info"),
                _badge(f"RMS={rms:.3f}", "warning"),
            ]
        elif solver == "occam2d":
            rms = result.final_rms
            items += [
                _badge(f"iter={result.n_iterations}", "info"),
                _badge(f"RMS={rms:.3f}", "warning"),
            ]
        elif solver == "mare2dem":
            rms_val = result.final_rms
            rms_str = f"{rms_val:.3f}" if rms_val is not None else "?"
            items += [
                _badge(f"iter={result.n_iterations}", "info"),
                _badge(f"RMS={rms_str}", "warning"),
            ]
    except Exception:
        pass
    return items


# ── Shared helper for building the controls dict ──────────────────────────────

def _build_controls(*args) -> dict:
    (which, iteration, rho_min, rho_max, depth_max, cmap,
     sect_mode, sect_dir, sect_offset,
     sect_start_lat, sect_start_lon, sect_end_lat, sect_end_lon,
     sect_n_samp, sect_orig_lat, sect_orig_lon,
     dm_depths, dm_ncols, dm_lat, dm_lon,
     ap_dir, ap_n, ap_offsets, ap_names, ap_ncols,
     cov_cmap, resp_station, resp_comp, ps_comp,
     show_sta, show_names, sta_tol) = args
    return dict(
        which=which, iteration=iteration,
        rho_min=rho_min, rho_max=rho_max,
        depth_max=depth_max, cmap=cmap,
        sect_mode=sect_mode, sect_dir=sect_dir, sect_offset=sect_offset,
        sect_start_lat=sect_start_lat, sect_start_lon=sect_start_lon,
        sect_end_lat=sect_end_lat,     sect_end_lon=sect_end_lon,
        sect_n_samp=sect_n_samp,
        sect_orig_lat=sect_orig_lat,   sect_orig_lon=sect_orig_lon,
        dm_depths=dm_depths, dm_ncols=dm_ncols,
        dm_lat=dm_lat, dm_lon=dm_lon,
        ap_dir=ap_dir, ap_n=ap_n, ap_offsets=ap_offsets,
        ap_names=ap_names, ap_ncols=ap_ncols,
        cov_cmap=cov_cmap,
        resp_station=resp_station, resp_comp=resp_comp, ps_comp=ps_comp,
        show_sta=bool(show_sta), show_names=bool(show_names),
        sta_tol=sta_tol,
    )


# ── Common State list (re-used in both Generate and Export) ───────────────────

def _control_states():
    return [
        State(IDs.INVR_WHICH,          "value"),
        State(IDs.INVR_ITERATION,      "value"),
        State(IDs.INVR_RHO_MIN,        "value"),
        State(IDs.INVR_RHO_MAX,        "value"),
        State(IDs.INVR_DEPTH_MAX,      "value"),
        State(IDs.INVR_CMAP,           "value"),
        State(IDs.INVR_SECT_MODE,      "value"),
        State(IDs.INVR_SECT_DIR,       "value"),
        State(IDs.INVR_SECT_OFFSET,    "value"),
        State(IDs.INVR_SECT_START_LAT, "value"),
        State(IDs.INVR_SECT_START_LON, "value"),
        State(IDs.INVR_SECT_END_LAT,   "value"),
        State(IDs.INVR_SECT_END_LON,   "value"),
        State(IDs.INVR_SECT_N_SAMP,    "value"),
        State(IDs.INVR_SECT_ORIG_LAT,  "value"),
        State(IDs.INVR_SECT_ORIG_LON,  "value"),
        State(IDs.INVR_DM_DEPTHS,      "value"),
        State(IDs.INVR_DM_NCOLS,       "value"),
        State(IDs.INVR_DM_LAT,         "value"),
        State(IDs.INVR_DM_LON,         "value"),
        State(IDs.INVR_AP_DIR,         "value"),
        State(IDs.INVR_AP_N,           "value"),
        State(IDs.INVR_AP_OFFSETS,     "value"),
        State(IDs.INVR_AP_NAMES,       "value"),
        State(IDs.INVR_AP_NCOLS,       "value"),
        State(IDs.INVR_COV_CMAP,       "value"),
        State(IDs.INVR_RESP_STATION,   "value"),
        State(IDs.INVR_RESP_COMP,      "value"),
        State(IDs.INVR_PS_COMP,        "value"),
        State(IDs.INVR_SHOW_STA,       "value"),
        State(IDs.INVR_SHOW_NAMES,     "value"),
        State(IDs.INVR_STA_TOL,        "value"),
    ]


# ── Register all callbacks on the app ────────────────────────────────────────

def register_inv_results(app) -> None:

    # ── 1. Tab button clicks → update active tab store ────────────────────────
    # Uses arguments() in JS so we don't need to name all 12 inputs explicitly.
    clientside_callback(
        """
        function() {
            var triggered = window.dash_clientside.callback_context.triggered;
            if (!triggered || !triggered.length) {
                return window.dash_clientside.no_update;
            }
            var prop_id = triggered[0].prop_id;
            var match = prop_id.match(/invr-tab-btn-([^.]+)/);
            return match ? match[1] : window.dash_clientside.no_update;
        }
        """,
        Output(IDs.INVR_ACTIVE_TAB, "data"),
        [Input(f"invr-tab-btn-{s}", "n_clicks") for s in _TAB_SLUGS],
        prevent_initial_call=True,
    )

    # ── 2. Active tab → update tab button classes ─────────────────────────────
    clientside_callback(
        """
        function(active_tab) {
            var slugs = """ + str(_TAB_SLUGS).replace("'", '"') + """;
            return slugs.map(function(s) {
                return s === active_tab ? "prof-tab-btn active" : "prof-tab-btn";
            });
        }
        """,
        [Output(f"invr-tab-btn-{s}", "className") for s in _TAB_SLUGS],
        Input(IDs.INVR_ACTIVE_TAB, "data"),
        prevent_initial_call=False,
    )

    # ── 3. Active tab → show/hide context panels ──────────────────────────────
    clientside_callback(
        """
        function(active_tab) {
            var slugs = """ + str(_TAB_SLUGS).replace("'", '"') + """;
            return slugs.map(function(s) {
                return s === active_tab ? {display:"block"} : {display:"none"};
            });
        }
        """,
        [Output(f"invr-ctx-{s}", "style") for s in _TAB_SLUGS],
        Input(IDs.INVR_ACTIVE_TAB, "data"),
        prevent_initial_call=False,
    )

    # ── 4. Section mode radio → show/hide axial vs arbitrary sub-panels ───────
    clientside_callback(
        """
        function(mode) {
            if (mode === "arbitrary") {
                return [{display:"none"}, {display:"block"}];
            }
            return [{display:"block"}, {display:"none"}];
        }
        """,
        [Output(IDs.INVR_SECT_AX_PANEL,  "style"),
         Output(IDs.INVR_SECT_ARB_PANEL, "style")],
        Input(IDs.INVR_SECT_MODE, "value"),
        prevent_initial_call=False,
    )

    # ── 5. Spinner show on load-button click (clientside) ─────────────────────
    clientside_callback(
        """
        function(n) {
            return n ? {display:"flex", alignItems:"center", justifyContent:"center"}
                     : {display:"none"};
        }
        """,
        Output("invr-load-spinner", "style"),
        Input(IDs.BTN_INVR_LOAD, "n_clicks"),
        prevent_initial_call=False,
    )

    # ── 6. Solver panel show/hide (clientside) ─────────────────────────────────
    clientside_callback(
        """
        function(solver_radio, meta) {
            var solver = (meta && meta.solver) ? meta.solver : (solver_radio || "");
            return [
                solver === "modem"    ? {display:"block"} : {display:"none"},
                solver === "occam2d"  ? {display:"block"} : {display:"none"},
                solver === "mare2dem" ? {display:"block"} : {display:"none"},
            ];
        }
        """,
        [Output("invr-acc-modem-wrap",   "style"),
         Output("invr-acc-occam2d-wrap",  "style"),
         Output("invr-acc-mare2dem-wrap", "style")],
        Input(IDs.INVR_SOLVER_RADIO, "value"),
        Input(IDs.INVR_STORE,        "data"),
        prevent_initial_call=False,
    )

    # ── 7. Solver tabs show/hide (clientside) ──────────────────────────────────
    _solver_tabs_js = {k: list(v) for k, v in _SOLVER_TABS.items()}
    clientside_callback(
        f"""
        function(meta) {{
            var solver_tabs = {json.dumps(_solver_tabs_js)};
            var all_slugs   = {json.dumps(_TAB_SLUGS)};
            var solver = meta ? meta.solver : null;
            var valid  = solver ? (solver_tabs[solver] || all_slugs) : all_slugs;
            return all_slugs.map(function(s) {{
                return valid.indexOf(s) >= 0
                       ? {{display:"inline-flex"}}
                       : {{display:"none"}};
            }});
        }}
        """,
        [Output(f"invr-tab-btn-{s}", "style") for s in _TAB_SLUGS],
        Input(IDs.INVR_STORE, "data"),
        prevent_initial_call=False,
    )

    # ── 8. Folder browser — open / cancel modal ─────────────────────────────────
    @app.callback(
        Output("invr-browser-modal", "is_open"),
        Output("invr-browser-path",  "data"),
        Input("btn-invr-browse",         "n_clicks"),
        Input("btn-invr-browser-cancel", "n_clicks"),
        State("invr-browser-modal",      "is_open"),
        State("invr-browser-path",       "data"),
        State(IDs.INVR_PATH_INPUT,       "value"),
        prevent_initial_call=True,
    )
    def _browser_toggle(n_browse, n_cancel, is_open, cur_path, typed_path):
        if ctx.triggered_id == "btn-invr-browse" and n_browse:
            start = os.path.expanduser("~")
            if typed_path:
                candidate = typed_path.strip()
                candidate = candidate if os.path.isdir(candidate) \
                            else os.path.dirname(candidate)
                if os.path.isdir(candidate):
                    start = candidate
            elif cur_path and os.path.isdir(cur_path):
                start = cur_path
            return True, start
        return False, cur_path

    # ── 9. Folder browser — navigate into subfolder ─────────────────────────────
    @app.callback(
        Output("invr-browser-path", "data", allow_duplicate=True),
        Input({"type": "invr-dir-item", "index": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def _browser_navigate(n_clicks_list):
        if not any(n for n in n_clicks_list if n):
            raise PreventUpdate
        triggered = ctx.triggered_id
        if triggered and isinstance(triggered, dict):
            new_path = triggered.get("index", "")
            if new_path and os.path.isdir(new_path):
                return new_path
        raise PreventUpdate

    # ── 10. Folder browser — go up one level ───────────────────────────────────
    @app.callback(
        Output("invr-browser-path", "data", allow_duplicate=True),
        Input("btn-invr-browser-up", "n_clicks"),
        State("invr-browser-path",   "data"),
        prevent_initial_call=True,
    )
    def _browser_up(n_clicks, cur_path):
        if not n_clicks or not cur_path:
            raise PreventUpdate
        parent = os.path.dirname(cur_path)
        if parent and parent != cur_path:
            return parent
        raise PreventUpdate

    # ── 11. Folder browser — refresh listing when path changes ─────────────────
    @app.callback(
        Output("invr-browser-list",    "children"),
        Output("invr-browser-cwd",     "children"),
        Input("invr-browser-path",     "data"),
    )
    def _browser_refresh(path_str):
        path = path_str or os.path.expanduser("~")
        return _dir_listing(path), path

    # ── 12. Folder browser — select current folder ──────────────────────────────
    @app.callback(
        Output(IDs.INVR_PATH_INPUT,      "value",   allow_duplicate=True),
        Output("invr-browser-modal",     "is_open", allow_duplicate=True),
        Input("btn-invr-browser-select", "n_clicks"),
        State("invr-browser-path",       "data"),
        prevent_initial_call=True,
    )
    def _browser_select(n_clicks, path_str):
        if not n_clicks:
            raise PreventUpdate
        return (path_str or no_update), False

    # ── 8. Load callback ───────────────────────────────────────────────────────
    @app.callback(
        Output(IDs.INVR_STORE,        "data"),
        Output(IDs.INVR_STATUS,       "children"),
        Output(IDs.INVR_INFO_STRIP,   "children"),
        Output(IDs.INVR_RESP_STATION, "options"),
        Output(IDs.INVR_ACTIVE_TAB,   "data",        allow_duplicate=True),
        Output("invr-load-spinner",   "style",        allow_duplicate=True),
        Input(IDs.BTN_INVR_LOAD,      "n_clicks"),
        State(IDs.INVR_PATH_INPUT,    "value"),
        State(IDs.INVR_SOLVER_RADIO,  "value"),
        State(IDs.SESSION_ID,         "data"),
        prevent_initial_call=True,
    )
    def _load_result(n_clicks, path_str, solver_hint, session_id):
        from dash import html

        _spinner_off = {"display": "none"}

        if not n_clicks or not path_str:
            raise PreventUpdate

        path = Path(path_str.strip())
        if not path.exists():
            return (
                no_update,
                [html.I(className="bi bi-x-circle me-1 text-danger"),
                 f"Path not found: {path_str}"],
                no_update, [], no_update, _spinner_off,
            )

        solver = solver_hint if solver_hint != "auto" else _detect_solver(path)
        if solver == "unknown" or solver not in _LOADERS:
            return (
                no_update,
                [html.I(className="bi bi-question-circle me-1 text-warning"),
                 f"Cannot detect solver in: {path.name}"],
                no_update, [], no_update, _spinner_off,
            )

        try:
            result = _LOADERS[solver](path)
        except Exception as exc:
            return (
                no_update,
                [html.I(className="bi bi-exclamation-triangle me-1 text-danger"),
                 f"Load error: {exc}"],
                no_update, [], no_update, _spinner_off,
            )

        _store_result(session_id or "default", result)

        # Station dropdown options
        station_opts = []
        try:
            if solver == "modem":
                data = result.data_obs or result.data_pred
                if data is not None and hasattr(data, "stations"):
                    station_opts = [{"label": s, "value": s}
                                    for s in sorted(data.stations)]
            elif solver == "occam2d":
                if hasattr(result, "response") and result.response is not None:
                    stas = getattr(result.response, "stations", [])
                    station_opts = [{"label": s, "value": s}
                                    for s in sorted(stas)]
        except Exception:
            pass

        n_iter  = int(getattr(result, "n_iter",
                      getattr(result, "n_iterations", 0)) or 0)
        rms_val = getattr(result, "final_rms", float("nan"))
        try:
            rms_f = float(rms_val or float("nan"))
        except (TypeError, ValueError):
            rms_f = float("nan")

        meta = {
            "path":   str(path),
            "solver": solver,
            "n_iter": n_iter,
            "rms":    rms_f,
        }
        rms_str = f"{rms_f:.3f}" if not __import__("math").isnan(rms_f) else "?"
        status  = [
            html.I(className="bi bi-check-circle me-1 text-success"),
            f"{solver.upper()} loaded · {n_iter} iter · RMS={rms_str}",
        ]
        info = _build_info_strip(result, solver, str(path))
        tab  = _DEFAULT_TABS.get(solver, "conv")

        return meta, status, info, station_opts, tab, _spinner_off

    # ── 13. Generate plot (auto on tab switch or manual Generate click) ──────────
    @app.callback(
        Output(IDs.INVR_CANVAS, "src"),
        Input(IDs.BTN_INVR_GEN,    "n_clicks"),
        Input(IDs.INVR_ACTIVE_TAB, "data"),    # ← auto-fire on tab change
        State(IDs.INVR_STORE,      "data"),
        State(IDs.SESSION_ID,      "data"),
        *_control_states(),
        prevent_initial_call=True,
    )
    def _generate_plot(n_clicks, active_tab, meta, session_id, *ctrl_args):
        if not meta:
            raise PreventUpdate

        result = _fetch_result(session_id or "default")
        if result is None:
            raise PreventUpdate

        solver     = meta.get("solver", "modem")
        tab        = active_tab or "conv"
        controls   = _build_controls(*ctrl_args)
        dispatcher = _DISPATCHERS.get(solver)
        if dispatcher is None:
            raise PreventUpdate

        try:
            fig = dispatcher(result, tab, controls)
            src = fig_to_src(fig, dpi=130)
            plt.close(fig)
            return src
        except Exception as exc:
            plt.close("all")
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.text(0.5, 0.5, f"Plot error:\n{exc}",
                    ha="center", va="center", transform=ax.transAxes,
                    color="red", fontsize=10, wrap=True)
            ax.axis("off")
            src = fig_to_src(fig, dpi=90)
            plt.close(fig)
            return src

    # ── 10. Export plot as 300 dpi PNG ────────────────────────────────────────
    @app.callback(
        Output(IDs.INVR_DOWNLOAD, "data"),
        Input(IDs.BTN_INVR_EXPORT, "n_clicks"),
        State(IDs.INVR_STORE,      "data"),
        State(IDs.INVR_ACTIVE_TAB, "data"),
        State(IDs.SESSION_ID,      "data"),
        *_control_states(),
        prevent_initial_call=True,
    )
    def _export_plot(n_clicks, meta, active_tab, session_id, *ctrl_args):
        if not n_clicks or not meta:
            raise PreventUpdate

        result = _fetch_result(session_id or "default")
        if result is None:
            raise PreventUpdate

        solver     = meta.get("solver", "modem")
        tab        = active_tab or "conv"
        controls   = _build_controls(*ctrl_args)
        dispatcher = _DISPATCHERS.get(solver)
        if dispatcher is None:
            raise PreventUpdate

        try:
            # Re-generate with publication-safe white background
            with matplotlib.rc_context(PUB_RC):
                fig = dispatcher(result, tab, controls)
            buf = io.BytesIO()
            fig.savefig(
                buf,
                format="png",
                dpi=300,
                bbox_inches="tight",
                facecolor="white",
                edgecolor="none",
            )
            plt.close(fig)
            buf.seek(0)
            b64      = base64.b64encode(buf.read()).decode()
            filename = f"invr_{solver}_{tab}.png"
            return dict(content=b64, filename=filename, base64=True)
        except Exception:
            plt.close("all")
            raise PreventUpdate
