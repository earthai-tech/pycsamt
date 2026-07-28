# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Callbacks for the Advanced Plots page.

Tab structure (tab-id → plot group)
------------------------------------
strike    → Strike Analysis    (STRIKE_PLOTS)
pt        → Phase Tensor       (PHASE_TENSOR_PLOTS)
induction → Induction / Tipper (INDUCTION_PLOTS)
impedance → Impedance / Z      (IMPEDANCE_PLOTS)
depth     → Depth Imaging      (DEPTH_PLOTS)
survey    → Survey Tools       (SURVEY_PLOTS)

Each tab keeps its own persistent figure (img-adv-{tab}).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dash import (
    Input,
    Output,
    State,
    clientside_callback,
    ctx,
    html,
    no_update,
)
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.advanced_controller import (
    ADVANCED_PLOT_DESCRIPTIONS,
    DEPTH_PLOTS,
    IMPEDANCE_PLOTS,
    INDUCTION_PLOTS,
    PHASE_TENSOR_PLOTS,
    STRIKE_PLOTS,
    SURVEY_PLOTS,
    AdvancedController,
)
from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    fig_to_src,
)

# Singleton controller
_CTRL = AdvancedController()

_ADV_TAB_IDS = ["strike", "pt", "induction", "impedance", "depth", "survey"]
_ADV_TAB_ORDER = {t: i for i, t in enumerate(_ADV_TAB_IDS)}

_TAB_TO_PLOTS = {
    "strike": STRIKE_PLOTS,
    "pt": PHASE_TENSOR_PLOTS,
    "induction": INDUCTION_PLOTS,
    "impedance": IMPEDANCE_PLOTS,
    "depth": DEPTH_PLOTS,
    "survey": SURVEY_PLOTS,
}

_TAB_TO_GROUP = {
    "strike": "Strike Analysis",
    "pt": "Phase Tensor",
    "induction": "Induction / Tipper",
    "impedance": "Impedance / Z",
    "depth": "Depth Imaging",
    "survey": "Survey Tools",
}

_TAB_ICONS = {
    "strike": "bi-compass",
    "pt": "bi-hexagon",
    "induction": "bi-arrows-expand",
    "impedance": "bi-lightning-charge",
    "depth": "bi-layers",
    "survey": "bi-clipboard-data",
}

# Build fn → (label, has_ax) lookup from all groups
_ALL_PLOTS: list = (
    STRIKE_PLOTS
    + PHASE_TENSOR_PLOTS
    + INDUCTION_PLOTS
    + IMPEDANCE_PLOTS
    + DEPTH_PLOTS
    + SURVEY_PLOTS
)
_HAS_AX: dict[str, bool] = {fn: has_ax for _, fn, has_ax in _ALL_PLOTS}
_FN_TO_LBL: dict[str, str] = {fn: lbl for lbl, fn, _ in _ALL_PLOTS}

# Plots that expose a period-index parameter (ip=)
_PERIOD_FNS = {
    "plot_phase_tensor_map",
    "plot_induction_arrows",
    "plot_induction_map",
    "plot_xyyx_crossover_map",
    "plot_dim_map",
    "plot_dim_occupancy_area",
    "plot_anisotropy",
    "plot_ellipticity_psection",
    "plot_depth_section",
    "plot_apparent_depth_psection",
    "plot_gradient_section",
    "plot_induction_section",
    "plot_theta_vs_period",
    "plot_theta_stability_stripe",
    "plot_strike_ribbon",
    "plot_strike_mapsticks",
}

_FIGSIZE_MAP = {
    "8x4": (8, 4),
    "10x6": (10, 6),
    "14x6": (14, 6),
    "10x8": (10, 8),
}

# Per-tab image IDs (same order as _ADV_TAB_IDS)
_TAB_IMG_IDS = [
    IDs.IMG_ADV_STRIKE,
    IDs.IMG_ADV_PT,
    IDs.IMG_ADV_INDUCTION,
    IDs.IMG_ADV_IMPEDANCE,
    IDs.IMG_ADV_DEPTH,
    IDs.IMG_ADV_SURVEY,
]


def _filter_sites(sites, station_filter, store_data, active_lines_store):
    """Return Sites limited to requested stations (priority: explicit > active lines > all)."""
    if sites is None:
        return None
    records = (store_data or {}).get("station_records", [])
    active_lines = (active_lines_store or {}).get("active", [])
    if station_filter:
        wanted = set(station_filter)
    elif active_lines:
        active_set = set(active_lines)
        wanted = {
            r["ID"]
            for r in records
            if r.get("Line", "") in active_set or not r.get("Line")
        }
    else:
        return sites
    try:
        from pycsamt.site.base import Sites

        filtered = [e for e in sites.edic if getattr(e, "station", None) in wanted]
        return Sites(filtered) if filtered else sites
    except Exception:
        return sites


def _pt_grid_profiles(store_data, active_lines_store, per_line: int) -> dict:
    """Group ``station_records`` by ``Line`` (active lines only) into
    ``{line_name: [representative_station, ...]}`` for
    ``plot_phase_tensor_strip_grid``."""
    from pycsamt.site.lines import (
        pick_representative_stations,
    )

    records = (store_data or {}).get("station_records", [])
    active_lines = (active_lines_store or {}).get("active", [])
    active_set = set(active_lines) if active_lines else None
    lines: dict[str, list] = {}
    for r in records:
        ln = r.get("Line", "")
        if not ln or (active_set is not None and ln not in active_set):
            continue
        lines.setdefault(ln, []).append(r["ID"])
    return {
        ln: pick_representative_stations(names, per_line) for ln, names in lines.items()
    }


def register_advanced(app) -> None:

    # ── 1. Tab button clicks → ADV_ACTIVE_TAB store ───────────────────────────
    @app.callback(
        Output(IDs.ADV_ACTIVE_TAB, "data"),
        [Input(f"adv-tab-{t}", "n_clicks") for t in _ADV_TAB_IDS],
        prevent_initial_call=True,
    )
    def switch_adv_tab(*_clicks):
        tid = ctx.triggered_id
        if not tid:
            return no_update
        return tid.replace("adv-tab-", "")

    # ── 2. Store → panel visibility + button classes (clientside) ─────────────
    clientside_callback(
        """
        function(active_tab) {
            if (!active_tab) return Array(12).fill(window.dash_clientside.no_update);
            var tabs = ["strike","pt","induction","impedance","depth","survey"];
            var flex = {display:"flex",flexDirection:"column",
                        height:"100%",minHeight:"0",flex:"1"};
            var none = {display:"none"};
            var styles  = tabs.map(function(t){ return t===active_tab ? flex : none; });
            var classes = tabs.map(function(t){
                return "prof-tab-btn" + (t===active_tab ? " active" : "");
            });
            return styles.concat(classes);
        }
        """,
        Output("adv-panel-strike", "style"),
        Output("adv-panel-pt", "style"),
        Output("adv-panel-induction", "style"),
        Output("adv-panel-impedance", "style"),
        Output("adv-panel-depth", "style"),
        Output("adv-panel-survey", "style"),
        Output("adv-tab-strike", "className"),
        Output("adv-tab-pt", "className"),
        Output("adv-tab-induction", "className"),
        Output("adv-tab-impedance", "className"),
        Output("adv-tab-depth", "className"),
        Output("adv-tab-survey", "className"),
        Input(IDs.ADV_ACTIVE_TAB, "data"),
        prevent_initial_call=False,
    )

    # ── 3. Tab change → plot dropdown (replaces old ADV_GROUP sync) ───────────
    @app.callback(
        Output(IDs.ADV_PLOT, "options"),
        Output(IDs.ADV_PLOT, "value"),
        Input(IDs.ADV_ACTIVE_TAB, "data"),
    )
    def sync_plots(active_tab):
        if not active_tab:
            raise PreventUpdate
        plots = _TAB_TO_PLOTS.get(active_tab, STRIKE_PLOTS)
        opts = [{"label": lbl, "value": fn} for lbl, fn, _ in plots]
        return opts, (opts[0]["value"] if opts else None)

    # ── 4. Plot selection → description + period/ATOM/PT-strip visibility ────
    @app.callback(
        Output(IDs.ADV_DESCRIPTION, "children"),
        Output("adv-period-section", "style"),
        Output("adv-atom-section", "style"),
        Output("adv-pt-station-section", "style"),
        Output("adv-pt-grid-section", "style"),
        Input(IDs.ADV_PLOT, "value"),
    )
    def update_param_panel(fn_name):
        _hidden, _shown = {"display": "none"}, {"display": "block"}
        if not fn_name:
            return "", _hidden, _hidden, _hidden, _hidden
        desc = ADVANCED_PLOT_DESCRIPTIONS.get(fn_name, "")
        period_style = _shown if fn_name in _PERIOD_FNS else _hidden
        atom_style = _shown if fn_name == "plot_atom_psection" else _hidden
        pt_sta_style = _shown if fn_name == "plot_phase_tensor_strip" else _hidden
        pt_grid_style = _shown if fn_name == "plot_phase_tensor_strip_grid" else _hidden
        return desc, period_style, atom_style, pt_sta_style, pt_grid_style

    # ── 5. Context info bar ───────────────────────────────────────────────────
    @app.callback(
        Output(IDs.ADV_INFO_BAR, "children"),
        Input(IDs.ADV_ACTIVE_TAB, "data"),
        Input(IDs.ADV_PLOT, "value"),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=False,
    )
    def update_info_bar(active_tab, fn_name, store_data):
        group = _TAB_TO_GROUP.get(active_tab or "strike", "Strike Analysis")
        icon = _TAB_ICONS.get(active_tab or "strike", "bi-graph-up")
        label = _FN_TO_LBL.get(fn_name or "", "")
        n_stn = (store_data or {}).get("n_stations", 0)

        parts = [
            html.Span(
                [html.I(className=f"bi {icon} me-1"), group],
                className="adv-info-group",
            ),
        ]
        if label:
            parts += [
                html.Span("·", className="adv-info-sep"),
                html.Span(label, className="adv-info-plot"),
            ]
        if n_stn:
            parts += [
                html.Span("·", className="adv-info-sep"),
                html.Span(
                    f"{n_stn} station{'s' if n_stn != 1 else ''}",
                    className="adv-info-stat",
                ),
            ]
        if not label:
            parts.append(
                html.Span(
                    " · Select a plot and click Generate",
                    style={"color": "var(--sub0)", "fontSize": "10.5px"},
                ),
            )
        return parts

    # ── 6. STORE_DATA + STORE_ACTIVE_LINES → filters + data bar ──────────────
    @app.callback(
        Output(IDs.ADV_LINE_FILTER, "options"),
        Output(IDs.ADV_STN_FILTER, "options"),
        Output(IDs.ADV_DATA_BAR, "children"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=False,
    )
    def update_filters(store_data, active_lines_store):
        if not store_data or not store_data.get("n_stations"):
            bar = html.Div("No survey loaded", className="adv-bar-empty")
            return [], [], bar

        records = store_data.get("station_records", [])
        active_lines = (active_lines_store or {}).get("active", [])
        if not active_lines:
            active_lines = list({r.get("Line", "") for r in records if r.get("Line")})

        line_opts = [{"label": ln, "value": ln} for ln in active_lines if ln]

        active_set = set(active_lines)
        stn_opts = [
            {"label": r["ID"], "value": r["ID"]}
            for r in records
            if r.get("Line", "") in active_set or not r.get("Line")
        ]

        n_s = len(stn_opts)
        n_l = len(active_lines)
        bar = html.Div(
            [
                html.Span(
                    f"{n_s} station{'s' if n_s != 1 else ''}",
                    className="adv-bar-chip",
                ),
                html.Span(
                    f"{n_l} active line{'s' if n_l != 1 else ''}",
                    className="adv-bar-chip",
                ),
            ],
            className="adv-bar",
        )
        return line_opts, stn_opts, bar

    # ── 7. Line filter → station scope ───────────────────────────────────────
    @app.callback(
        Output(IDs.ADV_STN_FILTER, "value", allow_duplicate=True),
        Output(IDs.ADV_STN_FILTER, "options", allow_duplicate=True),
        Input(IDs.ADV_LINE_FILTER, "value"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def line_to_stations(selected_lines, store_data, active_lines_store):
        if not store_data:
            return no_update, no_update
        records = store_data.get("station_records", [])
        active_lines = (active_lines_store or {}).get("active", [])
        if not active_lines:
            active_lines = list({r.get("Line", "") for r in records if r.get("Line")})
        active_set = set(active_lines)
        stn_opts = [
            {"label": r["ID"], "value": r["ID"]}
            for r in records
            if r.get("Line", "") in active_set or not r.get("Line")
        ]
        if not selected_lines:
            return None, stn_opts
        sel_set = set(selected_lines)
        selected = [r["ID"] for r in records if r.get("Line", "") in sel_set]
        return selected or None, stn_opts

    # ── 7b. STORE_DATA / active lines → PT-strip single-station options ──────
    @app.callback(
        Output(IDs.ADV_PT_STATION, "options"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=False,
    )
    def update_pt_station_options(store_data, active_lines_store):
        if not store_data:
            return []
        records = store_data.get("station_records", [])
        active_lines = (active_lines_store or {}).get("active", [])
        active_set = set(active_lines) if active_lines else None
        return [
            {"label": r["ID"], "value": r["ID"]}
            for r in records
            if active_set is None
            or r.get("Line", "") in active_set
            or not r.get("Line")
        ]

    # ── 8. Station All / None buttons ─────────────────────────────────────────
    @app.callback(
        Output(IDs.ADV_STN_FILTER, "value", allow_duplicate=True),
        Input(IDs.ADV_STN_ALL, "n_clicks"),
        State(IDs.ADV_STN_FILTER, "options"),
        prevent_initial_call=True,
    )
    def station_select_all(_n, opts):
        if not _n or not opts:
            return no_update
        return [o["value"] for o in opts]

    @app.callback(
        Output(IDs.ADV_STN_FILTER, "value", allow_duplicate=True),
        Input(IDs.ADV_STN_NONE, "n_clicks"),
        prevent_initial_call=True,
    )
    def station_deselect_all(_n):
        return [] if _n else no_update

    # ── 9. Train ATOM ─────────────────────────────────────────────────────────
    @app.callback(
        Output(IDs.ADV_TRAIN_SPINNER, "children"),
        Output(IDs.ADV_TRAIN_STATUS, "children"),
        Input(IDs.ADV_BTN_TRAIN, "n_clicks"),
        State(IDs.SESSION_ID, "data"),
        State("adv-atom-atoms", "value"),
        State("adv-atom-iter", "value"),
        State(IDs.ADV_STN_FILTER, "value"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def train_atom(
        n_clicks,
        session_id,
        n_atoms,
        n_iter,
        station_filter,
        store_data,
        active_lines_store,
    ):
        if not n_clicks:
            raise PreventUpdate
        try:
            sites = cache_get(session_id)
            filtered = _filter_sites(
                sites, station_filter, store_data, active_lines_store
            )
            if filtered is None:
                return "", html.Span("⚠ No data loaded", className="adv-train-warn")
            _CTRL.set_sites(filtered)
            _CTRL.train_dim_model(
                n_atoms=int(n_atoms or 6),
                n_iter=int(n_iter or 40),
            )
            return "", html.Span(
                f"✓ Trained — {int(n_atoms or 6)} atoms, {int(n_iter or 40)} iter",
                className="adv-train-ok",
            )
        except Exception as exc:
            return "", html.Span(f"✗ {exc}", className="adv-train-warn")

    # ── 10. Generate → render to active tab's figure ─────────────────────────
    @app.callback(
        Output(IDs.IMG_ADV_STRIKE, "src"),
        Output(IDs.IMG_ADV_PT, "src"),
        Output(IDs.IMG_ADV_INDUCTION, "src"),
        Output(IDs.IMG_ADV_IMPEDANCE, "src"),
        Output(IDs.IMG_ADV_DEPTH, "src"),
        Output(IDs.IMG_ADV_SURVEY, "src"),
        Output(IDs.ADV_SPINNER, "children"),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Input(IDs.BTN_ADV_RUN, "n_clicks"),
        State(IDs.ADV_ACTIVE_TAB, "data"),
        State(IDs.ADV_PLOT, "value"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.ADV_STN_FILTER, "value"),
        State("adv-period-idx", "value"),
        State("adv-figsize", "value"),
        State(IDs.ADV_PT_STATION, "value"),
        State(IDs.ADV_PT_PER_LINE, "value"),
        prevent_initial_call=True,
    )
    def run_advanced(
        n_clicks,
        active_tab,
        fn_name,
        session_id,
        store_data,
        active_lines_store,
        station_filter,
        period_idx,
        figsize,
        pt_station,
        pt_per_line,
    ):
        if not n_clicks or not fn_name:
            raise PreventUpdate

        try:
            sites = cache_get(session_id)
            filtered = _filter_sites(
                sites, station_filter, store_data, active_lines_store
            )

            apply_web_dark_theme()
            w, h = _FIGSIZE_MAP.get(figsize or "10x6", (10, 6))

            # ── Special-case the two PT functions needing extra required
            # kwargs (station= / profiles=) that the generic catalogue
            # dispatch (AdvancedController.draw) cannot supply.
            if fn_name == "plot_phase_tensor_strip":
                if not pt_station:
                    return (
                        *([no_update] * 6),
                        "",
                        True,
                        "Pick a station for the ellipse strip.",
                    )
                from pycsamt.emtools.tensor import (
                    plot_phase_tensor_strip,
                )

                fig = plt.figure(figsize=(w, h))
                ax = fig.add_subplot(111)
                plot_phase_tensor_strip(
                    sites,
                    station=pt_station,
                    ax=ax,
                    verbose=0,
                )
                render_fig = fig
            elif fn_name == "plot_phase_tensor_strip_grid":
                profiles = _pt_grid_profiles(
                    store_data,
                    active_lines_store,
                    int(pt_per_line or 4),
                )
                if not profiles:
                    return (
                        *([no_update] * 6),
                        "",
                        True,
                        "No survey lines detected for the strip grid.",
                    )
                from pycsamt.emtools.tensor import (
                    plot_phase_tensor_strip_grid,
                )

                render_fig = plot_phase_tensor_strip_grid(
                    sites,
                    profiles=profiles,
                    verbose=0,
                )
            else:
                _CTRL.set_sites(filtered)
                kwargs = {}
                if fn_name in _PERIOD_FNS and period_idx is not None:
                    kwargs["ip"] = max(0, int(period_idx))
                has_ax = _HAS_AX.get(fn_name, True)
                fig = plt.figure(figsize=(w, h))
                result = _CTRL.draw(fn_name, has_ax=has_ax, fig=fig, **kwargs)
                render_fig = result if result is not None else fig

            src = fig_to_src(render_fig)
            plt.close("all")

            # Write src only to the active tab's image; leave others unchanged
            outputs = [no_update] * 6
            tab_i = _ADV_TAB_ORDER.get(active_tab, 0)
            outputs[tab_i] = src

            return (*outputs, "", False, "")

        except Exception as exc:
            plt.close("all")
            return (*([no_update] * 6), "", True, f"{fn_name}: {exc}")
