# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/qc.py — QC page v2 callbacks.

Wiring
------
1. set_active_group  — pill click → QC_GROUP_SEL store + button classes
2. tab_switch        — clientside: Plot|Overview tab toggle
3. sync_qc_plots     — group change → RadioItems options+value
4. populate_qc_lines — STORE_DATA → QC_LINE_SEL options
5. populate_qc_stns  — sel_lines  → QC_STN_SEL options
6. run_qc            — main render: BTN_QC_RUN | auto (plot/group/data change)
7. run_qc_overview   — Overview quicklook: auto on nav to QC page or data load
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from dash import ALL, Input, Output, State, no_update
from dash import callback_context as ctx

from pycsamt.app.desktop.controllers.qc_controller import (
    ALL_GROUPS,
    QCController,
)
from pycsamt.app.web.cache import HAS_BG_MANAGER, cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    apply_web_light_theme,
    empty_src,
    fig_to_src,
    filter_sites_by_lines,
    no_active_lines_src,
)

# ── Catalogue helpers ─────────────────────────────────────────────────────────

_GROUP_MAP: dict = {gname: plots for gname, plots in ALL_GROUPS}
_GROUP_NAMES: list = [g for g, _ in ALL_GROUPS]

_FIGSIZE_MAP: dict = {
    "compact": (9, 5),
    "wide": (14, 5),
    "tall": (9, 8),
    "pub": (12, 7),
}


def _plot_options(group_name: str) -> list[dict]:
    plots = _GROUP_MAP.get(group_name, [])
    return [{"label": lbl, "value": fn} for lbl, fn, _ in plots]


def _has_ax(fn_name: str) -> bool:
    for _g, plots in ALL_GROUPS:
        for _lbl, fn, ha in plots:
            if fn == fn_name:
                return ha
    return True


# ── Registration ──────────────────────────────────────────────────────────────


def register_qc(app) -> None:

    # ── 1. Group pill click → QC_GROUP_SEL store + active classes ────────────
    @app.callback(
        Output(IDs.QC_GROUP_SEL, "data"),
        Output({"type": "qc-grp-btn", "index": ALL}, "className"),
        Input({"type": "qc-grp-btn", "index": ALL}, "n_clicks"),
        State(IDs.QC_GROUP_SEL, "data"),
        prevent_initial_call=True,
    )
    def set_active_group(n_clicks_list, current_group):
        triggered = ctx.triggered_id
        if not triggered or not isinstance(triggered, dict):
            return no_update, no_update
        new_group = triggered.get("index", current_group)
        classes = [
            "qc-grp-btn active" if g == new_group else "qc-grp-btn"
            for g in _GROUP_NAMES
        ]
        return new_group, classes

    # ── 2. Clientside tab switcher ────────────────────────────────────────────
    app.clientside_callback(
        """
        function(n_plot, n_ovl, cur_tab) {
            const tid = dash_clientside.callback_context.triggered_id;
            let tab = cur_tab || 'plot';
            if (tid === 'qc-tab-btn-plot')     tab = 'plot';
            if (tid === 'qc-tab-btn-overview') tab = 'overview';
            const isPlot = (tab === 'plot');
            return [
                tab,
                isPlot  ? {display: 'flex'} : {display: 'none'},
                !isPlot ? {display: 'flex'} : {display: 'none'},
                isPlot  ? 'qc-tab-btn active' : 'qc-tab-btn',
                !isPlot ? 'qc-tab-btn active' : 'qc-tab-btn',
            ];
        }
        """,
        Output(IDs.QC_ACTIVE_TAB, "data"),
        Output("qc-panel-plot", "style"),
        Output("qc-panel-overview", "style"),
        Output("qc-tab-btn-plot", "className"),
        Output("qc-tab-btn-overview", "className"),
        Input("qc-tab-btn-plot", "n_clicks"),
        Input("qc-tab-btn-overview", "n_clicks"),
        State(IDs.QC_ACTIVE_TAB, "data"),
        prevent_initial_call=True,
    )

    # ── 3. Sync plot RadioItems when group changes ────────────────────────────
    @app.callback(
        Output(IDs.QC_PLOT_SELECT, "options"),
        Output(IDs.QC_PLOT_SELECT, "value"),
        Input(IDs.QC_GROUP_SEL, "data"),
        prevent_initial_call=True,
    )
    def sync_qc_plots(group_name):
        opts = _plot_options(group_name or "Overview")
        val = opts[0]["value"] if opts else None
        return opts, val

    # ── 4. Populate Lines dropdown ────────────────────────────────────────────
    @app.callback(
        Output(IDs.QC_LINE_SEL, "options"),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def populate_qc_lines(store_data):
        records = (store_data or {}).get("station_records", [])
        lines = sorted({r.get("line", "") for r in records if r.get("line")})
        return [{"label": ln, "value": ln} for ln in lines]

    # ── 5. Populate Stations dropdown from selected lines ─────────────────────
    @app.callback(
        Output(IDs.QC_STN_SEL, "options"),
        Input(IDs.QC_LINE_SEL, "value"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def populate_qc_stns(sel_lines, store_data):
        records = (store_data or {}).get("station_records", [])
        if sel_lines:
            records = [r for r in records if r.get("line") in sel_lines]
        stns = sorted({r.get("station", "") for r in records if r.get("station")})
        return [{"label": s, "value": s} for s in stns]

    # ── 6. Main QC render ─────────────────────────────────────────────────────
    _bg = HAS_BG_MANAGER
    _cb_kw: dict = dict(prevent_initial_call=True)
    if _bg:
        _cb_kw["background"] = True
        _cb_kw["running"] = [
            (Output(IDs.BTN_QC_RUN, "disabled"), True, False),
            (Output(IDs.QC_SPINNER, "children"), " ", ""),
        ]

    @app.callback(
        Output(IDs.IMG_QC_PLOT, "src"),
        Output(IDs.QC_ACTIVE_TAB, "data", allow_duplicate=True),
        Output(IDs.QC_FEEDBACK, "children"),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        # ── Triggers ──────────────────────────────────────────────────────────
        Input(IDs.BTN_QC_RUN, "n_clicks"),
        Input(IDs.QC_PLOT_SELECT, "value"),
        Input(IDs.NAV_SECTION, "data"),
        Input(IDs.STORE_DATA, "data"),
        # ── Data filters ───────────────────────────────────────────────────────
        State(IDs.QC_LINE_SEL, "value"),
        State(IDs.QC_STN_SEL, "value"),
        State(IDs.QC_COMP_SEL, "value"),
        # ── Contextual params ──────────────────────────────────────────────────
        State(IDs.QC_PARAM_METRIC, "value"),
        State(IDs.QC_PARAM_METHOD, "value"),
        State(IDs.QC_PARAM_THRESH, "value"),
        State(IDs.QC_PARAM_NMAX, "value"),
        State(IDs.QC_PARAM_MAINS, "value"),
        # ── Fig size + auto toggle ─────────────────────────────────────────────
        State(IDs.QC_FIGSIZE, "value"),
        State(IDs.QC_AUTO_TOGGLE, "value"),
        # ── Session + theme ───────────────────────────────────────────────────
        State(IDs.STORE_THEME, "data"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        **_cb_kw,
    )
    def run_qc(
        n_run,
        fn_name,
        nav_section,
        store_data,
        sel_lines,
        sel_stns,
        comp,
        metric,
        method,
        thresh,
        n_max,
        mains_hz,
        figsize_key,
        auto_toggle,
        theme,
        session_id,
        active_lines_store,
    ):

        triggered = ctx.triggered_id
        manual_run = triggered == IDs.BTN_QC_RUN
        on_qc_page = nav_section == "qc"
        auto_on = bool(auto_toggle and "auto" in auto_toggle)
        dark = (theme or "dark") == "dark"
        _nu5 = (no_update,) * 5

        # Only render when on QC page (or manual override)
        if not on_qc_page and not manual_run:
            return _nu5

        # Auto guard: skip auto-render when Auto is toggled off
        if not manual_run and not auto_on:
            return _nu5

        if not fn_name or not session_id:
            if manual_run:
                return (
                    no_update,
                    no_update,
                    "Select a plot first.",
                    True,
                    "Select a plot first.",
                )
            return _nu5

        sites = cache_get(session_id)
        if sites is None:
            if manual_run:
                return (
                    empty_src(dark=dark),
                    no_update,
                    "No data loaded.",
                    True,
                    "Load survey data first.",
                )
            return _nu5

        # ── Filter by selected lines ──────────────────────────────────────────
        records = (store_data or {}).get("station_records", [])
        _als = active_lines_store or {}
        active_lines = sel_lines or _als.get("active", _als.get("all"))
        if active_lines:
            sites = filter_sites_by_lines(sites, records, active_lines)
            if sites is None:
                return (
                    no_active_lines_src(dark),
                    no_update,
                    "No matching lines.",
                    False,
                    "",
                )

        # ── Filter by selected stations ───────────────────────────────────────
        if sel_stns:
            _ctrl_tmp = QCController()
            _ctrl_tmp.set_sites(sites)
            _ctrl_tmp.filter_sites(sel_stns)
            sites = _ctrl_tmp._sites

        # ── Render ────────────────────────────────────────────────────────────
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        try:
            import matplotlib.pyplot as plt

            ctrl = QCController()
            ctrl.set_sites(sites)
            ctrl.dark = dark

            figsize = _FIGSIZE_MAP.get(figsize_key or "compact", (9, 5))
            has_ax = _has_ax(fn_name)
            fig = plt.figure(figsize=figsize)

            draw_kwargs = dict(
                metric=metric,
                method=method,
                thresh=thresh,
                n_max=n_max,
                mains_hz=mains_hz,
                component=comp,
            )

            alt = ctrl.draw(fn_name, has_ax, fig, figsize=None, **draw_kwargs)
            target = alt if alt is not None else fig
            src = fig_to_src(target)
            plt.close("all")
            return src, "plot", "", False, ""

        except Exception as exc:
            msg = f"{fn_name}: {exc}"
            return (
                empty_src(dark=dark),
                no_update,
                msg[:90],
                True,
                f"QC error: {exc}",
            )

    # ── 7. Overview quicklook — auto on nav / data-load ───────────────────────
    @app.callback(
        Output(IDs.IMG_QC_OVERVIEW, "src"),
        Input(IDs.NAV_SECTION, "data"),
        Input(IDs.STORE_DATA, "data"),
        State(IDs.STORE_THEME, "data"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def run_qc_overview(nav_section, store_data, theme, session_id, active_lines_store):
        if nav_section != "qc":
            return no_update

        sites = cache_get(session_id)
        if sites is None:
            return no_update

        dark = (theme or "dark") == "dark"
        records = (store_data or {}).get("station_records", [])
        _als = active_lines_store or {}
        active_lines = _als.get("active", _als.get("all"))
        if active_lines:
            sites = filter_sites_by_lines(sites, records, active_lines)
            if sites is None:
                return no_update

        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        try:
            import matplotlib.pyplot as plt

            ctrl = QCController()
            ctrl.set_sites(sites)
            ctrl.dark = dark
            fig = plt.figure(figsize=(9, 5))
            alt = ctrl.draw("plot_qc_quicklook", has_ax=False, fig=fig)
            target = alt if alt is not None else fig
            src = fig_to_src(target)
            plt.close("all")
            return src
        except Exception:
            return no_update
