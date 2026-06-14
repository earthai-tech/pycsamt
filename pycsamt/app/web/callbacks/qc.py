# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/qc.py — QC tab: group/plot dropdowns and background render.

Two callbacks:
  1. sync_qc_plots   — updates QC_PLOT_SELECT options when QC_GROUP_SELECT changes
  2. run_qc          — background callback, renders one QC plot into IMG_QC
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

from dash import Input, Output, State, no_update

from pycsamt.app.desktop.controllers.qc_controller import ALL_GROUPS, QCController
from pycsamt.app.web.cache import cache_get, has_diskcache
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    apply_web_light_theme,
    empty_src,
    fig_to_src,
)

# Build a quick lookup: group_name → [(label, fn_name, has_ax), ...]
_GROUP_MAP: dict = {gname: plots for gname, plots in ALL_GROUPS}


def _plot_options(group_name: str) -> list[dict]:
    plots = _GROUP_MAP.get(group_name, [])
    return [{"label": lbl, "value": fn} for lbl, fn, _ in plots]


def _has_ax(group_name: str, fn_name: str) -> bool:
    for lbl, fn, ha in _GROUP_MAP.get(group_name, []):
        if fn == fn_name:
            return ha
    return True


def register_qc(app) -> None:
    # ── 1. Sync plot dropdown when group changes ──────────────────────────────
    @app.callback(
        Output(IDs.QC_PLOT_SELECT, "options"),
        Output(IDs.QC_PLOT_SELECT, "value"),
        Input(IDs.QC_GROUP_SELECT, "value"),
        prevent_initial_call=True,
    )
    def sync_qc_plots(group_name):
        opts = _plot_options(group_name or "Overview")
        val  = opts[0]["value"] if opts else None
        return opts, val

    # ── 2. QC render (background when diskcache available, else sync) ─────────
    _bg = has_diskcache()
    _cb_kwargs: dict = dict(prevent_initial_call=True)
    if _bg:
        _cb_kwargs["background"] = True
        _cb_kwargs["running"] = [
            (Output(IDs.BTN_QC_RUN, "disabled"), True,  False),
            (Output(IDs.QC_SPINNER, "children"), " ",   ""),
        ]

    @app.callback(
        Output(IDs.IMG_QC,       "src"),
        Output(IDs.TOAST_ERROR,  "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,   "children", allow_duplicate=True),
        Input(IDs.BTN_QC_RUN,   "n_clicks"),
        State(IDs.QC_GROUP_SELECT, "value"),
        State(IDs.QC_PLOT_SELECT,  "value"),
        State(IDs.STORE_THEME,     "data"),
        State(IDs.SESSION_ID,      "data"),
        **_cb_kwargs,
    )
    def run_qc(_n, group_name, fn_name, theme, session_id):
        if not fn_name:
            return no_update, False, ""

        if not session_id:
            msg = "Session not initialised — refresh the page."
            return empty_src(dark=True), True, msg

        sites = cache_get(session_id)
        if sites is None:
            msg = "Load survey data before running a QC plot."
            return empty_src(dark=True), True, msg

        dark = (theme or "dark") == "dark"
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        try:
            import matplotlib.pyplot as plt

            ctrl = QCController()
            ctrl.set_sites(sites)
            ctrl.dark = dark

            has_ax = _has_ax(group_name or "Overview", fn_name)
            fig    = plt.figure(figsize=(9, 5))
            alt    = ctrl.draw(fn_name, has_ax, fig)

            target = alt if alt is not None else fig
            src    = fig_to_src(target)
            plt.close("all")
            return src, False, ""

        except Exception as exc:
            msg = f"QC plot failed ({fn_name}): {exc}"
            return empty_src(dark=True), True, msg
