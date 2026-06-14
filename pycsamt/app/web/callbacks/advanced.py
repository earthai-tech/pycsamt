# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the Advanced Plots page."""
from __future__ import annotations

import traceback

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dash import Input, Output, State, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.advanced_controller import (
    AdvancedController, ADVANCED_GROUPS,
)
from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import apply_web_dark_theme, empty_src, fig_to_src

_CTRL = AdvancedController()

# Filtered groups (skip Topography / Conversion stubs with empty plot lists)
_GROUPS = [(g, plots) for g, plots in ADVANCED_GROUPS if plots]


def register_advanced(app) -> None:

    @app.callback(
        Output(IDs.ADV_PLOT, "options"),
        Output(IDs.ADV_PLOT, "value"),
        Input(IDs.ADV_GROUP, "value"),
    )
    def sync_plots(group):
        if not group:
            raise PreventUpdate
        for g, plots in _GROUPS:
            if g == group:
                opts = [{"label": lbl, "value": fn} for lbl, fn, *_ in plots]
                return opts, opts[0]["value"] if opts else None
        return no_update, no_update

    @app.callback(
        Output(IDs.IMG_ADV,       "src"),
        Output(IDs.ADV_SPINNER,   "children"),
        Output(IDs.TOAST_ERROR,   "is_open",   allow_duplicate=True),
        Output(IDs.TOAST_BODY,    "children",  allow_duplicate=True),
        Input(IDs.BTN_ADV_RUN,    "n_clicks"),
        State(IDs.ADV_PLOT,       "value"),
        State(IDs.SESSION_ID,     "data"),
        prevent_initial_call=True,
    )
    def run_advanced(n_clicks, fn_name, session_id):
        if not n_clicks or not fn_name:
            raise PreventUpdate
        try:
            sites = cache_get(session_id)
            _CTRL.set_sites(sites)
            apply_web_dark_theme()
            fig = plt.figure(figsize=(10, 5))
            _CTRL.draw(fn_name, has_ax=True, fig=fig)
            src = fig_to_src(fig)
            return src, "", False, ""
        except Exception as exc:
            return empty_src(dark=True), "", True, f"{fn_name}: {exc}"
