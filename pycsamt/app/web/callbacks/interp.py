# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the Geological Interpretation page."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dash import Input, Output, State, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.interp_controller import (
    InterpController, WORKFLOW_CATALOGUE, CATEGORIES,
)
from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import apply_web_dark_theme, empty_src, fig_to_src

_CTRL = InterpController()


def register_interp(app) -> None:

    @app.callback(
        Output(IDs.INTERP_PLOT,  "options"),
        Output(IDs.INTERP_PLOT,  "value"),
        Output(IDs.INTERP_DESC,  "children"),
        Input(IDs.INTERP_CAT,    "value"),
    )
    def sync_interp_plots(cat):
        if not cat:
            raise PreventUpdate
        plots = WORKFLOW_CATALOGUE.get(cat, [])
        opts = [{"label": lbl, "value": fn} for lbl, fn, _ in plots]
        desc = plots[0][2] if plots else ""
        return opts, opts[0]["value"] if opts else None, desc

    @app.callback(
        Output(IDs.INTERP_DESC, "children", allow_duplicate=True),
        Input(IDs.INTERP_PLOT,  "value"),
        State(IDs.INTERP_CAT,   "value"),
        prevent_initial_call=True,
    )
    def update_desc(fn_name, cat):
        if not fn_name or not cat:
            raise PreventUpdate
        for lbl, fn, desc in WORKFLOW_CATALOGUE.get(cat, []):
            if fn == fn_name:
                return desc
        return no_update

    @app.callback(
        Output(IDs.IMG_INTERP,      "src"),
        Output(IDs.INTERP_SPINNER,  "children"),
        Output(IDs.TOAST_ERROR,     "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,      "children", allow_duplicate=True),
        Input(IDs.BTN_INTERP_RUN,   "n_clicks"),
        State(IDs.INTERP_PLOT,      "value"),
        State(IDs.SESSION_ID,       "data"),
        prevent_initial_call=True,
    )
    def run_interp(n_clicks, fn_name, session_id):
        if not n_clicks or not fn_name:
            raise PreventUpdate
        try:
            sites = cache_get(session_id)
            _CTRL.set_sites(sites)
            apply_web_dark_theme()
            fig = plt.figure(figsize=(11, 5))
            method = getattr(_CTRL, fn_name, None)
            if method is None:
                return empty_src(dark=True), "", True, f"Unknown method: {fn_name}"
            ax = fig.add_subplot(111)
            method(ax)
            src = fig_to_src(fig)
            return src, "", False, ""
        except Exception as exc:
            return empty_src(dark=True), "", True, f"{fn_name}: {exc}"
