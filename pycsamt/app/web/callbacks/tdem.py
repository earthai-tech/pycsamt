# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the TDEM analysis page."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dash import Input, Output, State, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.tdem_controller import (
    TDEMController, TDEM_GROUPS,
)
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import apply_web_dark_theme, empty_src, fig_to_src

_CTRL = TDEMController()


def register_tdem(app) -> None:

    @app.callback(
        Output(IDs.TDEM_INFO,       "children"),
        Output(IDs.TOAST_ERROR,     "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,      "children", allow_duplicate=True),
        Input(IDs.BTN_TDEM_LOAD,    "n_clicks"),
        State(IDs.TDEM_FOLDER,      "value"),
        prevent_initial_call=True,
    )
    def load_tdem(n_clicks, folder):
        if not n_clicks or not folder:
            raise PreventUpdate
        try:
            ok = _CTRL.load_folder(folder.strip())
            if ok:
                return _CTRL.summary, False, ""
            return "No TDEM files found in that folder.", True, "TDEM: folder empty"
        except Exception as exc:
            return f"Error: {exc}", True, f"TDEM load: {exc}"

    @app.callback(
        Output(IDs.TDEM_PLOT,   "options"),
        Output(IDs.TDEM_PLOT,   "value"),
        Input(IDs.TDEM_CAT,     "value"),
    )
    def sync_tdem_plots(cat):
        if not cat:
            raise PreventUpdate
        for group, plots in TDEM_GROUPS:
            if group == cat:
                opts = [{"label": lbl, "value": cls} for lbl, cls, *_ in plots]
                return opts, opts[0]["value"] if opts else None
        return no_update, no_update

    @app.callback(
        Output(IDs.IMG_TDEM,       "src"),
        Output(IDs.TDEM_SPINNER,   "children"),
        Output(IDs.TOAST_ERROR,    "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,     "children", allow_duplicate=True),
        Input(IDs.BTN_TDEM_RUN,    "n_clicks"),
        State(IDs.TDEM_PLOT,       "value"),
        prevent_initial_call=True,
    )
    def run_tdem(n_clicks, fn_name):
        if not n_clicks or not fn_name:
            raise PreventUpdate
        if not _CTRL.is_loaded():
            return empty_src(dark=True), "", True, "TDEM: load a folder first"
        try:
            apply_web_dark_theme()
            fig = plt.figure(figsize=(10, 5))
            _CTRL.draw(fn_name, has_ax=True, fig=fig)
            src = fig_to_src(fig)
            return src, "", False, ""
        except Exception as exc:
            return empty_src(dark=True), "", True, f"{fn_name}: {exc}"
