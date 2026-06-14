# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the Data Correction page."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.correction_controller import (
    CorrectionController, CATALOGUE, CATEGORIES,
)
from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import apply_web_dark_theme, empty_src, fig_to_src

_CTRL = CorrectionController()


def _sites_to_src(sites) -> str:
    """Render the Sites pseudo-section as a PNG src."""
    if sites is None:
        return empty_src(dark=True)
    try:
        apply_web_dark_theme()
        fig = plt.figure(figsize=(9, 4))
        ax = fig.add_subplot(111)
        try:
            sites.plot_pseudosection(ax=ax, component="res")
        except Exception:
            ax.text(0.5, 0.5, "No pseudo-section available",
                    ha="center", va="center", transform=ax.transAxes,
                    color="#a6adc8", fontsize=11)
        return fig_to_src(fig)
    except Exception:
        return empty_src(dark=True)


def register_correction(app) -> None:

    @app.callback(
        Output(IDs.CORR_METHOD, "options"),
        Output(IDs.CORR_METHOD, "value"),
        Input(IDs.CORR_CAT,     "value"),
    )
    def sync_methods(cat):
        if not cat:
            raise PreventUpdate
        methods = CATALOGUE.get(cat, {})
        opts = [{"label": k, "value": k} for k in methods]
        return opts, opts[0]["value"] if opts else None

    @app.callback(
        Output(IDs.IMG_CORR_LEFT,    "src"),
        Output(IDs.IMG_CORR_RIGHT,   "src"),
        Output(IDs.CORR_FEEDBACK,    "children"),
        Output(IDs.CORR_SPINNER,     "children"),
        Output(IDs.TOAST_ERROR,      "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,       "children", allow_duplicate=True),
        Input(IDs.BTN_CORR_PREVIEW,  "n_clicks"),
        State(IDs.CORR_CAT,          "value"),
        State(IDs.CORR_METHOD,       "value"),
        State(IDs.SESSION_ID,        "data"),
        prevent_initial_call=True,
    )
    def preview_correction(n_clicks, cat, method, session_id):
        if not n_clicks:
            raise PreventUpdate
        try:
            sites = cache_get(session_id)
            if sites is None:
                return (empty_src(dark=True), empty_src(dark=True),
                        "No data loaded.", "", True, "Correction: no data")
            _CTRL.set_raw_sites(sites)
            fn_name = CATALOGUE.get(cat, {}).get(method, {}).get("fn")
            kwargs  = {}
            before_src = _sites_to_src(_CTRL.current_sites)
            result = _CTRL.preview(fn_name, kwargs) if fn_name else None
            after_src  = _sites_to_src(result if result is not None else _CTRL.current_sites)
            return before_src, after_src, f"Preview: {method}", "", False, ""
        except Exception as exc:
            return (empty_src(dark=True), empty_src(dark=True),
                    f"Error: {exc}", "", True, str(exc))

    @app.callback(
        Output(IDs.CORR_STACK,      "children"),
        Output(IDs.CORR_FEEDBACK,   "children",  allow_duplicate=True),
        Output(IDs.TOAST_ERROR,     "is_open",   allow_duplicate=True),
        Output(IDs.TOAST_BODY,      "children",  allow_duplicate=True),
        Input(IDs.BTN_CORR_APPLY,   "n_clicks"),
        State(IDs.CORR_CAT,         "value"),
        State(IDs.CORR_METHOD,      "value"),
        State(IDs.SESSION_ID,       "data"),
        prevent_initial_call=True,
    )
    def apply_correction(n_clicks, cat, method, session_id):
        if not n_clicks:
            raise PreventUpdate
        try:
            sites = cache_get(session_id)
            if sites is None:
                return no_update, "No data loaded.", True, "Correction: no data"
            _CTRL.set_raw_sites(sites)
            fn_name = CATALOGUE.get(cat, {}).get(method, {}).get("fn")
            kwargs  = {}
            if fn_name:
                _CTRL.apply(fn_name, kwargs, label=method)
            items = [html.Li(s.label, className="text-muted small",
                             style={"fontSize": "11px"})
                     for s in _CTRL._stack]
            stack_el = html.Ul(items, style={"paddingLeft": "14px", "margin": "0"}) \
                if items else html.Span("None yet", className="text-muted small")
            return stack_el, f"Applied: {method}", False, ""
        except Exception as exc:
            return no_update, f"Error: {exc}", True, str(exc)

    @app.callback(
        Output(IDs.CORR_FEEDBACK,   "children",  allow_duplicate=True),
        Input(IDs.BTN_CORR_UNDO,    "n_clicks"),
        prevent_initial_call=True,
    )
    def undo_correction(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        _CTRL.undo_last()
        return f"Undone — {len(_CTRL._stack)} corrections remaining."

    @app.callback(
        Output(IDs.CORR_STACK,      "children",  allow_duplicate=True),
        Output(IDs.CORR_FEEDBACK,   "children",  allow_duplicate=True),
        Input(IDs.BTN_CORR_RESET,   "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_corrections(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        _CTRL.revert_all()
        return html.Span("None yet", className="text-muted small"), "Reset to raw."
