# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the Forward Modelling page."""
from __future__ import annotations

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dash import Input, Output, State, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.forward_controller import (
    ForwardController,
)
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import apply_web_dark_theme, empty_src, fig_to_src

_CTRL = ForwardController()


def register_forward(app) -> None:

    @app.callback(
        Output(IDs.FWD_LAYER_TABLE, "data"),
        Output(IDs.FWD_FEEDBACK,    "children"),
        Input(IDs.BTN_FWD_PRESET,   "n_clicks"),
        State(IDs.FWD_PRESET,       "value"),
        prevent_initial_call=True,
    )
    def load_preset(n_clicks, preset_name):
        if not n_clicks or not preset_name:
            raise PreventUpdate
        try:
            record = ForwardController.build_preset_1d(preset_name)
            rho  = record["resistivity"]
            thk  = record["thickness"]
            rows = [{"layer": i + 1, "resistivity": float(r), "thickness": float(t)}
                    for i, (r, t) in enumerate(zip(rho, thk))]
            return rows, f"Loaded preset: {preset_name}"
        except Exception as exc:
            return no_update, f"Error: {exc}"

    @app.callback(
        Output(IDs.IMG_FWD,        "src"),
        Output(IDs.FWD_SPINNER,    "children"),
        Output(IDs.FWD_FEEDBACK,   "children",  allow_duplicate=True),
        Output(IDs.TOAST_ERROR,    "is_open",   allow_duplicate=True),
        Output(IDs.TOAST_BODY,     "children",  allow_duplicate=True),
        Input(IDs.BTN_FWD_RUN,     "n_clicks"),
        State(IDs.FWD_LAYER_TABLE, "data"),
        State(IDs.FWD_DIM,         "value"),
        prevent_initial_call=True,
    )
    def run_forward(n_clicks, table_data, dim):
        if not n_clicks or not table_data:
            raise PreventUpdate
        try:
            rho = np.array([float(row["resistivity"]) for row in table_data])
            thk = np.array([float(row["thickness"])   for row in table_data])

            apply_web_dark_theme()
            fig, axes = plt.subplots(1, 2, figsize=(11, 5))

            from pycsamt.forward import LayeredModel
            model = LayeredModel(resistivity=rho, thickness=thk)
            model.compute_response()

            axes[0].step(np.log10(rho), -np.cumsum(thk) / 1e3, where="post",
                         color="#89b4fa", lw=2)
            axes[0].set_xlabel("log₁₀(ρ) [Ω·m]")
            axes[0].set_ylabel("Depth [km]")
            axes[0].set_title("Layer model")

            if hasattr(model, "periods_") and hasattr(model, "app_res_"):
                axes[1].loglog(model.periods_, model.app_res_, color="#a6e3a1", lw=2)
                axes[1].set_xlabel("Period [s]")
                axes[1].set_ylabel("ρₐ [Ω·m]")
                axes[1].set_title("Apparent resistivity")
            else:
                axes[1].text(0.5, 0.5, "Forward response\nnot available",
                             ha="center", va="center", transform=axes[1].transAxes,
                             color="#a6adc8")

            fig.tight_layout()
            src = fig_to_src(fig)
            return src, "", f"Done  ({len(rho)} layers, dim={dim})", False, ""
        except Exception as exc:
            return empty_src(dark=True), "", f"Error: {exc}", True, str(exc)

    @app.callback(
        Output(IDs.FWD_FEEDBACK,    "children",  allow_duplicate=True),
        Input(IDs.BTN_FWD_SAVE,     "n_clicks"),
        State(IDs.FWD_MODEL_NAME,   "value"),
        State(IDs.FWD_LAYER_TABLE,  "data"),
        State(IDs.FWD_DIM,          "value"),
        prevent_initial_call=True,
    )
    def save_model(n_clicks, name, table_data, dim):
        if not n_clicks:
            raise PreventUpdate
        if not name:
            return "Enter a model name first."
        if not table_data:
            return "No layer data to save."
        try:
            rho = np.array([float(r["resistivity"]) for r in table_data])
            thk = np.array([float(r["thickness"])   for r in table_data])
            record = ForwardController.model_to_record(name, dim, rho, thk)
            _CTRL.save_model(name, record)
            return f"Saved '{name}' to library ({len(_CTRL.model_names)} models)."
        except Exception as exc:
            return f"Save error: {exc}"

    @app.callback(
        Output(IDs.FWD_LAYER_TABLE, "data",     allow_duplicate=True),
        Input(IDs.BTN_FWD_ADD_LAYER, "n_clicks"),
        State(IDs.FWD_LAYER_TABLE,  "data"),
        prevent_initial_call=True,
    )
    def add_layer(n_clicks, current_data):
        if not n_clicks:
            raise PreventUpdate
        data = list(current_data or [])
        data.append({"layer": len(data) + 1, "resistivity": 100, "thickness": 50})
        return data
