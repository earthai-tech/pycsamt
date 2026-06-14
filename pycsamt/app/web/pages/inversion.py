# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MT Inversion page — Occam2D / ModEM / AI solvers."""
from __future__ import annotations

from dash import html
import dash_bootstrap_components as dbc

from pycsamt.app.web.layout import IDs, _icon
from pycsamt.app.web.utils import empty_src

PAGE_ID = "inversion"

_SOLVERS = {
    "Traditional": [
        {"label": "Occam2D  (smooth-model 2-D)",  "value": "occam2d"},
        {"label": "ModEM    (3-D anisotropic)",    "value": "modem"},
        {"label": "MARE2DEM (marine / land 2-D)",  "value": "mare2dem"},
    ],
    "AI / Neural": [
        {"label": "EMInverter1D  (LSTM 1-D)",      "value": "em1d"},
        {"label": "EMInverter2D  (U-Net 2-D)",     "value": "em2d"},
        {"label": "GCNInverter3D (GCN 3-D)",       "value": "gcn3d"},
    ],
}


def layout() -> html.Div:
    solver_radios = []
    for group_label, opts in _SOLVERS.items():
        solver_radios.append(html.Div(group_label, className="ctrl-label mt-2"))
        for opt in opts:
            solver_radios.append(
                html.Div(dbc.RadioItems(
                    options=[opt], value=None, id={"type": "solver-radio", "value": opt["value"]},
                    inline=False,
                ), className="solver-radio")
            )

    controls = html.Div([
        html.Div([
            html.Div("Inversion Dimension", className="ctrl-label"),
            dbc.RadioItems(
                id=IDs.INV_DIM,
                options=[{"label": "1-D (sounding)", "value": "1d"},
                         {"label": "2-D (profile)",  "value": "2d"},
                         {"label": "3-D (volume)",   "value": "3d"}],
                value="2d", inline=True,
                className="text-muted small",
            ),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Solver", className="ctrl-label"),
            dbc.Select(
                id=IDs.INV_SOLVER,
                options=[
                    {"label": "── Traditional ──", "value": "", "disabled": True},
                    {"label": "Occam2D  (smooth 2-D)",   "value": "occam2d"},
                    {"label": "ModEM    (3-D aniso.)",   "value": "modem"},
                    {"label": "MARE2DEM (marine 2-D)",   "value": "mare2dem"},
                    {"label": "── AI / Neural ──",       "value": "", "disabled": True},
                    {"label": "EMInverter1D  (LSTM)",    "value": "em1d"},
                    {"label": "EMInverter2D  (U-Net)",   "value": "em2d"},
                    {"label": "GCNInverter3D (GCN)",     "value": "gcn3d"},
                ],
                value="occam2d", size="sm",
            ),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Data Directory", className="ctrl-label"),
            dbc.Input(id=IDs.INV_DATA_DIR,
                      placeholder="/path/to/EDI  (leave blank for loaded data)",
                      size="sm"),
        ], className="ctrl-card"),
        html.Div([
            dbc.Button([html.I(className="bi bi-play-fill me-1"), "Launch Inversion"],
                       id=IDs.BTN_INV_RUN, color="danger", size="sm",
                       className="w-100 mb-2"),
            dbc.Spinner(html.Div(id=IDs.INV_SPINNER), size="sm", color="danger"),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Info", className="ctrl-label"),
            html.Div([
                html.P("Occam2D/ModEM require the respective binaries to be "
                       "installed and on PATH.", className="text-muted small"),
                html.P("AI solvers run entirely in Python (no external binary).",
                       className="text-muted small"),
            ]),
        ], className="ctrl-card"),
    ], className="analysis-controls")

    output = html.Div([
        dbc.Tabs([
            dbc.Tab(
                html.Pre("", id=IDs.INV_LOG, className="log-area",
                         style={"height": "500px", "margin": "8px"}),
                label="Solver Log", tab_id="inv-tab-log",
            ),
            dbc.Tab(
                html.Div(
                    html.Img(id=IDs.IMG_INV, src=empty_src(dark=True),
                             style={"width": "100%", "maxHeight": "100%",
                                    "objectFit": "contain"}),
                    className="fig-area",
                    style={"minHeight": "480px", "margin": "8px"},
                ),
                label="Result Section", tab_id="inv-tab-result",
            ),
        ], active_tab="inv-tab-log"),
    ], className="analysis-output")

    from pycsamt.app.web.layout import _command_bar
    return html.Div([
        _command_bar("inversion", "MT Inversion",
                     "1-D / 2-D / 3-D · Occam2D · ModEM · AI neural solvers"),
        html.Div([controls, output], className="analysis-layout", style={"flex": "1"}),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # callbacks/inversion.py (stub — launchers are external binaries)
