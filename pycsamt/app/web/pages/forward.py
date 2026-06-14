# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Forward modelling page (ForwardController)."""
from __future__ import annotations

from dash import dash_table, dcc, html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.controllers.forward_controller import (
    ForwardController, GEOLOGY_PRESET_NAMES,
)
from pycsamt.app.web.layout import IDs, _icon
from pycsamt.app.web.utils import empty_src

PAGE_ID = "forward"

_DIM_OPTS    = [{"label": "1-D Layered Earth", "value": "1d"},
                {"label": "2-D Section",        "value": "2d"},
                {"label": "3-D Volume",          "value": "3d"}]
_PRESET_OPTS = [{"label": p.replace("_", " ").title(), "value": p}
                for p in GEOLOGY_PRESET_NAMES]


def _layer_table():
    return dash_table.DataTable(
        id=IDs.FWD_LAYER_TABLE,
        columns=[
            {"name": "Layer", "id": "layer", "editable": False},
            {"name": "Resistivity (Ω·m)", "id": "resistivity", "editable": True,
             "type": "numeric"},
            {"name": "Thickness (m)", "id": "thickness", "editable": True,
             "type": "numeric"},
        ],
        data=[{"layer": i + 1, "resistivity": 100, "thickness": 50}
              for i in range(5)],
        editable=True,
        row_deletable=True,
        style_table={"overflowY": "auto", "maxHeight": "240px"},
        style_header={"backgroundColor": "#11111b", "color": "#6c7086",
                      "fontSize": "11px", "border": "none",
                      "textTransform": "uppercase", "fontWeight": "600"},
        style_cell={"backgroundColor": "#181825", "color": "#cdd6f4",
                    "border": "1px solid #313244", "fontSize": "11px",
                    "padding": "4px 8px"},
        style_data_conditional=[{"if": {"state": "active"},
                                 "backgroundColor": "#313244",
                                 "border": "1px solid #89b4fa"}],
    )


def layout() -> html.Div:
    controls = html.Div([
        html.Div([
            html.Div("Dimension", className="ctrl-label"),
            dbc.RadioItems(id=IDs.FWD_DIM, options=_DIM_OPTS,
                           value="1d", inline=False,
                           className="solver-radio"),
        ], className="ctrl-card"),
        html.Div([
            html.Div("Geological Preset", className="ctrl-label"),
            dbc.Select(id=IDs.FWD_PRESET, options=_PRESET_OPTS,
                       value=_PRESET_OPTS[0]["value"], size="sm"),
            dbc.Button([html.I(className="bi bi-magic me-1"), "Load Preset"],
                       id=IDs.BTN_FWD_PRESET, color="secondary",
                       size="sm", className="w-100 mt-2"),
        ], className="ctrl-card"),
        html.Div([
            dbc.Button([html.I(className="bi bi-play-fill me-1"), "Run Forward"],
                       id=IDs.BTN_FWD_RUN, color="primary", size="sm",
                       className="w-100 mb-2"),
            dbc.Row([
                dbc.Col(dbc.Input(id=IDs.FWD_MODEL_NAME,
                                  placeholder="Model name…", size="sm")),
                dbc.Col(dbc.Button([_icon("save", size=12, cls=""), " Save"],
                                   id=IDs.BTN_FWD_SAVE, color="secondary",
                                   size="sm"), width="auto"),
            ], className="g-1"),
            html.Div(
                dbc.Spinner(html.Div(id=IDs.FWD_SPINNER), size="sm", color="primary"),
                className="mt-2",
            ),
            html.Div(id=IDs.FWD_FEEDBACK, className="mt-2 text-muted small"),
        ], className="ctrl-card"),
    ], className="analysis-controls")

    output = html.Div([
        html.Div("Layer Model Editor", className="ctrl-label",
                 style={"padding": "8px 8px 4px"}),
        _layer_table(),
        dbc.Button([html.I(className="bi bi-plus-lg me-1"), "Add Layer"],
                   id=IDs.BTN_FWD_ADD_LAYER, color="secondary",
                   size="sm", className="mt-2 ms-1 mb-3"),
        html.Div([
            html.Img(id=IDs.IMG_FWD, src=empty_src(dark=True),
                     style={"width": "100%", "maxHeight": "100%",
                            "objectFit": "contain"}),
        ], className="fig-area", style={"flex": "1", "minHeight": "280px"}),
    ], className="analysis-output", style={"flexDirection": "column"})

    from pycsamt.app.web.layout import _command_bar
    return html.Div([
        _command_bar("forward", "Forward Modelling",
                     f"1-D / 2-D / 3-D · {len(GEOLOGY_PRESET_NAMES)} geological presets"),
        html.Div([controls, output], className="analysis-layout", style={"flex": "1"}),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # handled by callbacks/forward.py
