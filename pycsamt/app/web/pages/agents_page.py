# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""AI Agents page — full-page runner with card-grid agent browser."""
from __future__ import annotations

import json

from dash import ALL, Input, Output, dcc, html
import dash_bootstrap_components as dbc

from pycsamt.app.desktop.agent_registry import (
    agents_by_category, agent_names, default_params,
)
from pycsamt.app.web.layout import IDs, _command_bar
from pycsamt.app.web.utils import empty_src

PAGE_ID = "agents"

_BY_CAT    = agents_by_category()
_ALL_NAMES = agent_names()

_CAT_OPTS = [{"label": "All categories", "value": ""}] + [
    {"label": cat, "value": cat} for cat in _BY_CAT
]

# Map agent name → category for card labels
_NAME_TO_CAT: dict[str, str] = {
    name: cat
    for cat, names in _BY_CAT.items()
    for name in names
}

# Bootstrap icon per category (falls back to robot)
_CAT_ICON: dict[str, str] = {
    "Processing":   "bi-sliders",
    "Inversion":    "bi-layers",
    "QC":           "bi-shield-check",
    "Reporting":    "bi-file-earmark-text",
    "LLM":          "bi-chat-dots",
    "Interpretation": "bi-lightbulb",
}
_DEFAULT_ICON = "bi-robot"


def _params_json(name: str) -> str:
    params = default_params(name)
    return json.dumps(params, indent=2) if params else "{}"


def _agent_grid(names: list[str]) -> html.Div:
    """Card grid browser for all registered agents."""
    cards = []
    for i, name in enumerate(names):
        cat  = _NAME_TO_CAT.get(name, "")
        icon = _CAT_ICON.get(cat, _DEFAULT_ICON)
        cards.append(
            html.Button(
                [
                    html.I(className=f"bi {icon} agent-card-icon"),
                    html.Span(name, className="agent-card-name"),
                    html.Span(cat or "General", className="agent-card-cat"),
                ],
                id={"type": "agent-card", "index": name},
                className="agent-card" + (" selected" if i == 0 else ""),
                n_clicks=0,
                title=name,
            )
        )
    return html.Div(cards, className="agent-grid")


def layout() -> html.Div:
    total = len(_ALL_NAMES)
    cats  = len(_BY_CAT)

    controls = html.Div([
        # Category filter
        html.Div([
            html.Div([html.I(className="bi bi-funnel me-1"), "Category"],
                     className="ctrl-label"),
            dbc.Select(
                id=IDs.AGENTS_CAT,
                options=_CAT_OPTS,
                value="",
                size="sm",
            ),
        ], className="ctrl-card"),

        # Agent card grid
        html.Div([
            html.Div([html.I(className="bi bi-grid me-1"), "Select Agent"],
                     className="ctrl-label"),
            _agent_grid(_ALL_NAMES),
            # Hidden select keeps existing callbacks working unchanged
            dbc.Select(
                id=IDs.AGENTS_NAME,
                options=[{"label": n, "value": n} for n in _ALL_NAMES],
                value=_ALL_NAMES[0] if _ALL_NAMES else None,
                size="sm",
                style={"display": "none"},
            ),
        ], className="ctrl-card"),

        # Parameters
        html.Div([
            html.Div([html.I(className="bi bi-braces me-1"), "Parameters (JSON)"],
                     className="ctrl-label"),
            dbc.Textarea(
                id=IDs.AGENTS_PARAMS,
                value=_params_json(_ALL_NAMES[0]) if _ALL_NAMES else "{}",
                rows=7,
                style={
                    "fontFamily": "'JetBrains Mono','Fira Code','Consolas',monospace",
                    "fontSize": "11px",
                    "backgroundColor": "#181825",
                    "color": "#cdd6f4",
                    "border": "1px solid #313244",
                    "borderRadius": "4px",
                    "resize": "vertical",
                },
            ),
        ], className="ctrl-card"),

        # Run button
        html.Div([
            dbc.Button(
                [html.I(className="bi bi-play-fill me-1"), "Run Agent"],
                id=IDs.BTN_AGENTS_RUN, color="warning", size="sm",
                className="w-100 mb-2",
            ),
            dbc.Spinner(html.Div(id=IDs.AGENTS_SPINNER),
                        size="sm", color="warning"),
        ], className="ctrl-card"),

        dcc.Store(id=IDs.AGENTS_STORE, storage_type="memory"),
    ], className="analysis-controls")

    output = html.Div([
        dbc.Tabs([
            dbc.Tab(
                html.Pre(
                    id=IDs.AGENTS_OUT,
                    children="Agent output will appear here…",
                    className="log-area",
                    style={"height": "520px", "margin": "8px",
                           "whiteSpace": "pre-wrap", "wordBreak": "break-word"},
                ),
                label="Output", tab_id="agent-tab-out",
            ),
            dbc.Tab(
                html.Div(
                    html.Img(
                        id=IDs.IMG_AGENTS,
                        src=empty_src(dark=True),
                        style={"width": "100%", "maxHeight": "100%",
                               "objectFit": "contain"},
                    ),
                    className="fig-area",
                    style={"minHeight": "500px", "margin": "8px"},
                ),
                label="Figure", tab_id="agent-tab-fig",
            ),
        ], active_tab="agent-tab-out"),
    ], className="analysis-output")

    return html.Div([
        _command_bar(
            "agents", "AI Agents",
            f"{total} agents · {cats} categories · LLM · Processing · Inversion · Reporting",
        ),
        html.Div([controls, output], className="analysis-layout", style={"flex": "1"}),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    """Wire agent card clicks → update hidden select + highlight selected card."""
    from dash import callback, ctx, no_update

    @callback(
        Output(IDs.AGENTS_NAME, "value"),
        Input({"type": "agent-card", "index": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def _card_to_select(n_clicks_list):
        if not ctx.triggered_id:
            return no_update
        triggered = ctx.triggered_id
        if isinstance(triggered, dict) and "index" in triggered:
            return triggered["index"]
        return no_update

    @callback(
        Output({"type": "agent-card", "index": ALL}, "className"),
        Input(IDs.AGENTS_NAME, "value"),
        prevent_initial_call=True,
    )
    def _highlight_card(selected_name):
        if not selected_name:
            return no_update
        return [
            "agent-card selected" if name == selected_name else "agent-card"
            for name in _ALL_NAMES
        ]
