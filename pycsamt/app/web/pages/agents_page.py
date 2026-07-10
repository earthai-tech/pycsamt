# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""AI Agents page — 3-column professional layout with run history."""
from __future__ import annotations

import json
import uuid
from datetime import datetime

import dash_bootstrap_components as dbc
from dash import (
    ALL,
    Input,
    Output,
    State,
    ctx,
    dcc,
    html,
    no_update,
)

from pycsamt.app.desktop.agent_registry import (
    AGENT_REGISTRY,
    agent_names,
    agents_by_category,
    default_params,
    get_entry,
    processing_agents,
)
from pycsamt.app.web.layout import IDs, _command_bar
from pycsamt.app.web.utils import empty_src

# ── page-scoped IDs (not added to global IDs class) ───────────────────────
# Agent runner panel
_ID_SEARCH          = "agents-search"
_ID_LIST_WRAP       = "agents-list-wrap"
_ID_DETAIL_NAME     = "agents-detail-name"
_ID_HISTORY         = "agents-history-list"
_ID_SUMMARY         = "agents-result-summary"
_ID_DATA_BAR        = "agents-data-bar"
_ID_STATION_FILTER  = "agents-station-filter"
# View toggle
_ID_VIEW_RUNNER     = "agents-view-runner"
_ID_VIEW_CHAT       = "agents-view-chat"
_ID_TAB_RUNNER      = "agents-tab-btn-runner"
_ID_TAB_CHAT        = "agents-tab-btn-chat"
# Chat panel
_ID_CHAT_MESSAGES   = "agents-chat-messages"
_ID_CHAT_INPUT      = "agents-chat-input"
_ID_CHAT_SEND       = "agents-chat-send"
_ID_CHAT_STORE      = "agents-chat-store"
# Line/station filter
_ID_LINE_FILTER     = "agents-line-filter"
_ID_LINE_SEL_ALL    = "agents-line-sel-all"
_ID_LINE_DESEL_ALL  = "agents-line-desel-all"
_ID_STN_SEL_ALL     = "agents-stn-sel-all"
_ID_STN_DESEL_ALL   = "agents-stn-desel-all"

PAGE_ID = "agents"

# ── Module-level constants ─────────────────────────────────────────────────

_BY_CAT    = agents_by_category()          # {category: [names]}  (LLM only)
_ALL_NAMES = agent_names()                 # all registered agent names

_CAT_OPTS = (
    [{"label": "All categories", "value": ""}]
    + [{"label": cat, "value": cat} for cat in _BY_CAT]
    + [{"label": "⚡ Processing (no API key)", "value": "__processing__"}]
)

# LLM agent name → category label
_NAME_TO_CAT: dict[str, str] = {
    name: cat
    for cat, names in _BY_CAT.items()
    for name in names
}

_CAT_ICON: dict[str, str] = {
    "Processing":                 "bi-sliders",
    "Inversion":                  "bi-layers",
    "QC":                         "bi-shield-check",
    "Quality Control":            "bi-shield-check",
    "Reporting":                  "bi-file-earmark-text",
    "Interpretation & Reporting": "bi-lightbulb",
    "LLM":                        "bi-chat-dots",
    "Interpretation":             "bi-lightbulb",
    "Data Loading":               "bi-folder2-open",
    "Workflow":                   "bi-diagram-3",
}
_DEFAULT_ICON = "bi-robot"

_TYPE_META: dict[str, tuple[str, str, str]] = {
    "llm":        ("bi-stars",    "agent-badge-llm",        "LLM"),
    "processing": ("bi-cpu-fill", "agent-badge-processing", "Processing"),
}

# ── Chat constants ─────────────────────────────────────────────────────────

# (chip_key, display_label, prefill_text)
_QUICK_ACTIONS: list[tuple[str, str, str]] = [
    ("qc",         "Run full QC",          "Run quality control on all stations and identify bad SNR"),
    ("static",     "Correct static shift", "Detect and correct static shift across the survey"),
    ("dimension",  "Dimensionality",       "Classify 1-D / 2-D / 3-D structure from phase tensor"),
    ("inversion",  "Prep inversion",       "Prepare Occam2D inversion files for all profiles"),
    ("interpret",  "Interpret results",    "Interpret resistivity sections and map lithology"),
]

# keyword → [candidate agent names] (most specific first)
_KEYWORD_MAP: list[tuple[list[str], list[str]]] = [
    (["qc", "quality", "snr", "noise", "skew", "dead"],   ["QC Analysis",   "QC Quicklook"]),
    (["static", "shift"],                                  ["Static Shift",  "Static Shift (fast)"]),
    (["dimension", "1d", "2d", "3d", "strike"],            ["Dimensionality","Phase Analysis"]),
    (["occam", "inversion", "resistivity", "model"],       ["Inversion Prep","Occam2D"]),
    (["load", "import", "edi", "read", "file"],            ["MT Loader"]),
    (["interpret", "litholog", "correlat", "geolog"],      ["Interpretation"]),
    (["report", "export", "pdf", "summary", "docx"],       ["Report"]),
    (["pipeline", "workflow", "chain", "full process"],    ["Pipeline"]),
    (["denoise", "filter", "rpca", "emap"],                ["Denoising"]),
    (["tipper", "phase", "tensor", "mohr", "argand"],      ["Phase Analysis"]),
    (["forward", "synthetic", "model", "simulation"],      ["Forward Model"]),
    (["batch", "multiple", "profiles", "parallel"],        ["Batch Survey"]),
]


def _smart_route(message: str) -> dict:
    """Keyword router — no API key required. Returns a plan dict."""
    msg = message.lower()
    for keywords, candidates in _KEYWORD_MAP:
        if any(k in msg for k in keywords):
            valid = [a for a in candidates if a in AGENT_REGISTRY]
            if valid:
                return {
                    "agents":      valid,
                    "via":         "keyword",
                    "description": f"Matched by keywords → {', '.join(valid)}",
                }
    return {
        "agents":      [],
        "via":         "keyword",
        "description": (
            "I couldn't identify a specific workflow. "
            "Try: 'run QC', 'correct static shift', 'prep inversion', "
            "'classify dimensionality', or 'interpret results'."
        ),
    }


# ── Param form ─────────────────────────────────────────────────────────────

def _param_widget(key: str, spec: dict):
    ptype   = spec["type"]
    default = spec.get("default")
    wid     = {"type": "agent-param", "index": key}

    if ptype == "combo":
        return dbc.Select(
            id=wid,
            value=str(default) if default is not None else spec["options"][0],
            options=[{"label": o, "value": o} for o in spec["options"]],
            size="sm", className="param-input",
        )
    if ptype == "bool":
        return dbc.Checkbox(id=wid, value=bool(default), className="param-check mt-1")
    if ptype in ("float", "int"):
        lo, hi = spec.get("range", (None, None))
        return dbc.Input(
            id=wid, type="number", value=default,
            min=lo, max=hi, step=spec.get("step"),
            debounce=True, size="sm", className="param-input",
        )
    return dbc.Input(
        id=wid, type="text",
        value=str(default) if default is not None else "",
        debounce=True, size="sm", className="param-input",
    )


def _build_param_form(name: str | None) -> html.Div:
    entry = get_entry(name) if name else None
    if entry is None:
        return html.Div("No agent selected.", className="param-label")

    params = entry.get("params", {})
    if not params:
        return html.Div(
            [html.I(className="bi bi-info-circle me-1"),
             "No configurable parameters."],
            className="agent-no-params",
        )

    rows = []
    for key, spec in params.items():
        label = spec.get("label", key)
        tip   = spec.get("tip", "")
        lo, hi = (
            spec.get("range", (None, None))
            if spec["type"] in ("float", "int") else (None, None)
        )
        range_hint = f"  [{lo}–{hi}]" if lo is not None else ""
        rows.append(
            html.Div([
                html.Label(
                    label + range_hint,
                    title=tip or label,
                    className="param-label",
                    style={"flex": "0 0 120px", "cursor": "help" if tip else "default"},
                ),
                _param_widget(key, spec),
            ], className="param-row"),
        )
    return html.Div(rows, className="params-list")


# ── Agent info block ───────────────────────────────────────────────────────

def _agent_info_block(name: str | None) -> list:
    entry = get_entry(name) if name else None
    if entry is None:
        return [html.Span("No agent selected.", className="param-label")]

    atype = entry.get("type", "processing")
    icon_cls, badge_cls, badge_lbl = _TYPE_META.get(
        atype, ("bi-robot", "agent-badge-processing", atype.title())
    )
    cat  = entry.get("category", _NAME_TO_CAT.get(name or "", ""))
    desc = entry.get("description", "")

    pills = [
        html.Span(
            [html.I(className=f"bi {icon_cls} me-1"), badge_lbl],
            className=f"agent-type-badge {badge_cls}",
        ),
    ]
    if cat:
        pills.append(html.Span(cat, className="agent-cat-pill"))

    return [
        html.Div(pills, className="agent-info-row"),
        *([] if not desc else [html.P(desc, className="agent-desc-text")]),
    ]


# ── Agent list (left column) ───────────────────────────────────────────────

def _grouped_names(names: list[str]) -> list[tuple[str, list[str]]]:
    """Return [(section_label, [names])] preserving registry insertion order."""
    groups: dict[str, list[str]] = {}
    proc: list[str] = []
    for name in names:
        entry = AGENT_REGISTRY.get(name, {})
        if entry.get("type") == "processing":
            proc.append(name)
        else:
            cat = entry.get("category", "Other")
            groups.setdefault(cat, []).append(name)
    result = [(cat, nms) for cat, nms in groups.items() if nms]
    if proc:
        result.append(("⚡ Processing", proc))
    return result


def _agent_list(names: list[str], selected: str | None) -> html.Div:
    selected = selected or (names[0] if names else None)
    grouped  = _grouped_names(names)
    items: list = []
    for section, group_names in grouped:
        items.append(html.Div(section, className="agent-list-section"))
        for name in group_names:
            atype    = (AGENT_REGISTRY.get(name) or {}).get("type", "processing")
            is_proc  = atype == "processing"
            badge_cl = "agent-list-badge-proc" if is_proc else "agent-list-badge-llm"
            badge_lb = "⚡" if is_proc else "LLM"
            items.append(
                html.Button(
                    [
                        html.Span(name, className="agent-list-name"),
                        html.Span(badge_lb, className=badge_cl),
                    ],
                    id={"type": "agent-card", "index": name},
                    className="agent-list-row" + (" selected" if name == selected else ""),
                    n_clicks=0,
                    title=name,
                )
            )
    if not items:
        return html.Div("No agents found.", className="agent-no-params")
    return html.Div(items)


# ── Run history ────────────────────────────────────────────────────────────

def _render_history(history: list | None) -> list:
    if not history:
        return [html.Div("No runs yet.", className="agents-history-empty")]
    items = []
    for i, run in enumerate(history[:15]):
        ok      = run.get("status") == "success"
        icon    = "bi-check-circle-fill" if ok else "bi-x-circle-fill"
        icon_cl = "history-success" if ok else "history-error"
        items.append(
            html.Button(
                [
                    html.I(className=f"bi {icon} {icon_cl} me-1"),
                    html.Span(run.get("name", ""), className="history-name"),
                    html.Span(run.get("timestamp", ""), className="history-time"),
                ],
                id={"type": "history-item", "index": i},
                className="agents-history-item",
                n_clicks=0,
                title=run.get("name", ""),
            )
        )
    return items


# ── Page layout ───────────────────────────────────────────────────────────

# ── Data status bar ────────────────────────────────────────────────────────

def _build_data_bar(store_data: dict | None) -> list:
    """Return children for the data status bar from lightweight STORE_DATA."""
    if not store_data or not store_data.get("n_stations"):
        return [
            html.Span(className="agents-dot agents-dot-empty"),
            html.Span(
                "No survey loaded — use the Load Data dialog to enable agents.",
                className="agents-data-msg agents-data-msg-empty",
            ),
        ]

    n_sta    = store_data["n_stations"]
    n_lines  = store_data.get("n_lines", 1)
    data_dir = store_data.get("data_dir", "")
    line_counts: dict = store_data.get("line_counts", {})

    chips: list = [html.Span(className="agents-dot agents-dot-ok")]
    chips.append(html.Span(f"{n_sta} station{'s' if n_sta != 1 else ''}",
                            className="agents-data-chip"))
    chips.append(html.Span("·", className="agents-data-sep"))
    chips.append(html.Span(f"{n_lines} profile{'s' if n_lines != 1 else ''}",
                            className="agents-data-chip"))

    if line_counts:
        summary = ", ".join(
            f"{ln}: {cnt}" for ln, cnt in list(line_counts.items())[:3]
        )
        if len(line_counts) > 3:
            summary += f" +{len(line_counts) - 3} more"
        chips.append(html.Span("·", className="agents-data-sep"))
        chips.append(html.Span(summary, className="agents-data-chip agents-data-lines",
                               title=str(line_counts)))

    if data_dir and data_dir not in ("[uploaded]", "[browsed]"):
        short = ("…" + data_dir[-38:]) if len(data_dir) > 40 else data_dir
        chips.append(html.Span("·", className="agents-data-sep"))
        chips.append(html.Span(short, className="agents-data-chip agents-data-path",
                               title=data_dir))

    return chips


# ── Chat helpers ───────────────────────────────────────────────────────────

def _welcome_message() -> dict:
    return {
        "id":       "welcome",
        "role":     "assistant",
        "content": (
            "Hello! I'm your MT workflow assistant. "
            "I can run QC, correct static shift, prepare inversion files, "
            "classify dimensionality, or interpret resistivity results.\n\n"
            "Describe what you'd like to do with your loaded survey — "
            "or pick a quick action below."
        ),
        "plan":           None,
        "plan_status":    None,
        "result_log":     None,
        "result_src":     None,
        "result_summary": None,
        "timestamp":      "",
    }


def _render_chat(messages: list) -> list:
    """Render a list of message dicts into Dash components."""
    bubbles: list = []
    for msg in messages:
        role    = msg["role"]
        content = msg.get("content", "")
        ts      = msg.get("timestamp", "")
        plan    = msg.get("plan")
        p_stat  = msg.get("plan_status")
        r_log   = msg.get("result_log")
        r_src   = msg.get("result_src")
        r_sum   = msg.get("result_summary")
        msg_id  = msg.get("id", "")

        if role == "user":
            bubbles.append(
                html.Div([
                    html.Div([
                        html.I(className="bi bi-person-fill me-1"),
                        html.Span(ts, className="chat-ts"),
                    ], className="chat-bubble-meta"),
                    html.P(content, className="chat-bubble-text"),
                ], className="chat-bubble chat-bubble-user")
            )
            continue

        # ── assistant bubble ──────────────────────────────────────────────
        inner: list = [
            html.Div([
                html.I(className="bi bi-robot me-1 chat-robot-icon"),
                html.Span("Assistant", className="chat-assistant-label"),
                html.Span(ts, className="chat-ts"),
            ], className="chat-bubble-meta"),
        ]

        if content:
            # Preserve newlines
            parts = content.split("\n\n")
            inner.extend(html.P(p, className="chat-bubble-text") for p in parts if p.strip())

        # Plan card
        if plan and plan.get("agents") and p_stat in (None, "pending"):
            agent_items = [
                html.Li([
                    html.Span(a, className="chat-plan-agent"),
                    html.Span(
                        " LLM" if (AGENT_REGISTRY.get(a) or {}).get("type") == "llm" else " ⚡",
                        className=(
                            "chat-plan-badge-llm"
                            if (AGENT_REGISTRY.get(a) or {}).get("type") == "llm"
                            else "chat-plan-badge-proc"
                        ),
                    ),
                ]) for a in plan["agents"]
            ]
            via_note = (
                html.Span(" · AI-matched", className="chat-via-llm")
                if plan.get("via") == "llm"
                else html.Span(" · keyword-matched", className="chat-via-kw")
            )
            inner.append(
                html.Div([
                    html.Div(
                        [html.I(className="bi bi-list-check me-1"), "Proposed plan", via_note],
                        className="chat-plan-header",
                    ),
                    html.Ol(agent_items, className="chat-plan-list"),
                    html.Div([
                        html.Button(
                            [html.I(className="bi bi-play-fill me-1"), "Confirm & Run"],
                            id={"type": "chat-confirm", "index": msg_id},
                            className="chat-confirm-btn",
                            n_clicks=0,
                        ),
                        html.Button(
                            "Cancel",
                            id={"type": "chat-cancel", "index": msg_id},
                            className="chat-cancel-btn",
                            n_clicks=0,
                        ),
                    ], className="chat-plan-actions"),
                ], className="chat-plan-card")
            )

        # Cancelled plan note
        if p_stat == "cancelled":
            inner.append(html.P("Plan cancelled.", className="chat-cancelled-note"))

        # Running indicator
        if p_stat == "running":
            inner.append(
                html.Div([
                    dbc.Spinner(size="sm", color="warning"),
                    html.Span(" Running agents…", className="chat-running-note ms-2"),
                ], className="d-flex align-items-center mt-2")
            )

        # Inline results
        if r_log or r_sum or r_src:
            result_children: list = []
            if r_sum:
                result_children.append(html.P(r_sum, className="chat-result-summary"))
            if r_log:
                # open=True: log is always visible after execution —
                # html.Details collapses on every Dash re-render otherwise.
                result_children.append(
                    html.Details([
                        html.Summary("Execution log", className="chat-log-summary"),
                        html.Pre(r_log[:2000], className="chat-result-log"),
                    ], className="chat-result-details", open=True)
                )
            if r_src and r_src != empty_src(dark=True):
                result_children.append(
                    html.Img(src=r_src, className="chat-result-img")
                )
            inner.append(html.Div(result_children, className="chat-result-block"))

        bubbles.append(
            html.Div(inner, className="chat-bubble chat-bubble-assistant")
        )

    return bubbles


def _chat_layout() -> html.Div:
    """Return the full chat panel (message area + input + quick actions)."""
    return html.Div([
        # Scrollable conversation
        html.Div(
            _render_chat([_welcome_message()]),
            id=_ID_CHAT_MESSAGES,
            className="agents-chat-messages",
        ),

        # Input row
        html.Div([
            dbc.Textarea(
                id=_ID_CHAT_INPUT,
                placeholder="Describe what you want to do with your survey data…",
                rows=2,
                className="agents-chat-textarea",
                debounce=False,
                n_blur=0,
            ),
            dbc.Button(
                [html.I(className="bi bi-send-fill me-1"), "Send"],
                id=_ID_CHAT_SEND,
                color="primary",
                size="sm",
                className="agents-chat-send-btn",
                n_clicks=0,
            ),
        ], className="agents-chat-input-row"),

        # Quick action chips
        html.Div([
            html.Span("Quick actions:", className="agents-quick-label"),
            *[
                html.Button(
                    label,
                    id={"type": "quick-action", "index": key},
                    className="agents-quick-chip",
                    n_clicks=0,
                    title=prefill,
                )
                for key, label, prefill in _QUICK_ACTIONS
            ],
        ], className="agents-quick-row"),

        dcc.Store(id=_ID_CHAT_STORE, storage_type="session"),
    ], className="agents-chat-layout")


# ── Page layout ───────────────────────────────────────────────────────────

def layout() -> html.Div:
    total = len(_ALL_NAMES)
    cats  = len(_BY_CAT)
    first = _ALL_NAMES[0] if _ALL_NAMES else None

    # ── Left column — searchable, filterable agent list ────────────────
    col_list = html.Div([
        html.Div([
            dbc.Input(
                id=_ID_SEARCH,
                placeholder="Search agents…",
                type="search",
                size="sm",
                debounce=True,
                className="agents-search-input mb-1",
            ),
            dbc.Select(
                id=IDs.AGENTS_CAT,
                options=_CAT_OPTS,
                value="",
                size="sm",
            ),
        ], className="agents-col-filters"),

        html.Div(
            _agent_list(_ALL_NAMES, first),
            id=_ID_LIST_WRAP,
            className="agents-list-scroll",
        ),

        # Hidden select keeps AGENTS_NAME as the stable selection signal
        dbc.Select(
            id=IDs.AGENTS_NAME,
            options=[{"label": n, "value": n} for n in _ALL_NAMES],
            value=first,
            size="sm",
            style={"display": "none"},
        ),
    ], className="agents-col-list")

    # ── Middle column — agent detail + params + run button ─────────────
    col_detail = html.Div([
        # Header: agent name + type badge + description
        html.Div([
            html.H6(first or "Select an agent",
                    id=_ID_DETAIL_NAME,
                    className="agents-detail-name"),
            html.Div(_agent_info_block(first), id=IDs.AGENTS_DESC),
        ], className="agents-detail-header"),

        # Scrollable body: line filter → station filter → divider → parameters
        html.Div([
            # ── Line filter (active lines only) ───────────────────────────
            html.Div([
                html.Div(
                    [html.I(className="bi bi-layers me-1"), "Filter by Line"],
                    className="agents-section-label",
                ),
                html.Div([
                    html.Button("All",  id=_ID_LINE_SEL_ALL,   className="agents-mini-btn", n_clicks=0),
                    html.Button("None", id=_ID_LINE_DESEL_ALL,  className="agents-mini-btn", n_clicks=0),
                ], className="agents-filter-btn-row"),
            ], className="agents-filter-header-row"),
            dcc.Dropdown(
                id=_ID_LINE_FILTER,
                options=[],
                value=None,
                multi=True,
                placeholder="All active lines",
                className="agents-station-filter mb-1",
                optionHeight=28,
            ),

            # ── Station filter (filtered by selected lines) ───────────────
            html.Div([
                html.Div(
                    [html.I(className="bi bi-geo-alt me-1"), "Station Filter"],
                    className="agents-section-label",
                ),
                html.Div([
                    html.Button("All",  id=_ID_STN_SEL_ALL,   className="agents-mini-btn", n_clicks=0),
                    html.Button("None", id=_ID_STN_DESEL_ALL,  className="agents-mini-btn", n_clicks=0),
                ], className="agents-filter-btn-row"),
            ], className="agents-filter-header-row"),
            dcc.Dropdown(
                id=_ID_STATION_FILTER,
                options=[],
                value=None,
                multi=True,
                placeholder="All stations from active lines",
                className="agents-station-filter",
                optionHeight=28,
            ),
            html.Hr(className="agents-divider-sm"),

            # Parameters form
            html.Div(
                [html.I(className="bi bi-sliders me-1"), "Parameters"],
                className="agents-section-label",
            ),
            html.Div(_build_param_form(first), id=IDs.AGENTS_PARAM_FORM),
        ], className="agents-detail-params"),

        # Run button pinned to bottom
        html.Div([
            dbc.Button(
                [html.I(className="bi bi-play-fill me-1"), "Run Agent"],
                id=IDs.BTN_AGENTS_RUN,
                color="warning",
                size="sm",
                className="w-100 mb-2",
                disabled=True,          # enabled once data is loaded
            ),
            dbc.Spinner(
                html.Div(id=IDs.AGENTS_SPINNER),
                size="sm",
                color="warning",
            ),
        ], className="agents-detail-run"),
    ], className="agents-col-detail")

    # ── Right column — run history + result tabs ────────────────────────
    col_output = html.Div([
        # Run history strip
        html.Div([
            html.Div(
                [html.I(className="bi bi-clock-history me-1"), "Run History"],
                className="agents-section-label",
            ),
            html.Div(
                _render_history(None),
                id=_ID_HISTORY,
                className="agents-history-list",
            ),
        ], className="agents-output-history"),

        # Result tabs
        html.Div([
            html.Div(
                [html.I(className="bi bi-terminal me-1"), "Last Result"],
                className="agents-section-label",
            ),
            dbc.Tabs([
                dbc.Tab(
                    html.Pre(
                        id=IDs.AGENTS_OUT,
                        children="Agent output will appear here…",
                        className="log-area",
                        style={
                            "height": "370px", "margin": "8px",
                            "whiteSpace": "pre-wrap", "wordBreak": "break-word",
                        },
                    ),
                    label="Log", tab_id="agent-tab-log",
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
                        style={"minHeight": "350px", "margin": "8px"},
                    ),
                    label="Figure", tab_id="agent-tab-fig",
                ),
                dbc.Tab(
                    html.Pre(
                        id=_ID_SUMMARY,
                        children="Summary will appear here…",
                        className="log-area",
                        style={
                            "height": "370px", "margin": "8px",
                            "whiteSpace": "pre-wrap", "wordBreak": "break-word",
                        },
                    ),
                    label="Summary", tab_id="agent-tab-sum",
                ),
            ], active_tab="agent-tab-log"),
        ], className="agents-output-result"),

        dcc.Store(id=IDs.AGENTS_STORE, storage_type="memory"),
    ], className="agents-col-output")

    # ── View toggle bar ────────────────────────────────────────────────
    view_toggle = html.Div([
        html.Button(
            [html.I(className="bi bi-grid-3x3-gap me-1"), "Agent Runner"],
            id=_ID_TAB_RUNNER,
            className="agents-view-tab active",
            n_clicks=0,
        ),
        html.Button(
            [html.I(className="bi bi-chat-text me-1"), "Chat"],
            id=_ID_TAB_CHAT,
            className="agents-view-tab",
            n_clicks=0,
        ),
    ], className="agents-view-toggle")

    return html.Div([
        _command_bar(
            "agents", "AI Agents",
            f"{total} agents · {cats} categories · LLM · Processing · Inversion",
        ),
        view_toggle,
        # Data status bar — populated by _update_data_status callback
        html.Div(
            _build_data_bar(None),
            id=_ID_DATA_BAR,
            className="agents-data-bar",
        ),
        # ── Agent Runner view (default) ──────────────────────────────────
        html.Div(
            [col_list, col_detail, col_output],
            id=_ID_VIEW_RUNNER,
            className="agents-3col-layout",
        ),
        # ── Chat view (hidden until tab click) ───────────────────────────
        html.Div(
            _chat_layout(),
            id=_ID_VIEW_CHAT,
            style={"display": "none", "flex": "1", "overflow": "hidden"},
        ),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


# ── Callbacks ─────────────────────────────────────────────────────────────

def register_callbacks(app) -> None:
    from dash import callback

    # 1. Card/list-row click → update hidden select
    @callback(
        Output(IDs.AGENTS_NAME, "value"),
        Input({"type": "agent-card", "index": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def _card_to_select(_n_clicks_list):
        triggered = ctx.triggered_id
        if isinstance(triggered, dict) and "index" in triggered:
            return triggered["index"]
        return no_update

    # 2. Selected agent → highlight list row + refresh info + rebuild param form
    @callback(
        Output({"type": "agent-card", "index": ALL}, "className"),
        Output(IDs.AGENTS_DESC,        "children"),
        Output(IDs.AGENTS_PARAM_FORM,  "children"),
        Output(_ID_DETAIL_NAME,        "children"),
        Input(IDs.AGENTS_NAME, "value"),
        prevent_initial_call=True,
    )
    def _on_agent_selected(name):
        if not name:
            return no_update, no_update, no_update, no_update
        # ctx.outputs_list[0] contains the actual matched DOM elements
        # so this is correct even when the list is filtered to a subset
        classes = [
            "agent-list-row selected" if item["id"]["index"] == name
            else "agent-list-row"
            for item in ctx.outputs_list[0]
        ]
        return (
            classes,
            _agent_info_block(name),
            _build_param_form(name),
            name,
        )

    # 3. Category filter + search → rebuild list + reset selection
    @callback(
        Output(_ID_LIST_WRAP,       "children"),
        Output(IDs.AGENTS_NAME,     "value",   allow_duplicate=True),
        Input(IDs.AGENTS_CAT,  "value"),
        Input(_ID_SEARCH,      "value"),
        prevent_initial_call=True,
    )
    def _filter_list(cat, search):
        if cat == "__processing__":
            filtered = processing_agents()
        elif cat:
            filtered = _BY_CAT.get(cat, [])
        else:
            filtered = _ALL_NAMES

        if search:
            q = search.strip().lower()
            filtered = [n for n in filtered if q in n.lower()]

        new_sel = filtered[0] if filtered else None
        return _agent_list(filtered, new_sel), new_sel

    # 4. STORE_DATA / STORE_ACTIVE_LINES → data bar + line opts + station opts + run button
    @callback(
        Output(_ID_DATA_BAR,         "children"),
        Output(_ID_LINE_FILTER,      "options"),
        Output(_ID_STATION_FILTER,   "options"),
        Output(IDs.BTN_AGENTS_RUN,   "disabled"),
        Input(IDs.STORE_DATA,        "data"),
        Input(IDs.STORE_ACTIVE_LINES,"data"),
        prevent_initial_call=False,
    )
    def _update_data_status(store_data, active_lines_store):
        if not store_data or not store_data.get("n_stations"):
            return _build_data_bar(None), [], [], True

        records      = store_data.get("station_records", [])
        active_lines = (active_lines_store or {}).get("active", [])

        # If active_lines is empty (store not yet initialised), show everything
        if not active_lines:
            active_lines = list({r.get("Line", "") for r in records if r.get("Line")})

        # Build line options — only active lines
        line_opts = [{"label": ln, "value": ln} for ln in active_lines if ln]

        # Build station options — only stations belonging to active lines
        if active_lines:
            active_set = set(active_lines)
            stn_opts = [
                {"label": r["ID"], "value": r["ID"]}
                for r in records
                if r.get("Line", "") in active_set or not r.get("Line")
            ]
        else:
            stn_opts = [{"label": r["ID"], "value": r["ID"]} for r in records]

        return _build_data_bar(store_data), line_opts, stn_opts, False

    # 4b. Line "All" button → select all line options
    @callback(
        Output(_ID_LINE_FILTER, "value", allow_duplicate=True),
        Input(_ID_LINE_SEL_ALL, "n_clicks"),
        State(_ID_LINE_FILTER,  "options"),
        prevent_initial_call=True,
    )
    def _line_select_all(_n, opts):
        if not _n or not opts:
            return no_update
        return [o["value"] for o in opts]

    # 4c. Line "None" button → clear line selection
    @callback(
        Output(_ID_LINE_FILTER, "value", allow_duplicate=True),
        Input(_ID_LINE_DESEL_ALL, "n_clicks"),
        prevent_initial_call=True,
    )
    def _line_deselect_all(_n):
        return [] if _n else no_update

    # 4d. Selected lines → auto-fill matching stations in station filter
    @callback(
        Output(_ID_STATION_FILTER, "value",   allow_duplicate=True),
        Output(_ID_STATION_FILTER, "options", allow_duplicate=True),
        Input(_ID_LINE_FILTER,     "value"),
        State(IDs.STORE_DATA,      "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def _line_to_stations(selected_lines, store_data, active_lines_store):
        if not store_data:
            return no_update, no_update

        records      = store_data.get("station_records", [])
        active_lines = (active_lines_store or {}).get("active", [])
        if not active_lines:
            active_lines = list({r.get("Line", "") for r in records if r.get("Line")})
        active_set   = set(active_lines)

        # Always limit visible station options to active-line stations
        stn_opts = [
            {"label": r["ID"], "value": r["ID"]}
            for r in records
            if r.get("Line", "") in active_set or not r.get("Line")
        ]

        if not selected_lines:
            # No line filter → clear station selection (all stations implicitly)
            return None, stn_opts

        # Select every station belonging to any of the selected lines
        sel_set  = set(selected_lines)
        selected = [
            r["ID"] for r in records
            if r.get("Line", "") in sel_set
        ]
        return selected or None, stn_opts

    # 4e. Station "All" button → select all station options
    @callback(
        Output(_ID_STATION_FILTER, "value", allow_duplicate=True),
        Input(_ID_STN_SEL_ALL,     "n_clicks"),
        State(_ID_STATION_FILTER,  "options"),
        prevent_initial_call=True,
    )
    def _station_select_all(_n, opts):
        if not _n or not opts:
            return no_update
        return [o["value"] for o in opts]

    # 4f. Station "None" button → clear station selection
    @callback(
        Output(_ID_STATION_FILTER, "value", allow_duplicate=True),
        Input(_ID_STN_DESEL_ALL,   "n_clicks"),
        prevent_initial_call=True,
    )
    def _station_deselect_all(_n):
        return [] if _n else no_update

    # 5. Run Agent → collect form params + station filter → execute → history + output
    @callback(
        Output(IDs.AGENTS_OUT,     "children"),
        Output(IDs.IMG_AGENTS,     "src"),
        Output(IDs.AGENTS_SPINNER, "children"),
        Output(IDs.AGENTS_STORE,   "data"),
        Output(_ID_HISTORY,        "children"),
        Output(_ID_SUMMARY,        "children"),
        Input(IDs.BTN_AGENTS_RUN,  "n_clicks"),
        State(IDs.AGENTS_NAME,     "value"),
        State({"type": "agent-param", "index": ALL}, "value"),
        State(_ID_STATION_FILTER,  "value"),
        State(IDs.SESSION_ID,      "data"),
        State(IDs.AGENTS_STORE,    "data"),
        prevent_initial_call=True,
    )
    def _run_agent(_n, agent_name, _param_vals, station_filter, session_id, history):
        from pycsamt.app.web.callbacks.agents import (
            _exec_agent,
        )

        if not agent_name:
            return "No agent selected.", no_update, "", no_update, no_update, no_update
        if not session_id:
            return "⚠ Session not initialised — refresh.", no_update, "", no_update, no_update, no_update

        # Reconstruct {param_key: value} from pattern-matched States (index 1).
        # States are indexed in argument order: 0=AGENTS_NAME, 1=ALL params,
        # 2=station_filter, 3=SESSION_ID, 4=AGENTS_STORE.
        params: dict = {}
        if ctx.states_list and len(ctx.states_list) > 1:
            for item in ctx.states_list[1]:
                key = item["id"]["index"]
                val = item.get("value")
                if val is not None:
                    params[key] = val

        # Fill any remaining keys from registry defaults
        for k, v in default_params(agent_name).items():
            params.setdefault(k, v)

        # Type-cast per registry spec
        entry = get_entry(agent_name)
        for k, spec in (entry or {}).get("params", {}).items():
            if k not in params or params[k] is None:
                continue
            try:
                ptype = spec["type"]
                if ptype == "int":
                    params[k] = int(params[k])
                elif ptype == "float":
                    params[k] = float(params[k])
                elif ptype == "bool":
                    params[k] = bool(params[k])
            except (ValueError, TypeError):
                params[k] = spec.get("default")

        # Pass station subset (None = all stations)
        stations = station_filter or None

        log, src, summary, err = _exec_agent(agent_name, params, session_id,
                                             stations=stations)

        record = {
            "name":      agent_name,
            "timestamp": datetime.now().strftime("%H:%M:%S"),
            "status":    "error" if err else "success",
            "log":       log,
            "src":       src,
            "summary":   summary,
            "params":    json.dumps(params, indent=2, default=str),
            "stations":  stations,
        }
        updated = [record] + (history or [])[:19]

        return log, src, "", updated, _render_history(updated), summary

    # 7. Click a history item → restore that run's outputs
    @callback(
        Output(IDs.AGENTS_OUT,  "children",    allow_duplicate=True),
        Output(IDs.IMG_AGENTS,  "src",         allow_duplicate=True),
        Output(_ID_SUMMARY,     "children",    allow_duplicate=True),
        Input({"type": "history-item", "index": ALL}, "n_clicks"),
        State(IDs.AGENTS_STORE, "data"),
        prevent_initial_call=True,
    )
    def _load_history(_clicks, history):
        triggered = ctx.triggered_id
        if not isinstance(triggered, dict) or not history:
            return no_update, no_update, no_update
        idx = triggered.get("index", -1)
        if not (0 <= idx < len(history)):
            return no_update, no_update, no_update
        run = history[idx]
        return (
            run.get("log", ""),
            run.get("src", empty_src(dark=True)),
            run.get("summary", ""),
        )

    # ── Chat callbacks ─────────────────────────────────────────────────────

    # 8. Tab toggle: Runner ↔ Chat
    @callback(
        Output(_ID_VIEW_RUNNER, "style"),
        Output(_ID_VIEW_CHAT,   "style"),
        Output(_ID_TAB_RUNNER,  "className"),
        Output(_ID_TAB_CHAT,    "className"),
        Input(_ID_TAB_RUNNER,   "n_clicks"),
        Input(_ID_TAB_CHAT,     "n_clicks"),
        prevent_initial_call=True,
    )
    def _toggle_view(_nr, _nc):
        runner_visible = {"display": "grid",  "flex": "1", "overflow": "hidden", "minHeight": "0"}
        chat_visible   = {"display": "flex",  "flex": "1", "overflow": "hidden", "flexDirection": "column"}
        hidden         = {"display": "none"}
        if ctx.triggered_id == _ID_TAB_CHAT:
            return hidden, chat_visible, "agents-view-tab", "agents-view-tab active"
        return runner_visible, hidden, "agents-view-tab active", "agents-view-tab"

    # 9. Quick action chip → prefill the chat input
    @callback(
        Output(_ID_CHAT_INPUT, "value"),
        Input({"type": "quick-action", "index": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def _quick_action_fill(_clicks):
        triggered = ctx.triggered_id
        if not isinstance(triggered, dict):
            return no_update
        key = triggered.get("index", "")
        for k, _label, prefill in _QUICK_ACTIONS:
            if k == key:
                return prefill
        return no_update

    # 10. Send message → parse intent → propose plan → render chat
    @callback(
        Output(_ID_CHAT_MESSAGES, "children"),
        Output(_ID_CHAT_STORE,    "data"),
        Output(_ID_CHAT_INPUT,    "value"),
        Input(_ID_CHAT_SEND,      "n_clicks"),
        State(_ID_CHAT_INPUT,     "value"),
        State(_ID_CHAT_STORE,     "data"),
        State(IDs.SESSION_ID,     "data"),
        State(IDs.STORE_DATA,     "data"),
        prevent_initial_call=True,
    )
    def _send_message(_n, user_text, chat_history, session_id, store_data):
        if not user_text or not user_text.strip():
            return no_update, no_update, no_update

        history: list = chat_history or [_welcome_message()]
        ts = datetime.now().strftime("%H:%M:%S")

        # Add user message
        history = history + [{
            "id":             str(uuid.uuid4())[:8],
            "role":           "user",
            "content":        user_text.strip(),
            "plan":           None,
            "plan_status":    None,
            "result_log":     None,
            "result_src":     None,
            "result_summary": None,
            "timestamp":      ts,
        }]

        # Try LLM context agent, fall back to keyword router
        plan: dict = {}
        llm_content = ""
        has_data = bool(store_data and store_data.get("n_stations"))

        try:
            if has_data:
                from pycsamt.agents import ContextInputAgent
                from pycsamt.app.web.cache import cache_get
                sites = cache_get(session_id)
                if sites is not None:
                    agent  = ContextInputAgent()
                    result = agent.execute({"sites": sites, "request": user_text.strip()})
                    if result and getattr(result, "status", "") == "ok":
                        raw = getattr(result, "summary", "") or ""
                        # Extract agent names mentioned in the LLM response
                        matched = [
                            n for n in _ALL_NAMES
                            if n.lower() in raw.lower()
                        ]
                        if matched:
                            plan = {"agents": matched[:4], "via": "llm",
                                    "description": raw[:200]}
                            llm_content = "Based on your request, here's what I'll run:"
                        else:
                            plan = _smart_route(user_text)
                    else:
                        plan = _smart_route(user_text)
                else:
                    plan = _smart_route(user_text)
            else:
                plan = _smart_route(user_text)
        except Exception:
            plan = _smart_route(user_text)

        # Build assistant reply
        if plan.get("agents"):
            content = llm_content or "Based on your request, here's what I'll run:"
            plan_status = "pending"
        else:
            content = plan.get("description", "I couldn't identify a workflow.")
            plan_status = None

        if not has_data:
            content = (
                "No survey data is loaded. "
                "Please load EDI files first, then I can run agents on your data."
            )
            plan = {}
            plan_status = None

        history = history + [{
            "id":             str(uuid.uuid4())[:8],
            "role":           "assistant",
            "content":        content,
            "plan":           plan if plan.get("agents") else None,
            "plan_status":    plan_status,
            "result_log":     None,
            "result_src":     None,
            "result_summary": None,
            "timestamp":      datetime.now().strftime("%H:%M:%S"),
        }]

        return _render_chat(history), history, ""

    # 11. Confirm plan → execute agents → embed results in chat bubble
    @callback(
        Output(_ID_CHAT_MESSAGES, "children",  allow_duplicate=True),
        Output(_ID_CHAT_STORE,    "data",      allow_duplicate=True),
        Input({"type": "chat-confirm", "index": ALL}, "n_clicks"),
        State(_ID_CHAT_STORE,    "data"),
        State(IDs.SESSION_ID,    "data"),
        State(_ID_STATION_FILTER,"value"),
        prevent_initial_call=True,
    )
    def _confirm_plan(_clicks, chat_history, session_id, station_filter):
        from pycsamt.app.web.callbacks.agents import (
            _exec_agent,
        )

        # Pattern-match ALL fires for every newly-added button (n_clicks=0).
        # Only proceed when a real click occurred.
        if not _clicks or not any(_clicks):
            return no_update, no_update

        triggered = ctx.triggered_id
        if not isinstance(triggered, dict) or not chat_history:
            return no_update, no_update

        target_id = triggered.get("index", "")
        history   = list(chat_history)

        # Locate the message that owns this confirm button
        msg_idx = next(
            (i for i, m in enumerate(history) if m.get("id") == target_id),
            -1,
        )
        if msg_idx < 0:
            return no_update, no_update

        msg  = dict(history[msg_idx])
        plan = msg.get("plan") or {}
        agents_to_run = plan.get("agents", [])

        if not agents_to_run:
            return no_update, no_update

        # Mark as running
        msg["plan_status"] = "running"
        history[msg_idx]   = msg

        stations = station_filter or None
        all_logs: list[str] = []
        last_src  = empty_src(dark=True)
        summaries: list[str] = []
        errors:    list[str] = []

        for agent_name in agents_to_run:
            if agent_name not in AGENT_REGISTRY:
                all_logs.append(f"⚠ Agent '{agent_name}' not found — skipped.")
                continue
            params = default_params(agent_name)
            log, src, summary, err = _exec_agent(
                agent_name, params, session_id, stations=stations
            )
            all_logs.append(f"── {agent_name} ──\n{log}")
            if src and src != empty_src(dark=True):
                last_src = src
            if summary:
                summaries.append(f"{agent_name}: {summary}")
            if err:
                errors.append(f"{agent_name}: {err}")

        combined_log = "\n\n".join(all_logs)
        combined_sum = "\n".join(summaries)
        status_icon  = "✓" if not errors else "✗"
        result_note  = (
            f"{status_icon} Ran {len(agents_to_run)} agent(s)."
            + (f" Errors: {len(errors)}." if errors else "")
        )

        msg["plan_status"]    = "done"
        msg["result_log"]     = combined_log
        msg["result_src"]     = last_src
        msg["result_summary"] = result_note + ("\n\n" + combined_sum if combined_sum else "")
        history[msg_idx]      = msg

        # Compact follow-up — only shown if there are agent-level errors to highlight
        if errors:
            history = history + [{
                "id":             str(uuid.uuid4())[:8],
                "role":           "assistant",
                "content":        f"⚠ {len(errors)} agent(s) reported errors (see log above).",
                "plan":           None,
                "plan_status":    None,
                "result_log":     None,
                "result_src":     None,
                "result_summary": None,
                "timestamp":      datetime.now().strftime("%H:%M:%S"),
            }]

        return _render_chat(history), history

    # 12. Cancel plan → mark cancelled → re-render
    @callback(
        Output(_ID_CHAT_MESSAGES, "children",  allow_duplicate=True),
        Output(_ID_CHAT_STORE,    "data",      allow_duplicate=True),
        Input({"type": "chat-cancel", "index": ALL}, "n_clicks"),
        State(_ID_CHAT_STORE, "data"),
        prevent_initial_call=True,
    )
    def _cancel_plan(_clicks, chat_history):
        # Guard against spurious fires from newly-mounted cancel buttons (n_clicks=0).
        if not _clicks or not any(_clicks):
            return no_update, no_update

        triggered = ctx.triggered_id
        if not isinstance(triggered, dict) or not chat_history:
            return no_update, no_update

        target_id = triggered.get("index", "")
        history   = list(chat_history)
        for i, m in enumerate(history):
            if m.get("id") == target_id:
                msg = dict(m)
                msg["plan_status"] = "cancelled"
                history[i] = msg
                break

        return _render_chat(history), history
