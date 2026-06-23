# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Chat dispatch callbacks.

Flow
----
1. User types + clicks Send.
2. Message stored; thinking bubble shown.
3. Background thread runs agent chain.
4. dcc.Interval polls shared _JOBS dict.
5. On completion: append agent bubble,
   disable interval.

Figures embedded in agent output are stored
in STORE_FIGS keyed by UUID; each rendered
as an am-fig-card with view + export buttons.
"""

from __future__ import annotations

import base64
import io
import re
import threading
import time
import uuid
from datetime import datetime
from typing import Any

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
# Silence plt.show() globally — agent code
# must never open OS windows in a web app.
plt.show = lambda *a, **kw: None

from contextlib import nullcontext as _nullctx

from dash import ALL, Input, Output, State
from dash import ctx, html, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from .._ids import IDs

# ── shared job registry ────────────────────────────
_JOBS: dict[str, dict] = {}
_JOBS_LOCK = threading.Lock()

# Corrected-sites cache keyed by job ID.
# Populated after correction workflows so the
# post-processing modal can export EDI files.
_CORR_CACHE: dict[str, Any] = {}

# Workflows that produce corrected_sites data.
_CORRECTION_WFLOWS = frozenset({
    "static_shift",
    "denoise",
    "qc",
    "pre_inversion",
    "full",
})

# Plotting tasks handled by the lightweight PlotAgent (not the orchestrator):
# rho/phi sounding curves, scalar phase pseudo-section, and phase-tensor (Φ)
# ellipse pseudo-section. Each maps to a PlotAgent "kind".
_PLOT_KIND = {
    "rhophi":         "rhophi",
    "phase_psection": "phase_psection",
    "pt_psection":    "pt_psection",
    "tipper_plot":    "tipper",
    "phase_tensor_map": "pt_map",
    "station_response": "station_response",
    "strike_profile":  "strike_profile",
}
_PLOT_WORKFLOWS = frozenset(_PLOT_KIND)

# Analysis tools handled by the lightweight ToolAgent (table + figure).
_TOOL_KIND = {
    "strike":         "strike",
    "dimensionality": "dimensionality",
    "validator":      "validator",
    "coords":         "coords",
    "elevation":      "elevation",
    "converter":      "converter",
    "batch_export":   "batch_export",
}
_TOOL_WORKFLOWS = frozenset(_TOOL_KIND)

# ── response kinds ─────────────────────────────────
# Each agent response declares its KIND so the chat
# bubble renders the right shape instead of guessing
# from which fields happen to be populated.
KIND_ANSWER   = "answer"     # Q&A about the package
KIND_CODE     = "code"       # generated script
KIND_WORKFLOW = "workflow"   # pipeline result (+figs)
KIND_CLARIFY  = "clarify"    # needs more info
KIND_META     = "meta"       # capabilities / chitchat
KIND_ERROR    = "error"      # could not proceed

# Per-kind header chip: (icon class, label, css colour var).
_KIND_HEADER: dict[str, tuple[str, str, str]] = {
    KIND_ANSWER:   ("bi-chat-left-text", "Answer", "var(--blue)"),
    KIND_CODE:     ("bi-code-slash", "Generated code", "var(--green)"),
    KIND_WORKFLOW: ("bi-diagram-3", "Workflow result", "var(--blue)"),
    KIND_CLARIFY:  ("bi-question-circle", "Needs input", "var(--yellow)"),
    KIND_META:     ("bi-stars", "pyCSAMT assistant", "var(--blue)"),
    KIND_ERROR:    ("bi-exclamation-triangle", "Couldn't proceed", "var(--red)"),
}


def _new_job() -> str:
    jid = str(uuid.uuid4())
    with _JOBS_LOCK:
        _JOBS[jid] = {
            "status": "running",
            "steps": [],
            "result": None,
            "figs": {},
            "error": None,
            "kind": None,
            "started": time.time(),
            "workflow": None,
        }
    return jid


# Friendly labels for the executing-message header.
_WF_RUNNING_LABEL: dict[str, str] = {
    "qc":                 "quality control",
    "static_shift":       "static-shift correction",
    "phase_analysis":     "phase tensor analysis",
    "denoise":            "denoising",
    "tipper":             "tipper analysis",
    "sensitivity":        "sensitivity / DOI analysis",
    "rotation":           "tensor rotation",
    "freq_decimation":    "frequency decimation",
    "ai_inversion":       "1-D AI inversion",
    "inv1d":              "1-D AI inversion",
    "inv2d":              "2-D U-Net inversion",
    "inv3d":              "3-D GCN inversion",
    "pinn_inversion":     "PINN inversion",
    "hybrid_inversion":   "hybrid inversion",
    "ensemble_inversion": "ensemble inversion",
    "joint_inversion":    "joint inversion",
    "pre_inversion":      "inversion preparation",
    "modem":              "ModEM preparation",
    "report":             "survey report",
    "interpretation":     "geological interpretation",
    "interpret":          "geological interpretation",
    "forward":            "forward modelling",
    "full":               "full pipeline",
    "comparison":         "inversion comparison",
    "batch":              "batch processing",
    "code_gen":           "code generation",
}


def _fmt_elapsed(seconds: float) -> str:
    """Format seconds as M:SS for the executing header."""
    s = int(max(0, seconds))
    return f"{s // 60}:{s % 60:02d}"


# Phrases that mark the *action* ("generate code") rather than the
# *subject* ("static shift") of a code request. Stripping them lets the
# workflow classifier see the real target.
_CODE_ACTION_PHRASES = (
    "generate code for", "write code for", "give me code for",
    "show me code for", "code example for", "sample code for",
    "produce code for", "python script for", "create a notebook for",
    "notebook for", "write a script for", "script for", "code for",
    "generate code", "write code", "give me code", "show me code",
    "python script", "create notebook", "write a script",
    "produce code", "a script to", "script to", "code to",
)


def _code_target_workflow(text: str) -> str | None:
    """Classify the workflow a code request is *about*.

    "generate code for static shift" classifies to ``code_gen`` (the
    action) under the normal keyword table; strip the action phrases so
    the remaining subject ("static shift") routes to ``static_shift``.
    """
    from pycsamt.agents._workflows import classify_workflow
    t = " " + text.lower() + " "
    for p in _CODE_ACTION_PHRASES:
        t = t.replace(p, " ")
    wf = classify_workflow(t)
    return wf if wf and wf != "code_gen" else None


def _drop_workflow(inv_config: dict | None) -> dict:
    """Return *inv_config* without a persisted ``workflow`` key.

    STORE_INV_CONFIG (localStorage) holds inversion *parameters*, but a
    param-modal run also stamps the chosen ``workflow`` into it. Carrying
    that over would force every later request through the same workflow
    (e.g. all requests routed to ai_inversion). Normal/line-picker runs
    must classify the workflow from the request text, so they strip it.
    """
    return {
        k: v for k, v in (inv_config or {}).items()
        if k != "workflow"
    }


def _update_job(
    jid: str, **kw: Any
) -> None:
    with _JOBS_LOCK:
        if jid in _JOBS:
            _JOBS[jid].update(kw)


def _get_job(
    jid: str,
) -> dict | None:
    with _JOBS_LOCK:
        return _JOBS.get(jid)


# ── figure helpers ─────────────────────────────────

def _fig_to_b64(fig: Any) -> str:
    buf = io.BytesIO()
    fig.savefig(
        buf,
        format="png",
        dpi=150,
        bbox_inches="tight",
    )
    buf.seek(0)
    return base64.b64encode(
        buf.read()
    ).decode()


def _fig_thumb_item(
    fig_key: str,
    title: str,
    b64: str,
) -> html.Div:
    """Compact thumbnail tile inside the accordion."""
    short = (
        title if len(title) <= 28
        else title[:26] + "..."
    )
    return html.Div(
        [
            html.Img(
                src=(
                    f"data:image/png;base64,{b64}"
                ),
                className="am-fig-thumb",
                id={
                    "type": "am-fig-img",
                    "key": fig_key,
                },
                title=title,
            ),
            html.Div(
                short,
                className="am-fig-thumb-label",
                title=title,
            ),
            html.Button(
                [
                    html.I(
                        className=(
                            "bi bi-arrows-fullscreen"
                            " me-1"
                        )
                    ),
                    "View",
                ],
                className=(
                    "am-fig-btn am-fig-thumb-btn"
                ),
                id={
                    "type": "am-fig-open",
                    "key": fig_key,
                },
                n_clicks=0,
            ),
        ],
        className="am-fig-thumb-item",
    )


def _fig_accordion(figs: dict) -> html.Div:
    """
    Collapsible accordion for all figures.

    Shows a compact thumbnail grid when expanded.
    Collapsed by default to save space.
    """
    n = len(figs)
    thumbs = [
        _fig_thumb_item(k, v["title"], v["b64"])
        for k, v in figs.items()
    ]
    header = html.Span(
        [
            html.I(
                className=(
                    "bi bi-bar-chart-fill me-2"
                ),
                style={"color": "var(--blue)"},
            ),
            f"Figures ({n})",
            html.Span(
                f"{n}",
                className="am-fig-badge",
            ),
        ],
        className="am-fig-acc-title",
    )
    return html.Div(
        dbc.Accordion(
            dbc.AccordionItem(
                html.Div(
                    thumbs,
                    className="am-fig-grid",
                ),
                title=header,
                item_id="figs",
            ),
            start_collapsed=True,
            flush=True,
            className="am-fig-accordion",
        ),
        className="am-fig-accordion-wrap",
    )


# ── message bubble builders ────────────────────────

def _ts() -> str:
    return datetime.now().strftime("%H:%M")


def _mid() -> str:
    """Stable per-message id used for pinning + scroll-to."""
    return f"am-msg-{uuid.uuid4().hex[:8]}"


def _pin_button(mid: str) -> html.Button:
    """Pin/unpin toggle for a message toolbar."""
    return html.Button(
        html.I(className="bi bi-pin-angle"),
        className="am-msg-action am-pin-btn",
        id={"type": "am-pin-btn", "mid": mid},
        title="Pin message",
        n_clicks=0,
    )


def _user_bubble(text: str, mid: str | None = None) -> html.Div:
    toolbar_btns = [
        html.Button(
            html.I(className="bi bi-clipboard"),
            className="am-msg-action am-copy-btn",
            title="Copy",
            n_clicks=0,
        ),
        html.Button(
            html.I(
                className="bi bi-folder2-open"
            ),
            className="am-msg-action am-edi-msg-btn",
            title="Load EDI",
            n_clicks=0,
        ),
    ]
    if mid:
        toolbar_btns.append(_pin_button(mid))
    return html.Div(
        [
            html.Div(
                html.Div(
                    [
                        html.Div(
                            text,
                            className=(
                                "am-bubble user"
                            ),
                        ),
                        html.Div(
                            _ts(),
                            className="am-ts",
                        ),
                        html.Div(
                            toolbar_btns,
                            className=(
                                "am-msg-toolbar"
                            ),
                        ),
                    ]
                ),
                style={"maxWidth": "100%"},
            ),
            html.Div(
                html.I(
                    className="bi bi-person-fill"
                ),
                className="am-avatar user",
            ),
        ],
        className="am-msg-row user",
        **({"id": mid} if mid else {}),
    )


def _exec_step_row(label: str, status: str) -> html.Div:
    """One row of the executing timeline (rail dot + label)."""
    if status == "done":
        dot = html.Span(
            html.I(className="bi bi-check-lg"),
            className="am-tl-dot",
        )
    elif status == "error":
        dot = html.Span(
            html.I(className="bi bi-exclamation-lg"),
            className="am-tl-dot",
        )
    elif status == "running":
        dot = html.Span(
            html.I(
                className="bi bi-arrow-repeat am-tl-spin"
            ),
            className="am-tl-dot",
        )
    else:  # waiting / pending
        dot = html.Span(className="am-tl-dot")
    return html.Div(
        [dot, html.Span(label, className="am-tl-lbl")],
        className=f"am-tl-step {status}",
    )


def _thinking_bubble(
    steps: list[dict],
    workflow: str | None = None,
    elapsed: float | None = None,
) -> html.Div:
    steps = steps or []
    n_total = len(steps)
    n_done = sum(
        1 for s in steps if s.get("status") == "done"
    )

    # Header: "Running <workflow>" + live elapsed.
    name = (
        _WF_RUNNING_LABEL.get(workflow, workflow)
        if workflow
        else None
    )
    title = (
        [html.Span("Running "), html.Span(
            name, className="am-exec-name")]
        if name
        else [html.Span("Working")]
    )
    header_children = [
        html.I(
            className="bi bi-arrow-repeat am-tl-spin"
            " am-exec-spin"
        ),
        html.Span(title, className="am-exec-title"),
    ]
    if elapsed is not None:
        header_children.append(
            html.Span(
                _fmt_elapsed(elapsed),
                className="am-exec-elapsed",
            )
        )
    header = html.Div(
        header_children, className="am-exec-header"
    )

    # Indeterminate progress sweep.
    bar = html.Div(
        html.I(), className="am-exec-bar"
    )

    # Step timeline.
    timeline = html.Div(
        [
            _exec_step_row(s["label"], s["status"])
            for s in steps
        ],
        className="am-tl",
    ) if steps else html.Div()

    footer = (
        html.Div(
            f"Step {min(n_done + 1, n_total)} of {n_total}",
            className="am-exec-footer",
        )
        if n_total
        else html.Div()
    )

    return html.Div(
        [
            html.Div(
                html.I(className="bi bi-robot"),
                className="am-avatar agent",
            ),
            html.Div(
                [header, bar, timeline, footer],
                className=(
                    "am-bubble agent am-exec-bubble"
                ),
            ),
        ],
        className="am-msg-row",
        id="am-thinking-bubble",
    )


def _strip_thinking(msgs: list | None) -> list:
    """Return *msgs* without the live thinking bubble.

    Used when a running task is stopped: the in-progress bubble is
    removed so a final notice can take its place.
    """
    out = []
    for child in (msgs or []):
        if (
            isinstance(child, dict)
            and child.get("props", {}).get("id")
            == "am-thinking-bubble"
        ):
            continue
        out.append(child)
    return out


def _stop_job_response(
    current_msgs: list | None,
    stored_messages: list | None,
    job_store: dict | None,
):
    """Cancel the active job and build the send_message return tuple.

    Marks the job cancelled, drops the live thinking bubble, appends a
    "Task stopped" notice, clears the job store and disables polling so
    a late background result is never displayed. The user's typed input
    is preserved (``no_update``).
    """
    jid = (job_store or {}).get("jid")
    if jid:
        _update_job(jid, status="cancelled")
    msgs = _strip_thinking(current_msgs)
    msgs.append(_agent_bubble("Task stopped.", kind=KIND_ERROR))
    new_stored = list(stored_messages or [])
    new_stored.append({
        "role": "assistant",
        "content": "Task stopped by user.",
        "ts": _ts(),
    })
    return msgs, {}, True, no_update, new_stored, {}


# Pure markdown parsing lives in a dash-free module so it can be
# unit-tested without importing Dash / the GUI package.
from .._markdown import (
    split_inline_bold as _split_inline_bold,
    parse_markdown as _parse_markdown,
)


def _code_block(code: str) -> html.Div:
    """Render a collapsible, syntax-highlighted Python code block.

    Uses a native ``<details>`` element so the user can toggle the code
    open/closed without a callback. The Copy button sits outside the
    ``<summary>`` so clicking it copies (and flashes "Copied") instead of
    toggling the accordion.
    """
    import uuid as _uuid
    copy_id = f"am-copy-{_uuid.uuid4().hex[:8]}"
    n_lines = code.count("\n") + 1
    return html.Div(
        [
            html.Details(
                [
                    html.Summary(
                        [
                            html.I(
                                className=(
                                    "bi bi-chevron-right"
                                    " am-code-chevron"
                                )
                            ),
                            html.I(
                                className=(
                                    "bi bi-code-slash me-1"
                                ),
                                style={"color": "#61afef"},
                            ),
                            html.Span(
                                "python",
                                className="am-code-lang",
                            ),
                            html.Span(
                                f"{n_lines} lines",
                                className="am-code-meta",
                            ),
                        ],
                        className="am-code-summary",
                    ),
                    # code body — hljs highlights this
                    html.Pre(
                        html.Code(
                            code,
                            className="language-python",
                        ),
                    ),
                ],
                open=True,
                className="am-code-details",
            ),
            html.Button(
                [
                    html.I(
                        className="bi bi-clipboard me-1"
                    ),
                    "Copy",
                ],
                id=copy_id,
                className="am-code-copy-btn",
                title="Copy to clipboard",
                **{"data-code": code},
                n_clicks=0,
            ),
        ],
        className="am-code-block",
    )


def _render_inline(text: str) -> list:
    """Render inline ``**bold**`` markup into html spans."""
    return [
        html.Strong(chunk) if is_bold else html.Span(chunk)
        for is_bold, chunk in _split_inline_bold(text)
    ]


def _render_markdown(text: str) -> list:
    """Render lightweight markdown (headings, bullets, code, bold)."""
    children: list = []
    for tok in _parse_markdown(text):
        if tok[0] == "code":
            children.append(_code_block(tok[2]))
        elif tok[0] == "heading":
            children.append(
                html.P(
                    html.Strong(tok[1]),
                    className="am-md-h",
                )
            )
        elif tok[0] == "bullet":
            children.append(
                html.Li(
                    _render_inline(tok[1]),
                    style={"marginLeft": "12px"},
                )
            )
        elif tok[0] == "para":
            children.append(html.P(_render_inline(tok[1])))
        # "blank" tokens are skipped (paragraph spacing is via CSS)
    return children


def _kind_header(kind: str | None) -> html.Div | None:
    """Build the small per-kind label chip, or None for plain replies."""
    spec = _KIND_HEADER.get(kind or "")
    if not spec:
        return None
    icon, label, colour = spec
    return html.Div(
        [
            html.I(
                className=f"bi {icon} me-2",
                style={"color": colour},
            ),
            html.Span(label),
        ],
        className=f"am-kind-header am-kind-{kind}",
    )


def _agent_bubble(
    text: str,
    steps: list[dict] | None = None,
    figs: dict | None = None,
    code: str = "",
    kind: str | None = None,
    mid: str | None = None,
) -> html.Div:
    children: list = []

    # per-kind label chip (answer / code / clarify / …)
    header = _kind_header(kind)
    if header is not None:
        children.append(header)

    # body text — full lightweight-markdown rendering
    children.extend(_render_markdown(text))

    # step summary
    if steps:
        step_divs = [
            html.Div(
                [
                    html.I(
                        className=(
                            "bi bi-check-circle-fill"
                            " am-step-icon"
                        )
                    ),
                    html.Span(s["label"]),
                ],
                className="am-step done",
            )
            for s in steps
        ]
        children.append(
            html.Div(
                step_divs,
                className="am-steps",
            )
        )

    # generated code block
    if code and code.strip():
        children.append(_code_block(code))

    # figures go to sidebar only; show a compact
    # chip in the chat bubble so the user knows
    # results are ready without a full preview.
    if figs:
        n = len(figs)
        children.append(
            html.Div(
                [
                    html.I(
                        className=(
                            "bi bi-bar-chart-fill"
                            " me-2"
                        ),
                        style={
                            "color": "var(--blue)"
                        },
                    ),
                    html.Span(
                        f"{n} figure"
                        f"{'s' if n != 1 else ''}"
                        " generated — open the"
                        " Figures panel to view.",
                        style={
                            "fontSize": "12px",
                            "color": "var(--sub1)",
                        },
                    ),
                ],
                className="am-fig-notice",
            )
        )

    return html.Div(
        [
            html.Div(
                html.I(className="bi bi-robot"),
                className="am-avatar agent",
            ),
            html.Div(
                [
                    html.Div(children),
                    html.Div(
                        _ts(), className="am-ts"
                    ),
                    html.Div(
                        [
                            html.Button(
                                html.I(
                                    className="bi bi-clipboard"
                                ),
                                className=(
                                    "am-msg-action am-copy-btn"
                                ),
                                title="Copy answer",
                                n_clicks=0,
                                **{"data-copy": text or ""},
                            ),
                            _pin_button(mid),
                        ],
                        className="am-msg-toolbar",
                    ) if mid else html.Div(),
                ],
                className=(
                    "am-bubble agent"
                    + (f" am-bubble-{kind}" if kind else "")
                ),
            ),
        ],
        className="am-msg-row",
        **({"id": mid} if mid else {}),
    )


# ── Workflow figure filter ────────────────────────
# Maps workflow type → step names whose figures
# should appear in the agent bubble response.
# Steps NOT in this set are prerequisite steps
# (load, qc, denoise) whose figures are suppressed
# unless the user explicitly asked for them.
# None means "show figures from all steps."
_WORKFLOW_FIGURE_STEPS: dict[
    str, set[str] | None
] = {
    "qc":               {"qc", "static_shift", "report"},
    "static_shift":     {"static_shift"},
    "phase_analysis":   {"phase_analysis", "report"},
    "ai_inversion":     {"ai_inv", "interpret", "report"},
    "inv1d":            {"ai_inv", "interpret", "report"},
    "inv2d":            {"inv2d", "interpret", "report"},
    "inv3d":            {"inv3d", "interpret", "report"},
    "ensemble_inversion": {
        "ensemble", "interpret", "report",
    },
    "pinn_inversion":   {
        "pinn_inv", "interpret", "report",
    },
    "hybrid_inversion": {
        "hybrid_inv", "interpret", "report",
    },
    "joint_inversion":  {
        "joint", "interpret", "report",
    },
    "tipper":           {"tipper", "report"},
    "modem":            {"modem", "report"},
    "occam2d":          {"occam2d", "report"},
    "pre_inversion":    {
        "phase_analysis", "occam2d", "report",
    },
    "sensitivity":      {"sensitivity", "report"},
    "rotation":         {"rotate"},
    "freq_decimation":  {"decimate", "report"},
    "batch":            {"batch", "report"},
    "comparison":       {"compare", "report"},
    "full":             None,  # show everything
}
# Default set for unlisted workflows: skip
# only load and denoise steps.
_SKIP_ALWAYS: frozenset[str] = frozenset(
    {"load", "denoise"}
)


# ── Smart param detection ─────────────────────────

# Workflows that require user parameters
# before the job can start.
# Workflows that REQUIRE user parameters before
# starting.  Simple data-processing workflows
# (qc, static_shift, phase_analysis, tipper,
# rotation) run with sensible defaults and do
# NOT need a param modal.  Only inversion and
# a few analytical workflows need user input.
_NEEDS_PARAMS: frozenset[str] = frozenset({
    # Inversion workflows (model params needed)
    "ai_inversion", "inv1d",
    "inv2d", "inv3d",
    "pinn_inversion", "hybrid_inversion",
    "ensemble_inversion",
    "pre_inversion", "modem",
    # Analysis workflows (optional params)
    "denoise",
    "sensitivity",
    "interpret", "interpretation",
    "report",
    "code_gen",
    # Plotting tasks (station / component / period / publication)
    "rhophi", "phase_psection", "pt_psection", "tipper_plot",
    "phase_tensor_map", "station_response", "strike_profile",
    "strike", "dimensionality", "validator",
    # Data / IO tools
    "coords", "elevation", "converter", "batch_export",
})

_WF_LABELS: dict[str, str] = {
    "ai_inversion":       "1-D AI inversion",
    "inv1d":              "1-D AI inversion",
    "inv2d":              "2-D U-Net inversion",
    "inv3d":              "3-D GCN inversion",
    "pinn_inversion":     "PINN inversion",
    "hybrid_inversion":   "hybrid inversion",
    "ensemble_inversion": "ensemble inversion",
    "pre_inversion":      "pre-inversion setup",
    "modem":              "ModEM preparation",
    "qc":                 "QC pipeline",
    "phase_analysis":     "phase tensor analysis",
    "static_shift":    "static shift correction",
    "tipper":             "tipper analysis",
    "rotation":           "data rotation",
    "interpret":     "geological interpretation",
    "interpretation": "geological interpretation",
    "report":             "survey report",
    "code_gen":           "code generation",
    "denoise":            "data denoising",
    "sensitivity":    "sensitivity / DOI analysis",
    "rhophi":             "rho/phi sounding curves",
    "phase_psection":     "phase pseudo-section",
    "pt_psection":        "phase-tensor pseudo-section",
    "tipper_plot":        "tipper plot",
    "phase_tensor_map":   "phase-tensor map",
    "station_response":   "station response",
    "strike_profile":     "strike profile",
    "strike":             "strike analyzer",
    "dimensionality":     "dimensionality classifier",
    "validator":          "EDI validator",
    "coords":             "coordinate transformer",
    "elevation":          "elevation enrichment",
    "converter":          "format converter",
    "batch_export":       "batch plot export",
}


def _quick_workflow(text: str) -> str | None:
    """
    Fast regex classification.

    Returns workflow type string if user
    parameters are needed before the job
    can start, else None.
    """
    try:
        from pycsamt.agents.context import (
            _regex_extract,
            _normalise_config,
        )
        cfg = _regex_extract(text)
        cfg = _normalise_config(cfg, text)
        wf = cfg.get("workflow", "")
        return wf if wf in _NEEDS_PARAMS else None
    except Exception:
        return None


def _waiting_bubble(wf: str) -> html.Div:
    label = _WF_LABELS.get(
        wf, wf.replace("_", " ")
    )
    return html.Div(
        [
            html.Div(
                html.I(className="bi bi-robot"),
                className="am-avatar agent",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.I(
                                className=(
                                    "bi bi-sliders2"
                                    " me-2"
                                ),
                                style={
                                    "color": (
                                        "var(--blue)"
                                    )
                                },
                            ),
                            html.Span(
                                f"To run {label},"
                                " I need a few"
                                " parameters."
                            ),
                            html.Span(
                                " Please fill in"
                                " the form above"
                                " and click"
                                " Run Workflow.",
                                style={
                                    "opacity": ".7"
                                },
                            ),
                        ]
                    ),
                    html.Div(
                        _ts(), className="am-ts"
                    ),
                ],
                className="am-bubble agent",
            ),
        ],
        className="am-msg-row",
        id="am-waiting-bubble",
    )


def _line_waiting_bubble() -> html.Div:
    """Bubble shown while the line picker is open."""
    return html.Div(
        [
            html.Div(
                html.I(className="bi bi-robot"),
                className="am-avatar agent",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.I(
                                className=(
                                    "bi bi-layers"
                                    " me-2"
                                ),
                                style={
                                    "color": (
                                        "var(--blue)"
                                    )
                                },
                            ),
                            html.Span(
                                "Multiple survey"
                                " lines are loaded."
                                " Please select"
                                " which line(s) to"
                                " process using"
                                " the panel above.",
                            ),
                        ]
                    ),
                    html.Div(
                        _ts(), className="am-ts"
                    ),
                ],
                className="am-bubble agent",
            ),
        ],
        className="am-msg-row",
        id="am-line-waiting-bubble",
    )


def _extract_line_ref(
    text: str, groups: dict
) -> str | None:
    """Detect a line/profile reference in text.

    Checks known group names first (direct match),
    then 'line N' / 'profile N' patterns, then
    ordinal words ('first line', etc.).
    Returns the matched token or None.
    """
    t = text.lower()
    for key in groups:
        if key.lower() in t:
            return key
    m = re.search(
        r"\b(?:line|profile)\s+(\S+)", t
    )
    if m:
        return m.group(1)
    _ordinals = {
        "first": "1", "second": "2",
        "third": "3", "fourth": "4",
        "fifth": "5", "sixth": "6",
        "seventh": "7", "eighth": "8",
        "ninth": "9", "tenth": "10",
    }
    for word, num in _ordinals.items():
        if (
            f"{word} line" in t
            or f"{word} profile" in t
        ):
            return num
    return None


def _match_group(
    ref: str, groups: dict
) -> str | None:
    """Return the group key matching ref exactly
    or case-insensitively. No ordinal resolution —
    numeric refs always trigger the line picker.
    """
    if ref in groups:
        return ref
    ref_l = ref.lower()
    for key in groups:
        if key.lower() == ref_l:
            return key
    return None


# ── PINN / Hybrid keyword detection ───────────────

_PINN_KWS = frozenset({
    "pinn", "physics-informed",
    "physics informed", "no training data",
    "gradient descent",
})
_HYBRID_KWS = frozenset({
    "hybrid", "two-stage", "two stage",
    "ai + physics", "warm start", "warmstart",
})


def _is_pinn_or_hybrid(text: str) -> bool:
    t = text.lower()
    return (
        any(kw in t for kw in _PINN_KWS)
        or any(kw in t for kw in _HYBRID_KWS)
    )


# ── Web App launcher ──────────────────────────────

_WEB_APP_KWS: frozenset[str] = frozenset({
    "open web app",
    "launch web app",
    "open the web",
    "web interface",
    "full interface",
    "full app",
    "go to web app",
    "pycsamt web",
    "web application",
    "launch app",
    "full web",
    "open web",
    "start web app",
    "open full app",
})

# Tasks too complex for chat → redirect
_COMPLEX_VIZ_KWS: frozenset[str] = frozenset({
    "3d map",
    "station map",
    "explore results",
    "browse data",
    "browse edis",
    "pseudosection viewer",
    "phase tensor map",
    "full visualization",
    "map view",
    "3d visualization",
    "interactive plot",
    "interactive map",
    "interactive pseudosection",
    "full pipeline editor",
})

_webapp_state: dict = {
    "url": None,
    "running": False,
}
_webapp_lock = threading.Lock()


def _is_web_app_request(text: str) -> bool:
    t = text.lower()
    return any(kw in t for kw in _WEB_APP_KWS)


def _is_complex_viz(text: str) -> bool:
    t = text.lower()
    return any(
        kw in t for kw in _COMPLEX_VIZ_KWS
    )


def _free_port(preferred: int = 8051) -> int:
    import socket as _socket
    try:
        with _socket.socket(
            _socket.AF_INET,
            _socket.SOCK_STREAM,
        ) as s:
            s.setsockopt(
                _socket.SOL_SOCKET,
                _socket.SO_REUSEADDR,
                1,
            )
            s.bind(("127.0.0.1", preferred))
        return preferred
    except OSError:
        with _socket.socket(
            _socket.AF_INET,
            _socket.SOCK_STREAM,
        ) as s:
            s.bind(("127.0.0.1", 0))
            return int(s.getsockname()[1])


def _ensure_web_app() -> str:
    """
    Start the full pyCSAMT web app in a
    background thread if not already running.
    Returns the URL string.
    """
    global _webapp_state
    with _webapp_lock:
        if _webapp_state["running"]:
            return _webapp_state["url"]

        port = _free_port(8051)
        url = f"http://127.0.0.1:{port}"

        def _run_webapp():
            try:
                from pycsamt.app.web.app import (
                    create_app,
                )
                wa = create_app(debug=False)
                wa.run(
                    host="127.0.0.1",
                    port=port,
                    debug=False,
                    use_reloader=False,
                )
            except Exception:
                pass

        t = threading.Thread(
            target=_run_webapp,
            daemon=True,
        )
        t.start()
        _webapp_state["url"] = url
        _webapp_state["running"] = True
        return url


def _web_app_bubble(
    url: str,
    reason: str = "",
) -> html.Div:
    desc = reason or (
        "The full pyCSAMT web application"
        " provides interactive station maps,"
        " pseudosection viewers, phase-tensor"
        " tools, pipeline editor, inversion"
        " pages, and more — beyond what the"
        " agent chat can display."
    )
    return html.Div(
        [
            html.Div(
                html.I(className="bi bi-robot"),
                className="am-avatar agent",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.I(
                                        className=(
                                            "bi bi-window"
                                            "-fullscreen"
                                            " me-2"
                                        ),
                                        style={
                                            "color": (
                                                "var(--green)"
                                            )
                                        },
                                    ),
                                    html.Strong(
                                        "Launching"
                                        " pyCSAMT"
                                        " Web App"
                                    ),
                                ],
                                className=(
                                    "am-webapp-hdr"
                                ),
                            ),
                            html.P(
                                desc,
                                className=(
                                    "am-webapp-desc"
                                ),
                            ),
                            html.A(
                                [
                                    html.I(
                                        className=(
                                            "bi bi-box"
                                            "-arrow-up"
                                            "-right me-2"
                                        )
                                    ),
                                    url,
                                ],
                                href=url,
                                target="_blank",
                                rel="noopener",
                                className=(
                                    "am-webapp-link"
                                ),
                            ),
                            html.Div(
                                [
                                    html.I(
                                        className=(
                                            "bi bi-info"
                                            "-circle me-1"
                                        )
                                    ),
                                    "Server is"
                                    " starting —"
                                    " the link"
                                    " will be ready"
                                    " in a few"
                                    " seconds.",
                                ],
                                className=(
                                    "am-webapp-note"
                                ),
                            ),
                        ],
                        className="am-webapp-card",
                    ),
                    html.Div(
                        _ts(), className="am-ts"
                    ),
                ],
                className="am-bubble agent",
            ),
        ],
        className="am-msg-row",
    )


# ── intent dispatch helpers ────────────────────────

def _api_key_hint() -> str:
    """Reusable note nudging the user to set an API key for online mode."""
    return (
        "**Want richer, more natural answers?** Open **Settings** (top-right)"
        " and add an API key for any one provider — **Anthropic (Claude)**,"
        " **OpenAI**, **Google Gemini**, or **DeepSeek**. With a key I switch"
        " to online mode, which most improves:\n"
        "- **Questions** about pyCSAMT — fluent, synthesised answers"
        " (grounded in the package via RAG) instead of the offline summary;\n"
        "- **Code generation** — complete, tailored scripts;\n"
        "- **Request understanding** — better routing of free-form requests.\n"
        "Running workflows works offline too — no key required."
    )


def _capability_text() -> str:
    """Static capability summary for META / greeting / 'list tasks' intents."""
    return (
        "I'm the **pyCSAMT v2 assistant**. Here's what I can do for you:\n\n"
        "**Run processing workflows** on your loaded EDI data:\n"
        "- `qc` — quality control & per-station scan\n"
        "- `static_shift` — static-shift detection & AMA/LOESS correction\n"
        "- `denoise` — RPCA / Hampel / AI denoising\n"
        "- `phase_analysis` — phase-tensor, strike & dimensionality\n"
        "- `tipper` — tipper / induction arrows\n"
        "- `rotation` — tensor rotation to a strike frame\n"
        "- `freq_decimation` — period selection / decimation\n"
        "- `sensitivity` — sensitivity kernels / depth of investigation\n"
        "- `pre_inversion` — full pre-inversion preparation\n"
        "- `ai_inversion`, `inv2d`, `inv3d`, `pinn_inversion`,"
        " `hybrid_inversion`, `ensemble_inversion`, `joint_inversion` —"
        " inversions\n"
        "- `occam2d`, `modem` — inversion input preparation\n"
        "- `report` — generate a survey report\n\n"
        "**Make plots** from your data (I'll ask for stations, components,"
        " period range and publication styling):\n"
        "- rho/φ sounding curves\n"
        "- phase pseudo-section\n"
        "- phase-tensor (Φ) ellipse pseudo-section\n"
        "- phase-tensor map (ellipses at a chosen period)\n"
        "- station response inspector (Bode ρa/φ curves)\n"
        "- strike profile (strike vs station position)\n"
        "- tipper components / induction arrows (when tipper data exists)\n\n"
        "**Analyze** your data (results as a table + figure):\n"
        "- strike analyzer (geoelectric strike per station)\n"
        "- dimensionality classifier (1-D / 2-D / 3-D)\n"
        "- EDI validator (per-station quality checklist)\n\n"
        "**Data & I/O tools** (I'll ask for the options first):\n"
        "- coordinate transformer (station lat/lon → UTM)\n"
        "- elevation enrichment (fetch elevation from an open web API)\n"
        "- format converter (re-export the survey to CSV / JSON / EDI)\n"
        "- batch plot export (render & save a bundle of figures)\n\n"
        "**Answer questions** about pyCSAMT — classes, functions, the Sites"
        " data model, and which method to use.\n\n"
        "**Generate Python code** that reproduces a pyCSAMT workflow.\n\n"
        "To run a workflow, load an EDI dataset first with **Load EDI**"
        " (top-left). Questions and code requests need no data.\n\n"
        "For interactive maps, pseudosection viewers and the visual pipeline"
        " editor, use the full **pyCSAMT web application**.\n\n"
        + _api_key_hint()
    )


def _unknown_task_text(text: str) -> str:
    """Friendly 'I don't recognise that task' reply — used instead of
    silently defaulting an unrecognised request to the QC workflow."""
    snippet = (text or "").strip().replace("\n", " ")
    if len(snippet) > 80:
        snippet = snippet[:77] + "…"
    quoted = f' — "{snippet}"' if snippet else ""
    return (
        f"I'm not sure how to handle that as a task{quoted}.\n\n"
        "I didn't recognise a workflow I can run for that request, so I"
        " won't guess. Here's what I can actually do:\n"
        "- **Run workflows** on loaded EDI data: QC, static-shift,"
        " phase-tensor analysis, denoising, tipper, rotation, frequency"
        " decimation, sensitivity/DOI, inversions (AI 1-D/2-D/3-D, PINN,"
        " hybrid, ensemble, joint), ModEM/Occam2D prep, and reports.\n"
        "- **Answer questions** about pyCSAMT (classes, functions, methods).\n"
        "- **Generate Python code** for a pyCSAMT workflow.\n"
        "- **List my capabilities** — just ask \"what can you do?\".\n\n"
        "If you need maps, interactive pseudosection viewers, or a feature"
        " that isn't here, try the full **pyCSAMT web application**.\n\n"
        "Tip: phrase it as an action + target, e.g. \"run static shift\","
        " \"denoise the data\", or \"run AI inversion\".\n\n"
        + _api_key_hint()
    )


def _dispatch_question(
    jid: str,
    text: str,
    *,
    llm_prov: str,
    api_key: str | None,
    sel_model: str | None,
    offline: bool,
    history: list[dict] | None,
    step,
) -> None:
    """Answer a package question via PackageQAAgent."""
    from pycsamt.agents.package_qa import PackageQAAgent
    from pycsamt.api.agents import AGENT_CONFIG

    step("Answering question...", "running")
    ctx_str = ""
    if history:
        recent = [
            f"{m.get('role', 'user')}: {m.get('content', '')}"
            for m in history[-4:]
        ]
        ctx_str = "\n".join(recent)

    with (
        AGENT_CONFIG.offline() if offline else _nullctx()
    ):
        qa = PackageQAAgent(
            llm_provider=llm_prov,
            api_key=api_key,
            model=sel_model,
        )
        res = qa.execute(
            {"question": text, "context": ctx_str}
        )

    answer = (
        res.get("answer")
        or res.summary
        or "I couldn't find an answer in the pyCSAMT reference."
    )
    # When running without a key, the answer is the deterministic offline
    # composition — nudge the user that an API key unlocks a fluent reply.
    if offline or (res.get("source") == "rag_offline"):
        answer = (
            answer
            + "\n\n---\n*Offline answer composed from the pyCSAMT reference."
            " For a fuller, synthesised response, add an API key (Claude,"
            " OpenAI, Gemini or DeepSeek) in **Settings**.*"
        )
    step("Answer ready", "done")
    _update_job(
        jid,
        status="done",
        result=answer,
        steps=_JOBS[jid]["steps"],
        kind=KIND_ANSWER,
    )


def _dispatch_plot(
    jid: str,
    edi_path: Any,
    *,
    kind: str,
    params: dict,
    step,
    label: str = "",
) -> None:
    """Render a plotting task (rho/phi, phase / phase-tensor pseudo-section,
    or tipper) via the PlotAgent and publish the figure(s) to the chat +
    Figures panel. ``label`` names the dataset/line for user-facing messages."""
    import matplotlib.pyplot as plt
    from pycsamt.agents.plotting import PlotAgent

    _labels = {
        "rhophi": "rho/phi sounding curves",
        "phase_psection": "phase pseudo-section",
        "pt_psection": "phase-tensor pseudo-section",
        "tipper": "tipper plot",
    }
    where = f" for {label}" if label else ""
    step(f"Rendering {_labels.get(kind, kind)}...", "running")

    agent_input = {"path": edi_path, "kind": kind, **(params or {})}
    res = PlotAgent().execute(agent_input)

    if res.status == "failed":
        step("Plot unavailable", "done")
        # Tipper is frequently absent in CSAMT/AMT data — answer plainly.
        if res.get("reason") == "no_tipper":
            msg = (
                f"Tipper data is not available{where}. This dataset has no "
                "vertical magnetic transfer function (Tx/Ty), so induction "
                "arrows and tipper components can't be plotted. The other "
                "plots (rho/φ curves, phase and phase-tensor pseudo-sections) "
                "work fine — want one of those instead?"
            )
        else:
            msg = res.summary + (
                "\n\nHint: " + (res.get("hint") or "")
                if res.get("hint") else ""
            )
        _update_job(
            jid,
            status="done",
            result=msg,
            steps=_JOBS[jid]["steps"],
            kind=KIND_META if res.get("reason") == "no_tipper" else KIND_ERROR,
        )
        return

    figs: dict = {}
    for fname, fig in (res.data.get("figures") or {}).items():
        if isinstance(fig, plt.Figure):
            figs[str(uuid.uuid4())] = {
                "title": fname,
                "b64": _fig_to_b64(fig),
            }
            plt.close(fig)

    summary = res.summary
    if res.warnings:
        summary += "\n\n" + "\n".join(f"⚠ {w}" for w in res.warnings[:3])
    if not figs:
        summary += "\n\n(No figure was produced — check the parameters.)"

    _record_run(
        workflow=kind,
        path=str(edi_path),
        output_dir="",
        status=res.status,
        summary=res.summary,
        n_figures=len(figs),
    )
    step("Figure ready", "done")
    _update_job(
        jid,
        status="done",
        result=summary,
        steps=_JOBS[jid]["steps"],
        figs=figs,
        kind=KIND_WORKFLOW,
    )


def _dispatch_tool(
    jid: str,
    edi_path: Any,
    *,
    kind: str,
    params: dict,
    step,
    label: str = "",
) -> None:
    """Run an analysis tool (strike / dimensionality / validator) via the
    ToolAgent and publish its table + optional figure to the chat."""
    import matplotlib.pyplot as plt
    from pycsamt.agents.tooling import ToolAgent

    _labels = {
        "strike": "strike analysis",
        "dimensionality": "dimensionality classification",
        "validator": "EDI validation",
        "coords": "coordinate transform",
        "elevation": "elevation enrichment",
        "converter": "format conversion",
        "batch_export": "batch plot export",
    }
    where = f" for {label}" if label else ""
    step(f"Running {_labels.get(kind, kind)}...", "running")

    res = ToolAgent().execute(
        {"path": edi_path, "kind": kind, **(params or {})}
    )
    if res.status == "failed":
        step("Analysis failed", "done")
        _update_job(
            jid, status="done",
            result=res.summary + ("." if not res.summary.endswith(".") else "")
            + f" (dataset{where or ''})",
            steps=_JOBS[jid]["steps"], kind=KIND_ERROR,
        )
        return

    figs: dict = {}
    for fname, fig in (res.data.get("figures") or {}).items():
        if isinstance(fig, plt.Figure):
            figs[str(uuid.uuid4())] = {
                "title": fname, "b64": _fig_to_b64(fig),
            }
            plt.close(fig)

    table = (res.data or {}).get("table_text", "")
    result = res.summary
    if table:
        # a fenced block renders monospaced (the chat markdown has no tables)
        result += "\n\n```\n" + table + "\n```"
    if res.warnings:
        result += "\n\n" + "\n".join(f"⚠ {w}" for w in res.warnings[:3])

    _record_run(
        workflow=kind, path=str(edi_path), output_dir="",
        status=res.status, summary=res.summary, n_figures=len(figs),
    )
    step("Done", "done")
    _update_job(
        jid, status="done", result=result, steps=_JOBS[jid]["steps"],
        figs=figs, kind=KIND_WORKFLOW if figs else KIND_ANSWER,
    )


def _dispatch_code(
    jid: str,
    text: str,
    edi_store: dict,
    settings: dict,
    *,
    workflow: str | None,
    llm_prov: str,
    api_key: str | None,
    sel_model: str | None,
    offline: bool,
    step,
) -> None:
    """Generate a standalone pyCSAMT script via CodeGenerationAgent."""
    from pycsamt.agents.context import ContextInputAgent
    from pycsamt.agents.code_gen import CodeGenerationAgent
    from pycsamt.api.agents import AGENT_CONFIG

    step("Extracting configuration...", "done")
    with (
        AGENT_CONFIG.offline() if offline else _nullctx()
    ):
        ctx_agent = ContextInputAgent(
            llm_provider=llm_prov,
            api_key=api_key,
            model=sel_model,
        )
        ctx_res = ctx_agent.execute({"request": text})
    cfg = (
        ctx_res.data.get("config", {})
        if ctx_res and ctx_res.data
        else {}
    )
    # Pick the workflow the code should be ABOUT (the subject), not the
    # "code_gen" action. Priority: router slot (if specific) → subject
    # extracted from the text → sensible default.
    target = (
        workflow
        if (workflow and workflow != "code_gen")
        else _code_target_workflow(text)
    )
    if target:
        cfg["workflow"] = target
    elif cfg.get("workflow", "") in ("", "code_gen"):
        # plain "write me a script" with no subject → default to qc
        cfg["workflow"] = "qc"

    # ── RAG grounding: resolve a named survey line to its real path and
    # retrieve real symbols/recipe so the generated code is accurate. ──
    rag_text = ""
    resolved_line = None
    try:
        from pycsamt.assistant.rag.context_builder import (
            default_context_builder,
        )
        builder = default_context_builder()
        if builder is not None:
            ac = builder.build(text)
            rag_text = ac.context_text
            pc = ac.project_context
            if pc.get("exists") and pc.get("edi_dir"):
                resolved_line = pc.get("line")
                cfg["data_path"] = pc["edi_dir"]
    except Exception:  # noqa: BLE001 — RAG is best-effort
        rag_text = ""

    # Use a loaded EDI path when present (and no line was resolved) so the
    # script is immediately runnable; otherwise code_gen inserts a
    # /path/to/EDIs placeholder.
    if not resolved_line:
        edi_path = (edi_store or {}).get("path", "") or cfg.get(
            "data_path", ""
        )
        if edi_path:
            cfg["data_path"] = edi_path

    output_dir = (
        settings.get("output_dir") or ""
    ).strip() or "pycsamt_workflow_output"

    _update_job(jid, workflow="code_gen")
    step("Generating code...", "running")
    with (
        AGENT_CONFIG.offline() if offline else _nullctx()
    ):
        cg = CodeGenerationAgent(
            llm_provider=llm_prov,
            api_key=api_key,
            model=sel_model,
        )
        res = cg.execute(
            {
                "workflow_config": cfg,
                "results": {},
                "output_dir": output_dir,
                "rag_context": rag_text,
            }
        )

    code = res.get("code", "") if res else ""

    # Validate the generated script (deterministic): catch syntax errors
    # and any hallucinated pyCSAMT symbols before the user runs it.
    _valid_note = ""
    try:
        from pycsamt.assistant.tools.validation_tools import (
            validate_generated_code,
        )
        rep = validate_generated_code(code)
        if rep["ok"]:
            _valid_note = (
                "\n\n✓ Validated: syntax OK and all pyCSAMT imports"
                " resolve to real symbols."
            )
        elif not rep["syntax_ok"]:
            _valid_note = (
                "\n\n⚠ Validation: the script has a syntax error — "
                + "; ".join(rep["errors"][:2])
            )
        else:
            _valid_note = (
                "\n\n⚠ Validation: some symbols could not be verified — "
                + "; ".join(rep["errors"][:3])
            )
    except Exception:  # noqa: BLE001 — validation is best-effort
        _valid_note = ""

    _line_note = (
        f" for line {resolved_line}" if resolved_line else ""
    )
    summary = (
        "Here is a standalone pyCSAMT script that"
        f" reproduces the {cfg.get('workflow', 'qc')}"
        f" workflow{_line_note}. Copy it from the code block"
        " below — edit the data path if needed."
        + _valid_note
    )
    step("Code ready", "done")
    _update_job(
        jid,
        status="done",
        result=summary,
        code=code,
        steps=_JOBS[jid]["steps"],
        kind=KIND_CODE,
    )


# ── agent runner ───────────────────────────────────

def _run_agent(
    jid: str,
    text: str,
    edi_store: dict,
    settings: dict,
    inv_config: dict | None = None,
    history: list[dict] | None = None,
) -> None:
    """
    Execute in a background thread.
    Writes steps + result to _JOBS[jid].
    """
    def _step(label: str, status: str = "done"):
        with _JOBS_LOCK:
            _JOBS[jid]["steps"].append(
                {"label": label, "status": status}
            )

    try:
        _step("Parsing request...", "done")

        # configure provider + api key
        provider = settings.get(
            "provider", "offline"
        )
        key_map = {
            "claude": "ANTHROPIC_API_KEY",
            "openai": "OPENAI_API_KEY",
            "gemini": "GOOGLE_API_KEY",
            "deepseek": "DEEPSEEK_API_KEY",
        }
        api_key: str | None = None
        if provider in key_map:
            api_key = (
                settings.get(
                    f"key_{provider}", ""
                ) or None
            )
            if api_key:
                import os
                os.environ[
                    key_map[provider]
                ] = api_key

        # "offline" maps to "claude" provider
        # name so BaseAgent validates it, but
        # api_key stays None → regex fallback.
        llm_prov = (
            provider
            if provider != "offline"
            else "claude"
        )
        sel_model = settings.get(
            f"model_{llm_prov}"
        ) or None

        # import agents lazily
        from pycsamt.agents.context import (
            ContextInputAgent,
        )
        from pycsamt.agents.orchestrator import (
            WorkflowOrchestratorAgent,
        )
        from pycsamt.api.agents import AGENT_CONFIG

        # When truly offline, use AGENT_CONFIG.offline()
        # around every agent creation + execution so
        # that env-based keys (e.g. ANTHROPIC_API_KEY
        # from .env.local) are never picked up.
        # _nullctx is used for online providers so the
        # same code path works for both cases.
        _offline = provider == "offline"

        # ── top-level intent routing ──────────────
        # Decide WHAT KIND of request this is before
        # assuming it is a workflow to execute. This
        # is the master dispatch: questions, code, and
        # capability requests never touch the workflow
        # pipeline (and never need an EDI dataset).
        from pycsamt.agents.router import (
            IntentRouter,
            QUESTION as _I_QUESTION,
            CODE as _I_CODE,
            META as _I_META,
            CLARIFY as _I_CLARIFY,
        )

        with (
            AGENT_CONFIG.offline()
            if _offline else _nullctx()
        ):
            _router = IntentRouter(
                llm_provider=llm_prov,
                api_key=api_key,
                model=sel_model,
            )
            decision = _router.route(text, history=history)

        _step(f"Intent: {decision.intent}", "done")

        if decision.intent == _I_META:
            _update_job(
                jid,
                status="done",
                result=_capability_text(),
                steps=_JOBS[jid]["steps"],
                kind=KIND_META,
            )
            return
        if decision.intent == _I_CLARIFY:
            _update_job(
                jid,
                status="done",
                result=(
                    decision.clarification
                    or "Could you clarify what you'd"
                    " like me to do — answer a question,"
                    " generate code, or run a workflow?"
                ),
                steps=_JOBS[jid]["steps"],
                kind=KIND_CLARIFY,
            )
            return
        if decision.intent == _I_QUESTION:
            _dispatch_question(
                jid, text,
                llm_prov=llm_prov,
                api_key=api_key,
                sel_model=sel_model,
                offline=_offline,
                history=history,
                step=_step,
            )
            return
        if decision.intent == _I_CODE:
            _dispatch_code(
                jid, text, edi_store, settings,
                workflow=decision.workflow,
                llm_prov=llm_prov,
                api_key=api_key,
                sel_model=sel_model,
                offline=_offline,
                step=_step,
            )
            return

        # ── WORKFLOW / PLOT → run the pipeline ─────
        _step("Classifying workflow...", "done")

        with (
            AGENT_CONFIG.offline()
            if _offline else _nullctx()
        ):
            ctx_agent = ContextInputAgent(
                llm_provider=llm_prov,
                api_key=api_key,
                model=sel_model,
            )
            ctx_result = ctx_agent.execute(
                {"request": text}
            )

        if ctx_result.status == "failed":
            _update_job(
                jid,
                status="done",
                result=(
                    "I could not classify your "
                    "request. Please try rephrasing "
                    "or load an EDI dataset first."
                ),
                steps=_JOBS[jid]["steps"],
                kind=KIND_ERROR,
            )
            return

        cfg = ctx_result.data.get("config", {})

        # inject workflow-specific params
        _ic = inv_config or {}
        # Workflow-type authority (separation of concerns):
        #   1. param-modal selection (explicit user form), else
        #   2. ContextInputAgent classification — the shared, ordered
        #      keyword registry (regex-anchored, LLM-enriched). This is
        #      the authority on *which* workflow to run.
        # The IntentRouter decides the *intent* (question / code /
        # workflow / …) ONLY. It deliberately does not override the
        # workflow type: a terse-prompt LLM router can mislabel e.g.
        # "phase tensor and dimensionality analysis" as an inversion,
        # whereas the regex registry classifies it correctly.
        # Resolve the concrete workflow with an explicit "unknown" path so an
        # unrecognised request is NOT silently run as QC (unprofessional and
        # confusing). Authority order: param-modal form → registry keyword
        # match → LLM router slot → LLM ContextInputAgent (only when it gave a
        # genuine, non-default classification).
        from pycsamt.agents._workflows import classify_workflow as _cwf
        from pycsamt.agents.orchestrator import (
            _WORKFLOW_STEPS as _WF_STEPS,
        )
        _explicit_wf = _ic.get("workflow")
        _kw_wf = _cwf(text, default=None)
        _router_wf = (
            decision.workflow
            if (decision.workflow in _WF_STEPS
                or decision.workflow in _PLOT_WORKFLOWS)
            else None
        )
        _ctx_wf = cfg.get("workflow")
        _resolved_wf = (
            _explicit_wf
            or _kw_wf
            or _router_wf
            or (_ctx_wf if (_ctx_wf and _ctx_wf != "qc") else None)
        )
        # The user genuinely wants QC only when QC keywords actually matched.
        if _resolved_wf is None and _kw_wf == "qc":
            _resolved_wf = "qc"
        _known_wf = (
            _resolved_wf in _WF_STEPS
            or _resolved_wf in _PLOT_WORKFLOWS
            or _resolved_wf in _TOOL_WORKFLOWS
        )
        if not _resolved_wf or not _known_wf:
            _update_job(
                jid,
                status="done",
                result=_unknown_task_text(text),
                steps=_JOBS[jid]["steps"],
                kind=KIND_META,
            )
            return

        cfg["workflow"] = _resolved_wf
        wtype = _resolved_wf
        _update_job(jid, workflow=wtype)
        _step(f"Workflow: {wtype}", "done")
        if wtype in (
            "pinn_inversion", "hybrid_inversion"
        ):
            _pi = {
                "dim": _ic.get("dim", 1),
                "n_layers": _ic.get(
                    "n_layers", 10
                ),
                "depth_max": _ic.get(
                    "depth_max", 2000.0
                ),
                "epochs": _ic.get("epochs", 500),
                "lr": _ic.get("lr", 0.01),
                "smoothness_weight": _ic.get(
                    "smoothness_weight", 0.01
                ),
                "lateral_weight": _ic.get(
                    "lateral_weight", 0.005
                ),
                "graph_weight": _ic.get(
                    "graph_weight", 0.005
                ),
                "radius": _ic.get(
                    "radius", 5000.0
                ),
                "solver": _ic.get(
                    "solver", "mt1d"
                ),
            }
            cfg["pinn_params"] = _pi
            cfg["hybrid_params"] = _pi
            cfg["checkpoint"] = _ic.get(
                "checkpoint", ""
            )
        elif wtype in (
            "ai_inversion", "inv1d",
            "inv2d", "inv3d",
            "ensemble_inversion",
        ):
            cfg["ai_inv_params"] = {
                "n_layers": int(
                    _ic.get("n_layers", 10)
                ),
                "depth_max": float(
                    _ic.get("depth_max", 2000.0)
                ),
                "epochs": int(
                    _ic.get("epochs", 500)
                ),
                "lr": float(
                    _ic.get("lr", 0.01)
                ),
                "lateral_weight": float(
                    _ic.get(
                        "lateral_weight", 0.005
                    )
                ),
                "graph_weight": float(
                    _ic.get(
                        "graph_weight", 0.005
                    )
                ),
                "radius": float(
                    _ic.get("radius", 5000.0)
                ),
            }
            cfg["checkpoint"] = _ic.get(
                "checkpoint", ""
            )

        # Pass pipeline step params into cfg
        _step_p = _ic.get("step_params")
        if _step_p:
            cfg["step_params"] = _step_p

        # build orchestrator input_data
        edi_path = (
            edi_store.get("path", "")
            if edi_store
            else ""
        ) or cfg.get("data_path", "")

        # Filter to selected lines if set
        sel_lines = (edi_store or {}).get(
            "selected_lines", []
        )
        if sel_lines:
            grp = (edi_store or {}).get(
                "groups", {}
            )
            file_list: list[str] = []
            for ln in sel_lines:
                file_list.extend(
                    grp.get(ln, [])
                )
            if file_list:
                edi_path = file_list

        # Fall back to YAML line registry when no
        # EDI is loaded but user names a survey line.
        if not edi_path:
            _reg_yaml = (settings or {}).get(
                "line_registry", ""
            )
            if _reg_yaml:
                try:
                    import yaml as _yaml
                    _reg = (
                        _yaml.safe_load(_reg_yaml)
                        or {}
                    )
                    _tl = text.lower()
                    for _ln, _lp in _reg.items():
                        if (
                            str(_ln).lower() in _tl
                        ):
                            edi_path = str(_lp)
                            break
                except Exception:
                    pass

        # Fall back to the assistant project registry: a named survey
        # line ("run static shift on L22PLT") resolves to its real EDI
        # directory, so workflows run without loading data manually.
        _resolved_line = None
        if not edi_path:
            try:
                from pycsamt.assistant.tools.project_registry import (
                    ProjectRegistry,
                )
                _reg = ProjectRegistry.from_default()
                if _reg is not None:
                    _ln = _reg.find_line_in_text(text)
                    if _ln:
                        _info = _reg.resolve_line(_ln)
                        if _info.get("exists"):
                            edi_path = _info["edi_dir"]
                            _resolved_line = _ln
            except Exception:  # noqa: BLE001
                pass

        # Session fallback: a follow-up that names no line and loads no
        # data inherits the last-used dataset ("now run phase analysis").
        if not edi_path:
            _sess = _session()
            if _sess is not None and _sess.edi_path:
                edi_path = _sess.edi_path
                _resolved_line = _resolved_line or _sess.line

        # Remember the active dataset / line for later turns.
        _sess = _session()
        if _sess is not None and edi_path:
            _sess.set_data(
                edi_path=edi_path,
                line=_resolved_line or _sess.line,
            )
            _sess.record_workflow(wtype)

        output_dir = (
            settings.get("output_dir") or ""
        ).strip() or "pycsamt_workflow_output"

        # Guard: workflows that need EDI data
        # must have a valid path.
        _EDI_REQUIRED = frozenset({
            "qc", "static_shift",
            "phase_analysis", "ai_inversion",
            "inv1d", "inv2d", "inv3d",
            "pinn_inversion", "hybrid_inversion",
            "ensemble_inversion", "pre_inversion",
            "modem", "tipper", "rotation",
            "denoise", "sensitivity",
            "freq_decimation",
            "rhophi", "phase_psection", "pt_psection", "tipper_plot",
            "phase_tensor_map", "station_response", "strike_profile",
            "strike", "dimensionality", "validator",
            "coords", "elevation", "converter", "batch_export",
        })
        if (
            wtype in _EDI_REQUIRED
            and not edi_path
        ):
            _update_job(
                jid,
                status="done",
                result=(
                    "No EDI data loaded. "
                    "Please load an EDI dataset "
                    "first using the Load EDI "
                    "button, then retry."
                ),
                steps=_JOBS[jid]["steps"],
                kind=KIND_ERROR,
            )
            return

        # Friendly label for plot / tool messages (line, selected lines, dir).
        import os as _os
        _task_label = _resolved_line or (
            ", ".join(sel_lines) if sel_lines else ""
        )
        if not _task_label and isinstance(edi_path, str) and edi_path:
            _task_label = _os.path.basename(edi_path.rstrip("/\\"))

        # ── plotting tasks → lightweight PlotAgent (no orchestrator) ──
        if wtype in _PLOT_WORKFLOWS:
            _dispatch_plot(
                jid,
                edi_path,
                kind=_PLOT_KIND[wtype],
                params=(inv_config or {}),
                step=_step,
                label=_task_label,
            )
            return

        # ── analysis tools → lightweight ToolAgent (table + figure) ──
        if wtype in _TOOL_WORKFLOWS:
            _dispatch_tool(
                jid,
                edi_path,
                kind=_TOOL_KIND[wtype],
                params=(inv_config or {}),
                step=_step,
                label=_task_label,
            )
            return

        orch_input = {
            "config":     cfg,
            "request":    text,
            "data_path":  edi_path,
            "output_dir": output_dir,
        }

        with (
            AGENT_CONFIG.offline()
            if _offline else _nullctx()
        ):
            orch = WorkflowOrchestratorAgent(
                llm_provider=llm_prov,
                api_key=api_key,
                model=sel_model,
            )
            _step(
                f"Executing {wtype}...",
                "running",
            )
            result = orch.execute(orch_input)

        _step(
            f"Completed {wtype}",
            "done"
            if result.status != "failed"
            else "error",
        )

        # collect figures and generated code
        # result.data["result"] = coordinator AgentResult
        # coordinator AgentResult.data = {step: AgentResult}
        # each step AgentResult.data may have "figures"
        figs: dict = {}
        generated_code: str = ""
        step_results: dict = {}
        if result.status != "failed":
            exec_res = (
                result.data or {}
            ).get("result")
            step_results = (
                exec_res.data
                if exec_res and exec_res.data
                else {}
            )
            _fig_steps = (
                _WORKFLOW_FIGURE_STEPS.get(wtype)
            )
            # None means "show all steps";
            # a set means only those steps.
            # For unlisted workflows, skip
            # load/denoise by default.
            if (
                _fig_steps is None
                and wtype
                not in _WORKFLOW_FIGURE_STEPS
            ):
                _fig_steps = None  # show all
            for sname, sres in step_results.items():
                if not hasattr(sres, "data"):
                    continue
                # Skip prerequisite step figures
                # unless the workflow explicitly
                # keeps them (None = keep all).
                if _fig_steps is not None:
                    if sname not in _fig_steps:
                        # close any open figs for
                        # this suppressed step
                        _sf = (
                            sres.data or {}
                        ).get("figures", {})
                        for _f in (
                            _sf or {}
                        ).values():
                            if isinstance(
                                _f, plt.Figure
                            ):
                                plt.close(_f)
                        continue
                elif sname in _SKIP_ALWAYS:
                    continue
                step_figs = (
                    sres.data or {}
                ).get("figures", {})
                for fname, fig in (
                    step_figs or {}
                ).items():
                    if isinstance(fig, plt.Figure):
                        b64 = _fig_to_b64(fig)
                        plt.close(fig)
                        figs[str(uuid.uuid4())] = {
                            "title": (
                                f"[{sname}] {fname}"
                            ),
                            "b64": b64,
                        }

            # extract generated code (code_gen step)
            code_res = step_results.get("code_gen")
            if code_res and hasattr(
                code_res, "data"
            ):
                generated_code = (
                    (code_res.data or {}).get(
                        "code", ""
                    ) or ""
                )

        # cache corrected sites for post-proc modal
        if (
            result.status != "failed"
            and wtype in _CORRECTION_WFLOWS
        ):
            for _sn, _sr in (
                step_results.items()
            ):
                _d = (
                    getattr(_sr, "data", None)
                    or {}
                )
                _corr = _d.get(
                    "corrected_sites"
                )
                if _corr is not None:
                    _CORR_CACHE[jid] = _corr
                    _update_job(
                        jid,
                        postproc={
                            "jid": jid,
                            "workflow": wtype,
                            "output_dir": output_dir,
                        },
                    )
                    break

        # AgentResult.summary is the text field
        summary = (
            result.summary
            or result.error
            or "Workflow completed."
        )

        # Trace the run to the workflow history (observability +
        # the sidebar "Recent runs" view). Best-effort.
        _record_run(
            workflow=wtype,
            path=edi_path,
            output_dir=output_dir,
            status=result.status,
            summary=summary,
            n_figures=len(figs),
        )

        _update_job(
            jid,
            status="done",
            result=summary,
            steps=_JOBS[jid]["steps"],
            figs=figs,
            code=generated_code,
            kind=(
                KIND_ERROR
                if result.status == "failed"
                else KIND_WORKFLOW
            ),
        )

    except Exception as exc:  # noqa: BLE001
        _update_job(
            jid,
            status="error",
            error=str(exc),
            result=(
                f"An error occurred: {exc}\n\n"
                "Check that the EDI path is set "
                "and your API key is configured "
                "in Settings if needed."
            ),
            steps=_JOBS[jid]["steps"],
            kind=KIND_ERROR,
        )


# ── pinned-messages helpers ────────────────────────

def _pin_snippet(text: str, limit: int = 60) -> str:
    """One-line snippet for the sidebar pin entry."""
    line = " ".join((text or "").split())
    return line if len(line) <= limit else line[:limit - 1] + "…"


def _apply_pin_toggle(
    pins: list | None,
    mid: str,
    messages: list | None,
) -> list:
    """Toggle *mid* in *pins*: remove if present, else add from *messages*.

    Returns the new pins list. Raises KeyError if the message isn't found
    when adding (caller treats this as "no change").
    """
    pins = list(pins or [])
    if any(p.get("mid") == mid for p in pins):
        return [p for p in pins if p.get("mid") != mid]
    msg = next(
        (m for m in (messages or []) if m.get("mid") == mid),
        None,
    )
    if not msg:
        raise KeyError(mid)
    pins.append({
        "mid": mid,
        "role": msg.get("role", "assistant"),
        "snippet": _pin_snippet(msg.get("content", "")),
        "ts": msg.get("ts", _ts()),
    })
    return pins


def _remove_pin(pins: list | None, mid: str) -> list:
    """Return *pins* without the entry for *mid*."""
    return [
        p for p in (pins or []) if p.get("mid") != mid
    ]


def _pin_item(pin: dict) -> html.Div:
    """Render one pinned-message row in the sidebar."""
    mid = pin.get("mid", "")
    role = pin.get("role", "assistant")
    icon = (
        "bi-person-fill" if role == "user"
        else "bi-robot"
    )
    return html.Div(
        [
            html.Button(
                [
                    html.I(
                        className=f"bi {icon} am-pin-role"
                    ),
                    html.Span(
                        pin.get("snippet", ""),
                        className="am-pin-text",
                    ),
                ],
                id={"type": "am-pin-jump", "mid": mid},
                className="am-pin-jump",
                title="Jump to message",
                n_clicks=0,
            ),
            html.Button(
                html.I(className="bi bi-x"),
                id={"type": "am-unpin", "mid": mid},
                className="am-pin-remove",
                title="Unpin",
                n_clicks=0,
            ),
        ],
        className="am-pin-item",
    )


# ── workflow trace (Recent runs) ───────────────────

def _record_run(
    *,
    workflow: str,
    path: str,
    output_dir: str,
    status: str,
    summary: str,
    n_figures: int,
) -> None:
    """Append a completed workflow to the persistent trace (best-effort)."""
    try:
        from pycsamt.assistant.memory.workflow_history import (
            WorkflowHistory,
            WorkflowRun,
        )
        WorkflowHistory.default().record(
            WorkflowRun(
                workflow=workflow,
                status=status,
                path=path or None,
                output_dir=output_dir,
                summary=summary,
                n_figures=n_figures,
            )
        )
    except Exception:  # noqa: BLE001 — tracing must never break a job
        pass


# Module-level assistant session (the Agent Master is a local,
# single-user app — consistent with _JOBS / _CORR_CACHE singletons).
# Tracks the active data path / line / last workflow so follow-ups like
# "now run phase analysis" inherit context.
_SESSION: Any = None


def _session() -> Any:
    """Lazily create the per-process SessionState (or None if unavailable)."""
    global _SESSION
    if _SESSION is None:
        try:
            from pycsamt.assistant.memory import SessionState
            _SESSION = SessionState()
        except Exception:  # noqa: BLE001
            _SESSION = False  # mark as tried-and-unavailable
    return _SESSION or None


def _session_has_data() -> bool:
    s = _session()
    return bool(s and s.edi_path)


def _reset_session() -> None:
    """Clear the active session (e.g. on New Chat)."""
    global _SESSION
    s = _session()
    if s is not None:
        s.edi_path = None
        s.line = None
        s.last_workflow = None
        s.facts = {}


def _names_registry_line(text: str) -> bool:
    """True when *text* names a survey line that resolves to real data.

    Lets the no-EDI guard pass for "run X on line L22PLT" — _run_agent
    then resolves the line via the project registry.
    """
    try:
        from pycsamt.assistant.tools.project_registry import (
            ProjectRegistry,
        )
        reg = ProjectRegistry.from_default()
        if reg is None:
            return False
        line = reg.find_line_in_text(text)
        if not line:
            return False
        return bool(reg.resolve_line(line).get("exists"))
    except Exception:  # noqa: BLE001
        return False


def _recent_runs(n: int = 8) -> list[dict]:
    """Most-recent workflow runs (newest first), or []."""
    try:
        from pycsamt.assistant.memory.workflow_history import (
            WorkflowHistory,
        )
        return list(reversed(WorkflowHistory.default().recent(n)))
    except Exception:  # noqa: BLE001
        return []


def _run_item(run: dict) -> html.Div:
    """Render one Recent-runs row in the sidebar."""
    wf = run.get("workflow", "?")
    status = run.get("status", "")
    ok = status != "failed"
    ts = (run.get("timestamp", "") or "").replace("T", " ")[5:16]
    icon = (
        "bi-check-circle-fill" if ok
        else "bi-exclamation-triangle-fill"
    )
    color = "var(--green)" if ok else "var(--red)"
    label = _WF_RUNNING_LABEL.get(wf, wf.replace("_", " "))
    return html.Div(
        [
            html.I(
                className=f"bi {icon} am-run-icon",
                style={"color": color},
            ),
            html.Div(
                [
                    html.Div(label, className="am-run-wf"),
                    html.Div(ts, className="am-run-ts"),
                ],
                className="am-run-meta",
            ),
        ],
        className="am-run-item",
        title=run.get("summary", ""),
    )


# ── callbacks ──────────────────────────────────────

def register_chat(app) -> None:

    # 1. Send message → add user bubble,
    #    detect if params needed or start job
    @app.callback(
        Output(IDs.CHAT_WINDOW, "children"),
        Output(IDs.STORE_JOB, "data"),
        Output(
            IDs.INTERVAL_POLL, "disabled"
        ),
        Output(IDs.INPUT, "value"),
        Output(IDs.STORE_MESSAGES, "data"),
        Output(IDs.STORE_PENDING, "data"),
        Input(IDs.BTN_SEND, "n_clicks"),
        Input(
            {"type": "am-chip", "index": ALL},
            "n_clicks",
        ),
        State(IDs.INPUT, "value"),
        State(IDs.CHAT_WINDOW, "children"),
        State(IDs.STORE_EDI, "data"),
        State(IDs.STORE_SETTINGS, "data"),
        State(
            IDs.STORE_INV_CONFIG, "data"
        ),
        State(IDs.STORE_MESSAGES, "data"),
        State(IDs.INTERVAL_POLL, "disabled"),
        State(IDs.STORE_JOB, "data"),
        prevent_initial_call=True,
    )
    def send_message(
        n_send,
        chip_clicks,
        text,
        current_msgs,
        edi_store,
        settings,
        inv_config,
        stored_messages,
        poll_disabled,
        job_store,
    ):
        # ── STOP: button is in "stop" mode while a job runs ──
        # The interval is enabled (disabled is False) only while a
        # background job is active. A click on the send/stop button in
        # that state cancels the job rather than sending a new message.
        # Checked before ctx.triggered_id so the stop path is self
        # contained (and unit-testable without a callback context).
        if poll_disabled is False:
            return _stop_job_response(
                current_msgs, stored_messages, job_store
            )

        triggered = ctx.triggered_id
        if (
            isinstance(triggered, dict)
            and triggered.get("type")
            == "am-chip"
        ):
            idx = triggered["index"]
            from ..layout import _PROMPT_CHIPS
            text = _PROMPT_CHIPS[idx]

        if not text or not text.strip():
            raise PreventUpdate

        text = text.strip()

        _user_mid = _mid()
        new_stored = list(
            stored_messages or []
        )
        new_stored.append({
            "role": "user",
            "content": text,
            "ts": _ts(),
            "mid": _user_mid,
        })

        msgs = [
            c for c in (current_msgs or [])
            if not (
                isinstance(c, dict)
                and c.get("props", {}).get("id")
                == IDs.WELCOME
            )
        ]
        msgs.append(_user_bubble(text, mid=_user_mid))

        # ── web app redirect ──────────────────
        # Explicit launch request, or a task
        # too complex for the chat interface
        # (interactive viz, 3-D maps, etc.).
        # Both bypass the EDI guard because
        # the web app loads its own data.
        _is_web = _is_web_app_request(text)
        _is_viz = (
            not _is_web
            and _is_complex_viz(text)
        )
        if _is_web or _is_viz:
            reason = (
                ""
                if _is_web
                else (
                    "This task requires"
                    " interactive visualization"
                    " tools that go beyond the"
                    " agent chat. Launching the"
                    " full pyCSAMT web app"
                    " instead."
                )
            )
            wa_url = _ensure_web_app()
            wb = _web_app_bubble(wa_url, reason)
            msgs.append(wb)
            _wa_msg = (
                f"Redirecting to web app:"
                f" {wa_url}"
            )
            new_stored.append({
                "role": "assistant",
                "content": _wa_msg,
                "ts": _ts(),
            })
            return (
                msgs, no_update, True,
                "", new_stored, {},
            )

        # Strip any persisted "workflow" from the inversion-config
        # store: the workflow type must be classified per-request from
        # the text, never carried over from a previous param-modal run
        # (STORE_INV_CONFIG is localStorage). Only the param-modal
        # submit (params.py) sets the workflow for its own run.
        _inv_clean = _drop_workflow(inv_config)

        # ── quick intent gate ─────────────────
        # Questions, code, and capability requests
        # don't need data — skip the EDI guard,
        # line picker, and param modal and let the
        # router (inside _run_agent) dispatch them.
        from pycsamt.agents.router import (
            classify_intent_offline,
            NO_DATA_INTENTS,
        )
        _qi, _ = classify_intent_offline(text)
        if _qi in NO_DATA_INTENTS:
            msgs.append(
                _thinking_bubble([{
                    "label": "Parsing request...",
                    "status": "running",
                }])
            )
            jid = _new_job()
            t = threading.Thread(
                target=_run_agent,
                args=(
                    jid, text,
                    dict(edi_store or {}),
                    settings or {},
                    _inv_clean,
                    new_stored,
                ),
                daemon=True,
            )
            t.start()
            return (
                msgs, {"jid": jid}, False,
                "", new_stored, {},
            )

        # ── guard: no EDI loaded ──────────────
        # Skip the guard when the request names a known survey line —
        # _run_agent resolves it via the project registry, so the user
        # need not load EDIs manually for "run X on line L22PLT".
        edi_path = (edi_store or {}).get(
            "path", ""
        )
        if (
            not edi_path
            and not _names_registry_line(text)
            and not _session_has_data()
        ):
            _no_edi = (
                "No EDI dataset is loaded.\n"
                "Please click Load EDI "
                "(top-left) to select your "
                "files, then confirm."
            )
            msgs.append(_agent_bubble(_no_edi))
            new_stored.append({
                "role": "assistant",
                "content": _no_edi,
                "ts": _ts(),
            })
            return (
                msgs, no_update, True,
                "", new_stored, {},
            )

        # ── line disambiguation ───────────────
        # Detect if the user named a specific line
        # that can't be resolved to a group key.
        # If ambiguous: show the line picker modal.
        # If exact match: pre-filter edi_store.
        _groups = (edi_store or {}).get(
            "groups", {}
        )
        _sel_lines: list[str] = []
        if _groups and len(_groups) > 1:
            _ref = _extract_line_ref(
                text, _groups
            )
            if _ref is not None:
                _exact = _match_group(
                    _ref, _groups
                )
                if _exact is None:
                    # Ambiguous → show picker
                    msgs.append(
                        _line_waiting_bubble()
                    )
                    pending = {
                        "disambiguation": "lines",
                        "text": text,
                        "groups": {
                            k: list(v)
                            for k, v in
                            _groups.items()
                        },
                    }
                    return (
                        msgs, {}, True,
                        "", new_stored, pending,
                    )
                else:
                    _sel_lines = [_exact]

        # ── param detection ───────────────────
        wf = _quick_workflow(text)
        if wf:
            msgs.append(_waiting_bubble(wf))
            pending = {
                "workflow": wf,
                "text": text,
                # snapshot edi_path so _submit_params
                # can fall back if STORE_EDI is None
                "edi_path": (
                    (edi_store or {}).get(
                        "path", ""
                    ) or ""
                ),
                "edi_groups": (
                    (edi_store or {}).get(
                        "groups", {}
                    )
                ),
            }
            if _sel_lines:
                pending["selected_lines"] = (
                    _sel_lines
                )
            return (
                msgs, {}, True,
                "", new_stored, pending,
            )

        # ── normal flow: start job ────────────
        msgs.append(
            _thinking_bubble([{
                "label": "Parsing request...",
                "status": "running",
            }])
        )
        _edi_use = dict(edi_store or {})
        if _sel_lines:
            _edi_use["selected_lines"] = (
                _sel_lines
            )
        jid = _new_job()
        t = threading.Thread(
            target=_run_agent,
            args=(
                jid, text,
                _edi_use,
                settings or {},
                _inv_clean,
                new_stored,
            ),
            daemon=True,
        )
        t.start()

        return (
            msgs, {"jid": jid}, False,
            "", new_stored, {},
        )

    # 2. Poll → update thinking / finish
    @app.callback(
        Output(
            IDs.CHAT_WINDOW,
            "children",
            allow_duplicate=True,
        ),
        Output(
            IDs.INTERVAL_POLL,
            "disabled",
            allow_duplicate=True,
        ),
        Output(IDs.STORE_FIGS, "data"),
        Output(
            IDs.STORE_MESSAGES,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_POSTPROC, "data"
        ),
        Input(IDs.INTERVAL_POLL, "n_intervals"),
        State(IDs.STORE_JOB, "data"),
        State(IDs.CHAT_WINDOW, "children"),
        State(IDs.STORE_FIGS, "data"),
        State(IDs.STORE_MESSAGES, "data"),
        prevent_initial_call=True,
    )
    def poll_job(
        _n,
        job_store,
        current_msgs,
        fig_store,
        stored_messages,
    ):
        if not job_store:
            raise PreventUpdate
        jid = job_store.get("jid")
        if not jid:
            raise PreventUpdate

        job = _get_job(jid)
        if not job:
            raise PreventUpdate

        steps = job.get("steps", [])
        status = job.get("status", "running")

        # Job was stopped by the user — halt polling, render nothing.
        # send_message already replaced the thinking bubble.
        if status == "cancelled":
            return (
                no_update, True, no_update,
                no_update, no_update,
            )

        # replace last thinking bubble with
        # updated thinking bubble
        msgs = list(current_msgs or [])
        thinking_idx = None
        for i, child in enumerate(msgs):
            if isinstance(child, dict):
                cid = child.get(
                    "props", {}
                ).get("id", "")
                if cid == "am-thinking-bubble":
                    thinking_idx = i
                    break

        if status == "running":
            _started = job.get("started")
            _elapsed = (
                time.time() - _started
                if _started else None
            )
            new_thinking = _thinking_bubble(
                steps,
                workflow=job.get("workflow"),
                elapsed=_elapsed,
            )
            if thinking_idx is not None:
                msgs[thinking_idx] = new_thinking
            return (
                msgs, False, no_update,
                no_update, no_update,
            )

        # job done / error
        result_text = (
            job.get("result")
            or job.get("error")
            or "Done."
        )
        figs = job.get("figs", {})
        code = job.get("code", "")
        kind = job.get("kind")

        # merge figs into store
        new_fig_store = dict(fig_store or {})
        new_fig_store.update(figs)

        _agent_mid = _mid()
        agent_bub = _agent_bubble(
            result_text, steps, figs,
            code=code, kind=kind, mid=_agent_mid,
        )

        # replace thinking with agent bubble
        if thinking_idx is not None:
            msgs[thinking_idx] = agent_bub
        else:
            msgs.append(agent_bub)

        postproc = job.get("postproc")

        # clean up registry
        with _JOBS_LOCK:
            _JOBS.pop(jid, None)

        new_stored = list(
            stored_messages or []
        )
        new_stored.append({
            "role": "assistant",
            "content": result_text,
            "ts": _ts(),
            "mid": _agent_mid,
        })
        return (
            msgs, True, new_fig_store,
            new_stored,
            postproc if postproc else no_update,
        )

    # 3. Auto-scroll chat to bottom on update
    app.clientside_callback(
        """
        function(children) {
            setTimeout(function() {
                var el = document.getElementById(
                    'am-chat-window'
                );
                if (el) {
                    el.scrollTop = el.scrollHeight;
                }
            }, 60);
            return window.dash_clientside
                .no_update;
        }
        """,
        Output("am-scroll-dummy", "data"),
        Input(IDs.CHAT_WINDOW, "children"),
        prevent_initial_call=True,
    )

    # 3b. Toggle Send ⇄ Stop while a task runs.
    # The polling interval is enabled (disabled is False) only while a
    # background job is active, so it is the single source of truth for
    # "is the assistant busy".
    @app.callback(
        Output(IDs.BTN_SEND, "children"),
        Output(IDs.BTN_SEND, "className"),
        Output(IDs.BTN_SEND, "title"),
        Input(IDs.INTERVAL_POLL, "disabled"),
    )
    def toggle_send_stop(poll_disabled):
        if poll_disabled is False:
            # task running → Stop button
            return (
                html.I(className="bi bi-stop-fill"),
                "am-send-stop",
                "Stop task",
            )
        # idle → Send button
        return (
            html.I(className="bi bi-arrow-up"),
            "",
            "Send",
        )

    # 4. Figure open → modal
    @app.callback(
        Output(IDs.MODAL_FIG, "is_open"),
        Output(IDs.MODAL_FIG_IMG, "src"),
        Output(IDs.MODAL_FIG_TITLE, "children"),
        Output(IDs.MODAL_FIG_KEY, "data"),
        Input(
            {"type": "am-fig-open", "key": ALL},
            "n_clicks",
        ),
        State(IDs.STORE_FIGS, "data"),
        prevent_initial_call=True,
    )
    def open_fig_modal(n_clicks, fig_store):
        # Require the triggered button itself to
        # have a positive click count — prevents
        # auto-open when a new figure card is
        # dynamically added while an older button
        # still has n_clicks > 0.
        if not ctx.triggered:
            raise PreventUpdate
        if not ctx.triggered[0].get("value"):
            raise PreventUpdate
        triggered = ctx.triggered_id
        if not triggered:
            raise PreventUpdate
        fig_key = triggered.get("key")
        if not fig_store or fig_key not in (
            fig_store
        ):
            raise PreventUpdate
        info = fig_store[fig_key]
        src = (
            f"data:image/png;base64,"
            f"{info['b64']}"
        )
        return (
            True,
            src,
            info.get("title", "Figure"),
            fig_key,
        )

    # 5. Pin / unpin a message → STORE_PINS
    @app.callback(
        Output(IDs.STORE_PINS, "data"),
        Input(
            {"type": "am-pin-btn", "mid": ALL},
            "n_clicks",
        ),
        State(IDs.STORE_PINS, "data"),
        State(IDs.STORE_MESSAGES, "data"),
        prevent_initial_call=True,
    )
    def toggle_pin(_clicks, pins, messages):
        if not ctx.triggered or not ctx.triggered[0].get(
            "value"
        ):
            raise PreventUpdate
        trig = ctx.triggered_id
        if not isinstance(trig, dict):
            raise PreventUpdate
        mid = trig.get("mid")
        try:
            return _apply_pin_toggle(pins, mid, messages)
        except KeyError:
            raise PreventUpdate

    # 6. Unpin from the sidebar → STORE_PINS
    @app.callback(
        Output(
            IDs.STORE_PINS, "data",
            allow_duplicate=True,
        ),
        Input(
            {"type": "am-unpin", "mid": ALL},
            "n_clicks",
        ),
        State(IDs.STORE_PINS, "data"),
        prevent_initial_call=True,
    )
    def unpin(_clicks, pins):
        if not ctx.triggered or not ctx.triggered[0].get(
            "value"
        ):
            raise PreventUpdate
        trig = ctx.triggered_id
        if not isinstance(trig, dict):
            raise PreventUpdate
        return _remove_pin(pins, trig.get("mid"))

    # 7. Render the sidebar Pinned section
    @app.callback(
        Output(IDs.SIDEBAR_PINS, "children"),
        Input(IDs.STORE_PINS, "data"),
    )
    def render_pins(pins):
        if not pins:
            return html.Div(
                "No pinned messages yet.",
                className="am-sidebar-empty",
            )
        return [_pin_item(p) for p in pins]

    # 8. Click a pinned item → scroll to the message
    app.clientside_callback(
        """
        function(n_clicks) {
            const t = window.dash_clientside.callback_context
                .triggered;
            if (!t || !t.length || !t[0].value) {
                return window.dash_clientside.no_update;
            }
            let mid;
            try {
                mid = JSON.parse(t[0].prop_id.split('.')[0])
                    .mid;
            } catch (e) {
                return window.dash_clientside.no_update;
            }
            const el = document.getElementById(mid);
            if (el) {
                el.scrollIntoView({
                    behavior: 'smooth', block: 'center',
                });
                el.classList.add('am-msg-flash');
                setTimeout(function () {
                    el.classList.remove('am-msg-flash');
                }, 1600);
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output(IDs.PIN_SCROLL_DUMMY, "data"),
        Input(
            {"type": "am-pin-jump", "mid": ALL},
            "n_clicks",
        ),
        prevent_initial_call=True,
    )

    # 9. Recent runs → render from the workflow trace, refreshed on
    #    every chat update (a completed job updates the chat window).
    @app.callback(
        Output(IDs.SIDEBAR_RUNS, "children"),
        Input(IDs.CHAT_WINDOW, "children"),
    )
    def render_recent_runs(_children):
        runs = _recent_runs()
        if not runs:
            return html.Div(
                "No workflows run yet.",
                className="am-sidebar-empty",
            )
        return [_run_item(r) for r in runs]
