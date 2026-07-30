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
from dash.exceptions import PreventUpdate

from .._ids import IDs

# ── shared job registry ────────────────────────────
_JOBS: dict[str, dict] = {}
_JOBS_LOCK = threading.Lock()

# Corrected-sites cache keyed by job ID.
# Populated after correction workflows so the
# post-processing modal can export EDI files.
_CORR_CACHE: dict[str, Any] = {}

# Workflows that produce corrected_sites data.
_CORRECTION_WFLOWS = frozenset(
    {
        "static_shift",
        "denoise",
        "qc",
        "pre_inversion",
        "full",
    }
)

# Plotting tasks handled by the lightweight PlotAgent (not the orchestrator):
# rho/phi sounding curves, scalar phase pseudo-section, and phase-tensor (Φ)
# ellipse pseudo-section. Each maps to a PlotAgent "kind".
_PLOT_KIND = {
    "rhophi": "rhophi",
    "phase_psection": "phase_psection",
    "pt_psection": "pt_psection",
    "tipper_plot": "tipper",
    "phase_tensor_map": "pt_map",
    "station_response": "station_response",
    "strike_profile": "strike_profile",
    "pt_strip": "pt_strip",
    "pt_strip_grid": "pt_strip_grid",
}
_PLOT_WORKFLOWS = frozenset(_PLOT_KIND)

# Analysis tools handled by the lightweight ToolAgent (table + figure).
_TOOL_KIND = {
    "strike": "strike",
    "dimensionality": "dimensionality",
    "validator": "validator",
    "coords": "coords",
    "elevation": "elevation",
    "converter": "converter",
    "batch_export": "batch_export",
    "freq_editor": "freq_editor",
    "layered_model": "layered_model",
}
# Correction methods (Static Shift, Noise Removal, Tensor Rotation, …) all run
# through the single parameterised ``correction`` ToolAgent kind; the registry
# is the single source of truth (see pycsamt.agents._corrections).
from pycsamt.agents._corrections import (
    CORRECTION_METHODS as _CORR_METHODS,
)

_TOOL_KIND.update({_wf: "correction" for _wf in _CORR_METHODS})
_TOOL_WORKFLOWS = frozenset(_TOOL_KIND)

# Tool tasks that need no loaded EDI dataset (synthetic).
_NO_DATA_WORKFLOWS = frozenset({"layered_model"})

# ── response kinds ─────────────────────────────────
# Each agent response declares its KIND so the chat
# bubble renders the right shape instead of guessing
# from which fields happen to be populated.
KIND_ANSWER = "answer"  # Q&A about the package
KIND_CODE = "code"  # generated script
KIND_WORKFLOW = "workflow"  # pipeline result (+figs)
KIND_CLARIFY = "clarify"  # needs more info
KIND_META = "meta"  # capabilities / chitchat
KIND_ERROR = "error"  # could not proceed

# Per-kind header chip: (icon class, label, css colour var).
_KIND_HEADER: dict[str, tuple[str, str, str]] = {
    KIND_ANSWER: ("bi-chat-left-text", "Answer", "var(--blue)"),
    KIND_CODE: ("bi-code-slash", "Generated code", "var(--green)"),
    KIND_WORKFLOW: ("bi-diagram-3", "Workflow result", "var(--blue)"),
    KIND_CLARIFY: ("bi-question-circle", "Needs input", "var(--yellow)"),
    KIND_META: ("bi-stars", "pyCSAMT assistant", "var(--blue)"),
    KIND_ERROR: ("bi-exclamation-triangle", "Couldn't proceed", "var(--red)"),
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
    "qc": "quality control",
    "static_shift": "static-shift correction",
    "phase_analysis": "phase tensor analysis",
    "denoise": "denoising",
    "tipper": "tipper analysis",
    "sensitivity": "sensitivity / DOI analysis",
    "rotation": "tensor rotation",
    "freq_decimation": "frequency decimation",
    "ai_inversion": "1-D AI inversion",
    "inv1d": "1-D AI inversion",
    "inv2d": "2-D U-Net inversion",
    "inv3d": "3-D GCN inversion",
    "pinn_inversion": "PINN inversion",
    "hybrid_inversion": "hybrid inversion",
    "ensemble_inversion": "ensemble inversion",
    "joint_inversion": "joint inversion",
    "pre_inversion": "inversion preparation",
    "modem": "ModEM preparation",
    "report": "survey report",
    "interpretation": "geological interpretation",
    "interpret": "geological interpretation",
    "forward": "forward modelling",
    "full": "full pipeline",
    "comparison": "inversion comparison",
    "batch": "batch processing",
    "code_gen": "code generation",
    "data_overview": "survey data reading",
}

# Correction-method running labels, merged from the central registry.
_WF_RUNNING_LABEL.update(
    {
        _wf: _meta.get("running_label", "data correction")
        for _wf, _meta in _CORR_METHODS.items()
    }
)


def _fmt_elapsed(seconds: float) -> str:
    """Format seconds as M:SS for the executing header."""
    s = int(max(0, seconds))
    return f"{s // 60}:{s % 60:02d}"


# Phrases that mark the *action* ("generate code") rather than the
# *subject* ("static shift") of a code request. Stripping them lets the
# workflow classifier see the real target.
_CODE_ACTION_PHRASES = (
    "generate code for",
    "write code for",
    "give me code for",
    "show me code for",
    "code example for",
    "sample code for",
    "produce code for",
    "python script for",
    "create a notebook for",
    "notebook for",
    "write a script for",
    "script for",
    "code for",
    "generate code",
    "write code",
    "give me code",
    "show me code",
    "python script",
    "create notebook",
    "write a script",
    "produce code",
    "a script to",
    "script to",
    "code to",
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
    return {k: v for k, v in (inv_config or {}).items() if k != "workflow"}


def _update_job(jid: str, **kw: Any) -> None:
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
    return base64.b64encode(buf.read()).decode()


def _fig_thumb_item(
    fig_key: str,
    title: str,
    b64: str,
) -> html.Div:
    """Compact thumbnail tile inside the accordion."""
    short = title if len(title) <= 28 else title[:26] + "..."
    return html.Div(
        [
            html.Img(
                src=(f"data:image/png;base64,{b64}"),
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
                    html.I(className=("bi bi-arrows-fullscreen me-1")),
                    "View",
                ],
                className=("am-fig-btn am-fig-thumb-btn"),
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
        _fig_thumb_item(k, v["title"], v["b64"]) for k, v in figs.items()
    ]
    header = html.Span(
        [
            html.I(
                className=("bi bi-bar-chart-fill me-2"),
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
            html.I(className="bi bi-folder2-open"),
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
                            className=("am-bubble user"),
                        ),
                        html.Div(
                            _ts(),
                            className="am-ts",
                        ),
                        html.Div(
                            toolbar_btns,
                            className=("am-msg-toolbar"),
                        ),
                    ]
                ),
                style={"maxWidth": "100%"},
            ),
            html.Div(
                html.I(className="bi bi-person-fill"),
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
            html.I(className="bi bi-arrow-repeat am-tl-spin"),
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
    """Claude-style thinking line.

    Instead of a boxed timeline of every step, a single quiet line shows
    the *current* step with a text shimmer; when the agent advances, the
    old label fades away above ("ghost") and the new one slides in.  A
    2 px hairline tracks n_done/n_total, and hovering the line reveals
    the full step timeline as a floating panel (CSS :hover survives the
    600 ms poll re-renders, unlike a <details> open state).
    """
    steps = steps or []
    n_total = len(steps)
    n_done = sum(1 for s in steps if s.get("status") == "done")
    name = _WF_RUNNING_LABEL.get(workflow, workflow) if workflow else None

    # Current step = last non-done step, else last step, else a
    # workflow-level fallback while the first step is being created.
    cur_idx = None
    for i in range(n_total - 1, -1, -1):
        if steps[i].get("status") != "done":
            cur_idx = i
            break
    if cur_idx is None and n_total:
        cur_idx = n_total - 1
    if cur_idx is not None:
        cur_label = steps[cur_idx].get("label", "...")
        cur_status = steps[cur_idx].get("status", "running")
    else:
        cur_label = f"Running {name}..." if name else "Thinking..."
        cur_status = "running"

    # Ghost: the most recent completed step other than the current one.
    ghost = None
    for i in range(n_total - 1, -1, -1):
        if i != cur_idx and steps[i].get("status") == "done":
            ghost = (i, steps[i].get("label", ""))
            break

    # The ids embed the step index: across 600 ms polls React keeps the
    # same DOM node (animations keep running); when the step changes the
    # id changes, the node is re-mounted, and the enter/fade animations
    # replay exactly once per transition.
    meta_bits = []
    if n_total:
        meta_bits.append(f"{n_done}/{n_total}")
    if elapsed is not None:
        meta_bits.append(_fmt_elapsed(elapsed))
    line = html.Div(
        [
            html.I(
                className=(
                    "bi bi-stars am-think-glyph"
                    + (" error" if cur_status == "error" else "")
                )
            ),
            html.Span(
                cur_label,
                className="am-think-lbl",
                id=(f"am-think-lbl-{cur_idx if cur_idx is not None else 'x'}"),
            ),
            html.Span(
                " · ".join(meta_bits),
                className="am-think-meta",
            ),
        ],
        className="am-think-line",
    )

    ghost_el = (
        html.Div(
            ghost[1],
            className="am-think-ghost",
            id=f"am-think-ghost-{ghost[0]}",
        )
        if ghost
        else None
    )

    # Hairline progress: determinate once steps exist, sweep before.
    if n_total:
        pct = int(round(100.0 * n_done / n_total))
        track = html.Div(
            html.I(style={"width": f"{max(pct, 4)}%"}),
            className="am-think-track",
        )
    else:
        track = html.Div(html.I(), className="am-think-track indet")

    # Full timeline, revealed on hover (reuses the rail/dot rows).
    panel_children = []
    if name:
        panel_children.append(
            html.Div(
                f"Running {name}",
                className="am-think-panel-title",
            )
        )
    panel_children.append(
        html.Div(
            [
                _exec_step_row(
                    s.get("label", ""),
                    s.get("status", ""),
                )
                for s in steps
            ],
            className="am-tl",
        )
        if steps
        else html.Div(
            "Warming up...",
            className="am-think-panel-empty",
        )
    )
    panel = html.Div(panel_children, className="am-think-panel")

    body = [ghost_el] if ghost_el is not None else []
    body += [line, track, panel]

    return html.Div(
        html.Div(body, className="am-think"),
        className="am-msg-row am-msg-row--think",
        id="am-thinking-bubble",
    )


def _strip_thinking(msgs: list | None) -> list:
    """Return *msgs* without the live thinking bubble.

    Used when a running task is stopped: the in-progress bubble is
    removed so a final notice can take its place.
    """
    out = []
    for child in msgs or []:
        if (
            isinstance(child, dict)
            and child.get("props", {}).get("id") == "am-thinking-bubble"
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
    new_stored.append(
        {
            "role": "assistant",
            "content": "Task stopped by user.",
            "ts": _ts(),
        }
    )
    return msgs, {}, True, no_update, new_stored, {}


# Pure markdown parsing lives in a dash-free module so it can be
# unit-tested without importing Dash / the GUI package.
from .._markdown import (
    split_inline_bold as _split_inline_bold,
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
                                    "bi bi-chevron-right am-code-chevron"
                                )
                            ),
                            html.I(
                                className=("bi bi-code-slash me-1"),
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
                    html.I(className="bi bi-clipboard me-1"),
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
    """Render the reply as full GitHub-flavoured markdown.

    Uses :class:`dash.dcc.Markdown` so bold, italic, headings, nested lists,
    inline & fenced code, blockquotes, tables and links all format themselves
    — exposed like ChatGPT/Claude rather than a custom token subset. Colours
    come from the themed ``.am-md`` CSS, so light and dark both work."""
    if not (text or "").strip():
        return []
    return [
        dcc.Markdown(
            text,
            className="am-md",
            link_target="_blank",
        )
    ]


# Kind label chips ("Needs input", "Couldn't proceed") are retired:
# every reply — errors and clarifications included — renders as plain
# flat text in the main chat, exactly like the Claude / ChatGPT reply
# stream. Add a kind back here to restore its chip.
_HEADER_KINDS: frozenset = frozenset()


def _kind_header(kind: str | None) -> html.Div | None:
    """Build the small per-kind label chip, or None for plain replies."""
    if kind not in _HEADER_KINDS:
        return None
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
    card: dict | None = None,
) -> html.Div:
    children: list = []

    # per-kind label chip (answer / code / clarify / …)
    header = _kind_header(kind)
    if header is not None:
        children.append(header)

    # body text — full lightweight-markdown rendering
    children.extend(_render_markdown(text))

    # structured payloads (e.g. the survey data overview)
    if card:
        children.append(_data_overview_card(card))

    # step summary — collapsed by default once the request is done. The live
    # "thinking" bubble shows the full timeline; here we mask it behind a
    # compact, clickable <details> ("✓ N steps · trace") so the chat stays
    # clean, with a side button to surface the run in the sidebar.
    if steps:
        n_steps = len(steps)
        step_divs = [
            html.Div(
                [
                    html.I(className=("bi bi-check-circle-fill am-step-icon")),
                    html.Span(s["label"]),
                ],
                className="am-step done",
            )
            for s in steps
        ]
        _trace_idx = mid or uuid.uuid4().hex
        children.append(
            html.Div(
                [
                    html.Details(
                        [
                            html.Summary(
                                [
                                    html.I(
                                        className=(
                                            "bi bi-check2-circle am-trace-tick"
                                        )
                                    ),
                                    html.Span(
                                        f"{n_steps} step"
                                        f"{'s' if n_steps != 1 else ''}",
                                        className="am-trace-count",
                                    ),
                                    html.Span(
                                        "trace",
                                        className="am-trace-label",
                                    ),
                                    html.I(
                                        className=(
                                            "bi bi-chevron-right am-trace-chev"
                                        )
                                    ),
                                ],
                                className="am-trace-summary",
                            ),
                            html.Div(step_divs, className="am-steps"),
                        ],
                        className="am-trace",
                    ),
                    html.Button(
                        html.I(className="bi bi-box-arrow-up-right"),
                        id={
                            "type": "am-trace-open",
                            "index": _trace_idx,
                        },
                        className="am-trace-open-btn",
                        title="Open trace in the sidebar",
                        n_clicks=0,
                    ),
                ],
                className="am-trace-row",
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
                        className=("bi bi-bar-chart-fill me-2"),
                        style={"color": "var(--blue)"},
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
                    html.Div(_ts(), className="am-ts"),
                    html.Div(
                        [
                            html.Button(
                                html.I(className="bi bi-clipboard"),
                                className=("am-msg-action am-copy-btn"),
                                title="Copy answer",
                                n_clicks=0,
                                **{"data-copy": text or ""},
                            ),
                            _pin_button(mid),
                        ],
                        className="am-msg-toolbar",
                    )
                    if mid
                    else html.Div(),
                ],
                className=(
                    "am-bubble agent am-bubble-flat"
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
_WORKFLOW_FIGURE_STEPS: dict[str, set[str] | None] = {
    "qc": {"qc", "static_shift", "report"},
    "static_shift": {"static_shift"},
    "phase_analysis": {"phase_analysis", "report"},
    "ai_inversion": {"ai_inv", "interpret", "report"},
    "inv1d": {"ai_inv", "interpret", "report"},
    "inv2d": {"inv2d", "interpret", "report"},
    "inv3d": {"inv3d", "interpret", "report"},
    "ensemble_inversion": {
        "ensemble",
        "interpret",
        "report",
    },
    "pinn_inversion": {
        "pinn_inv",
        "interpret",
        "report",
    },
    "hybrid_inversion": {
        "hybrid_inv",
        "interpret",
        "report",
    },
    "joint_inversion": {
        "joint",
        "interpret",
        "report",
    },
    "tipper": {"tipper", "report"},
    "modem": {"modem", "report"},
    "mare2dem": {"mare2dem", "report"},
    "occam2d": {"occam2d", "report"},
    "pre_inversion": {
        "phase_analysis",
        "occam2d",
        "report",
    },
    "sensitivity": {"sensitivity", "report"},
    "rotation": {"rotate"},
    "freq_decimation": {"decimate", "report"},
    "batch": {"batch", "report"},
    "comparison": {"compare", "report"},
    "full": None,  # show everything
}
# Default set for unlisted workflows: skip
# only load and denoise steps.
_SKIP_ALWAYS: frozenset[str] = frozenset({"load", "denoise"})


# ── Smart param detection ─────────────────────────

# Workflows that require user parameters
# before the job can start.
# Workflows that REQUIRE user parameters before
# starting.  Simple data-processing workflows
# (qc, static_shift, phase_analysis, tipper,
# rotation) run with sensible defaults and do
# NOT need a param modal.  Only inversion and
# a few analytical workflows need user input.
_NEEDS_PARAMS: frozenset[str] = frozenset(
    {
        # Inversion workflows (model params needed)
        "ai_inversion",
        "inv1d",
        "inv2d",
        "inv3d",
        "pinn_inversion",
        "hybrid_inversion",
        "ensemble_inversion",
        "pre_inversion",
        "modem",
        "mare2dem",
        # Analysis workflows (optional params)
        "denoise",
        "sensitivity",
        "interpret",
        "interpretation",
        "report",
        "code_gen",
        # Plotting tasks (station / component / period / publication)
        "rhophi",
        "phase_psection",
        "pt_psection",
        "tipper_plot",
        "phase_tensor_map",
        "station_response",
        "strike_profile",
        "pt_strip",
        "pt_strip_grid",
        "strike",
        "dimensionality",
        "validator",
        # Data / IO tools
        "coords",
        "elevation",
        "converter",
        "batch_export",
        # Stateful tools
        "freq_editor",
        "layered_model",
    }
)

# Correction methods always open the parameter modal (the parameter set IS the
# "different correction with different control").
_NEEDS_PARAMS = _NEEDS_PARAMS | frozenset(_CORR_METHODS)

# Inversion-preparation workflows: their deliverable is a per-line input
# file set on disk (Occam2D / ModEM / MARE2DEM).  With several survey
# lines loaded they run once per selected line, and with no line named
# the line picker asks the user first ("Run all" builds every line).
_PREP_WORKFLOWS: frozenset[str] = frozenset(
    {
        "pre_inversion",
        "modem",
        "mare2dem",
    }
)

_WF_LABELS: dict[str, str] = {
    "ai_inversion": "1-D AI inversion",
    "inv1d": "1-D AI inversion",
    "inv2d": "2-D U-Net inversion",
    "inv3d": "3-D GCN inversion",
    "pinn_inversion": "PINN inversion",
    "hybrid_inversion": "hybrid inversion",
    "ensemble_inversion": "ensemble inversion",
    "pre_inversion": "pre-inversion setup",
    "modem": "ModEM preparation",
    "mare2dem": "MARE2DEM preparation",
    "qc": "QC pipeline",
    "phase_analysis": "phase tensor analysis",
    "static_shift": "static shift correction",
    "tipper": "tipper analysis",
    "rotation": "data rotation",
    "interpret": "geological interpretation",
    "interpretation": "geological interpretation",
    "report": "survey report",
    "code_gen": "code generation",
    "denoise": "data denoising",
    "sensitivity": "sensitivity / DOI analysis",
    "rhophi": "rho/phi sounding curves",
    "phase_psection": "phase pseudo-section",
    "pt_psection": "phase-tensor pseudo-section",
    "pt_strip": "phase-tensor ellipse strip",
    "pt_strip_grid": "phase-tensor strip grid (by line)",
    "tipper_plot": "tipper plot",
    "phase_tensor_map": "phase-tensor map",
    "station_response": "station response",
    "strike_profile": "strike profile",
    "strike": "strike analyzer",
    "dimensionality": "dimensionality classifier",
    "validator": "EDI validator",
    "coords": "coordinate transformer",
    "elevation": "elevation enrichment",
    "converter": "format converter",
    "batch_export": "batch plot export",
    "freq_editor": "frequency editor",
    "layered_model": "layered model builder",
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
            _normalise_config,
            _regex_extract,
        )

        cfg = _regex_extract(text)
        cfg = _normalise_config(cfg, text)
        wf = cfg.get("workflow", "")
        return wf if wf in _NEEDS_PARAMS else None
    except Exception:
        return None


def _waiting_bubble(wf: str) -> html.Div:
    label = _WF_LABELS.get(wf, wf.replace("_", " "))
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
                                className=("bi bi-sliders2 me-2"),
                                style={"color": ("var(--blue)")},
                            ),
                            html.Span(
                                f"To run {label}, I need a few parameters."
                            ),
                            html.Span(
                                " Please fill in"
                                " the form above"
                                " and click"
                                " Run Workflow.",
                                style={"opacity": ".7"},
                            ),
                        ]
                    ),
                    html.Div(_ts(), className="am-ts"),
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
                                className=("bi bi-layers me-2"),
                                style={"color": ("var(--blue)")},
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
                    html.Div(_ts(), className="am-ts"),
                ],
                className="am-bubble agent",
            ),
        ],
        className="am-msg-row",
        id="am-line-waiting-bubble",
    )


def _extract_line_ref(text: str, groups: dict) -> str | None:
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
    m = re.search(r"\b(?:line|profile)\s+(\S+)", t)
    if m:
        tok = m.group(1).strip(".,;:!?\"'”’")
        if tok:
            return tok
    _ordinals = {
        "first": "1",
        "second": "2",
        "third": "3",
        "fourth": "4",
        "fifth": "5",
        "sixth": "6",
        "seventh": "7",
        "eighth": "8",
        "ninth": "9",
        "tenth": "10",
    }
    for word, num in _ordinals.items():
        if f"{word} line" in t or f"{word} profile" in t:
            return num
    return None


def _match_group(ref: str, groups: dict) -> str | None:
    """Return the group key matching ref.

    Exact first, then case-insensitive,
    then — so "line 22" finds L22PLT — a
    unique embedded-number match for
    numeric refs, or a unique substring
    match for named refs ("l22").
    Ambiguous or unmatched refs return
    None and go to the line picker; an
    ordinal like "line 2" stays
    picker-bound unless a group actually
    carries that number.
    """
    if ref in groups:
        return ref
    ref_l = ref.lower().strip()
    for key in groups:
        if key.lower() == ref_l:
            return key
    if ref_l.isdigit():
        n = int(ref_l)
        hits = [
            key
            for key in groups
            if n in [int(run) for run in re.findall(r"\d+", key)]
        ]
    else:
        hits = [key for key in groups if ref_l in key.lower()]
    if len(hits) == 1:
        return hits[0]
    return None


# ── PINN / Hybrid keyword detection ───────────────

_PINN_KWS = frozenset(
    {
        "pinn",
        "physics-informed",
        "physics informed",
        "no training data",
        "gradient descent",
    }
)
_HYBRID_KWS = frozenset(
    {
        "hybrid",
        "two-stage",
        "two stage",
        "ai + physics",
        "warm start",
        "warmstart",
    }
)


def _is_pinn_or_hybrid(text: str) -> bool:
    t = text.lower()
    return any(kw in t for kw in _PINN_KWS) or any(
        kw in t for kw in _HYBRID_KWS
    )


# ── Application launchers (web / mapview / desktop) ──────────────

_WEB_APP_KWS: frozenset[str] = frozenset(
    {
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
    }
)

# Tasks too complex for chat → redirect to the web app.
# (Map-flavoured phrases route to MapView below; "phase tensor map"
# is deliberately absent — it is a chat plot workflow.)
_COMPLEX_VIZ_KWS: frozenset[str] = frozenset(
    {
        "explore results",
        "browse data",
        "browse edis",
        "full visualization",
        "interactive plot",
        "full pipeline editor",
    }
)

# Explicit MapView requests.
_MAPVIEW_KWS: frozenset[str] = frozenset(
    {
        "mapview",
        "map view",
        "map workbench",
        "open the map",
        "launch the map",
        "open map",
        "launch map",
        "start mapview",
        "map platform",
        "station map viewer",
    }
)

# Map-flavoured visualisation → MapView is its home turf.
_MAPVIEW_VIZ_KWS: frozenset[str] = frozenset(
    {
        "3d map",
        "station map",
        "interactive map",
        "3d visualization",
        "pseudosection viewer",
        "interactive pseudosection",
    }
)

# Native desktop application.
_DESKTOP_KWS: frozenset[str] = frozenset(
    {
        "desktop app",
        "desktop application",
        "desktop gui",
        "open the gui",
        "launch the gui",
        "start the gui",
        "gui app",
        "launch desktop",
        "open desktop",
        "start the desktop",
        "native app",
        "pycsamt gui",
        "pycsamt desktop",
    }
)

_webapp_state: dict = {
    "url": None,
    "running": False,
}
_webapp_lock = threading.Lock()

# Subprocess-backed apps (MapView Dash server, Qt desktop app).
_EXT_APPS: dict[str, dict] = {
    "mapview": {"proc": None, "url": None},
    "desktop": {"proc": None},
}
_ext_lock = threading.Lock()

_MAPVIEW_PORT = 8770


def _is_web_app_request(text: str) -> bool:
    t = text.lower()
    return any(kw in t for kw in _WEB_APP_KWS)


def _is_complex_viz(text: str) -> bool:
    t = text.lower()
    return any(kw in t for kw in _COMPLEX_VIZ_KWS)


_VIZ_REDIRECT_REASON = (
    "This task needs interactive visualization tools that go "
    "beyond the agent chat — launching the right app for it."
)


def _detect_app_request(text: str) -> tuple[str, str] | None:
    """Return ``(app_kind, reason)`` when *text* asks to open an
    application (or needs one): ``desktop`` / ``mapview`` / ``web``."""
    t = (text or "").lower()
    if any(kw in t for kw in _DESKTOP_KWS):
        return "desktop", ""
    if any(kw in t for kw in _MAPVIEW_KWS):
        return "mapview", ""
    if any(kw in t for kw in _MAPVIEW_VIZ_KWS):
        return "mapview", _VIZ_REDIRECT_REASON
    if _is_web_app_request(text):
        return "web", ""
    if _is_complex_viz(text):
        return "web", _VIZ_REDIRECT_REASON
    return None


def _port_in_use(port: int) -> bool:
    import socket as _socket

    with _socket.socket(_socket.AF_INET, _socket.SOCK_STREAM) as s:
        s.settimeout(0.4)
        return s.connect_ex(("127.0.0.1", port)) == 0


def _ensure_mapview() -> str:
    """Start the MapView workbench as a subprocess (reusing a live
    one) and return its URL."""
    import subprocess
    import sys as _sys

    with _ext_lock:
        st = _EXT_APPS["mapview"]
        if st["proc"] is not None and st["proc"].poll() is None and st["url"]:
            return st["url"]
        if _port_in_use(_MAPVIEW_PORT):
            # A MapView (or a previous launch) already owns the
            # default port — just point at it.
            st["url"] = f"http://127.0.0.1:{_MAPVIEW_PORT}"
            return st["url"]
        port = _free_port(_MAPVIEW_PORT)
        st["proc"] = subprocess.Popen(
            [
                _sys.executable,
                "-m",
                "pycsamt.app.mapview",
                "--port",
                str(port),
                "--no-browser",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        st["url"] = f"http://127.0.0.1:{port}"
        return st["url"]


def _ensure_desktop() -> tuple[bool, str]:
    """Launch the native desktop app as a subprocess.

    Returns ``(ok, note)`` — the note explains what to expect or
    what went wrong."""
    import subprocess
    import sys as _sys

    with _ext_lock:
        st = _EXT_APPS["desktop"]
        if st["proc"] is not None and st["proc"].poll() is None:
            return True, ("It's already running — check your open windows.")
        try:
            st["proc"] = subprocess.Popen(
                [_sys.executable, "-m", "pycsamt.app.desktop"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except Exception as exc:  # noqa: BLE001
            return False, f"It could not be started ({exc})."
        proc = st["proc"]
    # Catch immediate exits (typically PySide6 not installed).
    time.sleep(1.2)
    if proc.poll() is not None:
        return False, (
            "It exited immediately — PySide6 is probably not "
            "installed. Install the desktop extra with "
            '`pip install "pycsamt[desktop]"` and try again.'
        )
    return True, "A native pyCSAMT window should appear shortly."


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


# Per-app launch card copy: title, icon, description.
_APP_CARD: dict[str, tuple[str, str, str]] = {
    "web": (
        "Launching pyCSAMT Web App",
        "bi bi-window-fullscreen",
        "The full pyCSAMT web application provides interactive"
        " station maps, pseudosection viewers, phase-tensor tools,"
        " pipeline editor, inversion pages, and more — beyond what"
        " the agent chat can display.",
    ),
    "mapview": (
        "Launching MapView",
        "bi bi-map",
        "The MapView workbench renders interactive station maps,"
        " pseudosections, overlays, and 3-D quick-look views of"
        " your survey lines.",
    ),
    "desktop": (
        "Launching pyCSAMT Desktop",
        "bi bi-display",
        "The native desktop application gives point-and-click"
        " access to loading, processing, corrections, maps, and"
        " inversion workflows.",
    ),
}


def _launch_bubble(
    app_kind: str,
    url: str | None = None,
    reason: str = "",
    note: str = "",
    ok: bool = True,
) -> html.Div:
    """Chat card announcing an application launch (link when the
    app is a local web server, plain note for the desktop app)."""
    title, icon, default_desc = _APP_CARD[app_kind]
    desc = reason or default_desc

    card: list = [
        html.Div(
            [
                html.I(
                    className=f"{icon} me-2",
                    style={"color": ("var(--green)" if ok else "var(--red)")},
                ),
                html.Strong(
                    title
                    if ok
                    else title.replace("Launching", "Could not launch")
                ),
            ],
            className="am-webapp-hdr",
        ),
        html.P(desc, className="am-webapp-desc"),
    ]
    if url:
        card.append(
            html.A(
                [
                    html.I(className=("bi bi-box-arrow-up-right me-2")),
                    url,
                ],
                href=url,
                target="_blank",
                rel="noopener",
                className="am-webapp-link",
            )
        )
    default_note = (
        "Server is starting — the link will be ready in a few seconds."
        if url
        else ""
    )
    note = note or default_note
    if note:
        card.append(
            html.Div(
                [
                    html.I(className="bi bi-info-circle me-1"),
                    note,
                ],
                className="am-webapp-note",
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
                    html.Div(card, className="am-webapp-card"),
                    html.Div(_ts(), className="am-ts"),
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
        " **OpenAI**, **Google Gemini**, **DeepSeek**, or **MiniMax**. With a key I switch"
        " to online mode, which most improves:\n"
        "- **Questions** about pyCSAMT — fluent, synthesised answers"
        " (grounded in the package via RAG) instead of the offline summary;\n"
        "- **Code generation** — complete, tailored scripts;\n"
        "- **Request understanding** — better routing of free-form requests.\n"
        "Running workflows works offline too — no key required."
    )


def _correction_capability_block() -> str:
    """Markdown bullet list of correction categories, built from the registry.

    Keeps the capability summary in sync with whatever correction methods are
    registered (one bullet per catalogue category, with a few example methods).
    """
    from collections import OrderedDict

    by_cat: OrderedDict[str, list[str]] = OrderedDict()
    for meta in _CORR_METHODS.values():
        by_cat.setdefault(meta.get("category", "Correction"), []).append(
            meta.get("label", "")
        )
    if not by_cat:
        return ""
    lines = [
        "**Correct & condition** your data with full parameter control"
        " (I'll open a parameter form, then you can apply the result to"
        " the session or export it):"
    ]
    for cat, labels in by_cat.items():
        shown = ", ".join(labels[:4])
        more = f", +{len(labels) - 4} more" if len(labels) > 4 else ""
        lines.append(f"- **{cat}** — {shown}{more}")
    return "\n".join(lines) + "\n\n"


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
        "- `occam2d`, `modem`, `mare2dem` — inversion input preparation"
        " (per survey line; I'll ask which line(s) to build)\n"
        "- `report` — generate a survey report\n\n"
        "**Make plots** from your data (I'll ask for stations, components,"
        " period range and publication styling):\n"
        "- rho/φ sounding curves\n"
        "- phase pseudo-section\n"
        "- phase-tensor (Φ) ellipse pseudo-section\n"
        "- phase-tensor map (ellipses at a chosen period)\n"
        "- phase-tensor ellipse strip (single station vs period)\n"
        "- phase-tensor strip grid (ellipse strips tiled by survey line)\n"
        "- station response inspector (Bode ρa/φ curves)\n"
        "- strike profile (strike vs station position)\n"
        "- tipper components / induction arrows (when tipper data exists)\n\n"
        "**Analyze** your data (results as a table + figure):\n"
        "- strike analyzer (geoelectric strike per station)\n"
        "- dimensionality classifier (1-D / 2-D / 3-D)\n"
        "- EDI validator (per-station quality checklist)\n\n"
        + _correction_capability_block()
        + "**Data & I/O tools** (I'll ask for the options first):\n"
        "- coordinate transformer (station lat/lon → UTM)\n"
        "- elevation enrichment (fetch elevation from an open web API)\n"
        "- format converter (re-export the survey to CSV / JSON / EDI)\n"
        "- batch plot export (render & save a bundle of figures)\n\n"
        "**Edit & model**:\n"
        "- frequency editor (confidence QC: drop / mask / recover periods,"
        " then apply or export the edited survey)\n"
        "- layered model builder (preview a synthetic 1-D resistivity model;"
        " no data needed)\n\n"
        "**Tell you values** about a line — just ask and I'll compute and"
        " answer inline (for one line or all lines):\n"
        "- strike, azimuth/bearing, dimensionality (1-D/2-D/3-D), skew\n"
        "- station count, period & frequency range, coordinates & length\n"
        "- data-quality score, and a one-line summary"
        ' ("tell me about L22PLT")\n\n'
        "**Answer questions** about pyCSAMT — classes, functions, the Sites"
        " data model, and which method to use.\n\n"
        "**Generate Python code** that reproduces a pyCSAMT workflow.\n\n"
        "To run a workflow, load an EDI dataset first with **Load EDI**"
        " (top-left). Questions and code requests need no data.\n\n"
        "**Launch the other pyCSAMT apps** — just say the word and I'll"
        " start them for you:\n"
        "- *“open the map view”* — MapView workbench (station maps,"
        " pseudosections, 3-D quick looks);\n"
        "- *“open the web app”* — the full web application (visual pipeline"
        " editor, inversion pages);\n"
        "- *“open the desktop app”* — the native desktop GUI.\n\n"
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
        " hybrid, ensemble, joint), ModEM/Occam2D/MARE2DEM prep, and"
        " reports.\n"
        "- **Answer questions** about pyCSAMT (classes, functions, methods).\n"
        "- **Generate Python code** for a pyCSAMT workflow.\n"
        '- **List my capabilities** — just ask "what can you do?".\n\n'
        "If you need maps, interactive pseudosection viewers, or a feature"
        " that isn't here, say *“open the map view”*, *“open the web"
        " app”*, or *“open the desktop app”* and I'll launch it.\n\n"
        'Tip: phrase it as an action + target, e.g. "run static shift",'
        ' "denoise the data", or "run AI inversion".\n\n' + _api_key_hint()
    )


_PLOT_VERB_RE = re.compile(
    r"\b(plot|draw|display|render|graph|visuali[sz]e|show|view)\b"
)

# Generic analysis verbs — a request like "analyse only the line 22"
# names no concrete workflow; it resolves to the flagship phase-tensor
# analysis instead of bouncing an "unknown task" reply at the user.
_GENERIC_ANALYSIS_RE = re.compile(
    r"\b(?:re-?)?analy(?:se|ze|sis|ses|sed|zed|sing|zing)\b",
    re.IGNORECASE,
)

# (label, example) — the examples use the workflow-registry trigger
# phrases, so the suggested wording routes straight to the right plot.
_PLOT_MENU: tuple[tuple[str, str], ...] = (
    ("ρa/φ sounding curves", "plot the sounding curves"),
    ("phase pseudo-section", "plot the phase pseudosection"),
    ("phase-tensor ellipse section", "plot the phase tensor section"),
    ("phase-tensor map", "plot the phase tensor map"),
    ("phase-tensor ellipse strip", "plot the phase tensor ellipse strip"),
    ("phase-tensor strip grid (by line)", "plot the phase tensor strip grid"),
    ("station response (Bode)", "plot the station response"),
    ("strike profile", "plot the strike profile"),
    ("tipper / induction arrows", "plot the tipper"),
)


def _lines_bullets(groups: dict) -> str:
    return "\n".join(
        f"- **{ln}** — {len(groups[ln])} EDI file"
        f"{'s' if len(groups[ln]) != 1 else ''}"
        for ln in sorted(groups)
    )


def _smart_unknown_reply(text: str, edi_store: dict) -> str | None:
    """Context-aware reply for an unrecognised request.

    Handles the two frequent cases professionally instead of the
    generic capability dump: a request naming a line that isn't
    loaded (propose the loaded lines) and a plot request without a
    recognisable plot kind (propose the plot menu for the resolved
    line). Returns ``None`` when neither case applies.
    """
    t = (text or "").lower()
    plotish = bool(_PLOT_VERB_RE.search(t))
    groups = (edi_store or {}).get("groups", {}) or {}

    line_lbl: str | None = None
    ordinal_note = ""
    sel = [
        str(s)
        for s in ((edi_store or {}).get("selected_lines") or [])
        if not groups or s in groups
    ]
    if sel:
        # Line(s) already chosen via the line picker — that settles
        # the line question, so a textual ref ("line 22") is never
        # re-litigated into a "couldn't find that line" bounce.
        line_lbl = ", ".join(sel)
    elif groups:
        ref = _extract_line_ref(text, groups)
        if ref is not None:
            hit = _match_group(ref, groups)
            if hit is not None:
                line_lbl = hit
            else:
                keys = sorted(groups)
                if ref.isdigit() and 1 <= int(ref) <= len(keys):
                    # "line 2" → the 2nd loaded line
                    line_lbl = keys[int(ref) - 1]
                    ordinal_note = f" (your “line {ref}”)"
                else:
                    return (
                        f"I couldn't find a line called "
                        f"**“{ref}”**. Here's what's loaded:\n\n"
                        f"{_lines_bullets(groups)}\n\n"
                        "Name one of these — e.g. "
                        f"*“plot the phase pseudosection of "
                        f"{sorted(groups)[0]}”*."
                    )

    if not plotish:
        return None

    where = f" for **{line_lbl}**{ordinal_note}" if line_lbl else ""
    suffix = f" of {line_lbl}" if line_lbl else ""
    menu = "\n".join(
        f"- **{lbl}** — *“{ex}{suffix}”*" for lbl, ex in _PLOT_MENU
    )
    return (
        f"Happy to plot something{where} — which figure would "
        "you like?\n\n"
        f"{menu}\n\n"
        "I can also run analyses that end in a figure: the "
        "*strike analyzer*, the *dimensionality classifier*, or a "
        "full *phase-tensor analysis*."
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

    # Conversational context so a subject-less follow-up ("now do that",
    # "same for line 3") inherits the active workflow / line for retrieval.
    sess = _session()
    session_ctx = {
        "last_workflow": getattr(sess, "last_workflow", None)
        if sess
        else None,
        "last_line": getattr(sess, "line", None) if sess else None,
        "recent_turns": history[-6:] if history else None,
    }

    with AGENT_CONFIG.offline() if offline else _nullctx():
        qa = PackageQAAgent(
            llm_provider=llm_prov,
            api_key=api_key,
            model=sel_model,
        )
        res = qa.execute(
            {"question": text, "context": ctx_str, "session": session_ctx}
        )

    answer = (
        res.get("answer")
        or res.summary
        or "I couldn't find an answer in the pyCSAMT reference."
    )
    # When running without a key, the answer is the deterministic offline
    # composition — nudge the user that an API key unlocks a fluent reply.
    # (A clarify prompt is already complete — don't tack the nudge on.)
    if res.get("source") != "rag_clarify" and (
        offline or (res.get("source") == "rag_offline")
    ):
        answer = (
            answer
            + "\n\n---\n*Offline answer composed from the pyCSAMT reference."
            " For a fuller, synthesised response, add an API key (Claude,"
            " OpenAI, Gemini, DeepSeek or MiniMax) in **Settings**.*"
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
        "pt_strip": "phase-tensor ellipse strip",
        "pt_strip_grid": "phase-tensor strip grid",
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
                if res.get("hint")
                else ""
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
        "freq_editor": "frequency editing",
        "layered_model": "layered-model build",
        "correction": "data correction",
    }
    where = f" for {label}" if label else ""
    step(f"Running {_labels.get(kind, kind)}...", "running")

    res = ToolAgent().execute(
        {"path": edi_path, "kind": kind, **(params or {})}
    )
    if res.status == "failed":
        step("Analysis failed", "done")
        _update_job(
            jid,
            status="done",
            result=res.summary
            + ("." if not res.summary.endswith(".") else "")
            + f" (dataset{where or ''})",
            steps=_JOBS[jid]["steps"],
            kind=KIND_ERROR,
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

    table = (res.data or {}).get("table_text", "")
    result = res.summary
    if table:
        # a fenced block renders monospaced (the chat markdown has no tables)
        result += "\n\n```\n" + table + "\n```"
    if res.warnings:
        result += "\n\n" + "\n".join(f"⚠ {w}" for w in res.warnings[:3])

    _record_run(
        workflow=kind,
        path=str(edi_path),
        output_dir="",
        status=res.status,
        summary=res.summary,
        n_figures=len(figs),
    )
    step("Done", "done")

    # Mutating tools (freq_editor) hand back an edited survey: cache it and
    # raise the post-processing modal so the user can apply it to the session
    # or export it — the same pathway used by correction workflows.
    _corr = (res.data or {}).get("corrected_sites")
    _postproc = None
    if _corr is not None:
        _CORR_CACHE[jid] = _corr
        _postproc = {"jid": jid, "workflow": kind, "output_dir": ""}

    _update_job(
        jid,
        status="done",
        result=result,
        steps=_JOBS[jid]["steps"],
        figs=figs,
        kind=KIND_WORKFLOW if figs else KIND_ANSWER,
        **({"postproc": _postproc} if _postproc else {}),
    )


def _resolve_metric_targets(
    text: str,
    edi_store: dict,
    settings: dict,
    all_lines: bool,
) -> list[tuple[str, Any]]:
    """Resolve which line(s) a metric question is about.

    Returns ``[(label, src), …]`` where ``src`` is an EDI path, a list of EDI
    files, or a directory the MetricsAgent can load. Empty when no data is
    available."""
    import os

    groups = (edi_store or {}).get("groups", {}) or {}
    edi_path = (edi_store or {}).get("path", "") or ""
    sel_lines = (edi_store or {}).get("selected_lines", []) or []

    # An explicit line selection (from the line picker, or a pre-filter) wins
    # for a specific-line question: compute on exactly those line(s) and label
    # them by line name, never the whole upload folder. An explicit "all lines"
    # request is handled below instead.
    if sel_lines and groups and not all_lines:
        picked = [(ln, list(groups[ln])) for ln in sel_lines if groups.get(ln)]
        if picked:
            return picked

    if all_lines:
        if groups:
            return [(ln, list(fs)) for ln, fs in groups.items() if fs]
        if edi_path:
            return [("the survey", edi_path)]
        return []

    # A specific line named among the loaded groups.
    if groups:
        ref = _extract_line_ref(text, groups)
        if ref:
            m = _match_group(ref, groups)
            if m and groups.get(m):
                return [(m, list(groups[m]))]

    # A known survey line from the project registry ("strike of L22PLT").
    try:
        from pycsamt.assistant.tools.project_registry import (
            ProjectRegistry,
        )

        reg = ProjectRegistry.from_default()
        if reg is not None:
            ln = reg.find_line_in_text(text)
            if ln:
                info = reg.resolve_line(ln)
                if info.get("exists"):
                    return [(ln, info["edi_dir"])]
    except Exception:  # noqa: BLE001
        pass

    # The whole loaded dataset.
    if edi_path:
        label = os.path.basename(str(edi_path).rstrip("/\\")) or "the survey"
        return [(label, edi_path)]

    # Session fallback (a follow-up that inherits the last dataset).
    _sess = _session()
    if _sess is not None and getattr(_sess, "edi_path", None):
        return [(getattr(_sess, "line", "") or "the survey", _sess.edi_path)]
    return []


def _dispatch_metrics(
    jid: str,
    text: str,
    edi_store: dict,
    settings: dict,
    *,
    step,
) -> None:
    """Answer a question about computed line value(s) inline via MetricsAgent."""
    from pycsamt.agents.metrics import (
        MetricsAgent,
        parse_metric_request,
    )

    kinds, all_lines = parse_metric_request(text)
    if not kinds:
        kinds = ["summary"]

    step("Resolving line(s)...", "done")
    targets = _resolve_metric_targets(text, edi_store, settings, all_lines)
    if not targets:
        _update_job(
            jid,
            status="done",
            result=(
                "I don't have any survey data to read values from yet. "
                "Load an EDI dataset with **Load EDI** (top-left), or name a "
                "known survey line, then ask again."
            ),
            steps=_JOBS[jid]["steps"],
            kind=KIND_META,
        )
        return

    step("Computing values...", "running")
    warnings: list[str] = []
    if len(targets) == 1:
        label, src = targets[0]
        res = MetricsAgent().execute(
            {"sites": src, "kinds": kinds, "label": label}
        )
        result_text = res.summary
        warnings = list(res.warnings or [])
    else:
        # All lines: one compact line per survey line.
        out_lines = []
        for label, src in targets:
            res = MetricsAgent().execute(
                {"sites": src, "kinds": kinds, "label": label}
            )
            if res.status != "success":
                out_lines.append(f"- **{label}**: {res.summary}")
                continue
            vals = (res.data or {}).get("values", {})
            if len(kinds) == 1:
                out_lines.append(f"- **{label}**: {vals.get(kinds[0], 'n/a')}")
            else:
                sub = "; ".join(f"{k}: {vals.get(k, 'n/a')}" for k in kinds)
                out_lines.append(f"- **{label}**: {sub}")
        header = (
            f"Here's the {kinds[0]} for each line:"
            if len(kinds) == 1
            else "Here's what I found per line:"
        )
        result_text = header + "\n" + "\n".join(out_lines)

    if warnings:
        result_text += "\n\n" + "\n".join(f"⚠ {w}" for w in warnings[:3])

    _record_run(
        workflow="metrics",
        path="",
        output_dir="",
        status="success",
        summary=result_text[:200],
        n_figures=0,
    )
    step("Done", "done")
    _update_job(
        jid,
        status="done",
        result=result_text,
        steps=_JOBS[jid]["steps"],
        kind=KIND_ANSWER,
    )


# ── inversion preparation: per-line file-set builds ──────────────────────
# The deliverable of pre_inversion / modem / mare2dem is a solver input
# file set on disk.  Unlike the generic pipeline path (which merges all
# selected lines into one dataset), each survey line is built separately —
# a 2-D/2.5-D profile inversion is defined per line — and the reply lists
# the written artefacts per line, grounded in the package reference (RAG).

_PREP_STEP_NAME: dict[str, str] = {
    "pre_inversion": "occam2d",
    "modem": "modem",
    "mare2dem": "mare2dem",
}

_PREP_CODE_LABEL: dict[str, str] = {
    "pre_inversion": "Occam2D",
    "modem": "ModEM",
    "mare2dem": "MARE2DEM",
}

# inv-config (param modal) keys each prep agent reads from its input.
_PREP_PARAM_KEYS: dict[str, tuple] = {
    "pre_inversion": ("modes", "period_range", "title"),
    "modem": ("error_floor", "period_range", "component_types"),
    "mare2dem": (
        "error_floor",
        "output_modes",
        "initial_rho",
        "target_rms",
        "max_iterations",
        "topo",
    ),
}

# Output data keys that hold written file paths in the prep AgentResult.
_PREP_FILE_KEYS = (
    "data_path",
    "mesh_path",
    "model_path",
    "startup_path",
    "cov_path",
    "ctrl_path",
    "resistivity_path",
    "settings_path",
)

_PREP_NEXT_STEPS: dict[str, str] = {
    "pre_inversion": (
        "`OccamDataFile.dat` holds the observed responses, `Occam2DMesh` "
        "the finite-element mesh, `Occam2DModel` the mesh→parameter "
        "mapping, and `OccamStartup` the iteration-zero controls. Run the "
        "folder with `pycsamt.models.occam2d.OccamRunner` (or the "
        "Occam2DMT binary directly)."
    ),
    "modem": (
        "`ModEM_Data.dat` holds the impedances, `m0.ws` the half-space "
        "starting model, `ModEM.cov` the model covariance, and `ModEM.inv` "
        "the NLCG inversion controls. Launch the run with "
        "`pycsamt.models.modem.ModEmRunner` or the ModEM binary."
    ),
    "mare2dem": (
        "`mare2dem.emdata` holds the TE/TM responses projected onto the "
        "profile, `mare2dem.resistivity` the starting half-space model, "
        "and `mare2dem.settings` the parallel decomposition. Run them "
        "with `pycsamt.models.mare2dem.Mare2DEMRunner` — the MARE2DEM "
        "binary must be compiled once per machine via `SourceManager`."
    ),
}


def _prep_all_lines_requested(text: str) -> bool:
    """True when the request names all survey lines explicitly."""
    return bool(
        re.search(
            r"\b(?:all|every|each)\s+(?:the\s+)?(?:survey\s+)?"
            r"(?:lines?|profiles?)\b",
            (text or "").lower(),
        )
    )


def _safe_dirname(label: str) -> str:
    """Filesystem-safe per-line folder name."""
    safe = "".join(c if c.isalnum() or c in "-_ " else "_" for c in str(label))
    return safe.strip().replace(" ", "_") or "line"


def _fmt_file_size(path) -> str:
    import os

    try:
        kb = os.path.getsize(str(path)) / 1024.0
    except OSError:
        return ""
    return f"{kb / 1024.0:.1f} MB" if kb >= 1024 else f"{max(kb, 1):.0f} KB"


def _collect_prep_files(result, prep_step: str):
    """Extract written files + stats from one orchestrator run.

    Returns ``(files, stats, warns, step_results)`` where *files* is a
    list of existing paths written by the prep step."""
    import os

    exec_res = (result.data or {}).get("result")
    step_results = (
        exec_res.data if exec_res is not None and exec_res.data else {}
    )
    sres = step_results.get(prep_step)
    d = getattr(sres, "data", None) or {}
    files = []
    for key in _PREP_FILE_KEYS:
        p = d.get(key)
        if p is None:
            continue
        try:
            if os.path.exists(str(p)):
                files.append(str(p))
        except (OSError, ValueError):
            continue
    stats = {
        k: d.get(k) for k in ("n_stations", "n_periods", "n_data") if d.get(k)
    }
    warns = list(getattr(sres, "warnings", None) or [])
    return files, stats, warns, step_results


def _prep_rag_note(query: str) -> str:
    """Best-effort RAG grounding line for the prep reply.

    Retrieves the package symbols behind the build so the closing
    guidance cites real, verifiable APIs. Degrades silently to an empty
    string when the RAG index is unavailable."""
    try:
        from pycsamt.assistant.rag.context_builder import (
            default_context_builder,
        )

        builder = default_context_builder()
        if builder is None:
            return ""
        ctx_r = builder.build(query)
        if ctx_r.is_empty() or not ctx_r.citations:
            return ""
        seen: set[str] = set()
        syms: list[str] = []
        for c in ctx_r.citations:
            s = c.get("symbol") or c.get("source_path")
            if s and s not in seen:
                seen.add(s)
                syms.append(str(s))
        if not syms:
            return ""
        return (
            "*Grounded in the pyCSAMT reference: "
            + " · ".join(f"`{s}`" for s in syms[:4])
            + "*"
        )
    except Exception:  # noqa: BLE001 — RAG is best-effort
        return ""


def _dispatch_inversion_prep(
    jid: str,
    wtype: str,
    text: str,
    edi_store: dict,
    settings: dict,
    cfg: dict,
    *,
    inv_config: dict,
    llm_prov: str,
    api_key: str | None,
    sel_model: str | None,
    offline: bool,
    step,
) -> None:
    """Build the inversion input file set, one orchestrator run per line."""
    import os

    from pycsamt.agents.orchestrator import (
        WorkflowOrchestratorAgent,
    )
    from pycsamt.api.agents import AGENT_CONFIG

    code = _PREP_CODE_LABEL.get(wtype, wtype)
    prep_step = _PREP_STEP_NAME[wtype]

    step("Resolving survey line(s)...", "done")
    targets = _resolve_metric_targets(
        text, edi_store, settings, _prep_all_lines_requested(text)
    )
    if not targets:
        _update_job(
            jid,
            status="done",
            result=(
                "I don't have survey data to prepare inversion files "
                "from yet. Load an EDI dataset with **Load EDI** "
                "(top-left), or name a known survey line, then ask again."
            ),
            steps=_JOBS[jid]["steps"],
            kind=KIND_META,
        )
        return

    root = (
        (settings or {}).get("output_dir") or ""
    ).strip() or "pycsamt_inversion_prep"

    # Forward the parameters the prep agent understands (from the param
    # modal) via the orchestrator's per-step injection.
    _ic = inv_config or {}
    prep_params = {
        k: _ic[k]
        for k in _PREP_PARAM_KEYS.get(wtype, ())
        if _ic.get(k) not in (None, "", [])
    }
    if prep_params:
        cfg.setdefault("step_params", {}).setdefault(prep_step, {}).update(
            prep_params
        )

    multi = len(targets) > 1
    figs: dict = {}
    sections: list[str] = []
    n_ok = 0

    for label, src in targets:
        step(f"Building {code} files for {label}...", "running")
        out_dir = os.path.join(root, _safe_dirname(label)) if multi else root
        try:
            with AGENT_CONFIG.offline() if offline else _nullctx():
                orch = WorkflowOrchestratorAgent(
                    llm_provider=llm_prov,
                    api_key=api_key,
                    model=sel_model,
                )
                result = orch.execute(
                    {
                        "config": dict(cfg),
                        "request": text,
                        "data_path": src,
                        "output_dir": out_dir,
                    }
                )
        except Exception as exc:  # noqa: BLE001
            sections.append(f"**{label}** — ⚠ {exc}")
            step(f"{label} failed", "error")
            continue

        ok = result.status != "failed"
        files, stats, warns, step_results = _collect_prep_files(
            result, prep_step
        )

        # keep only the prep/report step figures, labelled per line
        _fig_steps = _WORKFLOW_FIGURE_STEPS.get(wtype) or set()
        for sname, sres in step_results.items():
            sfigs = (getattr(sres, "data", None) or {}).get(
                "figures", {}
            ) or {}
            for fname, fig in sfigs.items():
                if not isinstance(fig, plt.Figure):
                    continue
                if ok and sname in _fig_steps:
                    figs[str(uuid.uuid4())] = {
                        "title": f"[{label} · {sname}] {fname}",
                        "b64": _fig_to_b64(fig),
                    }
                plt.close(fig)

        if ok and files:
            n_ok += 1
            head = f"**{label}**"
            if stats.get("n_stations"):
                head += f" — {stats['n_stations']} station(s)"
            if stats.get("n_periods"):
                head += f" × {stats['n_periods']} period(s)"
            rows = [
                f"- `{os.path.basename(p)}`"
                + (f" ({_fmt_file_size(p)})" if _fmt_file_size(p) else "")
                + f" — `{p}`"
                for p in files
            ]
            sec = head + "\n" + "\n".join(rows)
        else:
            err = (
                result.error
                or "; ".join(warns[:1])
                or result.summary
                or "no files were written"
            )
            sec = f"**{label}** — ⚠ {err}"
        for w in warns[:2]:
            sec += f"\n  - ⚠ {w}"
        sections.append(sec)
        step(f"{label} done", "done" if ok and files else "error")

    if multi:
        header = (
            f"{code} inversion preparation — {n_ok}/{len(targets)} "
            f"line(s) built under `{root}`."
        )
    elif n_ok:
        header = f"{code} inversion preparation complete."
    else:
        header = f"{code} inversion preparation failed."
    parts = [header, "\n\n".join(sections)]
    guide = _PREP_NEXT_STEPS.get(wtype, "")
    if n_ok and guide:
        parts.append(f"**Next steps.** {guide}")
    rag_note = _prep_rag_note(f"{code} inversion input files runner prepare")
    if n_ok and rag_note:
        parts.append(rag_note)

    _record_run(
        workflow=wtype,
        path=", ".join(str(lbl) for lbl, _ in targets),
        output_dir=root,
        status="success" if n_ok else "failed",
        summary=header,
        n_figures=len(figs),
    )
    _update_job(
        jid,
        status="done",
        result="\n\n".join(parts),
        steps=_JOBS[jid]["steps"],
        figs=figs,
        kind=KIND_WORKFLOW if n_ok else KIND_ERROR,
    )


# ── data overview: "read the EDI data / stations / sites" ───────────────
# Deterministic, offline-first task: read the survey already stored in the
# session, compute statistics, and answer with a rich overview card — no
# LLM call and no pipeline run.
_DATA_READ_VERBS = (
    "read",
    "show",
    "list",
    "describe",
    "summar",
    "inspect",
    "check",
    "overview",
    "stat",
)
_DATA_READ_NOUNS = (
    "edi",
    "data",
    "dataset",
    "station",
    "site",
    "survey",
    "line",
)
# Phrasing that belongs to another intent (plots, workflows, code, docs).
_DATA_READ_BLOCK = (
    "plot",
    "map",
    "pseudosection",
    "pseudo-section",
    "section",
    "figure",
    "chart",
    "rose",
    "diagram",
    "code",
    "script",
    "run ",
    "invert",
    "inversion",
    "correct",
    "filter",
    "denoise",
    "export",
    "save",
    "write",
    "how do i",
    "how to",
    "what is",
    "what does",
    "why",
)


_LINES_QUERY_PATTERNS = (
    "what are the lines",
    "what lines",
    "which lines",
    "list lines",
    "list the lines",
    "available lines",
    "loaded lines",
    "how many lines",
    "the lines?",
    "what are the profiles",
    "which profiles",
    "list profiles",
    "loaded profiles",
    "how many profiles",
)


def _looks_like_lines_query(text: str) -> bool:
    """True when *text* asks which survey lines are loaded."""
    t = (text or "").lower()
    return any(p in t for p in _LINES_QUERY_PATTERNS)


_NO_DATA_GUIDANCE = (
    "**No survey data is stored yet.**\n\n"
    "Load a dataset first and I'll read it and report the "
    "statistics:\n\n"
    "- **Load EDI** — use the ⊕ menu (top-left) and "
    "pick your EDI folder;\n"
    "- or name a known survey line, e.g. "
    "*“read line L22PLT”*.\n\n"
    "Once loaded, ask me to *“read the data”* again."
)


def _looks_like_data_read(text: str) -> bool:
    """True when *text* asks to read/summarise the loaded survey data."""
    t = (text or "").lower()
    if _looks_like_lines_query(text):
        return True
    if any(b in t for b in _DATA_READ_BLOCK):
        return False
    if not any(n in t for n in _DATA_READ_NOUNS):
        return False
    if any(v in t for v in _DATA_READ_VERBS):
        return True
    return any(
        p in t
        for p in (
            "what data",
            "which data",
            "what's loaded",
            "whats loaded",
        )
    )


def _line_stats(label: str, sites, warnings: list[str]) -> dict:
    """Per-line statistics shown on the overview card."""
    import numpy as np

    from pycsamt.agents.metrics import (
        MetricsAgent,
        _has_tipper,
        _ll_to_utm,
        _station_coords,
    )

    rows = MetricsAgent._scan(sites)
    names = [str(r.get("station", "?")) for r in rows]
    tmins = [r["t_min_s"] for r in rows if r.get("t_min_s")]
    tmaxs = [r["t_max_s"] for r in rows if r.get("t_max_s")]
    nfr = [r["n_freq"] for r in rows if r.get("n_freq")]

    freq = None
    if tmins and tmaxs:
        freq = (1.0 / max(tmaxs), 1.0 / min(tmins))

    scores = [
        float(r["qc_score"]) for r in rows if r.get("qc_score") is not None
    ]
    qc = float(np.nanmean(scores)) if scores else None
    flagged = sum(
        1
        for r in rows
        if not r.get("has_z")
        or not r.get("has_coords")
        or (r.get("qc_score") or 0) < 50
    )
    if flagged:
        warnings.append(
            f"{label}: {flagged} of {len(rows)} station(s) flagged "
            "(low QC, missing Z, or no coordinates)"
        )

    length_km = None
    try:
        coords = [
            (la, lo)
            for _, la, lo in _station_coords(sites)
            if la is not None and lo is not None
        ]
        if len(coords) >= 2:
            es, ns = [], []
            for la, lo in coords:
                e, n, _ = _ll_to_utm(la, lo, None, "N", "WGS84")
                es.append(e)
                ns.append(n)
            es, ns = np.asarray(es), np.asarray(ns)
            length_km = (
                float(np.hypot(es.max() - es.min(), ns.max() - ns.min()))
                / 1000.0
            )
    except Exception:  # noqa: BLE001
        pass

    tipper = False
    try:
        tipper = bool(_has_tipper(sites))
    except Exception:  # noqa: BLE001
        pass

    return {
        "label": label,
        "n_stations": len(rows),
        "stations": names,
        "freq": freq,
        "max_nfreq": max(nfr) if nfr else None,
        "qc": qc,
        "flagged": flagged,
        "length_km": length_km,
        "tipper": tipper,
    }


def _fmt_hz(v: float) -> str:
    """Human frequency: 0.125 Hz, 384 Hz, 8.19 kHz — never 8.19e+03."""
    if v >= 1000.0:
        return f"{v / 1000.0:.3g} kHz"
    return f"{v:.3g} Hz"


def _fmt_freq(freq) -> str:
    if not freq:
        return "n/a"
    return f"{_fmt_hz(freq[0])} – {_fmt_hz(freq[1])}"


def _overview_payload(lines: list[dict], warnings: list[str]) -> dict:
    """JSON-safe payload for the survey overview card."""
    import numpy as np

    n_st = sum(ln["n_stations"] for ln in lines)
    flo = [ln["freq"][0] for ln in lines if ln.get("freq")]
    fhi = [ln["freq"][1] for ln in lines if ln.get("freq")]
    qcs = [ln["qc"] for ln in lines if ln.get("qc") is not None]
    lens = [ln["length_km"] for ln in lines if ln.get("length_km")]
    tiles = [
        {"value": f"{n_st}", "label": "stations"},
        {
            "value": (_fmt_freq((min(flo), max(fhi))) if flo else "n/a"),
            "label": "frequency range",
        },
        {
            "value": (f"≈ {sum(lens):.1f} km" if lens else "n/a"),
            "label": ("profile length" if len(lines) == 1 else "total length"),
        },
        {
            "value": (f"{float(np.mean(qcs)):.0f}/100" if qcs else "n/a"),
            "label": "mean QC score",
        },
    ]
    scope_bits = [f"{len(lines)} line{'s' if len(lines) != 1 else ''}"]
    if any(ln.get("tipper") for ln in lines):
        scope_bits.append("tipper present")
    return {
        "tiles": tiles,
        "scope": " · ".join(scope_bits),
        "lines": [
            {
                "label": ln["label"],
                "n_stations": ln["n_stations"],
                "freq": _fmt_freq(ln.get("freq")),
                "qc": (
                    f"{ln['qc']:.0f}/100"
                    if ln.get("qc") is not None
                    else "n/a"
                ),
                "flagged": ln.get("flagged", 0),
            }
            for ln in lines
        ],
        "stations": (lines[0]["stations"] if len(lines) == 1 else []),
        "warnings": warnings,
    }


def _data_overview_card(card: dict) -> html.Div:
    """Rich 'survey at a glance' card rendered under the reply text."""
    children: list = [
        html.Div(
            [
                html.I(className="bi bi-database-check"),
                html.Span("Survey at a glance"),
                html.Span(
                    card.get("scope", ""),
                    className="am-ov-scope",
                ),
            ],
            className="am-ov-head",
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.B(t.get("value", "")),
                        html.Span(t.get("label", "")),
                    ],
                    className="am-ov-tile",
                )
                for t in card.get("tiles", [])
            ],
            className="am-ov-tiles",
        ),
    ]

    rows = card.get("lines", [])
    if len(rows) > 1:
        children.append(
            html.Table(
                [
                    html.Thead(
                        html.Tr(
                            [
                                html.Th("Line"),
                                html.Th("Stations"),
                                html.Th("Frequency range"),
                                html.Th("QC"),
                            ]
                        )
                    ),
                    html.Tbody(
                        [
                            html.Tr(
                                [
                                    html.Td(r["label"]),
                                    html.Td(str(r["n_stations"])),
                                    html.Td(r["freq"]),
                                    html.Td(
                                        r["qc"]
                                        + (
                                            f" · {r['flagged']} flagged"
                                            if r.get("flagged")
                                            else ""
                                        )
                                    ),
                                ]
                            )
                            for r in rows
                        ]
                    ),
                ],
                className="am-ov-table",
            )
        )

    chips = card.get("stations", [])
    if chips:
        shown = chips[:12]
        extra = len(chips) - len(shown)
        children.append(
            html.Div(
                [html.Span(c, className="am-ov-chip") for c in shown]
                + (
                    [
                        html.Span(
                            f"+{extra} more",
                            className="am-ov-chip more",
                        )
                    ]
                    if extra > 0
                    else []
                ),
                className="am-ov-chips",
            )
        )

    for w in card.get("warnings", [])[:3]:
        children.append(
            html.Div(
                [
                    html.I(className="bi bi-exclamation-triangle"),
                    html.Span(w),
                ],
                className="am-ov-warn",
            )
        )

    children.append(
        html.Div(
            [
                html.I(className="bi bi-lightbulb"),
                html.Span(
                    "Next: try “run qc”, “plot a "
                    "station map”, or “prepare an "
                    "Occam2D inversion”.",
                ),
            ],
            className="am-ov-hint",
        )
    )
    return html.Div(children, className="am-ov-card")


def _dispatch_data_overview(
    jid: str,
    text: str,
    edi_store: dict,
    settings: dict,
    *,
    step,
) -> None:
    """Read the stored survey inline and reply with statistics."""
    _update_job(jid, workflow="data_overview")
    step("Checking stored survey data...", "done")

    groups = (edi_store or {}).get("groups", {}) or {}

    # "what are the lines?" — answer instantly from the store, no
    # file reading needed.
    if _looks_like_lines_query(text):
        step("Listing survey lines...", "done")
        if groups:
            n = len(groups)
            keys = sorted(groups)
            rows = "\n".join(
                f"- **{ln}** — {len(groups[ln])} EDI file"
                f"{'s' if len(groups[ln]) != 1 else ''}"
                for ln in keys
            )
            result = (
                f"**{n} survey line"
                f"{'s are' if n != 1 else ' is'} loaded:**\n\n"
                f"{rows}\n\n"
                f"Say *“read line {keys[0]}”* for statistics, or "
                f"*“plot the phase pseudosection of {keys[0]}”* "
                "for a figure."
            )
            kind = KIND_ANSWER
        elif (edi_store or {}).get("path"):
            import os as _os

            label = (
                _os.path.basename(str(edi_store["path"]).rstrip("/\\"))
                or "the survey"
            )
            result = (
                f"One dataset is loaded (**{label}**) with no "
                "separate line grouping. Say *“read the data”* "
                "for statistics, or *“plot the sounding "
                "curves”* for a figure."
            )
            kind = KIND_ANSWER
        else:
            result = _NO_DATA_GUIDANCE
            kind = KIND_META
        _update_job(
            jid,
            status="done",
            result=result,
            steps=_JOBS[jid]["steps"],
            kind=kind,
        )
        return

    # A named line reads that line only; otherwise read every line.
    named = bool(groups) and bool(_extract_line_ref(text, groups))
    targets = _resolve_metric_targets(
        text,
        edi_store,
        settings,
        all_lines=not named,
    )
    if not targets:
        _update_job(
            jid,
            status="done",
            result=_NO_DATA_GUIDANCE,
            steps=_JOBS[jid]["steps"],
            kind=KIND_META,
        )
        return

    from pycsamt.emtools._core import ensure_sites

    lines: list[dict] = []
    warnings: list[str] = []
    for label, src in targets:
        step(f"Reading {label}...", "running")
        try:
            sites = ensure_sites(
                src,
                recursive=True,
                strict=False,
                verbose=0,
            )
            lines.append(_line_stats(label, sites, warnings))
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"{label}: could not be read ({exc})")

    if not lines:
        _update_job(
            jid,
            status="done",
            result=(
                "I found stored survey entries but could not read "
                "any of them.\n\n" + "\n".join(f"- {w}" for w in warnings[:4])
            ),
            steps=_JOBS[jid]["steps"],
            kind=KIND_ERROR,
        )
        return

    step("Computing statistics...", "running")
    card = _overview_payload(lines, warnings)
    n_st = sum(ln["n_stations"] for ln in lines)
    if len(lines) == 1:
        lead = (
            f"Read **{lines[0]['label']}** — "
            f"{n_st} station{'s' if n_st != 1 else ''}. "
            "Here's the survey at a glance."
        )
    else:
        lead = (
            f"Read **{n_st} station"
            f"{'s' if n_st != 1 else ''}** across "
            f"**{len(lines)} lines** from the stored survey. "
            "Here's the data at a glance."
        )

    _record_run(
        workflow="data_overview",
        path="",
        output_dir="",
        status="success",
        summary=lead[:200],
        n_figures=0,
    )
    step("Overview ready", "done")
    _update_job(
        jid,
        status="done",
        result=lead,
        card=card,
        steps=_JOBS[jid]["steps"],
        kind=KIND_ANSWER,
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
    from pycsamt.agents.code_gen import CodeGenerationAgent
    from pycsamt.agents.context import ContextInputAgent
    from pycsamt.api.agents import AGENT_CONFIG

    step("Extracting configuration...", "done")
    with AGENT_CONFIG.offline() if offline else _nullctx():
        ctx_agent = ContextInputAgent(
            llm_provider=llm_prov,
            api_key=api_key,
            model=sel_model,
        )
        ctx_res = ctx_agent.execute({"request": text})
    cfg = ctx_res.data.get("config", {}) if ctx_res and ctx_res.data else {}
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
    with AGENT_CONFIG.offline() if offline else _nullctx():
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

    _line_note = f" for line {resolved_line}" if resolved_line else ""
    summary = (
        "Here is a standalone pyCSAMT script that"
        f" reproduces the {cfg.get('workflow', 'qc')}"
        f" workflow{_line_note}. Copy it from the code block"
        " below — edit the data path if needed." + _valid_note
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
            _JOBS[jid]["steps"].append({"label": label, "status": status})

    try:
        _step("Parsing request...", "done")

        # configure provider + api key
        provider = settings.get("provider", "offline")
        key_map = {
            "claude": "ANTHROPIC_API_KEY",
            "openai": "OPENAI_API_KEY",
            "gemini": "GOOGLE_API_KEY",
            "deepseek": "DEEPSEEK_API_KEY",
            "minimax": "MINIMAX_API_KEY",
        }
        api_key: str | None = None
        if provider in key_map:
            api_key = settings.get(f"key_{provider}", "") or None
            if api_key:
                import os

                os.environ[key_map[provider]] = api_key

        # "offline" maps to "claude" provider
        # name so BaseAgent validates it, but
        # api_key stays None → regex fallback.
        llm_prov = provider if provider != "offline" else "claude"
        sel_model = settings.get(f"model_{llm_prov}") or None

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
            CLARIFY as _I_CLARIFY,
        )
        from pycsamt.agents.router import (
            CODE as _I_CODE,
        )
        from pycsamt.agents.router import (
            META as _I_META,
        )
        from pycsamt.agents.router import (
            METRICS as _I_METRICS,
        )
        from pycsamt.agents.router import (
            QUESTION as _I_QUESTION,
        )
        from pycsamt.agents.router import (
            IntentRouter,
        )

        # Deterministic data-overview gate: "read the EDI data /
        # stations / sites" is answered inline from the stored
        # survey — no LLM router call, no pipeline.
        if _looks_like_data_read(text):
            _step("Intent: data overview", "done")
            _dispatch_data_overview(
                jid,
                text,
                edi_store,
                settings,
                step=_step,
            )
            return

        with AGENT_CONFIG.offline() if _offline else _nullctx():
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
        if decision.intent == _I_METRICS:
            _dispatch_metrics(
                jid,
                text,
                edi_store,
                settings,
                step=_step,
            )
            return
        if decision.intent == _I_QUESTION:
            _dispatch_question(
                jid,
                text,
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
                jid,
                text,
                edi_store,
                settings,
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

        with AGENT_CONFIG.offline() if _offline else _nullctx():
            ctx_agent = ContextInputAgent(
                llm_provider=llm_prov,
                api_key=api_key,
                model=sel_model,
            )
            ctx_result = ctx_agent.execute({"request": text})

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
        from pycsamt.agents._workflows import (
            classify_workflow as _cwf,
        )
        from pycsamt.agents.orchestrator import (
            _WORKFLOW_STEPS as _WF_STEPS,
        )

        _explicit_wf = _ic.get("workflow")
        _kw_wf = _cwf(text, default=None)
        _router_wf = (
            decision.workflow
            if (
                decision.workflow in _WF_STEPS
                or decision.workflow in _PLOT_WORKFLOWS
            )
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
        # A generic "analyse …" that names no concrete workflow (e.g.
        # "analyse only the line 22", typically straight after the line
        # picker) runs the flagship phase-tensor analysis on the
        # loaded/selected data instead of an "unknown task" bounce.
        if (
            _resolved_wf is None
            and _GENERIC_ANALYSIS_RE.search(text or "")
            and (
                (edi_store or {}).get("path")
                or (edi_store or {}).get("groups")
            )
        ):
            _resolved_wf = "phase_analysis"
        _known_wf = (
            _resolved_wf in _WF_STEPS
            or _resolved_wf in _PLOT_WORKFLOWS
            or _resolved_wf in _TOOL_WORKFLOWS
        )
        if not _resolved_wf or not _known_wf:
            # Try a context-aware reply first (unknown line →
            # propose the loaded lines; plot request without a
            # plot kind → propose the plot menu).
            _smart = _smart_unknown_reply(text, edi_store)
            _update_job(
                jid,
                status="done",
                result=_smart or _unknown_task_text(text),
                steps=_JOBS[jid]["steps"],
                kind=KIND_CLARIFY if _smart else KIND_META,
            )
            return

        cfg["workflow"] = _resolved_wf
        wtype = _resolved_wf
        # The Pre-Inversion form lets the user pick the target code
        # (Occam2D / ModEM / MARE2DEM); honour it by switching to that
        # code's dedicated chain.
        if wtype == "pre_inversion":
            _code = str(
                _ic.get("inversion_code") or cfg.get("inversion_code") or ""
            ).lower()
            if _code in ("modem", "mare2dem"):
                wtype = _code
                cfg["workflow"] = _code
        _update_job(jid, workflow=wtype)
        _step(f"Workflow: {wtype}", "done")
        if wtype in ("pinn_inversion", "hybrid_inversion"):
            _pi = {
                "dim": _ic.get("dim", 1),
                "n_layers": _ic.get("n_layers", 10),
                "depth_max": _ic.get("depth_max", 2000.0),
                "epochs": _ic.get("epochs", 500),
                "lr": _ic.get("lr", 0.01),
                "smoothness_weight": _ic.get("smoothness_weight", 0.01),
                "lateral_weight": _ic.get("lateral_weight", 0.005),
                "graph_weight": _ic.get("graph_weight", 0.005),
                "radius": _ic.get("radius", 5000.0),
                "solver": _ic.get("solver", "mt1d"),
            }
            cfg["pinn_params"] = _pi
            cfg["hybrid_params"] = _pi
            cfg["checkpoint"] = _ic.get("checkpoint", "")
        elif wtype in (
            "ai_inversion",
            "inv1d",
            "inv2d",
            "inv3d",
            "ensemble_inversion",
        ):
            cfg["ai_inv_params"] = {
                "n_layers": int(_ic.get("n_layers", 10)),
                "depth_max": float(_ic.get("depth_max", 2000.0)),
                "epochs": int(_ic.get("epochs", 500)),
                "lr": float(_ic.get("lr", 0.01)),
                "lateral_weight": float(_ic.get("lateral_weight", 0.005)),
                "graph_weight": float(_ic.get("graph_weight", 0.005)),
                "radius": float(_ic.get("radius", 5000.0)),
            }
            cfg["checkpoint"] = _ic.get("checkpoint", "")

        # Pass pipeline step params into cfg
        _step_p = _ic.get("step_params")
        if _step_p:
            cfg["step_params"] = _step_p

        # build orchestrator input_data
        edi_path = (edi_store.get("path", "") if edi_store else "") or cfg.get(
            "data_path", ""
        )

        # Filter to selected lines if set
        sel_lines = (edi_store or {}).get("selected_lines", [])
        if sel_lines:
            grp = (edi_store or {}).get("groups", {})
            file_list: list[str] = []
            for ln in sel_lines:
                file_list.extend(grp.get(ln, []))
            if file_list:
                edi_path = file_list

        # Fall back to YAML line registry when no
        # EDI is loaded but user names a survey line.
        if not edi_path:
            _reg_yaml = (settings or {}).get("line_registry", "")
            if _reg_yaml:
                try:
                    import yaml as _yaml

                    _reg = _yaml.safe_load(_reg_yaml) or {}
                    _tl = text.lower()
                    for _ln, _lp in _reg.items():
                        if str(_ln).lower() in _tl:
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
        _EDI_REQUIRED = frozenset(
            {
                "qc",
                "static_shift",
                "phase_analysis",
                "ai_inversion",
                "inv1d",
                "inv2d",
                "inv3d",
                "pinn_inversion",
                "hybrid_inversion",
                "ensemble_inversion",
                "pre_inversion",
                "modem",
                "mare2dem",
                "tipper",
                "rotation",
                "denoise",
                "sensitivity",
                "freq_decimation",
                "rhophi",
                "phase_psection",
                "pt_psection",
                "tipper_plot",
                "phase_tensor_map",
                "station_response",
                "strike_profile",
                "pt_strip",
                "pt_strip_grid",
                "strike",
                "dimensionality",
                "validator",
                "coords",
                "elevation",
                "converter",
                "batch_export",
                "freq_editor",  # layered_model is synthetic → no data needed
            }
        ) | frozenset(_CORR_METHODS)  # all correction methods need data
        if wtype in _EDI_REQUIRED and not edi_path:
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

        # ── inversion preparation → per-line file-set builds ──────────
        if wtype in _PREP_WORKFLOWS:
            _dispatch_inversion_prep(
                jid,
                wtype,
                text,
                edi_store,
                settings,
                cfg,
                inv_config=(inv_config or {}),
                llm_prov=llm_prov,
                api_key=api_key,
                sel_model=sel_model,
                offline=_offline,
                step=_step,
            )
            return

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
            _tool_params = dict(inv_config or {})
            # Correction workflows share one ToolAgent kind; tell the agent
            # which catalogue method to apply (and let it coerce the params).
            if wtype in _CORR_METHODS:
                _tool_params["corr_wf"] = wtype
                _tool_params["fn_name"] = _CORR_METHODS[wtype]["fn"]
            _dispatch_tool(
                jid,
                edi_path,
                kind=_TOOL_KIND[wtype],
                params=_tool_params,
                step=_step,
                label=_task_label,
            )
            return

        orch_input = {
            "config": cfg,
            "request": text,
            "data_path": edi_path,
            "output_dir": output_dir,
        }

        with AGENT_CONFIG.offline() if _offline else _nullctx():
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
            "done" if result.status != "failed" else "error",
        )

        # collect figures and generated code
        # result.data["result"] = coordinator AgentResult
        # coordinator AgentResult.data = {step: AgentResult}
        # each step AgentResult.data may have "figures"
        figs: dict = {}
        generated_code: str = ""
        step_results: dict = {}
        if result.status != "failed":
            exec_res = (result.data or {}).get("result")
            step_results = exec_res.data if exec_res and exec_res.data else {}
            _fig_steps = _WORKFLOW_FIGURE_STEPS.get(wtype)
            # None means "show all steps";
            # a set means only those steps.
            # For unlisted workflows, skip
            # load/denoise by default.
            if _fig_steps is None and wtype not in _WORKFLOW_FIGURE_STEPS:
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
                        _sf = (sres.data or {}).get("figures", {})
                        for _f in (_sf or {}).values():
                            if isinstance(_f, plt.Figure):
                                plt.close(_f)
                        continue
                elif sname in _SKIP_ALWAYS:
                    continue
                step_figs = (sres.data or {}).get("figures", {})
                for fname, fig in (step_figs or {}).items():
                    if isinstance(fig, plt.Figure):
                        b64 = _fig_to_b64(fig)
                        plt.close(fig)
                        figs[str(uuid.uuid4())] = {
                            "title": (f"[{sname}] {fname}"),
                            "b64": b64,
                        }

            # extract generated code (code_gen step)
            code_res = step_results.get("code_gen")
            if code_res and hasattr(code_res, "data"):
                generated_code = (code_res.data or {}).get("code", "") or ""

        # cache corrected sites for post-proc modal
        if result.status != "failed" and wtype in _CORRECTION_WFLOWS:
            for _sn, _sr in step_results.items():
                _d = getattr(_sr, "data", None) or {}
                _corr = _d.get("corrected_sites")
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
        summary = result.summary or result.error or "Workflow completed."

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
            kind=(KIND_ERROR if result.status == "failed" else KIND_WORKFLOW),
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
    return line if len(line) <= limit else line[: limit - 1] + "…"


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
    pins.append(
        {
            "mid": mid,
            "role": msg.get("role", "assistant"),
            "snippet": _pin_snippet(msg.get("content", "")),
            "ts": msg.get("ts", _ts()),
        }
    )
    return pins


def _remove_pin(pins: list | None, mid: str) -> list:
    """Return *pins* without the entry for *mid*."""
    return [p for p in (pins or []) if p.get("mid") != mid]


def _pin_item(pin: dict) -> html.Div:
    """Render one pinned-message row in the sidebar."""
    mid = pin.get("mid", "")
    role = pin.get("role", "assistant")
    icon = "bi-person-fill" if role == "user" else "bi-robot"
    return html.Div(
        [
            html.Button(
                [
                    html.I(className=f"bi {icon} am-pin-role"),
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
    icon = "bi-check-circle-fill" if ok else "bi-exclamation-triangle-fill"
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
        Output(IDs.INTERVAL_POLL, "disabled"),
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
        State(IDs.STORE_INV_CONFIG, "data"),
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
            return _stop_job_response(current_msgs, stored_messages, job_store)

        triggered = ctx.triggered_id
        if isinstance(triggered, dict) and triggered.get("type") == "am-chip":
            idx = triggered["index"]
            from ..layout import _PROMPT_CHIPS

            text = _PROMPT_CHIPS[idx]

        if not text or not text.strip():
            raise PreventUpdate

        text = text.strip()

        _user_mid = _mid()
        new_stored = list(stored_messages or [])
        new_stored.append(
            {
                "role": "user",
                "content": text,
                "ts": _ts(),
                "mid": _user_mid,
            }
        )

        msgs = [
            c
            for c in (current_msgs or [])
            if not (
                isinstance(c, dict)
                and c.get("props", {}).get("id") == IDs.WELCOME
            )
        ]
        msgs.append(_user_bubble(text, mid=_user_mid))

        # ── application launch ────────────────
        # "open the map view / web app / desktop app", or a task
        # too complex for the chat (interactive viz → MapView or
        # the web app). All bypass the EDI guard: the apps load
        # their own data.
        _app_req = _detect_app_request(text)
        if _app_req is not None:
            _app_kind, _reason = _app_req
            if _app_kind == "desktop":
                _ok, _note = _ensure_desktop()
                wb = _launch_bubble(
                    "desktop",
                    reason=_reason,
                    note=_note,
                    ok=_ok,
                )
                _log = (
                    "Launched the pyCSAMT desktop app."
                    if _ok
                    else f"Desktop app launch failed: {_note}"
                )
            elif _app_kind == "mapview":
                _mv_url = _ensure_mapview()
                wb = _launch_bubble(
                    "mapview",
                    url=_mv_url,
                    reason=_reason,
                )
                _log = f"Launched MapView: {_mv_url}"
            else:
                _wa_url = _ensure_web_app()
                wb = _launch_bubble(
                    "web",
                    url=_wa_url,
                    reason=_reason,
                )
                _log = f"Redirecting to web app: {_wa_url}"
            msgs.append(wb)
            new_stored.append(
                {
                    "role": "assistant",
                    "content": _log,
                    "ts": _ts(),
                }
            )
            return (
                msgs,
                no_update,
                True,
                "",
                new_stored,
                {},
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
            NO_DATA_INTENTS,
            classify_intent_offline,
        )

        _qi, _ = classify_intent_offline(text)
        # Data-overview requests skip the EDI guard too: with no data
        # stored, the dispatcher replies with load instructions instead
        # of the terse guard message.
        if _qi in NO_DATA_INTENTS or _looks_like_data_read(text):
            msgs.append(
                _thinking_bubble(
                    [
                        {
                            "label": "Parsing request...",
                            "status": "running",
                        }
                    ]
                )
            )
            jid = _new_job()
            t = threading.Thread(
                target=_run_agent,
                args=(
                    jid,
                    text,
                    dict(edi_store or {}),
                    settings or {},
                    _inv_clean,
                    new_stored,
                ),
                daemon=True,
            )
            t.start()
            return (
                msgs,
                {"jid": jid},
                False,
                "",
                new_stored,
                {},
            )

        # ── guard: no EDI loaded ──────────────
        # Skip the guard when the request names a known survey line —
        # _run_agent resolves it via the project registry, so the user
        # need not load EDIs manually for "run X on line L22PLT".
        edi_path = (edi_store or {}).get("path", "")
        # Synthetic tools (e.g. layered model) need no dataset — let them
        # through the no-EDI guard straight to the param modal.
        from pycsamt.agents._workflows import (
            classify_workflow as _cwf_guard,
        )

        _dataless = _cwf_guard(text) in _NO_DATA_WORKFLOWS
        if (
            not edi_path
            and not _dataless
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
            new_stored.append(
                {
                    "role": "assistant",
                    "content": _no_edi,
                    "ts": _ts(),
                }
            )
            return (
                msgs,
                no_update,
                True,
                "",
                new_stored,
                {},
            )

        # ── line disambiguation ───────────────
        # Detect if the user named a specific line
        # that can't be resolved to a group key.
        # If ambiguous: show the line picker modal.
        # If exact match: pre-filter edi_store.
        _groups = (edi_store or {}).get("groups", {})
        _sel_lines: list[str] = []
        if _groups and len(_groups) > 1:
            _ref = _extract_line_ref(text, _groups)
            if _ref is not None:
                _exact = _match_group(_ref, _groups)
                if _exact is None:
                    # Ambiguous → show picker
                    msgs.append(_line_waiting_bubble())
                    pending = {
                        "disambiguation": "lines",
                        "text": text,
                        "groups": {k: list(v) for k, v in _groups.items()},
                    }
                    return (
                        msgs,
                        {},
                        True,
                        "",
                        new_stored,
                        pending,
                    )
                else:
                    _sel_lines = [_exact]
            elif not _prep_all_lines_requested(text):
                # Inversion-prep builds are per
                # line: with several lines loaded
                # and none named, ask which line(s)
                # to build — "Run all" builds every
                # line. ("all lines" in the text
                # skips the picker.)
                from pycsamt.agents._workflows import (
                    classify_workflow as _cwf_prep,
                )

                if _cwf_prep(text, default=None) in _PREP_WORKFLOWS:
                    msgs.append(_line_waiting_bubble())
                    pending = {
                        "disambiguation": "lines",
                        "text": text,
                        "groups": {k: list(v) for k, v in _groups.items()},
                    }
                    return (
                        msgs,
                        {},
                        True,
                        "",
                        new_stored,
                        pending,
                    )

        # ── param detection ───────────────────
        wf = _quick_workflow(text)
        if wf:
            msgs.append(_waiting_bubble(wf))
            pending = {
                "workflow": wf,
                "text": text,
                # snapshot edi_path so _submit_params
                # can fall back if STORE_EDI is None
                "edi_path": ((edi_store or {}).get("path", "") or ""),
                "edi_groups": ((edi_store or {}).get("groups", {})),
            }
            if _sel_lines:
                pending["selected_lines"] = _sel_lines
            return (
                msgs,
                {},
                True,
                "",
                new_stored,
                pending,
            )

        # ── normal flow: start job ────────────
        msgs.append(
            _thinking_bubble(
                [
                    {
                        "label": "Parsing request...",
                        "status": "running",
                    }
                ]
            )
        )
        _edi_use = dict(edi_store or {})
        if _sel_lines:
            _edi_use["selected_lines"] = _sel_lines
        jid = _new_job()
        t = threading.Thread(
            target=_run_agent,
            args=(
                jid,
                text,
                _edi_use,
                settings or {},
                _inv_clean,
                new_stored,
            ),
            daemon=True,
        )
        t.start()

        return (
            msgs,
            {"jid": jid},
            False,
            "",
            new_stored,
            {},
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
        Output(IDs.STORE_POSTPROC, "data"),
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
                no_update,
                True,
                no_update,
                no_update,
                no_update,
            )

        # replace last thinking bubble with
        # updated thinking bubble
        msgs = list(current_msgs or [])
        thinking_idx = None
        for i, child in enumerate(msgs):
            if isinstance(child, dict):
                cid = child.get("props", {}).get("id", "")
                if cid == "am-thinking-bubble":
                    thinking_idx = i
                    break

        if status == "running":
            _started = job.get("started")
            _elapsed = time.time() - _started if _started else None
            new_thinking = _thinking_bubble(
                steps,
                workflow=job.get("workflow"),
                elapsed=_elapsed,
            )
            if thinking_idx is not None:
                msgs[thinking_idx] = new_thinking
            return (
                msgs,
                False,
                no_update,
                no_update,
                no_update,
            )

        # job done / error
        result_text = job.get("result") or job.get("error") or "Done."
        figs = job.get("figs", {})
        code = job.get("code", "")
        kind = job.get("kind")
        card = job.get("card")

        # merge figs into store
        new_fig_store = dict(fig_store or {})
        new_fig_store.update(figs)

        _agent_mid = _mid()
        agent_bub = _agent_bubble(
            result_text,
            steps,
            figs,
            code=code,
            kind=kind,
            mid=_agent_mid,
            card=card,
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

        new_stored = list(stored_messages or [])
        new_stored.append(
            {
                "role": "assistant",
                "content": result_text,
                "ts": _ts(),
                "mid": _agent_mid,
            }
        )
        return (
            msgs,
            True,
            new_fig_store,
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
        if not fig_store or fig_key not in (fig_store):
            raise PreventUpdate
        info = fig_store[fig_key]
        src = f"data:image/png;base64,{info['b64']}"
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
        if not ctx.triggered or not ctx.triggered[0].get("value"):
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
            IDs.STORE_PINS,
            "data",
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
        if not ctx.triggered or not ctx.triggered[0].get("value"):
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
