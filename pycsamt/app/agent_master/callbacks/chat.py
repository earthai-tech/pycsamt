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

from dash import ALL, Input, Output, State
from dash import ctx, html, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from .._ids import IDs

# ── shared job registry ────────────────────────────
_JOBS: dict[str, dict] = {}
_JOBS_LOCK = threading.Lock()


def _new_job() -> str:
    jid = str(uuid.uuid4())
    with _JOBS_LOCK:
        _JOBS[jid] = {
            "status": "running",
            "steps": [],
            "result": None,
            "figs": {},
            "error": None,
        }
    return jid


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


def _fig_card(
    fig_key: str,
    title: str,
    b64: str,
) -> html.Div:
    return html.Div(
        [
            html.Div(
                [
                    html.Span(
                        title,
                        className="am-fig-title",
                    ),
                    html.Button(
                        [
                            html.I(
                                className=(
                                    "bi bi-arrows-"
                                    "fullscreen"
                                )
                            ),
                            " View",
                        ],
                        className="am-fig-btn",
                        id={
                            "type": "am-fig-open",
                            "key": fig_key,
                        },
                        n_clicks=0,
                    ),
                ],
                className="am-fig-toolbar",
            ),
            html.Div(
                html.Img(
                    src=f"data:image/png;base64,"
                    f"{b64}",
                    style={"cursor": "zoom-in"},
                    id={
                        "type": "am-fig-img",
                        "key": fig_key,
                    },
                ),
                className="am-fig-preview",
            ),
        ],
        className="am-fig-card",
    )


# ── message bubble builders ────────────────────────

def _ts() -> str:
    return datetime.now().strftime("%H:%M")


def _user_bubble(text: str) -> html.Div:
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
                            [
                                html.Button(
                                    html.I(
                                        className=(
                                            "bi "
                                            "bi-clipboard"
                                        )
                                    ),
                                    className=(
                                        "am-msg-action"
                                        " am-copy-btn"
                                    ),
                                    title="Copy",
                                    n_clicks=0,
                                ),
                                html.Button(
                                    html.I(
                                        className=(
                                            "bi "
                                            "bi-folder2"
                                            "-open"
                                        )
                                    ),
                                    className=(
                                        "am-msg-action"
                                        " am-edi-msg-btn"
                                    ),
                                    title=(
                                        "Load EDI"
                                    ),
                                    n_clicks=0,
                                ),
                            ],
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
    )


def _thinking_bubble(steps: list[dict]) -> html.Div:
    step_els = [
        html.Div(
            [
                html.I(
                    className=(
                        "bi bi-check-circle-fill"
                        " am-step-icon"
                        if s["status"] == "done"
                        else (
                            "bi bi-exclamation-triangle"
                            "-fill am-step-icon"
                            if s["status"] == "error"
                            else
                            "bi bi-arrow-repeat"
                            " am-step-icon"
                        )
                    )
                ),
                html.Span(s["label"]),
            ],
            className=(
                f"am-step {s['status']}"
            ),
        )
        for s in steps
    ]
    dots = html.Div(
        [
            html.Div(
                [
                    html.I(
                        className=(
                            "bi bi-robot me-2"
                        )
                    ),
                    html.Span("Thinking"),
                    html.Div(
                        [
                            html.Span(),
                            html.Span(),
                            html.Span(),
                        ],
                        className="am-dots",
                    ),
                ],
                className="am-thinking",
            ),
            html.Div(
                step_els,
                className="am-steps",
            )
            if step_els
            else html.Div(),
        ]
    )
    return html.Div(
        [
            html.Div(
                html.I(className="bi bi-robot"),
                className="am-avatar agent",
            ),
            html.Div(
                dots,
                className="am-bubble agent",
            ),
        ],
        className="am-msg-row",
        id="am-thinking-bubble",
    )


def _agent_bubble(
    text: str,
    steps: list[dict] | None = None,
    figs: dict | None = None,
) -> html.Div:
    children: list = []
    # plain text / markdown (rendered as plain)
    for line in (text or "").split("\n"):
        if line.startswith("**") and line.endswith(
            "**"
        ):
            children.append(
                html.Strong(line[2:-2])
            )
        elif line.startswith("- "):
            children.append(
                html.Li(
                    line[2:],
                    style={"marginLeft": "12px"},
                )
            )
        elif line.startswith("```"):
            pass
        else:
            children.append(html.P(line))

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

    # figure cards
    if figs:
        for fig_key, fig_info in figs.items():
            children.append(
                _fig_card(
                    fig_key,
                    fig_info.get("title", "Figure"),
                    fig_info["b64"],
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
                ],
                className="am-bubble agent",
            ),
        ],
        className="am-msg-row",
    )


# ── Smart param detection ─────────────────────────

# Workflows that require user parameters
# before the job can start.
_NEEDS_PARAMS: frozenset[str] = frozenset({
    # Inversion workflows
    "ai_inversion", "inv1d",
    "inv2d", "inv3d",
    "pinn_inversion", "hybrid_inversion",
    "ensemble_inversion",
    "pre_inversion", "modem",
    # Pipeline-only workflows
    "qc", "phase_analysis",
    "static_shift", "tipper", "rotation",
    # Agent-focused workflows
    "interpret", "interpretation",
    "report",
    "code_gen",
    "denoise",
    "sensitivity",
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


# ── agent runner ───────────────────────────────────

def _run_agent(
    jid: str,
    text: str,
    edi_store: dict,
    settings: dict,
    inv_config: dict | None = None,
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

        _step("Classifying workflow...", "done")

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
            )
            return

        cfg = ctx_result.data.get("config", {})

        # inject workflow-specific params
        _ic = inv_config or {}
        # param-modal workflow takes precedence
        # over re-classification from text alone
        if _ic.get("workflow"):
            cfg["workflow"] = _ic["workflow"]
        wtype = cfg.get("workflow", "qc")
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
        orch_input = {
            "config":    cfg,
            "request":   text,
            "data_path": edi_path,
        }

        orch = WorkflowOrchestratorAgent(
            llm_provider=llm_prov,
            api_key=api_key,
            model=sel_model,
        )
        _step(f"Executing {wtype}...", "running")

        result = orch.execute(orch_input)

        _step(
            f"Completed {wtype}",
            "done"
            if result.status != "failed"
            else "error",
        )

        # collect figures from all workflow steps
        # result.data["result"] = coordinator AgentResult
        # coordinator AgentResult.data = {step: AgentResult}
        # each step AgentResult.data may have "figures"
        figs: dict = {}
        if result.status != "failed":
            exec_res = (
                result.data or {}
            ).get("result")
            step_results = (
                exec_res.data
                if exec_res and exec_res.data
                else {}
            )
            for sname, sres in step_results.items():
                if not hasattr(sres, "data"):
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

        # AgentResult.summary is the text field
        summary = (
            result.summary
            or result.error
            or "Workflow completed."
        )

        _update_job(
            jid,
            status="done",
            result=summary,
            steps=_JOBS[jid]["steps"],
            figs=figs,
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
    ):
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

        new_stored = list(
            stored_messages or []
        )
        new_stored.append({
            "role": "user",
            "content": text,
            "ts": _ts(),
        })

        msgs = [
            c for c in (current_msgs or [])
            if not (
                isinstance(c, dict)
                and c.get("props", {}).get("id")
                == IDs.WELCOME
            )
        ]
        msgs.append(_user_bubble(text))

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

        # ── guard: no EDI loaded ──────────────
        edi_path = (edi_store or {}).get(
            "path", ""
        )
        if not edi_path:
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

        # ── param detection ───────────────────
        wf = _quick_workflow(text)
        if wf:
            msgs.append(_waiting_bubble(wf))
            pending = {
                "workflow": wf, "text": text,
            }
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
        jid = _new_job()
        t = threading.Thread(
            target=_run_agent,
            args=(
                jid, text,
                edi_store or {},
                settings or {},
                inv_config or {},
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
            new_thinking = _thinking_bubble(
                steps
            )
            if thinking_idx is not None:
                msgs[thinking_idx] = new_thinking
            return (
                msgs, False, no_update, no_update
            )

        # job done / error
        result_text = (
            job.get("result")
            or job.get("error")
            or "Done."
        )
        figs = job.get("figs", {})

        # merge figs into store
        new_fig_store = dict(fig_store or {})
        new_fig_store.update(figs)

        agent_bub = _agent_bubble(
            result_text, steps, figs
        )

        # replace thinking with agent bubble
        if thinking_idx is not None:
            msgs[thinking_idx] = agent_bub
        else:
            msgs.append(agent_bub)

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
        })
        return (
            msgs, True, new_fig_store,
            new_stored,
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
