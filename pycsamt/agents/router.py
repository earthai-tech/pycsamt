# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents.router
=====================

:class:`IntentRouter` — the true "master" entry point.

The orchestrator classifies *which geophysical workflow* a request maps to.
But that is only the right question once we already know the user wants to
**run a workflow** at all.  Many messages are not workflow requests:

* a **question** about the package ("what does StaticShiftAgent do?"),
* a **code** request ("write a script to load EDIs"),
* a **plot** request ("show me the pseudosection"),
* a **meta** request ("what can you do?", "hi"),

The :class:`IntentRouter` sits above the workflow classifier and decides the
*intent* first, then lets each intent dispatch to the right specialist agent.

Design
------
* **LLM-first when a key is available** — a single structured tool-style call
  returns ``{intent, workflow, confidence, clarification, reasoning}``.
* **Deterministic offline fallback** — :func:`classify_intent_offline` is a
  fast, pure-function heuristic used when no key is configured *and* as the
  cheap pre-check the GUI uses to decide whether a message even needs an EDI
  dataset loaded.
* **Low-confidence → clarify** — rather than guessing, the router can return
  ``intent == CLARIFY`` so the caller asks one disambiguating question.

The router never executes anything; it only decides *where a message goes*.
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass
from typing import Any

from ._base import AgentResult, BaseAgent

# ── intent constants ───────────────────────────────────────────────────────────

QUESTION = "question"  # answer about pycsamt → PackageQAAgent
CODE = "code"  # emit a script        → CodeGenerationAgent
PLOT = "plot"  # produce a figure     → workflow / web app
WORKFLOW = "workflow"  # run a pipeline       → Context → Orchestrator
META = "meta"  # capabilities / chitchat (no agent, static reply)
CLARIFY = "clarify"  # too ambiguous → ask the user one question
METRICS = "metrics"  # computed VALUE of a line → MetricsAgent (inline)

INTENTS = frozenset({QUESTION, CODE, PLOT, WORKFLOW, META, CLARIFY, METRICS})

# Intents that do NOT require an EDI dataset to be loaded.
# METRICS needs data (a line) — deliberately excluded.
NO_DATA_INTENTS = frozenset({QUESTION, CODE, META, CLARIFY})

# Confidence below which the LLM path downgrades to CLARIFY.
_CLARIFY_THRESHOLD = 0.45


# ── decision record ────────────────────────────────────────────────────────────


@dataclass
class RouterDecision:
    """Structured routing decision returned by :meth:`IntentRouter.route`.

    Attributes
    ----------
    intent : str
        One of :data:`INTENTS`.
    workflow : str or None
        Workflow slot, only meaningful when ``intent`` is
        :data:`WORKFLOW` or :data:`PLOT`.
    confidence : float
        ``0.0``–``1.0`` self-reported confidence.
    clarification : str or None
        A question to ask the user when ``intent == CLARIFY``.
    reasoning : str
        One-sentence rationale (useful for logs / debugging).
    source : str
        ``"llm"`` or ``"offline"`` — which path produced the decision.
    """

    intent: str
    workflow: str | None = None
    confidence: float = 0.0
    clarification: str | None = None
    reasoning: str = ""
    source: str = "offline"

    @property
    def needs_data(self) -> bool:
        """True when this intent requires a loaded EDI dataset."""
        return self.intent not in NO_DATA_INTENTS


# ── offline heuristic ──────────────────────────────────────────────────────────

# Greetings / capability probes → META.
_META_PHRASES = (
    "what can you do",
    "what can i ask",
    "who are you",
    "what are you",
    "your capabilities",
    "what do you do",
    "how do you work",
    "what is pycsamt agent",
    "help me get started",
    "getting started",
    "what can you help",
    "what you can do",
    "things you can do",
    "what you can perform",
    "introduce yourself",
    "introduce yourself.",
    "tell me about yourself",
    "about yourself",
    "who is this",
    "what is this assistant",
)

# Capability / catalogue listing → META ("list the agents/tasks/workflows you
# can run", "which workflows are available", "show me what you can do"). These
# must NOT fall through to WORKFLOW (which would blindly default to QC). The
# pattern deliberately needs a *listing* verb next to a *catalogue noun* so a
# genuine question like "what does the static-shift agent do" stays a QUESTION.
# Plural-only catalogue nouns: a question about ONE named thing ("what does
# the static-shift agent do") stays a QUESTION; only a plural list ("agents",
# "tasks", "workflows") trips the capability route.
_CAP_NOUN = (
    r"(agents|tasks|workflows|pipelines|capabilities|"
    r"commands|features|functions)\b"
)
_CAPABILITY_RE = re.compile(
    r"\b(list|show|give|get|enumerate|provide|display|see)\b.{0,30}"
    + _CAP_NOUN
    + r"|\b(available|supported|all the|all your)\b.{0,20}"
    + _CAP_NOUN
    + r"|\bwhich\b.{0,20}"
    + _CAP_NOUN
    + r"|\bwhat\b.{0,20}"
    + _CAP_NOUN
    + r".{0,20}\b(can|do|are|available|support|perform|run|have|there)\b"
    + r"|\b"
    + _CAP_NOUN
    + r" (you|i) can\b"
    + r"|\bwhat (can|do) you (do|perform|run|offer|help)\b",
    re.IGNORECASE,
)
_GREETINGS = frozenset(
    {
        "hi",
        "hello",
        "hey",
        "yo",
        "hiya",
        "thanks",
        "thank you",
        "good morning",
        "good afternoon",
        "good evening",
    }
)

# Explicit code-generation phrases → CODE.
_CODE_PHRASES = (
    "write code",
    "write python code",
    "write a python script",
    "generate code",
    "generate a script",
    "python script",
    "code for",
    "script for",
    "script that",
    "script to",
    "write a script",
    "write function",
    "write a function",
    "write class",
    "create notebook",
    "notebook for",
    "give me code",
    "show me code",
    "code example",
    "code snippet",
    "sample code",
    "example code",
    "snippet",
    "how to code",
    "produce code",
)

# Plot / visualise verbs → PLOT.
_PLOT_PHRASES = (
    "plot",
    "visualise",
    "visualize",
    "draw",
    "chart",
    "graph the",
    "show me a",
    "display the",
    "render the",
)

# Imperative verbs that signal "run this on my data" → WORKFLOW.
_WORKFLOW_VERBS = (
    "run",
    "perform",
    "execute",
    "compute",
    "invert",
    "inversion",
    "process",
    "correct",
    "denoise",
    "rotate",
    "decimate",
    "prepare",
    "analyse",
    "analyze",
    "evaluate",
    "load",
    "apply",
    "do ",
    "generate report",
    "static shift",
    "phase tensor",
    "quality control",
)

# Leading tokens / phrases that mark a genuine question → QUESTION.
_QUESTION_STARTS = (
    "what",
    "how",
    "why",
    "which",
    "who",
    "where",
    "when",
    "explain",
    "describe",
    "define",
    "tell me",
    "list the",
    "difference between",
    "can you explain",
    "what's",
    "whats",
    "is there",
    "are there",
    "does pycsamt",
    "do i need",
)

_PATH_RE = re.compile(r"[/~][\w/\\.\-]{3,}")
_LINE_RE = re.compile(r"\bL\d{2}PLT\b", re.I)

_AMBIGUOUS_COMMANDS = frozenset(
    {
        "run the workflow",
        "process my data",
        "fix the line",
        "make the plot",
        "do the inversion",
        "analyse this survey",
        "analyze this survey",
        "correct it",
        "prepare files",
        "run everything",
        "use the checkpoint",
    }
)


def classify_intent_offline(text: str) -> tuple[str, float]:
    """Classify *text* into an intent without any LLM call.

    Returns ``(intent, confidence)``.  Fast and deterministic; used both as
    the offline router backend and as the GUI pre-check that decides whether a
    message needs an EDI dataset (so questions/code/help never trip the data
    guard).

    Examples
    --------
    >>> classify_intent_offline("what does StaticShiftAgent do?")[0]
    'question'
    >>> classify_intent_offline("run AI inversion on /data/willy")[0]
    'workflow'
    >>> classify_intent_offline("write a script to load EDIs")[0]
    'code'
    >>> classify_intent_offline("hi")[0]
    'meta'
    """
    t = text.lower().strip()
    if not t:
        return META, 0.5

    # ── META: greetings & capability probes ────────────────────────────────
    if t in _GREETINGS or any(p in t for p in _META_PHRASES):
        return META, 0.9
    # capability / catalogue listing ("list the agents", "which workflows…")
    if _CAPABILITY_RE.search(t):
        return META, 0.88

    # ── CODE: explicit script requests ─────────────────────────────────────
    if any(p in t for p in _CODE_PHRASES):
        return CODE, 0.85

    # ── METRICS: a question about a computed VALUE of a line ────────────────
    # ("what's the strike of L22PLT?", "azimuth of all lines"). Checked before
    # QUESTION so value questions compute the number rather than fall through
    # to the package-concept Q&A path.
    from .metrics import looks_like_metric_query

    if looks_like_metric_query(text):
        return METRICS, 0.82

    has_path = bool(_PATH_RE.search(text))
    has_line = bool(_LINE_RE.search(text))
    has_wf_verb = any(v in t for v in _WORKFLOW_VERBS)

    if t in _AMBIGUOUS_COMMANDS:
        return CLARIFY, 0.72

    # ── QUESTION: genuine "what/how/why…" or trailing "?" ──────────────────
    is_question = t.endswith("?") or t.startswith(_QUESTION_STARTS)
    if is_question:
        # An imperative *with* a data path is a workflow phrased politely
        # ("can you run QC on /data/x?") — keep it a workflow.
        if has_wf_verb and (has_path or has_line):
            return WORKFLOW, 0.7
        return QUESTION, 0.8

    # ── PLOT: explicit visualisation verbs ─────────────────────────────────
    if any(p in t for p in _PLOT_PHRASES):
        return PLOT, 0.7

    # ── WORKFLOW: imperative verbs or a bare data path ─────────────────────
    if has_wf_verb or has_path:
        return WORKFLOW, 0.65

    # Nothing matched: most pycsamt chat is workflow-oriented, but with low
    # confidence so the LLM path (when available) can override.
    return WORKFLOW, 0.4


# ── LLM routing prompt ─────────────────────────────────────────────────────────

_KNOWN_WF = (
    "qc, static_shift, phase_analysis, forward, pre_inversion, "
    "inversion_eval, interpretation, report, full, ai_inversion, "
    "inv2d, inv3d, ensemble_inversion, joint_inversion, "
    "pinn_inversion, hybrid_inversion, modem, occam2d, tipper, "
    "sensitivity, rotation, freq_decimation, batch, comparison, "
    "code_gen, denoise, rhophi, phase_psection, pt_psection, tipper_plot, "
    "phase_tensor_map, pt_strip, pt_strip_grid, "
    "station_response, strike_profile, "
    "strike, dimensionality, validator, "
    "coords, elevation, converter, batch_export, "
    "freq_editor, layered_model"
)

# Append correction-method workflow ids from the central registry so the LLM
# router can name them directly (e.g. "corr_ss_ama").
try:  # pragma: no cover - registry is the single source of truth
    from ._corrections import (
        CORRECTION_METHODS as _CORR_METHODS,
    )

    if _CORR_METHODS:
        _KNOWN_WF = _KNOWN_WF + ", " + ", ".join(_CORR_METHODS)
except Exception:  # noqa: BLE001
    pass

_ROUTER_SYSTEM = f"""\
You are the top-level intent router for the pycsamt v2 magnetotelluric (MT)
assistant. Classify the user's message into exactly one INTENT.

INTENTS:
- "question": the user asks ABOUT pycsamt — concepts, what a class/function
  does, how something works, which method to use, definitions. They want an
  explanation, not execution.
- "code": the user wants a Python script / function / notebook generated.
- "plot": the user wants a figure/plot/visualisation produced from data.
- "workflow": the user wants to RUN a processing pipeline on their data
  (QC, static-shift, phase analysis, an inversion, report, etc.).
- "meta": greetings, "what can you do", "introduce yourself", capability
  questions about the assistant itself, and requests to LIST/ENUMERATE the
  agents, tasks or workflows the assistant can run ("list the agents", "which
  workflows are available", "what tasks can you perform"). A question about
  what ONE named agent/function does is a "question", not "meta".
- "metrics": the user asks for a COMPUTED VALUE of their survey line(s) and
  wants the number back inline — strike, azimuth/bearing, dimensionality,
  skew, station count, period/frequency range, coordinates/length, quality
  score, or a one-line summary ("what's the strike of L22PLT?", "azimuth of
  all lines", "how many stations", "tell me about this line"). This is NOT a
  plot/figure request and NOT "run an analysis".

Return ONLY a JSON object:
{{
  "intent": one of question|code|plot|workflow|meta|metrics,
  "workflow": one of [{_KNOWN_WF}] or null (only for workflow/plot/code),
  "confidence": float 0..1,
  "clarification": a single question to ask IF the request is too ambiguous
     to route, else null,
  "reasoning": one short sentence
}}

Rules:
- "How do I run an inversion?" is a QUESTION (asking for guidance).
- "Run an inversion on /data/x" is a WORKFLOW (asking to execute).
- "Write code to run an inversion" is CODE.
- Prefer "question" when the user clearly wants to learn, not execute.
- Set a low confidence (<0.5) and provide "clarification" when genuinely
  ambiguous.
"""


# ── agent ───────────────────────────────────────────────────────────────────────


class IntentRouter(BaseAgent):
    r"""Classify a chat message into a :class:`RouterDecision`.

    Parameters
    ----------
    api_key : str or None
        LLM key. ``None`` → pure offline heuristic (:func:`classify_intent_offline`).
    model, llm_provider : str
        Passed to :class:`~pycsamt.agents._base.BaseAgent`.

    Examples
    --------
    Offline::

        router = IntentRouter()
        d = router.route("what does StaticShiftAgent do?")
        d.intent          # 'question'
        d.needs_data      # False

    Online::

        router = IntentRouter(api_key="sk-ant-...", llm_provider="claude")
        d = router.route("run QC on /data/willy")
        d.intent          # 'workflow'
    """

    SYSTEM_PROMPT = _ROUTER_SYSTEM

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
    ) -> None:
        self._caller_key = api_key
        super().__init__(
            "IntentRouter",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
        )

    # The router is most ergonomic via .route(); execute() is provided so it
    # also satisfies the BaseAgent contract and can be chained if needed.
    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        text = (
            input_data.get("request")
            or input_data.get("text")
            or input_data.get("question")
            or ""
        )
        history = input_data.get("history")
        t0 = time.time()
        decision = self.route(text, history=history)
        return AgentResult(
            status="success",
            summary=(
                f"intent={decision.intent}"
                + (
                    f", workflow={decision.workflow}"
                    if decision.workflow
                    else ""
                )
            ),
            data={"decision": decision},
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=self._last_cost,
        )

    def route(
        self,
        text: str,
        *,
        history: list[dict] | None = None,
    ) -> RouterDecision:
        """Return a :class:`RouterDecision` for *text*.

        Uses the LLM when a key is configured, falling back to the offline
        heuristic on any failure or when offline.
        """
        text = (text or "").strip()
        if not text:
            return RouterDecision(
                intent=META,
                confidence=0.5,
                reasoning="empty message",
                source="offline",
            )

        # offline path
        if not self._caller_key:
            return self._offline_decision(text)

        # online path — structured JSON
        self._last_cost = 0.0
        user_msg = text
        if history:
            recent = history[-6:]
            convo = "\n".join(
                f"{m.get('role', 'user')}: {m.get('content', '')}"
                for m in recent
            )
            user_msg = (
                f"Recent conversation:\n{convo}\n\nCurrent message: {text}"
            )
        try:
            raw = self.query_llm(
                user_msg,
                system_message=self.SYSTEM_PROMPT,
                max_tokens=200,
                temperature=0.0,
            )
            parsed = self.extract_json(raw) if raw else None
        except Exception as exc:  # noqa: BLE001
            self._log.warning("Router LLM failed: %s", exc)
            parsed = None

        if not isinstance(parsed, dict):
            return self._offline_decision(text)

        intent = str(parsed.get("intent", "")).lower().strip()
        if intent not in INTENTS:
            return self._offline_decision(text)

        try:
            conf = float(parsed.get("confidence", 0.5))
        except (TypeError, ValueError):
            conf = 0.5
        wf = parsed.get("workflow") or None
        clar = parsed.get("clarification") or None

        # Low confidence with a usable clarification → ask instead of guess.
        if conf < _CLARIFY_THRESHOLD and clar:
            intent = CLARIFY

        return RouterDecision(
            intent=intent,
            workflow=str(wf) if wf else None,
            confidence=conf,
            clarification=clar,
            reasoning=str(parsed.get("reasoning", "")),
            source="llm",
        )

    # ── helpers ────────────────────────────────────────────────────────────
    @staticmethod
    def _offline_decision(text: str) -> RouterDecision:
        intent, conf = classify_intent_offline(text)
        return RouterDecision(
            intent=intent,
            confidence=conf,
            reasoning="offline heuristic",
            source="offline",
        )


__all__ = [
    "IntentRouter",
    "RouterDecision",
    "classify_intent_offline",
    "QUESTION",
    "CODE",
    "PLOT",
    "WORKFLOW",
    "META",
    "CLARIFY",
    "METRICS",
    "INTENTS",
    "NO_DATA_INTENTS",
]
