# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents.package_qa
==========================
:class:`PackageQAAgent` — answer free-form questions
about the pycsamt v2 package.

This agent solves the "LLM does not know pycsamt"
problem by:

1. Extracting the full pycsamt v2 API reference
   (docstrings + signatures + parameter tables)
   from the live package at import time.

2. Using *query-adaptive injection*: instead of
   always dumping the entire reference into the
   prompt, it selects the most relevant tiers
   based on keywords in the user's question.

   Example:
     "How do I access impedance Z values?"
     -> injects TIER_SITES + TIER_EXAMPLES only
        (not the 15-agent class dump)

   "What does StaticShiftAgent do?"
     -> injects TIER_AGENTS + TIER_EXAMPLES

   "What workflows are supported?"
     -> injects TIER_CORE only

   This keeps the LLM's signal-to-noise ratio
   high — the relevant information is at the top,
   not buried in unrelated content.

Offline mode (no api_key)
--------------------------
Uses keyword-matching over extracted docstrings
to return a best-effort answer. No LLM call, no
cost, CI-safe.

Online mode (api_key provided)
--------------------------------
Sends question + selected context tiers to the
LLM provider and returns the grounded answer.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from ._base import AgentResult, BaseAgent
from ._package_context import (
    PACKAGE_CONTEXT,
    TIER_AGENTS,
    TIER_CORE,
    TIER_EMTOOLS,
    TIER_EXAMPLES,
    TIER_SITES,
    _collect_agents,
    _collect_core,
    _first_para as _short_doc,
)

if TYPE_CHECKING:
    pass


# ── system prompt template ─────────────────────────────────────────────────────

_SYSTEM_TMPL = (
    "You are a helpful expert on the pycsamt v2"
    " Python library for magnetotelluric data"
    " processing.\n\n"
    "The pycsamt v2 API reference extracted from"
    " the live package is provided below."
    " Use it as your *primary source of truth*."
    " If an answer is not in the reference, say"
    " so — do not invent class names, parameters,"
    " or behaviours that are not listed.\n\n"
    "---\n{context}\n---\n\n"
    "Answer in clear technical English."
    " Include short code examples when helpful."
)

# Full-dump system prompt (used as class attribute
# so tests can inspect it without instantiating)
_FULL_SYSTEM_PROMPT = _SYSTEM_TMPL.format(
    context=PACKAGE_CONTEXT
)


# ── query → tier selection ─────────────────────────────────────────────────────

# Keywords that indicate the question is about
# data access, Sites, EDI, or the data model.
_SITES_KW = frozenset(
    "sites site edi impedance tensor z freq"
    " frequency resistivity phase rho"
    " station load loader ensure_sites"
    " write_sites edicollection edif edifile"
    " write export access data attribute"
    " from_any select edit closest profile".split()
)

# Keywords that indicate the question is about
# a specific agent class.
_AGENT_KW = frozenset(
    "agent staticshift dataqa dataqc qc"
    " denois denoise phaseanaly phase"
    " forward inversion invert pinn"
    " hybrid neural ai inv inv2d inv3d"
    " report loader orchestrat interpret"
    " codegen code_gen correct correction"
    " shift galvanic method parameter"
    " parameter checkpoint epoch".split()
)

# Keywords for workflow / routing questions.
_WORKFLOW_KW = frozenset(
    "workflow workflows supported support"
    " pipeline route routing list available"
    " types kind kinds step steps".split()
)

# Keywords that clearly need code examples.
_EXAMPLE_KW = frozenset(
    "example how show code snippet use"
    " usage run execute load apply"
    " create instantiate import".split()
)

# Keywords indicating emtools.ss / static-shift
# function-level questions (not just the agent).
_EMTOOLS_KW = frozenset(
    "ama estimate_ss correct_ss read_edis"
    " emtools plot_ss fac_z fac_rho"
    " delta_log10 half_window sounding"
    " pseudosection logrho survey"
    " collection shift_factor"
    " correct_ss_ama estimate_ss_ama"
    " plot_ss_summary plot_ss_1d"
    " static galvanic shift_correction"
    " write_code write_script script"
    " logRho logRho_before".split()
)

_STOPWORDS = frozenset(
    "a an the is are was were be been"
    " do does did can could would how"
    " what which who where when why"
    " i you we they it in of to for"
    " on at by with about from".split()
)


def _keywords(text: str) -> list[str]:
    """Extract meaningful lowercase tokens."""
    words = re.findall(
        r"[a-z][a-z_0-9]{2,}", text.lower()
    )
    return [w for w in words if w not in _STOPWORDS]


def _select_tiers(question: str) -> str:
    r"""
    Choose which context tiers to inject.

    Returns a string containing 1-5 tiers,
    selected by keyword overlap with the question.
    Always includes TIER_CORE (it's small and
    provides the essential overview).
    """
    kw = set(_keywords(question))

    need_sites    = bool(kw & _SITES_KW)
    need_agents   = bool(kw & _AGENT_KW)
    need_workflow = bool(kw & _WORKFLOW_KW)
    need_examples = bool(kw & _EXAMPLE_KW)
    need_emtools  = bool(kw & _EMTOOLS_KW)

    # If nothing matched, include everything
    # (better to over-inform than under-inform)
    nothing_matched = not any(
        [need_sites, need_agents,
         need_workflow, need_examples,
         need_emtools]
    )

    parts = [TIER_CORE]
    if need_agents or nothing_matched:
        parts.append(TIER_AGENTS)
    if need_sites or nothing_matched:
        parts.append(TIER_SITES)
    if need_examples or nothing_matched:
        parts.append(TIER_EXAMPLES)
    if need_emtools or nothing_matched:
        parts.append(TIER_EMTOOLS)

    return "\n".join(parts)


# ── offline lookup ─────────────────────────────────────────────────────────────

_WORKFLOW_PREFIXES = (
    "qc", "static_shift",
    "phase_analysis", "forward",
    "inversion", "ai_inv",
    "pinn", "hybrid", "full",
    "report", "interp", "modem",
    "occam", "tipper",
)

_WORKFLOW_TRIGGERS = frozenset(
    "workflow workflows support supports"
    " available list".split()
)


def _workflow_lines() -> list[str]:
    return [
        ln
        for ln in PACKAGE_CONTEXT.split("\n")
        if ln.strip() and any(
            ln.strip().startswith(k)
            for k in _WORKFLOW_PREFIXES
        )
    ]


def _offline_answer(question: str) -> dict:
    r"""
    Answer a question without an LLM.

    Strategy:
    1. Workflow-list questions -> return workflow table
    2. Other questions -> score agents/classes by
       keyword overlap, return top-3 docstrings
    """
    kw = set(_keywords(question))

    # workflow-listing shortcut
    if kw & _WORKFLOW_TRIGGERS:
        lines = _workflow_lines()
        answer = (
            "pycsamt v2 supports the following"
            " workflows:\n\n"
            + "\n".join("  " + ln for ln in lines)
        )
        return {
            "answer": answer,
            "source": "docstring_lookup",
            "excerpts": [],
        }

    # class-level scoring
    agents = _collect_agents()
    core   = _collect_core()
    all_cls = {**agents, **core}

    scored: list[tuple[int, str, str]] = []
    for name, cls in all_cls.items():
        doc  = _short_doc(
            inspect_full_doc(cls), max_chars=600
        )
        hits = kw & set(_keywords(name + " " + doc))
        scored.append((len(hits), name, doc))

    scored.sort(key=lambda t: -t[0])
    top = [t for t in scored if t[0] > 0][:3]

    if not top:
        lines = _workflow_lines()
        answer = (
            "I could not find a specific match"
            " in the pycsamt API reference."
            " Supported workflows:\n\n"
            + "\n".join("  " + ln for ln in lines)
        )
        return {
            "answer": answer,
            "source": "docstring_lookup",
            "excerpts": [],
        }

    excerpts = [
        {"class": t[1], "excerpt": t[2]}
        for t in top
    ]
    answer = (
        "Based on the pycsamt v2 API reference:\n\n"
        + "\n\n".join(
            f"**{t[1]}**: {t[2]}" for t in top
        )
    )
    return {
        "answer": answer,
        "source": "docstring_lookup",
        "excerpts": excerpts,
    }


def inspect_full_doc(obj: Any) -> str:
    """Re-export: return full docstring."""
    import inspect as _i
    return _i.getdoc(obj) or ""


# keep Any importable at runtime
from typing import Any  # noqa: E402


# ── agent class ────────────────────────────────────────────────────────────────

class PackageQAAgent(BaseAgent):
    r"""Answer free-form questions about pycsamt v2.

    Parameters
    ----------
    api_key : str or None
        LLM key.  ``None`` activates offline
        docstring-lookup mode (no cost, no network).
    model, llm_provider : str
        Passed to BaseAgent.

    Notes
    -----
    Uses *query-adaptive context injection*: selects
    the most relevant documentation tiers (agents,
    Sites data model, usage examples) based on
    keywords in the question, rather than always
    dumping the full reference. This keeps the LLM
    focused on what's relevant.

    Examples
    --------
    Offline (no LLM)::

        agent = PackageQAAgent()
        r = agent.execute({
            "question": (
                "What does StaticShiftAgent do?"
            )
        })
        print(r["answer"])

    Online (Claude)::

        agent = PackageQAAgent(
            api_key="sk-ant-...",
            llm_provider="claude",
        )
        r = agent.execute({
            "question": (
                "How do I access impedance Z?"
            )
        })
        print(r["answer"])
    """

    # Full-dump system prompt stored as class attr
    # so tests can inspect it without instantiating.
    SYSTEM_PROMPT = _FULL_SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        use_rag: bool = True,
    ) -> None:
        # capture before AGENT_CONFIG resolves key
        self._caller_key = api_key
        self.use_rag = use_rag
        super().__init__(
            "PackageQAAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
        )

    def _build_rag(self, question: str):
        """Return an AssembledContext for *question*, or ``None``.

        Lazily uses the RAG layer (:mod:`pycsamt.assistant.rag`); any
        failure (assistant not installed, no source tree) degrades
        silently so QA still works without retrieval.
        """
        if not self.use_rag:
            return None
        try:
            from pycsamt.assistant.rag.context_builder import (
                default_context_builder,
            )
            builder = default_context_builder()
            if builder is None:
                return None
            ctx = builder.build(question)
            return None if ctx.is_empty() else ctx
        except Exception:  # noqa: BLE001 — RAG is best-effort
            return None

    def execute(
        self, input_data: dict
    ) -> AgentResult:
        question = (
            input_data.get("question")
            or input_data.get("request")
            or ""
        ).strip()

        if not question:
            return AgentResult(
                status="failed",
                summary="No question provided.",
                error=(
                    "No 'question' in input_data."
                ),
                data={},
            )

        extra_ctx = input_data.get("context", "")
        rag = self._build_rag(question)
        citations = rag.citations if rag else []

        # offline path
        if self._caller_key is None:
            if rag is not None:
                # RAG-composed answer beats the docstring keyword lookup
                answer = rag.compose_offline_answer()
                return AgentResult(
                    status="success",
                    summary=answer[:120],
                    data={
                        "answer": answer,
                        "source": "rag_offline",
                        "citations": citations,
                        "excerpts": [],
                    },
                )
            data = _offline_answer(question)
            return AgentResult(
                status="success",
                summary=data["answer"][:120],
                data=data,
            )

        # online path — select relevant tiers
        selected_ctx = _select_tiers(question)
        system_prompt = _SYSTEM_TMPL.format(
            context=selected_ctx
        )

        # Ground the LLM in retrieved, citable package facts.
        msg_parts: list[str] = []
        if rag is not None and rag.context_text:
            msg_parts.append(
                "Retrieved pyCSAMT context (prefer these real symbols; "
                "cite sources as [n] when used):\n"
                f"{rag.context_text}"
            )
        if extra_ctx:
            msg_parts.append(f"Session context: {extra_ctx}")
        msg_parts.append(f"Question: {question}")
        user_msg = "\n\n".join(msg_parts)

        try:
            llm_answer = self.query_llm(
                user_msg,
                max_tokens=1024,
                system_message=system_prompt,
            )
        except Exception as exc:
            return AgentResult(
                status="failed",
                summary="LLM call failed.",
                error=str(exc),
                data={},
            )

        data = {
            "answer": llm_answer or "(no answer)",
            "source": "llm+rag" if rag is not None else "llm",
            "excerpts": [],
            "citations": citations,
            "tiers_used": _tiers_used(question),
        }
        return AgentResult(
            status="success",
            summary=(llm_answer or "")[:120],
            data=data,
        )


def _tiers_used(question: str) -> list[str]:
    """Return names of tiers selected for question."""
    kw = set(_keywords(question))
    out = ["core"]
    if kw & _AGENT_KW:
        out.append("agents")
    if kw & _SITES_KW:
        out.append("sites")
    if kw & _EXAMPLE_KW:
        out.append("examples")
    if kw & _EMTOOLS_KW:
        out.append("emtools")
    if len(out) == 1:
        out = [
            "core", "agents", "sites",
            "examples", "emtools",
        ]
    return out
