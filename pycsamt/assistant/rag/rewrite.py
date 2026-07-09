# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.rewrite
=============================

Deterministic **conversational query rewriting** for retrieval.

Retrieval is single-turn by nature: it sees only the words of the current
message. But chat is multi-turn — *"now do it on line 3"*, *"plot that
again"*, *"why did that fail"* carry their real subject in the **previous**
turn. Without folding that context in, the retriever matches on
``it`` / ``that`` / ``again`` (i.e. nothing) and drifts.

This module rewrites the *retrieval* query only (never the user-facing
text) by inheriting the session's last workflow / active line when — and
only when — the current message is a subject-less follow-up. It is
deliberately deterministic (no LLM): an explicit subject in the message
always wins, so a fresh, self-contained question is never altered.
"""

from __future__ import annotations

import re

from .config import WORKFLOW_KEYWORDS, infer_workflow

__all__ = ["is_follow_up", "rewrite_query"]

# Cues that a message continues a previous one rather than introducing a
# new subject ("run it again", "same for line 3", "now do that").
_FOLLOWUP_CUES = (
    "again", "same", "as well", "as before", "like before", "redo",
    "repeat", "once more", "also do", "then do", "do that", "do this",
    "do it", "run it", "that one", "this one", "the other one",
    "same thing", "same as", "now do", "and now",
)
# Bare anaphora words that, standing in for the subject, signal a follow-up.
_ANAPHORA = frozenset(
    "it that this them those these one".split()
)
# Generic references to "the line" without naming one.
_LINE_REF = re.compile(
    r"\b(that|this|the|same|it|its)\s+line\b|\bsame\s+data\b", re.IGNORECASE
)

_WORD = re.compile(r"[a-z0-9]+")


def is_follow_up(query: str) -> bool:
    """True when *query* reads as a continuation lacking its own subject.

    Requires a follow-up cue (or a bare anaphor) *and* the absence of an
    explicit workflow keyword — an on-topic message that happens to say
    "again" (e.g. "run static shift again") names its own subject and is
    not treated as context-dependent.
    """
    ql = (query or "").lower()
    if infer_workflow(ql):
        return False
    if any(cue in ql for cue in _FOLLOWUP_CUES):
        return True
    words = _WORD.findall(ql)
    # A short message that is essentially just an anaphor ("that", "it too").
    content = [w for w in words if w not in _ANAPHORA]
    return bool(words) and len(content) <= 1 and any(
        w in _ANAPHORA for w in words
    )


def _workflow_terms(workflow: str) -> str:
    """A couple of human-readable keywords for *workflow* (for the query)."""
    kws = WORKFLOW_KEYWORDS.get(workflow, ())
    # Prefer multi-word natural phrasings over bare identifiers.
    natural = [k for k in kws if " " in k] or list(kws)
    return " ".join(natural[:2])


def rewrite_query(
    query: str,
    *,
    last_workflow: str | None = None,
    last_line: str | None = None,
    recent_turns: list[dict] | None = None,
) -> str:
    """Return a retrieval query enriched with inherited session context.

    When *query* is a subject-less follow-up, append the *last_workflow*'s
    keywords (so the retriever recovers the topic) and the *last_line*
    token if the message refers to "the line" without naming it. A
    self-contained query is returned unchanged.
    """
    q = (query or "").strip()
    if not is_follow_up(q):
        return q

    # Fall back to the most recent workflow named in the transcript.
    wf = last_workflow
    if wf is None and recent_turns:
        for turn in reversed(recent_turns):
            wf = infer_workflow(turn.get("content", ""))
            if wf:
                break

    extra: list[str] = []
    if wf:
        terms = _workflow_terms(wf)
        if terms:
            extra.append(terms)
    if last_line and _LINE_REF.search(q) and last_line.lower() not in q.lower():
        extra.append(last_line)

    return f"{q} {' '.join(extra)}".strip() if extra else q
