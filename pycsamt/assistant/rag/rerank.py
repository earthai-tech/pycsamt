# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.rerank
============================

Optional **LLM re-ranking** of retrieved chunks — a precision pass on top
of hybrid retrieval.

Retrieval (BM25 + expansion + optional dense) is tuned for *recall*: get
the right chunk into the top ~20. A cross-encoder / LLM re-ranker then
reads the query against each candidate and reorders for *precision*,
which is where lexical/dense retrievers drift (a word collision ranks an
unrelated helper above the real answer).

This is deliberately dependency-free and model-agnostic: the caller
passes a ``rank_fn(prompt, system) -> str`` (typically a thin wrapper over
:meth:`~pycsamt.agents._base.BaseAgent.query_llm`), so re-ranking is only
ever active when an API key is configured. Any failure or unparseable
reply falls back to the original order — re-ranking can improve results
but never break them.
"""

from __future__ import annotations

import re
from typing import Callable

from .schemas import RAGChunk

__all__ = [
    "build_rerank_prompt",
    "parse_ranking",
    "llm_rerank",
    "RANK_SYSTEM",
]

RANK_SYSTEM = (
    "You are a precise retrieval re-ranker for the pyCSAMT geophysics "
    "library. Given a user query and numbered candidate snippets, return "
    "the candidate numbers ordered from most to least relevant to the "
    "query. Reply with ONLY a comma-separated list of numbers (e.g. "
    "'3, 1, 5'). Omit candidates that are irrelevant."
)

_SNIPPET_CHARS = 240


def _label(chunk: RAGChunk) -> str:
    return chunk.symbol or chunk.title or chunk.source_path


def build_rerank_prompt(query: str, chunks: list[RAGChunk]) -> str:
    """Compose the candidate list shown to the re-ranker."""
    lines = [f"Query: {query}", "", "Candidates:"]
    for i, c in enumerate(chunks, start=1):
        snippet = " ".join((c.text or "").split())[:_SNIPPET_CHARS]
        lines.append(f"[{i}] ({c.kind}) {_label(c)} — {snippet}")
    lines.append("")
    lines.append("Return the candidate numbers, most relevant first, comma-separated.")
    return "\n".join(lines)


def parse_ranking(reply: str, n: int) -> list[int]:
    """Parse ``'3, 1, 5'`` into 0-based indices, de-duplicated and bounded.

    Numbers outside ``1..n`` are ignored. Returns an empty list when
    nothing parseable is found (caller then keeps the original order).
    """
    order: list[int] = []
    seen: set[int] = set()
    for tok in re.findall(r"\d+", reply or ""):
        v = int(tok) - 1
        if 0 <= v < n and v not in seen:
            seen.add(v)
            order.append(v)
    return order


def llm_rerank(
    query: str,
    chunks: list[RAGChunk],
    *,
    rank_fn: Callable[[str, str], str],
    top_n: int | None = None,
) -> list[RAGChunk]:
    """Reorder *chunks* by LLM-judged relevance to *query*.

    Parameters
    ----------
    rank_fn : callable ``(prompt, system) -> str``
        Runs the LLM; usually wraps ``agent.query_llm``. Only the first
        *top_n* candidates (default: all) are sent, to bound cost.
    top_n : int, optional
        Cap on candidates to re-rank.

    Any candidate the model omits is appended after the ranked ones in
    its original order, so nothing is silently dropped. Falls back to the
    input order on any error.
    """
    if not chunks:
        return chunks
    pool = chunks[:top_n] if top_n else chunks
    try:
        reply = rank_fn(build_rerank_prompt(query, pool), RANK_SYSTEM)
        order = parse_ranking(reply, len(pool))
    except Exception:  # noqa: BLE001 — rerank is best-effort
        return chunks
    if not order:
        return chunks
    ranked = [pool[i] for i in order]
    ranked_set = set(order)
    leftovers = [pool[i] for i in range(len(pool)) if i not in ranked_set]
    tail = chunks[len(pool) :] if top_n else []
    return ranked + leftovers + tail
