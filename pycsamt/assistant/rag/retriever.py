# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.retriever
===============================

Hybrid lexical retriever over the RAG corpus — pure Python, no heavy
dependencies (no faiss / sentence-transformers).

It combines:

* **BM25** ranking over tokenised chunk text (Okapi BM25), with an
  identifier-aware tokenizer so ``estimate_ss_ama`` /
  ``StaticShiftAgent`` split into searchable terms; and
* **structural boosts** — chunk priority, workflow match (inferred from
  the query or supplied), recipe chunks, and exact symbol mentions.

The interface is deliberately backend-agnostic: a future vector / rerank
stage can be layered behind :meth:`Retriever.search` without changing
callers. Results are returned as a
:class:`~pycsamt.assistant.rag.schemas.RetrievedContext`.
"""

from __future__ import annotations

import math
import re
from typing import Iterable

from .config import infer_workflow
from .schemas import RAGChunk, RetrievedContext

__all__ = ["tokenize", "BM25", "Retriever", "build_retriever"]

_STOPWORDS = frozenset(
    "a an the of to for and or in on at by with from is are be do does did "
    "how what which when where why can could would should i you we my our "
    "me it this that these those use using used run write get make show "
    "please code script function class method".split()
)

_WORD = re.compile(r"[A-Za-z0-9]+")
_CAMEL = re.compile(r"([a-z0-9])([A-Z])")


def _split_identifier(token: str) -> list[str]:
    """Split ``StaticShiftAgent`` / ``estimate_ss_ama`` into word parts."""
    spaced = _CAMEL.sub(r"\1 \2", token).replace("_", " ")
    return spaced.split()


def tokenize(text: str) -> list[str]:
    """Identifier-aware tokenizer: lowercase word + sub-word terms."""
    out: list[str] = []
    for raw in _WORD.findall(text or ""):
        low = raw.lower()
        parts = _split_identifier(raw)
        # keep the whole identifier (helps exact symbol matches) …
        if len(low) >= 2 and low not in _STOPWORDS:
            out.append(low)
        # … and its sub-words
        for p in parts:
            pl = p.lower()
            if pl != low and len(pl) >= 2 and pl not in _STOPWORDS:
                out.append(pl)
    return out


class BM25:
    """Okapi BM25 over pre-tokenised documents (inverted index)."""

    def __init__(
        self,
        corpus_tokens: list[list[str]],
        *,
        k1: float = 1.5,
        b: float = 0.75,
    ) -> None:
        self.k1 = k1
        self.b = b
        self.n_docs = len(corpus_tokens)
        self.doc_len = [len(d) for d in corpus_tokens]
        self.avgdl = (
            sum(self.doc_len) / self.n_docs if self.n_docs else 0.0
        )
        postings: dict[str, list[tuple[int, int]]] = {}
        df: dict[str, int] = {}
        for i, toks in enumerate(corpus_tokens):
            counts: dict[str, int] = {}
            for t in toks:
                counts[t] = counts.get(t, 0) + 1
            for t, f in counts.items():
                postings.setdefault(t, []).append((i, f))
                df[t] = df.get(t, 0) + 1
        self.postings = postings
        self.idf = {
            t: math.log(1.0 + (self.n_docs - n + 0.5) / (n + 0.5))
            for t, n in df.items()
        }

    def scores(self, query_tokens: Iterable[str]) -> list[float]:
        """BM25 score for every document (0.0 when no term matches)."""
        out = [0.0] * self.n_docs
        if not self.avgdl:
            return out
        for t in set(query_tokens):
            idf = self.idf.get(t)
            if idf is None:
                continue
            for i, f in self.postings[t]:
                dl = self.doc_len[i]
                denom = f + self.k1 * (
                    1.0 - self.b + self.b * dl / self.avgdl
                )
                out[i] += idf * (f * (self.k1 + 1.0)) / denom
        return out


class Retriever:
    """Rank :class:`RAGChunk` objects for a query (BM25 + boosts)."""

    # boost factors
    _PRIORITY_STEP = 0.12     # per priority level above 1
    _WORKFLOW_BOOST = 1.6
    _RECIPE_BOOST = 1.25
    _SYMBOL_BOOST = 1.5

    def __init__(self, chunks: list[RAGChunk]) -> None:
        self.chunks = chunks
        self._bm25 = BM25([tokenize(c.text) for c in chunks])

    def search(
        self,
        query: str,
        *,
        k: int = 8,
        kinds: Iterable[str] | None = None,
        workflow: str | None = None,
    ) -> RetrievedContext:
        """Return the top-*k* chunks for *query* as a RetrievedContext.

        Parameters
        ----------
        k : int
            Max chunks to return.
        kinds : iterable of str, optional
            Hard filter on chunk kind (e.g. ``{"recipe", "python_symbol"}``).
        workflow : str, optional
            Workflow boost hint; if omitted it is inferred from the query.
        """
        q_tokens = tokenize(query)
        base = self._bm25.scores(q_tokens)
        wf = workflow or infer_workflow(query)
        ql = query.lower()
        kinds_set = set(kinds) if kinds else None

        scored: list[tuple[float, int]] = []
        for i, chunk in enumerate(self.chunks):
            s = base[i]
            if s <= 0.0:
                continue
            if kinds_set and chunk.kind not in kinds_set:
                continue
            s *= 1.0 + self._PRIORITY_STEP * (chunk.priority - 1)
            if wf and chunk.workflow == wf:
                s *= self._WORKFLOW_BOOST
            if chunk.kind == "recipe":
                s *= self._RECIPE_BOOST
            if chunk.symbol:
                leaf = chunk.symbol.rsplit(".", 1)[-1].lower()
                if len(leaf) >= 3 and leaf in ql:
                    s *= self._SYMBOL_BOOST
            scored.append((s, i))

        scored.sort(key=lambda t: -t[0])
        top = [self.chunks[i] for _, i in scored[:k]]

        symbols = [
            {
                "symbol": c.symbol,
                "kind": c.kind,
                "workflow": c.workflow,
                "signature": c.metadata.get("signature"),
                "source_path": c.source_path,
            }
            for c in top
            if c.symbol
            and c.kind in ("python_symbol", "python_method", "module_doc")
        ]
        return RetrievedContext(
            query=query, chunks=top, symbols=symbols
        )


# ── convenience builder with a tiny cache ───────────────────────────────────────
_CACHE: dict[str, Retriever] = {}


def build_retriever(
    root=None,
    *,
    chunks: list[RAGChunk] | None = None,
    use_cache: bool = True,
) -> Retriever:
    """Build (or fetch a cached) :class:`Retriever`.

    With *chunks* given, builds directly from them (no disk scan).
    Otherwise ingests the repo via
    :func:`pycsamt.assistant.rag.ingest.build_chunks` and caches the
    result keyed by root.
    """
    if chunks is not None:
        return Retriever(chunks)

    from .ingest import build_chunks, repo_root

    key = str(root or repo_root())
    if use_cache and key in _CACHE:
        return _CACHE[key]
    retr = Retriever(build_chunks(root))
    if use_cache:
        _CACHE[key] = retr
    return retr
