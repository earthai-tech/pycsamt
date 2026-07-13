# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.feedback
==============================

Turn user **thumbs up / down** into a retrieval signal — the assistant
that learns from use.

When a user marks an answer helpful (or not), we record ``(query,
target, verdict)``. On a later, *similar* query, :meth:`FeedbackStore.
adjustments` returns a per-target score delta: chunks that were helpful
on near-identical questions are boosted, chunks that were rejected
become **hard negatives** and are pushed down. Similarity is token
Jaccard over the two queries, and every delta is scaled by it, so only
genuinely related history moves the ranking.

A *target* is a chunk's stable key — its fully-qualified symbol when it
has one, else its chunk id (see :func:`chunk_key`). Keying on the symbol
alone would silently ignore documentation and recipe chunks, which are
often the top hit.

The store is a plain JSONL file next to the index (``.pycsamt_rag/
feedback.jsonl``); it survives index rebuilds and is safe to delete. All
of this is deterministic and offline — no model, no network.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

__all__ = ["FeedbackStore", "chunk_key", "record_feedback"]


def chunk_key(chunk: Any) -> str:
    """Stable feedback key for a chunk: its symbol, else its chunk id.

    Code chunks are keyed by symbol so the verdict survives edits that
    shift line numbers; doc / recipe chunks fall back to their id.
    """
    return getattr(chunk, "symbol", None) or getattr(chunk, "id", "")


class FeedbackStore:
    """Append-only thumbs up/down log with a query-similarity adjuster."""

    FILENAME = "feedback.jsonl"
    _BOOST = 0.5  # max upward nudge for a helpful symbol
    _PENALTY = 0.6  # max downward nudge for a rejected symbol
    _MIN_OVERLAP = 0.34  # ignore history less similar than this (Jaccard)

    def __init__(
        self,
        path: Path | str | None = None,
        *,
        root: Path | str | None = None,
    ) -> None:
        if path is None:
            from .index_store import default_index_dir

            path = default_index_dir(root) / self.FILENAME
        self.path = Path(path)
        self._cache: list[dict] = []
        self._cache_mtime: float | None = None

    # ── write ───────────────────────────────────────────────────────────────
    def record(self, query: str, target: str, helpful: bool) -> None:
        """Append one verdict: *target* was/ wasn't helpful for *query*.

        *target* is a chunk key from :func:`chunk_key`.
        """
        if not (query and target):
            return
        rec = {
            "query": query,
            "target": target,
            "helpful": bool(helpful),
            "ts": time.time(),
        }
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

    # ── read (mtime-cached so the hot path doesn't re-read every query) ──────
    def load(self) -> list[dict]:
        try:
            mtime = self.path.stat().st_mtime
        except OSError:
            return []
        if self._cache_mtime == mtime:
            return self._cache
        out: list[dict] = []
        try:
            with self.path.open(encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line:
                        try:
                            out.append(json.loads(line))
                        except json.JSONDecodeError:
                            continue
        except OSError:
            return []
        self._cache, self._cache_mtime = out, mtime
        return out

    # ── the retrieval signal ────────────────────────────────────────────────
    def adjustments(
        self, query: str, *, min_overlap: float | None = None
    ) -> dict[str, float]:
        """Return ``{target: delta}`` learned from queries similar to *query*.

        ``delta > 0`` boosts a previously-helpful chunk; ``delta < 0``
        demotes a rejected one. Each contribution is weighted by the
        Jaccard similarity of the stored query to *query*, so unrelated
        history has no effect.
        """
        from .retriever import tokenize

        floor = self._MIN_OVERLAP if min_overlap is None else min_overlap
        q = set(tokenize(query))
        if not q:
            return {}
        deltas: dict[str, float] = {}
        for rec in self.load():
            target = rec.get("target")
            if not target:
                continue
            r = set(tokenize(rec.get("query", "")))
            if not r:
                continue
            jac = len(q & r) / len(q | r)
            if jac < floor:
                continue
            base = self._BOOST if rec.get("helpful") else -self._PENALTY
            deltas[target] = deltas.get(target, 0.0) + base * jac
        return deltas


def record_feedback(
    query: str, target: str, helpful: bool, *, root: Path | str | None = None
) -> None:
    """Convenience: record one verdict into the default feedback store.

    *target* is a chunk key — pass ``chunk_key(chunk)``.
    """
    FeedbackStore(root=root).record(query, target, helpful)
