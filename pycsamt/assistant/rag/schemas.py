# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.schemas
=============================

Core data types for the RAG layer.

A :class:`RAGChunk` is one retrievable unit — a Python symbol
(function/class/method), a documentation section, a recipe, etc. —
carrying enough metadata to rank, filter and cite it. :class:`RetrievedContext`
is what the retriever hands to the context/prompt builder.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any

__all__ = ["RAGChunk", "RetrievedContext"]


@dataclass
class RAGChunk:
    """A single retrievable unit of knowledge.

    Attributes
    ----------
    id : str
        Stable identifier, e.g. ``"pycsamt/emtools/ss.py:120:180"``.
    text : str
        The emb<eddable / searchable text (docstring + code, or doc prose).
    source_path : str
        Repo-relative path the chunk came from.
    kind : str
        ``"python_symbol"`` | ``"python_method"`` | ``"module_doc"`` |
        ``"doc_section"`` | ``"recipe"`` …
    symbol : str or None
        Fully-qualified symbol, e.g. ``"pycsamt.emtools.ss.estimate_ss_ama"``.
    workflow : str or None
        Associated workflow id (``"static_shift"``, ``"qc"``, …) when known.
    module : str or None
        Dotted module path for code chunks.
    title : str or None
        Heading/name for doc or recipe chunks.
    priority : int
        Static importance (higher = boosted at retrieval time).
    metadata : dict
        Free-form extras (signature, line range, parameters, …).
    """

    id: str
    text: str
    source_path: str
    kind: str
    symbol: str | None = None
    workflow: str | None = None
    module: str | None = None
    title: str | None = None
    priority: int = 1
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """JSON-serialisable representation (for chunks.jsonl)."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "RAGChunk":
        """Rebuild a chunk from :meth:`to_dict` output."""
        known = {f for f in cls.__dataclass_fields__}  # noqa: PLC0208
        return cls(**{k: v for k, v in d.items() if k in known})


@dataclass
class RetrievedContext:
    """Result of a retrieval call, ready for context assembly.

    Attributes
    ----------
    query : str
        The (possibly rewritten) query that produced these results.
    chunks : list of RAGChunk
        Ranked chunks, most relevant first.
    symbols : list of dict
        Symbol-graph hits (API entries) when available.
    project_context : dict
        Resolved project/line context (from the project registry).
    """

    query: str
    chunks: list[RAGChunk] = field(default_factory=list)
    symbols: list[dict[str, Any]] = field(default_factory=list)
    project_context: dict[str, Any] = field(default_factory=dict)
