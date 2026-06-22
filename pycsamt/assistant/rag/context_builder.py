# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.context_builder
=====================================

Turn a user query into a compact, **cited** context block for the LLM —
the bridge between retrieval (:mod:`~pycsamt.assistant.rag.retriever`),
the project registry
(:class:`~pycsamt.assistant.tools.project_registry.ProjectRegistry`) and
the answering agent (:class:`~pycsamt.agents.package_qa.PackageQAAgent`).

It also composes a deterministic **offline answer** from the retrieved
chunks, so the assistant is useful even without an LLM key.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .retriever import Retriever, build_retriever
from .schemas import RAGChunk

__all__ = [
    "AssembledContext",
    "ContextBuilder",
    "default_context_builder",
]

_SNIPPET_CHARS = 700


def _snippet(chunk: RAGChunk, max_len: int = _SNIPPET_CHARS) -> str:
    """Most informative slice of a chunk (docstring/body over raw code)."""
    text = chunk.text or ""
    if "Docstring:" in text:
        after = text.split("Docstring:", 1)[1]
        body = after.split("\nCode:", 1)[0].strip()
        if body:
            text = body
    elif text.startswith("Document:"):
        # drop the "Document: <path>" first line
        text = text.split("\n", 1)[-1].strip()
    text = text.strip()
    return text if len(text) <= max_len else text[:max_len] + " …"


@dataclass
class AssembledContext:
    """Everything the answerer needs for one query."""

    query: str
    context_text: str = ""
    citations: list[dict[str, Any]] = field(default_factory=list)
    project_context: dict[str, Any] = field(default_factory=dict)
    symbols: list[dict[str, Any]] = field(default_factory=list)
    chunks: list[RAGChunk] = field(default_factory=list)

    def is_empty(self) -> bool:
        return not self.chunks and not self.project_context

    def compose_offline_answer(self, top: int = 3) -> str:
        """A readable, no-LLM answer assembled from the top chunks."""
        if self.is_empty():
            return (
                "I couldn't find anything relevant in the pyCSAMT "
                "reference for that. Try naming a class, function or "
                "workflow (e.g. static shift, phase tensor, AI inversion)."
            )
        parts: list[str] = []
        if self.project_context:
            pc = self.project_context
            parts.append(
                f"**Project line {pc.get('line')}** → "
                f"{pc.get('n_edi_files', 0)} EDI file(s) at "
                f"`{pc.get('edi_dir')}`"
                + ("" if pc.get("exists") else "  (path not found)")
            )
        parts.append("Based on the pyCSAMT reference:")
        for c in self.chunks[:top]:
            label = c.symbol or c.title or c.source_path
            parts.append(f"\n**{label}**\n{_snippet(c, 360)}")
        if self.citations:
            cites = "  ".join(
                f"[{c['n']}] {c['source_path']}"
                for c in self.citations[:top]
            )
            parts.append(f"\nSources: {cites}")
        return "\n".join(parts)


class ContextBuilder:
    """Assemble :class:`AssembledContext` from a retriever + registry."""

    def __init__(self, retriever: Retriever, registry: Any | None = None):
        self.retriever = retriever
        self.registry = registry

    def build(
        self,
        query: str,
        *,
        k: int = 6,
        max_chars: int = 6000,
    ) -> AssembledContext:
        ctx = self.retriever.search(query, k=k)

        project_context: dict[str, Any] = {}
        if self.registry is not None:
            try:
                line = self.registry.find_line_in_text(query)
                if line:
                    project_context = self.registry.resolve_line(line)
            except Exception:  # noqa: BLE001
                project_context = {}

        blocks: list[str] = []
        citations: list[dict[str, Any]] = []

        if project_context:
            pc = project_context
            blocks.append(
                "Project context:\n"
                f"  line: {pc.get('line')}\n"
                f"  edi_dir: {pc.get('edi_dir')} "
                f"(exists={pc.get('exists')}, "
                f"n_edi={pc.get('n_edi_files')})\n"
                f"  sort_by: {pc.get('sort_by')}, "
                f"output_root: {pc.get('output_root')}\n"
                f"  default_workflows: {pc.get('default_workflows')}\n"
                f"  static_shift_defaults: {pc.get('static_shift')}"
            )

        used = len("\n\n".join(blocks))
        for i, c in enumerate(ctx.chunks, start=1):
            label = c.symbol or c.title or c.source_path
            sig = c.metadata.get("signature")
            head = f"[{i}] {label}"
            if sig:
                head += f"  —  {sig}"
            block = (
                f"{head}\n"
                f"    (source: {c.source_path})\n"
                f"{_snippet(c)}"
            )
            if used + len(block) > max_chars and i > 1:
                break
            blocks.append(block)
            used += len(block)
            citations.append(
                {
                    "n": i,
                    "source_path": c.source_path,
                    "symbol": c.symbol,
                    "kind": c.kind,
                    "title": c.title,
                }
            )

        return AssembledContext(
            query=query,
            context_text="\n\n".join(blocks).strip(),
            citations=citations,
            project_context=project_context,
            symbols=ctx.symbols,
            chunks=ctx.chunks,
        )


# ── cached default builder (retriever + registry) ───────────────────────────────
_BUILDER_CACHE: dict[str, "ContextBuilder | None"] = {}


def default_context_builder(root=None) -> "ContextBuilder | None":
    """Build (and cache) a :class:`ContextBuilder` from the repo + registry.

    Returns ``None`` if the corpus can't be built (e.g. running from a
    wheel install with no source tree) — callers should degrade
    gracefully.
    """
    from .ingest import repo_root

    key = str(root or repo_root())
    if key in _BUILDER_CACHE:
        return _BUILDER_CACHE[key]

    builder: ContextBuilder | None
    try:
        retriever = build_retriever(root)
        if not retriever.chunks:
            builder = None
        else:
            registry = None
            try:
                from ..tools.project_registry import ProjectRegistry
                registry = ProjectRegistry.from_default(root)
            except Exception:  # noqa: BLE001
                registry = None
            builder = ContextBuilder(retriever, registry)
    except Exception:  # noqa: BLE001
        builder = None

    _BUILDER_CACHE[key] = builder
    return builder
