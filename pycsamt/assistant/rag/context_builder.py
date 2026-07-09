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

import re
from dataclasses import dataclass, field
from typing import Any

from .retriever import Retriever, build_retriever
from .schemas import RAGChunk, RetrievedContext

__all__ = [
    "AssembledContext",
    "ContextBuilder",
    "default_context_builder",
    "needs_clarification",
]

_SNIPPET_CHARS = 700


# NumPy-doc section headers: everything from the first one on is API
# reference detail, not answer prose.
_NUMPYDOC_SECTION_RE = re.compile(
    r"\n\s*(Parameters|Returns|Yields|Raises|See Also|Notes|Examples"
    r"|Attributes|References|Warnings)\s*\n\s*-{3,}",
)


def _snippet(chunk: RAGChunk, max_len: int = _SNIPPET_CHARS) -> str:
    """Most informative prose slice of a chunk.

    For code chunks: the docstring prose (cut before the first numpydoc
    section). For doc/recipe chunks: the body with the (separately-shown)
    heading and embedded code fences stripped, so the offline answer
    reads cleanly instead of dumping raw markdown.
    """
    text = chunk.text or ""
    if "Docstring:" in text:
        body = text.split("Docstring:", 1)[1].split("\nCode:", 1)[0]
        body = _NUMPYDOC_SECTION_RE.split(body, 1)[0]
        if body.strip():
            text = body
    elif text.startswith("Document:"):
        text = text.split("\n", 1)[-1]  # drop "Document: <path>"
        # drop a leading heading line equal to the chunk title (it is
        # already shown as the block label) and any fenced code blocks
        title = (chunk.title or "").strip()
        lines = text.splitlines()
        if lines and lines[0].strip() == title:
            lines = lines[1:]
        text = "\n".join(lines)
        text = re.sub(r"```.*?```", "", text, flags=re.S)
    # Drop indexer metadata lines ("Module: …", "File: …", "Class: …",
    # "Signature: …") — they are chunk headers, not prose, and make the
    # offline answer read like a raw index dump.
    text = "\n".join(
        ln for ln in text.splitlines()
        if not re.match(
            r"\s*(Module|File|Class|Function|Method|Signature)\s*:\s",
            ln,
        )
    )
    text = re.sub(r"\n{3,}", "\n\n", text).strip()
    return text if len(text) <= max_len else text[:max_len].rstrip() + " …"


def _first_sentence(text: str, max_len: int = 180) -> tuple[str, str]:
    """Split *text* into (first sentence, remainder) for lead-in style."""
    text = (text or "").strip()
    if not text:
        return "", ""
    first_para = text.split("\n\n", 1)[0].replace("\n", " ").strip()
    m = re.search(r"(?<=[.!?])\s", first_para)
    first = first_para[: m.start()] if m else first_para
    if len(first) > max_len:
        first = first[: max_len - 1].rstrip() + "…"
    rest = text[len(first):].lstrip(" .\n") if text.startswith(first) else (
        text.split("\n\n", 1)[1].strip() if "\n\n" in text else ""
    )
    return first, rest


def _chunk_label(chunk: RAGChunk) -> str:
    """Readable label: the symbol path, or ``docname: Section`` for
    doc/recipe chunks (a bare section heading lacks context)."""
    if chunk.symbol:
        return chunk.symbol
    title = (chunk.title or "").strip()
    src = (chunk.source_path or "").replace("\\", "/")
    stem = src.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    if title and stem and title.lower() != stem.lower():
        return f"{stem}: {title}"
    return title or src


# Boosted-BM25 score below which retrieval is treated as weak. Calibrated
# on the pyCSAMT corpus: on-topic / paraphrased queries score ~30–80,
# contentless follow-ups ("the thing", "do that again") score <20.
_CLARIFY_SCORE_FLOOR = 25.0


@dataclass
class AssembledContext:
    """Everything the answerer needs for one query."""

    query: str
    context_text: str = ""
    citations: list[dict[str, Any]] = field(default_factory=list)
    project_context: dict[str, Any] = field(default_factory=dict)
    symbols: list[dict[str, Any]] = field(default_factory=list)
    chunks: list[RAGChunk] = field(default_factory=list)
    top_score: float = 0.0

    def is_empty(self) -> bool:
        return not self.chunks and not self.project_context

    def is_confident(self, floor: float = _CLARIFY_SCORE_FLOOR) -> bool:
        """True when a project line resolved or retrieval scored well."""
        return bool(self.project_context) or self.top_score >= floor

    def _lead_chunk(self) -> RAGChunk | None:
        """Pick the best 'definition' chunk to open the answer.

        A code symbol leads only when the query actually names it
        (e.g. "StaticShiftAgent"); otherwise the retriever's top-ranked
        chunk leads — for how-to questions that is typically a recipe
        or documentation section, not an internal helper.
        """
        ql = self.query.lower()
        code_kinds = ("python_symbol", "python_method", "module_doc")
        for c in self.chunks:
            if c.symbol and c.kind in code_kinds:
                leaf = c.symbol.rsplit(".", 1)[-1].lower()
                if len(leaf) >= 3 and leaf in ql:
                    return c
        return self.chunks[0] if self.chunks else None

    def compose_offline_answer(self, top: int = 3) -> str:
        """A readable, no-LLM answer assembled from the top chunks.

        Leads with the most relevant symbol's docstring, then adds a
        couple of supporting snippets and the key/related symbols. (With
        an API key the LLM synthesises a fluent answer from the same
        retrieved context instead.)
        """
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
                f"**Survey line {pc.get('line')}** — "
                f"{pc.get('n_edi_files', 0)} EDI file(s) at "
                f"`{pc.get('edi_dir')}`"
                + ("" if pc.get("exists") else "  (path not found)")
            )

        seen: set[str] = set()

        # Definition-style opening: the best symbol, dash, its first
        # docstring sentence — then the rest of that docstring as a
        # normal paragraph.
        lead = self._lead_chunk()
        if lead is not None:
            label = _chunk_label(lead)
            seen.add(label)
            body = _snippet(lead, 460)
            if body.lstrip().startswith(("-", "*", "1.")):
                # A list-shaped chunk (e.g. a recipe's step list) reads
                # best kept as markdown, not flattened to one line.
                opening = f"**`{label}`**\n\n{body}"
            else:
                first, rest = _first_sentence(body)
                opening = f"**`{label}`**"
                if first:
                    opening += f" — {first}"
                if rest:
                    opening += f"\n\n{rest}"
            parts.append(opening)

        # Supporting chunks become one-line bullets instead of
        # further docstring dumps.
        bullets: list[str] = []
        for c in self.chunks:
            if len(bullets) >= max(0, top - 1):
                break
            label = _chunk_label(c)
            if label in seen:
                continue
            first, _ = _first_sentence(_snippet(c, 220))
            if not first:
                continue
            # A lead-in to a list ("Emit a tidy header section:")
            # squashes badly on one line — keep only the lead-in.
            if ":" in first:
                first = first.split(":", 1)[0].rstrip() + "."
            seen.add(label)
            bullets.append(f"- **`{label}`** — {first}")
        if bullets:
            parts.append("Also relevant:\n" + "\n".join(bullets))

        # Related API symbols (dedup, skip shown + test fixtures).
        rel: list[str] = []
        for s in self.symbols:
            sym = s.get("symbol")
            if not sym or sym in rel or sym in seen:
                continue
            if ".tests." in sym or ".test_" in sym:
                continue
            rel.append(sym)
        if rel:
            parts.append(
                "See also: "
                + ", ".join(f"`{s}`" for s in rel[:4])
            )

        if self.citations:
            paths: list[str] = []
            for c in self.citations:
                p = c["source_path"]
                if p not in paths:
                    paths.append(p)
            cites = " · ".join(f"`{p}`" for p in paths[:top])
            parts.append(f"*Sources: {cites}*")
        return "\n\n".join(parts)


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
        rerank_fn: Any | None = None,
        rerank_top_n: int = 12,
        session: dict[str, Any] | None = None,
    ) -> AssembledContext:
        """Assemble context for *query*.

        When *rerank_fn* (a ``(prompt, system) -> str`` LLM callable) is
        given, retrieval fetches a deeper candidate pool (*rerank_top_n*),
        an LLM re-ranks it for precision, and the top *k* are kept. Without
        it, the top *k* from hybrid retrieval are used directly.

        When *session* (``{last_workflow, last_line, recent_turns}``) is
        given, a subject-less follow-up query is rewritten for *retrieval*
        so it inherits the conversation's topic — the original *query* is
        still used for project-line resolution and is echoed back.
        """
        retrieval_query = query
        if session:
            from .rewrite import rewrite_query

            retrieval_query = rewrite_query(
                query,
                last_workflow=session.get("last_workflow"),
                last_line=session.get("last_line"),
                recent_turns=session.get("recent_turns"),
            )

        if rerank_fn is not None:
            pool = self.retriever.search(retrieval_query, k=max(k, rerank_top_n))
            from .rerank import llm_rerank

            reranked = llm_rerank(
                retrieval_query, pool.chunks, rank_fn=rerank_fn,
                top_n=rerank_top_n,
            )
            ctx = RetrievedContext(
                query=retrieval_query,
                chunks=reranked[:k],
                symbols=pool.symbols,
                top_score=pool.top_score,
            )
        else:
            ctx = self.retriever.search(retrieval_query, k=k)

        project_context: dict[str, Any] = {}
        if self.registry is not None:
            try:
                line = self.registry.find_line_in_text(query)
                if not line and session and session.get("last_line"):
                    # A follow-up like "now on that line" inherits the line.
                    from .rewrite import is_follow_up

                    if is_follow_up(query):
                        line = session["last_line"]
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
            top_score=ctx.top_score,
        )


# ── retrieval-confidence → clarification ────────────────────────────────────────

_CLARIFY_QUESTION = (
    "I want to point this at the right thing before I answer. Could you say "
    "which **dataset or survey line** (e.g. L22PLT) and which **workflow** "
    "you mean — QC, static shift, phase-tensor analysis, denoising, "
    "sensitivity, or an inversion?"
)


def needs_clarification(
    assembled: AssembledContext,
    *,
    floor: float = _CLARIFY_SCORE_FLOOR,
) -> str | None:
    """Return a clarifying question when answering would be a guess.

    Fires only for an **anaphoric follow-up** ("do that again", "same for
    this one") that stayed weak — i.e. retrieval is not confident
    (:meth:`AssembledContext.is_confident`) *and* the query refers back to
    something the assistant can't identify (no session resolved it). A
    self-contained question, even a rare low-scoring one ("the Sites
    class"), is answered best-effort rather than bounced.
    """
    if assembled.is_confident(floor):
        return None
    from .rewrite import is_follow_up

    if not is_follow_up(assembled.query):
        return None
    return _CLARIFY_QUESTION


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
