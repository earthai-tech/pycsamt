# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.doc_indexer
=================================

Turn documentation (``.rst`` / ``.md``) into section-level
:class:`~pycsamt.assistant.rag.schemas.RAGChunk` objects.

Sectioning (rather than whole-file chunks) keeps retrieval focused: a
question about static-shift should pull the static-shift *section*, not
an entire tutorial. A lightweight splitter handles both:

* reStructuredText — a title line followed by a run of underline
  punctuation (``===``, ``---``, ``~~~`` …);
* Markdown — ATX headings (``#``, ``##`` …).
"""

from __future__ import annotations

import re
from pathlib import Path

from .config import (
    infer_workflow,
    priority_for,
    workflow_for_path,
)
from .schemas import RAGChunk

__all__ = ["index_doc_file", "split_sections"]

_MD_HEADING = re.compile(r"^(#{1,6})\s+(.*\S)\s*$")
_RST_UNDERLINE = re.compile(r"^[=\-~^\"'`#*+.]{3,}\s*$")
_MIN_SECTION_CHARS = 24


def split_sections(text: str, suffix: str) -> list[tuple[str, str]]:
    """Split *text* into ``(title, body)`` sections.

    Content before the first heading is returned under an empty title.
    """
    lines = text.splitlines()
    sections: list[tuple[str, list[str]]] = [("", [])]

    if suffix == ".md":
        for line in lines:
            m = _MD_HEADING.match(line)
            if m:
                sections.append((m.group(2).strip(), []))
            else:
                sections[-1][1].append(line)
    else:  # reStructuredText
        i = 0
        n = len(lines)
        while i < n:
            line = lines[i]
            nxt = lines[i + 1] if i + 1 < n else ""
            if (
                line.strip()
                and _RST_UNDERLINE.match(nxt)
                and len(nxt.strip()) >= len(line.strip()) - 2
            ):
                sections.append((line.strip(), []))
                i += 2
                continue
            sections[-1][1].append(line)
            i += 1

    out: list[tuple[str, str]] = []
    for title, body in sections:
        joined = "\n".join(body).strip()
        if joined or title:
            out.append((title, joined))
    return out


def index_doc_file(
    path: Path, root: Path, *, kind: str = "doc_section"
) -> list[RAGChunk]:
    """Index one ``.rst`` / ``.md`` file into section chunks.

    *kind* lets callers tag recipes (``"recipe"``) distinctly from
    ordinary documentation so the retriever can boost them.
    """
    rel_path = path.relative_to(root).as_posix()
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return []

    suffix = path.suffix.lower()
    prio = priority_for(rel_path)
    chunks: list[RAGChunk] = []

    for idx, (title, body) in enumerate(split_sections(text, suffix)):
        if len(body) < _MIN_SECTION_CHARS and not title:
            continue
        header = f"{title}\n\n" if title else ""
        chunk_text = f"Document: {rel_path}\n{header}{body}".strip()
        chunks.append(
            RAGChunk(
                id=f"{rel_path}:sec{idx}",
                text=chunk_text,
                source_path=rel_path,
                kind=kind,
                title=title or None,
                workflow=infer_workflow(rel_path, title, body[:400])
                or workflow_for_path(rel_path),
                priority=prio,
                metadata={"section_index": idx},
            )
        )

    return chunks
