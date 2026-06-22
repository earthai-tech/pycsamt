# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.ingest
============================

Walk the repository and build the RAG corpus — a flat list of
:class:`~pycsamt.assistant.rag.schemas.RAGChunk` — by dispatching each
indexable file to the AST or doc indexer.

This is the entry point for Step 1. Later steps consume the chunks to
build BM25 / vector indexes and a symbol graph; here we only produce the
chunks and (optionally) persist them as JSONL.

Examples
--------
>>> from pycsamt.assistant.rag.ingest import build_chunks, repo_root
>>> chunks = build_chunks(repo_root())          # doctest: +SKIP
>>> sum(c.kind == "python_symbol" for c in chunks)   # doctest: +SKIP
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Iterator

from .ast_indexer import index_python_file
from .config import should_index
from .doc_indexer import index_doc_file
from .schemas import RAGChunk

__all__ = [
    "repo_root",
    "iter_index_files",
    "build_chunks",
    "save_chunks",
    "load_chunks",
    "corpus_stats",
]


def repo_root() -> Path:
    """Locate the repository root (the dir that contains ``pycsamt/``).

    Walks up from this file; falls back to four levels up.
    """
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "pycsamt").is_dir() and (
            parent / "pyproject.toml"
        ).exists():
            return parent
    return here.parents[3]


def iter_index_files(root: Path) -> Iterator[Path]:
    """Yield absolute paths of files that should be indexed under *root*."""
    root = Path(root)
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() not in (".py", ".rst", ".md"):
            continue
        rel = path.relative_to(root).as_posix()
        if should_index(rel):
            yield path


def build_chunks(
    root: Path | str | None = None,
    *,
    files: Iterable[Path] | None = None,
) -> list[RAGChunk]:
    """Build the full chunk corpus.

    Parameters
    ----------
    root : path-like, optional
        Repo root; defaults to :func:`repo_root`.
    files : iterable of Path, optional
        Restrict ingestion to these files (useful for tests / incremental
        builds). When ``None``, all indexable files under *root* are used.
    """
    root = Path(root) if root is not None else repo_root()
    paths = list(files) if files is not None else iter_index_files(root)

    chunks: list[RAGChunk] = []
    for path in paths:
        suffix = path.suffix.lower()
        try:
            if suffix == ".py":
                chunks.extend(index_python_file(path, root))
            elif suffix in (".rst", ".md"):
                chunks.extend(index_doc_file(path, root))
        except Exception:  # noqa: BLE001 — ingestion is best-effort
            continue
    return chunks


def save_chunks(chunks: list[RAGChunk], out_path: Path | str) -> Path:
    """Write *chunks* as JSON Lines to *out_path* (parent dirs created)."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        for c in chunks:
            f.write(json.dumps(c.to_dict(), ensure_ascii=False) + "\n")
    return out_path


def load_chunks(in_path: Path | str) -> list[RAGChunk]:
    """Read chunks back from a JSONL file written by :func:`save_chunks`."""
    in_path = Path(in_path)
    out: list[RAGChunk] = []
    with in_path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                out.append(RAGChunk.from_dict(json.loads(line)))
    return out


def corpus_stats(chunks: list[RAGChunk]) -> dict[str, int]:
    """Return simple counts (total + per kind + per workflow) for reporting."""
    stats: dict[str, int] = {"total": len(chunks)}
    for c in chunks:
        stats[f"kind:{c.kind}"] = stats.get(f"kind:{c.kind}", 0) + 1
        if c.workflow:
            stats[f"wf:{c.workflow}"] = stats.get(f"wf:{c.workflow}", 0) + 1
    return stats
