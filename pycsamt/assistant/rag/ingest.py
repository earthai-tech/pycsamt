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
>>> chunks = build_chunks(repo_root())  # doctest: +SKIP
>>> sum(c.kind == "python_symbol" for c in chunks)  # doctest: +SKIP
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Iterator
from pathlib import Path

from .ast_indexer import index_python_file
from .config import INDEX_ROOTS, ROOT_DOCS, should_index
from .doc_indexer import index_doc_file
from .schemas import RAGChunk

__all__ = [
    "repo_root",
    "iter_index_files",
    "build_chunks",
    "save_chunks",
    "load_chunks",
    "corpus_stats",
    "source_fingerprint",
]


def repo_root() -> Path:
    """Locate the repository root (the dir that contains ``pycsamt/``).

    Walks up from this file; falls back to four levels up.
    """
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "pycsamt").is_dir() and (parent / "pyproject.toml").exists():
            return parent
    return here.parents[3]


def iter_index_files(root: Path) -> Iterator[Path]:
    """Yield absolute paths of files that should be indexed under *root*.

    Walks only the :data:`~pycsamt.assistant.rag.config.INDEX_ROOTS` (plus
    the root-level docs) rather than the whole repository — ``data/``,
    ``results/`` and ``.git/`` hold tens of thousands of files that could
    never be indexed, and scanning them dominated both ingestion and the
    freshness check. The yielded set is unchanged.
    """
    root = Path(root)
    for name in ROOT_DOCS:
        path = root / name
        if path.is_file() and should_index(name):
            yield path
    for index_root in INDEX_ROOTS:
        base = root / index_root
        if not base.is_dir():
            continue
        for path in base.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix.lower() not in (".py", ".rst", ".md"):
                continue
            if should_index(path.relative_to(root).as_posix()):
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
                rel = path.relative_to(root).as_posix()
                kind = (
                    "recipe" if rel.startswith("assistant_recipes/") else "doc_section"
                )
                chunks.extend(index_doc_file(path, root, kind=kind))
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


def source_fingerprint(
    root: Path | str | None = None,
    *,
    files: Iterable[Path] | None = None,
) -> str:
    """Content fingerprint of the indexable source tree.

    Hashes each indexable file's ``(relpath, sha256(content))`` — no
    parsing — so it changes exactly when indexed content changes. Callers
    use it to detect a *stale* persisted index without re-ingesting.

    Content, not mtime: merely running the test suite (or any tool that
    rewrites a file byte-identically) bumps mtimes and would otherwise
    report a false "stale index". Hashing costs no more than the directory
    walk it replaces.
    """
    import hashlib

    root = Path(root) if root is not None else repo_root()
    paths = sorted(files) if files is not None else sorted(iter_index_files(root))
    h = hashlib.sha256()
    for p in paths:
        try:
            data = p.read_bytes()
        except OSError:
            continue
        rel = p.relative_to(root).as_posix()
        h.update(rel.encode("utf-8"))
        h.update(b":")
        h.update(hashlib.sha256(data).digest())
    return h.hexdigest()
