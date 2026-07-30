# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.index_store
=================================

Persist the RAG corpus to disk so the retriever loads in milliseconds
instead of re-ingesting the repository (~seconds) on every process.

Layout (under ``<root>/.pycsamt_rag/`` by default)::

    .pycsamt_rag/
    ├── manifest.json     # version, timestamp, counts, stats
    └── chunks.jsonl      # one RAGChunk per line

Build it once (``python -m pycsamt.assistant.rag build``) or let the
builders create it lazily; refresh it whenever the package source
changes.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from .ingest import (
    build_chunks,
    corpus_stats,
    load_chunks,
    repo_root,
    save_chunks,
    source_fingerprint,
)
from .schemas import RAGChunk

__all__ = [
    "DEFAULT_INDEX_DIRNAME",
    "default_index_dir",
    "build_index",
    "load_index",
    "index_exists",
    "index_is_stale",
    "read_manifest",
]

DEFAULT_INDEX_DIRNAME = ".pycsamt_rag"
# v2: source_fingerprint switched from (mtime, size) to a content hash.
_MANIFEST_VERSION = 2


def default_index_dir(root: Path | str | None = None) -> Path:
    """``<root>/.pycsamt_rag`` (root defaults to the repo root)."""
    base = Path(root) if root is not None else repo_root()
    return base / DEFAULT_INDEX_DIRNAME


def index_exists(
    out_dir: Path | str | None = None, root: Path | str | None = None
) -> bool:
    """True when a persisted ``chunks.jsonl`` is present."""
    d = Path(out_dir) if out_dir is not None else default_index_dir(root)
    return (d / "chunks.jsonl").is_file()


def read_manifest(
    out_dir: Path | str | None = None, root: Path | str | None = None
) -> dict | None:
    """Return the manifest dict, or ``None`` if absent/unreadable."""
    d = Path(out_dir) if out_dir is not None else default_index_dir(root)
    mf = d / "manifest.json"
    if not mf.is_file():
        return None
    try:
        return json.loads(mf.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def index_is_stale(
    out_dir: Path | str | None = None, root: Path | str | None = None
) -> bool:
    """True when the persisted index no longer matches the source tree.

    Compares the manifest's stored :func:`~pycsamt.assistant.rag.ingest.
    source_fingerprint` against a freshly-computed one. A missing manifest
    or fingerprint (e.g. an index built before this field existed) counts
    as stale, so callers rebuild rather than trust an unknown state.
    """
    mf = read_manifest(out_dir, root)
    if not mf or "source_fingerprint" not in mf:
        return True
    if mf.get("version") != _MANIFEST_VERSION:
        return True  # fingerprint scheme changed → not comparable
    return mf["source_fingerprint"] != source_fingerprint(root)


def build_index(
    root: Path | str | None = None,
    out_dir: Path | str | None = None,
    *,
    save: bool = True,
    embed: bool = False,
    embed_api_key: str | None = None,
    embed_provider: str | None = None,
    embed_model: str | None = None,
) -> dict:
    """Ingest the repo and (optionally) persist chunks + manifest.

    When *embed* is true and an embedding backend resolves (an API key is
    supplied and its client imports), every chunk is embedded and the
    vectors are saved next to ``chunks.jsonl`` for dense retrieval. This
    is strictly opt-in — the default build is text-only and offline.

    Returns the manifest dict (with ``out_dir``, ``stats`` and, when
    embedded, ``embedded``/``embed_model``/``embed_dim``).
    """
    root = Path(root) if root is not None else repo_root()
    out_dir = Path(out_dir) if out_dir is not None else default_index_dir(root)

    chunks = build_chunks(root)
    stats = corpus_stats(chunks)
    manifest = {
        "version": _MANIFEST_VERSION,
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "root": str(root),
        "n_chunks": len(chunks),
        "stats": stats,
        "source_fingerprint": source_fingerprint(root),
        "embedded": False,
        "out_dir": str(out_dir),
    }

    vectors = None
    if embed:
        from .embeddings import resolve_embedding_backend

        backend = resolve_embedding_backend(
            api_key=embed_api_key,
            provider=embed_provider,
            model=embed_model,
        )
        if backend is None:
            raise RuntimeError(
                "Embedding requested but no backend resolved — pass an "
                "embed_api_key (and ensure the provider client is "
                "installed)."
            )
        vectors = backend.embed([c.text for c in chunks])
        manifest["embedded"] = True
        manifest["embed_model"] = backend.name
        manifest["embed_dim"] = int(vectors.shape[1])

    if save:
        out_dir.mkdir(parents=True, exist_ok=True)
        save_chunks(chunks, out_dir / "chunks.jsonl")
        if vectors is not None:
            from .embeddings import (
                VECTOR_FILENAME,
                save_vectors,
            )

            save_vectors(
                out_dir / VECTOR_FILENAME, [c.id for c in chunks], vectors
            )
        (out_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8"
        )
    return manifest


def load_index(
    out_dir: Path | str | None = None,
    root: Path | str | None = None,
) -> list[RAGChunk] | None:
    """Load persisted chunks, or ``None`` when no index is present."""
    d = Path(out_dir) if out_dir is not None else default_index_dir(root)
    chunks_path = d / "chunks.jsonl"
    if not chunks_path.is_file():
        return None
    try:
        return load_chunks(chunks_path)
    except (OSError, json.JSONDecodeError):
        return None
