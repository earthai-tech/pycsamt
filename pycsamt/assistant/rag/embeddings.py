# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.embeddings
================================

Optional **dense-retrieval** stage — embeddings + a persisted vector
store + fusion helpers — layered *behind* the BM25 retriever.

Design contract (why this module is small and dependency-light):

* **Never required.** The assistant works fully offline on BM25 +
  :mod:`~pycsamt.assistant.rag.expansion`. Everything here activates only
  when a caller has built corpus vectors *and* an embedding backend is
  available (an API key, or a local model someone plugs in).
* **Backend-agnostic.** :class:`EmbeddingBackend` is a tiny protocol;
  :class:`OpenAIEmbeddingBackend` is the default because the ``openai``
  client ships in the working env (Anthropic has no embeddings endpoint,
  so Claude users still embed via OpenAI or a local model). A local
  ``sentence-transformers`` backend can be dropped in without touching
  callers.
* **Rank fusion, not score fusion.** BM25 and cosine scores live on
  different scales, so we combine them with **Reciprocal Rank Fusion**
  (:func:`rrf_fuse`) — robust and normalisation-free.

The heavy import (``numpy``) is deferred to call time so importing this
module is free on installs that never use dense retrieval.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Protocol, Sequence, runtime_checkable

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np

__all__ = [
    "EmbeddingBackend",
    "OpenAIEmbeddingBackend",
    "resolve_embedding_backend",
    "save_vectors",
    "load_vectors",
    "cosine_scores",
    "rrf_fuse",
    "VECTOR_FILENAME",
]

VECTOR_FILENAME = "embeddings.npz"


@runtime_checkable
class EmbeddingBackend(Protocol):
    """Anything that turns texts into L2-normalised row vectors."""

    name: str

    def embed(self, texts: Sequence[str]) -> "np.ndarray":
        """Return a ``(len(texts), dim)`` float32 array, rows unit-norm."""
        ...


def _l2_normalize(mat: "np.ndarray") -> "np.ndarray":
    import numpy as np

    mat = np.asarray(mat, dtype=np.float32)
    norms = np.linalg.norm(mat, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    return mat / norms


class OpenAIEmbeddingBackend:
    """Embed via the OpenAI embeddings endpoint (default dense backend).

    Parameters
    ----------
    api_key : str
        OpenAI API key.
    model : str
        Embedding model id (default ``text-embedding-3-small`` — cheap,
        1536-dim, ample for code/doc chunks).
    batch_size : int
        Texts per request (the endpoint accepts many; we batch to stay
        well under token limits).
    """

    def __init__(
        self,
        api_key: str,
        *,
        model: str = "text-embedding-3-small",
        batch_size: int = 256,
    ) -> None:
        self.api_key = api_key
        self.model = model
        self.batch_size = batch_size
        self.name = f"openai:{model}"

    def embed(self, texts: Sequence[str]) -> "np.ndarray":
        import numpy as np
        import openai

        client = openai.OpenAI(api_key=self.api_key)
        vecs: list[list[float]] = []
        for i in range(0, len(texts), self.batch_size):
            batch = [t or " " for t in texts[i : i + self.batch_size]]
            resp = client.embeddings.create(model=self.model, input=batch)
            vecs.extend(d.embedding for d in resp.data)
        return _l2_normalize(np.asarray(vecs, dtype=np.float32))


def resolve_embedding_backend(
    *,
    api_key: str | None = None,
    provider: str | None = None,
    model: str | None = None,
) -> EmbeddingBackend | None:
    """Return an available :class:`EmbeddingBackend`, or ``None``.

    Currently resolves an OpenAI backend when an API key is supplied and
    the ``openai`` client imports. Returns ``None`` (dense retrieval
    stays off) otherwise — callers must degrade to BM25.
    """
    if not api_key:
        return None
    prov = (provider or "openai").lower()
    if prov in ("openai", "oai", ""):
        try:
            import openai  # noqa: F401
        except Exception:  # noqa: BLE001
            return None
        return OpenAIEmbeddingBackend(
            api_key, model=model or "text-embedding-3-small"
        )
    return None


# ── persisted vector store ──────────────────────────────────────────────────────

def save_vectors(
    path: Path | str, ids: Sequence[str], vectors: "np.ndarray"
) -> Path:
    """Persist ``ids`` + a ``(n, dim)`` matrix to ``path`` (``.npz``)."""
    import numpy as np

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        ids=np.asarray(list(ids), dtype=object),
        vectors=np.asarray(vectors, dtype=np.float32),
    )
    # np.savez appends .npz if absent; normalise the returned path
    return path if path.suffix == ".npz" else path.with_suffix(".npz")


def load_vectors(
    path: Path | str,
) -> "tuple[list[str], np.ndarray] | None":
    """Load ``(ids, vectors)`` written by :func:`save_vectors`, or ``None``."""
    import numpy as np

    path = Path(path)
    if not path.is_file():
        return None
    try:
        data = np.load(path, allow_pickle=True)
        ids = [str(x) for x in data["ids"].tolist()]
        vectors = np.asarray(data["vectors"], dtype=np.float32)
    except Exception:  # noqa: BLE001 — a corrupt cache must not crash startup
        return None
    return ids, vectors


# ── fusion ──────────────────────────────────────────────────────────────────────

def cosine_scores(qvec: "np.ndarray", matrix: "np.ndarray") -> "np.ndarray":
    """Cosine similarity of a query vector against every row of *matrix*.

    Assumes rows (and *qvec*) are L2-normalised, so cosine reduces to a
    dot product.
    """
    import numpy as np

    q = np.asarray(qvec, dtype=np.float32).ravel()
    n = np.linalg.norm(q)
    if n:
        q = q / n
    return np.asarray(matrix, dtype=np.float32) @ q


def rrf_fuse(
    rankings: Sequence[Sequence[int]],
    *,
    k: int = 60,
    weights: Sequence[float] | None = None,
) -> dict[int, float]:
    """Reciprocal-Rank-Fusion of several ranked id-lists.

    Each element of *rankings* is a list of item indices in descending
    relevance. An item's fused score is ``Σ_i w_i / (k + rank_i)`` (rank
    0-based). Items absent from a ranking simply contribute nothing from
    it. Returns ``{index: fused_score}``.
    """
    if weights is None:
        weights = [1.0] * len(rankings)
    fused: dict[int, float] = {}
    for ranking, w in zip(rankings, weights):
        for rank, idx in enumerate(ranking):
            fused[idx] = fused.get(idx, 0.0) + w / (k + rank)
    return fused
