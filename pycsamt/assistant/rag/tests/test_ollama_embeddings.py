# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Local Ollama embedding backend (no network: the HTTP call is faked)."""
from __future__ import annotations

import io
import json

import numpy as np
import pytest

from pycsamt.assistant.rag import embeddings as emb


class _Resp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_ollama_backend_batches_and_normalises(monkeypatch):
    calls = []

    def fake_urlopen(req, timeout=None):
        body = json.loads(req.data)
        calls.append(body)
        return _Resp(json.dumps(
            {"embeddings": [[3.0, 4.0] for _ in body["input"]]}).encode())

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    be = emb.OllamaEmbeddingBackend(model="nomic-embed-text", batch_size=2)
    out = be.embed(["a", "", "c"])
    assert be.name == "ollama:nomic-embed-text"
    assert [len(c["input"]) for c in calls] == [2, 1]
    assert calls[0]["input"][1] == " "  # empty text is not sent empty
    np.testing.assert_allclose(out, [[0.6, 0.8]] * 3, rtol=1e-6)


def test_ollama_backend_count_mismatch_raises(monkeypatch):
    monkeypatch.setattr(
        "urllib.request.urlopen",
        lambda req, timeout=None: _Resp(json.dumps({"embeddings": []}).encode()),
    )
    with pytest.raises(RuntimeError, match="returned 0 embeddings"):
        emb.OllamaEmbeddingBackend().embed(["x"])


def test_resolver_returns_ollama_without_key():
    be = emb.resolve_embedding_backend(provider="ollama", model="bge-m3")
    assert isinstance(be, emb.OllamaEmbeddingBackend)
    assert be.name == "ollama:bge-m3"
