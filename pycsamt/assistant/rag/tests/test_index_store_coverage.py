# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Targeted coverage for :mod:`pycsamt.assistant.rag.index_store`.

``test_index_store.py`` covers the text-only build/load/staleness happy
paths. This file adds: ``read_manifest`` on a directory with no manifest
and with a corrupt one, ``index_is_stale`` with no manifest at all,
``load_index`` on a corrupt ``chunks.jsonl``, and the ``embed=True``
branches of :func:`build_index` (both the "no backend resolved" error and
a successful embed via a monkeypatched backend, without any network
access or a real API key).
"""

from __future__ import annotations

from pathlib import Path
from tempfile import mkdtemp

import numpy as np
import pytest

from pycsamt.assistant.rag import embeddings as emb_mod
from pycsamt.assistant.rag.index_store import (
    build_index,
    index_is_stale,
    load_index,
    read_manifest,
)

_PY = '''\
"""Static-shift helpers."""

def estimate_ss_ama(sites, half_window=3):
    """Estimate AMA static-shift factors."""
    return sites
'''


def _tree() -> Path:
    tmp = Path(mkdtemp())
    (tmp / "pycsamt" / "emtools").mkdir(parents=True)
    (tmp / "pycsamt" / "emtools" / "ss.py").write_text(
        _PY, encoding="utf-8"
    )
    (tmp / "README.md").write_text(
        "# pyCSAMT\nProcessing suite.\n", encoding="utf-8"
    )
    return tmp


class _FakeBackend:
    name = "fake:unit-test"

    def embed(self, texts):
        return np.ones((len(texts), 4), dtype=np.float32)


def test_read_manifest_missing_returns_none(tmp_path):
    assert read_manifest(root=tmp_path) is None


def test_read_manifest_corrupt_json_returns_none(tmp_path):
    idx = tmp_path / ".pycsamt_rag"
    idx.mkdir()
    (idx / "manifest.json").write_text("{not valid json", encoding="utf-8")
    assert read_manifest(root=tmp_path) is None


def test_index_is_stale_true_when_no_manifest_at_all(tmp_path):
    assert index_is_stale(root=tmp_path) is True


def test_load_index_corrupt_chunks_returns_none(tmp_path):
    idx = tmp_path / ".pycsamt_rag"
    idx.mkdir()
    (idx / "chunks.jsonl").write_text("{not json\n", encoding="utf-8")
    assert load_index(root=tmp_path) is None


def test_build_index_embed_without_resolvable_backend_raises(monkeypatch):
    monkeypatch.setattr(
        emb_mod, "resolve_embedding_backend", lambda **kw: None
    )
    root = _tree()
    with pytest.raises(RuntimeError, match="no backend resolved"):
        build_index(root=root, embed=True, embed_api_key=None)


def test_build_index_embed_success_persists_vectors(monkeypatch, tmp_path):
    monkeypatch.setattr(
        emb_mod, "resolve_embedding_backend", lambda **kw: _FakeBackend()
    )
    root = _tree()
    out_dir = tmp_path / "idx"
    mf = build_index(
        root=root, out_dir=out_dir, embed=True, embed_api_key="unused"
    )
    assert mf["embedded"] is True
    assert mf["embed_model"] == "fake:unit-test"
    assert mf["embed_dim"] == 4
    assert (out_dir / emb_mod.VECTOR_FILENAME).is_file()
    assert (out_dir / "manifest.json").is_file()


def test_build_index_embed_no_save_skips_vector_persist(monkeypatch, tmp_path):
    monkeypatch.setattr(
        emb_mod, "resolve_embedding_backend", lambda **kw: _FakeBackend()
    )
    root = _tree()
    out_dir = tmp_path / "idx2"
    mf = build_index(
        root=root,
        out_dir=out_dir,
        embed=True,
        embed_api_key="unused",
        save=False,
    )
    assert mf["embedded"] is True
    assert not out_dir.exists()
