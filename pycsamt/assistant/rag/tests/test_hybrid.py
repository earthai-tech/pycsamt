# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Retriever Tier 1: expansion integration + dense fusion (stub backend).

Hermetic — a deterministic hashing stub stands in for a real embedding
service, so the fusion mechanism is proven with no network or heavy deps.
"""
from __future__ import annotations

import unittest

import numpy as np

from pycsamt.assistant.rag.retriever import Retriever
from pycsamt.assistant.rag.schemas import RAGChunk


class _StubBackend:
    """Deterministic bag-of-words hashing embedder (proves the mechanism)."""

    name = "stub"
    dim = 64

    def _vec(self, text: str) -> np.ndarray:
        v = np.zeros(self.dim, dtype=np.float32)
        for w in text.lower().split():
            v[hash(w) % self.dim] += 1.0
        n = np.linalg.norm(v)
        return v / n if n else v

    def embed(self, texts):
        return np.vstack([self._vec(t) for t in texts])


class _BoomBackend:
    name = "boom"

    def embed(self, texts):  # noqa: D401
        raise RuntimeError("embedding service down")


def _chunks() -> list[RAGChunk]:
    return [
        RAGChunk(
            id="ss", text="static shift galvanic ama correction of impedance",
            source_path="pycsamt/emtools/ss.py", kind="python_symbol",
            workflow="static_shift",
            symbol="pycsamt.emtools.ss.correct_ss_ama",
        ),
        RAGChunk(
            id="dn", text="rpca hampel outlier spatial denoising removal",
            source_path="pycsamt/agents/denoising.py", kind="python_symbol",
            workflow="denoise",
            symbol="pycsamt.agents.denoising.DenoisingAgent",
        ),
        RAGChunk(
            id="ph", text="phase tensor strike dimensionality skew analysis",
            source_path="pycsamt/agents/phase_analysis.py",
            kind="python_symbol", workflow="phase_analysis",
            symbol="pycsamt.agents.phase_analysis.PhaseAnalysisAgent",
        ),
    ]


class TestExpansionIntegration(unittest.TestCase):

    def test_expansion_surfaces_chunk_with_no_lexical_overlap(self):
        # "spikes" shares NO word with the denoise chunk (rpca/hampel/…);
        # only query expansion can bridge them.
        r = Retriever(_chunks())
        ctx = r.search("clean up spikes", k=1)
        self.assertTrue(ctx.chunks)
        self.assertEqual(ctx.chunks[0].id, "dn")

    def test_expansion_off_reverts_to_lexical(self):
        r = Retriever(_chunks(), use_expansion=False)
        ctx = r.search("clean up spikes", k=3)
        # No lexical overlap and no expansion → nothing retrieved.
        self.assertEqual(ctx.chunks, [])


class TestDenseFusion(unittest.TestCase):

    def setUp(self):
        self.chunks = _chunks()
        self.mat = _StubBackend().embed([c.text for c in self.chunks])

    def test_dense_enabled_flag(self):
        self.assertFalse(Retriever(self.chunks).dense_enabled)
        r = Retriever(self.chunks, vectors=self.mat, embed_backend=_StubBackend())
        self.assertTrue(r.dense_enabled)

    def test_fusion_returns_relevant_first(self):
        r = Retriever(self.chunks, vectors=self.mat, embed_backend=_StubBackend())
        ctx = r.search("galvanic static shift correction", k=2)
        self.assertEqual(ctx.chunks[0].id, "ss")

    def test_kind_filter_respected_under_fusion(self):
        r = Retriever(self.chunks, vectors=self.mat, embed_backend=_StubBackend())
        ctx = r.search("phase tensor", k=3, kinds={"python_symbol"})
        self.assertTrue(all(c.kind == "python_symbol" for c in ctx.chunks))

    def test_embedding_failure_falls_back_to_lexical(self):
        # A broken backend must not break search — lexical order survives.
        r = Retriever(self.chunks, vectors=self.mat, embed_backend=_BoomBackend())
        ctx = r.search("galvanic static shift", k=2)
        self.assertTrue(ctx.chunks)
        self.assertEqual(ctx.chunks[0].id, "ss")


if __name__ == "__main__":
    unittest.main(verbosity=2)
