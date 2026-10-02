# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the optional dense stage: fusion, cosine, vector store (Tier 1).

Uses numpy (present in the env) but never any embedding service — the
backend resolution path is checked to *decline* without a key.
"""

from __future__ import annotations

import sys
import tempfile
import types
import unittest
import unittest.mock
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pycsamt.assistant.rag.embeddings import (
    OpenAIEmbeddingBackend,
    _l2_normalize,
    cosine_scores,
    load_vectors,
    resolve_embedding_backend,
    rrf_fuse,
    save_vectors,
)


class TestRRF(unittest.TestCase):
    def test_fuses_and_rewards_agreement(self):
        # Item 1 is rank-0 in both lists → must score highest.
        fused = rrf_fuse([[1, 0, 2], [1, 2, 0]])
        best = max(fused, key=lambda i: fused[i])
        self.assertEqual(best, 1)

    def test_weights_bias_a_ranker(self):
        base = rrf_fuse([[0], [1]])
        self.assertAlmostEqual(base[0], base[1])  # symmetric
        weighted = rrf_fuse([[0], [1]], weights=[3.0, 1.0])
        self.assertGreater(weighted[0], weighted[1])  # ranker 0 wins

    def test_absent_items_contribute_nothing(self):
        fused = rrf_fuse([[0, 1]])
        self.assertNotIn(2, fused)


class TestCosine(unittest.TestCase):
    def test_identical_vectors_score_one(self):
        mat = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
        s = cosine_scores(np.array([1.0, 0.0]), mat)
        self.assertAlmostEqual(float(s[0]), 1.0, places=5)
        self.assertAlmostEqual(float(s[1]), 0.0, places=5)

    def test_handles_unnormalised_query(self):
        mat = np.array([[1.0, 0.0]], dtype=np.float32)
        s = cosine_scores(np.array([5.0, 0.0]), mat)  # not unit-norm
        self.assertAlmostEqual(float(s[0]), 1.0, places=5)

    def test_zero_query_vector_skips_normalization(self):
        mat = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
        s = cosine_scores(np.array([0.0, 0.0]), mat)
        self.assertTrue(np.allclose(s, [0.0, 0.0]))


class TestL2Normalize(unittest.TestCase):
    def test_normalizes_rows_and_guards_zero_rows(self):
        mat = np.array([[3.0, 4.0], [0.0, 0.0]], dtype=np.float32)
        out = _l2_normalize(mat)
        self.assertAlmostEqual(float(np.linalg.norm(out[0])), 1.0, places=5)
        # zero row is left as all-zero, not divided by zero
        self.assertTrue(np.allclose(out[1], [0.0, 0.0]))


class TestVectorStore(unittest.TestCase):
    def test_roundtrip(self):
        mat = np.random.RandomState(0).rand(4, 8).astype(np.float32)
        p = Path(tempfile.mkdtemp()) / "e.npz"
        save_vectors(p, ["a", "b", "c", "d"], mat)
        got = load_vectors(p)
        self.assertIsNotNone(got)
        ids, m2 = got
        self.assertEqual(ids, ["a", "b", "c", "d"])
        self.assertTrue(np.allclose(m2, mat))

    def test_missing_file_returns_none(self):
        self.assertIsNone(load_vectors(Path(tempfile.mkdtemp()) / "nope.npz"))

    def test_corrupt_file_returns_none(self):
        p = Path(tempfile.mkdtemp()) / "bad.npz"
        p.write_bytes(b"not a valid npz payload")
        self.assertIsNone(load_vectors(p))


class TestBackendResolution(unittest.TestCase):
    def test_no_key_declines(self):
        # Without a key, dense retrieval must stay off (returns None).
        self.assertIsNone(resolve_embedding_backend(api_key=None))

    def test_unknown_provider_declines(self):
        self.assertIsNone(resolve_embedding_backend(api_key="x", provider="nonesuch"))

    def test_declines_when_openai_import_fails(self):
        with unittest.mock.patch.dict(sys.modules, {"openai": None}):
            self.assertIsNone(
                resolve_embedding_backend(api_key="x", provider="openai")
            )

    def test_resolves_openai_backend_when_importable(self):
        fake_module = types.ModuleType("openai")
        fake_module.OpenAI = object  # never instantiated on this path
        with unittest.mock.patch.dict(sys.modules, {"openai": fake_module}):
            backend = resolve_embedding_backend(
                api_key="key123", provider=None, model="my-model"
            )
        self.assertIsInstance(backend, OpenAIEmbeddingBackend)
        self.assertEqual(backend.api_key, "key123")
        self.assertEqual(backend.model, "my-model")
        self.assertEqual(backend.name, "openai:my-model")


class _FakeEmbeddingsAPI:
    def __init__(self, dim=2):
        self.dim = dim
        self.calls: list[list[str]] = []

    def create(self, model, input):
        self.calls.append(list(input))
        return SimpleNamespace(
            data=[SimpleNamespace(embedding=[1.0] * self.dim) for _ in input]
        )


class _FakeOpenAIClient:
    def __init__(self, api_key=None):
        self.api_key = api_key
        self.embeddings = _FakeEmbeddingsAPI()


class TestOpenAIEmbeddingBackend(unittest.TestCase):
    def test_embed_batches_replaces_empty_strings_and_normalizes(self):
        fake_module = types.ModuleType("openai")
        fake_module.OpenAI = _FakeOpenAIClient
        with unittest.mock.patch.dict(sys.modules, {"openai": fake_module}):
            backend = OpenAIEmbeddingBackend(
                "key123", model="test-model", batch_size=2
            )
            self.assertEqual(backend.name, "openai:test-model")
            vecs = backend.embed(["a", "", "c"])

        self.assertEqual(vecs.shape, (3, 2))
        # rows [1.0, 1.0] normalized to unit length
        expected = 1.0 / np.sqrt(2.0)
        self.assertTrue(np.allclose(vecs, expected))


if __name__ == "__main__":
    unittest.main(verbosity=2)
