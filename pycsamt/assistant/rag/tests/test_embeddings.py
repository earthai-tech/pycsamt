# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the optional dense stage: fusion, cosine, vector store (Tier 1).

Uses numpy (present in the env) but never any embedding service — the
backend resolution path is checked to *decline* without a key.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from pycsamt.assistant.rag.embeddings import (
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


class TestBackendResolution(unittest.TestCase):
    def test_no_key_declines(self):
        # Without a key, dense retrieval must stay off (returns None).
        self.assertIsNone(resolve_embedding_backend(api_key=None))

    def test_unknown_provider_declines(self):
        self.assertIsNone(
            resolve_embedding_backend(api_key="x", provider="nonesuch")
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
