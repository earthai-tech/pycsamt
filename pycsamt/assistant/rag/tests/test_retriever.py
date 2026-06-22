# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the hybrid BM25 retriever and recipe ingestion (RAG Step 3).

Hermetic: builds a small in-memory chunk set; no repo scan, no deps.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from pycsamt.assistant.rag.ingest import build_chunks
from pycsamt.assistant.rag.retriever import BM25, Retriever, tokenize
from pycsamt.assistant.rag.schemas import RAGChunk


class TestTokenize(unittest.TestCase):

    def test_identifier_splitting(self):
        toks = set(tokenize("StaticShiftAgent estimate_ss_ama"))
        # camelCase + snake_case split into sub-words
        for t in ("static", "shift", "agent", "estimate", "ss", "ama"):
            self.assertIn(t, toks)
        # whole identifiers retained too
        self.assertIn("staticshiftagent", toks)

    def test_stopwords_and_short_dropped(self):
        toks = tokenize("how do I run the QC")
        self.assertNotIn("the", toks)
        self.assertNotIn("do", toks)
        self.assertIn("qc", toks)


class TestBM25(unittest.TestCase):

    def test_ranks_matching_doc_first(self):
        corpus = [
            tokenize("static shift galvanic correction ama"),
            tokenize("phase tensor strike dimensionality"),
            tokenize("load edi files into sites"),
        ]
        bm = BM25(corpus)
        scores = bm.scores(tokenize("galvanic static shift"))
        self.assertGreater(scores[0], scores[1])
        self.assertGreater(scores[0], scores[2])

    def test_no_match_zero(self):
        bm = BM25([tokenize("alpha beta")])
        self.assertEqual(bm.scores(tokenize("gamma")), [0.0])


def _chunks() -> list[RAGChunk]:
    return [
        RAGChunk(
            id="1", text="Static shift correction with AMA. "
            "estimate_ss_ama, correct_ss_ama, galvanic distortion.",
            source_path="assistant_recipes/static_shift.md",
            kind="recipe", workflow="static_shift", priority=3,
            title="Static-shift recipe",
        ),
        RAGChunk(
            id="2", text="StaticShiftAgent detects and corrects galvanic "
            "static shift in MT data.",
            source_path="pycsamt/agents/static_shift.py",
            kind="python_symbol", workflow="static_shift", priority=3,
            symbol="pycsamt.agents.static_shift.StaticShiftAgent",
            metadata={"signature": "class StaticShiftAgent"},
        ),
        RAGChunk(
            id="3", text="Phase tensor, strike and dimensionality analysis.",
            source_path="pycsamt/agents/phase_analysis.py",
            kind="python_symbol", workflow="phase_analysis", priority=3,
            symbol="pycsamt.agents.phase_analysis.PhaseAnalysisAgent",
        ),
        RAGChunk(
            id="4", text="Load EDI files into a Sites collection.",
            source_path="pycsamt/emtools/_core.py",
            kind="python_symbol", workflow=None, priority=3,
            symbol="pycsamt.emtools._core.ensure_sites",
        ),
    ]


class TestRetriever(unittest.TestCase):

    def setUp(self):
        self.r = Retriever(_chunks())

    def test_relevant_first(self):
        ctx = self.r.search("correct galvanic static shift", k=3)
        self.assertTrue(ctx.chunks)
        self.assertEqual(ctx.chunks[0].workflow, "static_shift")

    def test_symbol_boost(self):
        ctx = self.r.search("how do I use ensure_sites", k=2)
        self.assertEqual(
            ctx.chunks[0].symbol, "pycsamt.emtools._core.ensure_sites"
        )

    def test_kinds_filter(self):
        ctx = self.r.search(
            "static shift", k=5, kinds={"recipe"}
        )
        self.assertTrue(ctx.chunks)
        self.assertTrue(all(c.kind == "recipe" for c in ctx.chunks))

    def test_symbols_extracted(self):
        ctx = self.r.search("StaticShiftAgent", k=3)
        syms = [s["symbol"] for s in ctx.symbols]
        self.assertIn(
            "pycsamt.agents.static_shift.StaticShiftAgent", syms
        )

    def test_no_results_for_irrelevant(self):
        ctx = self.r.search("completely unrelated xyzzy query", k=3)
        self.assertEqual(ctx.chunks, [])


class TestRecipeIngestion(unittest.TestCase):

    def test_recipes_tagged_recipe_kind(self):
        tmp = Path(tempfile.mkdtemp())
        rec = tmp / "assistant_recipes"
        rec.mkdir(parents=True)
        (rec / "static_shift.md").write_text(
            "# Static shift\nCorrect galvanic static shift.\n",
            encoding="utf-8",
        )
        chunks = build_chunks(tmp)
        self.assertTrue(chunks)
        self.assertTrue(all(c.kind == "recipe" for c in chunks))
        self.assertTrue(
            any(c.workflow == "static_shift" for c in chunks)
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
