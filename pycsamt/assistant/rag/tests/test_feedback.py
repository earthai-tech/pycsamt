# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Thumbs up/down feedback → retrieval hard negatives (RAG Tier 3)."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from pycsamt.assistant.rag.feedback import FeedbackStore, chunk_key
from pycsamt.assistant.rag.retriever import Retriever
from pycsamt.assistant.rag.schemas import RAGChunk

_QUERY = "correct galvanic static shift with ama"
_GOOD = "pycsamt.emtools.ss.correct_ss_ama"
_BAD = "pycsamt.zonge.ops.ASTATIC"


def _store() -> FeedbackStore:
    return FeedbackStore(Path(tempfile.mkdtemp()) / "feedback.jsonl")


def _chunks() -> list[RAGChunk]:
    # Both mention the query terms; only feedback should separate them.
    return [
        RAGChunk(
            id="bad", text="ASTATIC galvanic static shift ama correction",
            source_path="pycsamt/zonge/ops.py", kind="python_symbol",
            symbol=_BAD,
        ),
        RAGChunk(
            id="good", text="correct_ss_ama galvanic static shift ama correction",
            source_path="pycsamt/emtools/ss.py", kind="python_symbol",
            symbol=_GOOD,
        ),
    ]


class TestFeedbackStore(unittest.TestCase):

    def test_record_and_load_roundtrip(self):
        s = _store()
        s.record(_QUERY, _GOOD, True)
        s.record(_QUERY, _BAD, False)
        recs = s.load()
        self.assertEqual(len(recs), 2)
        self.assertTrue(recs[0]["helpful"])
        self.assertFalse(recs[1]["helpful"])

    def test_missing_file_loads_empty(self):
        self.assertEqual(_store().load(), [])
        self.assertEqual(_store().adjustments(_QUERY), {})

    def test_adjustments_boost_and_penalty(self):
        s = _store()
        s.record(_QUERY, _GOOD, True)
        s.record(_QUERY, _BAD, False)
        adj = s.adjustments(_QUERY)
        self.assertGreater(adj[_GOOD], 0.0)
        self.assertLess(adj[_BAD], 0.0)

    def test_dissimilar_query_gets_no_adjustment(self):
        s = _store()
        s.record(_QUERY, _GOOD, True)
        self.assertEqual(s.adjustments("plot the phase tensor map"), {})

    def test_repeated_verdicts_accumulate(self):
        s = _store()
        s.record(_QUERY, _BAD, False)
        one = s.adjustments(_QUERY)[_BAD]
        s.record(_QUERY, _BAD, False)
        two = s.adjustments(_QUERY)[_BAD]
        self.assertLess(two, one)  # more negative

    def test_blank_input_ignored(self):
        s = _store()
        s.record("", _GOOD, True)
        s.record(_QUERY, "", True)
        self.assertEqual(s.load(), [])


class TestChunkKey(unittest.TestCase):

    def test_code_chunk_keys_on_symbol(self):
        c = RAGChunk(id="x", text="", source_path="a.py",
                     kind="python_symbol", symbol=_GOOD)
        self.assertEqual(chunk_key(c), _GOOD)

    def test_doc_chunk_falls_back_to_id(self):
        # Doc/recipe chunks have no symbol; they must still be votable.
        c = RAGChunk(id="docs/x.rst:sec3", text="", source_path="docs/x.rst",
                     kind="doc_section")
        self.assertEqual(chunk_key(c), "docs/x.rst:sec3")


class TestRetrieverFeedbackIntegration(unittest.TestCase):
    """A rejected chunk is demoted on a similar query — but never erased."""

    def test_hard_negative_demotes_symbol(self):
        chunks = _chunks()
        baseline = Retriever(chunks).search(_QUERY, k=2)
        self.assertEqual(len(baseline.chunks), 2)

        s = _store()
        s.record(_QUERY, chunk_key(baseline.chunks[0]), False)  # reject winner
        s.record(_QUERY, chunk_key(baseline.chunks[1]), True)   # promote other

        tuned = Retriever(chunks, feedback_adjust=s.adjustments).search(_QUERY, k=2)
        self.assertEqual(tuned.chunks[0].symbol, baseline.chunks[1].symbol)
        # demoted, not removed
        self.assertEqual(len(tuned.chunks), 2)
        self.assertIn(baseline.chunks[0].symbol, [c.symbol for c in tuned.chunks])

    def test_doc_chunk_can_be_demoted(self):
        # A symbol-less doc section is the top hit; rejecting it must work.
        docs = [
            RAGChunk(id="doc:sec1", text="galvanic static shift ama overview",
                     source_path="docs/ss.rst", kind="doc_section", priority=3),
            RAGChunk(id="code", text="galvanic static shift ama correction",
                     source_path="pycsamt/emtools/ss.py",
                     kind="python_symbol", symbol=_GOOD),
        ]
        baseline = Retriever(docs).search(_QUERY, k=2)
        top = baseline.chunks[0]
        self.assertIsNone(top.symbol)  # the doc section leads

        s = _store()
        s.record(_QUERY, chunk_key(top), False)
        tuned = Retriever(docs, feedback_adjust=s.adjustments).search(_QUERY, k=2)
        self.assertNotEqual(tuned.chunks[0].id, top.id)

    def test_no_feedback_leaves_order_unchanged(self):
        chunks = _chunks()
        a = Retriever(chunks).search(_QUERY, k=2)
        b = Retriever(chunks, feedback_adjust=_store().adjustments).search(_QUERY, k=2)
        self.assertEqual([c.id for c in a.chunks], [c.id for c in b.chunks])

    def test_broken_adjuster_does_not_break_search(self):
        def boom(_q):
            raise RuntimeError("feedback store corrupt")

        ctx = Retriever(_chunks(), feedback_adjust=boom).search(_QUERY, k=2)
        self.assertEqual(len(ctx.chunks), 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
