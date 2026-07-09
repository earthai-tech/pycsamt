# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the optional LLM re-ranker (RAG Tier 1). No network."""
from __future__ import annotations

import unittest

from pycsamt.assistant.rag.rerank import (
    llm_rerank,
    parse_ranking,
    build_rerank_prompt,
)
from pycsamt.assistant.rag.schemas import RAGChunk


def _chunks() -> list[RAGChunk]:
    return [
        RAGChunk(id="a", text="alpha", source_path="a.py", kind="python_symbol"),
        RAGChunk(id="b", text="bravo", source_path="b.py", kind="python_symbol"),
        RAGChunk(id="c", text="charlie", source_path="c.py", kind="python_symbol"),
    ]


class TestParseRanking(unittest.TestCase):

    def test_parses_and_bounds(self):
        self.assertEqual(parse_ranking("3, 1, 5", 3), [2, 0])  # 5 out of range

    def test_dedupes(self):
        self.assertEqual(parse_ranking("2 2 1", 3), [1, 0])

    def test_garbage_is_empty(self):
        self.assertEqual(parse_ranking("no numbers here", 3), [])


class TestLLMRerank(unittest.TestCase):

    def test_reorders_by_model_reply(self):
        # Model says candidate 3 then 1 → chunk c, then a; b appended.
        reply = "3, 1"
        out = llm_rerank(
            "q", _chunks(), rank_fn=lambda p, s: reply
        )
        self.assertEqual([c.id for c in out], ["c", "a", "b"])

    def test_prompt_lists_candidates(self):
        prompt = build_rerank_prompt("my query", _chunks())
        self.assertIn("my query", prompt)
        self.assertIn("[1]", prompt)
        self.assertIn("[3]", prompt)

    def test_unparseable_reply_keeps_order(self):
        out = llm_rerank("q", _chunks(), rank_fn=lambda p, s: "???")
        self.assertEqual([c.id for c in out], ["a", "b", "c"])

    def test_exception_falls_back(self):
        def boom(p, s):
            raise RuntimeError("llm down")

        out = llm_rerank("q", _chunks(), rank_fn=boom)
        self.assertEqual([c.id for c in out], ["a", "b", "c"])

    def test_empty_input(self):
        self.assertEqual(llm_rerank("q", [], rank_fn=lambda p, s: "1"), [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
