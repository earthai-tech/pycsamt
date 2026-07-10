# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the context builder and its wiring into PackageQAAgent
(RAG Step 4). Hermetic — small in-memory corpus, fake registry,
monkeypatched default builder.
"""
from __future__ import annotations

import unittest

from pycsamt.assistant.rag import context_builder as cb
from pycsamt.assistant.rag.context_builder import (
    AssembledContext,
    ContextBuilder,
    needs_clarification,
)
from pycsamt.assistant.rag.retriever import Retriever
from pycsamt.assistant.rag.schemas import RAGChunk


def _chunks() -> list[RAGChunk]:
    return [
        RAGChunk(
            id="1",
            text="Symbol: pycsamt.emtools.ss.estimate_ss_ama\n"
            "Signature: estimate_ss_ama(sites, half_window=3)\n"
            "File: pycsamt/emtools/ss.py\n\n"
            "Docstring:\nEstimate AMA static-shift factors per station.\n\n"
            "Code:\ndef estimate_ss_ama(...): ...",
            source_path="pycsamt/emtools/ss.py",
            kind="python_symbol", workflow="static_shift", priority=3,
            symbol="pycsamt.emtools.ss.estimate_ss_ama",
            metadata={"signature": "estimate_ss_ama(sites, half_window=3)"},
        ),
        RAGChunk(
            id="2",
            text="Static shift correction recipe. Use StaticShiftAgent or "
            "estimate_ss_ama / correct_ss_ama. galvanic distortion.",
            source_path="assistant_recipes/static_shift.md",
            kind="recipe", workflow="static_shift", priority=3,
            title="Static-shift recipe",
        ),
    ]


class _FakeRegistry:
    def find_line_in_text(self, text):
        return "L22PLT" if "l22plt" in text.lower() else None

    def resolve_line(self, name):
        return {
            "line": "L22PLT", "edi_dir": "/data/L22PLT", "exists": True,
            "n_edi_files": 25, "sort_by": "lon",
            "output_root": "results/willy/L22PLT",
            "default_workflows": ["load", "qc", "static_shift"],
            "static_shift": {"method": "ama", "half_window": 3},
            "plot": {"dpi": 150},
        }


class TestContextBuilder(unittest.TestCase):

    def setUp(self):
        self.builder = ContextBuilder(
            Retriever(_chunks()), _FakeRegistry()
        )

    def test_rerank_fn_reorders_chunks(self):
        # A rerank_fn that names candidates "2, 1" reverses the two-chunk
        # pool, whatever order retrieval produced.
        default_ids = [c.id for c in self.builder.build("static shift").chunks]
        self.assertEqual(len(default_ids), 2)
        reranked = self.builder.build(
            "static shift", rerank_fn=lambda prompt, system: "2, 1"
        )
        self.assertEqual(
            [c.id for c in reranked.chunks], list(reversed(default_ids))
        )

    def test_rerank_absent_leaves_default_path(self):
        # No rerank_fn → identical to a plain build (no crash, same lead).
        a = self.builder.build("static shift")
        b = self.builder.build("static shift", rerank_fn=None)
        self.assertEqual(
            [c.id for c in a.chunks], [c.id for c in b.chunks]
        )

    def test_session_rewrites_followup_for_retrieval(self):
        # A subject-less follow-up retrieves nothing on its own, but with a
        # session it inherits "static_shift" and finds the recipe/symbol.
        cold = self.builder.build("run that again")
        self.assertTrue(cold.is_empty())
        warm = self.builder.build(
            "run that again", session={"last_workflow": "static_shift"}
        )
        self.assertTrue(warm.chunks)
        # The user-facing query is preserved (only retrieval was rewritten).
        self.assertEqual(warm.query, "run that again")

    def test_top_score_populated(self):
        ac = self.builder.build("static shift correction")
        self.assertGreater(ac.top_score, 0.0)

    def test_build_has_context_and_citations(self):
        ac = self.builder.build("how to correct static shift")
        self.assertFalse(ac.is_empty())
        self.assertIn("estimate_ss_ama", ac.context_text)
        self.assertTrue(ac.citations)
        self.assertEqual(ac.citations[0]["n"], 1)

    def test_project_context_resolved_from_line(self):
        ac = self.builder.build("process line L22PLT for static shift")
        self.assertEqual(ac.project_context.get("line"), "L22PLT")
        self.assertIn("L22PLT", ac.context_text)
        self.assertIn("edi_dir", ac.context_text)

    def test_no_project_context_without_line(self):
        ac = self.builder.build("how to correct static shift")
        self.assertEqual(ac.project_context, {})

    def test_offline_answer_composed(self):
        ac = self.builder.build("static shift")
        ans = ac.compose_offline_answer()
        # leads with a real symbol's docstring, ends with citations,
        # and has no empty code fences
        self.assertIn("estimate_ss_ama", ans)
        self.assertIn("Sources:", ans)
        self.assertNotIn("```", ans)

    def test_empty_context(self):
        empty = AssembledContext(query="x")
        self.assertTrue(empty.is_empty())
        self.assertIn("couldn't find", empty.compose_offline_answer())


class TestClarification(unittest.TestCase):

    def test_confident_context_needs_no_clarification(self):
        ac = AssembledContext(query="static shift", top_score=50.0)
        self.assertIsNone(needs_clarification(ac))

    def test_project_context_is_confident(self):
        ac = AssembledContext(
            query="do that again", top_score=1.0,
            project_context={"line": "L22PLT"},
        )
        self.assertIsNone(needs_clarification(ac))

    def test_anaphoric_followup_clarifies(self):
        ac = AssembledContext(query="do that again", top_score=5.0)
        msg = needs_clarification(ac)
        self.assertIsNotNone(msg)
        self.assertIn("workflow", msg.lower())

    def test_selfcontained_lowscore_query_does_not_clarify(self):
        # A real (if rare) question is answered best-effort, never bounced —
        # even at a low score, because it names its own subject.
        ac = AssembledContext(query="the Sites class", top_score=5.0)
        self.assertIsNone(needs_clarification(ac))


class TestApiCardsAndGraph(unittest.TestCase):
    """Tier 3: signature-aware cards + symbol cross-references in context."""

    def _chunks(self):
        return [
            RAGChunk(
                id="1",
                text="Symbol: pycsamt.emtools.ss.correct_ss_ama\n\n"
                "Docstring:\nCorrect static shift by AMA.",
                source_path="pycsamt/emtools/ss.py",
                kind="python_symbol", workflow="static_shift", priority=3,
                symbol="pycsamt.emtools.ss.correct_ss_ama",
                module="pycsamt.emtools.ss",
                metadata={
                    "signature": "correct_ss_ama(sites, half_window=3) -> dict",
                    "params": [
                        {"name": "sites", "annotation": None, "default": None},
                        {"name": "half_window", "annotation": "int",
                         "default": "3"},
                    ],
                    "returns": "dict",
                    "refs": ["ensure_sites"],
                },
            ),
            RAGChunk(
                id="2",
                text="Symbol: pycsamt.emtools.ss.ensure_sites\n\n"
                "Docstring:\nLoad EDIs into Sites.",
                source_path="pycsamt/emtools/ss.py",
                kind="python_symbol", priority=3,
                symbol="pycsamt.emtools.ss.ensure_sites",
                module="pycsamt.emtools.ss",
                metadata={"signature": "ensure_sites(path)", "params": [],
                          "refs": []},
            ),
        ]

    def setUp(self):
        from pycsamt.assistant.rag.graph import SymbolGraph

        chunks = self._chunks()
        self.builder = ContextBuilder(
            Retriever(chunks), None, SymbolGraph(chunks)
        )

    def test_api_card_rendered_for_lead_symbol(self):
        ac = self.builder.build("correct static shift ama")
        self.assertIn("API: correct_ss_ama(sites, half_window=3)", ac.context_text)
        self.assertIn("half_window: int = 3", ac.context_text)

    def test_related_symbols_from_graph(self):
        ac = self.builder.build("correct static shift ama")
        self.assertIn("pycsamt.emtools.ss.ensure_sites", ac.related_symbols)
        self.assertIn("Related API", ac.context_text)

    def test_api_cards_structured(self):
        ac = self.builder.build("correct static shift ama")
        cards = {c["symbol"]: c for c in ac.api_cards()}
        card = cards["pycsamt.emtools.ss.correct_ss_ama"]
        self.assertEqual(card["returns"], "dict")
        self.assertEqual(
            [p["name"] for p in card["params"]], ["sites", "half_window"]
        )

    def test_related_symbols_lead_the_see_also(self):
        ac = self.builder.build("correct static shift ama")
        ans = ac.compose_offline_answer()
        self.assertIn("ensure_sites", ans)

    def test_no_graph_means_no_related(self):
        plain = ContextBuilder(Retriever(self._chunks()), None)
        ac = plain.build("correct static shift ama")
        self.assertEqual(ac.related_symbols, [])
        self.assertNotIn("Related API", ac.context_text)


class TestPackageQAWithRAG(unittest.TestCase):
    """PackageQAAgent offline uses the RAG-composed answer when available."""

    def setUp(self):
        self._orig = cb.default_context_builder
        builder = ContextBuilder(Retriever(_chunks()), _FakeRegistry())
        cb.default_context_builder = lambda root=None: builder

    def tearDown(self):
        cb.default_context_builder = self._orig

    def test_offline_uses_rag(self):
        from pycsamt.agents.package_qa import PackageQAAgent
        qa = PackageQAAgent()  # no key → offline
        r = qa.execute({"question": "how to correct static shift"})
        self.assertEqual(r.status, "success")
        self.assertEqual(r.data["source"], "rag_offline")
        self.assertTrue(r.data["citations"])
        self.assertIn("static", r.data["answer"].lower())

    def test_use_rag_false_falls_back(self):
        from pycsamt.agents.package_qa import PackageQAAgent
        qa = PackageQAAgent(use_rag=False)
        r = qa.execute({"question": "what workflows are supported?"})
        # falls back to the docstring lookup, not rag_offline
        self.assertNotEqual(r.data["source"], "rag_offline")


if __name__ == "__main__":
    unittest.main(verbosity=2)
