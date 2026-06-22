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
