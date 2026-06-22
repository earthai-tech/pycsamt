# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
RAG eval suite gate (Step 8).

Runs the bundled suite against the real retriever + registry and asserts
minimum quality thresholds, so regressions in intent/workflow/line
classification or retrieval recall fail CI. Hermetic harness tests use a
tiny in-memory retriever.
"""
from __future__ import annotations

import unittest

from pycsamt.assistant.evals.harness import (
    EvalReport,
    evaluate,
    load_suite,
    suites_dir,
)
from pycsamt.assistant.rag.retriever import Retriever
from pycsamt.assistant.rag.schemas import RAGChunk


# ── thresholds (kept below observed: intent 90 / wf 100 / line 100 /
#    symbol-recall 100) so honest dips don't false-fail, real ones do ──
_MIN_INTENT = 0.80
_MIN_WORKFLOW = 0.90
_MIN_LINE = 0.90
_MIN_SYMBOL_RECALL = 0.70


class TestEvalSuiteGate(unittest.TestCase):
    """End-to-end gate over the real corpus + project registry."""

    @classmethod
    def setUpClass(cls):
        cls.records = load_suite(suites_dir() / "rag_questions.jsonl")
        cls.report = evaluate(cls.records)

    def test_suite_nonempty(self):
        self.assertGreaterEqual(len(self.records), 5)

    def test_intent_accuracy(self):
        self.assertGreaterEqual(
            self.report.metrics["intent_accuracy"], _MIN_INTENT,
            msg=self.report.summary(),
        )

    def test_workflow_accuracy(self):
        self.assertGreaterEqual(
            self.report.metrics["workflow_accuracy"], _MIN_WORKFLOW,
            msg=self.report.summary(),
        )

    def test_line_accuracy(self):
        self.assertGreaterEqual(
            self.report.metrics["line_accuracy"], _MIN_LINE,
            msg=self.report.summary(),
        )

    def test_symbol_recall(self):
        self.assertGreaterEqual(
            self.report.metrics["symbol_recall"], _MIN_SYMBOL_RECALL,
            msg=self.report.summary(),
        )

    def test_no_hallucination_violations(self):
        self.assertEqual(
            self.report.violations, [], msg=str(self.report.violations)
        )


class TestHarnessHermetic(unittest.TestCase):
    """Harness mechanics on a tiny controlled corpus."""

    def _retriever(self):
        return Retriever([
            RAGChunk(
                id="1", text="estimate_ss_ama AMA static shift factors",
                source_path="pycsamt/emtools/ss.py", kind="python_symbol",
                workflow="static_shift", priority=3,
                symbol="pycsamt.emtools.ss.estimate_ss_ama",
            ),
        ])

    def test_symbol_recall_and_violation(self):
        recs = [
            {
                "query": "estimate_ss_ama static shift",
                "expected_intent": "workflow",
                "expected_symbols": ["pycsamt.emtools.ss.estimate_ss_ama"],
                "must_not_contain": ["from pycsamt import static_shift"],
            }
        ]
        rep = evaluate(recs, retriever=self._retriever(), registry=None)
        self.assertEqual(rep.metrics["symbol_recall"], 1.0)
        self.assertEqual(rep.violations, [])

    def test_missing_symbol_lowers_recall(self):
        recs = [
            {
                "query": "estimate_ss_ama",
                "expected_symbols": [
                    "pycsamt.emtools.ss.estimate_ss_ama",
                    "pycsamt.does.not.Exist",
                ],
            }
        ]
        rep = evaluate(recs, retriever=self._retriever(), registry=None)
        self.assertEqual(rep.metrics["symbol_recall"], 0.5)

    def test_report_is_dataclass(self):
        rep = evaluate([], retriever=self._retriever(), registry=None)
        self.assertIsInstance(rep, EvalReport)
        self.assertEqual(rep.n, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
