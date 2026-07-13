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
    load_all_suites,
)
from pycsamt.assistant.rag.retriever import Retriever
from pycsamt.assistant.rag.schemas import RAGChunk

# ── thresholds (kept below observed over all suites: intent ~92 /
#    wf ~91 / line 100 / symbol-recall 100) so honest dips don't
#    false-fail, real regressions do ──
_MIN_INTENT = 0.85
_MIN_WORKFLOW = 0.85
_MIN_LINE = 0.90
_MIN_SYMBOL_RECALL = 0.80
# Retrieval-quality guards (Tier 2). Observed over all suites incl. the
# adversarial paraphrases: no_test 100 / wf_in_topk 100 / nonempty ~96
# (only the "what can you do?" META probes retrieve nothing, by design).
_MIN_WF_IN_TOPK = 0.90
_MIN_NONEMPTY = 0.85


class TestEvalSuiteGate(unittest.TestCase):
    """End-to-end gate over all suites + the real corpus + registry."""

    @classmethod
    def setUpClass(cls):
        cls.records = load_all_suites()
        cls.report = evaluate(cls.records)

    def test_suite_nonempty(self):
        self.assertGreaterEqual(len(self.records), 20)

    def test_intent_accuracy(self):
        self.assertGreaterEqual(
            self.report.metrics["intent_accuracy"],
            _MIN_INTENT,
            msg=self.report.summary(),
        )

    def test_workflow_accuracy(self):
        self.assertGreaterEqual(
            self.report.metrics["workflow_accuracy"],
            _MIN_WORKFLOW,
            msg=self.report.summary(),
        )

    def test_line_accuracy(self):
        self.assertGreaterEqual(
            self.report.metrics["line_accuracy"],
            _MIN_LINE,
            msg=self.report.summary(),
        )

    def test_symbol_recall(self):
        self.assertGreaterEqual(
            self.report.metrics["symbol_recall"],
            _MIN_SYMBOL_RECALL,
            msg=self.report.summary(),
        )

    def test_no_hallucination_violations(self):
        self.assertEqual(
            self.report.violations, [], msg=str(self.report.violations)
        )

    def test_no_test_file_pollution(self):
        # Tier 0/2 guard: a unit-test module must never be retrieved.
        self.assertEqual(
            self.report.test_pollution,
            [],
            msg=str(self.report.test_pollution),
        )

    def test_workflow_in_topk(self):
        # Tier 1/2 guard: paraphrased queries surface the right workflow.
        self.assertGreaterEqual(
            self.report.metrics["workflow_in_topk"],
            _MIN_WF_IN_TOPK,
            msg=self.report.summary(),
        )

    def test_retrieval_nonempty(self):
        self.assertGreaterEqual(
            self.report.metrics["retrieval_nonempty"],
            _MIN_NONEMPTY,
            msg=self.report.summary(),
        )


class TestHarnessHermetic(unittest.TestCase):
    """Harness mechanics on a tiny controlled corpus."""

    def _retriever(self):
        return Retriever(
            [
                RAGChunk(
                    id="1",
                    text="estimate_ss_ama AMA static shift factors",
                    source_path="pycsamt/emtools/ss.py",
                    kind="python_symbol",
                    workflow="static_shift",
                    priority=3,
                    symbol="pycsamt.emtools.ss.estimate_ss_ama",
                ),
            ]
        )

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

    def test_detects_test_file_pollution(self):
        # A retriever that surfaces a test module must be flagged.
        retr = Retriever(
            [
                RAGChunk(
                    id="t",
                    text="estimate_ss_ama static shift unit test",
                    source_path="pycsamt/emtools/tests/test_ss.py",
                    kind="python_symbol",
                    workflow="static_shift",
                    symbol="pycsamt.emtools.tests.test_ss.test_ama",
                ),
            ]
        )
        rep = evaluate(
            [{"query": "estimate_ss_ama static shift"}],
            retriever=retr,
            registry=None,
        )
        self.assertEqual(len(rep.test_pollution), 1)
        self.assertEqual(rep.metrics["no_test_in_topk"], 0.0)

    def test_workflow_in_topk_metric(self):
        recs = [
            {
                "query": "estimate_ss_ama static shift",
                "expected_retrieval_workflow": "static_shift",
            }
        ]
        rep = evaluate(recs, retriever=self._retriever(), registry=None)
        self.assertEqual(rep.metrics["workflow_in_topk"], 1.0)
        self.assertEqual(rep.metrics["no_test_in_topk"], 1.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
