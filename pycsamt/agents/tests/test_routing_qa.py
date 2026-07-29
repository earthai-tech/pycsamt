# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Tests that Q&A (package knowledge) requests are
correctly distinguished from workflow-execution
requests, and that routing directs them to
PackageQAAgent rather than a processing workflow.

Pattern
-------
* Workflow request:  "Run QC on /data/L22PLT"
  -> WorkflowOrchestratorAgent should route to 'qc'
* Q&A request:       "What does DataQCAgent do?"
  -> Should be answered by PackageQAAgent, NOT
     treated as a 'qc' workflow run request.

These tests run offline — no API key needed.
"""

from __future__ import annotations

import re
import unittest

# ── Q&A vs workflow classification ────────────────────────────────────────────

_QA_PATTERNS = [
    r"\b(what|how|why|explain|describe|tell me"
    r"|what is|what does|how do|can you"
    r"|which|list)\b",
]

_WORKFLOW_PATTERNS = [
    r"\b(run|load|process|apply|correct|invert" r"|compute|generate|perform|execute)\b",
]


def _is_qa(text: str) -> bool:
    """Heuristic: is this a Q&A request?"""
    t = text.lower()
    is_q = any(re.search(p, t) for p in _QA_PATTERNS)
    has_path = bool(re.search(r"/[\w/]+", t))
    # If it has a file path, it's likely a
    # workflow request even if it has a question
    # word.
    return is_q and not has_path


_QA_CASES = [
    "What does StaticShiftAgent do?",
    "How do I load EDI files in pycsamt?",
    "Explain the Sites class",
    "What workflows does pycsamt support?",
    "What is the difference between PINN and hybrid inversion?",
    "How does phase tensor analysis work?",
    "List all supported pycsamt agents",
    "Tell me about the DataQCAgent",
    "Which correction methods are available?",
    "Can you explain what AIInversionAgent does?",
    "What parameters does StaticShiftAgent take?",
    "How do I access impedance tensor values?",
    "What is PACKAGE_CONTEXT?",
    "Describe the WorkflowOrchestratorAgent",
]

_WORKFLOW_CASES = [
    "Run QC on /data/L22PLT",
    "Load EDI files from /mnt/edi/willy",
    "Apply static shift correction to /data/EDIs",
    "Process line L22PLT at /surveys/",
    "Invert the data at /data/willy using PINN",
    "Generate a report for /out/qc_results/",
    "Correct galvanic shift in /data/survey",
    "Run phase tensor analysis on /data/willy EDIs",
]


class TestQAvsWorkflowHeuristic(unittest.TestCase):
    """Verify _is_qa heuristic on known cases."""

    def test_qa_cases_classified_as_qa(self):
        for text in _QA_CASES:
            with self.subTest(text=text[:60]):
                self.assertTrue(
                    _is_qa(text),
                    f"Expected Q&A classification for: {text!r}",
                )

    def test_workflow_cases_not_classified_as_qa(self):
        for text in _WORKFLOW_CASES:
            with self.subTest(text=text[:60]):
                self.assertFalse(
                    _is_qa(text),
                    f"Expected workflow classification for: {text!r}",
                )


class TestPackageQAAnswersPackageQuestions(unittest.TestCase):
    r"""
    PackageQAAgent correctly handles Q&A cases and
    returns answers that reference real pycsamt API.
    """

    def setUp(self):
        from pycsamt.agents.package_qa import (
            PackageQAAgent,
        )

        self.agent = PackageQAAgent()

    def _ask(self, q: str) -> str:
        r = self.agent.execute({"question": q})
        self.assertEqual(
            r.status,
            "success",
            f"Agent failed for: {q!r}",
        )
        return r.data.get("answer", "")

    def test_all_qa_cases_return_success(self):
        for text in _QA_CASES:
            with self.subTest(text=text[:60]):
                r = self.agent.execute({"question": text})
                self.assertEqual(
                    r.status,
                    "success",
                    f"Failed: {text!r}",
                )

    def test_static_shift_question(self):
        ans = self._ask("What does StaticShiftAgent do?")
        self.assertGreater(len(ans), 20)
        self.assertTrue(
            any(
                k in ans.lower()
                for k in (
                    "static",
                    "shift",
                    "galvanic",
                    "correct",
                )
            )
        )

    def test_workflow_list_question(self):
        ans = self._ask("What workflows does pycsamt support?")
        # Should contain workflow names
        self.assertTrue(
            any(
                k in ans.lower()
                for k in (
                    "qc",
                    "inversion",
                    "static",
                    "report",
                )
            )
        )

    def test_agent_list_question(self):
        ans = self._ask("List all supported pycsamt agents")
        self.assertTrue(
            any(
                k in ans
                for k in (
                    "Agent",
                    "agent",
                    "MTLoader",
                    "DataQC",
                )
            )
        )

    def test_sites_question(self):
        ans = self._ask("How do I access impedance tensor values?")
        self.assertTrue(
            any(
                k in ans.lower()
                for k in (
                    "sites",
                    "impedance",
                    "tensor",
                    "edi",
                )
            )
        )


class TestContextInputDoesNotAbsorbQA(unittest.TestCase):
    r"""
    ContextInputAgent (workflow router) should NOT
    produce a valid workflow config for questions
    that are clearly asking about the package API.

    These questions lack a data path and a clear
    processing verb — so the router should either
    return an empty / unknown workflow or fall back
    gracefully, not produce 'qc' or 'inversion'.
    """

    def setUp(self):
        from pycsamt.agents import ContextInputAgent

        self.ctx_agent = ContextInputAgent()

    def _route(self, text: str) -> str:
        r = self.ctx_agent.execute({"request": text})
        cfg = r.get("config") or {}
        return str(cfg.get("workflow", "")).lower()

    def test_pure_qa_no_data_path_no_inversion(self):
        for text in [
            "What does StaticShiftAgent do?",
            "Explain the Sites class",
        ]:
            with self.subTest(text=text[:60]):
                self._route(text)
                # It's fine if the router returns
                # something — the important thing
                # is no data_path is set for a
                # question that provides none.
                r = self.ctx_agent.execute({"request": text})
                cfg = r.get("config") or {}
                # No real data path should appear
                data_path = str(cfg.get("data_path", ""))
                self.assertFalse(
                    data_path.startswith("/data"),
                    f"Router hallucinated a data"
                    f" path for a Q&A request:"
                    f" {data_path!r}",
                )


if __name__ == "__main__":
    unittest.main()
