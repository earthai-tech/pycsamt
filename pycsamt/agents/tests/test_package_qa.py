# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Tests for PackageQAAgent — offline and live modes.

Offline tests (class TestPackageQAOffline)
-------------------------------------------
* Run without any LLM API key — no cost, no network.
* Verify the agent returns structured answers that
  reference real pycsamt classes and workflows.
* These are the tests that verify the package is
  "understandable" by the AI assistant even without
  an LLM call — i.e. the context is correctly built.

Live LLM tests (class TestPackageQALive)
-----------------------------------------
* Require an API key in the environment variable
  PYCSAMT_TEST_API_KEY (or ANTHROPIC_API_KEY).
* Marked ``pytest.mark.live`` and skipped in CI.
* Verify that the LLM returns factually correct
  answers about pycsamt classes and workflows.

Run offline tests::

    pytest pycsamt/agents/tests/test_package_qa.py \
        -v -k "not live"

Run live tests (requires API key)::

    PYCSAMT_TEST_API_KEY=sk-ant-... \
    pytest pycsamt/agents/tests/test_package_qa.py \
        -v -m live
"""

from __future__ import annotations

import os
import unittest

import pytest

# ── fixtures / helpers ────────────────────────────────────────────────────────


def _api_key() -> str | None:
    return (
        os.environ.get("PYCSAMT_TEST_API_KEY")
        or os.environ.get("ANTHROPIC_API_KEY")
        or os.environ.get("OPENAI_API_KEY")
    )


def _provider() -> str:
    if os.environ.get("OPENAI_API_KEY"):
        return "openai"
    return "claude"


# ── offline tests ─────────────────────────────────────────────────────────────


class TestPackageQAOffline(unittest.TestCase):
    """All tests without LLM calls (offline mode)."""

    def setUp(self):
        from pycsamt.agents.package_qa import (
            PackageQAAgent,
        )

        self.agent = PackageQAAgent()  # no api_key

    # ── basic contract ────────────────────────────

    def test_execute_returns_success(self):
        r = self.agent.execute(
            {"question": ("What does StaticShiftAgent do?")}
        )
        self.assertEqual(r.status, "success")

    def test_answer_key_present(self):
        r = self.agent.execute({"question": "What is DataQCAgent?"})
        self.assertIn("answer", r.data)
        self.assertIsInstance(r.data["answer"], str)

    def test_answer_nonempty(self):
        r = self.agent.execute({"question": "How do I load EDI files?"})
        self.assertGreater(len(r.data.get("answer", "")), 10)

    def test_source_is_docstring_lookup(self):
        # Offline QA now prefers the RAG-composed answer when the corpus
        # is available; the docstring lookup remains the fallback and is
        # tested explicitly with RAG disabled.
        from pycsamt.agents.package_qa import PackageQAAgent

        agent = PackageQAAgent(use_rag=False)
        r = agent.execute({"question": "What is the Sites class?"})
        self.assertEqual(
            r.data.get("source"),
            "docstring_lookup",
        )

    def test_no_answer_for_empty_question(self):
        r = self.agent.execute({"question": ""})
        self.assertEqual(r.status, "failed")
        self.assertIn("question", r.error.lower())

    def test_uses_request_key_as_fallback(self):
        r = self.agent.execute(
            {"request": ("Tell me about StaticShiftAgent")}
        )
        self.assertEqual(r.status, "success")

    # ── content correctness (offline) ────────────

    def test_static_shift_mentions_correction(self):
        r = self.agent.execute(
            {"question": ("What does StaticShiftAgent do?")}
        )
        ans = r.data.get("answer", "").lower()
        self.assertTrue(
            any(
                kw in ans
                for kw in (
                    "shift",
                    "correct",
                    "static",
                    "galvanic",
                )
            ),
            f"Expected shift/correct in answer, got: {ans[:200]}",
        )

    def test_qc_mentions_quality(self):
        r = self.agent.execute(
            {"question": ("What is DataQCAgent used for?")}
        )
        ans = r.data.get("answer", "").lower()
        self.assertTrue(
            any(
                kw in ans
                for kw in (
                    "qc",
                    "quality",
                    "flag",
                    "station",
                    "data",
                )
            ),
            f"Expected qc/quality in answer, got: {ans[:200]}",
        )

    def test_sites_class_question(self):
        r = self.agent.execute({"question": "What is the Sites class?"})
        ans = r.data.get("answer", "").lower()
        self.assertTrue(
            any(
                kw in ans
                for kw in (
                    "site",
                    "edi",
                    "station",
                    "impedance",
                    "sites",
                )
            ),
            f"Expected sites/edi/station in answer, got: {ans[:200]}",
        )

    def test_loader_question(self):
        r = self.agent.execute(
            {"question": ("How do I load EDI files in pycsamt?")}
        )
        ans = r.data.get("answer", "").lower()
        self.assertTrue(
            any(
                kw in ans
                for kw in (
                    "mtloaderagent",
                    "loader",
                    "edi",
                    "load",
                )
            ),
            f"Expected loader/edi in answer, got: {ans[:200]}",
        )

    def test_workflow_question_mentions_workflow(self):
        r = self.agent.execute(
            {"question": ("What workflows does pycsamt support?")}
        )
        ans = r.data.get("answer", "").lower()
        # at least one workflow keyword present
        self.assertTrue(
            any(
                kw in ans
                for kw in (
                    "qc",
                    "inversion",
                    "static",
                    "phase",
                )
            ),
            f"No workflow keywords in answer: {ans[:200]}",
        )

    def test_inversion_question(self):
        r = self.agent.execute({"question": ("What is AIInversionAgent?")})
        ans = r.data.get("answer", "").lower()
        self.assertTrue(
            any(
                kw in ans
                for kw in (
                    "inver",
                    "neural",
                    "ai",
                    "network",
                    "model",
                )
            ),
            f"Expected inversion terms in answer, got: {ans[:200]}",
        )

    def test_pinn_question(self):
        r = self.agent.execute(
            {"question": ("How does PINN inversion work in pycsamt?")}
        )
        ans = r.data.get("answer", "").lower()
        self.assertTrue(
            any(
                kw in ans
                for kw in (
                    "pinn",
                    "physics",
                    "inver",
                    "neural",
                )
            ),
            f"Expected pinn/physics in answer, got: {ans[:200]}",
        )

    def test_excerpts_returned(self):
        r = self.agent.execute(
            {"question": ("What does StaticShiftAgent do?")}
        )
        excerpts = r.data.get("excerpts", [])
        self.assertIsInstance(excerpts, list)

    def test_unknown_question_graceful(self):
        r = self.agent.execute(
            {"question": ("How do I order pizza with pycsamt?")}
        )
        # Should not crash; returns success with
        # a fallback answer
        self.assertEqual(r.status, "success")
        self.assertGreater(len(r.data.get("answer", "")), 5)

    # ── system prompt contract ─────────────────────

    def test_system_prompt_contains_pycsamt(self):
        from pycsamt.agents.package_qa import (
            PackageQAAgent,
        )

        prompt = PackageQAAgent.SYSTEM_PROMPT
        self.assertIn("pycsamt", prompt)

    def test_system_prompt_contains_agent_names(self):
        from pycsamt.agents.package_qa import (
            PackageQAAgent,
        )

        prompt = PackageQAAgent.SYSTEM_PROMPT
        self.assertIn("StaticShiftAgent", prompt)
        self.assertIn("DataQCAgent", prompt)

    def test_system_prompt_contains_workflows(self):
        from pycsamt.agents.package_qa import (
            PackageQAAgent,
        )

        prompt = PackageQAAgent.SYSTEM_PROMPT
        self.assertIn("qc", prompt)
        self.assertIn("static_shift", prompt)

    def test_system_prompt_contains_examples(self):
        from pycsamt.agents.package_qa import (
            PackageQAAgent,
        )

        prompt = PackageQAAgent.SYSTEM_PROMPT
        self.assertIn("from pycsamt", prompt)


# ── live LLM tests ────────────────────────────────────────────────────────────


@pytest.mark.live
class TestPackageQALive(unittest.TestCase):
    r"""
    LLM-powered Q&A tests.

    Skipped automatically when no API key is set.
    Run with::

        ANTHROPIC_API_KEY=sk-ant-... pytest \
            pycsamt/agents/tests/test_package_qa.py \
            -m live -v
    """

    @classmethod
    def setUpClass(cls):
        key = _api_key()
        if not key:
            raise unittest.SkipTest(
                "No API key found in environment."
                " Set PYCSAMT_TEST_API_KEY or"
                " ANTHROPIC_API_KEY."
            )
        from pycsamt.agents.package_qa import (
            PackageQAAgent,
        )

        cls.agent = PackageQAAgent(
            api_key=key,
            llm_provider=_provider(),
        )

    def _ask(self, question: str) -> str:
        r = self.agent.execute({"question": question})
        self.assertEqual(
            r.status,
            "success",
            f"LLM call failed: {r.error}",
        )
        return r.data.get("answer", "")

    def test_static_shift_answer_is_grounded(self):
        ans = self._ask(
            "What does StaticShiftAgent do and"
            " what correction methods does it"
            " support?"
        )
        ans_lower = ans.lower()
        self.assertTrue(
            any(
                kw in ans_lower
                for kw in (
                    "ama",
                    "loess",
                    "median",
                    "static",
                    "galvanic",
                )
            ),
            f"Answer missing expected terms: {ans[:300]}",
        )

    def test_sites_class_described_correctly(self):
        ans = self._ask(
            "What is the Sites class and how do"
            " I access impedance data from it?"
        )
        ans_lower = ans.lower()
        self.assertTrue(
            any(
                kw in ans_lower
                for kw in (
                    "impedance",
                    "frequenc",
                    "resistiv",
                    "edi",
                )
            ),
            f"Answer missing impedance/freq terms: {ans[:300]}",
        )

    def test_workflow_listing_accurate(self):
        ans = self._ask("List all supported pycsamt workflows.")
        for wf in (
            "qc",
            "static_shift",
            "ai_inversion",
            "report",
        ):
            self.assertIn(
                wf,
                ans.lower(),
                f"Workflow {wf!r} missing from answer: {ans[:400]}",
            )

    def test_code_example_is_valid_python(self):
        ans = self._ask(
            "Show me a Python code example of"
            " loading EDI files and running QC"
            " with pycsamt."
        )
        # Should include at least one import
        self.assertIn(
            "import",
            ans.lower(),
            f"Expected import in code example: {ans[:300]}",
        )
        self.assertTrue(
            any(
                kw in ans
                for kw in (
                    "MTLoaderAgent",
                    "DataQCAgent",
                    "pycsamt",
                )
            ),
            f"Expected pycsamt class in example: {ans[:300]}",
        )

    def test_unknown_class_not_hallucinated(self):
        ans = self._ask("What does the QuantumInverterAgent do in pycsamt?")
        ans_lower = ans.lower()
        # LLM should say it doesn't exist,
        # not hallucinate a description
        self.assertTrue(
            any(
                kw in ans_lower
                for kw in (
                    "not",
                    "no such",
                    "does not exist",
                    "not found",
                    "cannot find",
                )
            ),
            f"LLM may have hallucinated an answer: {ans[:400]}",
        )

    def test_source_is_llm(self):
        r = self.agent.execute({"question": ("What is StaticShiftAgent?")})
        self.assertEqual(r.data.get("source"), "llm")

    def test_answer_references_pycsamt(self):
        ans = self._ask("How does pycsamt handle phase tensor analysis?")
        self.assertTrue(
            any(
                kw in ans.lower()
                for kw in (
                    "pycsamt",
                    "phase",
                    "tensor",
                    "phaseanaly",
                )
            ),
            f"Answer doesn't reference pycsamt or phase: {ans[:300]}",
        )


if __name__ == "__main__":
    unittest.main()
