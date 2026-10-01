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
        r = self.agent.execute({"question": ("What does StaticShiftAgent do?")})
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
        r = self.agent.execute({"request": ("Tell me about StaticShiftAgent")})
        self.assertEqual(r.status, "success")

    # ── content correctness (offline) ────────────

    def test_static_shift_mentions_correction(self):
        r = self.agent.execute({"question": ("What does StaticShiftAgent do?")})
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
        r = self.agent.execute({"question": ("What is DataQCAgent used for?")})
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
        r = self.agent.execute({"question": ("How do I load EDI files in pycsamt?")})
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
        r = self.agent.execute({"question": ("What workflows does pycsamt support?")})
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
        r = self.agent.execute({"question": ("What does StaticShiftAgent do?")})
        excerpts = r.data.get("excerpts", [])
        self.assertIsInstance(excerpts, list)

    def test_unknown_question_graceful(self):
        r = self.agent.execute({"question": ("How do I order pizza with pycsamt?")})
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
            "What is the Sites class and how do" " I access impedance data from it?"
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
        # "llm+rag" is expected whenever retrieval finds relevant context
        # (the normal case); "llm" is the fallback with no RAG corpus hit.
        # Either confirms a real LLM call was made, as opposed to an
        # offline/docstring-only answer.
        r = self.agent.execute({"question": ("What is StaticShiftAgent?")})
        self.assertIn(r.data.get("source"), ("llm", "llm+rag"))

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


# ── tier selection (private helpers) ─────────────────────────────────────────


class TestSelectTiers(unittest.TestCase):
    """Direct coverage of ``_select_tiers`` / ``_tiers_used`` — these are
    normally reached only from the online (LLM) path, which the offline
    tests above never exercise."""

    def test_select_tiers_always_includes_core(self):
        from pycsamt.agents.package_qa import TIER_CORE, _select_tiers

        ctx = _select_tiers("What is DataQCAgent?")
        self.assertIn(TIER_CORE, ctx)

    def test_select_tiers_nothing_matched_includes_everything(self):
        from pycsamt.agents.package_qa import (
            TIER_AGENTS,
            TIER_EMTOOLS,
            TIER_EXAMPLES,
            TIER_SITES,
            _select_tiers,
        )

        ctx = _select_tiers("xyzzy qwerty plugh")
        for tier in (TIER_AGENTS, TIER_SITES, TIER_EXAMPLES, TIER_EMTOOLS):
            self.assertIn(tier, ctx)

    def test_select_tiers_sites_keyword_only(self):
        from pycsamt.agents.package_qa import TIER_SITES, _select_tiers

        ctx = _select_tiers("How do I access impedance tensor data?")
        self.assertIn(TIER_SITES, ctx)

    def test_select_tiers_agent_keyword(self):
        from pycsamt.agents.package_qa import TIER_AGENTS, _select_tiers

        ctx = _select_tiers("What does the StaticShiftAgent correction do?")
        self.assertIn(TIER_AGENTS, ctx)

    def test_select_tiers_workflow_keyword(self):
        from pycsamt.agents.package_qa import TIER_CORE, _select_tiers

        ctx = _select_tiers("What workflows are supported?")
        # workflow-only match still always carries TIER_CORE
        self.assertIn(TIER_CORE, ctx)

    def test_select_tiers_emtools_keyword(self):
        from pycsamt.agents.package_qa import TIER_EMTOOLS, _select_tiers

        ctx = _select_tiers("How does estimate_ss_ama compute the shift_factor?")
        self.assertIn(TIER_EMTOOLS, ctx)


class TestTiersUsed(unittest.TestCase):
    def test_tiers_used_nothing_matched_returns_all(self):
        from pycsamt.agents.package_qa import _tiers_used

        out = _tiers_used("xyzzy qwerty plugh")
        self.assertEqual(
            out, ["core", "agents", "sites", "examples", "emtools"]
        )

    def test_tiers_used_agent_keyword(self):
        from pycsamt.agents.package_qa import _tiers_used

        out = _tiers_used("What galvanic shift correction method is used?")
        self.assertIn("core", out)
        self.assertIn("agents", out)

    def test_tiers_used_sites_keyword(self):
        from pycsamt.agents.package_qa import _tiers_used

        out = _tiers_used("How do I access impedance from Sites?")
        self.assertIn("sites", out)

    def test_tiers_used_examples_keyword(self):
        from pycsamt.agents.package_qa import _tiers_used

        out = _tiers_used("Show me a code usage example.")
        self.assertIn("examples", out)

    def test_tiers_used_emtools_keyword(self):
        from pycsamt.agents.package_qa import _tiers_used

        out = _tiers_used("Explain the sounding pseudosection collection.")
        self.assertIn("emtools", out)


# ── offline "no match" fallback ───────────────────────────────────────────────


class TestOfflineNoMatch(unittest.TestCase):
    def test_offline_answer_no_match_returns_workflow_fallback(self):
        from pycsamt.agents.package_qa import PackageQAAgent

        agent = PackageQAAgent(use_rag=False)
        r = agent.execute({"question": "xyzzy qwerty plugh"})
        self.assertEqual(r.status, "success")
        self.assertIn(
            "could not find a specific match", r.data["answer"]
        )
        self.assertEqual(r.data["excerpts"], [])


# ── RAG wiring (_build_rag, needs_clarification, rag_offline) ────────────────


class _FakeRagContext:
    def __init__(self, context_text="ctx", citations=None, empty=False):
        self.context_text = context_text
        self.citations = citations if citations is not None else [{"n": 1}]
        self._empty = empty

    def is_empty(self):
        return self._empty

    def compose_offline_answer(self):
        return "RAG says: use StaticShiftAgent for static shift correction."


class _FakeBuilder:
    def __init__(self, ctx):
        self._ctx = ctx

    def build(self, question, session=None):
        return self._ctx


class TestRagWiring(unittest.TestCase):
    """These exercise ``_build_rag``'s try/except and the RAG-dependent
    branches of ``execute`` by monkeypatching
    ``pycsamt.assistant.rag.context_builder`` — the real corpus builder
    returns ``None`` in this environment (no bundled RAG index), so the
    non-trivial branches are otherwise unreachable."""

    def setUp(self):
        import pycsamt.assistant.rag.context_builder as cb

        self.cb = cb
        self._orig_builder = cb.default_context_builder
        self._orig_clarify = cb.needs_clarification

    def tearDown(self):
        self.cb.default_context_builder = self._orig_builder
        self.cb.needs_clarification = self._orig_clarify

    def test_build_rag_returns_context_when_builder_succeeds(self):
        from pycsamt.agents.package_qa import PackageQAAgent

        fake_ctx = _FakeRagContext()
        self.cb.default_context_builder = lambda *a, **k: _FakeBuilder(fake_ctx)

        agent = PackageQAAgent()
        rag = agent._build_rag("What is the Sites class?")
        self.assertIs(rag, fake_ctx)

    def test_build_rag_swallows_exceptions(self):
        from pycsamt.agents.package_qa import PackageQAAgent

        def _boom(*a, **k):
            raise RuntimeError("corpus unavailable")

        self.cb.default_context_builder = _boom

        agent = PackageQAAgent()
        rag = agent._build_rag("What is the Sites class?")
        self.assertIsNone(rag)

    def test_execute_returns_clarifying_question_when_rag_unsure(self):
        from pycsamt.agents.package_qa import PackageQAAgent

        fake_ctx = _FakeRagContext()
        self.cb.default_context_builder = lambda *a, **k: _FakeBuilder(fake_ctx)
        self.cb.needs_clarification = lambda *a, **k: "Which agent do you mean?"

        agent = PackageQAAgent()
        r = agent.execute({"question": "what about it?"})
        self.assertEqual(r.status, "success")
        self.assertEqual(r.data["source"], "rag_clarify")
        self.assertEqual(r.data["answer"], "Which agent do you mean?")
        self.assertEqual(r.data["citations"], [])

    def test_execute_offline_uses_rag_composed_answer(self):
        from pycsamt.agents.package_qa import PackageQAAgent

        fake_ctx = _FakeRagContext(citations=[{"n": 1, "source": "x.py"}])
        self.cb.default_context_builder = lambda *a, **k: _FakeBuilder(fake_ctx)
        self.cb.needs_clarification = lambda *a, **k: None

        agent = PackageQAAgent()  # no api_key -> offline path
        r = agent.execute({"question": "What does StaticShiftAgent do?"})
        self.assertEqual(r.status, "success")
        self.assertEqual(r.data["source"], "rag_offline")
        self.assertIn("StaticShiftAgent", r.data["answer"])
        self.assertEqual(r.data["citations"], fake_ctx.citations)


# ── online (LLM) path ─────────────────────────────────────────────────────────


class TestOnlinePath(unittest.TestCase):
    """The LLM call itself is mocked out on the agent instance so these
    stay network-free and CI-safe, while still exercising the real
    ``execute`` online branch (tier selection, prompt assembly, error
    handling)."""

    def test_execute_online_success_without_rag(self):
        from pycsamt.agents.package_qa import PackageQAAgent

        agent = PackageQAAgent(api_key="fake-key", use_rag=False)
        agent.query_llm = lambda *a, **k: "The Sites class holds EDI data."
        r = agent.execute({"question": "What is the Sites class?"})
        self.assertEqual(r.status, "success")
        self.assertEqual(r.data["source"], "llm")
        self.assertEqual(
            r.data["answer"], "The Sites class holds EDI data."
        )
        self.assertIn("tiers_used", r.data)

    def test_execute_online_success_with_rag_and_extra_context(self):
        import pycsamt.assistant.rag.context_builder as cb

        from pycsamt.agents.package_qa import PackageQAAgent

        orig_builder, orig_clarify = cb.default_context_builder, cb.needs_clarification
        try:
            fake_ctx = _FakeRagContext(context_text="StaticShiftAgent docs")
            cb.default_context_builder = lambda *a, **k: _FakeBuilder(fake_ctx)
            cb.needs_clarification = lambda *a, **k: None

            agent = PackageQAAgent(api_key="fake-key")
            agent.query_llm = lambda *a, **k: "grounded answer"
            r = agent.execute(
                {
                    "question": "What does StaticShiftAgent do?",
                    "context": "session note",
                }
            )
            self.assertEqual(r.status, "success")
            self.assertEqual(r.data["source"], "llm+rag")
            self.assertEqual(r.data["citations"], fake_ctx.citations)
        finally:
            cb.default_context_builder = orig_builder
            cb.needs_clarification = orig_clarify

    def test_execute_online_llm_failure_returns_failed_status(self):
        from pycsamt.agents.package_qa import PackageQAAgent

        agent = PackageQAAgent(api_key="fake-key", use_rag=False)

        def _boom(*a, **k):
            raise RuntimeError("provider unreachable")

        agent.query_llm = _boom
        r = agent.execute({"question": "What is the Sites class?"})
        self.assertEqual(r.status, "failed")
        self.assertIn("provider unreachable", r.error)

    def test_execute_online_no_answer_falls_back_to_placeholder(self):
        from pycsamt.agents.package_qa import PackageQAAgent

        agent = PackageQAAgent(api_key="fake-key", use_rag=False)
        agent.query_llm = lambda *a, **k: ""
        r = agent.execute({"question": "What is the Sites class?"})
        self.assertEqual(r.status, "success")
        self.assertEqual(r.data["answer"], "(no answer)")


if __name__ == "__main__":
    unittest.main()


class TestQuestionSubjects(unittest.TestCase):
    """Offline composition answers each compared or named subject."""

    def test_subjects(self):
        from pycsamt.agents.package_qa import _question_subjects

        self.assertEqual(
            _question_subjects(
                "What is the difference between quality control and denoising in pyCSAMT?"),
            ["quality control", "denoising"])
        self.assertEqual(_question_subjects("Occam2D vs ModEM?"), ["Occam2D", "ModEM"])
        self.assertEqual(
            _question_subjects("Explain estimate_ss_ama and correct_ss_ama."),
            ["estimate_ss_ama", "correct_ss_ama"])
        self.assertEqual(_question_subjects("How do I load EDI files?"), [])

    def test_compose_subjects_retrieves_each(self):
        from pycsamt.agents.package_qa import PackageQAAgent

        asked = []

        class Rag:
            def __init__(self, subject):
                self.subject = subject

            def is_empty(self):
                return False

            def compose_offline_answer(self, top=3):
                return f"about {self.subject}"

        agent = PackageQAAgent()
        agent._build_rag = lambda q, session=None: asked.append(q) or Rag(q)
        text = agent._compose_subjects("Occam2D vs ModEM?")
        self.assertEqual(asked, ["Occam2D", "ModEM"])
        self.assertIn("### Occam2D\n\nabout Occam2D", text)
        self.assertIsNone(agent._compose_subjects("How do I load EDI files?"))
