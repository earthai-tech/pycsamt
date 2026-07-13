# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Tests for _package_context.py — the runtime API context builder.

These tests verify that the context builder:
* Produces a non-trivial, non-empty string
* Contains entries for the main pycsamt agents
* Contains entries for core data classes
* Contains workflow keywords
* Contains usage examples
* Can be rebuilt with custom parameters without error

No LLM calls — fully offline, CI-safe.
"""

from __future__ import annotations

import unittest

from pycsamt.agents._package_context import (
    PACKAGE_CONTEXT,
    _collect_agents,
    _collect_core,
    _short_doc,
    _sig,
    build_package_context,
)
from pycsamt.agents.package_qa import _keywords


class TestContextBuilder(unittest.TestCase):
    """PACKAGE_CONTEXT singleton."""

    def test_non_empty(self):
        self.assertGreater(len(PACKAGE_CONTEXT), 500)

    def test_contains_pycsamt_header(self):
        self.assertIn("pycsamt", PACKAGE_CONTEXT)

    def test_contains_workflow_keywords(self):
        for kw in (
            "qc",
            "static_shift",
            "phase_analysis",
            "ai_inversion",
            "pinn_inversion",
            "report",
        ):
            self.assertIn(
                kw,
                PACKAGE_CONTEXT,
                msg=f"Workflow keyword missing: {kw!r}",
            )

    def test_contains_agent_class_names(self):
        for name in (
            "StaticShiftAgent",
            "DataQCAgent",
            "WorkflowOrchestratorAgent",
            "MTLoaderAgent",
        ):
            self.assertIn(
                name,
                PACKAGE_CONTEXT,
                msg=f"Agent class missing: {name!r}",
            )

    def test_contains_usage_examples(self):
        self.assertIn("Usage examples", PACKAGE_CONTEXT)
        self.assertIn("from pycsamt", PACKAGE_CONTEXT)

    def test_rebuild_with_custom_limits(self):
        ctx = build_package_context(
            max_agent_doc=100,
            max_core_doc=200,
        )
        self.assertGreater(len(ctx), 200)

    def test_rebuild_idempotent(self):
        ctx1 = build_package_context()
        ctx2 = build_package_context()
        self.assertEqual(ctx1, ctx2)


class TestCollectors(unittest.TestCase):
    """_collect_agents and _collect_core helpers."""

    def test_collect_agents_returns_dict(self):
        agents = _collect_agents()
        self.assertIsInstance(agents, dict)

    def test_collect_agents_nonempty(self):
        agents = _collect_agents()
        self.assertGreater(
            len(agents),
            0,
            "No agents collected — import failed?",
        )

    def test_collect_agents_static_shift(self):
        agents = _collect_agents()
        self.assertIn(
            "StaticShiftAgent",
            agents,
            "StaticShiftAgent not collected",
        )

    def test_collect_agents_qc(self):
        agents = _collect_agents()
        self.assertIn(
            "DataQCAgent",
            agents,
            "DataQCAgent not collected",
        )

    def test_collect_core_returns_dict(self):
        core = _collect_core()
        self.assertIsInstance(core, dict)

    def test_collect_core_sites(self):
        core = _collect_core()
        self.assertIn(
            "Sites",
            core,
            "Sites class not collected",
        )


class TestHelpers(unittest.TestCase):
    """_keywords, _short_doc, _sig helpers."""

    def test_keywords_filters_stopwords(self):
        kw = _keywords("what is the best way")
        # 'what','is','the' are stopwords
        self.assertNotIn("what", kw)
        self.assertNotIn("the", kw)

    def test_keywords_extracts_terms(self):
        kw = _keywords("StaticShiftAgent corrects galvanic shift")
        self.assertIn("staticshiftagent", kw)
        self.assertIn("corrects", kw)
        self.assertIn("galvanic", kw)

    def test_short_doc_truncates(self):
        from pycsamt.agents.static_shift import (
            StaticShiftAgent,
        )

        doc = _short_doc(StaticShiftAgent, 50)
        self.assertLessEqual(len(doc), 60)

    def test_short_doc_nonempty(self):
        from pycsamt.agents.static_shift import (
            StaticShiftAgent,
        )

        doc = _short_doc(StaticShiftAgent)
        self.assertGreater(len(doc), 5)

    def test_sig_returns_string(self):
        from pycsamt.agents.static_shift import (
            StaticShiftAgent,
        )

        sig = _sig(StaticShiftAgent)
        self.assertIsInstance(sig, str)

    def test_sig_excludes_self(self):
        from pycsamt.agents.static_shift import (
            StaticShiftAgent,
        )

        sig = _sig(StaticShiftAgent)
        self.assertNotIn("self", sig)


if __name__ == "__main__":
    unittest.main()
