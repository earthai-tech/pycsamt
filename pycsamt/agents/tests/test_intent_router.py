# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Intent-routing benchmark for :class:`pycsamt.agents.router.IntentRouter`.

The router is the top-level dispatch that decides WHAT KIND of request a
message is (question / code / plot / workflow / meta) before assuming it is a
workflow to execute.  This suite measures the offline (no-LLM) classifier,
which is also the cheap pre-check the GUI uses to decide whether a message
needs an EDI dataset loaded.

Run::

    pytest pycsamt/agents/tests/test_intent_router.py -v
"""

from __future__ import annotations

import unittest

from pycsamt.agents.router import (
    CODE,
    META,
    NO_DATA_INTENTS,
    PLOT,
    QUESTION,
    WORKFLOW,
    IntentRouter,
    RouterDecision,
    classify_intent_offline,
)

# (request, expected_intent)
_CASES: list[tuple[str, str]] = [
    # ── questions (want an explanation, not execution) ──────────────────────
    ("what does StaticShiftAgent do?", QUESTION),
    ("what does the static shift agent do?", QUESTION),
    ("How do I access impedance Z values?", QUESTION),
    ("explain the difference between qc and denoise", QUESTION),
    ("which method should I use for static shift?", QUESTION),
    ("tell me about the Sites data model", QUESTION),
    ("what is a phase tensor?", QUESTION),
    # ── code generation ─────────────────────────────────────────────────────
    ("write a script to load EDIs", CODE),
    ("generate code for AI inversion", CODE),
    ("give me a python script for qc", CODE),
    ("create a notebook for phase tensor analysis", CODE),
    ("show me code to correct static shift", CODE),
    # ── meta / capability / greeting ────────────────────────────────────────
    ("hi", META),
    ("hello", META),
    ("what can you do?", META),
    ("who are you", META),
    ("how do you work", META),
    # capability / catalogue listing must be META, not a QC workflow
    ("get the list of agents", META),
    ("list agents", META),
    ("list the tasks you can perform", META),
    ("which workflows are available", META),
    ("what workflows are supported?", META),
    ("what tasks can you do", META),
    ("show me what you can do", META),
    ("available workflows", META),
    ("what agents do you have", META),
    # ── workflow execution (imperative on data) ─────────────────────────────
    ("run AI inversion on /data/willy", WORKFLOW),
    ("run QC on /data/L22PLT", WORKFLOW),
    ("static shift correction for /data/willy", WORKFLOW),
    ("perform phase tensor analysis on /data/x", WORKFLOW),
    ("denoise /surveys/line1", WORKFLOW),
    ("can you run qc on /data/x?", WORKFLOW),  # polite imperative + path
    ("invert /data/profile1 with the CNN", WORKFLOW),
    # ── plot ────────────────────────────────────────────────────────────────
    ("plot the pseudosection", PLOT),
    ("visualize the resistivity model", PLOT),
]

_MIN_ACCURACY = 0.85


class TestIntentRouterOffline(unittest.TestCase):
    """Offline intent classification accuracy + behaviour."""

    @classmethod
    def setUpClass(cls):
        cls.results = [
            {
                "text": txt,
                "expected": exp,
                "got": classify_intent_offline(txt)[0],
            }
            for txt, exp in _CASES
        ]

    def _check_group(self, expected: str) -> None:
        cases = [r for r in self.results if r["expected"] == expected]
        wrong = [c for c in cases if c["got"] != expected]
        self.assertEqual(
            wrong,
            [],
            msg=(
                f"Intent {expected!r}: {len(wrong)}/{len(cases)} wrong\n"
                + "\n".join(
                    f"  got [{c['got']}] for: {c['text']!r}" for c in wrong
                )
            ),
        )

    def test_question_group(self):
        self._check_group(QUESTION)

    def test_code_group(self):
        self._check_group(CODE)

    def test_meta_group(self):
        self._check_group(META)

    def test_workflow_group(self):
        self._check_group(WORKFLOW)

    def test_plot_group(self):
        self._check_group(PLOT)

    def test_overall_accuracy(self):
        n_ok = sum(1 for r in self.results if r["got"] == r["expected"])
        acc = n_ok / len(self.results)
        wrong = [r for r in self.results if r["got"] != r["expected"]]
        self.assertGreaterEqual(
            acc,
            _MIN_ACCURACY,
            msg=(
                f"Intent accuracy {acc:.1%} < {_MIN_ACCURACY:.0%}\n"
                + "\n".join(
                    f"  [{r['got']} != {r['expected']}] {r['text']!r}"
                    for r in wrong
                )
            ),
        )


class TestRouterDecision(unittest.TestCase):
    """RouterDecision semantics and the offline route() path."""

    def test_needs_data_flag(self):
        for intent in (QUESTION, CODE, META):
            d = RouterDecision(intent=intent)
            self.assertFalse(
                d.needs_data,
                msg=f"{intent} should not need data",
            )
        for intent in (WORKFLOW, PLOT):
            d = RouterDecision(intent=intent)
            self.assertTrue(
                d.needs_data,
                msg=f"{intent} should need data",
            )

    def test_no_data_intents_set(self):
        self.assertIn(QUESTION, NO_DATA_INTENTS)
        self.assertIn(CODE, NO_DATA_INTENTS)
        self.assertIn(META, NO_DATA_INTENTS)
        self.assertNotIn(WORKFLOW, NO_DATA_INTENTS)
        self.assertNotIn(PLOT, NO_DATA_INTENTS)

    def test_offline_route_returns_decision(self):
        router = IntentRouter()  # no key → offline
        d = router.route("what does StaticShiftAgent do?")
        self.assertIsInstance(d, RouterDecision)
        self.assertEqual(d.intent, QUESTION)
        self.assertEqual(d.source, "offline")
        self.assertFalse(d.needs_data)

    def test_offline_route_workflow(self):
        router = IntentRouter()
        d = router.route("run QC on /data/willy")
        self.assertEqual(d.intent, WORKFLOW)
        self.assertTrue(d.needs_data)

    def test_empty_message_is_meta(self):
        router = IntentRouter()
        d = router.route("   ")
        self.assertEqual(d.intent, META)

    def test_execute_contract(self):
        """execute() wraps route() in an AgentResult."""
        router = IntentRouter()
        res = router.execute({"request": "write a script for qc"})
        self.assertEqual(res.status, "success")
        self.assertIn("decision", res.data)
        self.assertEqual(res["decision"].intent, CODE)


class TestLazyExports(unittest.TestCase):
    """Router and QA agents are reachable from the package root."""

    def test_intent_router_exported(self):
        import pycsamt.agents as A

        self.assertTrue(hasattr(A, "IntentRouter"))

    def test_package_qa_exported(self):
        import pycsamt.agents as A

        self.assertTrue(hasattr(A, "PackageQAAgent"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
