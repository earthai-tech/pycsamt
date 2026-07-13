# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for _base.py — AgentResult, BaseAgent, and _pricing."""

import unittest


class TestAgentResult(unittest.TestCase):
    def _make(self, **kw):
        from pycsamt.agents._base import AgentResult

        defaults = dict(
            status="success",
            summary="ok",
            data={},
            warnings=[],
            llm_interpretation=None,
            elapsed_seconds=0.1,
            cost_estimate_usd=0.0,
            error=None,
            error_fix_hint=None,
        )
        defaults.update(kw)
        return AgentResult(**defaults)

    def test_status_success(self):
        r = self._make(status="success")
        self.assertEqual(r.status, "success")

    def test_data_access_via_getitem(self):
        r = self._make(data={"key": 42})
        self.assertEqual(r["key"], 42)

    def test_data_access_via_get(self):
        r = self._make(data={"x": 99})
        self.assertEqual(r.get("x"), 99)
        self.assertIsNone(r.get("missing"))
        self.assertEqual(r.get("missing", -1), -1)

    def test_failed_classmethod(self):
        from pycsamt.agents._base import AgentResult

        r = AgentResult.failed("something went wrong", elapsed=0.05)
        self.assertEqual(r.status, "failed")
        self.assertIsNotNone(r.error)
        self.assertIn("something went wrong", r.error)

    def test_failed_with_hint(self):
        from pycsamt.agents._base import AgentResult

        r = AgentResult.failed("bad", hint="try this fix", elapsed=0.0)
        self.assertEqual(r.error_fix_hint, "try this fix")

    def test_bool_true_for_success(self):
        r = self._make(status="success")
        self.assertTrue(bool(r))

    def test_bool_false_for_failed(self):
        from pycsamt.agents._base import AgentResult

        r = AgentResult.failed("err", elapsed=0.0)
        self.assertFalse(bool(r))

    def test_repr_contains_status(self):
        r = self._make(status="needs_review")
        self.assertIn("needs_review", repr(r))

    def test_warnings_list(self):
        r = self._make(warnings=["w1", "w2"])
        self.assertEqual(len(r.warnings), 2)


class TestPricing(unittest.TestCase):
    def test_estimate_cost_positive(self):
        from pycsamt.agents._pricing import estimate_cost

        c = estimate_cost(
            "claude",
            "claude-sonnet-4-6",
            input_tokens=1000,
            output_tokens=200,
        )
        self.assertGreater(c, 0.0)

    def test_estimate_cost_unknown_model_non_negative(self):
        from pycsamt.agents._pricing import estimate_cost

        c = estimate_cost(
            "claude",
            "claude-fake-model",
            input_tokens=1000,
            output_tokens=200,
        )
        self.assertGreaterEqual(c, 0.0)

    def test_format_cost_dollar_sign(self):
        from pycsamt.agents._pricing import format_cost

        s = format_cost(0.001234)
        self.assertIn("$", s)

    def test_get_rate_returns_dict(self):
        from pycsamt.agents._pricing import get_rate

        r = get_rate("claude", "claude-haiku-4-5")
        self.assertIn("input", r)
        self.assertIn("output", r)

    def test_provider_rates_structure(self):
        from pycsamt.agents._pricing import PROVIDER_RATES

        for provider in ("claude", "openai", "gemini"):
            self.assertIn(provider, PROVIDER_RATES)
        for models in PROVIDER_RATES.values():
            for rates in models.values():
                self.assertIn("input", rates)
                self.assertIn("output", rates)


class TestBaseAgentInstantiation(unittest.TestCase):
    def test_no_llm_repr(self):
        from pycsamt.agents import ContextInputAgent
        from pycsamt.api.agents import AGENT_CONFIG

        # offline(): ignore API keys from the developer environment
        # (.env.local) so the agent really has no LLM configured
        with AGENT_CONFIG.offline():
            ag = ContextInputAgent()
        self.assertIn("no-LLM", repr(ag))

    def test_name_attribute(self):
        from pycsamt.agents import DataQCAgent

        ag = DataQCAgent()
        self.assertEqual(ag.name, "DataQCAgent")

    def test_query_llm_without_key_returns_none(self):
        from pycsamt.agents import ContextInputAgent
        from pycsamt.api.agents import AGENT_CONFIG

        with AGENT_CONFIG.offline():
            ag = ContextInputAgent()
            result = ag.query_llm("any prompt", max_tokens=50)
        self.assertIsNone(result)

    def test_extract_json_valid(self):
        from pycsamt.agents import ContextInputAgent

        ag = ContextInputAgent()
        obj = ag.extract_json('{"a": 1, "b": "hello"}')
        self.assertEqual(obj["a"], 1)

    def test_extract_json_embedded(self):
        from pycsamt.agents import ContextInputAgent

        ag = ContextInputAgent()
        obj = ag.extract_json('Some text {"x": 99} trailing')
        self.assertEqual(obj["x"], 99)

    def test_extract_json_invalid_returns_none(self):
        from pycsamt.agents import ContextInputAgent

        ag = ContextInputAgent()
        obj = ag.extract_json("not json at all")
        self.assertIsNone(obj)


if __name__ == "__main__":
    unittest.main()
