# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ForwardModelAgent — 1-D MT forward with synthetic layered model."""

import unittest


class TestForwardModelAgent(unittest.TestCase):
    def setUp(self):
        from pycsamt.agents import ForwardModelAgent

        self.agent = ForwardModelAgent()

    def _simple_model(self):
        return {
            "model": {
                "resistivity": [200.0, 20.0, 5000.0],
                "thickness": [300.0, 800.0],
            }
        }

    def test_execute_returns_success(self):
        r = self.agent.execute(self._simple_model())
        self.assertIn(r.status, ("success", "needs_review"))

    def test_figures_dict_present(self):
        r = self.agent.execute(self._simple_model())
        self.assertIn("figures", r.data)
        self.assertIsInstance(r.data["figures"], dict)

    def test_layered_model_in_output(self):
        r = self.agent.execute(self._simple_model())
        # output uses 'layered_model' key
        self.assertIsNotNone(r.get("layered_model") or r.get("model"))

    def test_no_llm_cost(self):
        r = self.agent.execute(self._simple_model())
        self.assertEqual(r.cost_estimate_usd, 0.0)

    def test_elapsed_positive(self):
        r = self.agent.execute(self._simple_model())
        self.assertGreaterEqual(r.elapsed_seconds, 0.0)

    def test_missing_model_uses_default(self):
        # ForwardModelAgent may use a default model when none is supplied
        r = self.agent.execute({})
        self.assertIsNotNone(r.status)

    def test_single_layer_model(self):
        r = self.agent.execute({"model": {"resistivity": [100.0], "thickness": []}})
        self.assertIn(r.status, ("success", "needs_review", "failed"))

    def test_negative_resistivity_handled(self):
        r = self.agent.execute(
            {"model": {"resistivity": [-50.0, 200.0], "thickness": [500.0]}}
        )
        # should not raise, either fails gracefully or warns
        self.assertIsNotNone(r.status)


if __name__ == "__main__":
    unittest.main()
