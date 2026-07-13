# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for DataQCAgent and StaticShiftAgent."""

import unittest
from pathlib import Path

_EDI_DIR = Path(__file__).parents[4] / "data" / "3edis"
_HAS_EDIS = _EDI_DIR.exists() and any(_EDI_DIR.glob("*.edi"))


def _load_sites():
    from pycsamt.agents import MTLoaderAgent

    r = MTLoaderAgent().execute({"path": str(_EDI_DIR)})
    if r.status != "success":
        return None
    return r["sites"]


@unittest.skipUnless(_HAS_EDIS, "3edis dataset not found")
class TestDataQCAgent(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pycsamt.agents import DataQCAgent

        cls.sites = _load_sites()
        if cls.sites is None:
            raise unittest.SkipTest("Could not load sites.")
        cls.agent = DataQCAgent()
        cls.result = cls.agent.execute({"sites": cls.sites})

    def test_status_success_or_needs_review(self):
        self.assertIn(self.result.status, ("success", "needs_review"))

    def test_n_stations_present(self):
        n = self.result.get("n_stations")
        self.assertIsNotNone(n)
        self.assertGreater(n, 0)

    def test_qc_scores_dict(self):
        scores = self.result.get("qc_scores") or {}
        self.assertIsInstance(scores, dict)

    def test_n_flagged_non_negative(self):
        n = self.result.get("n_flagged", 0)
        self.assertGreaterEqual(n, 0)

    def test_sites_in_output(self):
        self.assertIsNotNone(self.result.get("sites"))

    def test_no_llm_cost(self):
        self.assertEqual(self.result.cost_estimate_usd, 0.0)


@unittest.skipUnless(_HAS_EDIS, "3edis dataset not found")
class TestStaticShiftAgent(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pycsamt.agents import StaticShiftAgent

        cls.sites = _load_sites()
        if cls.sites is None:
            raise unittest.SkipTest("Could not load sites.")
        cls.agent = StaticShiftAgent()
        cls.result = cls.agent.execute({"sites": cls.sites})

    def test_status_not_failed(self):
        self.assertNotEqual(self.result.status, "failed")

    def test_corrected_sites_present(self):
        corrected = self.result.get("corrected_sites")
        self.assertIsNotNone(corrected)

    def test_delta_stats_present(self):
        ds = self.result.get("delta_stats")
        self.assertIsNotNone(ds)

    def test_cost_zero(self):
        self.assertEqual(self.result.cost_estimate_usd, 0.0)

    def test_elapsed_positive(self):
        self.assertGreater(self.result.elapsed_seconds, 0.0)


if __name__ == "__main__":
    unittest.main()
