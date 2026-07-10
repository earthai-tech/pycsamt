# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for MTLoaderAgent — load EDI files, QC table, station count."""

import unittest
from pathlib import Path

_EDI_DIR = Path(__file__).parents[4] / "data" / "3edis"
_HAS_EDIS = _EDI_DIR.exists() and any(_EDI_DIR.glob("*.edi"))


@unittest.skipUnless(_HAS_EDIS, "3edis dataset not found")
class TestMTLoaderAgentWithRealData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from pycsamt.agents import MTLoaderAgent
        cls.agent = MTLoaderAgent()
        cls.result = cls.agent.execute({"path": str(_EDI_DIR)})

    def test_status_success(self):
        self.assertEqual(self.result.status, "success")

    def test_n_stations_matches_files(self):
        n_files = len(list(_EDI_DIR.glob("*.edi")))
        self.assertEqual(self.result["n_stations"], n_files)

    def test_sites_object_present(self):
        sites = self.result["sites"]
        self.assertIsNotNone(sites)

    def test_qc_table_present(self):
        qt = self.result.get("qc_table")
        self.assertIsNotNone(qt)

    def test_station_names_non_empty(self):
        names = self.result.get("station_names") or []
        self.assertGreater(len(names), 0)

    def test_cost_zero_without_llm(self):
        self.assertEqual(self.result.cost_estimate_usd, 0.0)

    def test_elapsed_positive(self):
        self.assertGreater(self.result.elapsed_seconds, 0.0)

    def test_no_error(self):
        self.assertIsNone(self.result.error)


class TestMTLoaderAgentMissingPath(unittest.TestCase):

    def test_nonexistent_path_fails(self):
        from pycsamt.agents import MTLoaderAgent
        r = MTLoaderAgent().execute({"path": "/nonexistent/path/to/edis"})
        self.assertEqual(r.status, "failed")
        self.assertIsNotNone(r.error)

    def test_no_path_fails(self):
        from pycsamt.agents import MTLoaderAgent
        r = MTLoaderAgent().execute({})
        self.assertEqual(r.status, "failed")


if __name__ == "__main__":
    unittest.main()
