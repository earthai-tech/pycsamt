# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ContextInputAgent — regex fallback, JSON extraction."""

import unittest


class TestContextInputAgentRegex(unittest.TestCase):
    """All tests run without an LLM key (regex-only mode)."""

    def setUp(self):
        from pycsamt.agents import ContextInputAgent

        self.agent = ContextInputAgent()  # no api_key

    def test_execute_returns_success(self):
        r = self.agent.execute({"request": "Load EDIs from /data/L22PLT"})
        self.assertEqual(r.status, "success")

    def test_config_key_present(self):
        r = self.agent.execute({"request": "Load /data/L22PLT, QC it"})
        self.assertIn("config", r.data)

    def test_data_path_extracted(self):
        r = self.agent.execute({"request": "Load EDIs from /data/WILLY_EDIs"})
        cfg = r.get("config") or {}
        self.assertIn("/data/WILLY_EDIs", str(cfg.get("data_path", "")))

    def test_period_range_extracted(self):
        r = self.agent.execute(
            {"request": "Load /data/EDIs, period range 0.001 to 100 s"}
        )
        cfg = r.get("config") or {}
        t_min = cfg.get("T_min") or cfg.get("period_range", [None])[0]
        self.assertIsNotNone(t_min)

    def test_workflow_type_detected_qc(self):
        r = self.agent.execute({"request": "QC the EDI files in /data/EDIs"})
        cfg = r.get("config") or {}
        wf = str(cfg.get("workflow", "")).lower()
        self.assertIn("qc", wf)

    def test_workflow_type_detected_phase_analysis(self):
        r = self.agent.execute({"request": "Run phase tensor analysis on /data/WILLY"})
        cfg = r.get("config") or {}
        wf = str(cfg.get("workflow", "")).lower()
        self.assertTrue(
            "phase" in wf or "analysis" in wf or "tensor" in wf,
            f"Expected phase/analysis/tensor in workflow, got: {wf!r}",
        )

    def test_empty_request_handled(self):
        r = self.agent.execute({"request": ""})
        # empty request may return failed — that is acceptable behaviour
        self.assertIsNotNone(r.status)

    def test_cost_zero_without_llm(self):
        r = self.agent.execute({"request": "Load /data/x"})
        self.assertEqual(r.cost_estimate_usd, 0.0)

    def test_llm_interpretation_none_without_key(self):
        r = self.agent.execute({"request": "Load /data/x"})
        self.assertIsNone(r.llm_interpretation)

    def test_elapsed_seconds_positive(self):
        r = self.agent.execute({"request": "test"})
        self.assertGreaterEqual(r.elapsed_seconds, 0.0)


class TestContextInputAgentComponents(unittest.TestCase):
    def test_extract_json_object(self):
        from pycsamt.agents import ContextInputAgent

        ag = ContextInputAgent()
        d = ag.extract_json('{"workflow": "qc", "value": 3.14}')
        self.assertEqual(d["workflow"], "qc")
        self.assertAlmostEqual(d["value"], 3.14)

    def test_extract_json_embedded_in_prose(self):
        from pycsamt.agents import ContextInputAgent

        ag = ContextInputAgent()
        obj = ag.extract_json('Here is the answer: {"x": 1} — done.')
        self.assertEqual(obj["x"], 1)

    def test_regex_extract_helper(self):
        from pycsamt.agents.context import _regex_extract

        cfg = _regex_extract("Load /data/L22 EDIs, period 0.01 to 10 s, QC")
        self.assertIn("data_path", cfg)


if __name__ == "__main__":
    unittest.main()
