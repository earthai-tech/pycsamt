# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for IoTFieldAgent."""

import unittest

import matplotlib

matplotlib.use("Agg")


def _sim_packets(n_stations=3, n_samples=256, seed=3):
    from pycsamt.iot import simulate_iot_network

    return simulate_iot_network(
        n_stations=n_stations,
        n_samples=n_samples,
        seed=seed,
    )


class TestIoTFieldAgentFromPackets(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pycsamt.agents import IoTFieldAgent

        cls.packets = _sim_packets()
        cls.agent = IoTFieldAgent()
        cls.result = cls.agent.execute(
            {
                "packets": cls.packets,
                "survey_id": "TEST_IOT",
                "manifest": True,
                "figures": True,
            }
        )

    def test_exported_and_lazy_loadable(self):
        import pycsamt.agents as A

        self.assertIn("IoTFieldAgent", A.__all__)
        self.assertIsNotNone(A.IoTFieldAgent)

    def test_status_not_failed(self):
        self.assertIn(self.result.status, ("success", "needs_review"))

    def test_session_resolved(self):
        session = self.result.get("session")
        self.assertIsNotNone(session)
        self.assertEqual(self.result.get("n_packets"), len(self.packets))
        self.assertGreater(self.result.get("n_stations"), 0)

    def test_monitoring_status_and_level(self):
        status = self.result.get("status")
        self.assertIsNotNone(status)
        self.assertIn(
            self.result.get("level"), ("ok", "warn", "critical", "unknown")
        )
        self.assertIsInstance(self.result.get("issues"), list)

    def test_tables_present(self):
        for key in (
            "status_table",
            "telemetry_summary",
            "packet_table",
            "station_table",
            "pipeline_input",
        ):
            self.assertIn(key, self.result, key)

    def test_manifest_built(self):
        self.assertIn("manifest", self.result)

    def test_figures_produced(self):
        figs = self.result.get("figures") or {}
        self.assertIn("dashboard", figs)
        self.assertGreaterEqual(len(figs), 1)

    def test_offline_no_cost(self):
        self.assertEqual(self.result.cost_estimate_usd, 0.0)
        self.assertIsNone(self.result.llm_interpretation)


class TestIoTFieldAgentEdgeCases(unittest.TestCase):
    def test_missing_source_fails_cleanly(self):
        from pycsamt.agents import IoTFieldAgent

        res = IoTFieldAgent().execute({})
        self.assertEqual(res.status, "failed")
        self.assertIsNotNone(res.error)

    def test_session_input_accepted(self):
        from pycsamt.agents import IoTFieldAgent
        from pycsamt.iot import FieldSession

        session = FieldSession("S1")
        session.add_packets(_sim_packets(n_stations=2, seed=5))
        res = IoTFieldAgent().execute({"session": session, "figures": False})
        self.assertNotEqual(res.status, "failed")
        self.assertIs(res.get("session"), session)

    def test_manifest_written_to_output_dir(self):
        import os
        import tempfile

        from pycsamt.agents import IoTFieldAgent

        outdir = tempfile.mkdtemp(prefix="iotagent_test_")
        res = IoTFieldAgent().execute(
            {
                "packets": _sim_packets(n_stations=2, seed=9),
                "survey_id": "WRITE_ME",
                "manifest": True,
                "figures": False,
                "output_dir": outdir,
            }
        )
        path = res.get("manifest_path")
        self.assertIsNotNone(path)
        self.assertTrue(os.path.isfile(path))


if __name__ == "__main__":
    unittest.main()
