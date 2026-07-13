# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ModelZooAgent — list action, metadata, error handling."""

import unittest


class TestModelZooAgentList(unittest.TestCase):
    def setUp(self):
        from pycsamt.agents import ModelZooAgent

        self.agent = ModelZooAgent()

    def test_list_status_success(self):
        r = self.agent.execute({"action": "list"})
        self.assertEqual(r.status, "success")

    def test_list_returns_models_dict(self):
        r = self.agent.execute({"action": "list"})
        models = r.get("models")
        self.assertIsInstance(models, dict)
        self.assertGreater(len(models), 0)

    def test_known_model_in_list(self):
        r = self.agent.execute({"action": "list"})
        self.assertIn("mt1d-resnet-5layer-v1", r["models"])

    def test_details_list_populated(self):
        r = self.agent.execute({"action": "list"})
        details = r.get("details") or []
        self.assertGreater(len(details), 0)
        for row in details:
            self.assertIn("name", row)
            self.assertIn("arch", row)
            self.assertIn("solver", row)

    def test_cost_zero(self):
        r = self.agent.execute({"action": "list"})
        self.assertEqual(r.cost_estimate_usd, 0.0)

    def test_no_llm_interpretation(self):
        r = self.agent.execute({"action": "list"})
        self.assertIsNone(r.llm_interpretation)


class TestModelZooAgentDownload(unittest.TestCase):
    def test_download_unknown_model_fails(self):
        from pycsamt.agents import ModelZooAgent

        ag = ModelZooAgent()
        r = ag.execute(
            {"action": "download", "model_name": "nonexistent-model-v99"}
        )
        self.assertEqual(r.status, "failed")

    def test_download_without_model_name_fails(self):
        from pycsamt.agents import ModelZooAgent

        ag = ModelZooAgent()
        r = ag.execute({"action": "download"})
        self.assertEqual(r.status, "failed")

    def test_download_known_model_graceful(self):
        # Weights may not be released yet — should succeed or needs_review, not crash
        from pycsamt.agents import ModelZooAgent

        ag = ModelZooAgent()
        r = ag.execute(
            {"action": "download", "model_name": "mt1d-resnet-5layer-v1"}
        )
        self.assertIn(r.status, ("success", "needs_review"))
        self.assertIn("model_info", r.data)
        self.assertEqual(r.data["model_info"]["arch"], "resnet")


class TestModelZooAgentPredictNoSites(unittest.TestCase):
    def test_predict_without_sites_fails(self):
        from pycsamt.agents import ModelZooAgent

        ag = ModelZooAgent()
        r = ag.execute(
            {
                "action": "predict",
                "model_name": "mt1d-resnet-5layer-v1",
                # no path / sites
            }
        )
        self.assertEqual(r.status, "failed")


class TestModelZooRegistry(unittest.TestCase):
    def test_list_pretrained_function(self):
        from pycsamt.ai._zoo import list_pretrained

        models = list_pretrained()
        self.assertIsInstance(models, dict)
        for name, desc in models.items():
            self.assertIsInstance(name, str)
            self.assertIsInstance(desc, str)

    def test_get_pretrained_info_known(self):
        from pycsamt.ai._zoo import get_pretrained_info

        info = get_pretrained_info("mt1d-cnn-5layer-v1")
        self.assertEqual(info["arch"], "cnn1d")
        self.assertEqual(info["n_layers"], 5)

    def test_get_pretrained_info_unknown_raises(self):
        from pycsamt.ai._zoo import get_pretrained_info

        with self.assertRaises(KeyError):
            get_pretrained_info("bogus-model-v0")

    def test_unknown_action_fails(self):
        from pycsamt.agents import ModelZooAgent

        ag = ModelZooAgent()
        r = ag.execute(
            {"action": "teleport", "model_name": "mt1d-resnet-5layer-v1"}
        )
        self.assertEqual(r.status, "failed")


if __name__ == "__main__":
    unittest.main()
