# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for WorkflowOrchestratorAgent — routing, dry-run, keyword classification."""

import unittest


class TestKeywordClassification(unittest.TestCase):
    def _classify(self, text):
        from pycsamt.agents.orchestrator import (
            _keyword_classify,
        )

        wf, reason = _keyword_classify(text)
        return wf

    def test_qc_keywords(self):
        self.assertEqual(self._classify("run QC on the data"), "qc")

    def test_phase_analysis_keywords(self):
        self.assertEqual(
            self._classify("compute phase tensor and strike analysis"),
            "phase_analysis",
        )

    def test_ai_inversion_keywords(self):
        wf = self._classify("run neural network CNN inversion")
        self.assertEqual(wf, "ai_inversion")

    def test_inv3d_keywords(self):
        wf = self._classify("GCN graph convolutional 3D spatial inversion")
        self.assertEqual(wf, "inv3d")

    def test_inv2d_keywords(self):
        wf = self._classify("U-Net 2D profile inversion")
        self.assertEqual(wf, "inv2d")

    def test_ensemble_keywords(self):
        wf = self._classify("ensemble uncertainty quantification")
        self.assertEqual(wf, "ensemble_inversion")

    def test_joint_keywords(self):
        wf = self._classify("joint inversion MT+TEM multi-modal")
        self.assertEqual(wf, "joint_inversion")

    def test_modem_keywords(self):
        wf = self._classify("prepare ModEM 3D inversion")
        self.assertEqual(wf, "modem")

    def test_pre_inversion_keywords(self):
        wf = self._classify("set up Occam2D mesh and startup file")
        self.assertEqual(wf, "pre_inversion")

    def test_full_keywords(self):
        wf = self._classify("run the full pipeline")
        self.assertEqual(wf, "full")

    def test_unknown_defaults_to_qc(self):
        wf = self._classify("xyzzyx random text with no keywords")
        self.assertEqual(wf, "qc")


class TestWorkflowStepsRegistry(unittest.TestCase):
    def test_all_phase4_workflows_registered(self):
        from pycsamt.agents.orchestrator import (
            _WORKFLOW_STEPS,
        )

        for wf in ("inv2d", "inv3d", "ensemble_inversion", "joint_inversion"):
            self.assertIn(wf, _WORKFLOW_STEPS)
            self.assertGreater(len(_WORKFLOW_STEPS[wf]), 0)

    def test_all_classic_workflows_registered(self):
        from pycsamt.agents.orchestrator import (
            _WORKFLOW_STEPS,
        )

        for wf in (
            "qc",
            "phase_analysis",
            "pre_inversion",
            "ai_inversion",
            "modem",
            "full",
        ):
            self.assertIn(wf, _WORKFLOW_STEPS)

    def test_step_tuples_have_4_or_5_elements(self):
        from pycsamt.agents.orchestrator import (
            _WORKFLOW_STEPS,
        )

        for wf, steps in _WORKFLOW_STEPS.items():
            for step in steps:
                self.assertIn(
                    len(step),
                    (4, 5),
                    f"{wf}: step {step[0]!r} has wrong length {len(step)}",
                )


class TestOrchestratorDryRun(unittest.TestCase):
    def test_dry_run_qc_workflow(self):
        from pycsamt.agents import WorkflowOrchestratorAgent

        ag = WorkflowOrchestratorAgent()
        r = ag.execute(
            {
                "request": "QC the data",
                "data_path": "/data/test",
                "dry_run": True,
            }
        )
        self.assertIn(r.status, ("success", "needs_review"))

    def test_workflow_type_returned(self):
        from pycsamt.agents import WorkflowOrchestratorAgent

        ag = WorkflowOrchestratorAgent()
        r = ag.execute(
            {
                "request": "phase tensor analysis",
                "data_path": "/data/test",
                "dry_run": True,
            }
        )
        self.assertIn("workflow_type", r.data)

    def test_steps_list_in_output(self):
        from pycsamt.agents import WorkflowOrchestratorAgent

        ag = WorkflowOrchestratorAgent()
        r = ag.execute(
            {
                "config": {"workflow": "qc"},
                "data_path": "/data/test",
                "dry_run": True,
            }
        )
        steps = r.get("steps") or []
        self.assertIsInstance(steps, list)

    def test_explicit_config_workflow_respected(self):
        from pycsamt.agents import WorkflowOrchestratorAgent

        ag = WorkflowOrchestratorAgent()
        r = ag.execute(
            {
                "config": {"workflow": "modem"},
                "data_path": "/data/test",
                "dry_run": True,
            }
        )
        self.assertEqual(r.get("workflow_type"), "modem")


if __name__ == "__main__":
    unittest.main()
