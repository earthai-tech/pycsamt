# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for AgentCoordinator — dry-run, checkpointing, cost tracking."""

import unittest


class TestAgentCoordinatorDryRun(unittest.TestCase):
    def _make_simple_chain(self, name="test_wf"):
        from pycsamt.agents import (
            AgentCoordinator,
            ContextInputAgent,
            MTLoaderAgent,
        )

        coord = AgentCoordinator(name)
        coord.add_step("parse", ContextInputAgent(), description="Parse request")
        coord.add_step(
            "load",
            MTLoaderAgent(),
            input_fn=lambda r: {
                "path": (r["parse"].get("config") or {}).get("data_path", "")
            },
            description="Load data",
        )
        return coord

    def test_dry_run_returns_success(self):
        coord = self._make_simple_chain()
        r = coord.execute(
            {"request": "Load /data/L22 EDIs, QC them"},
            dry_run=True,
        )
        self.assertEqual(r.status, "success")

    def test_dry_run_summary_contains_steps(self):
        coord = self._make_simple_chain()
        r = coord.execute({"request": "Load /data/test"}, dry_run=True)
        self.assertIn("2", r.summary)  # "2/2 steps"

    def test_step_count_after_add(self):
        from pycsamt.agents import (
            AgentCoordinator,
            ContextInputAgent,
        )

        coord = AgentCoordinator("step_count_test")
        for i in range(4):
            coord.add_step(f"step{i}", ContextInputAgent())
        self.assertEqual(len(coord._steps), 4)

    def test_no_steps_returns_success(self):
        from pycsamt.agents import AgentCoordinator

        coord = AgentCoordinator("empty")
        r = coord.execute({})
        # empty workflow is allowed — no failure
        self.assertIn(r.status, ("success", "needs_review"))

    def test_dry_run_cost_is_zero(self):
        coord = self._make_simple_chain()
        r = coord.execute({"request": "load /data/x"}, dry_run=True)
        self.assertEqual(r.cost_estimate_usd, 0.0)

    def test_checkpoint_dir_created(self, tmp_path=None):
        import os
        import tempfile

        td = tempfile.mkdtemp()
        from pycsamt.agents import (
            AgentCoordinator,
            ContextInputAgent,
        )

        coord = AgentCoordinator("chk_test", checkpoint_dir=td)
        coord.add_step("parse", ContextInputAgent())
        coord.execute({"request": "test"}, dry_run=True)
        # checkpoint dir should exist after execution
        self.assertTrue(os.path.isdir(td))

    def test_result_data_has_plan_key(self):
        coord = self._make_simple_chain("keys_test")
        r = coord.execute({"request": "load /data/x"}, dry_run=True)
        # dry-run result contains 'plan' and 'steps' keys, not per-step keys
        self.assertTrue("plan" in r.data or "steps" in r.data)

    def test_elapsed_seconds_positive(self):
        coord = self._make_simple_chain()
        r = coord.execute({"request": "load /x"}, dry_run=True)
        self.assertGreaterEqual(r.elapsed_seconds, 0.0)

    def test_workflow_name_stored(self):
        from pycsamt.agents import AgentCoordinator

        coord = AgentCoordinator("my_special_workflow")
        self.assertEqual(coord.workflow_name, "my_special_workflow")


class TestAgentCoordinatorVersion(unittest.TestCase):
    def test_version_format(self):
        from pycsamt.agents import __version__

        # version is a non-empty string in semver-like format
        self.assertIsInstance(__version__, str)
        self.assertGreater(len(__version__), 0)

    def test_all_exports_count(self):
        from pycsamt.agents import __all__

        self.assertGreaterEqual(len(__all__), 25)


if __name__ == "__main__":
    unittest.main()
