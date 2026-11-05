# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regression checks for the CAGEO audit's execution-contract failures.

Uses dummy agents, no model or solver, and ordinary test-owned directories.
Can run via unittest without pytest's platform-specific temporary fixtures.
"""
from __future__ import annotations

import json
from pathlib import Path
import unittest
from unittest.mock import patch
from uuid import uuid4

from pycsamt.agents._base import AgentResult
from pycsamt.agents._workflow_plan import WorkflowPlan
from pycsamt.agents.coordinator import AgentCoordinator
from pycsamt.agents._checkpoint_identity import workflow_identity


class DummyAgent:
    def __init__(self, status="success", scale=1):
        self.status = status
        self.scale = scale
        self.calls = []

    def execute(self, config):
        self.calls.append(dict(config))
        return AgentResult(status=self.status, summary="dummy result",
                           data={"input": dict(config), "call": len(self.calls), "scale": self.scale},
                           error="dummy failure" if self.status == "failed" else None)


class ExecutionSafeguards(unittest.TestCase):
    def setUp(self):
        # Ordinary mkdir avoids managed-Windows private-temp-directory ACLs.
        root = Path(__file__).resolve().parents[3]
        self.directory = root / ".tmp" / "execution_safeguards" / uuid4().hex
        self.directory.mkdir(parents=True)
        self.input = self.directory / "input.edi"
        self.input.write_text("original input", encoding="utf-8")

    def coordinator(self, name="test", agent=None):
        agent = agent or DummyAgent()
        coord = AgentCoordinator(name, checkpoint_dir=self.directory / name, verbose=False)
        coord.add_step("first", agent)
        return coord, agent

    def test_rejects_invalid_scientific_parameters(self):
        cases = [
            {"period_range": [-1, -2]}, {"period_range": [2, 1]},
            {"period_range": [1, 1]}, {"period_range": [True, 2]},
            {"period_range": [1, float("nan")]}, {"period_range": "1,2"},
            {"period_range": [1]}, {"period_range": ["bad", 2]},
            {"component": "invalid"}, {"component": ["xy"]},
            {"inversion_code": "unregistered"}, {"depth_max_km": 0},
            {"depth_max_km": float("inf")}, {"n_periods": 2.5},
            {"n_periods": True}, {"step_params": {"load": {"period_range": [5, 1]}}},
        ]
        for parameters in cases:
            with self.subTest(parameters=parameters):
                plan = WorkflowPlan("QC", "qc", data_path=str(self.input), parameters=parameters)
                self.assertFalse(plan.is_valid())
                self.assertTrue(plan.validation_errors())

    def test_valid_core_parameters(self):
        plan = WorkflowPlan("QC", "qc", data_path=str(self.input), parameters={
            "period_range": [0.001, 1.0], "component": "xy",
            "inversion_code": "occam2d", "depth_max_km": 2, "n_periods": 20})
        self.assertTrue(plan.is_valid())

    def test_normalization_preserves_invalid_intent(self):
        from pycsamt.agents.context import _normalise_config, _validate_config

        config = _normalise_config({"workflow": "qc", "data_path": str(self.input),
                                    "period_range": [10, 1], "component": "invalid"}, "QC")
        self.assertEqual(config["period_range"], [10, 1])
        self.assertEqual(config["component"], "invalid")
        self.assertTrue(_validate_config(config))

    def test_plan_preserves_review_flags_and_all_settings(self):
        config = {"workflow": "qc", "data_path": str(self.input), "requires_human_review": True,
                  "risk_flags": ["expert decision required"], "step_params": {"qc": {"threshold": 3}}}
        plan = WorkflowPlan.from_config(config, request="QC")
        self.assertTrue(plan.requires_human_review)
        self.assertIn("expert decision required", plan.risk_flags)
        self.assertEqual(plan.parameters["step_params"], config["step_params"])

    def test_required_review_stops_downstream(self):
        coord, first = self.coordinator(agent=DummyAgent("needs_review"))
        downstream = DummyAgent()
        coord.add_step("downstream", downstream)
        result = coord.execute({"path": str(self.input)})
        self.assertEqual(result.status, "needs_review")
        self.assertEqual(len(first.calls), 1)
        self.assertEqual(downstream.calls, [])
        self.assertEqual(set(result.data), {"first"})

    def test_optional_review_also_stops_downstream(self):
        coord = AgentCoordinator("optional_review", checkpoint_dir=self.directory / "optional", verbose=False)
        coord.add_step("review", DummyAgent("needs_review"), required=False)
        downstream = DummyAgent()
        coord.add_step("downstream", downstream)
        self.assertEqual(coord.execute({}).status, "needs_review")
        self.assertEqual(downstream.calls, [])

    def test_required_failure_and_optional_failure_contracts(self):
        for required in (True, False):
            with self.subTest(required=required):
                coord = AgentCoordinator(str(required), checkpoint_dir=self.directory / str(required), verbose=False)
                coord.add_step("failure", DummyAgent("failed"), required=required)
                downstream = DummyAgent()
                coord.add_step("downstream", downstream)
                result = coord.execute({})
                self.assertEqual(result.status, "failed" if required else "needs_review")
                self.assertEqual(len(downstream.calls), 0 if required else 1)

    def test_unchanged_success_can_resume(self):
        coord, agent = self.coordinator()
        config = {"path": str(self.input)}
        self.assertIsNotNone(workflow_identity(config, coord._steps, coord.workflow_name))
        coord.execute(config)
        result = coord.execute(config, resume=True)
        self.assertEqual(len(agent.calls), 1)
        self.assertEqual(result.status, "success")
        metadata = json.loads((coord._ckpt_dir / "first.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["schema_version"], 1)
        self.assertEqual(len(metadata["payload_sha256"]), 64)

    def test_changed_config_recomputes(self):
        coord, agent = self.coordinator()
        coord.execute({"dataset": "A"})
        result = coord.execute({"dataset": "B"}, resume=True)
        self.assertEqual(len(agent.calls), 2)
        self.assertEqual(result.data["first"].data["input"], {"dataset": "B"})

    def test_changed_input_bytes_at_same_path_recomputes(self):
        coord, agent = self.coordinator()
        config = {"path": str(self.input)}
        coord.execute(config)
        self.input.write_text("modified input", encoding="utf-8")
        coord.execute(config, resume=True)
        self.assertEqual(len(agent.calls), 2)

    def test_changed_agent_settings_recompute(self):
        coord, agent = self.coordinator()
        coord.execute({})
        agent.scale = 2
        result = coord.execute({}, resume=True)
        self.assertEqual(len(agent.calls), 2)
        self.assertEqual(result.data["first"].data["scale"], 2)

    def test_changed_adapter_closure_recomputes(self):
        coord, first = self.coordinator()
        settings = {"scale": 1}
        downstream = DummyAgent()
        coord.add_step("second", downstream, input_fn=lambda results: dict(settings))
        coord.execute({})
        settings["scale"] = 2
        coord.execute({}, resume=True)
        self.assertEqual(len(downstream.calls), 2)
        self.assertEqual(downstream.calls[-1], {"scale": 2})

    def test_changed_package_source_identity_recomputes(self):
        from pycsamt.agents import _checkpoint_identity as identity

        coord, agent = self.coordinator()
        coord.execute({})
        original = identity.file_hash
        with patch.object(identity, "file_hash", side_effect=lambda p: "changed" if p.name == "_workflow_plan.py" else original(p)):
            coord.execute({}, resume=True)
        self.assertEqual(len(agent.calls), 2)

    def test_changed_upstream_result_invalidates_downstream(self):
        coord, first = self.coordinator()
        downstream = DummyAgent()
        coord.add_step("second", downstream)
        coord.execute({})
        (coord._ckpt_dir / "first.json").unlink()
        coord.execute({}, resume=True)
        self.assertEqual(len(first.calls), 2)
        self.assertEqual(len(downstream.calls), 2)

    def test_failed_checkpoint_is_reexecuted(self):
        coord, agent = self.coordinator(agent=DummyAgent("failed"))
        self.assertEqual(coord.execute({}).status, "failed")
        self.assertEqual(coord.execute({}, resume=True).status, "failed")
        self.assertEqual(len(agent.calls), 2)

    def test_review_checkpoint_is_reexecuted(self):
        coord, agent = self.coordinator(agent=DummyAgent("needs_review"))
        coord.execute({})
        result = coord.execute({}, resume=True)
        self.assertEqual(result.status, "needs_review")
        self.assertEqual(len(agent.calls), 2)

    def test_corrupt_payload_is_not_deserialized(self):
        coord, agent = self.coordinator()
        coord.execute({})
        (coord._ckpt_dir / "first.pkl").write_bytes(b"corrupted")
        with patch("pycsamt.agents.coordinator.pickle.load") as unpickle:
            result = coord.execute({}, resume=True)
            unpickle.assert_not_called()
        self.assertEqual(result.status, "success")
        self.assertEqual(len(agent.calls), 2)

    def test_opaque_configuration_disables_reuse(self):
        coord, agent = self.coordinator()
        config = {"sites": object()}
        self.assertIsNone(workflow_identity(config, coord._steps, coord.workflow_name))
        coord.execute(config)
        coord.execute(config, resume=True)
        self.assertEqual(len(agent.calls), 2)

    def test_approval_binds_parameters_and_input_bytes(self):
        plan = WorkflowPlan("QC", "qc", data_path=str(self.input), requires_human_review=True)
        first = plan.approval_digest()
        plan.parameters["period_range"] = [0.001, 1]
        self.assertNotEqual(first, plan.approval_digest())
        second = plan.approval_digest()
        self.input.write_text("new input", encoding="utf-8")
        self.assertNotEqual(second, plan.approval_digest())

    def test_nested_input_identity_invalidates_approval(self):
        other = self.directory / "other.npz"
        other.write_bytes(b"weights")
        plan = WorkflowPlan("QC", "qc", data_path=str(self.input),
                            parameters={"step_params": {"invert": {"checkpoint": str(other)}}})
        first = plan.approval_digest()
        other.write_bytes(b"changed weights")
        self.assertNotEqual(first, plan.approval_digest())

    def test_orchestrator_requires_exact_approval_and_blocks_invalid_parameters(self):
        import pycsamt.agents.orchestrator as module
        from pycsamt.api.agents import AGENT_CONFIG

        with AGENT_CONFIG.offline():
            orchestrator = module.WorkflowOrchestratorAgent()
        dummy = DummyAgent()
        payload = {"request": "QC", "data_path": str(self.input),
                   "output_dir": str(self.directory / "run"),
                   "config": {"workflow": "qc", "requires_human_review": True}}
        with patch.object(module, "_WORKFLOW_STEPS", {"qc": [("first", "Dummy", None, "dummy")]}), \
             patch.object(module, "_build_registry", return_value=({"Dummy": dummy}, {})), \
             patch.object(orchestrator, "_retrieve_context", return_value=None), \
             patch.object(module, "_write_provenance"):
            pending = orchestrator.execute(payload)
            self.assertEqual(pending.status, "needs_review")
            self.assertEqual(dummy.calls, [])
            self.assertFalse(Path(payload["output_dir"]).exists())
            payload["approved_plan_hash"] = pending.data["approval_plan_hash"]
            completed = orchestrator.execute(payload)
            self.assertEqual(completed.status, "success")
            self.assertEqual(len(dummy.calls), 1)
            payload["config"]["period_range"] = [0.01, 1]
            self.assertEqual(orchestrator.execute(payload).status, "needs_review")
            self.assertEqual(len(dummy.calls), 1)
            payload["config"]["period_range"] = [-1, 1]
            self.assertEqual(orchestrator.execute(payload).status, "failed")
            self.assertEqual(len(dummy.calls), 1)


if __name__ == "__main__":
    unittest.main()


class ConventionalRoutingSafeguards(unittest.TestCase):
    def test_specific_execution_route_preserves_preparation(self):
        from pycsamt.agents._workflows import classify_workflow
        self.assertEqual(classify_workflow('Execute prepared Occam2D conventional inversion'),
                         'conventional_inversion')
        self.assertEqual(classify_workflow('prepare Occam2D files'), 'pre_inversion')

    def test_execution_review_cannot_be_disabled_by_config(self):
        plan = WorkflowPlan.from_config({'workflow': 'conventional_inversion',
            'data_path': '.', 'requires_human_review': False})
        self.assertTrue(plan.requires_human_review)

    def test_nested_executable_and_project_bytes_bind_approval(self):
        directory = Path(__file__).resolve().parents[3]/'.tmp'/('approval_native_'+uuid4().hex)
        project = directory/'project'
        project.mkdir(parents=True)
        binary = directory/'solver.exe'
        binary.write_bytes(b'first executable')
        startup = project/'Startup'
        startup.write_bytes(b'first input')
        plan = WorkflowPlan.from_config({'workflow': 'conventional_inversion',
            'data_path': str(project), 'step_params': {'conventional_inversion': {
                'backend_options': {'prepared_project': {'project_dir': str(project),
                    'binary_path': str(binary)}}}}})
        first = plan.approval_digest()
        binary.write_bytes(b'changed executable')
        second = plan.approval_digest()
        self.assertNotEqual(first, second)
        startup.write_bytes(b'changed native input')
        self.assertNotEqual(second, plan.approval_digest())


class FieldInversionSafeguards(unittest.TestCase):
    def test_field_request_routes_to_executed_inversion(self):
        from pycsamt.agents._workflows import classify_workflow
        self.assertEqual(
            classify_workflow('Run the Occam2D inversion of line L26PLT after AMA'),
            'field_inversion')
        # The exact paper request names static shift, rotation and inversion;
        # the compound-request rule used to reroute it to "full", which ran
        # an unrequested AI inversion without review and no Occam2D run.
        self.assertEqual(classify_workflow(
            'Run the Occam2D inversion of line L26PLT (Baohuashan AMT): AMA '
            'static-shift correction in the 10-1100 Hz band, rotation to the '
            '52.5 degree strike, joint TE+TM, 25 iterations.'), 'field_inversion')
        # A genuine compound request without an execution route stays "full".
        self.assertEqual(classify_workflow(
            'run quality control, prepare the inversion, run the AI inversion '
            'and generate a report'), 'full')
        # Preparing files is still preparation, not execution.
        self.assertEqual(classify_workflow('prepare Occam2D files'), 'pre_inversion')

    def test_field_inversion_always_requires_review(self):
        plan = WorkflowPlan.from_config({'workflow': 'field_inversion',
            'data_path': '.', 'requires_human_review': False})
        self.assertTrue(plan.requires_human_review)

    def test_missing_scientific_settings_refused_before_any_step(self):
        from pycsamt.agents.orchestrator import WorkflowOrchestratorAgent
        from pycsamt.api.agents import AGENT_CONFIG
        out = Path(__file__).resolve().parents[3]/'.tmp'/('field_refuse_'+uuid4().hex)
        payload = {'request': 'Run the Occam2D inversion of line L26PLT',
                   'output_dir': str(out),
                   'config': {'data_path': '.', 'workflow': 'field_inversion',
                              'step_params': {'inversion': {'backend': 'occam2d'}}}}
        with AGENT_CONFIG.offline():
            agent = WorkflowOrchestratorAgent()
            first = agent.execute(payload)
            token = (first.data or {}).get('approval_plan_hash')
            result = agent.execute(dict(payload, approved_plan_hash=token)) if token else first
        self.assertEqual(result.status, 'failed')
        self.assertIn('explicit settings', result.error or '')
        self.assertIn('strike_deg', result.error or '')
