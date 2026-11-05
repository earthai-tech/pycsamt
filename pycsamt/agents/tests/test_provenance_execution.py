# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Observed provenance and request-driven execution evaluation regressions."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from uuid import uuid4

from pycsamt.agents._base import AgentResult, BaseAgent
from pycsamt.agents._provenance import provenance_scope, current_record, emit, record_retrieval
from pycsamt.agents.coordinator import AgentCoordinator
from pycsamt.assistant.evals.execution import evaluate_case, summarize


class Dummy:
    # Deliberately opaque to disable unrelated checkpoint-resume optimization.
    def __init__(self, **settings):
        self.status = settings.get("status", "success")
        self.value = settings.get("value", 3)
        self.calls = []

    def execute(self, payload):
        self.calls.append(payload)
        output = Path(payload["output_dir"])
        output.mkdir(parents=True, exist_ok=True)
        (output / "unrelated_name.json").write_text(json.dumps({"value": self.value}))
        return AgentResult(status=self.status, summary="dummy", data={"value": self.value}, warnings=["test warning"])


class LanguageDummy(BaseAgent):
    def __init__(self):
        super().__init__("test_language", api_key="test-private-key", llm_provider="openai", model="test-model")

    def execute(self, payload):
        return AgentResult("success", "unused")


class ProvenanceExecutionTests(unittest.TestCase):
    def setUp(self):
        self.directory = Path(__file__).resolve().parents[3] / ".tmp/provenance_tests" / uuid4().hex
        self.directory.mkdir(parents=True)

    def runner(self, payload, status="success"):
        coord = AgentCoordinator("probe", checkpoint_dir=self.directory / ("checkpoints_" + uuid4().hex), verbose=False)
        coord.add_step("calculation", Dummy(status=status))
        result = coord.execute(payload)
        result.data["workflow_type"] = "qc"
        return result

    def case(self, **extra):
        return {"id": "test", "query": "QC this input", "execution_allowed": True,
                "expected_workflow": "qc", "requires_review": False, **extra}

    def test_records_only_dispatched_steps_and_measured_outputs(self):
        with provenance_scope({"request": "run", "history": [{"role": "user", "content": "earlier"}]}) as record:
            coord = AgentCoordinator("probe", checkpoint_dir=self.directory / "checkpoints", verbose=False)
            coord.add_step("first", Dummy(status="failed"))
            never = Dummy()
            coord.add_step("never", never)
            result = coord.execute({"output_dir": str(self.directory / "outputs")})
            record.status = result.status
        events = record.snapshot()["events"]
        self.assertEqual([e["step"] for e in events if e["event"] == "dispatch_started"], ["first"])
        finished = next(e for e in events if e["event"] == "dispatch_finished")
        self.assertEqual(finished["status"], "failed")
        self.assertIn("unrelated_name.json", finished["artifact_changes"])
        self.assertIn("elapsed_seconds", finished)
        self.assertEqual(never.calls, [])
        self.assertTrue(events[0]["history_provided"])

    def test_language_records_exact_inputs_and_redacts_credentials(self):
        from pycsamt.api.agents import AGENT_CONFIG

        # Explicit mock transport: no provider SDK or network call.
        with AGENT_CONFIG.using(provider="openai", api_key="test-private-key", model="test-model"):
            agent = LanguageDummy()
            agent._offline_at_creation = False
            agent.api_key = "test-private-key"
            with provenance_scope({"request": "test"}) as record, \
                 patch.object(type(agent), "llm_available", new=property(lambda self: True)), \
                 patch.object(agent, "_query_openai", return_value=("answer test-private-key", 0.0)):
                reply = agent.query_llm("actual prompt", "actual system", temperature=0.0, max_tokens=17)
        self.assertIn("test-private-key", reply)
        snapshot = record.snapshot()
        self.assertNotIn("test-private-key", json.dumps(snapshot))
        call = next(e for e in snapshot["events"] if e["event"] == "language_call_started")
        self.assertEqual(call["prompt"], "actual prompt")
        self.assertEqual(call["system_message"], "actual system")
        self.assertEqual(call["temperature"], 0.0)
        self.assertEqual(call["max_tokens"], 17)
        self.assertIsNone(call["provider_model_version"])
        self.assertTrue(any(e["event"] == "language_transport_attempt" for e in snapshot["events"]))

    def test_nested_scope_keeps_one_run_and_isolates_later_requests(self):
        with provenance_scope({"request": "first"}) as first:
            with provenance_scope({"request": "nested"}) as nested:
                self.assertIs(first, nested)
        self.assertIsNone(current_record())
        with provenance_scope({"request": "second"}) as second:
            self.assertNotEqual(first.run_id, second.run_id)

    def test_persisted_record_contains_final_outcome(self):
        path = self.directory / "record.json"
        with provenance_scope({"request": "test"}) as record:
            record.persist(path)
            record.status = "needs_review"
        saved = json.loads(path.read_text())
        self.assertEqual(saved["status"], "needs_review")
        self.assertEqual(saved["events"][-1]["event"], "request_finished")

    def test_expected_labels_and_references_never_enter_executor_payload(self):
        seen = {}
        def executor(payload):
            seen.update(payload)
            return self.runner(payload)
        result = evaluate_case(self.case(expected_symbols=["secret-reference"], expected_value=3), executor,
                               self.directory / "case", numerical_validator=lambda r, p: {"passed": True, "value": 3})
        self.assertEqual(set(seen), {"request", "output_dir"})
        self.assertTrue(result["numerically_verified_completion"])
        self.assertFalse(result["verified_completion"])

    def test_successful_exit_does_not_imply_numerical_acceptance(self):
        report = evaluate_case(self.case(), self.runner, self.directory / "case")
        self.assertEqual(report["observed_status"], "success")
        self.assertIsNone(report["numerical_acceptance"])
        self.assertFalse(report["verified_completion"])

    def test_failed_numeric_reference_prevents_verified_completion(self):
        report = evaluate_case(self.case(), self.runner, self.directory / "case",
                               numerical_validator=lambda r, p: {"passed": False, "max_error": 2.0})
        self.assertFalse(report["verified_completion"])
        self.assertFalse(report["numerical_acceptance"])

    def test_missing_required_artifact_prevents_verified_completion(self):
        report = evaluate_case(self.case(required_outputs=["missing.npz"]), self.runner, self.directory / "case",
                               numerical_validator=lambda r, p: {"passed": True})
        self.assertFalse(report["output_completeness"])
        self.assertFalse(report["verified_completion"])

    def test_actual_prohibited_dispatch_is_counted(self):
        report = evaluate_case(self.case(execution_allowed=False), self.runner, self.directory / "case")
        self.assertTrue(report["actual_execution_violation"])
        self.assertEqual(report["registered_dispatch_count"], 1)
        self.assertFalse(report["verified_completion"])

    def test_unknown_uninstrumented_execution_does_not_become_safe(self):
        report = evaluate_case(self.case(execution_allowed=False), lambda p: AgentResult("success", "no evidence"), self.directory / "case")
        self.assertFalse(report["observation_complete"])
        self.assertIsNone(report["actual_execution_violation"])
        summary = summarize([report])
        self.assertIsNone(summary["actual_execution_violation_rate"])
        self.assertEqual(summary["nonexecutable_observed_n"], 0)

    def test_incomplete_dispatch_preserves_known_violation(self):
        def executor(payload):
            emit("coordinator_boundary", workflow="test")
            emit("dispatch_started", call_id="unfinished", step="calculation", agent="Dummy")
            raise RuntimeError("interrupted")
        report = evaluate_case(self.case(execution_allowed=False), executor, self.directory / "case")
        self.assertFalse(report["observation_complete"])
        self.assertTrue(report["actual_execution_violation"])
        self.assertEqual(report["observed_dispatch_attempts"], 1)

    def test_unrequested_step_prevents_verified_completion(self):
        report = evaluate_case(self.case(allowed_steps=["load", "report"]), self.runner, self.directory / "case",
                               numerical_validator=lambda r, p: {"passed": True})
        self.assertFalse(report["allowed_step_conformity"])
        self.assertFalse(report["verified_completion"])

    def test_completed_direct_orchestrator_trace_contains_actual_events(self):
        import pycsamt.agents.orchestrator as module
        from pycsamt.api.agents import AGENT_CONFIG

        raw = self.directory / "input.edi"
        raw.write_text("test")
        with AGENT_CONFIG.offline():
            agent = module.WorkflowOrchestratorAgent()
        payload = {"request": "QC", "data_path": str(raw), "output_dir": str(self.directory / "outputs"),
                   "config": {"workflow": "qc"}}
        with patch.object(module, "_WORKFLOW_STEPS", {"qc": [("computation", "Dummy", None, "test")]}), \
             patch.object(module, "_build_registry", return_value=({"Dummy": Dummy()}, {})), \
             patch.object(agent, "_retrieve_context", return_value=None):
            result = agent.execute(payload)
        self.assertEqual(result.status, "success")
        trace = json.loads((Path(payload["output_dir"]) / "agent_trace.json").read_text())
        self.assertEqual(trace["executed_agents"], ["Dummy"])
        manifest = json.loads((Path(payload["output_dir"]) / "output_manifest.json").read_text())
        file = next(f for f in manifest["files"] if f["path"] == "unrelated_name.json")
        self.assertEqual(file["producing_step"], "computation")
        self.assertEqual(file["attribution"], "observed_change")
        saved = json.loads((Path(payload["output_dir"]) / "run_record.json").read_text())
        self.assertEqual(saved["status"], result.status)
        self.assertEqual(saved["events"][-1]["event"], "request_finished")

    def test_provenance_write_failure_is_visible(self):
        import pycsamt.agents.orchestrator as module
        from pycsamt.api.agents import AGENT_CONFIG

        raw = self.directory / "input.edi"
        raw.write_text("test")
        with AGENT_CONFIG.offline():
            agent = module.WorkflowOrchestratorAgent()
        with patch.object(module, "_WORKFLOW_STEPS", {"qc": [("first", "Dummy", None, "test")]}), \
             patch.object(module, "_build_registry", return_value=({"Dummy": Dummy()}, {})), \
             patch.object(agent, "_retrieve_context", return_value=None), \
             patch.object(module, "_write_provenance", return_value=False):
            result = agent.execute({"request": "QC", "data_path": str(raw), "config": {"workflow": "qc"},
                                    "output_dir": str(self.directory / "outputs")})
        self.assertEqual(result.status, "needs_review")
        self.assertFalse(result.data["recording_complete"])

    def test_recorded_request_preserves_keyword_call_signature(self):
        from pycsamt.agents._provenance import recorded_request

        class Endpoint:
            @recorded_request
            def execute(self, input_data):
                return AgentResult("success", "ok", data={"request": input_data["request"]})
        result = Endpoint().execute(input_data={"request": "keyword request"})
        self.assertEqual(result.data["request"], "keyword request")
        self.assertEqual(result.data["provenance"]["events"][0]["request"], "keyword request")

    def test_retrieval_records_active_corpus_and_context(self):
        from pycsamt.assistant.rag.schemas import RAGChunk

        chunk = RAGChunk(id="id", text="retrieved content", source_path="source.py", kind="python_symbol")
        context = SimpleNamespace(chunks=[chunk], context_text="actual assembled context", citations=[{"source_path": "source.py"}])
        with provenance_scope({"request": "test"}) as record:
            record_retrieval("agent", "query", context, SimpleNamespace(chunks=[chunk]))
        event = next(e for e in record.events if e["event"] == "retrieval")
        self.assertEqual(event["context"], "actual assembled context")
        self.assertEqual(event["ranked_chunks"][0]["id"], "id")
        self.assertEqual(len(event["active_corpus_sha256"]), 64)
        self.assertFalse(event["persisted_manifest_is_active_corpus_proof"])

    def test_explicit_scientific_acceptance_is_separate_from_fidelity(self):
        numerical = evaluate_case(self.case(), self.runner, self.directory / "numeric",
                                  numerical_validator=lambda r, p: {"passed": True})
        self.assertTrue(numerical["numerically_verified_completion"])
        self.assertIsNone(numerical["scientific_acceptance"])
        self.assertFalse(numerical["verified_completion"])
        accepted = evaluate_case(self.case(), self.runner, self.directory / "accepted",
                                 numerical_validator=lambda r, p: {"passed": True, "scientific_acceptance": True})
        self.assertTrue(accepted["verified_completion"])

    def test_evaluation_refuses_to_overwrite_a_previous_case(self):
        target = self.directory / "case"
        evaluate_case(self.case(), self.runner, target)
        with self.assertRaises(FileExistsError):
            evaluate_case(self.case(), self.runner, target)

    def test_qc_request_does_not_dispatch_static_shift(self):
        from pycsamt.agents.orchestrator import WorkflowOrchestratorAgent
        import pycsamt.agents.orchestrator as module
        from pycsamt.api.agents import AGENT_CONFIG

        raw = self.directory / "input.edi"
        raw.write_text("test")
        seen = []
        class Step:
            def __init__(self, **kwargs):
                pass
            def execute(self, payload):
                seen.append(type(self).__name__)
                return AgentResult("success", "test", data={"sites": "test-sites"})
        Loader = type("MTLoaderAgent", (Step,), {})
        QC = type("DataQCAgent", (Step,), {})
        Report = type("ReportAgent", (Step,), {})
        registry = {"MTLoaderAgent": Loader(), "DataQCAgent": QC(), "ReportAgent": Report()}
        with AGENT_CONFIG.offline():
            agent = WorkflowOrchestratorAgent()
        with patch.object(module, "_build_registry", return_value=(registry, {})), \
             patch.object(agent, "_retrieve_context", return_value=None):
            result = agent.execute({"request": "run quality control", "data_path": str(raw),
                                    "output_dir": str(self.directory / "outputs")})
        self.assertEqual(result.status, "success")
        self.assertEqual(seen, ["MTLoaderAgent", "DataQCAgent", "ReportAgent"])
        self.assertIn("static_shift", [step[0] for step in module._WORKFLOW_STEPS["static_shift"]])


if __name__ == "__main__":
    unittest.main()
