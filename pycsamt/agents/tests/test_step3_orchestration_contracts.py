from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.agents import ContextInputAgent, WorkflowOrchestratorAgent
from pycsamt.agents._base import AgentResult, BaseAgent
from pycsamt.agents._workflows import classify_workflow
from pycsamt.agents.coordinator import AgentCoordinator


class _FakeAgent(BaseAgent):
    def __init__(self, result: AgentResult, *, name: str = "fake") -> None:
        super().__init__(name)
        self.result = result
        self.calls: list[dict] = []

    def execute(self, input_data: dict) -> AgentResult:
        self.calls.append(dict(input_data))
        return self.result


def _success(**data) -> AgentResult:
    return AgentResult(status="success", summary="ok", data=data)


def _checkpoint_dir(name: str) -> Path:
    path = Path("test_tmp_step3") / name
    path.mkdir(parents=True, exist_ok=True)
    return path


def test_context_and_orchestrator_use_registry_for_specific_workflows():
    probes = [
        ("run QC on /data/edis", "qc"),
        ("correct static shift for /data/edis", "static_shift"),
        ("phase tensor strike analysis on /data/edis", "phase_analysis"),
        ("write code for static shift", "code_gen"),
        (
            "run quality control, prepare the inversion, run AI inversion "
            "and generate a report for /data/edis",
            "full",
        ),
    ]

    context = ContextInputAgent()
    orchestrator = WorkflowOrchestratorAgent()

    for request, expected in probes:
        assert classify_workflow(request, default="qc") == expected

        ctx_result = context.execute({"request": request})
        assert ctx_result.status in {"success", "needs_review"}
        assert ctx_result["config"]["workflow"] == expected

        orch_result = orchestrator.execute(
            {
                "request": request,
                "data_path": "/data/edis",
                "dry_run": True,
            }
        )
        assert orch_result.status in {"success", "needs_review"}
        assert orch_result["workflow_type"] == expected


def test_orchestrator_returns_agent_result_contract_for_dry_run():
    result = WorkflowOrchestratorAgent().execute(
        {
            "request": "QC the survey",
            "data_path": "/data/edis",
            "output_dir": str(_checkpoint_dir("orchestrator")),
            "dry_run": True,
        }
    )

    assert isinstance(result, AgentResult)
    assert result.status == "success"
    assert "orchestrator routed" in result.summary.lower()
    assert result["workflow_type"] == "qc"
    assert result["result"].status == "success"
    assert isinstance(result["steps"], list)
    assert result["steps"][0]["name"] == "load"
    assert result.cost_estimate_usd == pytest.approx(0.0)


def test_required_step_failure_aborts_and_preserves_partial_results():
    first = _FakeAgent(_success(sites="loaded"), name="load")
    failing = _FakeAgent(
        AgentResult.failed("qc failed", hint="inspect station errors"),
        name="qc",
    )
    never = _FakeAgent(_success(report="done"), name="report")
    coordinator = AgentCoordinator(
        "failure_contract",
        checkpoint_dir=_checkpoint_dir("required_failure"),
        verbose=False,
    )
    coordinator.add_step("load", first)
    coordinator.add_step(
        "qc",
        failing,
        input_fn=lambda results: {"sites": results["load"]["sites"]},
    )
    coordinator.add_step("report", never)

    result = coordinator.execute({"path": "/data/edis"})

    assert result.status == "failed"
    assert result.error == "qc failed"
    assert result.error_fix_hint == "inspect station errors"
    assert set(result.data) == {"load", "qc"}
    assert result.data["load"].status == "success"
    assert result.data["qc"].status == "failed"
    assert first.calls == [{"path": "/data/edis"}]
    assert failing.calls == [{"sites": "loaded"}]
    assert never.calls == []


def test_optional_step_failure_yields_needs_review_and_continues():
    optional = _FakeAgent(AgentResult.failed("optional plot failed"))
    final = _FakeAgent(_success(report="done"))
    coordinator = AgentCoordinator(
        "optional_contract",
        checkpoint_dir=_checkpoint_dir("optional_failure"),
        verbose=False,
    )
    coordinator.add_step("optional_plot", optional, required=False)
    coordinator.add_step("report", final)

    result = coordinator.execute({"path": "/data/edis"})

    assert result.status == "needs_review"
    assert result.data["optional_plot"].status == "failed"
    assert result.data["report"].status == "success"
    assert final.calls == [{"path": "/data/edis"}]
