# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Catalogue corrections as orchestrator workflows (2.6.6).

Each correction id is a registered two-step workflow (load, correction), so a
correction requested in the chat is validated, cross-checked, held at the
review gate when flagged and recorded like any other workflow.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

import pycsamt.agents.orchestrator as orch
from pycsamt.agents import _workflow_plan
from pycsamt.agents._checkpoint_identity import data_identity
from pycsamt.agents._corrections import CORRECTION_METHODS, fn_for
from pycsamt.agents._workflow_plan import WorkflowPlan
from pycsamt.agents.tooling import CorrectionStepAgent


@pytest.mark.parametrize("wf", sorted(CORRECTION_METHODS))
def test_every_correction_is_a_registered_workflow(wf):
    assert wf in _workflow_plan.VALID_WORKFLOWS
    chain = orch._WORKFLOW_STEPS[wf]
    assert [s[0] for s in chain] == ["load", "correction"]
    assert chain[1][1] == "CorrectionStepAgent"


def test_correction_step_agent_never_uses_an_llm():
    agent = CorrectionStepAgent()
    assert agent.llm_available is False
    res = agent.execute({"corr_wf": "corr_ss_ama"})
    assert res.status == "failed"  # no data supplied, nothing invented


def test_plan_accepts_a_list_of_existing_files(tmp_path):
    a, b = tmp_path / "a.edi", tmp_path / "b.edi"
    a.write_text("x")
    b.write_text("y")
    plan = WorkflowPlan("QC", "qc", data_path=[str(a), str(b)])
    assert plan.validation_errors() == []
    ident = data_identity([str(a), str(b)])
    assert [p["path"] for p in ident["paths"]] == [str(a.resolve()), str(b.resolve())]
    plan.approval_digest()  # a list input has a stable identity


def test_plan_reports_missing_files_in_a_list(tmp_path):
    a = tmp_path / "a.edi"
    a.write_text("x")
    plan = WorkflowPlan("QC", "qc", data_path=[str(a), str(tmp_path / "gone.edi")])
    errs = plan.validation_errors()
    assert errs and "gone.edi" in errs[0]


def _copy(edi_dir, tmp_path):
    data = tmp_path / "edi"
    shutil.copytree(edi_dir, data)
    return data


def test_correction_runs_through_the_orchestrator(edi_dir, tmp_path, monkeypatch):
    from pycsamt.api.agents import AGENT_CONFIG

    data = _copy(edi_dir, tmp_path)
    out = tmp_path / "run"
    with AGENT_CONFIG.offline():
        agent = orch.WorkflowOrchestratorAgent()
        monkeypatch.setattr(agent, "_retrieve_context", lambda q: None)
        res = agent.execute({
            "request": "apply ama static shift",
            "data_path": str(data),
            "output_dir": str(out),
            "config": {"workflow": "corr_ss_ama"},
        })
    assert res.status == "success", res.error
    plan = res.data["workflow_plan"]
    assert plan.workflow_type == "corr_ss_ama"
    steps = res.data["result"].data
    corr = steps["correction"]
    assert corr.data["corrected_sites"] is not None
    for name in ("workflow_plan.json", "agent_trace.json",
                 "run_record.json"):
        assert (out / name).exists(), name  # provenance trail
    # inputs are never modified by a correction
    for f in Path(edi_dir).glob("*.edi"):
        assert (data / f.name).read_bytes() == f.read_bytes()


def test_flagged_correction_waits_for_exact_approval(edi_dir, tmp_path, monkeypatch):
    from pycsamt.api.agents import AGENT_CONFIG

    data = _copy(edi_dir, tmp_path)
    payload = {
        "request": "apply ama static shift",
        "data_path": str(data),
        "output_dir": str(tmp_path / "run"),
        "config": {"workflow": "corr_ss_ama",
                   "risk_flags": ["routes disagree"]},
    }
    with AGENT_CONFIG.offline():
        agent = orch.WorkflowOrchestratorAgent()
        monkeypatch.setattr(agent, "_retrieve_context", lambda q: None)
        held = agent.execute(payload)
        assert held.status == "needs_review"
        h = held.data["approval_plan_hash"]
        wrong = agent.execute(dict(payload, approved_plan_hash="0" * 64))
        assert wrong.status == "needs_review"
        ok = agent.execute(dict(payload, approved_plan_hash=h))
    assert ok.status == "success", ok.error
    assert ok.data["result"].data["correction"].data["corrected_sites"] is not None


def test_correction_step_applies_the_workflow_method(monkeypatch, tmp_path):
    from pycsamt.api.agents import AGENT_CONFIG

    seen = {}

    class Loader:
        def execute(self, cfg):
            from pycsamt.agents._base import AgentResult
            return AgentResult("success", "loaded", {"sites": "SITES"})

    class Corr:
        llm_available = False
        llm_provider = None
        model = None

        def execute(self, cfg):
            from pycsamt.agents._base import AgentResult
            seen.update(cfg)
            return AgentResult("success", "corrected", {"corrected_sites": "C"})

    data = tmp_path / "edi"
    data.mkdir()
    monkeypatch.setattr(orch, "_build_registry", lambda **k: (
        {"MTLoaderAgent": Loader(), "CorrectionStepAgent": Corr()}, {}))
    monkeypatch.setattr(orch, "_write_provenance", lambda **k: True)
    with AGENT_CONFIG.offline():
        agent = orch.WorkflowOrchestratorAgent()
        monkeypatch.setattr(agent, "_retrieve_context", lambda q: None)
        res = agent.execute({
            "request": "loess static shift correction",
            "data_path": str(data),
            "output_dir": str(tmp_path / "run"),
            "config": {"workflow": "corr_ss_loess",
                       "step_params": {"correction": {"corr_wf": "corr_ss_ama",
                                                      "stations": ["S01"]}}},
        })
    assert res.status == "success", res.error
    assert seen["corr_wf"] == "corr_ss_loess"  # the plan's id, never another
    assert seen["fn_name"] == fn_for("corr_ss_loess")
    assert seen["stations"] == ["S01"]
    assert seen["sites"] == "SITES"
