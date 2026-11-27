# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""The Agent Master chat applies the paper's safeguards (2.6.6).

* the orchestrator cross-check runs on the workflow the chat picks;
* a plan held at the review gate is shown with an Approve button, and
  approving resends exactly that plan with its approval hash;
* catalogue corrections run through the orchestrator (see
  ``test_run_agent_dispatch``), so they are validated, gated and recorded.
"""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")

import pycsamt.agents.context as context_mod
import pycsamt.agents.orchestrator as orch_mod
import pycsamt.agents.router as router_mod
import pycsamt.agents._workflows as wf_mod
import pycsamt.app.agent_master.callbacks.chat as C
from pycsamt.agents._base import AgentResult
from pycsamt.agents.router import WORKFLOW, RouterDecision


class _SyncThread:
    def __init__(self, target=None, args=(), kwargs=None, daemon=None):
        self._t, self._a, self._k = target, args, kwargs or {}

    def start(self):
        self._t(*self._a, **self._k)


@pytest.fixture
def routed(monkeypatch):
    """Route every request to the workflow intent with no LLM involved."""
    monkeypatch.setattr(
        router_mod.IntentRouter, "route",
        lambda self, text, history=None: RouterDecision(
            intent=WORKFLOW, workflow=None, confidence=0.9, source="offline"),
    )
    monkeypatch.setattr(
        context_mod.ContextInputAgent, "execute",
        lambda self, d: AgentResult("success", "ctx", {"config": {}}),
    )
    monkeypatch.setattr(C.threading, "Thread", _SyncThread)
    C._PENDING_PLANS.clear()
    yield monkeypatch
    C._PENDING_PLANS.clear()


def _held(hash_="ab" * 32):
    from pycsamt.agents._workflow_plan import WorkflowPlan

    plan = WorkflowPlan("denoise", "qc", data_path="/tmp/edis",
                        risk_flags=["LLM chose workflow 'qc' but the keyword route reads 'denoise'"])
    return AgentResult("needs_review",
                       "Workflow requires approval of this exact plan before execution.",
                       {"workflow_type": "qc", "workflow_plan": plan,
                        "approval_plan_hash": hash_,
                        "steps": [{"name": "load"}, {"name": "qc"}, {"name": "report"}]})


def test_cross_check_flags_a_form_choice_that_contradicts_the_text(routed):
    calls = []

    def _fake(self, d):
        calls.append(d)
        return AgentResult("success", "done", {"result": AgentResult("success", "i", {})})

    routed.setattr(orch_mod.WorkflowOrchestratorAgent, "execute", _fake)
    jid = C._new_job()
    # the form says QC, the text asks for denoising
    C._run_agent(jid, "denoise the survey data", {"path": "/tmp/edis"},
                 {"provider": "offline"}, {"workflow": "qc"})
    flags = calls[0]["config"].get("risk_flags") or []
    assert any("keyword route" in f for f in flags)


def test_cross_check_adds_nothing_when_routes_agree(routed):
    calls = []
    routed.setattr(orch_mod.WorkflowOrchestratorAgent, "execute",
                   lambda self, d: calls.append(d) or AgentResult(
                       "success", "done", {"result": AgentResult("success", "i", {})}))
    jid = C._new_job()
    C._run_agent(jid, "run QC on the data", {"path": "/tmp/edis"}, {"provider": "offline"})
    assert not calls[0]["config"].get("risk_flags")


def test_held_plan_is_shown_with_an_approve_button(routed):
    routed.setattr(orch_mod.WorkflowOrchestratorAgent, "execute", lambda self, d: _held())
    jid = C._new_job()
    C._run_agent(jid, "run QC on the data", {"path": "/tmp/edis"}, {"provider": "offline"})
    job = C._get_job(jid)
    assert job["kind"] == C.KIND_CLARIFY
    assert "needs your approval" in job["result"]
    assert "keyword route" in job["result"]
    token = job["approval"]["token"]
    assert job["approval"]["hash"] == "ab" * 32
    assert token in C._PENDING_PLANS
    bubble = C._agent_bubble(job["result"], kind=job["kind"], approval=job["approval"])
    assert "am-approve-plan" in str(bubble) and token in str(bubble)


def test_approve_resends_the_same_plan_with_its_hash(routed):
    calls = []

    def _fake(self, d):
        calls.append(dict(d))
        if d.get("approved_plan_hash") == "ab" * 32:
            return AgentResult("success", "ran", {"result": AgentResult("success", "i", {})})
        return _held()

    routed.setattr(orch_mod.WorkflowOrchestratorAgent, "execute", _fake)
    jid = C._new_job()
    C._run_agent(jid, "run QC on the data", {"path": "/tmp/edis"}, {"provider": "offline"})
    token = C._get_job(jid)["approval"]["token"]

    msgs, job_store, disabled, stored = C._approve_plan_response(token, [], [])
    assert disabled is False and "jid" in job_store
    assert stored[-1]["approved_plan_hash"] == "ab" * 32
    second = calls[-1]
    assert second["approved_plan_hash"] == "ab" * 32
    first = {k: v for k, v in calls[0].items()}
    assert {k: v for k, v in second.items() if k != "approved_plan_hash"} == first
    done = C._get_job(job_store["jid"])
    assert done["kind"] == C.KIND_WORKFLOW
    assert token not in C._PENDING_PLANS  # one approval, one run


def test_unknown_or_used_approval_runs_nothing(routed):
    routed.setattr(orch_mod.WorkflowOrchestratorAgent, "execute",
                   lambda self, d: pytest.fail("nothing may run"))
    msgs, job_store, disabled, stored = C._approve_plan_response("nope", [], [])
    assert disabled is True
    assert "no longer valid" in stored[-1]["content"]


def _input_type(entry):
    import json

    if not entry["inputs"]:
        return None
    cid = entry["inputs"][0]["id"]
    if isinstance(cid, str) and cid.startswith("{"):
        cid = json.loads(cid)  # dash stores pattern ids as JSON text
    return cid.get("type") if isinstance(cid, dict) else None


def _approve_cb(agent_app):
    keys = [k for k, e in agent_app.callback_map.items()
            if _input_type(e) == "am-approve-plan"]
    assert len(keys) == 1, keys
    fn = agent_app.callback_map[keys[0]]["callback"]
    return getattr(fn, "__wrapped__", fn)


def _trigger(value):
    import json

    import dash._callback_context as cc
    from dash._utils import AttributeDict

    prop = json.dumps({"token": "t", "type": "am-approve-plan"}, sort_keys=True)
    cc.context_value.set(AttributeDict(
        triggered_inputs=[{"prop_id": prop + ".n_clicks", "value": value}]))


@pytest.mark.parametrize("value,poll_disabled", [(0, True), (1, False)])
def test_approve_callback_ignores_render_and_running_jobs(agent_app, value, poll_disabled):
    from dash.exceptions import PreventUpdate

    fn = _approve_cb(agent_app)
    _trigger(value)
    try:
        with pytest.raises(PreventUpdate):
            fn([value], [], [], poll_disabled)
    finally:
        import dash._callback_context as cc
        cc.context_value.set({})


def test_approve_callback_starts_the_held_plan(agent_app, monkeypatch):
    started = {}
    monkeypatch.setattr(C, "_approve_plan_response",
                        lambda token, msgs, stored: started.setdefault("token", token) and
                        ([], {"jid": "j"}, False, []))
    fn = _approve_cb(agent_app)
    _trigger(1)
    try:
        out = fn([1], [], [], True)
    finally:
        import dash._callback_context as cc
        cc.context_value.set({})
    assert started["token"] == "t"
    assert out[1] == {"jid": "j"} and out[2] is False
