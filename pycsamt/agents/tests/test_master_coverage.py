# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.master`.

Complements ``test_master.py`` (lazy import, provider aliases,
workflow/plot dispatch) with the remaining ``AgentMaster.run`` intent
branches (meta, clarify, question, code, metrics) and their dispatch
helpers (``_dispatch_code``, ``_dispatch_metrics``, ``_agent_result``).
``IntentRouter.route`` is stubbed at the class level so each intent is
exercised deterministically without any LLM call.
"""

from __future__ import annotations

from pycsamt.agents import AgentMaster
from pycsamt.agents._base import AgentResult
from pycsamt.agents.router import CLARIFY, CODE, META, METRICS, QUESTION
from pycsamt.agents.router import IntentRouter, RouterDecision


def _stub_route(monkeypatch, decision: RouterDecision):
    monkeypatch.setattr(IntentRouter, "route", lambda self, text, **k: decision)


def test_meta_intent_returns_capability_summary(monkeypatch):
    _stub_route(
        monkeypatch,
        RouterDecision(intent=META, confidence=0.9, source="offline"),
    )
    master = AgentMaster()
    result = master.run("hi there")
    assert result.status == "success"
    assert result.data["intent"] == META
    assert "pyCSAMT" in result.llm_interpretation


def test_clarify_intent_returns_needs_review_with_clarification(monkeypatch):
    _stub_route(
        monkeypatch,
        RouterDecision(
            intent=CLARIFY,
            confidence=0.5,
            clarification="Which line do you mean?",
            source="offline",
        ),
    )
    master = AgentMaster()
    result = master.run("fix the line")
    assert result.status == "needs_review"
    assert result.summary == "Which line do you mean?"


def test_clarify_intent_default_message_when_no_clarification(monkeypatch):
    _stub_route(
        monkeypatch,
        RouterDecision(intent=CLARIFY, confidence=0.5, source="offline"),
    )
    master = AgentMaster()
    result = master.run("do the thing")
    assert result.status == "needs_review"
    assert "clarify" in result.summary.lower()


def test_question_intent_dispatches_to_package_qa(monkeypatch):
    import pycsamt.agents.package_qa as pqa_mod

    _stub_route(
        monkeypatch,
        RouterDecision(intent=QUESTION, confidence=0.8, source="offline"),
    )
    monkeypatch.setattr(
        pqa_mod.PackageQAAgent,
        "execute",
        lambda self, input_data: AgentResult(
            status="success", summary="answer", data={}
        ),
    )
    master = AgentMaster()
    result = master.run("what does StaticShiftAgent do?")
    assert result.status == "success"
    assert result.data["intent"] == QUESTION


def test_code_intent_dispatches_to_code_gen(monkeypatch):
    import pycsamt.agents.code_gen as cg_mod
    import pycsamt.agents.context as ctx_mod

    _stub_route(
        monkeypatch,
        RouterDecision(intent=CODE, confidence=0.85, source="offline"),
    )
    monkeypatch.setattr(
        ctx_mod.ContextInputAgent,
        "execute",
        lambda self, input_data: AgentResult(
            status="success",
            summary="ctx",
            data={"config": {"workflow": "qc"}},
        ),
    )
    captured = {}

    def _fake_cg_execute(self, input_data):
        captured.update(input_data)
        return AgentResult(status="success", summary="code", data={})

    monkeypatch.setattr(
        cg_mod.CodeGenerationAgent, "execute", _fake_cg_execute
    )
    master = AgentMaster()
    result = master.run(
        "write a script to load EDIs",
        data_path="/data/edi",
        output_dir="/out",
    )
    assert result.status == "success"
    assert result.data["intent"] == "code"
    assert captured["workflow_config"]["data_path"] == "/data/edi"
    assert captured["output_dir"] == "/out"


def test_metrics_intent_dispatches_to_metrics_agent(monkeypatch):
    import pycsamt.agents.metrics as metrics_mod

    _stub_route(
        monkeypatch,
        RouterDecision(intent=METRICS, confidence=0.82, source="offline"),
    )
    captured = {}

    def _fake_execute(self, input_data):
        captured.update(input_data)
        return AgentResult(status="success", summary="metrics", data={})

    monkeypatch.setattr(metrics_mod.MetricsAgent, "execute", _fake_execute)
    master = AgentMaster()
    result = master.run(
        "what's the strike of L22PLT?", data_path="/data/edi"
    )
    assert result.status == "success"
    assert result.data["intent"] == "metrics"
    assert captured["data_path"] == "/data/edi"


def test_plan_is_a_dry_run_shortcut(monkeypatch):
    seen = {}

    class _FakeOrch:
        def execute(self, payload):
            seen.update(payload)
            return AgentResult(status="success", summary="planned")

    master = AgentMaster()
    master._orchestrator = _FakeOrch()
    _stub_route(
        monkeypatch,
        RouterDecision(intent="workflow", confidence=0.7, source="offline"),
    )
    master.plan("run QC")
    assert seen["dry_run"] is True
