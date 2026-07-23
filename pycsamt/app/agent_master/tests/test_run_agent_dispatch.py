# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the _run_agent routing tree in
pycsamt.app.agent_master.callbacks.chat — the core job-execution engine.

Strategy: _run_agent is a plain function (no Dash context needed), so it is
called directly. Heavy dependencies (IntentRouter, ContextInputAgent,
WorkflowOrchestratorAgent, and the smaller _dispatch_* helpers) are
monkeypatched at the class/module level so each test isolates one branch
of the routing tree without touching real LLM/agent execution.
"""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")

import pycsamt.agents.context as context_mod
import pycsamt.agents.orchestrator as orch_mod
import pycsamt.agents.router as router_mod
import pycsamt.app.agent_master.callbacks.chat as C
from pycsamt.agents._base import AgentResult
from pycsamt.agents.router import CLARIFY, CODE, META, METRICS, QUESTION, WORKFLOW


def _patch_router(monkeypatch, intent, workflow=None, clarification=None):
    from pycsamt.agents.router import RouterDecision

    monkeypatch.setattr(
        router_mod.IntentRouter,
        "route",
        lambda self, text, history=None: RouterDecision(
            intent=intent,
            workflow=workflow,
            confidence=0.9,
            clarification=clarification,
            source="offline",
        ),
    )


def _patch_context(monkeypatch, config=None, status="success"):
    monkeypatch.setattr(
        context_mod.ContextInputAgent,
        "execute",
        lambda self, input_data: AgentResult(
            status, "ctx", {"config": config or {}}
        ),
    )


def _patch_orchestrator(monkeypatch, result):
    monkeypatch.setattr(
        orch_mod.WorkflowOrchestratorAgent,
        "execute",
        lambda self, input_data: result,
    )


def _new_job():
    return C._new_job()


def test_data_read_text_dispatches_to_overview(monkeypatch):
    called = {}
    monkeypatch.setattr(
        C,
        "_dispatch_data_overview",
        lambda jid, text, edi_store, settings, step: called.setdefault(
            "hit", True
        ),
    )
    jid = _new_job()
    C._run_agent(jid, "read the loaded data", {}, {"provider": "offline"})
    assert called.get("hit") is True


class TestIntentRouting:
    def test_meta_intent_returns_capability_text(self, monkeypatch):
        _patch_router(monkeypatch, META)
        jid = _new_job()
        C._run_agent(jid, "what can you do", {}, {"provider": "offline"})
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_META
        assert job["status"] == "done"

    def test_clarify_intent_uses_clarification_text(self, monkeypatch):
        _patch_router(monkeypatch, CLARIFY, clarification="Which line?")
        jid = _new_job()
        C._run_agent(jid, "run it", {}, {"provider": "offline"})
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_CLARIFY
        assert job["result"] == "Which line?"

    def test_clarify_intent_falls_back_to_default_prompt(self, monkeypatch):
        _patch_router(monkeypatch, CLARIFY, clarification=None)
        jid = _new_job()
        C._run_agent(jid, "run it", {}, {"provider": "offline"})
        job = C._get_job(jid)
        assert "clarify" in job["result"].lower()

    def test_metrics_intent_dispatches(self, monkeypatch):
        called = {}
        monkeypatch.setattr(
            C,
            "_dispatch_metrics",
            lambda jid, text, edi_store, settings, step: called.setdefault(
                "hit", True
            ),
        )
        _patch_router(monkeypatch, METRICS)
        jid = _new_job()
        # avoid "stat..." substring words that trip the data-read gate
        C._run_agent(
            jid,
            "compute the mean resistivity for line 1",
            {},
            {"provider": "offline"},
        )
        assert called.get("hit") is True

    def test_question_intent_dispatches(self, monkeypatch):
        called = {}
        monkeypatch.setattr(
            C,
            "_dispatch_question",
            lambda jid, text, **kw: called.setdefault("hit", True),
        )
        _patch_router(monkeypatch, QUESTION)
        jid = _new_job()
        C._run_agent(
            jid, "what does StaticShiftAgent do", {}, {"provider": "offline"}
        )
        assert called.get("hit") is True

    def test_code_intent_dispatches(self, monkeypatch):
        called = {}
        monkeypatch.setattr(
            C,
            "_dispatch_code",
            lambda jid, text, edi_store, settings, **kw: called.setdefault(
                "hit", True
            ),
        )
        _patch_router(monkeypatch, CODE)
        jid = _new_job()
        C._run_agent(
            jid, "generate code for static shift", {}, {"provider": "offline"}
        )
        assert called.get("hit") is True


class TestContextClassificationFailure:
    def test_failed_context_result_reports_error(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, status="failed")
        jid = _new_job()
        C._run_agent(jid, "run qc pipeline", {}, {"provider": "offline"})
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ERROR
        assert "could not classify" in job["result"].lower()


class TestWorkflowResolutionAuthority:
    """explicit (param modal) > keyword registry > router slot > ctx classification."""

    def _run(self, monkeypatch, *, inv_config=None, router_wf=None, ctx_wf=None,
              kw_wf_return="_UNSET_"):
        _patch_router(monkeypatch, WORKFLOW, workflow=router_wf)
        _patch_context(monkeypatch, config=({"workflow": ctx_wf} if ctx_wf else {}))
        if kw_wf_return != "_UNSET_":
            import pycsamt.agents._workflows as wf_mod

            monkeypatch.setattr(
                wf_mod, "classify_workflow", lambda text, default=None: kw_wf_return
            )
        captured = {}
        _patch_orchestrator(
            monkeypatch,
            AgentResult(
                "success",
                "done",
                {"result": AgentResult("success", "inner", {})},
            ),
        )

        def _fake_execute(self, input_data):
            captured["cfg"] = dict(input_data.get("config", {}))
            return AgentResult(
                "success",
                "done",
                {"result": AgentResult("success", "inner", {})},
            )

        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent, "execute", _fake_execute
        )
        jid = _new_job()
        C._run_agent(
            jid,
            "run something",
            {"path": "/tmp/edis"},
            {"provider": "offline"},
            inv_config or {},
        )
        return captured, C._get_job(jid)

    def test_explicit_form_selection_wins(self, monkeypatch):
        captured, job = self._run(
            monkeypatch,
            inv_config={"workflow": "static_shift"},
            router_wf="ai_inversion",
            ctx_wf="qc",
            kw_wf_return=None,
        )
        assert captured["cfg"]["workflow"] == "static_shift"

    def test_keyword_registry_wins_over_router_and_ctx(self, monkeypatch):
        captured, job = self._run(
            monkeypatch,
            router_wf="ai_inversion",
            ctx_wf="phase_analysis",
            kw_wf_return="static_shift",
        )
        assert captured["cfg"]["workflow"] == "static_shift"

    def test_router_slot_wins_when_no_keyword_match(self, monkeypatch):
        captured, job = self._run(
            monkeypatch,
            router_wf="phase_analysis",
            ctx_wf="qc",
            kw_wf_return=None,
        )
        assert captured["cfg"]["workflow"] == "phase_analysis"

    def test_ctx_classification_used_as_last_resort(self, monkeypatch):
        captured, job = self._run(
            monkeypatch,
            router_wf=None,
            ctx_wf="phase_analysis",
            kw_wf_return=None,
        )
        assert captured["cfg"]["workflow"] == "phase_analysis"

    def test_qc_only_used_when_keyword_matched_qc(self, monkeypatch):
        # ctx_wf defaults to "qc" often; that alone must NOT resolve the
        # workflow — only an explicit qc keyword match should.
        captured, job = self._run(
            monkeypatch,
            router_wf=None,
            ctx_wf="qc",
            kw_wf_return="qc",
        )
        assert captured["cfg"]["workflow"] == "qc"


class TestUnknownWorkflow:
    def test_unresolved_workflow_returns_clarify_or_meta(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod, "classify_workflow", lambda text, default=None: None
        )
        jid = _new_job()
        C._run_agent(
            jid, "asdkjfhaslkdjfh gibberish nonsense", {}, {"provider": "offline"}
        )
        job = C._get_job(jid)
        assert job["status"] == "done"
        assert job["kind"] in (C.KIND_CLARIFY, C.KIND_META)


class TestPreInversionCodeSwitch:
    def test_modem_code_selection_switches_workflow(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "pre_inversion"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod,
            "classify_workflow",
            lambda text, default=None: "pre_inversion",
        )
        called = {}
        monkeypatch.setattr(
            C,
            "_dispatch_inversion_prep",
            lambda jid, wtype, *a, **kw: called.setdefault("wtype", wtype),
        )
        jid = _new_job()
        C._run_agent(
            jid,
            "prep modem inversion inputs",
            {"path": "/tmp/edis"},
            {"provider": "offline"},
            {"inversion_code": "modem"},
        )
        assert called.get("wtype") == "modem"


class TestParamInjection:
    def test_pinn_params_injected(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "pinn_inversion"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod,
            "classify_workflow",
            lambda text, default=None: "pinn_inversion",
        )
        captured = {}

        def _fake_execute(self, input_data):
            captured["cfg"] = dict(input_data.get("config", {}))
            return AgentResult(
                "success",
                "done",
                {"result": AgentResult("success", "inner", {})},
            )

        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent, "execute", _fake_execute
        )
        jid = _new_job()
        C._run_agent(
            jid,
            "run pinn inversion",
            {"path": "/tmp/edis"},
            {"provider": "offline"},
            {"n_layers": 12, "checkpoint": "ckpt.pt"},
        )
        assert captured["cfg"]["pinn_params"]["n_layers"] == 12
        assert captured["cfg"]["hybrid_params"]["n_layers"] == 12
        assert captured["cfg"]["checkpoint"] == "ckpt.pt"

    def test_ai_inversion_params_injected(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "ai_inversion"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod,
            "classify_workflow",
            lambda text, default=None: "ai_inversion",
        )
        captured = {}

        def _fake_execute(self, input_data):
            captured["cfg"] = dict(input_data.get("config", {}))
            return AgentResult(
                "success",
                "done",
                {"result": AgentResult("success", "inner", {})},
            )

        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent, "execute", _fake_execute
        )
        jid = _new_job()
        C._run_agent(
            jid,
            "run 1d ai inversion",
            {"path": "/tmp/edis"},
            {"provider": "offline"},
            {"n_layers": 8, "epochs": 200},
        )
        assert captured["cfg"]["ai_inv_params"]["n_layers"] == 8
        assert captured["cfg"]["ai_inv_params"]["epochs"] == 200

    def test_step_params_passed_through(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "qc"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod, "classify_workflow", lambda text, default=None: "qc"
        )
        captured = {}

        def _fake_execute(self, input_data):
            captured["cfg"] = dict(input_data.get("config", {}))
            return AgentResult(
                "success",
                "done",
                {"result": AgentResult("success", "inner", {})},
            )

        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent, "execute", _fake_execute
        )
        jid = _new_job()
        C._run_agent(
            jid,
            "run qc",
            {"path": "/tmp/edis"},
            {"provider": "offline"},
            {"step_params": {"load": {"period_min": 0.01}}},
        )
        assert captured["cfg"]["step_params"] == {
            "load": {"period_min": 0.01}
        }


class TestEdiPathResolution:
    def test_selected_lines_filters_edi_path(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "qc"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod, "classify_workflow", lambda text, default=None: "qc"
        )
        captured = {}

        def _fake_execute(self, input_data):
            captured["data_path"] = input_data.get("data_path")
            return AgentResult(
                "success",
                "done",
                {"result": AgentResult("success", "inner", {})},
            )

        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent, "execute", _fake_execute
        )
        jid = _new_job()
        C._run_agent(
            jid,
            "run qc",
            {
                "path": "/tmp/all_edis",
                "groups": {"L1": ["a.edi", "b.edi"], "L2": ["c.edi"]},
                "selected_lines": ["L1"],
            },
            {"provider": "offline"},
        )
        assert captured["data_path"] == ["a.edi", "b.edi"]

    def test_line_registry_yaml_fallback(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "qc"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod, "classify_workflow", lambda text, default=None: "qc"
        )
        captured = {}

        def _fake_execute(self, input_data):
            captured["data_path"] = input_data.get("data_path")
            return AgentResult(
                "success",
                "done",
                {"result": AgentResult("success", "inner", {})},
            )

        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent, "execute", _fake_execute
        )
        # avoid the project-registry / session fallbacks interfering
        monkeypatch.setattr(C, "_session", lambda: None)
        jid = _new_job()
        C._run_agent(
            jid,
            "run qc on l22plt",
            {},
            {
                "provider": "offline",
                "line_registry": "L22PLT: /data/L22PLT\n",
            },
        )
        assert captured["data_path"] == "/data/L22PLT"


class TestEdiRequiredGuard:
    def test_missing_edi_path_reports_error(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "qc"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod, "classify_workflow", lambda text, default=None: "qc"
        )
        monkeypatch.setattr(C, "_session", lambda: None)
        jid = _new_job()
        C._run_agent(jid, "run qc", {}, {"provider": "offline"})
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ERROR
        assert "No EDI data loaded" in job["result"]


class TestPlotAndToolDispatchRouting:
    def test_plot_workflow_routes_to_dispatch_plot(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "rhophi"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod, "classify_workflow", lambda text, default=None: "rhophi"
        )
        called = {}
        monkeypatch.setattr(
            C,
            "_dispatch_plot",
            lambda jid, edi_path, *, kind, params, step, label="": called.update(
                kind=kind, edi_path=edi_path
            ),
        )
        jid = _new_job()
        C._run_agent(
            jid,
            "plot rho phi curves",
            {"path": "/tmp/edis"},
            {"provider": "offline"},
        )
        assert called["kind"] == "rhophi"
        assert called["edi_path"] == "/tmp/edis"

    def test_tool_workflow_routes_to_dispatch_tool(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "strike"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod, "classify_workflow", lambda text, default=None: "strike"
        )
        called = {}
        monkeypatch.setattr(
            C,
            "_dispatch_tool",
            lambda jid, edi_path, *, kind, params, step, label="": called.update(
                kind=kind, params=params
            ),
        )
        jid = _new_job()
        C._run_agent(
            jid,
            "run strike analysis",
            {"path": "/tmp/edis"},
            {"provider": "offline"},
        )
        assert called["kind"] == "strike"

    def test_correction_workflow_sets_corr_wf_and_fn_name(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        corr_wf = next(iter(C._CORR_METHODS))
        _patch_context(monkeypatch, config={"workflow": corr_wf})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod,
            "classify_workflow",
            lambda text, default=None: corr_wf,
        )
        called = {}
        monkeypatch.setattr(
            C,
            "_dispatch_tool",
            lambda jid, edi_path, *, kind, params, step, label="": called.update(
                kind=kind, params=params
            ),
        )
        jid = _new_job()
        C._run_agent(
            jid,
            "apply the correction",
            {"path": "/tmp/edis"},
            {"provider": "offline"},
        )
        assert called["kind"] == "correction"
        assert called["params"]["corr_wf"] == corr_wf
        assert called["params"]["fn_name"] == C._CORR_METHODS[corr_wf]["fn"]


class TestOrchestratorResultHandling:
    def test_success_collects_figures_and_code(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "qc"})
        import matplotlib.pyplot as plt

        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod, "classify_workflow", lambda text, default=None: "qc"
        )

        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        step_result = AgentResult(
            "success",
            "ok",
            {"figures": {"qc_summary": fig}},
        )
        code_result = AgentResult(
            "success", "ok", {"code": "print('hi')"}
        )
        outer = AgentResult(
            "success",
            "All good",
            {
                "result": AgentResult(
                    "success",
                    "inner",
                    {"qc": step_result, "code_gen": code_result},
                )
            },
        )
        _patch_orchestrator(monkeypatch, outer)
        jid = _new_job()
        C._run_agent(
            jid, "run qc", {"path": "/tmp/edis"}, {"provider": "offline"}
        )
        job = C._get_job(jid)
        assert job["status"] == "done"
        assert job["kind"] == C.KIND_WORKFLOW
        assert len(job["figs"]) == 1
        assert job["code"] == "print('hi')"

    def test_failed_orchestrator_result_marks_error(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "qc"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod, "classify_workflow", lambda text, default=None: "qc"
        )
        _patch_orchestrator(
            monkeypatch, AgentResult("failed", "", {}, error="boom")
        )
        jid = _new_job()
        C._run_agent(
            jid, "run qc", {"path": "/tmp/edis"}, {"provider": "offline"}
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ERROR
        assert job["result"] == "boom"

    def test_correction_workflow_caches_corrected_sites(self, monkeypatch):
        _patch_router(monkeypatch, WORKFLOW, workflow=None)
        _patch_context(monkeypatch, config={"workflow": "static_shift"})
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod,
            "classify_workflow",
            lambda text, default=None: "static_shift",
        )
        sentinel_sites = object()
        step_result = AgentResult(
            "success", "ok", {"corrected_sites": sentinel_sites}
        )
        outer = AgentResult(
            "success",
            "Corrected.",
            {
                "result": AgentResult(
                    "success",
                    "inner",
                    {"static_shift": step_result},
                )
            },
        )
        _patch_orchestrator(monkeypatch, outer)
        jid = _new_job()
        C._run_agent(
            jid,
            "run static shift",
            {"path": "/tmp/edis"},
            {"provider": "offline"},
        )
        assert C._CORR_CACHE[jid] is sentinel_sites
        job = C._get_job(jid)
        assert job["postproc"]["workflow"] == "static_shift"


class TestRunAgentExceptionHandling:
    def test_unexpected_exception_marks_job_error(self, monkeypatch):
        def boom(self, text, history=None):
            raise RuntimeError("router exploded")

        monkeypatch.setattr(router_mod.IntentRouter, "route", boom)
        jid = _new_job()
        C._run_agent(jid, "run qc", {}, {"provider": "offline"})
        job = C._get_job(jid)
        assert job["status"] == "error"
        assert job["kind"] == C.KIND_ERROR
        assert "router exploded" in job["error"]
