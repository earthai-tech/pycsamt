# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the ``_dispatch_*`` task handlers in
pycsamt.app.agent_master.callbacks.chat — the leaf functions ``_run_agent``
routes into once an intent/workflow has been classified.

These handlers are plain functions (no Dash context needed): each is called
directly with a fresh job id, its heavy agent dependency (PackageQAAgent,
PlotAgent, ToolAgent, MetricsAgent, WorkflowOrchestratorAgent,
CodeGenerationAgent) monkeypatched at the class level so the test exercises
the dispatcher's own branching (result shaping, figure collection, warning
surfacing, error paths) without running real LLM/plotting/IO code.
"""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")

import pycsamt.agents.code_gen as codegen_mod
import pycsamt.agents.context as context_mod
import pycsamt.agents.metrics as metrics_mod
import pycsamt.agents.orchestrator as orch_mod
import pycsamt.agents.package_qa as qa_mod
import pycsamt.agents.plotting as plot_mod
import pycsamt.agents.tooling as tool_mod
import pycsamt.app.agent_master.callbacks.chat as C
from pycsamt.agents._base import AgentResult


def _noop_step(label, status="done"):
    pass


def _new_job():
    return C._new_job()


# ── _dispatch_question ────────────────────────────────────────────────────


class TestDispatchQuestion:
    def test_online_answer_no_offline_nudge(self, monkeypatch):
        monkeypatch.setattr(
            qa_mod.PackageQAAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success", "s", {"answer": "It does X.", "source": "llm"}
            ),
        )
        jid = _new_job()
        C._dispatch_question(
            jid,
            "what does StaticShiftAgent do",
            llm_prov="claude",
            api_key="sk-x",
            sel_model=None,
            offline=False,
            history=None,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ANSWER
        assert job["result"] == "It does X."
        assert "Offline answer" not in job["result"]

    def test_offline_answer_gets_nudge(self, monkeypatch):
        monkeypatch.setattr(
            qa_mod.PackageQAAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success", "s", {"answer": "Offline text.", "source": "rag_offline"}
            ),
        )
        jid = _new_job()
        C._dispatch_question(
            jid,
            "what is qc",
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            history=[{"role": "user", "content": "earlier"}],
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert "Offline answer composed" in job["result"]

    def test_clarify_source_skips_nudge_even_offline(self, monkeypatch):
        monkeypatch.setattr(
            qa_mod.PackageQAAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success", "s", {"answer": "Which line?", "source": "rag_clarify"}
            ),
        )
        jid = _new_job()
        C._dispatch_question(
            jid,
            "what about it",
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            history=None,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["result"] == "Which line?"

    def test_missing_answer_falls_back_to_summary(self, monkeypatch):
        monkeypatch.setattr(
            qa_mod.PackageQAAgent,
            "execute",
            lambda self, input_data: AgentResult("success", "the summary", {}),
        )
        jid = _new_job()
        C._dispatch_question(
            jid,
            "explain qc",
            llm_prov="claude",
            api_key="k",
            sel_model=None,
            offline=False,
            history=None,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert "the summary" in job["result"]

    def test_missing_answer_and_summary_uses_default_text(self, monkeypatch):
        monkeypatch.setattr(
            qa_mod.PackageQAAgent,
            "execute",
            lambda self, input_data: AgentResult("success", "", {}),
        )
        jid = _new_job()
        C._dispatch_question(
            jid,
            "explain qc",
            llm_prov="claude",
            api_key="k",
            sel_model=None,
            offline=False,
            history=None,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert "couldn't find an answer" in job["result"]


# ── _dispatch_plot ─────────────────────────────────────────────────────────


class TestDispatchPlot:
    def test_success_collects_figure_and_records_run(self, monkeypatch):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot([0, 1], [1, 0])
        monkeypatch.setattr(
            plot_mod.PlotAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success", "Plot ready.", {"figures": {"rhophi": fig}}
            ),
        )
        jid = _new_job()
        C._dispatch_plot(
            jid,
            "/tmp/edis",
            kind="rhophi",
            params={},
            step=_noop_step,
            label="L22",
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_WORKFLOW
        assert len(job["figs"]) == 1
        assert "Plot ready." in job["result"]

    def test_no_tipper_reason_is_meta_not_error(self, monkeypatch):
        monkeypatch.setattr(
            plot_mod.PlotAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "failed", "no tipper", {"reason": "no_tipper"}
            ),
        )
        jid = _new_job()
        C._dispatch_plot(
            jid, "/tmp/edis", kind="tipper", params={}, step=_noop_step, label="L1"
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_META
        assert "Tipper data is not available for L1" in job["result"]

    def test_generic_failure_uses_hint(self, monkeypatch):
        monkeypatch.setattr(
            plot_mod.PlotAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "failed", "could not plot", {"hint": "check the path"}
            ),
        )
        jid = _new_job()
        C._dispatch_plot(
            jid, "/tmp/edis", kind="rhophi", params={}, step=_noop_step
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ERROR
        assert "Hint: check the path" in job["result"]

    def test_no_figures_produced_notes_it(self, monkeypatch):
        monkeypatch.setattr(
            plot_mod.PlotAgent,
            "execute",
            lambda self, input_data: AgentResult("success", "Done.", {}),
        )
        jid = _new_job()
        C._dispatch_plot(
            jid, "/tmp/edis", kind="rhophi", params={}, step=_noop_step
        )
        job = C._get_job(jid)
        assert "No figure was produced" in job["result"]

    def test_warnings_surfaced_in_summary(self, monkeypatch):
        monkeypatch.setattr(
            plot_mod.PlotAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success", "Done.", {}, warnings=["low coherence"]
            ),
        )
        jid = _new_job()
        C._dispatch_plot(
            jid, "/tmp/edis", kind="rhophi", params={}, step=_noop_step
        )
        job = C._get_job(jid)
        assert "low coherence" in job["result"]


# ── _dispatch_tool ─────────────────────────────────────────────────────────


class TestDispatchTool:
    def test_success_with_table_and_figure(self, monkeypatch):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        monkeypatch.setattr(
            tool_mod.ToolAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success",
                "Strike computed.",
                {"figures": {"strike": fig}, "table_text": "az | conf\n10 | 0.9"},
            ),
        )
        jid = _new_job()
        C._dispatch_tool(
            jid, "/tmp/edis", kind="strike", params={}, step=_noop_step, label="L1"
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_WORKFLOW
        assert len(job["figs"]) == 1
        assert "```" in job["result"]
        assert "az | conf" in job["result"]

    def test_success_without_figures_is_answer_kind(self, monkeypatch):
        monkeypatch.setattr(
            tool_mod.ToolAgent,
            "execute",
            lambda self, input_data: AgentResult("success", "Validated.", {}),
        )
        jid = _new_job()
        C._dispatch_tool(
            jid, "/tmp/edis", kind="validator", params={}, step=_noop_step
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ANSWER
        assert job["figs"] == {}

    def test_failure_reports_dataset_label(self, monkeypatch):
        monkeypatch.setattr(
            tool_mod.ToolAgent,
            "execute",
            lambda self, input_data: AgentResult("failed", "could not validate", {}),
        )
        jid = _new_job()
        C._dispatch_tool(
            jid,
            "/tmp/edis",
            kind="validator",
            params={},
            step=_noop_step,
            label="L1",
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ERROR
        assert "(dataset for L1)" in job["result"]

    def test_corrected_sites_trigger_postproc_modal(self, monkeypatch):
        sentinel = object()
        monkeypatch.setattr(
            tool_mod.ToolAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success", "Edited.", {"corrected_sites": sentinel}
            ),
        )
        jid = _new_job()
        C._dispatch_tool(
            jid, "/tmp/edis", kind="freq_editor", params={}, step=_noop_step
        )
        job = C._get_job(jid)
        assert C._CORR_CACHE[jid] is sentinel
        assert job["postproc"]["workflow"] == "freq_editor"


# ── _dispatch_metrics ──────────────────────────────────────────────────────


class TestDispatchMetrics:
    def test_no_targets_reports_no_data(self):
        jid = _new_job()
        C._dispatch_metrics(
            jid, "compute mean resistivity", {}, {}, step=_noop_step
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_META
        assert "don't have any survey data" in job["result"]

    def test_single_target_uses_metrics_agent_summary(self, monkeypatch):
        monkeypatch.setattr(
            metrics_mod.MetricsAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success", "Mean rho: 120 Ohm.m", {"values": {"mean_rho": 120}}
            ),
        )
        jid = _new_job()
        C._dispatch_metrics(
            jid,
            "mean resistivity",
            {"path": "/tmp/edis"},
            {},
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ANSWER
        assert "Mean rho: 120" in job["result"]

    def test_multi_target_lists_per_line_single_kind(self, monkeypatch):
        def _fake_execute(self, input_data):
            label = input_data["label"]
            val = 100 if label == "L1" else 200
            return AgentResult(
                "success", f"{label} ok", {"values": {"mean_rho": val}}
            )

        monkeypatch.setattr(metrics_mod.MetricsAgent, "execute", _fake_execute)
        monkeypatch.setattr(
            metrics_mod, "parse_metric_request", lambda text: (["mean_rho"], True)
        )
        jid = _new_job()
        C._dispatch_metrics(
            jid,
            "mean resistivity for all lines",
            {"groups": {"L1": ["a.edi"], "L2": ["b.edi"]}},
            {},
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert "**L1**: 100" in job["result"]
        assert "**L2**: 200" in job["result"]
        assert "Here's the mean_rho for each line" in job["result"]

    def test_multi_target_multi_kind_and_failure_line(self, monkeypatch):
        def _fake_execute(self, input_data):
            if input_data["label"] == "L1":
                return AgentResult("failed", "could not read L1", {})
            return AgentResult(
                "success", "ok", {"values": {"a": 1, "b": 2}}
            )

        monkeypatch.setattr(metrics_mod.MetricsAgent, "execute", _fake_execute)
        monkeypatch.setattr(
            metrics_mod, "parse_metric_request", lambda text: (["a", "b"], True)
        )
        jid = _new_job()
        C._dispatch_metrics(
            jid,
            "a and b for all lines",
            {"groups": {"L1": ["a.edi"], "L2": ["b.edi"]}},
            {},
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert "**L1**: could not read L1" in job["result"]
        assert "**L2**: a: 1; b: 2" in job["result"]
        assert "Here's what I found per line" in job["result"]

    def test_warnings_appended(self, monkeypatch):
        monkeypatch.setattr(
            metrics_mod.MetricsAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success", "ok", {"values": {}}, warnings=["missing coords"]
            ),
        )
        jid = _new_job()
        C._dispatch_metrics(
            jid, "summary", {"path": "/tmp/edis"}, {}, step=_noop_step
        )
        job = C._get_job(jid)
        assert "missing coords" in job["result"]


# ── _line_stats ────────────────────────────────────────────────────────────


class TestLineStats:
    def test_computes_freq_qc_and_flags(self, monkeypatch):
        monkeypatch.setattr(
            metrics_mod.MetricsAgent, "_scan", staticmethod(lambda sites: self._rows_static())
        )
        monkeypatch.setattr(
            metrics_mod, "_station_coords", lambda sites: []
        )
        monkeypatch.setattr(metrics_mod, "_has_tipper", lambda sites: True)
        warnings: list = []
        out = C._line_stats("L1", object(), warnings)
        assert out["n_stations"] == 2
        assert out["qc"] == pytest.approx(55.0)
        assert out["flagged"] == 1
        assert out["tipper"] is True
        assert any("1 of 2 station(s) flagged" in w for w in warnings)

    @staticmethod
    def _rows_static():
        return [
            {
                "station": "S1",
                "t_min_s": 1.0,
                "t_max_s": 100.0,
                "n_freq": 40,
                "qc_score": 90.0,
                "has_z": True,
                "has_coords": True,
            },
            {
                "station": "S2",
                "t_min_s": 1.0,
                "t_max_s": 50.0,
                "n_freq": 30,
                "qc_score": 20.0,
                "has_z": False,
                "has_coords": True,
            },
        ]

    def test_length_km_computed_from_coords(self, monkeypatch):
        monkeypatch.setattr(
            metrics_mod.MetricsAgent, "_scan", staticmethod(lambda sites: [])
        )
        monkeypatch.setattr(
            metrics_mod,
            "_station_coords",
            lambda sites: [("S1", 38.0, -118.0), ("S2", 38.01, -118.0)],
        )
        monkeypatch.setattr(
            metrics_mod,
            "_ll_to_utm",
            lambda la, lo, zone, hemi, datum: (500000.0 + la * 1000, 4200000.0, 11),
        )
        monkeypatch.setattr(metrics_mod, "_has_tipper", lambda sites: False)
        out = C._line_stats("L1", object(), [])
        assert out["length_km"] is not None
        assert out["length_km"] > 0

    def test_coord_exception_leaves_length_none(self, monkeypatch):
        monkeypatch.setattr(
            metrics_mod.MetricsAgent, "_scan", staticmethod(lambda sites: [])
        )

        def _boom(sites):
            raise RuntimeError("no coords")

        monkeypatch.setattr(metrics_mod, "_station_coords", _boom)
        monkeypatch.setattr(metrics_mod, "_has_tipper", lambda sites: False)
        out = C._line_stats("L1", object(), [])
        assert out["length_km"] is None

    def test_tipper_exception_defaults_false(self, monkeypatch):
        monkeypatch.setattr(
            metrics_mod.MetricsAgent, "_scan", staticmethod(lambda sites: [])
        )
        monkeypatch.setattr(metrics_mod, "_station_coords", lambda sites: [])

        def _boom(sites):
            raise RuntimeError("no tipper info")

        monkeypatch.setattr(metrics_mod, "_has_tipper", _boom)
        out = C._line_stats("L1", object(), [])
        assert out["tipper"] is False


# ── _dispatch_data_overview: full "read the data" path ─────────────────────


class TestDispatchDataOverviewFullPath:
    def test_single_named_line_reads_only_that_line(self, monkeypatch):
        calls = []

        def _fake_line_stats(label, sites, warnings):
            calls.append((label, sites))
            return {
                "label": label,
                "n_stations": 3,
                "stations": ["a", "b", "c"],
                "freq": (1.0, 100.0),
                "max_nfreq": 20,
                "qc": 80.0,
                "flagged": 0,
                "length_km": 1.0,
                "tipper": False,
            }

        import pycsamt.emtools._core as core_mod

        monkeypatch.setattr(core_mod, "ensure_sites", lambda *a, **k: object())
        monkeypatch.setattr(C, "_line_stats", _fake_line_stats)
        jid = _new_job()
        C._dispatch_data_overview(
            jid,
            "read line L22PLT",
            {"groups": {"L22PLT": ["a.edi"], "K1": ["b.edi"]}},
            {},
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ANSWER
        assert len(calls) == 1
        assert calls[0][0] == "L22PLT"
        assert "Read **L22PLT**" in job["result"]
        assert job["card"]["tiles"][0]["value"] == "3"

    def test_all_lines_reads_every_group(self, monkeypatch):
        def _fake_line_stats(label, sites, warnings):
            return {
                "label": label,
                "n_stations": 1,
                "stations": [label],
                "freq": None,
                "max_nfreq": None,
                "qc": None,
                "flagged": 0,
                "length_km": None,
                "tipper": False,
            }

        import pycsamt.emtools._core as core_mod

        monkeypatch.setattr(core_mod, "ensure_sites", lambda *a, **k: object())
        monkeypatch.setattr(C, "_line_stats", _fake_line_stats)
        jid = _new_job()
        C._dispatch_data_overview(
            jid,
            "read the data",
            {"groups": {"L1": ["a.edi"], "L2": ["b.edi"]}},
            {},
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert "2 lines" in job["result"] or "2 lines" in job["card"]["scope"]
        assert len(job["card"]["lines"]) == 2

    def test_read_exception_is_collected_as_warning(self, monkeypatch):
        import pycsamt.emtools._core as core_mod

        def _boom(*a, **k):
            raise RuntimeError("bad EDI")

        monkeypatch.setattr(core_mod, "ensure_sites", _boom)
        jid = _new_job()
        C._dispatch_data_overview(
            jid,
            "read the data",
            {"path": "/tmp/edis"},
            {},
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ERROR
        assert "could not read any of them" in job["result"]
        assert "bad EDI" in job["result"]


# ── _dispatch_inversion_prep ────────────────────────────────────────────────


def _prep_result(prep_step, files_map, extra_warns=None, figures=None):
    """Build the nested AgentResult shape _collect_prep_files expects."""
    step_data = dict(files_map)
    if figures:
        step_data["figures"] = figures
    step_res = AgentResult(
        "success", "step ok", step_data, warnings=extra_warns or []
    )
    inner = AgentResult("success", "inner", {prep_step: step_res})
    return AgentResult("success", "All good.", {"result": inner})


class TestDispatchInversionPrep:
    def test_no_targets_reports_no_data(self):
        jid = _new_job()
        C._dispatch_inversion_prep(
            jid,
            "modem",
            "prep modem",
            {},
            {},
            {},
            inv_config={},
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_META
        assert "prepare inversion files" in job["result"]

    def test_single_line_success_lists_files(self, monkeypatch, tmp_path):
        data_file = tmp_path / "ModEM_Data.dat"
        data_file.write_text("x" * 100, encoding="utf-8")
        result = _prep_result("modem", {"data_path": str(data_file)})
        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent, "execute", lambda self, d: result
        )
        jid = _new_job()
        C._dispatch_inversion_prep(
            jid,
            "modem",
            "prep modem inversion",
            {"path": str(tmp_path)},
            {"output_dir": str(tmp_path / "out")},
            {},
            inv_config={},
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_WORKFLOW
        assert "ModEM_Data.dat" in job["result"]
        assert "ModEM inversion preparation complete." in job["result"]

    def test_multi_line_reports_per_line_and_header_count(self, monkeypatch, tmp_path):
        ok_file = tmp_path / "ModEM_Data.dat"
        ok_file.write_text("x", encoding="utf-8")

        calls = {"n": 0}

        def _fake_execute(self, input_data):
            calls["n"] += 1
            if calls["n"] == 1:
                return _prep_result("modem", {"data_path": str(ok_file)})
            return AgentResult("failed", "no data", {}, error="no data for L2")

        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent, "execute", _fake_execute
        )
        jid = _new_job()
        C._dispatch_inversion_prep(
            jid,
            "modem",
            "prep modem for all lines",
            {"groups": {"L1": ["a.edi"], "L2": ["b.edi"]}},
            {"output_dir": str(tmp_path / "out")},
            {},
            inv_config={},
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert "1/2" in job["result"]
        assert "no data for L2" in job["result"]
        assert job["kind"] == C.KIND_WORKFLOW  # n_ok > 0

    def test_all_lines_fail_marks_error(self, monkeypatch, tmp_path):
        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent,
            "execute",
            lambda self, d: AgentResult("failed", "boom", {}, error="boom"),
        )
        jid = _new_job()
        C._dispatch_inversion_prep(
            jid,
            "pre_inversion",
            "prep occam2d",
            {"path": str(tmp_path)},
            {"output_dir": str(tmp_path / "out")},
            {},
            inv_config={},
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ERROR
        assert "preparation failed" in job["result"]

    def test_orchestrator_exception_is_caught_per_line(self, monkeypatch, tmp_path):
        def _boom(self, input_data):
            raise RuntimeError("kaboom")

        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent, "execute", _boom
        )
        jid = _new_job()
        C._dispatch_inversion_prep(
            jid,
            "mare2dem",
            "prep mare2dem",
            {"path": str(tmp_path)},
            {"output_dir": str(tmp_path / "out")},
            {},
            inv_config={},
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_ERROR
        assert "kaboom" in job["result"]

    def test_prep_params_forwarded_into_step_params(self, monkeypatch, tmp_path):
        captured = {}
        ok_file = tmp_path / "ModEM_Data.dat"
        ok_file.write_text("x", encoding="utf-8")

        def _fake_execute(self, input_data):
            captured["cfg"] = dict(input_data["config"])
            return _prep_result("modem", {"data_path": str(ok_file)})

        monkeypatch.setattr(
            orch_mod.WorkflowOrchestratorAgent, "execute", _fake_execute
        )
        jid = _new_job()
        C._dispatch_inversion_prep(
            jid,
            "modem",
            "prep modem",
            {"path": str(tmp_path)},
            {"output_dir": str(tmp_path / "out")},
            {},
            inv_config={"error_floor": 5.0},
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        assert captured["cfg"]["step_params"]["modem"]["error_floor"] == 5.0


# ── _prep_rag_note: real citation extraction path ──────────────────────────


class TestPrepRagNoteSuccess:
    def test_citations_rendered_as_symbols(self, monkeypatch):
        from pycsamt.assistant.rag.context_builder import AssembledContext

        ctx = AssembledContext(
            query="q",
            context_text="...",
            citations=[
                {"symbol": "OccamRunner.run"},
                {"source_path": "models/occam2d/data.py"},
                {"symbol": "OccamRunner.run"},  # duplicate, deduped
            ],
            chunks=["chunk-placeholder"],
        )

        class _Builder:
            def build(self, query):
                return ctx

        import pycsamt.assistant.rag.context_builder as cb_mod

        monkeypatch.setattr(cb_mod, "default_context_builder", lambda: _Builder())
        note = C._prep_rag_note("build occam2d files")
        assert "OccamRunner.run" in note
        assert "models/occam2d/data.py" in note
        assert note.count("OccamRunner.run") == 1

    def test_builder_returns_none(self, monkeypatch):
        import pycsamt.assistant.rag.context_builder as cb_mod

        monkeypatch.setattr(cb_mod, "default_context_builder", lambda: None)
        assert C._prep_rag_note("q") == ""

    def test_empty_context_returns_blank(self, monkeypatch):
        from pycsamt.assistant.rag.context_builder import AssembledContext

        class _Builder:
            def build(self, query):
                return AssembledContext(query=query)

        import pycsamt.assistant.rag.context_builder as cb_mod

        monkeypatch.setattr(cb_mod, "default_context_builder", lambda: _Builder())
        assert C._prep_rag_note("q") == ""


# ── _dispatch_code: extra branches beyond the RAG-happy-path test ─────────


class TestDispatchCodeExtraBranches:
    def _patch_context_and_codegen(self, monkeypatch, cfg=None, code="print(1)"):
        monkeypatch.setattr(
            context_mod.ContextInputAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success", "ok", {"config": cfg or {}}
            ),
        )
        monkeypatch.setattr(
            codegen_mod.CodeGenerationAgent,
            "execute",
            lambda self, input_data: AgentResult(
                "success", "ok", {"code": code}
            ),
        )

    def test_unspecified_script_asks_for_task_instead_of_defaulting_to_qc(self, monkeypatch):
        self._patch_context_and_codegen(monkeypatch, cfg={"workflow": "code_gen"})
        import pycsamt.assistant.rag.context_builder as cb_mod

        monkeypatch.setattr(cb_mod, "default_context_builder", lambda: None)
        monkeypatch.setattr(
            "pycsamt.assistant.tools.validation_tools.validate_generated_code",
            lambda code: {"ok": True, "syntax_ok": True, "errors": []},
        )
        jid = _new_job()
        C._dispatch_code(
            jid,
            "write me a script",
            {},
            {},
            workflow=None,
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_CLARIFY
        assert "What should the script" in job["result"]
        assert "code" not in job

    def test_uses_loaded_edi_path_when_no_line_resolved(self, monkeypatch):
        self._patch_context_and_codegen(monkeypatch)
        import pycsamt.assistant.rag.context_builder as cb_mod

        monkeypatch.setattr(cb_mod, "default_context_builder", lambda: None)
        monkeypatch.setattr(
            "pycsamt.assistant.tools.validation_tools.validate_generated_code",
            lambda code: {"ok": True, "syntax_ok": True, "errors": []},
        )
        captured = {}
        orig = codegen_mod.CodeGenerationAgent.execute

        def _spy(self, input_data):
            captured["cfg"] = dict(input_data["workflow_config"])
            return AgentResult("success", "ok", {"code": "print(1)"})

        monkeypatch.setattr(codegen_mod.CodeGenerationAgent, "execute", _spy)
        jid = _new_job()
        C._dispatch_code(
            jid,
            "generate code for qc",
            {"path": "/data/loaded_edis"},
            {},
            workflow="qc",
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        assert captured["cfg"]["data_path"] == "/data/loaded_edis"

    def test_syntax_error_validation_branch(self, monkeypatch):
        self._patch_context_and_codegen(monkeypatch)
        import pycsamt.assistant.rag.context_builder as cb_mod

        monkeypatch.setattr(cb_mod, "default_context_builder", lambda: None)
        monkeypatch.setattr(
            "pycsamt.assistant.tools.validation_tools.validate_generated_code",
            lambda code: {
                "ok": False,
                "syntax_ok": False,
                "errors": ["SyntaxError: bad"],
            },
        )
        jid = _new_job()
        C._dispatch_code(
            jid,
            "generate code for qc",
            {},
            {},
            workflow="qc",
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert "syntax error" in job["result"]
        assert "SyntaxError: bad" in job["result"]

    def test_unresolved_symbol_validation_branch(self, monkeypatch):
        self._patch_context_and_codegen(monkeypatch)
        import pycsamt.assistant.rag.context_builder as cb_mod

        monkeypatch.setattr(cb_mod, "default_context_builder", lambda: None)
        monkeypatch.setattr(
            "pycsamt.assistant.tools.validation_tools.validate_generated_code",
            lambda code: {
                "ok": False,
                "syntax_ok": True,
                "errors": ["Unknown symbol: pycsamt.foo.Bar"],
            },
        )
        jid = _new_job()
        C._dispatch_code(
            jid,
            "generate code for qc",
            {},
            {},
            workflow="qc",
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert "some symbols could not be verified" in job["result"]

    def test_validation_failure_is_explicitly_unverifiable(self, monkeypatch):
        self._patch_context_and_codegen(monkeypatch)
        import pycsamt.assistant.rag.context_builder as cb_mod

        monkeypatch.setattr(cb_mod, "default_context_builder", lambda: None)

        def _boom(code):
            raise RuntimeError("validator unavailable")

        monkeypatch.setattr(
            "pycsamt.assistant.tools.validation_tools.validate_generated_code",
            _boom,
        )
        jid = _new_job()
        C._dispatch_code(
            jid,
            "generate code for qc",
            {},
            {},
            workflow="qc",
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_CODE
        assert "Validated" not in job["result"]
        assert "Validation" in job["result"]
        assert "unverifiable" in job["result"]
        assert "validator unavailable" in job["result"]

    def test_rag_exception_degrades_silently(self, monkeypatch):
        self._patch_context_and_codegen(monkeypatch)
        import pycsamt.assistant.rag.context_builder as cb_mod

        def _boom():
            raise RuntimeError("rag down")

        monkeypatch.setattr(cb_mod, "default_context_builder", _boom)
        monkeypatch.setattr(
            "pycsamt.assistant.tools.validation_tools.validate_generated_code",
            lambda code: {"ok": True, "syntax_ok": True, "errors": []},
        )
        jid = _new_job()
        C._dispatch_code(
            jid,
            "generate code for qc",
            {},
            {},
            workflow="qc",
            llm_prov="claude",
            api_key=None,
            sel_model=None,
            offline=True,
            step=_noop_step,
        )
        job = C._get_job(jid)
        assert job["kind"] == C.KIND_CODE
