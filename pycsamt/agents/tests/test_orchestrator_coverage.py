# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.orchestrator`.

``test_orchestrator.py`` only exercises ``dry_run=True`` requests, which
never reaches step registration edge cases, plan validation, real agent
execution, or provenance writing. This file adds:

* direct unit tests for the small ``_ifn_*`` input-function builders,
  ``_require``'s error branch, and ``_make_outdir_injector``'s closure;
* direct unit tests for ``_sha256_file``, ``_infer_producing_step``,
  ``_build_output_manifest``, and a ``_build_registry`` instantiation
  failure;
* ``execute()``-level tests for LLM-classification success, the two
  unknown-workflow failure branches, request-based data_path extraction,
  a missing root-step agent, checkpoint/step-param injection (via
  ``dry_run=True`` so setup runs without real training), and the plan
  validation gate; and
* one real, non-dry-run "report" workflow execution against the bundled
  3-EDI dataset, which exercises the entire provenance-writing path
  (``_write_provenance`` and everything it calls) for real.
"""

from __future__ import annotations

import json

import pytest

from pycsamt.agents.orchestrator import (
    WorkflowOrchestratorAgent,
    _build_output_manifest,
    _build_registry,
    _ifn_ai_inv_model,
    _ifn_codegen,
    _ifn_denoise_sites,
    _ifn_empty_model,
    _ifn_ensemble_model,
    _ifn_hybrid_model,
    _ifn_inv2d_model,
    _ifn_load_sites,
    _ifn_pa_corrected,
    _ifn_pinn_from_qc,
    _ifn_pinn_model,
    _ifn_qc_sites,
    _ifn_results,
    _ifn_rotate,
    _ifn_ss_corrected,
    _infer_producing_step,
    _make_outdir_injector,
    _require,
    _sha256_file,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


# ── _require and the _ifn_* input-function builders ─────────────────────────


def test_require_missing_step_raises():
    with pytest.raises(RuntimeError, match="missing"):
        _require({}, "load", "sites")


def test_require_present_step_returns_value():
    assert _require({"load": {"sites": "SITES"}}, "load", "sites") == "SITES"


def test_ifn_load_and_qc_and_ss_and_denoise():
    r = {
        "load": {"sites": "L"},
        "qc": {"sites": "Q"},
        "static_shift": {"corrected_sites": "SS"},
        "denoise": {"denoised_sites": "D"},
    }
    assert _ifn_load_sites(r) == {"sites": "L"}
    assert _ifn_qc_sites(r) == {"sites": "Q"}
    assert _ifn_ss_corrected(r) == {"sites": "SS"}
    assert _ifn_denoise_sites(r) == {"sites": "D"}
    assert _ifn_pa_corrected(r) == {"sites": "SS"}


def test_ifn_ai_inv_model_uses_best_model_when_present():
    r = {"ai_inv": {"best_model": {"resistivity": [1.0, 2.0]}}}
    assert _ifn_ai_inv_model(r) == {"model": {"resistivity": [1.0, 2.0]}}


def test_ifn_ai_inv_model_falls_back_to_first_prediction():
    r = {
        "ai_inv": {
            "best_model": {},
            "predictions": {"s1": [1.0, 2.0, 3.0]},
        }
    }
    out = _ifn_ai_inv_model(r)
    assert out["model"]["station"] == "s1"
    assert len(out["model"]["resistivity"]) == 3
    assert len(out["model"]["thickness"]) == 2


def test_ifn_ai_inv_model_empty_when_nothing_available():
    assert _ifn_ai_inv_model({"ai_inv": {}}) == {"model": None}


def test_ifn_ensemble_and_empty_and_inv2d_model():
    assert _ifn_ensemble_model(
        {"ensemble": {"best_model": {"a": 1}}}
    ) == {"model": {"a": 1}}
    assert _ifn_empty_model({}) == {"model": {}}
    out = _ifn_inv2d_model({"inv2d": {"pred_section": [[1, 2], [3, 4]]}})
    assert out == {"model": {"resistivity": [[1, 2], [3, 4]]}}


def test_ifn_inv2d_model_converts_ndarray():
    import numpy as np

    out = _ifn_inv2d_model(
        {"inv2d": {"pred_section": np.array([[1.0, 2.0]])}}
    )
    assert out["model"]["resistivity"] == [[1.0, 2.0]]


def test_ifn_results_and_codegen():
    r = {"a": 1}
    assert _ifn_results(r) == {"results": r}
    out = _ifn_codegen(r)
    assert out == {"workflow_config": {}, "results": r}


def test_ifn_rotate():
    r = {
        "qc": {"sites": "Q"},
        "phase_analysis": {"strike_consensus": 42.0},
    }
    assert _ifn_rotate(r) == {"sites": "Q", "strike_deg": 42.0}


def test_ifn_pinn_and_hybrid():
    r = {
        "qc": {"sites": "Q"},
        "pinn_inv": {"section": {"x": 1}},
        "hybrid_inv": {"section": {"y": 2}},
    }
    assert _ifn_pinn_from_qc(r) == {"sites": "Q"}
    assert _ifn_pinn_model(r) == {"model": {"x": 1}}
    assert _ifn_hybrid_model(r) == {"model": {"y": 2}}


def test_make_outdir_injector_sets_default_but_not_override():
    injected = _make_outdir_injector(lambda r: {"sites": "S"}, "/out/a")
    assert injected({}) == {"sites": "S", "output_dir": "/out/a"}

    injected2 = _make_outdir_injector(
        lambda r: {"sites": "S", "output_dir": "/already/set"}, "/out/a"
    )
    assert injected2({})["output_dir"] == "/already/set"


# ── provenance helpers ────────────────────────────────────────────────────────


def test_sha256_file(tmp_path):
    p = tmp_path / "f.txt"
    p.write_bytes(b"hello world")
    import hashlib

    expected = hashlib.sha256(b"hello world").hexdigest()
    assert _sha256_file(p) == expected


def test_infer_producing_step_matches_prep_subfolder_or_name():
    assert _infer_producing_step("pycsamt_occam2d/data.dat", ["occam2d"]) == (
        "occam2d"
    )
    assert _infer_producing_step("qc_report.md", ["qc", "report"]) == "qc"
    assert _infer_producing_step("unrelated.txt", ["qc"]) is None


def test_build_output_manifest_empty_dir_returns_empty_list(tmp_path):
    assert _build_output_manifest(str(tmp_path / "missing"), ["qc"]) == []


def test_build_output_manifest_lists_files(tmp_path):
    (tmp_path / "qc_report.md").write_text("hi")
    (tmp_path / "workflow_plan.json").write_text("{}")  # excluded
    entries = _build_output_manifest(str(tmp_path), ["qc"])
    paths = [e["path"] for e in entries]
    assert "qc_report.md" in paths
    assert "workflow_plan.json" not in paths
    entry = next(e for e in entries if e["path"] == "qc_report.md")
    assert entry["producing_step"] == "qc"
    assert "sha256" in entry


# ── _build_registry ───────────────────────────────────────────────────────────


def test_build_registry_records_instantiation_failures(monkeypatch):
    import pycsamt.agents.orchestrator as orch

    real_import = orch._import

    def _boom_import(module, cls):
        if module == "loader":
            raise RuntimeError("boom loader")
        return real_import(module, cls)

    monkeypatch.setattr(orch, "_import", _boom_import)
    registry, failures = _build_registry()
    assert "MTLoaderAgent" in failures
    assert "boom loader" in failures["MTLoaderAgent"]
    assert "MTLoaderAgent" not in registry
    # unrelated agents still instantiate fine
    assert "DataQCAgent" in registry


# ── execute()-level branch tests ─────────────────────────────────────────────


def test_llm_classification_success(monkeypatch):
    agent = WorkflowOrchestratorAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: (
        '{"workflow_type": "qc", "reasoning": "quality check requested"}'
    )
    result = agent.execute(
        {
            "request": "please clean up /data/x",
            "dry_run": True,
        }
    )
    assert result["workflow_type"] == "qc"
    assert result["reasoning"] == "quality check requested"


def test_unknown_workflow_with_valid_default_fails(monkeypatch):
    agent = WorkflowOrchestratorAgent(default_workflow="qc")
    result = agent.execute(
        {"config": {"workflow": "not_a_real_workflow"}, "dry_run": True}
    )
    assert result.status == "failed"
    assert "recognised by the classifier" in result.error


def test_unknown_workflow_with_invalid_default_fails(monkeypatch):
    agent = WorkflowOrchestratorAgent(default_workflow="also_not_real")
    result = agent.execute(
        {"config": {"workflow": "not_a_real_workflow"}, "dry_run": True}
    )
    assert result.status == "failed"
    assert "Unknown workflow" in result.error


def test_data_path_extracted_from_request(monkeypatch):
    agent = WorkflowOrchestratorAgent()
    result = agent.execute(
        {
            "request": "run QC on /data/some_survey please",
            "config": {"workflow": "qc"},
            "dry_run": True,
        }
    )
    # dry-run preview still builds the plan, which carries the resolved
    # data_path extracted from the free-text request.
    assert result["workflow_plan"].data_path == "/data/some_survey"


def test_missing_root_agent_fails(monkeypatch):
    import pycsamt.agents.orchestrator as orch

    monkeypatch.setattr(
        orch, "_build_registry", lambda **k: ({}, {"MTLoaderAgent": "boom"})
    )
    agent = WorkflowOrchestratorAgent()
    result = agent.execute(
        {
            "request": "run QC on /data/x",
            "config": {"workflow": "qc"},
            "dry_run": True,
        }
    )
    assert result.status == "failed"
    assert "MTLoaderAgent" in result.error
    assert "could not be loaded" in result.error


def test_checkpoint_injection_during_dry_run_setup():
    agent = WorkflowOrchestratorAgent()
    result = agent.execute(
        {
            "request": "run AI inversion on /data/x",
            "config": {"workflow": "ai_inversion", "checkpoint": "ckpt.npz"},
            "dry_run": True,
        }
    )
    assert result.status in ("success", "needs_review")
    assert result["workflow_type"] == "ai_inversion"


def test_step_params_injection_during_dry_run_setup():
    agent = WorkflowOrchestratorAgent()
    result = agent.execute(
        {
            "request": "run QC on /data/x",
            "config": {
                "workflow": "qc",
                "step_params": {"qc": {"snr_threshold": 5.0}},
            },
            "dry_run": True,
        }
    )
    assert result.status in ("success", "needs_review")


def test_plan_validation_failure_blocks_execution_without_dry_run():
    # An empty 'request' fails WorkflowPlan validation; since dry_run is
    # False the validation gate must block before any agent runs.
    agent = WorkflowOrchestratorAgent()
    result = agent.execute(
        {"config": {"workflow": "report", "data_path": "data/3edis"}}
    )
    assert result.status == "failed"
    assert "failed validation" in result.summary


# ── full non-dry-run execution: exercises provenance writing for real ───────


def test_real_report_workflow_writes_provenance(tmp_output):
    agent = WorkflowOrchestratorAgent()
    result = agent.execute(
        {
            "request": "Generate a QC report for /data/3edis",
            "config": {"workflow": "report", "data_path": "data/3edis"},
            "output_dir": str(tmp_output),
        }
    )
    assert result.status == "success"

    plan_path = tmp_output / "workflow_plan.json"
    trace_path = tmp_output / "agent_trace.json"
    env_path = tmp_output / "environment.json"
    manifest_path = tmp_output / "output_manifest.json"
    for p in (plan_path, trace_path, env_path, manifest_path):
        assert p.exists(), p

    trace = json.loads(trace_path.read_text())
    assert trace["parsed_workflow"] == "report"
    assert trace["executed_agents"] == [
        "MTLoaderAgent",
        "DataQCAgent",
        "ReportAgent",
    ]
    assert trace["exec_status"] in ("success", "needs_review")

    manifest = json.loads(manifest_path.read_text())
    assert manifest["files"]
    assert any(
        f["path"].endswith("survey_report.md") for f in manifest["files"]
    )
