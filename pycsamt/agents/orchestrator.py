"""
pycsamt.agents.orchestrator
============================

:class:`WorkflowOrchestratorAgent` — Intelligent workflow dispatcher.

Given a natural-language MT processing request the agent:

1. Uses the LLM (or a rule-based fallback) to classify the workflow type.
2. Builds an :class:`~pycsamt.agents.coordinator.AgentCoordinator` with the
   appropriate agent chain.
3. Executes (or previews) the full workflow, returning a consolidated result.

Supported workflow types
------------------------
``"qc"``
    QC-only pipeline: load → QC → static-shift → report

``"phase_analysis"``
    Tensor analysis: load → QC → static-shift → phase analysis → report

``"pre_inversion"``
    Pre-inversion: load → QC → static-shift → phase analysis → Occam2D prep → codegen

``"ai_inversion"``
    AI 1-D end-to-end: load → QC → denoising → AI inversion → interpretation → report

``"inv3d"``
    GCN 3-D spatial: load → QC → static-shift → Inv3DAgent → interpretation → report

``"inv2d"``
    U-Net 2-D: load → QC → denoising → Inv2DAgent → interpretation → report

``"ensemble_inversion"``
    Ensemble uncertainty: load → QC → denoising → EnsembleAgent → interpretation → report

``"joint_inversion"``
    Multi-modal DRCNN: load → QC → static-shift → JointInversionAgent → interpretation → report

``"modem"``
    ModEM prep: load → QC → static-shift → ModEM file → report

``"full"``
    Everything: load → QC → static-shift → phase analysis →
    denoising → AI inversion → Occam2D prep → report

``"denoise"``
    Denoising only: load → QC → denoising → report

``"static_shift"``
    Static-shift only: load → QC → static-shift → report

``"inversion_eval"``
    Evaluate result: load → QC → static-shift → eval → report

``"code_gen"``
    Code generation: load → QC → code generation → report

``"report"``
    Report only: load → QC → report

``"forward"``
    Forward modelling: load → QC → forward model → report

``"interpretation"``
    Interpretation: load → QC → static-shift →
    geological interpretation → report
"""

from __future__ import annotations

import json
import os
import platform
import sys
import time
from importlib.metadata import version as _pkg_version
from typing import Any

from ._base import AgentResult, BaseAgent
from ..api.agents import AGENT_CONFIG

_SYSTEM_PROMPT = """\
You are a pycsamt MT workflow routing expert.
Given a natural-language MT processing request,
return a JSON object with:
{
  "workflow_type": one of:
    "qc", "phase_analysis", "pre_inversion",
    "ai_inversion", "inv2d", "inv3d",
    "ensemble_inversion", "joint_inversion",
    "modem", "full", "full_ai_workflow",
    "pinn_inversion", "hybrid_inversion",
    "tipper", "sensitivity", "rotation",
    "freq_decimation", "batch", "comparison",
    "denoise", "static_shift", "inversion_eval",
    "code_gen", "report", "forward",
    "interpretation",
  "reasoning": one sentence explaining the choice
}

Rules:
- "qc" if quality, cleaning, or flagging only.
- "denoise" if denoising or noise removal only.
- "static_shift" if static-shift or galvanic
  distortion only (no full inversion).
- "phase_analysis" if phase tensor, strike,
  dimensionality, or Mohr circles.
- "pre_inversion" if Occam2D, mesh preparation.
- "inversion_eval" if evaluating an existing
  inversion result, RMS, misfit, or residuals.
- "ai_inversion" if 1-D AI or neural inversion.
- "inv3d" if GCN or 3-D AI inversion.
- "inv2d" if U-Net or 2-D profile AI inversion.
- "ensemble_inversion" if ensemble or uncertainty.
- "joint_inversion" if multi-modal or TEM+MT.
- "modem" if ModEM or 3-D conventional inversion.
- "pinn_inversion" if PINN or physics-informed.
- "hybrid_inversion" if two-stage or AI warm-start.
- "tipper" if tipper or induction arrows.
- "sensitivity" if sensitivity kernels or DOI.
- "rotation" if tensor rotation or strike frame.
- "freq_decimation" if period selection/decimation.
- "batch" if batch processing multiple profiles.
- "comparison" if comparing inversion results.
- "code_gen" if generating a Python script.
- "report" if generating a survey report only.
- "forward" if forward modelling or synthetic data.
- "interpretation" if geological interpretation.
- "full" if complete pipeline or multiple methods.
Default to "qc" when uncertain.
Return ONLY the JSON.
"""

# ── keyword-based fallback ────────────────────────────────────────────────────
# Workflow keyword table is the single source of truth in
# pycsamt.agents._workflows (shared with ContextInputAgent).  Kept under the
# old private name so any external references keep working.
from ._workflows import (  # noqa: E402
    WORKFLOW_KEYWORDS as _WORKFLOW_KEYWORDS,
    classify_workflow as _classify_workflow,
)

# ── input-function builders (replaces unsafe eval) ────────────────────────────
# Each function receives the running results dict and returns the
# input dict for the next agent step.  Keyed by (workflow, step_name).

def _require(r: dict, step: str, key: str) -> Any:
    """Return r[step][key], raising a clear RuntimeError
    if the step result is absent (agent was skipped)."""
    if step not in r:
        raise RuntimeError(
            f"Step {step!r} result is missing — the "
            f"responsible agent may have failed to "
            f"initialise. Check orchestrator logs."
        )
    return r[step][key]


def _ifn_load_sites(r):
    return {"sites": _require(r, "load", "sites")}

def _ifn_qc_sites(r):
    return {"sites": _require(r, "qc", "sites")}

def _ifn_ss_corrected(r):
    return {
        "sites": _require(
            r, "static_shift", "corrected_sites"
        )
    }

def _ifn_denoise_sites(r):
    return {
        "sites": _require(
            r, "denoise", "denoised_sites"
        )
    }

def _ifn_pa_corrected(r):
    return {
        "sites": _require(
            r, "static_shift", "corrected_sites"
        )
    }

def _ifn_ai_inv_model(r):
    ai_r = r["ai_inv"]
    bm = ai_r.get("best_model") or {}
    if not bm.get("resistivity"):
        # best_model empty — try first prediction
        preds = ai_r.get("predictions") or {}
        if preds:
            import numpy as _np
            nm, log_r = next(iter(preds.items()))
            n = len(log_r)
            ths = [
                10.0 * (1.5 ** i)
                for i in range(n - 1)
            ]
            bm = {
                "station":     nm,
                "resistivity": [
                    float(10.0 ** v)
                    for v in log_r
                ],
                "thickness": ths,
            }
    return {"model": bm or None}

def _ifn_ensemble_model(r):
    return {
        "model": r["ensemble"].get("best_model", {})
    }

def _ifn_empty_model(r):
    return {"model": {}}

def _ifn_inv2d_model(r):
    sec = r["inv2d"].get("pred_section", [])
    return {
        "model": {
            "resistivity": (
                sec.tolist()
                if hasattr(sec, "tolist")
                else sec
            )
        }
    }

def _ifn_results(r):
    return {"results": r}

def _ifn_codegen(r):
    return {"workflow_config": {}, "results": r}

def _ifn_rotate(r):
    return {
        "sites":      _require(r, "qc", "sites"),
        "strike_deg": r["phase_analysis"].get(
            "strike_consensus", 0.0
        ),
    }

def _ifn_pinn_from_qc(r):
    return {"sites": _require(r, "qc", "sites")}

def _ifn_pinn_model(r):
    return {
        "model": r["pinn_inv"].get("section", {})
    }

def _ifn_hybrid_model(r):
    return {
        "model": r["hybrid_inv"].get("section", {})
    }


_WORKFLOW_STEPS = {
    "qc": [
        ("load",         "MTLoaderAgent",     None,
         "Load EDI files and scan per-station quality"),
        ("qc",           "DataQCAgent",       _ifn_load_sites,
         "Frequency confidence assessment and QC flags"),
        ("static_shift", "StaticShiftAgent",  _ifn_qc_sites,
         "Static-shift detection and AMA correction"),
        ("report",       "ReportAgent",       _ifn_results,
         "Generate QC report"),
    ],
    "phase_analysis": [
        ("load",          "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("static_shift",  "StaticShiftAgent", _ifn_qc_sites,
         "Static-shift correction"),
        ("phase_analysis","PhaseAnalysisAgent",_ifn_ss_corrected,
         "Phase tensor, strike, and dimensionality analysis"),
        ("report",        "ReportAgent",      _ifn_results,
         "Generate survey report"),
    ],
    "pre_inversion": [
        ("load",          "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("static_shift",  "StaticShiftAgent", _ifn_qc_sites,
         "Static-shift correction"),
        ("phase_analysis","PhaseAnalysisAgent",_ifn_ss_corrected,
         "Phase tensor and strike analysis"),
        ("occam2d",       "Occam2DAgent",     _ifn_ss_corrected,
         "Write Occam2D data + mesh + startup files"),
        ("code_gen",      "CodeGenerationAgent",_ifn_codegen,
         "Generate reproducible Python script"),
    ],
    "ai_inversion": [
        ("load",      "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",        "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("denoise",   "DenoisingAgent",   _ifn_qc_sites,
         "RPCA denoising + optional AI denoising"),
        ("ai_inv",    "AIInversionAgent", _ifn_denoise_sites,
         "AI 1-D inversion (EMInverter1D)"),
        ("interpret", "InterpretationAgent",
         _ifn_ai_inv_model,
         "Geological interpretation of predicted models",
         False),  # optional: model may be empty
        ("report",    "ReportAgent",      _ifn_results,
         "Generate AI inversion report"),
    ],
    "inv1d": [  # alias for ai_inversion
        ("load",      "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",        "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("denoise",   "DenoisingAgent",   _ifn_qc_sites,
         "RPCA denoising"),
        ("ai_inv",    "AIInversionAgent", _ifn_denoise_sites,
         "AI 1-D inversion (EMInverter1D)"),
        ("interpret", "InterpretationAgent",
         _ifn_ai_inv_model,
         "Geological interpretation",
         False),
        ("report",    "ReportAgent",      _ifn_results,
         "Generate AI inversion report"),
    ],
    "modem": [
        ("load",         "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",           "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("static_shift", "StaticShiftAgent", _ifn_qc_sites,
         "Static-shift correction"),
        ("modem",        "ModEmAgent",       _ifn_ss_corrected,
         "Write ModEM3D data file"),
        ("report",       "ReportAgent",      _ifn_results,
         "Generate report"),
    ],
    "inv3d": [
        ("load",         "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",           "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("static_shift", "StaticShiftAgent", _ifn_qc_sites,
         "Static-shift correction"),
        ("inv3d",        "Inv3DAgent",       _ifn_ss_corrected,
         "GCN 3-D spatial AI inversion"),
        ("interpret",    "InterpretationAgent",
         _ifn_empty_model,
         "Geological interpretation of 3-D volume",
         False),
        ("report",       "ReportAgent",      _ifn_results,
         "Generate 3-D inversion report"),
    ],
    "inv2d": [
        ("load",      "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",        "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("denoise",   "DenoisingAgent",   _ifn_qc_sites,
         "RPCA denoising"),
        ("inv2d",     "Inv2DAgent",       _ifn_denoise_sites,
         "U-Net 2-D profile AI inversion"),
        ("interpret", "InterpretationAgent",
         _ifn_inv2d_model,
         "Geological interpretation of 2-D section",
         False),
        ("report",    "ReportAgent",      _ifn_results,
         "Generate 2-D inversion report"),
    ],
    "ensemble_inversion": [
        ("load",      "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",        "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("denoise",   "DenoisingAgent",   _ifn_qc_sites,
         "RPCA denoising"),
        ("ensemble",  "EnsembleAgent",    _ifn_denoise_sites,
         "Ensemble 1-D inversion with uncertainty bands"),
        ("interpret", "InterpretationAgent",
         _ifn_ensemble_model,
         "Geological interpretation with uncertainty",
         False),
        ("report",    "ReportAgent",      _ifn_results,
         "Generate ensemble inversion + uncertainty report"),
    ],
    "joint_inversion": [
        ("load",         "MTLoaderAgent",    None,
         "Load primary MT EDI files"),
        ("qc",           "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("static_shift", "StaticShiftAgent", _ifn_qc_sites,
         "Static-shift correction"),
        ("joint",        "JointInversionAgent",
         _ifn_ss_corrected,
         "DRCNN multi-modal joint inversion"),
        ("interpret",    "InterpretationAgent",
         _ifn_empty_model,
         "Geological interpretation of joint section",
         False),
        ("report",       "ReportAgent",      _ifn_results,
         "Generate joint inversion report"),
    ],
    "tipper": [
        ("load",   "MTLoaderAgent",      None,
         "Load EDI files"),
        ("tipper", "TipperAnalysisAgent",_ifn_load_sites,
         "Tipper analysis — induction arrows and amplitude maps"),
        ("report", "ReportAgent",        _ifn_results,
         "Generate tipper report"),
    ],
    "sensitivity": [
        ("load",        "MTLoaderAgent",   None,
         "Load EDI files"),
        ("qc",          "DataQCAgent",     _ifn_load_sites,
         "Data quality control"),
        ("sensitivity", "SensitivityAgent",_ifn_qc_sites,
         "Bostick sensitivity kernels and DOI analysis"),
        ("report",      "ReportAgent",     _ifn_results,
         "Generate sensitivity report"),
    ],
    "rotation": [
        ("load",          "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("phase_analysis","PhaseAnalysisAgent",_ifn_qc_sites,
         "Strike estimation from phase tensor"),
        ("rotate",        "TensorRotationAgent", _ifn_rotate,
         "Rotate tensors and write corrected EDIs"),
    ],
    "freq_decimation": [
        ("load",       "MTLoaderAgent",          None,
         "Load EDI files"),
        ("qc",         "DataQCAgent",            _ifn_load_sites,
         "Data quality control"),
        ("decimate",   "FrequencyDecimationAgent",_ifn_qc_sites,
         "Frequency decimation / period selection"),
        ("report",     "ReportAgent",            _ifn_results,
         "Generate decimation report"),
    ],
    "batch": [
        ("load",   "MTLoaderAgent",  None,
         "Load all survey EDIs in batch"),
        ("batch",  "BatchSurveyAgent",_ifn_load_sites,
         "Batch processing across profiles"),
        ("report", "ReportAgent",    _ifn_results,
         "Generate batch report"),
    ],
    "comparison": [
        ("load",     "MTLoaderAgent",           None,
         "Load EDI files"),
        ("qc",       "DataQCAgent",             _ifn_load_sites,
         "Data quality control"),
        ("compare",  "InversionComparisonAgent",_ifn_qc_sites,
         "Compare inversion results"),
        ("report",   "ReportAgent",             _ifn_results,
         "Generate comparison report"),
    ],
    # ── standalone single-agent workflows ────────────
    "denoise": [
        ("load",   "MTLoaderAgent",  None,
         "Load EDI files"),
        ("qc",     "DataQCAgent",    _ifn_load_sites,
         "Data quality control"),
        ("denoise","DenoisingAgent", _ifn_qc_sites,
         "RPCA denoising and noise filtering"),
        ("report", "ReportAgent",    _ifn_results,
         "Generate denoising report"),
    ],
    "static_shift": [
        ("load",         "MTLoaderAgent",   None,
         "Load EDI files"),
        ("qc",           "DataQCAgent",     _ifn_load_sites,
         "Data quality control"),
        ("static_shift", "StaticShiftAgent",_ifn_qc_sites,
         "Static-shift detection and AMA correction"),
        ("report",       "ReportAgent",     _ifn_results,
         "Generate static-shift report"),
    ],
    "inversion_eval": [
        ("load",         "MTLoaderAgent",          None,
         "Load EDI files"),
        ("qc",           "DataQCAgent",            _ifn_load_sites,
         "Data quality control"),
        ("static_shift", "StaticShiftAgent",       _ifn_qc_sites,
         "Static-shift correction"),
        ("inversion_eval","InversionEvaluationAgent",
         _ifn_ss_corrected,
         "Evaluate inversion: RMS, misfit, residuals"),
        ("report",       "ReportAgent",            _ifn_results,
         "Generate evaluation report"),
    ],
    "code_gen": [
        ("load",    "MTLoaderAgent",       None,
         "Load EDI files"),
        ("qc",      "DataQCAgent",         _ifn_load_sites,
         "Data quality control"),
        ("code_gen","CodeGenerationAgent", _ifn_codegen,
         "Generate reproducible Python script"),
        ("report",  "ReportAgent",         _ifn_results,
         "Generate report"),
    ],
    "report": [
        ("load",   "MTLoaderAgent", None,
         "Load EDI files"),
        ("qc",     "DataQCAgent",   _ifn_load_sites,
         "Data quality control"),
        ("report", "ReportAgent",   _ifn_results,
         "Generate survey report"),
    ],
    "forward": [
        ("load",    "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",      "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("forward", "ForwardModelAgent",_ifn_qc_sites,
         "MT forward modelling and synthetic response"),
        ("report",  "ReportAgent",      _ifn_results,
         "Generate forward model report"),
    ],
    "interpretation": [
        ("load",         "MTLoaderAgent",       None,
         "Load EDI files"),
        ("qc",           "DataQCAgent",         _ifn_load_sites,
         "Data quality control"),
        ("static_shift", "StaticShiftAgent",    _ifn_qc_sites,
         "Static-shift correction"),
        ("interpret",    "InterpretationAgent", _ifn_empty_model,
         "Geological and lithological interpretation",
         False),
        ("report",       "ReportAgent",         _ifn_results,
         "Generate interpretation report"),
    ],
    "full": [
        ("load",          "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("static_shift",  "StaticShiftAgent", _ifn_qc_sites,
         "Static-shift correction"),
        ("phase_analysis","PhaseAnalysisAgent",_ifn_ss_corrected,
         "Phase tensor and strike analysis"),
        ("denoise",       "DenoisingAgent",   _ifn_ss_corrected,
         "RPCA denoising"),
        ("ai_inv",        "AIInversionAgent", _ifn_denoise_sites,
         "AI 1-D inversion"),
        ("occam2d",       "Occam2DAgent",     _ifn_ss_corrected,
         "Write Occam2D input files"),
        ("report",        "ReportAgent",      _ifn_results,
         "Generate full survey report"),
    ],
    "full_ai_workflow": [
        ("load",          "MTLoaderAgent",    None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",      _ifn_load_sites,
         "Data quality control"),
        ("static_shift",  "StaticShiftAgent", _ifn_qc_sites,
         "Static-shift correction"),
        ("phase_analysis","PhaseAnalysisAgent",_ifn_ss_corrected,
         "Phase tensor and dimensionality"),
        ("denoise",       "DenoisingAgent",   _ifn_ss_corrected,
         "RPCA denoising"),
        ("ai_inv",        "AIInversionAgent", _ifn_denoise_sites,
         "AI 1-D inversion (EMInverter1D)"),
        ("inv2d",         "Inv2DAgent",       _ifn_denoise_sites,
         "U-Net 2-D profile inversion"),
        ("code_gen",      "CodeGenerationAgent", _ifn_codegen,
         "Generate reproducible Python script"),
        ("report",        "ReportAgent",      _ifn_results,
         "Generate full AI workflow report"),
    ],
    "pinn_inversion": [
        ("load",     "MTLoaderAgent",        None,
         "Load EDI files"),
        ("qc",       "DataQCAgent",          _ifn_load_sites,
         "Data quality control"),
        ("pinn_inv", "PINNInversionAgent",   _ifn_pinn_from_qc,
         "Physics-informed MT inversion"),
        ("interpret","InterpretationAgent",  _ifn_pinn_model,
         "Geological interpretation"),
        ("report",   "ReportAgent",          _ifn_results,
         "Generate PINN inversion report"),
    ],
    "hybrid_inversion": [
        ("load",       "MTLoaderAgent",          None,
         "Load EDI files"),
        ("qc",         "DataQCAgent",            _ifn_load_sites,
         "Data quality control"),
        ("hybrid_inv", "HybridInversionAgent",   _ifn_pinn_from_qc,
         "Two-stage AI + physics inversion"),
        ("interpret",  "InterpretationAgent",    _ifn_hybrid_model,
         "Geological interpretation"),
        ("report",     "ReportAgent",            _ifn_results,
         "Generate hybrid inversion report"),
    ],
}


class WorkflowOrchestratorAgent(BaseAgent):
    """Intelligently route an NL request to the correct agent chain.

    Parameters
    ----------
    api_key, model, llm_provider : str
    default_workflow : str
        Fallback when NL classification is ambiguous (default ``"qc"``).

    Input keys
    ----------
    ``request`` : str — natural-language processing request
    ``config`` : dict, optional — pre-built config (skips NL parsing)
    ``dry_run`` : bool — preview without executing (default False)
    ``output_dir`` : str
    ``data_path`` : str — EDI path (overrides extracted path)

    Output data keys
    ----------------
    ``workflow_type``   str
    ``reasoning``       str
    ``coordinator``     AgentCoordinator instance
    ``result``          AgentResult from the coordinator
    ``steps``           list of step metadata

    Examples
    --------
    Dry-run preview::

        agent = WorkflowOrchestratorAgent()
        r = agent.execute({
            "request": "Load L22PLT EDIs, run full phase tensor analysis",
            "dry_run": True,
        })
        print(r["workflow_type"])  # "phase_analysis"

    Full run with LLM::

        agent = WorkflowOrchestratorAgent(api_key="sk-ant-…")
        r = agent.execute({
            "request": "Clean and denoise the WILLY data, then run AI inversion",
            "data_path": "/data/WILLY_DATA",
        })
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        default_workflow: str = "qc",
    ) -> None:
        super().__init__(
            "WorkflowOrchestratorAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
        )
        self.default_workflow = default_workflow

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        request    = str(input_data.get("request", ""))
        config     = input_data.get("config") or {}
        dry_run    = bool(input_data.get("dry_run", False))
        output_dir = input_data.get("output_dir", "pycsamt_workflow_output")
        data_path  = input_data.get("data_path") or config.get("data_path", "")

        # ── classify workflow ─────────────────────────────────────────────────
        workflow_type = config.get("workflow", "")
        reasoning = ""

        if not workflow_type and request:
            if self.api_key:
                raw = self.query_llm(request, max_tokens=120)
                parsed = self.extract_json(raw or "")
                if parsed and isinstance(parsed, dict):
                    workflow_type = str(parsed.get("workflow_type", ""))
                    reasoning     = str(parsed.get("reasoning", ""))

            if not workflow_type:
                workflow_type, reasoning = _keyword_classify(request)

        available = sorted(_WORKFLOW_STEPS)
        if workflow_type not in _WORKFLOW_STEPS:
            if self.default_workflow in _WORKFLOW_STEPS:
                return AgentResult.failed(
                    f"Workflow {workflow_type!r} is recognised by "
                    f"the classifier but has no registered step "
                    f"chain.  Available workflows: "
                    f"{available}",
                    elapsed=time.time() - t0,
                )
            return AgentResult.failed(
                f"Unknown workflow {workflow_type!r}.  "
                f"Available: {available}",
                elapsed=time.time() - t0,
            )

        steps_spec = _WORKFLOW_STEPS[workflow_type]

        # ── if no data_path yet, try to extract from request ─────────────────
        if not data_path and request:
            from .context import _regex_extract
            extracted = _regex_extract(request)
            data_path = extracted.get("data_path", "")

        # ── build coordinator ─────────────────────────────────────────────────
        from .coordinator import AgentCoordinator
        coord = AgentCoordinator(
            f"orchestrated_{workflow_type}",
            checkpoint_dir=f"{output_dir}/.checkpoints",
            verbose=True,
        )

        # Extract per-workflow constructor params
        pinn_init = _extract_inv_init(
            config.get("pinn_params", {}),
            ["dim", "n_layers", "depth_max",
             "smoothness_weight",
             "lateral_weight", "graph_weight",
             "radius", "epochs", "lr",
             "solver", "comp"],
        )
        hybrid_init = _extract_inv_init(
            config.get(
                "hybrid_params",
                config.get("pinn_params", {}),
            ),
            ["dim", "max_iter",
             "smoothness_weight",
             "lateral_weight", "graph_weight",
             "radius", "lr", "solver",
             "comp", "n_freqs"],
        )
        ai_inv_init = _extract_inv_init(
            config.get("ai_inv_params", {}),
            ["n_layers", "depth_max",
             "epochs", "lr"],
        )
        checkpoint = config.get(
            "checkpoint", ""
        )

        # Push resolved LLM config into AGENT_CONFIG so every
        # sub-agent in the registry inherits it automatically.
        with AGENT_CONFIG.using(
            provider=self.llm_provider if self.api_key else None,
            api_key=self.api_key,
            model=self.model if self.api_key else None,
        ):
            agent_registry, reg_failures = _build_registry(
                pinn_init=pinn_init,
                hybrid_init=hybrid_init,
                ai_inv_init=ai_inv_init,
            )

        # Per-step params from the smart modal
        _step_params_cfg = config.get(
            "step_params", {}
        )

        step_meta = []
        for _item in steps_spec:
            step_name      = _item[0]
            agent_class_name = _item[1]
            input_fn       = _item[2]
            desc           = _item[3]
            step_required  = (
                _item[4] if len(_item) > 4 else True
            )
            agent_obj = agent_registry.get(agent_class_name)
            if agent_obj is None:
                # Root steps (input_fn=None) have no
                # upstream dependency — every downstream
                # step depends on their output.  If a root
                # step agent is missing, abort immediately
                # with a diagnostic rather than letting a
                # cryptic KeyError surface later.
                if input_fn is None:
                    fail_reason = reg_failures.get(
                        agent_class_name,
                        "unknown import error",
                    )
                    return AgentResult.failed(
                        f"Agent {agent_class_name!r} "
                        f"could not be loaded and is "
                        f"required as the first step of "
                        f"workflow {workflow_type!r}. "
                        f"Reason: {fail_reason}",
                        elapsed=time.time() - t0,
                    )
                warnings.append(
                    f"Agent {agent_class_name!r} not "
                    f"available; skipping step."
                )
                continue

            # Inject checkpoint into any
            # AI/PINN/Hybrid inversion agent
            # step via a closure over input_fn.
            _CKPT_AGENTS = {
                "HybridInversionAgent",
                "AIInversionAgent",
                "PINNInversionAgent",
                "Inv2DAgent",
                "Inv3DAgent",
            }
            step_fn = input_fn
            if (
                agent_class_name in _CKPT_AGENTS
                and checkpoint
            ):
                def _make_ckpt_fn(fn, ckpt):
                    def _wrapped(r):
                        base = fn(r) if fn else {}
                        base.setdefault(
                            "checkpoint", ckpt
                        )
                        return base
                    return _wrapped
                step_fn = _make_ckpt_fn(
                    input_fn, checkpoint
                )

            # Inject per-step params from modal
            # for non-root steps only (root steps
            # receive exec_config directly below).
            _s_extra = _step_params_cfg.get(
                step_name
            )
            if _s_extra and step_fn is not None:
                def _make_step_injector(fn, ex):
                    def _injected(r):
                        base = fn(r)
                        base.update(ex)
                        return base
                    return _injected
                step_fn = _make_step_injector(
                    step_fn, _s_extra
                )

            coord.add_step(
                step_name, agent_obj,
                input_fn=step_fn,
                description=desc,
                required=step_required,
            )
            step_meta.append({
                "name": step_name,
                "agent": agent_class_name,
                "description": desc,
            })

        # ── build and validate WorkflowPlan ──────────────────────
        from ._workflow_plan import WorkflowPlan
        plan_config = dict(config)
        plan_config.setdefault("workflow", workflow_type)
        plan_config.setdefault("data_path", data_path)
        plan_config.setdefault("output_dir", output_dir)
        plan = WorkflowPlan.from_config(
            plan_config,
            request=request,
            provider=(
                self.llm_provider if self.api_key else "offline"
            ),
        )
        warnings.extend(plan.risk_flags)

        # ── run (or preview) ──────────────────────────────────────
        exec_config = {
            "path":       data_path,
            "output_dir": output_dir,
            "request":    request,
        }
        # Merge "load" step params into exec_config
        # so the root load agent sees period range,
        # component, etc. without an input_fn.
        _load_p = _step_params_cfg.get("load", {})
        if _load_p:
            exec_config.update(_load_p)
        exec_result = coord.execute(
            exec_config, dry_run=dry_run
        )

        elapsed = time.time() - t0

        # ── write provenance trail ────────────────────────────────
        if not dry_run:
            _write_provenance(
                output_dir=output_dir,
                plan=plan,
                step_meta=step_meta,
                exec_result=exec_result,
                elapsed=elapsed,
                all_warnings=warnings + exec_result.warnings,
            )

        return AgentResult(
            status=exec_result.status,
            summary=(
                f"Orchestrator routed to {workflow_type!r} "
                f"workflow ({len(step_meta)} steps). "
                f"{exec_result.summary}"
            ),
            data={
                "workflow_type":  workflow_type,
                "reasoning":      reasoning,
                "workflow_plan":  plan,
                "coordinator":    coord,
                "result":         exec_result,
                "steps":          step_meta,
            },
            warnings=warnings + exec_result.warnings,
            elapsed_seconds=elapsed,
            cost_estimate_usd=(
                self._last_cost + exec_result.cost_estimate_usd
            ),
        )


# ── helpers ───────────────────────────────────────────────────────────────────

def _write_provenance(
    output_dir: str,
    plan: Any,
    step_meta: list,
    exec_result: Any,
    elapsed: float,
    all_warnings: list,
) -> None:
    """
    Write provenance artefacts to *output_dir*.

    Produces:
    * ``workflow_plan.json``   — validated :class:`WorkflowPlan`
    * ``agent_trace.json``     — execution trace
    * ``environment.json``     — Python / package versions
    """
    try:
        os.makedirs(output_dir, exist_ok=True)

        # workflow_plan.json
        plan_path = os.path.join(output_dir, "workflow_plan.json")
        plan.save(plan_path)

        # environment.json
        env: dict[str, Any] = {
            "python": sys.version,
            "platform": platform.platform(),
            "packages": {},
        }
        for pkg in ("pycsamt", "numpy", "scipy",
                    "matplotlib", "torch", "tensorflow"):
            try:
                env["packages"][pkg] = _pkg_version(pkg)
            except Exception:
                env["packages"][pkg] = "not installed"
        env_path = os.path.join(output_dir, "environment.json")
        with open(env_path, "w", encoding="utf-8") as fh:
            json.dump(env, fh, indent=2)

        # agent_trace.json
        trace: dict[str, Any] = {
            "user_request":     plan.request,
            "parsed_workflow":  plan.workflow_type,
            "provider":         plan.provider,
            "validation_status": (
                "passed" if plan.is_valid() else "failed"
            ),
            "executed_agents":  [s["agent"] for s in step_meta],
            "steps":            step_meta,
            "parameters":       plan.parameters,
            "warnings":         all_warnings,
            "human_review_required": plan.requires_human_review,
            "expected_outputs": plan.expected_outputs,
            "elapsed_seconds":  round(elapsed, 3),
            "exec_status":      exec_result.status,
        }
        trace_path = os.path.join(output_dir, "agent_trace.json")
        with open(trace_path, "w", encoding="utf-8") as fh:
            json.dump(trace, fh, indent=2, default=str)

    except Exception as exc:
        import logging
        logging.getLogger("pycsamt.agents.orchestrator").warning(
            "Could not write provenance: %s", exc,
        )


def _keyword_classify(text: str) -> tuple[str, str]:
    """Rule-based workflow classification from *text*.

    Delegates to the shared, ordered keyword table in
    :mod:`pycsamt.agents._workflows` so the orchestrator and
    :class:`ContextInputAgent` can never drift apart.
    """
    wf = _classify_workflow(text, default=None)
    if wf:
        return wf, f"Matched keywords for workflow {wf!r}."
    return "qc", "No specific keywords matched; defaulted to QC."


def _build_registry(
    pinn_init: "dict | None" = None,
    hybrid_init: "dict | None" = None,
    ai_inv_init: "dict | None" = None,
) -> "tuple[dict[str, Any], dict[str, str]]":
    r"""Instantiate all known agent classes.

    Each agent resolves its LLM config via
    :data:`~pycsamt.api.agents.AGENT_CONFIG` automatically.
    Call inside an ``AGENT_CONFIG.using()`` context.

    Parameters
    ----------
    pinn_init : dict, optional
        Extra keyword arguments forwarded to
        :class:`~pycsamt.agents.PINNInversionAgent`.
    hybrid_init : dict, optional
        Extra keyword arguments forwarded to
        :class:`~pycsamt.agents.HybridInversionAgent`.

    Returns
    -------
    registry : dict[str, agent]
    failures : dict[str, str]
        Agents that could not be instantiated, mapped
        to their exception message.
    """
    registry: dict[str, Any] = {}
    failures: dict[str, str] = {}
    _pi = pinn_init or {}
    _hi = hybrid_init or {}
    _ai = ai_inv_init or {}

    def _try(name: str, factory):
        try:
            registry[name] = factory()
        except Exception as exc:
            import logging
            failures[name] = str(exc)
            logging.getLogger(
                "pycsamt.agents.orchestrator"
            ).warning(
                "Could not instantiate %s: %s",
                name, exc,
            )

    _try("ContextInputAgent",         lambda: _import("context",          "ContextInputAgent")())
    _try("MTLoaderAgent",             lambda: _import("loader",            "MTLoaderAgent")())
    _try("DataQCAgent",               lambda: _import("qc",                "DataQCAgent")())
    _try("StaticShiftAgent",          lambda: _import("static_shift",      "StaticShiftAgent")())
    _try("PhaseAnalysisAgent",        lambda: _import("phase_analysis",    "PhaseAnalysisAgent")())
    _try("ForwardModelAgent",         lambda: _import("forward",           "ForwardModelAgent")())
    _try("InversionPrepAgent",        lambda: _import("inversion_prep",    "InversionPrepAgent")())
    _try("InversionEvaluationAgent",  lambda: _import("inversion_eval",    "InversionEvaluationAgent")())
    _try("InterpretationAgent",       lambda: _import("interpretation",    "InterpretationAgent")())
    _try("ReportAgent",               lambda: _import("report",            "ReportAgent")())
    _try("CodeGenerationAgent",       lambda: _import("code_gen",          "CodeGenerationAgent")())
    _try("DenoisingAgent",            lambda: _import("denoising",         "DenoisingAgent")())
    _try("AIInversionAgent",
         lambda: _import(
             "ai_inversion", "AIInversionAgent"
         )(
             n_layers=int(_ai.get("n_layers", 5)),
             epochs=int(_ai.get("epochs", 30)),
         ))
    _try("Occam2DAgent",              lambda: _import("occam2d_agent",     "Occam2DAgent")())
    _try("ModEmAgent",                lambda: _import("modem_agent",       "ModEmAgent")())
    _try("AnomalyDetectionAgent",     lambda: _import("anomaly_agent",     "AnomalyDetectionAgent")())
    _try("Inv3DAgent",                lambda: _import("inv3d_agent",       "Inv3DAgent")())
    _try("Inv2DAgent",                lambda: _import("inv2d_agent",       "Inv2DAgent")())
    _try("EnsembleAgent",             lambda: _import("ensemble_agent",    "EnsembleAgent")())
    _try("JointInversionAgent",       lambda: _import("joint_agent",       "JointInversionAgent")())
    _try("ModelZooAgent",             lambda: _import("model_zoo_agent",   "ModelZooAgent")())
    _try("TensorRotationAgent",       lambda: _import("tensor_rotation",   "TensorRotationAgent")())
    _try("EDIExportAgent",            lambda: _import("edi_export",        "EDIExportAgent")())
    _try("TipperAnalysisAgent",       lambda: _import("tipper_analysis",   "TipperAnalysisAgent")())
    _try("SensitivityAgent",          lambda: _import("sensitivity",       "SensitivityAgent")())
    _try("FrequencyDecimationAgent",  lambda: _import("freq_decimation",   "FrequencyDecimationAgent")())
    _try("InversionComparisonAgent",  lambda: _import("inversion_comparison", "InversionComparisonAgent")())
    _try("ResistivityMapAgent",       lambda: _import("resistivity_map",   "ResistivityMapAgent")())
    _try("BatchSurveyAgent",
         lambda: _import("batch_survey",      "BatchSurveyAgent")())
    _try("InversionBackendAgent",
         lambda: _import("inversion_backend", "InversionBackendAgent")())
    _try("PINNInversionAgent",
         lambda: _import("pinn_agent",   "PINNInversionAgent")(**_pi))
    _try("HybridInversionAgent",
         lambda: _import("hybrid_agent", "HybridInversionAgent")(**_hi))

    return registry, failures


def _extract_inv_init(
    params: dict, allowed: list
) -> dict:
    """Return only the *allowed* keys from *params*."""
    return {k: v for k, v in params.items() if k in allowed}


def _import(module: str, cls: str) -> type:
    import importlib
    mod = importlib.import_module(
        f".{module}", package="pycsamt.agents"
    )
    return getattr(mod, cls)


__all__ = ["WorkflowOrchestratorAgent"]
