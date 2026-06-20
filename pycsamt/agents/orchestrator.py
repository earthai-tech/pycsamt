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
"""

from __future__ import annotations

import time
from typing import Any

from ._base import AgentResult, BaseAgent
from ..api.agents import AGENT_CONFIG

_SYSTEM_PROMPT = """\
You are a pycsamt MT workflow routing expert.
Given a natural-language MT processing request, return a JSON object with:
{
  "workflow_type": one of: "qc", "phase_analysis", "pre_inversion",
                           "ai_inversion", "inv2d", "inv3d",
                           "ensemble_inversion", "joint_inversion",
                           "modem", "full",
  "reasoning": one sentence explaining the choice
}

Rules:
- "qc" if the request mentions quality, cleaning, flagging, or noise removal only.
- "phase_analysis" if phase tensor, strike, dimensionality, or Mohr circles are mentioned.
- "pre_inversion" if Occam2D, 2D inversion, or mesh preparation is mentioned.
- "ai_inversion" if neural network, CNN, 1D deep learning, or 1D AI inversion is mentioned.
- "inv3d" if GCN, graph convolutional, 3-D AI, or spatial graph inversion is mentioned.
- "inv2d" if U-Net, 2-D AI, 2-D neural, or profile deep learning inversion is mentioned.
- "ensemble_inversion" if ensemble, uncertainty, confidence interval, or calibrated inversion is mentioned.
- "joint_inversion" if joint, multi-modal, multi-physics, TEM+MT, or combined modality is mentioned.
- "modem" if ModEM, 3D inversion, or 3D model is mentioned.
- "full" if "complete", "full pipeline", or multiple methods are mentioned.
Default to "qc" when uncertain.
Return ONLY the JSON.
"""

# ── keyword-based fallback ────────────────────────────────────────────────────

_WORKFLOW_KEYWORDS = {
    "full":               ["full", "complete", "everything", "all steps", "all pipeline"],
    "joint_inversion":    ["joint", "joint inversion", "multi-modal", "multi-physics",
                           "tem+mt", "mt+tem", "combined modality", "secondary modality"],
    "ensemble_inversion": ["ensemble", "uncertainty quantification", "confidence interval",
                           "calibrated inversion", "bayesian", "ensemble inverter"],
    "inv3d":              ["gcn", "graph convolutional", "3d ai", "3d gcn",
                           "spatial inversion", "3d neural", "inv3d",
                           "graph network inversion"],
    "inv2d":              ["unet", "u-net", "2d ai", "2d neural", "2d deep",
                           "profile inversion", "lateral continuity ai", "inv2d"],
    "ai_inversion":       ["neural", "cnn", "deep learning", "ai inversion",
                           "machine learning", "inverter", "inv1d"],
    "modem":              ["modem", "3d inversion", "3-d inversion", "3d model"],
    "pre_inversion":      ["occam", "occam2d", "2d inversion", "mesh", "startup",
                           "pre-inversion", "pre inversion", "inversion prep"],
    "phase_analysis":     ["phase tensor", "pt analysis", "strike", "dimensionality",
                           "mohr", "argand", "fingerprint", "skew analysis"],
    "qc":                 ["qc", "quality", "clean", "flag", "noise", "snr", "static shift"],
    "tipper":             ["tipper", "induction arrow", "wiese", "parkinson", "tippers"],
    "sensitivity":        ["sensitivity", "depth of investigation", "doi", "resolution kernel",
                           "bostick", "vertical resolution"],
    "freq_decimation":    ["decimate", "period selection", "frequency selection",
                           "dead band", "frequency decimation"],
    "batch":              ["batch", "multiple profiles", "all profiles", "survey batch",
                           "parallel processing"],
    "rotation":           ["rotate", "rotation", "strike rotation", "tensor rotation",
                           "coordinate frame", "principal axis"],
    "comparison":         ["compare", "comparison", "versus", "before after",
                           "difference section", "model comparison"],
}

_WORKFLOW_STEPS = {
    "qc": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files and scan per-station quality"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Frequency confidence assessment and QC flags"),
        ("static_shift",  "StaticShiftAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "Static-shift detection and AMA correction"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate QC report"),
    ],
    "phase_analysis": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("static_shift",  "StaticShiftAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "Static-shift correction"),
        ("phase_analysis","PhaseAnalysisAgent",
         "lambda r: {'sites': r['static_shift']['corrected_sites']}",
         "Phase tensor, strike, and dimensionality analysis"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate survey report"),
    ],
    "pre_inversion": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("static_shift",  "StaticShiftAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "Static-shift correction"),
        ("phase_analysis","PhaseAnalysisAgent",
         "lambda r: {'sites': r['static_shift']['corrected_sites']}",
         "Phase tensor and strike analysis"),
        ("occam2d",       "Occam2DAgent",
         "lambda r: {'sites': r['static_shift']['corrected_sites']}",
         "Write Occam2D data + mesh + startup files"),
        ("code_gen",      "CodeGenerationAgent",
         "lambda r: {'workflow_config': {}, 'results': r}",
         "Generate reproducible Python script"),
    ],
    "ai_inversion": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("denoise",       "DenoisingAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "RPCA denoising + optional AI denoising"),
        ("ai_inv",        "AIInversionAgent",
         "lambda r: {'sites': r['denoise']['denoised_sites']}",
         "AI 1-D inversion (EMInverter1D)"),
        ("interpret",     "InterpretationAgent",
         "lambda r: {'model': r['ai_inv'].get('best_model')}",
         "Geological interpretation of predicted models"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate AI inversion report"),
    ],
    "modem": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("static_shift",  "StaticShiftAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "Static-shift correction"),
        ("modem",         "ModEmAgent",
         "lambda r: {'sites': r['static_shift']['corrected_sites']}",
         "Write ModEM3D data file"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate report"),
    ],
    "inv3d": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("static_shift",  "StaticShiftAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "Static-shift correction"),
        ("inv3d",         "Inv3DAgent",
         "lambda r: {'sites': r['static_shift']['corrected_sites']}",
         "GCN 3-D spatial AI inversion"),
        ("interpret",     "InterpretationAgent",
         "lambda r: {'model': {}}",
         "Geological interpretation of 3-D volume"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate 3-D inversion report"),
    ],
    "inv2d": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("denoise",       "DenoisingAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "RPCA denoising"),
        ("inv2d",         "Inv2DAgent",
         "lambda r: {'sites': r['denoise']['denoised_sites']}",
         "U-Net 2-D profile AI inversion"),
        ("interpret",     "InterpretationAgent",
         "lambda r: {'model': {'resistivity': r['inv2d'].get('pred_section', []).tolist() "
         "if hasattr(r['inv2d'].get('pred_section', []), 'tolist') else []}}",
         "Geological interpretation of 2-D section"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate 2-D inversion report"),
    ],
    "ensemble_inversion": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("denoise",       "DenoisingAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "RPCA denoising"),
        ("ensemble",      "EnsembleAgent",
         "lambda r: {'sites': r['denoise']['denoised_sites']}",
         "Ensemble 1-D inversion with uncertainty bands"),
        ("interpret",     "InterpretationAgent",
         "lambda r: {'model': r['ensemble'].get('best_model', {})}",
         "Geological interpretation with uncertainty"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate ensemble inversion + uncertainty report"),
    ],
    "joint_inversion": [
        ("load",          "MTLoaderAgent",      None,
         "Load primary MT EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("static_shift",  "StaticShiftAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "Static-shift correction"),
        ("joint",         "JointInversionAgent",
         "lambda r: {'sites': r['static_shift']['corrected_sites']}",
         "DRCNN multi-modal joint inversion"),
        ("interpret",     "InterpretationAgent",
         "lambda r: {'model': {}}",
         "Geological interpretation of joint section"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate joint inversion report"),
    ],
    "tipper": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("tipper",        "TipperAnalysisAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Tipper analysis — induction arrows and amplitude maps"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate tipper report"),
    ],
    "sensitivity": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("sensitivity",   "SensitivityAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "Bostick sensitivity kernels and DOI analysis"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate sensitivity report"),
    ],
    "rotation": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("phase_analysis","PhaseAnalysisAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "Strike estimation from phase tensor"),
        ("rotate",        "TensorRotationAgent",
         "lambda r: {'sites': r['qc']['sites'], "
         "'strike_deg': r['phase_analysis'].get('strike_consensus', 0.0)}",
         "Rotate tensors and write corrected EDIs"),
    ],
    "full": [
        ("load",          "MTLoaderAgent",      None,
         "Load EDI files"),
        ("qc",            "DataQCAgent",
         "lambda r: {'sites': r['load']['sites']}",
         "Data quality control"),
        ("static_shift",  "StaticShiftAgent",
         "lambda r: {'sites': r['qc']['sites']}",
         "Static-shift correction"),
        ("phase_analysis","PhaseAnalysisAgent",
         "lambda r: {'sites': r['static_shift']['corrected_sites']}",
         "Phase tensor and strike analysis"),
        ("denoise",       "DenoisingAgent",
         "lambda r: {'sites': r['static_shift']['corrected_sites']}",
         "RPCA denoising"),
        ("ai_inv",        "AIInversionAgent",
         "lambda r: {'sites': r['denoise']['denoised_sites']}",
         "AI 1-D inversion"),
        ("occam2d",       "Occam2DAgent",
         "lambda r: {'sites': r['static_shift']['corrected_sites']}",
         "Write Occam2D input files"),
        ("report",        "ReportAgent",
         "lambda r: {'results': r}",
         "Generate full survey report"),
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

        if workflow_type not in _WORKFLOW_STEPS:
            candidate = self.default_workflow
            if candidate not in _WORKFLOW_STEPS:
                candidate = "qc"
            warnings.append(
                f"Workflow {workflow_type!r} is not recognised; "
                f"falling back to {candidate!r}."
            )
            workflow_type = candidate

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

        # Temporarily push the orchestrator's resolved LLM config into
        # AGENT_CONFIG so every sub-agent in the registry inherits it
        # naturally via AGENT_CONFIG.resolve() — no explicit propagation needed.
        with AGENT_CONFIG.using(
            provider=self.llm_provider if self.api_key else None,
            api_key=self.api_key,
            model=self.model if self.api_key else None,
        ):
            agent_registry = _build_registry()

        step_meta = []
        for step_name, agent_class_name, input_fn_str, desc in steps_spec:
            agent_obj = agent_registry.get(agent_class_name)
            if agent_obj is None:
                warnings.append(f"Agent {agent_class_name!r} not available; skipping.")
                continue

            # compile input_fn from string (safe eval in constrained scope)
            input_fn = None
            if input_fn_str:
                try:
                    input_fn = eval(input_fn_str)  # noqa: PGH001
                except Exception as exc:
                    warnings.append(f"input_fn for {step_name!r}: {exc}")

            coord.add_step(
                step_name, agent_obj,
                input_fn=input_fn,
                description=desc,
            )
            step_meta.append({
                "name": step_name,
                "agent": agent_class_name,
                "description": desc,
            })

        # ── run (or preview) ──────────────────────────────────────────────────
        exec_config = {
            "path":       data_path,
            "output_dir": output_dir,
            "request":    request,
        }
        exec_result = coord.execute(exec_config, dry_run=dry_run)

        elapsed = time.time() - t0
        return AgentResult(
            status=exec_result.status,
            summary=(
                f"Orchestrator routed to {workflow_type!r} workflow "
                f"({len(step_meta)} steps). {exec_result.summary}"
            ),
            data={
                "workflow_type": workflow_type,
                "reasoning":     reasoning,
                "coordinator":   coord,
                "result":        exec_result,
                "steps":         step_meta,
            },
            warnings=warnings + exec_result.warnings,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost + exec_result.cost_estimate_usd,
        )


# ── helpers ───────────────────────────────────────────────────────────────────

def _keyword_classify(text: str) -> tuple[str, str]:
    """Rule-based workflow classification from *text*."""
    t = text.lower()
    for wf in (
        "full", "joint_inversion", "ensemble_inversion", "inv3d", "inv2d",
        "ai_inversion", "modem", "pre_inversion", "phase_analysis",
        "tipper", "sensitivity", "rotation", "freq_decimation",
        "batch", "comparison", "qc",
    ):
        if any(kw in t for kw in _WORKFLOW_KEYWORDS.get(wf, [])):
            return wf, f"Matched keywords for workflow {wf!r}."
    return "qc", "No specific keywords matched; defaulted to QC."


def _build_registry() -> dict[str, Any]:
    """Instantiate all known agent classes.

    Each agent resolves its LLM config via :data:`~pycsamt.api.agents.AGENT_CONFIG`
    automatically — callers push the desired config via ``AGENT_CONFIG.using()``
    before calling this function.
    """
    registry: dict[str, Any] = {}

    def _try(name: str, factory):
        try:
            registry[name] = factory()
        except Exception as exc:
            import logging
            logging.getLogger("pycsamt.agents.orchestrator").warning(
                "Could not instantiate %s: %s", name, exc,
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
    _try("AIInversionAgent",          lambda: _import("ai_inversion",      "AIInversionAgent")())
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
    _try("BatchSurveyAgent",          lambda: _import("batch_survey",          "BatchSurveyAgent")(**common))
    _try("InversionBackendAgent",     lambda: _import("inversion_backend",     "InversionBackendAgent")(**common))

    return registry


def _import(module: str, cls: str) -> type:
    import importlib
    mod = importlib.import_module(f".{module}", package="pycsamt.agents")
    return getattr(mod, cls)


__all__ = ["WorkflowOrchestratorAgent"]
