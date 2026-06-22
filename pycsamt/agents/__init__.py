"""
pycsamt.agents
==============

AI-assisted MT/AMT/CSAMT workflow automation.

All agents are **lazy-loaded**: importing this module costs nothing unless
you actually instantiate an agent.  This keeps the base pycsamt import fast
even when LLM libraries (``anthropic``, ``openai``, ``google-generativeai``)
are not installed.

Agent catalogue
---------------
:class:`ContextInputAgent`
    Translate a natural-language request into a structured workflow config.
:class:`MTLoaderAgent`
    Load EDI / AVG / J files into a :class:`~pycsamt.site.Sites` object
    with a per-station QC report.
:class:`AgentCoordinator`
    Chain agents into named workflows with checkpointing and cost tracking.
:class:`DataQCAgent`
    SNR section, dead-band detection, per-station quality scores.
:class:`StaticShiftAgent`
    Static-shift detection and correction (AMA / LOESS / spatial median).
:class:`PhaseAnalysisAgent`
    Phase tensor, strike, dimensionality, Mohr circles, Argand diagrams.
:class:`ForwardModelAgent`
    1-D, 2-D, and 3-D MT forward modelling.
:class:`InversionPrepAgent` / :class:`Occam2DAgent` / :class:`ModEmAgent` / :class:`Mare2DEMAgent`
    Write Occam2D, ModEM3D, or MARE2DEM inversion input files; optionally run
    the binary and scan results.
:class:`InversionEvaluationAgent`
    Compute RMS, residual PT section, misfit pseudosection.
:class:`InterpretationAgent`
    Map resistivity ranges to lithology; correlate with borehole logs.
:class:`ReportAgent`
    Assemble all figures and tables into Markdown / HTML / PDF.
:class:`CodeGenerationAgent`
    Emit a standalone reproducible Python script from the workflow config.
:class:`WorkflowOrchestratorAgent`
    NL → classify workflow → build and run the correct agent chain.
:class:`DenoisingAgent`
    RPCA / Hampel / EMAP / AI-CAE denoising.
:class:`AnomalyDetectionAgent`
    Unsupervised CAE anomaly flagging per (station, frequency).
:class:`AIInversionAgent`
    End-to-end 1-D AI inversion (ResNet / CNN / FCN).
:class:`Inv2DAgent`
    U-Net 2-D profile inversion with lateral continuity.
:class:`Inv3DAgent`
    GCN 3-D spatial inversion using inter-station graph message-passing.
:class:`EnsembleAgent`
    Ensemble 1-D inversion with conformal uncertainty bands.
:class:`JointInversionAgent`
    DRCNN multi-modal joint inversion (MT + TEM / CSAMT / gravity).
:class:`ModelZooAgent`
    Browse, download, and deploy pre-trained EM inverter checkpoints.
:class:`PINNInversionAgent`
    Physics-informed 1-D/2-D/3-D MT inversion via Adam
    (no labelled training data required).
:class:`HybridInversionAgent`
    Two-stage AI warm-start + physics refinement for
    1-D, 2-D, and 3-D MT inversion.

Supporting classes
------------------
:class:`AgentResult`
    Standardised return type for every agent.
:class:`BaseAgent`
    Abstract base class; inherit to build custom agents.

Quick start
-----------
Without an LLM key (regex fallback, no cost)::

    from pycsamt.agents import ContextInputAgent, MTLoaderAgent, AgentCoordinator

    ctx    = ContextInputAgent()   # no api_key → pure regex
    loader = MTLoaderAgent()

    coord  = AgentCoordinator("mt_qc")
    coord.add_step("parse", ctx,
                   description="Parse request into config")
    coord.add_step("load",  loader,
                   input_fn=lambda r: {"path": r["parse"]["config"]["data_path"]},
                   description="Load EDI files and QC scan")

    result = coord.execute(
        {"request": "Load /data/EDIs, QC, period range 1e-4 to 1 s"},
        dry_run=True,
    )

With a Claude API key::

    ctx = ContextInputAgent(api_key="sk-ant-…", llm_provider="claude")
    res = ctx.execute({"request": "Load L22PLT EDIs, run PT analysis …"})
    print(res["config"])
    print(res.llm_interpretation)

LLM providers supported
-----------------------
* **Anthropic Claude** (default) — ``llm_provider="claude"``
* **OpenAI**                      — ``llm_provider="openai"``
* **Google Gemini**               — ``llm_provider="gemini"``
* **DeepSeek**                    — ``llm_provider="deepseek"``

Install the relevant client library before use:
``pip install anthropic`` / ``pip install openai`` /
``pip install google-generativeai``
DeepSeek reuses the ``openai`` package with a custom
``base_url``; no extra install is needed.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

# ── always-available, zero-cost exports ──────────────────────────────────────
from ._base import AgentResult, BaseAgent
from ._pricing import (
    estimate_cost, format_cost, get_rate, PROVIDER_RATES,
)
from ._workflow_plan import (
    WorkflowPlan,
    validate_workflow_plan,
    VALID_WORKFLOWS,
)
from .coordinator import AgentCoordinator, WorkflowStep
from ..api.agents import (
    AGENT_CONFIG, AgentConfig, BudgetExceededError,
    configure_agents, reset_agents,
)

# ── lazy agent map ────────────────────────────────────────────────────────────
_LAZY: dict[str, str] = {
    "ContextInputAgent":          ".context",
    "IntentRouter":               ".router",
    "PackageQAAgent":             ".package_qa",
    "MTLoaderAgent":              ".loader",
    "DataQCAgent":                ".qc",
    "StaticShiftAgent":           ".static_shift",
    "PhaseAnalysisAgent":         ".phase_analysis",
    "ForwardModelAgent":          ".forward",
    "InversionPrepAgent":         ".inversion_prep",
    "InversionEvaluationAgent":   ".inversion_eval",
    "InterpretationAgent":        ".interpretation",
    "ReportAgent":                ".report",
    "CodeGenerationAgent":        ".code_gen",
    "WorkflowOrchestratorAgent":  ".orchestrator",
    "DenoisingAgent":             ".denoising",
    "AIInversionAgent":           ".ai_inversion",
    "Occam2DAgent":               ".occam2d_agent",
    "ModEmAgent":                 ".modem_agent",
    "AnomalyDetectionAgent":      ".anomaly_agent",
    "Inv2DAgent":                 ".inv2d_agent",
    "Inv3DAgent":                 ".inv3d_agent",
    "EnsembleAgent":              ".ensemble_agent",
    "JointInversionAgent":        ".joint_agent",
    "ModelZooAgent":              ".model_zoo_agent",
    # new agents
    "TensorRotationAgent":        ".tensor_rotation",
    "EDIExportAgent":             ".edi_export",
    "TipperAnalysisAgent":        ".tipper_analysis",
    "SensitivityAgent":           ".sensitivity",
    "FrequencyDecimationAgent":   ".freq_decimation",
    "InversionComparisonAgent":   ".inversion_comparison",
    "ResistivityMapAgent":        ".resistivity_map",
    "BatchSurveyAgent":           ".batch_survey",
    "InversionBackendAgent":      ".inversion_backend",
    "PipelineAgent":              ".pipeline_agent",
    "Mare2DEMAgent":              ".mare2dem_agent",
    "PINNInversionAgent":         ".pinn_agent",
    "HybridInversionAgent":       ".hybrid_agent",
}


def __getattr__(name: str):
    if name in _LAZY:
        module = importlib.import_module(_LAZY[name], package=__name__)
        obj = getattr(module, name)
        globals()[name] = obj   # cache — skip import machinery on next access
        return obj
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


if TYPE_CHECKING:
    from .context          import ContextInputAgent
    from .router           import IntentRouter, RouterDecision
    from .package_qa       import PackageQAAgent
    from .loader           import MTLoaderAgent
    from .qc               import DataQCAgent
    from .static_shift     import StaticShiftAgent
    from .phase_analysis   import PhaseAnalysisAgent
    from .forward          import ForwardModelAgent
    from .inversion_prep   import InversionPrepAgent
    from .inversion_eval   import InversionEvaluationAgent
    from .interpretation   import InterpretationAgent
    from .report           import ReportAgent
    from .code_gen         import CodeGenerationAgent
    from .orchestrator     import WorkflowOrchestratorAgent
    from .denoising        import DenoisingAgent
    from .ai_inversion     import AIInversionAgent
    from .occam2d_agent    import Occam2DAgent
    from .modem_agent      import ModEmAgent
    from .anomaly_agent    import AnomalyDetectionAgent
    from .inv2d_agent      import Inv2DAgent
    from .inv3d_agent      import Inv3DAgent
    from .ensemble_agent   import EnsembleAgent
    from .joint_agent      import JointInversionAgent
    from .model_zoo_agent       import ModelZooAgent
    from .tensor_rotation       import TensorRotationAgent
    from .edi_export            import EDIExportAgent
    from .tipper_analysis       import TipperAnalysisAgent
    from .sensitivity           import SensitivityAgent
    from .freq_decimation       import FrequencyDecimationAgent
    from .inversion_comparison  import InversionComparisonAgent
    from .resistivity_map       import ResistivityMapAgent
    from .batch_survey          import BatchSurveyAgent
    from .inversion_backend     import InversionBackendAgent
    from .pipeline_agent        import PipelineAgent
    from .mare2dem_agent        import Mare2DEMAgent
    from .pinn_agent            import PINNInversionAgent
    from .hybrid_agent          import HybridInversionAgent


__version__ = "2.0.0"

__all__ = [
    # global LLM config
    "AGENT_CONFIG",
    "AgentConfig",
    "BudgetExceededError",
    "configure_agents",
    "reset_agents",
    # core infrastructure
    "AgentResult",
    "BaseAgent",
    "AgentCoordinator",
    "WorkflowStep",
    # pricing helpers
    "estimate_cost",
    "format_cost",
    "get_rate",
    "PROVIDER_RATES",
    # workflow plan
    "WorkflowPlan",
    "validate_workflow_plan",
    "VALID_WORKFLOWS",
    # agents — alphabetical within logical groups
    "ContextInputAgent",
    "IntentRouter",
    "PackageQAAgent",
    "MTLoaderAgent",
    "DataQCAgent",
    "StaticShiftAgent",
    "PhaseAnalysisAgent",
    "ForwardModelAgent",
    "InversionPrepAgent",
    "InversionEvaluationAgent",
    "InterpretationAgent",
    "ReportAgent",
    "CodeGenerationAgent",
    "WorkflowOrchestratorAgent",
    "DenoisingAgent",
    "AIInversionAgent",
    "AnomalyDetectionAgent",
    "Occam2DAgent",
    "ModEmAgent",
    "Inv2DAgent",
    "Inv3DAgent",
    "EnsembleAgent",
    "JointInversionAgent",
    "ModelZooAgent",
    "TensorRotationAgent",
    "EDIExportAgent",
    "TipperAnalysisAgent",
    "SensitivityAgent",
    "FrequencyDecimationAgent",
    "InversionComparisonAgent",
    "ResistivityMapAgent",
    "BatchSurveyAgent",
    "InversionBackendAgent",
    "PipelineAgent",
    "Mare2DEMAgent",
    "PINNInversionAgent",
    "HybridInversionAgent",
]
