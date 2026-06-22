# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents._workflows
=========================

**Single source of truth** for workflow identifiers, their human-readable
descriptions, the natural-language keyword → workflow mapping, and the alias
table.

Before this module the keyword table lived in *two* places — once in
:mod:`pycsamt.agents.context` (used by :class:`ContextInputAgent`) and again,
slightly drifted, in :mod:`pycsamt.agents.orchestrator` (used by
:func:`_keyword_classify`).  They had to be hand-synced, and any divergence
made classification inconsistent depending on which path ran first.

Everything related to *which workflow a phrase maps to* now lives here, and
both consumers import :func:`classify_workflow` / :data:`WORKFLOW_KEYWORDS`
from this module.

Ordering matters
----------------
:data:`WORKFLOW_KEYWORDS` is an **ordered** mapping: classification is
first-match-wins, so the most specific phrases must come first (e.g. ``"gcn"``
before bare ``"inversion"``).  Do not reorder casually — the routing
regression suite (``test_nl_routing.py``) is validated against this order.
"""

from __future__ import annotations

__all__ = [
    "WORKFLOW_DESCRIPTIONS",
    "WORKFLOW_KEYWORDS",
    "WORKFLOW_ALIASES",
    "classify_workflow",
    "normalise_workflow",
]


# ── canonical workflow descriptions ─────────────────────────────────────────────

WORKFLOW_DESCRIPTIONS: dict[str, str] = {
    # --- core processing ---
    "qc":                 "Data quality control and preprocessing",
    "static_shift":       "Static-shift detection and correction",
    "phase_analysis":     "Phase tensor, strike, dimensionality analysis",
    "forward":            "Forward modelling (1D / 2D / 3D)",
    "inversion_prep":     "Prepare Occam2D / ModEM inversion files",
    "pre_inversion":      "Prepare Occam2D / ModEM inversion files",
    "inversion_eval":     "Evaluate inversion result quality",
    "interpretation":     "Geological interpretation",
    "report":             "Generate survey report",
    "full":               "Full pipeline: QC to report",
    # --- AI / DL inversion ---
    "ai_inversion":       "AI 1-D inversion (EMInverter1D / CNN)",
    "inv1d":              "AI 1-D inversion (EMInverter1D / CNN)",
    "inv2d":              "U-Net 2-D profile AI inversion",
    "inv3d":              "GCN 3-D spatial AI inversion",
    "ensemble_inversion": "Ensemble 1-D inversion with uncertainty",
    "joint_inversion":    "Multi-modal DRCNN joint inversion",
    "pinn_inversion":     "Physics-informed (PINN) inversion",
    "hybrid_inversion":   "Two-stage AI warm-start + physics inversion",
    # --- external codes ---
    "modem":              "ModEM 3-D inversion file preparation",
    "occam2d":            "Occam2D 2-D inversion file preparation",
    # --- geophysical analysis ---
    "tipper":             "Tipper / induction arrow analysis",
    "sensitivity":        "Bostick sensitivity and DOI analysis",
    "rotation":           "Tensor rotation to principal axes",
    "freq_decimation":    "Frequency decimation / period selection",
    # --- survey-level ---
    "batch":              "Batch processing of multiple profiles",
    "comparison":         "Compare inversion results",
    "full_ai_workflow":   "Full AI-assisted pipeline",
    # --- agent-focused ---
    "code_gen":           "Generate a reproducible Python script",
    "denoise":            "Data denoising",
}


# ── natural-language keyword mapping (ordered, most-specific first) ──────────────
# First-match wins.  Narrow phrases precede bare keywords so that, e.g.,
# "gcn" does not get swallowed by "inversion".

WORKFLOW_KEYWORDS: dict[str, list[str]] = {
    "full_ai_workflow": [
        "full ai workflow",
        "full ai pipeline",
        "ai full pipeline",
    ],
    "full": [
        "full pipeline",
        "full workflow",
        "end to end",
        "all steps",
    ],
    "joint_inversion": [
        "joint inversion",
        "multi-modal",
        "multi-physics",
        "tem+mt", "mt+tem",
        "combined modality",
    ],
    "ensemble_inversion": [
        "ensemble inversion",
        "ensemble method",
        "ensemble",
        "uncertainty quantification",
        "calibrated inversion",
        "confidence interval",
        "bayesian",
    ],
    "inv3d": [
        "gcn",
        "graph convolutional",
        "3d ai", "3d gcn",
        "spatial inversion",
        "3d neural", "graph network",
        "inv3d",
    ],
    "inv2d": [
        "unet", "u-net",
        "2d ai", "2d neural", "2d deep",
        "profile inversion",
        "lateral continuity ai",
        "inv2d",
    ],
    "pinn_inversion": [
        "pinn inversion", "pinn",
        "physics-informed",
        "physics informed",
        "no training data",
        "gradient descent inversion",
    ],
    "hybrid_inversion": [
        "hybrid inversion", "two-stage",
        "two stage", "ai + physics",
        "warm start", "ai warmstart",
        "stage 1 stage 2",
    ],
    "comparison": [
        "compare",
        "comparison",
        "versus",
        "before after",
    ],
    # inversion_eval, modem, pre_inversion all before ai_inversion so specific
    # keywords ("modem", "occam2d", "inversion result") are not swallowed by
    # bare "inversion".
    "inversion_eval": [
        "evaluate inversion",
        "inversion result",
        "inversion quality",
        "check inversion",
        "rms", "misfit",
        "residual pt",
    ],
    "modem": [
        "modem",
        "3d inversion",
        "3d model",
    ],
    "pre_inversion": [
        "occam2d", "occam",
        "2d inversion", "mesh",
        "startup", "pre-inversion",
        "pre inversion", "inversion prep",
        "prepare inversion",
        "inversion data file",
        "inversion data",
        "data file for inversion",
    ],
    "tipper": [
        "tipper", "induction arrow",
        "wiese", "parkinson",
    ],
    "sensitivity": [
        "sensitivity",
        "depth of investigation",
        "bostick",
    ],
    "rotation": [
        "rotate", "rotation",
        "strike rotation",
        "tensor rotation",
    ],
    # freq_decimation before ai_inversion so "decimate/period range" never
    # routes to the bare "inversion" catch-all below.
    "freq_decimation": [
        "decimate", "decimation",
        "period selection",
        "frequency selection",
        "select period",
        "select frequency",
        "period range",
        "frequency range",
        "dead band",
    ],
    "batch": [
        "batch process",
        "batch mode",
        "batch run",
        "multiple profiles",
        "all profiles",
        "survey batch",
    ],
    # Action-verb workflows before content keywords so "write code for X"
    # routes to code_gen, not X's workflow.
    "code_gen": [
        "write code", "generate code",
        "python script", "code for",
        "script for", "write a script",
        "write function", "write class",
        "create notebook",
        "notebook for",
    ],
    "report": [
        "generate report",
        "write report",
        "survey report",
        "create report",
        "make report",
    ],
    "denoise": [
        "denoise", "denoising",
        "remove noise",
        "noise reduction",
        "filter noise",
    ],
    # Mostly multi-word phrases to avoid false positives on bare "mohr"
    # / "argand". Bare "strike" and "skew" are included: in MT they
    # almost always mean geoelectric strike / phase-tensor skew, and
    # rotation (checked earlier) already owns "strike rotation" /
    # "rotate", so rotation requests are not mis-routed here.
    "phase_analysis": [
        "phase tensor",
        "phase analysis",
        "pt analysis",
        "phase tensor analysis",
        "mohr circle",
        "argand diagram",
        "strike analysis",
        "geoelectric strike",
        "strike direction",
        "strike angle",
        "strike",
        "dimensionality analysis",
        "dimensionality",
        "bahr skew",
        "skew",
    ],
    "static_shift": [
        "static shift",
        "static-shift",
        "galvanic distortion",
        "galvanic",
    ],
    "forward": [
        "forward model",
        "forward modeling",
        "synthetic data",
        "simulate",
    ],
    # Catch bare "inversion/invert" only after all specific variants have been
    # checked.
    "ai_inversion": [
        "ai inversion",
        "neural", "cnn",
        "deep learning",
        "machine learning",
        "inverter", "inv1d", "1d ai",
        "1d neural", "1d inversion",
        "inversion", "invert", "inverting",
        "run inversion", "do inversion",
        "perform inversion",
        "start inversion",
    ],
    # Removed bare "geology", "resistor", "conductor" -- too broad for MT text.
    "interpretation": [
        "geological interpretation",
        "lithological interpretation",
        "geological context",
        "geological unit",
        "interpret geology",
        "interpret",
        "lithology",
        "geological",
        "lithological",
    ],
    "qc": [
        "qc", "quality control",
        "quality check",
        "data quality",
        "remove noisy",
        "remove bad data",
    ],
}


# ── alias normalisation ─────────────────────────────────────────────────────────

WORKFLOW_ALIASES: dict[str, str] = {
    "qc_preprocessing":        "qc",
    "preprocessing":           "qc",
    # galvanic distortion -> static shift
    "galvanic":                "static_shift",
    "pre-inversion":           "pre_inversion",
    "pre_inversion_prep":      "pre_inversion",
    "post-inversion":          "inversion_eval",
    "interpretation":          "interpretation",
    "report":                  "report",
    "full_pipeline":           "full",
    "1d_inversion":            "ai_inversion",
    "1d_ai":                   "ai_inversion",
    "ai_1d":                   "ai_inversion",
    "ai1d":                    "ai_inversion",
    "2d_inversion":            "inv2d",
    "ai_2d":                   "inv2d",
    "3d_inversion":            "inv3d",
    "ai_3d":                   "inv3d",
    "ensemble":                "ensemble_inversion",
    "joint":                   "joint_inversion",
    "occam":                   "pre_inversion",
    "occam2d":                 "pre_inversion",
    "full_ai":                 "full_ai_workflow",
    # PINN / Hybrid aliases
    "pinn":                    "pinn_inversion",
    "physics_informed":        "pinn_inversion",
    "physics-informed":        "pinn_inversion",
    "pinn_inv":                "pinn_inversion",
    "hybrid":                  "hybrid_inversion",
    "two_stage":               "hybrid_inversion",
    "hybrid_inv":              "hybrid_inversion",
    "ai_physics":              "hybrid_inversion",
    # Agent-focused workflow aliases
    "code_generation":         "code_gen",
    "codegen":                 "code_gen",
    "write_code":              "code_gen",
    "denoising":               "denoise",
    "noise_removal":           "denoise",
    "doi":                     "sensitivity",
    "depth_of_investigation":  "sensitivity",
    "frequency_decimation":    "freq_decimation",
    "freq_dec":                "freq_decimation",
    "period_decimation":       "freq_decimation",
}


# ── public API ──────────────────────────────────────────────────────────────────

def classify_workflow(
    text: str,
    *,
    default: str | None = None,
) -> str | None:
    """Classify free *text* into a workflow id via :data:`WORKFLOW_KEYWORDS`.

    First-match-wins over the ordered keyword table (case-insensitive).

    Parameters
    ----------
    text : str
        Natural-language request.
    default : str or None
        Returned when no keyword matches.  Callers that want the historical
        "leave unset, normalise later" behaviour pass ``None``; callers that
        want a hard default (e.g. the orchestrator) pass ``"qc"``.

    Returns
    -------
    str or None
        The matched workflow id, or *default*.

    Examples
    --------
    >>> classify_workflow("run a GCN 3D inversion")
    'inv3d'
    >>> classify_workflow("xyzzy", default="qc")
    'qc'
    >>> classify_workflow("xyzzy") is None
    True
    """
    t = (text or "").lower()
    for wf, kws in WORKFLOW_KEYWORDS.items():
        if any(kw in t for kw in kws):
            return wf
    return default


def normalise_workflow(workflow: str) -> str:
    """Normalise a workflow string: lowercase, separators, alias resolution.

    Returns ``"qc"`` for empty / falsy input.

    Examples
    --------
    >>> normalise_workflow("PINN")
    'pinn_inversion'
    >>> normalise_workflow("1d-inversion")
    'ai_inversion'
    >>> normalise_workflow("")
    'qc'
    """
    wf = (
        str(workflow or "")
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
    )
    return WORKFLOW_ALIASES.get(wf, wf) or "qc"
