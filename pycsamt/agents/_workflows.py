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
    "qc": "Data quality control and preprocessing",
    "static_shift": "Static-shift detection and correction",
    "phase_analysis": "Phase tensor, strike, dimensionality analysis",
    "forward": "Forward modelling (1D / 2D / 3D)",
    "inversion_prep": "Prepare Occam2D / ModEM inversion files",
    "pre_inversion": "Prepare Occam2D / ModEM inversion files",
    "inversion_eval": "Evaluate inversion result quality",
    "interpretation": "Geological interpretation",
    "report": "Generate survey report",
    "full": "Full pipeline: QC to report",
    # --- AI / DL inversion ---
    "ai_inversion": "AI 1-D inversion (EMInverter1D / CNN)",
    "inv1d": "AI 1-D inversion (EMInverter1D / CNN)",
    "inv2d": "U-Net 2-D profile AI inversion",
    "inv3d": "GCN 3-D spatial AI inversion",
    "ensemble_inversion": "Ensemble 1-D inversion with uncertainty",
    "joint_inversion": "Multi-modal DRCNN joint inversion",
    "pinn_inversion": "Physics-informed (PINN) inversion",
    "hybrid_inversion": "Two-stage AI warm-start + physics inversion",
    # --- external codes ---
    "modem": "ModEM 3-D inversion file preparation",
    "mare2dem": "MARE2DEM 2.5-D inversion file preparation",
    "occam2d": "Occam2D 2-D inversion file preparation",
    # --- geophysical analysis ---
    "tipper": "Tipper / induction arrow analysis",
    "sensitivity": "Bostick sensitivity and DOI analysis",
    "rotation": "Tensor rotation to principal axes",
    "freq_decimation": "Frequency decimation / period selection",
    # --- survey-level ---
    "batch": "Batch processing of multiple profiles",
    "comparison": "Compare inversion results",
    "full_ai_workflow": "Full AI-assisted pipeline",
    # --- agent-focused ---
    "code_gen": "Generate a reproducible Python script",
    "denoise": "Data denoising",
    # --- plotting ---
    "rhophi": "Apparent-resistivity / phase sounding curves",
    "phase_psection": "Scalar phase pseudo-section",
    "pt_psection": "Phase-tensor (Phi) ellipse pseudo-section",
    "pt_strip": "Single-station phase-tensor ellipse strip vs period",
    "pt_strip_grid": "Phase-tensor ellipse strips tiled by survey line",
    "tipper_plot": "Tipper components / induction-arrow plot",
    "phase_tensor_map": "Geographic map of phase-tensor ellipses",
    "station_response": "Per-station impedance response (Bode) curves",
    "strike_profile": "Geoelectric strike vs station position",
    "strike": "Geoelectric strike analyzer (per-station table)",
    "dimensionality": "1-D / 2-D / 3-D dimensionality classifier",
    "validator": "Per-station EDI quality checklist",
    # --- data / IO tools ---
    "coords": "Transform station lat/lon to UTM coordinates",
    "elevation": "Enrich stations with elevation from a web API",
    "converter": "Re-export the survey to CSV / JSON / EDI",
    "batch_export": "Render and save a bundle of standard plots",
    # --- stateful tools ---
    "freq_editor": "Confidence-based frequency QC (drop/mask/recover)",
    "layered_model": "Build & preview a synthetic 1-D resistivity model",
}


# ── natural-language keyword mapping (ordered, most-specific first) ──────────────
# First-match wins.  Narrow phrases precede bare keywords so that, e.g.,
# "gcn" does not get swallowed by "inversion".

WORKFLOW_KEYWORDS: dict[str, list[str]] = {
    # ── plotting tasks (highest priority: the phase-tensor pseudo-section
    # must win over the bare "phase tensor" → phase_analysis match) ──────────
    # pt_strip_grid / pt_strip come before pt_psection: "phase tensor ellipse
    # strip" would otherwise be swallowed by pt_psection's "phase tensor
    # ellipse" substring (first-match-wins over dict insertion order).
    "pt_strip_grid": [
        "phase tensor strip grid",
        "phase-tensor strip grid",
        "pt strip grid",
        "ellipse strip grid",
        "phase tensor by line",
        "phase tensor per line",
        "phase tensor across lines",
        "multi-profile phase tensor",
        "multi profile phase tensor",
        "phase tensor grid",
        "strip grid",
    ],
    "pt_strip": [
        "phase tensor strip",
        "phase-tensor strip",
        "pt strip",
        "ellipse strip",
        "ellipse timeseries",
        "ellipse time series",
        "phase tensor ellipse strip",
        "phase tensor for one station",
        "phase tensor for station",
        "single station phase tensor",
        "single-station phase tensor",
        "ellipse vs period",
        "ellipse versus period",
    ],
    "pt_psection": [
        "phase tensor pseudosection",
        "phase tensor pseudo-section",
        "phase tensor psection",
        "phase tensor section",
        "phase tensor ellipse",
        "pt pseudosection",
        "pt pseudo-section",
        "pt psection",
        "pt section",
        "phi pseudosection",
        "phi tensor",
        "ellipse pseudosection",
        "ellipse section",
    ],
    "phase_tensor_map": [
        "phase tensor map",
        "phase-tensor map",
        "pt map",
        "phase tensor ellipse map",
        "map of phase tensor",
        "phase tensor on a map",
        "phase tensor geographic",
        "ellipse map",
    ],
    "phase_psection": [
        "phase pseudosection",
        "phase pseudo-section",
        "phase psection",
        "phase section",
        "pseudosection phase",
        "pseudosection of phase",
        "phase pseudo section",
        "pseudosection",
        "pseudo-section",
        "psection",
    ],
    "rhophi": [
        "rho phi",
        "rho/phi",
        "rhophi",
        "rho-phi",
        "rho and phi",
        "resistivity and phase",
        "apparent resistivity and phase",
        "res phase",
        "resphase",
        "res/phase",
        "sounding curve",
        "sounding curves",
        "rho phase",
        "rho/phase",
        "rho phi plot",
        "rho phi block",
        "rho phi curve",
    ],
    "strike_profile": [
        "strike profile",
        "strike vs",
        "strike along",
        "strike line chart",
        "strike viewer",
        "strike versus position",
        "strike per station",
    ],
    "station_response": [
        "station response",
        "response inspector",
        "impedance response",
        "response curve",
        "inspect station",
        "station inspector",
        "bode plot",
        "bode diagram",
    ],
    # ── analysis tools (specific triggers; bare "strike"/"dimensionality"
    # stay with the full phase_analysis workflow) ────────────────────────────
    "strike": [
        "strike analyzer",
        "strike analyser",
        "analyze strike",
        "analyse strike",
        "estimate strike",
        "regional strike",
        "strike rose",
        "strike estimation",
        "strike table",
    ],
    "dimensionality": [
        "dimensionality classifier",
        "classify dimensionality",
        "classify dimension",
        "dimensionality classification",
        "1d 2d 3d",
        "1d/2d/3d",
        "dimensionality table",
    ],
    "validator": [
        "edi validator",
        "validator",
        "validate edi",
        "validate the data",
        "validate data",
        "validate stations",
        "data quality checklist",
        "quality checklist",
        "qc checklist",
    ],
    # ── data / IO tools (specific multi-word triggers so bare verbs like
    # "convert" / "export" / "batch" don't steal them) ───────────────────────
    "coords": [
        "coordinate transform",
        "coordinate transformer",
        "coordinate conversion",
        "transform coordinates",
        "transform the coordinates",
        "convert coordinates",
        "lat lon to utm",
        "latlon to utm",
        "lat/lon to utm",
        "ll to utm",
        "utm conversion",
        "to utm",
        "reproject coordinates",
        "reproject stations",
    ],
    "elevation": [
        "elevation enrichment",
        "elevation enrich",
        "enrich elevation",
        "enrich with elevation",
        "fetch elevation",
        "get elevation",
        "add elevation",
        "station elevation",
        "elevation from api",
        "elevation lookup",
        "dem lookup",
        "lookup elevation",
    ],
    "batch_export": [
        "batch export",
        "batch plot export",
        "batch-export",
        "export plots",
        "export all plots",
        "export figures",
        "export all figures",
        "save all plots",
        "save all figures",
        "save the plots",
        "save figures to",
    ],
    "converter": [
        "format converter",
        "format conversion",
        "convert format",
        "convert the format",
        "convert to csv",
        "convert to json",
        "convert to edi",
        "export to csv",
        "export to json",
        "export survey",
        "export the survey",
        "re-export edi",
        "reexport edi",
        "export station metadata",
        "convert survey",
    ],
    # ── stateful tools ───────────────────────────────────────────────────────
    # freq_editor before freq_decimation so "edit frequencies by confidence"
    # is not swallowed by the decimation route.
    "freq_editor": [
        "frequency editor",
        "frequency edit",
        "edit frequencies",
        "edit the frequencies",
        "edit frequency",
        "confidence qc",
        "confidence-based qc",
        "confidence based qc",
        "recover frequencies",
        "drop frequencies",
        "mask frequencies",
        "recover periods",
        "frequency confidence",
    ],
    "layered_model": [
        "layered model",
        "layer model",
        "layered earth",
        "model builder",
        "build a model",
        "build model",
        "resistivity model",
        "1d model",
        "1-d model",
        "synthetic model",
        "earth model",
        "preview model",
    ],
    "tipper_plot": [
        "tipper plot",
        "plot tipper",
        "plot the tipper",
        "tipper component",
        "tipper components",
        "tipper curve",
        "tipper curves",
        "induction arrow",
        "induction arrows",
        "tipper arrow",
        "tipper arrows",
        "tipper map",
        "show tipper",
        "view tipper",
        "plot induction",
    ],
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
        "tem+mt",
        "mt+tem",
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
        "3d ai",
        "3d gcn",
        "spatial inversion",
        "3d neural",
        "graph network",
        "inv3d",
    ],
    "inv2d": [
        "unet",
        "u-net",
        "2d ai",
        "2d neural",
        "2d deep",
        "profile inversion",
        "lateral continuity ai",
        "inv2d",
    ],
    "pinn_inversion": [
        "pinn inversion",
        "pinn",
        "physics-informed",
        "physics informed",
        "no training data",
        "gradient descent inversion",
    ],
    "hybrid_inversion": [
        "hybrid inversion",
        "two-stage",
        "two stage",
        "ai + physics",
        "warm start",
        "ai warmstart",
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
        "rms",
        "misfit",
        "residual pt",
        "evaluate the inversion",
        "evaluate the result",
        "evaluate the results",
        "assess the inversion",
        "assess the result",
    ],
    # mare2dem before modem / pre_inversion: "prepare the inversion for
    # mare2dem" must not be swallowed by the generic "prepare the
    # inversion" → pre_inversion phrases.
    "mare2dem": [
        "mare2dem",
        "mare 2dem",
        "mare2d",
        "2.5d inversion",
        "2.5-d inversion",
        "2.5d",
        "2.5-d",
    ],
    "modem": [
        "modem",
        "3d inversion",
        "3d model",
    ],
    "pre_inversion": [
        "occam2d",
        "occam",
        "2d inversion",
        "mesh",
        "startup",
        "pre-inversion",
        "pre inversion",
        "inversion prep",
        "prepare inversion",
        "prepare the inversion",
        "prepare an inversion",
        "inversion preparation",
        "set up the inversion",
        "setup the inversion",
        "inversion data file",
        "inversion data",
        "data file for inversion",
        "build inversion input",
        "build inversion inputs",
        "inversion input",
        "inversion inputs",
        "third party 2d inversion",
        "third-party 2d inversion",
        "third party 2-d inversion",
        "third-party 2-d inversion",
        "inversion export",
        "before inversion export",
        "inversion-preparation diagnostics",
        "pre inversion coverage",
    ],
    "tipper": [
        "tipper",
        "induction arrow",
        "induction arrows",
        "vertical magnetic transfer",
        "vertical magnetic transfer function",
        "vertical magnetic transfer functions",
        "vertical transfer",
        "magnetic transfer function",
        "magnetic transfer functions",
        "wiese",
        "parkinson",
    ],
    "sensitivity": [
        "sensitivity",
        "depth of investigation",
        "depth-of-investigation",
        "investigation depth",
        "how far down",
        "trust my model",
        "depth sensitivity",
        "sensitivity versus depth",
        "bostick",
    ],
    "rotation": [
        "rotate",
        "rotation",
        "strike rotation",
        "tensor rotation",
    ],
    # freq_decimation before ai_inversion so "decimate/period range" never
    # routes to the bare "inversion" catch-all below.
    "freq_decimation": [
        "decimate",
        "decimation",
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
        "write code",
        "write python code",
        "generate code",
        "generate a script",
        "generate script",
        "python script",
        "code for",
        "script for",
        "script that",
        "script to",
        "write a script",
        "write a python script",
        "write function",
        "write class",
        "create notebook",
        "notebook for",
    ],
    "report": [
        "generate report",
        "write report",
        "survey report",
        "create report",
        "make report",
        # natural articled / adjectived phrasings
        "generate a report",
        "write a report",
        "create a report",
        "make a report",
        "produce a report",
        "reproducible report",
        "pdf report",
        "final report",
        "full report",
    ],
    "denoise": [
        "denoise",
        "denoising",
        "remove noise",
        "noise reduction",
        "filter noise",
        "remove spikes",
        "spikes",
        "outliers",
        "robust filtering",
        "hampel",
        "glitches",
        "clean cultural noise",
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
        "neural",
        "cnn",
        "deep learning",
        "machine learning",
        "inverter",
        "trained inverter",
        "inv1d",
        "1d ai",
        "1d neural",
        "1d inversion",
        "1-d checkpoint",
        "1d checkpoint",
        "checkpoint validated",
        "forward checked",
        "forward-checked",
        "resistivity model",
        "inversion",
        "invert",
        "inverting",
        "run inversion",
        "do inversion",
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
        "qc",
        "quality control",
        "quality check",
        "data quality",
        "remove noisy",
        "remove bad data",
        "bad frequencies",
        "noisy stations",
        "low-confidence",
        "low confidence",
        "flag low confidence",
        "station rejected",
        "rejected during preprocessing",
        "qc report",
    ],
}


# ── alias normalisation ─────────────────────────────────────────────────────────

WORKFLOW_ALIASES: dict[str, str] = {
    "qc_preprocessing": "qc",
    "preprocessing": "qc",
    # galvanic distortion -> static shift
    "galvanic": "static_shift",
    "pre-inversion": "pre_inversion",
    "pre_inversion_prep": "pre_inversion",
    "post-inversion": "inversion_eval",
    "interpretation": "interpretation",
    "report": "report",
    "full_pipeline": "full",
    "1d_inversion": "ai_inversion",
    "1d_ai": "ai_inversion",
    "ai_1d": "ai_inversion",
    "ai1d": "ai_inversion",
    "2d_inversion": "inv2d",
    "ai_2d": "inv2d",
    "3d_inversion": "inv3d",
    "ai_3d": "inv3d",
    "ensemble": "ensemble_inversion",
    "joint": "joint_inversion",
    "occam": "pre_inversion",
    "occam2d": "pre_inversion",
    "mare_2dem": "mare2dem",
    "mare2d": "mare2dem",
    "full_ai": "full_ai_workflow",
    # PINN / Hybrid aliases
    "pinn": "pinn_inversion",
    "physics_informed": "pinn_inversion",
    "physics-informed": "pinn_inversion",
    "pinn_inv": "pinn_inversion",
    "hybrid": "hybrid_inversion",
    "two_stage": "hybrid_inversion",
    "hybrid_inv": "hybrid_inversion",
    "ai_physics": "hybrid_inversion",
    # Agent-focused workflow aliases
    "code_generation": "code_gen",
    "codegen": "code_gen",
    "write_code": "code_gen",
    "denoising": "denoise",
    "noise_removal": "denoise",
    "doi": "sensitivity",
    "depth_of_investigation": "sensitivity",
    "frequency_decimation": "freq_decimation",
    "freq_dec": "freq_decimation",
    "period_decimation": "freq_decimation",
}


# ── correction methods (from the central registry) ──────────────────────────────
# The desktop *correction section* methods become first-class workflow ids.
# Their keyword entries are PREPENDED so a specific correction phrase
# ("ama static shift") wins over the generic "static shift" → static_shift
# match (classification is first-match-wins). Keep keywords specific.
def _register_correction_workflows() -> None:
    global WORKFLOW_KEYWORDS
    try:
        from ._corrections import CORRECTION_METHODS
    except Exception:  # noqa: BLE001 — corrections are optional
        return
    corr_kw = {
        wf: list(meta.get("keywords", []))
        for wf, meta in CORRECTION_METHODS.items()
        if meta.get("keywords")
    }
    for wf, meta in CORRECTION_METHODS.items():
        WORKFLOW_DESCRIPTIONS.setdefault(
            wf, meta.get("description", meta.get("title", wf))
        )
    # corrections first, then the existing (more generic) table
    WORKFLOW_KEYWORDS = {**corr_kw, **WORKFLOW_KEYWORDS}


_register_correction_workflows()


# ── public API ──────────────────────────────────────────────────────────────────

# Pipeline-stage families for compound-request detection. A request that
# names stages from three or more distinct families ("run QC, prepare the
# inversion, invert with AI and write a report") describes the *full*
# pipeline — returning the first keyword that happens to match would
# silently execute only one stage of it. Plot / tool / correction ids are
# deliberately absent: they never make a request "compound".
_STAGE_FAMILY: dict[str, str] = {
    "qc": "qc",
    "static_shift": "preprocess",
    "denoise": "preprocess",
    "rotation": "preprocess",
    "freq_decimation": "preprocess",
    "phase_analysis": "analysis",
    "tipper": "analysis",
    "sensitivity": "analysis",
    "dimensionality": "analysis",
    "pre_inversion": "prep",
    "modem": "prep",
    "mare2dem": "prep",
    "ai_inversion": "inversion",
    "inv2d": "inversion",
    "inv3d": "inversion",
    "pinn_inversion": "inversion",
    "hybrid_inversion": "inversion",
    "ensemble_inversion": "inversion",
    "joint_inversion": "inversion",
    "full_ai_workflow": "inversion",
    "inversion_eval": "evaluate",
    "comparison": "evaluate",
    "report": "report",
    "interpretation": "report",
}

_COMPOUND_MIN_FAMILIES = 3


def classify_workflow(
    text: str,
    *,
    default: str | None = None,
) -> str | None:
    """Classify free *text* into a workflow id via :data:`WORKFLOW_KEYWORDS`.

    First-match-wins over the ordered keyword table (case-insensitive),
    with one exception: a *compound* request that names pipeline stages
    from :data:`_COMPOUND_MIN_FAMILIES` or more distinct families (QC,
    preprocessing, inversion prep, inversion, evaluation, report, …)
    classifies as ``"full"`` — the end-to-end pipeline — instead of the
    first stage that happens to match.

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
    >>> classify_workflow(
    ...     "run quality control, prepare the inversion, run the AI "
    ...     "inversion and generate a report"
    ... )
    'full'
    >>> classify_workflow("xyzzy", default="qc")
    'qc'
    >>> classify_workflow("xyzzy") is None
    True
    """
    t = (text or "").lower()
    matched = [
        wf for wf, kws in WORKFLOW_KEYWORDS.items() if any(kw in t for kw in kws)
    ]
    if not matched:
        return default
    families = {_STAGE_FAMILY[wf] for wf in matched if wf in _STAGE_FAMILY}
    if "full" not in matched and len(families) >= _COMPOUND_MIN_FAMILIES:
        return "full"
    return matched[0]


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
    wf = str(workflow or "").lower().replace("-", "_").replace(" ", "_")
    return WORKFLOW_ALIASES.get(wf, wf) or "qc"
