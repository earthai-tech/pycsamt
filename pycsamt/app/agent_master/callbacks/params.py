# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Smart parameter-collection modal callbacks."""

from __future__ import annotations

import threading
from typing import Any

from dash import ALL, Input, Output, State, ctx
from dash import html, dcc, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from .._ids import IDs
from .chat import (
    _new_job,
    _run_agent,
    _thinking_bubble,
    _ts,
)

# ── Workflow-level field IDs (pattern-matching)
# type="am-pf" lets _submit_params collect
# all rendered workflow fields via ALL.

_FID_N_LAYERS  = {
    "type": "am-pf", "key": "n_layers",
}
_FID_DEPTH_MAX = {
    "type": "am-pf", "key": "depth_max",
}
_FID_EPOCHS    = {
    "type": "am-pf", "key": "epochs",
}
_FID_LR        = {
    "type": "am-pf", "key": "lr",
}
_FID_SMOOTH    = {
    "type": "am-pf",
    "key": "smoothness_weight",
}
_FID_LATERAL   = {
    "type": "am-pf",
    "key": "lateral_weight",
}
_FID_GRAPH_W   = {
    "type": "am-pf", "key": "graph_weight",
}
_FID_RADIUS    = {
    "type": "am-pf", "key": "radius",
}
_FID_INV_CODE  = {
    "type": "am-pf",
    "key": "inversion_code",
}

# ── Pipeline step field IDs (pattern-matching)
# type="am-ps" for step-level fields.

_FID_PERIOD_MIN = {
    "type": "am-ps", "key": "period_min",
}
_FID_PERIOD_MAX = {
    "type": "am-ps", "key": "period_max",
}
_FID_COMPONENT  = {
    "type": "am-ps", "key": "component",
}
_FID_SNR        = {
    "type": "am-ps",
    "key": "snr_threshold",
}
_FID_SS_METHOD  = {
    "type": "am-ps", "key": "ss_method",
}
_FID_DENOISE    = {
    "type": "am-ps",
    "key": "denoise_method",
}
_FID_PLOT_STYLE = {
    "type": "am-ps", "key": "plot_style",
}
_FID_RPT_FMT    = {
    "type": "am-ps",
    "key": "report_format",
}

# ── Reusable step section definitions ────────

_STEP_SECTIONS: dict[str, dict] = {
    "load": {
        "title": "Load & Filter",
        "icon": "bi-folder2-open",
        "fields": [
            {
                "id": _FID_PERIOD_MIN,
                "key": "period_min",
                "label": "Min period (s)",
                "type": "number",
                "min": 1e-6,
                "max": 1000.0,
                "step": None,
                "default": 0.0001,
                "help": (
                    "Shortest period to keep."
                ),
            },
            {
                "id": _FID_PERIOD_MAX,
                "key": "period_max",
                "label": "Max period (s)",
                "type": "number",
                "min": 0.001,
                "max": 100000.0,
                "step": None,
                "default": 1.0,
                "help": (
                    "Longest period to keep."
                ),
            },
            {
                "id": _FID_COMPONENT,
                "key": "component",
                "label": "Component",
                "type": "select",
                "options": [
                    {
                        "label": "xy  (E-polarisation)",
                        "value": "xy",
                    },
                    {
                        "label": "yx  (H-polarisation)",
                        "value": "yx",
                    },
                    {
                        "label": "all  (both)",
                        "value": "all",
                    },
                ],
                "default": "xy",
            },
        ],
    },
    "qc": {
        "title": "Quality Control",
        "icon": "bi-shield-check",
        "fields": [
            {
                "id": _FID_SNR,
                "key": "snr_threshold",
                "label": "SNR threshold",
                "type": "slider",
                "min": 0.1,
                "max": 10.0,
                "step": 0.1,
                "default": 0.5,
                "help": (
                    "Min signal-to-noise"
                    " ratio to keep."
                ),
            },
        ],
    },
    "static_shift": {
        "title": "Static Shift",
        "icon": "bi-arrows-expand-vertical",
        "fields": [
            {
                "id": _FID_SS_METHOD,
                "key": "ss_method",
                "label": "Correction method",
                "type": "radio",
                "options": [
                    {
                        "label": "SPI (default)",
                        "value": "spi",
                    },
                    {
                        "label": "Distortion",
                        "value": "distortion",
                    },
                    {
                        "label": "Skip",
                        "value": "none",
                    },
                ],
                "default": "spi",
                "inline": True,
            },
        ],
    },
    "denoise": {
        "title": "Denoising",
        "icon": "bi-soundwave",
        "fields": [
            {
                "id": _FID_DENOISE,
                "key": "denoise_method",
                "label": "Method",
                "type": "radio",
                "options": [
                    {
                        "label": "Wavelet",
                        "value": "wavelet",
                    },
                    {
                        "label": "Median",
                        "value": "median",
                    },
                    {
                        "label": "None",
                        "value": "none",
                    },
                ],
                "default": "wavelet",
                "inline": True,
            },
        ],
    },
    "phase_analysis": {
        "title": "Phase Tensor Analysis",
        "icon": "bi-graph-up-arrow",
        "fields": [
            {
                "id": _FID_PLOT_STYLE,
                "key": "plot_style",
                "label": "Plot type",
                "type": "select",
                "options": [
                    {
                        "label": "Mohr circle",
                        "value": "mohr",
                    },
                    {
                        "label": "Argand diagram",
                        "value": "argand",
                    },
                    {
                        "label": "Polar diagram",
                        "value": "polar",
                    },
                    {
                        "label": "Strike rose",
                        "value": "rose",
                    },
                ],
                "default": "mohr",
            },
        ],
    },
    "report": {
        "title": "Report",
        "icon": "bi-file-earmark-text",
        "fields": [
            {
                "id": _FID_RPT_FMT,
                "key": "report_format",
                "label": "Output format",
                "type": "radio",
                "options": [
                    {
                        "label": "HTML",
                        "value": "html",
                    },
                    {
                        "label": "Text",
                        "value": "text",
                    },
                ],
                "default": "html",
                "inline": True,
            },
        ],
    },
}

# ── Workflow schemas ──────────────────────────

_SCHEMAS: dict[str, dict] = {
    "ai_inversion": {
        "title": "1-D AI Inversion",
        "icon": "bi-lightning-charge-fill",
        "color": "var(--mauve)",
        "desc": (
            "Neural-network 1-D inversion. "
            "Defaults work for most "
            "MT / AMT datasets."
        ),
        "fields": [
            {
                "id": _FID_N_LAYERS,
                "key": "n_layers",
                "label": "Number of layers",
                "type": "slider",
                "min": 2, "max": 50,
                "step": 1, "default": 10,
                "help": (
                    "Depth layers in the "
                    "resistivity model (2-50)."
                ),
            },
            {
                "id": _FID_DEPTH_MAX,
                "key": "depth_max",
                "label": "Max depth (m)",
                "type": "number",
                "min": 100, "max": 50000,
                "step": 100, "default": 2000.0,
                "help": (
                    "Bottom of the model "
                    "column in metres."
                ),
            },
            {
                "id": _FID_EPOCHS,
                "key": "epochs",
                "label": "Training epochs",
                "type": "slider",
                "min": 100, "max": 3000,
                "step": 100, "default": 500,
                "help": (
                    "More epochs improve fit "
                    "at the cost of speed."
                ),
            },
            {
                "id": _FID_LR,
                "key": "lr",
                "label": "Learning rate",
                "type": "select",
                "options": [
                    {
                        "label": "0.001  slow/stable",
                        "value": 0.001,
                    },
                    {
                        "label": "0.01   default",
                        "value": 0.01,
                    },
                    {
                        "label": "0.05   fast",
                        "value": 0.05,
                    },
                    {
                        "label": "0.1    aggressive",
                        "value": 0.1,
                    },
                ],
                "default": 0.01,
            },
        ],
        "steps": ["load", "qc", "denoise"],
    },
    "inv1d": None,
    "pinn_inversion": {
        "title": "PINN Inversion",
        "icon": "bi-cpu-fill",
        "color": "var(--blue)",
        "desc": (
            "Physics-informed neural-network "
            "inversion. No labelled training "
            "data required."
        ),
        "fields": [
            {
                "id": _FID_N_LAYERS,
                "key": "n_layers",
                "label": "Number of layers",
                "type": "slider",
                "min": 2, "max": 50,
                "step": 1, "default": 10,
            },
            {
                "id": _FID_DEPTH_MAX,
                "key": "depth_max",
                "label": "Max depth (m)",
                "type": "number",
                "min": 100, "max": 50000,
                "step": 100, "default": 2000.0,
            },
            {
                "id": _FID_EPOCHS,
                "key": "epochs",
                "label": "Optimisation epochs",
                "type": "slider",
                "min": 200, "max": 5000,
                "step": 200, "default": 1000,
                "help": (
                    "PINN needs more iterations "
                    "than pure ML."
                ),
            },
            {
                "id": _FID_LR,
                "key": "lr",
                "label": "Learning rate",
                "type": "select",
                "options": [
                    {
                        "label": "0.001  slow/stable",
                        "value": 0.001,
                    },
                    {
                        "label": "0.01   default",
                        "value": 0.01,
                    },
                    {
                        "label": "0.05   fast",
                        "value": 0.05,
                    },
                ],
                "default": 0.01,
            },
            {
                "id": _FID_SMOOTH,
                "key": "smoothness_weight",
                "label": "Smoothness weight",
                "type": "slider",
                "min": 0.0, "max": 0.1,
                "step": 0.005, "default": 0.01,
                "help": (
                    "Penalty for rapid "
                    "resistivity changes."
                ),
            },
        ],
        "steps": ["load", "qc"],
    },
    "hybrid_inversion": {
        "title": "Hybrid Inversion",
        "icon": "bi-diagram-3-fill",
        "color": "var(--green)",
        "desc": (
            "Two-stage: AI warm-start then "
            "PINN physics refinement."
        ),
        "fields": [
            {
                "id": _FID_N_LAYERS,
                "key": "n_layers",
                "label": "Number of layers",
                "type": "slider",
                "min": 2, "max": 50,
                "step": 1, "default": 10,
            },
            {
                "id": _FID_DEPTH_MAX,
                "key": "depth_max",
                "label": "Max depth (m)",
                "type": "number",
                "min": 100, "max": 50000,
                "step": 100, "default": 2000.0,
            },
            {
                "id": _FID_EPOCHS,
                "key": "epochs",
                "label": "PINN refinement epochs",
                "type": "slider",
                "min": 100, "max": 3000,
                "step": 100, "default": 500,
            },
            {
                "id": _FID_SMOOTH,
                "key": "smoothness_weight",
                "label": "Smoothness weight",
                "type": "slider",
                "min": 0.0, "max": 0.1,
                "step": 0.005, "default": 0.01,
            },
        ],
        "steps": ["load", "qc"],
    },
    "inv2d": {
        "title": "2-D AI Inversion (U-Net)",
        "icon": "bi-layers-fill",
        "color": "var(--yellow)",
        "desc": (
            "U-Net profile inversion. "
            "Requires >= 3 collinear sites."
        ),
        "fields": [
            {
                "id": _FID_N_LAYERS,
                "key": "n_layers",
                "label": "Depth layers",
                "type": "slider",
                "min": 5, "max": 100,
                "step": 5, "default": 30,
            },
            {
                "id": _FID_DEPTH_MAX,
                "key": "depth_max",
                "label": "Max depth (m)",
                "type": "number",
                "min": 500, "max": 50000,
                "step": 500, "default": 5000.0,
            },
            {
                "id": _FID_LATERAL,
                "key": "lateral_weight",
                "label": "Lateral smoothness",
                "type": "slider",
                "min": 0.0, "max": 0.1,
                "step": 0.005, "default": 0.005,
                "help": (
                    "Penalty for lateral "
                    "resistivity jumps."
                ),
            },
            {
                "id": _FID_EPOCHS,
                "key": "epochs",
                "label": "Training epochs",
                "type": "slider",
                "min": 100, "max": 3000,
                "step": 100, "default": 500,
            },
        ],
        "steps": ["load", "qc"],
    },
    "inv3d": {
        "title": "3-D GCN Inversion",
        "icon": "bi-boxes",
        "color": "var(--peach)",
        "desc": (
            "Graph-convolutional 3-D "
            "spatial inversion. "
            "All loaded profiles are used."
        ),
        "fields": [
            {
                "id": _FID_N_LAYERS,
                "key": "n_layers",
                "label": "Depth layers",
                "type": "slider",
                "min": 5, "max": 60,
                "step": 5, "default": 20,
            },
            {
                "id": _FID_DEPTH_MAX,
                "key": "depth_max",
                "label": "Max depth (m)",
                "type": "number",
                "min": 500, "max": 50000,
                "step": 500, "default": 5000.0,
            },
            {
                "id": _FID_RADIUS,
                "key": "radius",
                "label": "Graph radius (m)",
                "type": "number",
                "min": 500, "max": 50000,
                "step": 500, "default": 5000.0,
                "help": (
                    "Sites within this distance "
                    "are linked as graph nodes."
                ),
            },
            {
                "id": _FID_GRAPH_W,
                "key": "graph_weight",
                "label": "Graph smoothness",
                "type": "slider",
                "min": 0.0, "max": 0.05,
                "step": 0.001,
                "default": 0.005,
            },
        ],
        "steps": ["load", "qc"],
    },
    "ensemble_inversion": {
        "title": "Ensemble Inversion",
        "icon": "bi-collection-fill",
        "color": "var(--teal)",
        "desc": (
            "Multi-model ensemble with "
            "uncertainty quantification."
        ),
        "fields": [
            {
                "id": _FID_N_LAYERS,
                "key": "n_layers",
                "label": "Number of layers",
                "type": "slider",
                "min": 2, "max": 50,
                "step": 1, "default": 10,
            },
            {
                "id": _FID_DEPTH_MAX,
                "key": "depth_max",
                "label": "Max depth (m)",
                "type": "number",
                "min": 100, "max": 50000,
                "step": 100, "default": 2000.0,
            },
            {
                "id": _FID_EPOCHS,
                "key": "epochs",
                "label": "Training epochs",
                "type": "slider",
                "min": 100, "max": 3000,
                "step": 100, "default": 500,
            },
        ],
        "steps": ["load", "qc", "denoise"],
    },
    "pre_inversion": {
        "title": "Pre-Inversion Setup",
        "icon": "bi-file-earmark-code-fill",
        "color": "var(--sapphire)",
        "desc": (
            "Prepare input files for "
            "classical 2-D or 3-D codes."
        ),
        "fields": [
            {
                "id": _FID_INV_CODE,
                "key": "inversion_code",
                "label": "Inversion code",
                "type": "radio",
                "options": [
                    {"label": "Occam2D",
                     "value": "occam2d"},
                    {"label": "ModEM",
                     "value": "modem"},
                ],
                "default": "occam2d",
                "inline": True,
            },
            {
                "id": _FID_N_LAYERS,
                "key": "n_layers",
                "label": "Number of layers",
                "type": "slider",
                "min": 5, "max": 100,
                "step": 5, "default": 30,
            },
            {
                "id": _FID_DEPTH_MAX,
                "key": "depth_max",
                "label": "Max depth (m)",
                "type": "number",
                "min": 500, "max": 100000,
                "step": 500, "default": 10000.0,
            },
        ],
        "steps": ["load", "qc", "static_shift"],
    },
    "modem": {
        "title": "ModEM 3-D Setup",
        "icon": "bi-grid-3x3-gap-fill",
        "color": "var(--red)",
        "desc": (
            "Generate ModEM input data "
            "and mesh configuration."
        ),
        "fields": [
            {
                "id": _FID_N_LAYERS,
                "key": "n_layers",
                "label": "Vertical layers",
                "type": "slider",
                "min": 10, "max": 100,
                "step": 5, "default": 30,
            },
            {
                "id": _FID_DEPTH_MAX,
                "key": "depth_max",
                "label": "Max depth (m)",
                "type": "number",
                "min": 1000, "max": 200000,
                "step": 1000, "default": 50000.0,
            },
        ],
        "steps": ["load", "qc"],
    },
    # ── Pipeline-only workflows ───────────────
    "qc": {
        "title": "QC Pipeline",
        "icon": "bi-shield-check",
        "color": "var(--green)",
        "desc": (
            "Clean and flag MT data. Configure "
            "period range, noise thresholds, "
            "and static shift correction."
        ),
        "fields": [],
        "steps": [
            "load", "qc", "static_shift", "report"
        ],
    },
    "phase_analysis": {
        "title": "Phase Tensor Analysis",
        "icon": "bi-graph-up-arrow",
        "color": "var(--blue)",
        "desc": (
            "Compute phase tensors and impedance "
            "strike. Configure data range and "
            "plot style."
        ),
        "fields": [],
        "steps": [
            "load", "qc", "phase_analysis",
            "report",
        ],
    },
    "static_shift": {
        "title": "Static Shift Correction",
        "icon": "bi-arrows-expand-vertical",
        "color": "var(--yellow)",
        "desc": (
            "Detect and correct galvanic "
            "distortion from MT impedances."
        ),
        "fields": [],
        "steps": [
            "load", "qc", "static_shift", "report"
        ],
    },
    "tipper": {
        "title": "Tipper Analysis",
        "icon": "bi-arrow-up-right-circle",
        "color": "var(--teal)",
        "desc": (
            "Analyse the vertical magnetic "
            "transfer function (tipper). Choose "
            "the tipper component to plot."
        ),
        "fields": [
            {
                "id": {
                    "type": "am-pf",
                    "key": "tipper_component",
                },
                "key": "tipper_component",
                "label": "Tipper component",
                "type": "radio",
                "options": [
                    {
                        "label": "Tx  (x-direction)",
                        "value": "Tx",
                    },
                    {
                        "label": "Ty  (y-direction)",
                        "value": "Ty",
                    },
                    {
                        "label": "Both",
                        "value": "all",
                    },
                ],
                "default": "Tx",
                "inline": True,
            },
        ],
        "steps": ["load", "qc", "report"],
    },
    "rotation": {
        "title": "Data Rotation",
        "icon": "bi-arrow-repeat",
        "color": "var(--peach)",
        "desc": (
            "Rotate MT impedance tensors to a "
            "specified strike direction. Enter "
            "the angle measured from geographic "
            "North (positive = clockwise)."
        ),
        "fields": [
            {
                "id": {
                    "type": "am-pf",
                    "key": "strike_angle",
                },
                "key": "strike_angle",
                "label": "Strike angle (deg)",
                "type": "number",
                "min": -180, "max": 180,
                "step": 5, "default": 0.0,
                "help": (
                    "Rotation from N, positive"
                    " clockwise."
                ),
            },
        ],
        "steps": ["load", "qc", "report"],
    },
    # ── Agent-focused workflows ───────────────
    "interpret": {
        "title": "Geological Interpretation",
        "icon": "bi-globe-americas",
        "color": "var(--teal)",
        "desc": (
            "Interpret resistivity models. "
            "Set the depth focus and geological "
            "context to guide the analysis."
        ),
        "fields": [
            {
                "id": _FID_DEPTH_MAX,
                "key": "depth_max",
                "label": "Max depth (m)",
                "type": "number",
                "min": 100, "max": 20000,
                "step": 100, "default": 2000.0,
                "help": (
                    "Limit interpretation to "
                    "this depth."
                ),
            },
            {
                "id": {
                    "type": "am-pf",
                    "key": "geology",
                },
                "key": "geology",
                "label": "Geological context",
                "type": "select",
                "options": [
                    {
                        "label": "Generic",
                        "value": "generic",
                    },
                    {
                        "label": "Sedimentary basin",
                        "value": "sedimentary",
                    },
                    {
                        "label": "Basement / crystalline",
                        "value": "basement",
                    },
                    {
                        "label": "Volcanic / hydrothermal",
                        "value": "volcanic",
                    },
                    {
                        "label": "Ore deposit",
                        "value": "ore",
                    },
                ],
                "default": "generic",
                "help": (
                    "Constrains the vocabulary "
                    "and focus of the LLM."
                ),
            },
            {
                "id": {
                    "type": "am-pf",
                    "key": "style",
                },
                "key": "style",
                "label": "Output style",
                "type": "radio",
                "options": [
                    {
                        "label": "Brief",
                        "value": "brief",
                    },
                    {
                        "label": "Detailed",
                        "value": "detailed",
                    },
                ],
                "default": "detailed",
                "inline": True,
            },
        ],
        "steps": [],
    },
    "report": {
        "title": "Survey Report",
        "icon": "bi-file-earmark-text-fill",
        "color": "var(--sapphire)",
        "desc": (
            "Generate a structured survey report. "
            "Choose format and detail level."
        ),
        "fields": [
            {
                "id": {
                    "type": "am-pf",
                    "key": "report_format",
                },
                "key": "report_format",
                "label": "Output format",
                "type": "radio",
                "options": [
                    {
                        "label": "HTML (browser)",
                        "value": "html",
                    },
                    {
                        "label": "Plain text",
                        "value": "text",
                    },
                ],
                "default": "html",
                "inline": True,
            },
            {
                "id": {
                    "type": "am-pf",
                    "key": "detail_level",
                },
                "key": "detail_level",
                "label": "Detail level",
                "type": "radio",
                "options": [
                    {
                        "label": "Summary",
                        "value": "summary",
                    },
                    {
                        "label": "Full",
                        "value": "full",
                    },
                    {
                        "label": "Technical",
                        "value": "technical",
                    },
                ],
                "default": "full",
                "inline": True,
            },
            {
                "id": {
                    "type": "am-pf",
                    "key": "include_figures",
                },
                "key": "include_figures",
                "label": "Include figures",
                "type": "radio",
                "options": [
                    {
                        "label": "Yes",
                        "value": "yes",
                    },
                    {
                        "label": "No",
                        "value": "no",
                    },
                ],
                "default": "yes",
                "inline": True,
            },
        ],
        "steps": [],
    },
    "code_gen": {
        "title": "Code Generation",
        "icon": "bi-code-slash",
        "color": "var(--green)",
        "desc": (
            "Generate Python code for your "
            "MT data task. Choose the output "
            "type and preferred framework."
        ),
        "fields": [
            {
                "id": {
                    "type": "am-pf",
                    "key": "output_type",
                },
                "key": "output_type",
                "label": "Output type",
                "type": "radio",
                "options": [
                    {
                        "label": "Script (.py)",
                        "value": "script",
                    },
                    {
                        "label": "Function",
                        "value": "function",
                    },
                    {
                        "label": "Class",
                        "value": "class",
                    },
                    {
                        "label": "Notebook (.ipynb)",
                        "value": "notebook",
                    },
                ],
                "default": "script",
                "inline": False,
            },
            {
                "id": {
                    "type": "am-pf",
                    "key": "framework",
                },
                "key": "framework",
                "label": "Framework",
                "type": "select",
                "options": [
                    {
                        "label": "NumPy / SciPy",
                        "value": "numpy",
                    },
                    {
                        "label": "PyTorch",
                        "value": "torch",
                    },
                    {
                        "label": "TensorFlow",
                        "value": "tensorflow",
                    },
                    {
                        "label": "pycsamt only",
                        "value": "pycsamt",
                    },
                ],
                "default": "numpy",
            },
            {
                "id": {
                    "type": "am-pf",
                    "key": "verbosity",
                },
                "key": "verbosity",
                "label": "Code verbosity",
                "type": "radio",
                "options": [
                    {
                        "label": "Minimal",
                        "value": "minimal",
                    },
                    {
                        "label": "Normal",
                        "value": "normal",
                    },
                    {
                        "label": "Verbose (with docs)",
                        "value": "verbose",
                    },
                ],
                "default": "normal",
                "inline": False,
            },
        ],
        "steps": [],
    },
    "denoise": {
        "title": "Data Denoising",
        "icon": "bi-soundwave",
        "color": "var(--blue)",
        "desc": (
            "Remove noise from MT impedances. "
            "Choose the method and the minimum "
            "signal-to-noise threshold."
        ),
        "fields": [
            {
                "id": {
                    "type": "am-pf",
                    "key": "denoise_method",
                },
                "key": "denoise_method",
                "label": "Denoising method",
                "type": "radio",
                "options": [
                    {
                        "label": "Wavelet",
                        "value": "wavelet",
                    },
                    {
                        "label": "Median filter",
                        "value": "median",
                    },
                    {
                        "label": "Skip",
                        "value": "none",
                    },
                ],
                "default": "wavelet",
                "inline": True,
            },
            {
                "id": {
                    "type": "am-pf",
                    "key": "snr_threshold",
                },
                "key": "snr_threshold",
                "label": "SNR threshold",
                "type": "slider",
                "min": 0.1, "max": 10.0,
                "step": 0.1, "default": 0.5,
                "help": (
                    "Periods below this SNR "
                    "are flagged and removed."
                ),
            },
        ],
        "steps": ["load"],
    },
    "sensitivity": {
        "title": "Sensitivity / DOI Analysis",
        "icon": "bi-layers-half",
        "color": "var(--mauve)",
        "desc": (
            "Estimate depth of investigation and "
            "model sensitivity. Choose method and "
            "maximum depth."
        ),
        "fields": [
            {
                "id": {
                    "type": "am-pf",
                    "key": "doi_method",
                },
                "key": "doi_method",
                "label": "DOI method",
                "type": "radio",
                "options": [
                    {
                        "label": "Bostick (fast)",
                        "value": "bostick",
                    },
                    {
                        "label": "Sensitivity kernel",
                        "value": "sensitivity_kernel",
                    },
                    {
                        "label": "Jacobian",
                        "value": "jacobian",
                    },
                ],
                "default": "bostick",
                "inline": False,
            },
            {
                "id": _FID_DEPTH_MAX,
                "key": "depth_max",
                "label": "Max depth (m)",
                "type": "number",
                "min": 100, "max": 20000,
                "step": 100, "default": 2000.0,
                "help": (
                    "Compute DOI down to "
                    "this depth."
                ),
            },
        ],
        "steps": [],
    },
}
# aliases
_SCHEMAS["inv1d"] = _SCHEMAS["ai_inversion"]
_SCHEMAS["interpretation"] = _SCHEMAS["interpret"]


# ── Field & section renderers ─────────────────

def _field_el(f: dict, val: Any) -> Any:
    fid = f["id"]
    ft = f["type"]
    if ft == "slider":
        return dcc.Slider(
            id=fid,
            min=f["min"],
            max=f["max"],
            step=f["step"],
            value=val,
            marks=None,
            tooltip={
                "placement": "bottom",
                "always_visible": True,
            },
            className="am-pf-slider",
        )
    if ft == "number":
        return dbc.Input(
            id=fid,
            type="number",
            value=val,
            min=f.get("min"),
            max=f.get("max"),
            step=f.get("step"),
            className="am-pf-number",
        )
    if ft == "select":
        return dcc.Dropdown(
            id=fid,
            options=f["options"],
            value=val,
            clearable=False,
            className="am-pf-select",
        )
    if ft == "radio":
        return dbc.RadioItems(
            id=fid,
            options=f["options"],
            value=val,
            inline=f.get("inline", False),
            className="am-pf-radio",
        )
    return html.Div(id=fid)


def _render_form(
    fields: list[dict],
    inv_config: dict,
) -> list:
    rows = []
    for f in fields:
        val = inv_config.get(
            f["key"], f["default"]
        )
        help_txt = f.get("help", "")
        rows.append(
            html.Div(
                [
                    html.Div(
                        [
                            html.Label(
                                f["label"],
                                className=(
                                    "am-pf-label"
                                ),
                            ),
                            html.Small(
                                help_txt,
                                className=(
                                    "am-pf-help"
                                ),
                            ) if help_txt else None,
                        ],
                        className="am-pf-label-row",
                    ),
                    _field_el(f, val),
                ],
                className="am-pf-row",
            )
        )
    return [r for r in rows if r]


def _render_steps_accordion(
    step_names: list[str],
    config: dict,
) -> html.Div | None:
    items = []
    active = []
    for sname in step_names:
        sdef = _STEP_SECTIONS.get(sname)
        if not sdef:
            continue
        sc = config.get("step_params", {}).get(
            sname, {}
        )
        merged = {**config, **sc}
        fields_html = _render_form(
            sdef["fields"], merged
        )
        item_id = f"am-pstep-{sname}"
        active.append(item_id)
        items.append(
            dbc.AccordionItem(
                fields_html,
                title=html.Span([
                    html.I(
                        className=(
                            f"{sdef['icon']} me-2"
                        ),
                    ),
                    sdef["title"],
                ]),
                item_id=item_id,
                class_name="am-step-item",
            )
        )
    if not items:
        return None
    return html.Div(
        [
            html.H6(
                [
                    html.I(
                        className=(
                            "bi bi-diagram-3 me-2"
                        ),
                    ),
                    "Pipeline Steps",
                ],
                className="am-pf-section-hdr",
            ),
            dbc.Accordion(
                items,
                always_open=True,
                active_item=active,
                className="am-step-accordion",
            ),
        ],
    )


def _collect_params(
    schema: dict,
    field_vals: dict,
    inv_config: dict,
) -> dict:
    result = dict(inv_config or {})
    for f in schema.get("fields", []):
        key = f["key"]
        val = field_vals.get(key)
        if val is None:
            val = f["default"]
        result[key] = val
    return result


def _collect_step_params(
    schema: dict,
    field_vals: dict,
) -> dict:
    result: dict[str, dict] = {}
    for sname in schema.get("steps", []):
        sdef = _STEP_SECTIONS.get(sname)
        if not sdef:
            continue
        sdata: dict = {}
        for f in sdef["fields"]:
            v = field_vals.get(f["key"])
            sdata[f["key"]] = (
                v
                if v is not None
                else f["default"]
            )
        result[sname] = sdata
    return result


def _cancel_bubble() -> html.Div:
    return html.Div(
        [
            html.Div(
                html.I(className="bi bi-robot"),
                className="am-avatar agent",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.I(
                                className=(
                                    "bi bi-x-circle"
                                    " me-2"
                                ),
                                style={
                                    "color": (
                                        "var(--red)"
                                    )
                                },
                            ),
                            "Workflow cancelled.",
                        ]
                    ),
                    html.Div(
                        _ts(), className="am-ts"
                    ),
                ],
                className="am-bubble agent",
            ),
        ],
        className="am-msg-row",
    )


# ── Callbacks ─────────────────────────────────

def register_params(app) -> None:

    # Open modal when STORE_PENDING is set
    @app.callback(
        Output(IDs.MODAL_PARAMS, "is_open"),
        Output(
            IDs.PARAM_MODAL_TITLE, "children"
        ),
        Output(
            IDs.PARAM_MODAL_DESC, "children"
        ),
        Output(
            IDs.PARAM_FORM_BODY, "children"
        ),
        Input(IDs.STORE_PENDING, "data"),
        State(IDs.STORE_INV_CONFIG, "data"),
        prevent_initial_call=True,
    )
    def _open_param_modal(pending, inv_config):
        if not pending:
            raise PreventUpdate
        wf = pending.get("workflow")
        schema = _SCHEMAS.get(wf)
        if not schema:
            raise PreventUpdate
        ic = inv_config or {}
        title = [
            html.I(
                className=(
                    f"{schema['icon']} me-2"
                ),
                style={"color": schema["color"]},
            ),
            schema["title"],
        ]
        form: list = []
        wf_fields = schema.get("fields", [])
        if wf_fields:
            form.extend(
                _render_form(wf_fields, ic)
            )
        step_names = schema.get("steps", [])
        if step_names:
            if wf_fields:
                form.append(
                    html.Hr(
                        className="am-pf-divider"
                    )
                )
            acc = _render_steps_accordion(
                step_names, ic
            )
            if acc:
                form.append(acc)
        return (
            True, title, schema["desc"], form
        )

    # Submit: collect form → start job
    @app.callback(
        Output(
            IDs.CHAT_WINDOW,
            "children",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_JOB,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.INTERVAL_POLL,
            "disabled",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_MESSAGES,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_INV_CONFIG,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_PENDING,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.MODAL_PARAMS,
            "is_open",
            allow_duplicate=True,
        ),
        Input(IDs.BTN_PARAM_RUN, "n_clicks"),
        State(IDs.STORE_PENDING, "data"),
        State(IDs.CHAT_WINDOW, "children"),
        State(IDs.STORE_EDI, "data"),
        State(IDs.STORE_SETTINGS, "data"),
        State(IDs.STORE_MESSAGES, "data"),
        State(IDs.STORE_INV_CONFIG, "data"),
        State(
            {"type": "am-pf", "key": ALL},
            "value",
        ),
        State(
            {"type": "am-ps", "key": ALL},
            "value",
        ),
        prevent_initial_call=True,
    )
    def _submit_params(
        n_run,
        pending,
        current_msgs,
        edi_store,
        settings,
        stored_msgs,
        inv_config,
        _pf,
        _ps,
    ):
        if not n_run or not pending:
            raise PreventUpdate
        wf = pending.get("workflow", "")
        schema = _SCHEMAS.get(wf, {})

        pf_vals = {
            s["id"]["key"]: s["value"]
            for s in (
                ctx.states_list[6] or []
            )
        }
        ps_vals = {
            s["id"]["key"]: s["value"]
            for s in (
                ctx.states_list[7] or []
            )
        }
        new_ic = _collect_params(
            schema,
            pf_vals,
            inv_config or {},
        )
        new_ic["step_params"] = (
            _collect_step_params(
                schema, ps_vals
            )
        )
        new_ic["workflow"] = wf

        # Replace waiting bubble with thinking
        msgs = list(current_msgs or [])
        for i, c in enumerate(msgs):
            cid = (
                c.get("props", {})
                 .get("id", "")
                if isinstance(c, dict)
                else ""
            )
            if cid == "am-waiting-bubble":
                msgs[i] = _thinking_bubble([{
                    "label": "Starting...",
                    "status": "running",
                }])
                break

        jid = _new_job()
        threading.Thread(
            target=_run_agent,
            args=(
                jid,
                pending.get("text", ""),
                edi_store or {},
                settings or {},
                new_ic,
            ),
            daemon=True,
        ).start()

        return (
            msgs,
            {"jid": jid},
            False,
            stored_msgs or [],
            new_ic,
            {},
            False,
        )

    # Cancel: close modal, add cancellation msg
    @app.callback(
        Output(
            IDs.MODAL_PARAMS,
            "is_open",
            allow_duplicate=True,
        ),
        Output(
            IDs.STORE_PENDING,
            "data",
            allow_duplicate=True,
        ),
        Output(
            IDs.CHAT_WINDOW,
            "children",
            allow_duplicate=True,
        ),
        Input(
            IDs.BTN_PARAM_CANCEL, "n_clicks"
        ),
        State(IDs.CHAT_WINDOW, "children"),
        prevent_initial_call=True,
    )
    def _cancel_params(n, current_msgs):
        if not n:
            raise PreventUpdate
        msgs = list(current_msgs or [])
        replaced = False
        for i, c in enumerate(msgs):
            cid = (
                c.get("props", {})
                 .get("id", "")
                if isinstance(c, dict)
                else ""
            )
            if cid == "am-waiting-bubble":
                msgs[i] = _cancel_bubble()
                replaced = True
                break
        if not replaced:
            msgs.append(_cancel_bubble())
        return False, {}, msgs
