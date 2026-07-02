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
_FID_CHECKPOINT = {
    "type": "am-pf", "key": "checkpoint",
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
    "type": "am-ps", "key": "method",
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

# ── Plotting-task fields (workflow-level → config → PlotAgent) ──
_FID_PLOT_STATIONS  = {"type": "am-pf", "key": "stations"}
_FID_PLOT_COMPONENT = {"type": "am-pf", "key": "components"}
_FID_PLOT_PMIN      = {"type": "am-pf", "key": "period_min"}
_FID_PLOT_PMAX      = {"type": "am-pf", "key": "period_max"}
_FID_PLOT_PUB       = {"type": "am-pf", "key": "publication"}
_FID_PLOT_ERRBAR    = {"type": "am-pf", "key": "errorbar"}
_FID_PLOT_COLORBY   = {"type": "am-pf", "key": "color_by"}
_FID_PLOT_SCALE     = {"type": "am-pf", "key": "scale"}
_FID_PLOT_VIEW      = {"type": "am-pf", "key": "view"}
_FID_PLOT_PARTS     = {"type": "am-pf", "key": "parts"}
_FID_PLOT_CONV      = {"type": "am-pf", "key": "convention"}
_FID_PLOT_PERIOD    = {"type": "am-pf", "key": "period"}
_FID_PLOT_METHOD    = {"type": "am-pf", "key": "method"}
_FID_PLOT_SORTBY    = {"type": "am-pf", "key": "sort_by"}
_FID_SKEW_TH        = {"type": "am-pf", "key": "skew_th"}
_FID_ELLIPT_TH      = {"type": "am-pf", "key": "ellipt_th"}
# ── Data / IO tool fields (Wave C) ──────────────────────────────
_FID_DATUM          = {"type": "am-pf", "key": "datum"}
_FID_ZONE           = {"type": "am-pf", "key": "zone"}
_FID_HEMISPHERE     = {"type": "am-pf", "key": "hemisphere"}
_FID_API            = {"type": "am-pf", "key": "api"}
_FID_FORMAT         = {"type": "am-pf", "key": "format"}
_FID_OUTDIR         = {"type": "am-pf", "key": "output_dir"}
_FID_PLOTS          = {"type": "am-pf", "key": "plots"}
_FID_DPI            = {"type": "am-pf", "key": "dpi"}
# ── Stateful tool fields (Wave D) ───────────────────────────────
_FID_MODE           = {"type": "am-pf", "key": "mode"}
_FID_THRESHOLD      = {"type": "am-pf", "key": "threshold"}
_FID_CI_HI          = {"type": "am-pf", "key": "ci_hi"}
_FID_CI_LO          = {"type": "am-pf", "key": "ci_lo"}
_FID_ALSO           = {"type": "am-pf", "key": "also"}
_FID_REJECT         = {"type": "am-pf", "key": "reject"}
_FID_PRESET         = {"type": "am-pf", "key": "preset"}
_FID_RESISTIVITIES  = {"type": "am-pf", "key": "resistivities"}
_FID_THICKNESSES    = {"type": "am-pf", "key": "thicknesses"}

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
                "key": "method",
                "label": "Correction method",
                "type": "radio",
                "options": [
                    {
                        "label": "AMA (default)",
                        "value": "ama",
                    },
                    {
                        "label": "Loess",
                        "value": "loess",
                    },
                    {
                        "label": "Ref-median",
                        "value": "refmedian",
                    },
                    {
                        "label": "Skip",
                        "value": "none",
                    },
                ],
                "default": "ama",
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
            {
                "id": _FID_CHECKPOINT,
                "key": "checkpoint",
                "label": "Model / checkpoint path",
                "type": "text",
                "default": "",
                "placeholder": (
                    "Leave blank to train"
                    " from scratch"
                ),
                "help": (
                    "Optional: path to a "
                    "pre-trained model file "
                    "or checkpoint folder to "
                    "skip full training."
                ),
            },
        ],
        "steps": [],
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
            {
                "id": _FID_CHECKPOINT,
                "key": "checkpoint",
                "label": (
                    "Checkpoint path (optional)"
                ),
                "type": "text",
                "default": "",
                "placeholder": (
                    "Leave blank to train"
                    " from scratch"
                ),
                "help": (
                    "Optional: resume PINN "
                    "training from a saved "
                    "checkpoint file or folder."
                ),
            },
        ],
        "steps": [],
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
            {
                "id": _FID_CHECKPOINT,
                "key": "checkpoint",
                "label": (
                    "AI model / checkpoint path"
                ),
                "type": "text",
                "default": "",
                "placeholder": (
                    "/path/to/model.pt or"
                    " checkpoint/"
                ),
                "help": (
                    "Required: path to a trained"
                    " AI inverter or checkpoint"
                    " folder. Run ai_inversion"
                    " first if none is available."
                ),
            },
        ],
        "steps": [],
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
            {
                "id": _FID_CHECKPOINT,
                "key": "checkpoint",
                "label": (
                    "U-Net checkpoint (optional)"
                ),
                "type": "text",
                "default": "",
                "placeholder": (
                    "Leave blank to train"
                    " from scratch"
                ),
                "help": (
                    "Optional: path to a "
                    "pre-trained U-Net model "
                    "file or checkpoint folder."
                ),
            },
        ],
        "steps": [],
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
            {
                "id": _FID_CHECKPOINT,
                "key": "checkpoint",
                "label": (
                    "GCN checkpoint (optional)"
                ),
                "type": "text",
                "default": "",
                "placeholder": (
                    "Leave blank to train"
                    " from scratch"
                ),
                "help": (
                    "Optional: path to a "
                    "pre-trained GCN model "
                    "file or checkpoint folder."
                ),
            },
        ],
        "steps": [],
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
        "steps": [],
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
                    {"label": "MARE2DEM",
                     "value": "mare2dem"},
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
        "steps": [],
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
            {
                "id": {"type": "am-pf",
                       "key": "error_floor"},
                "key": "error_floor",
                "label": "Error floor (rel.)",
                "type": "number",
                "min": 0.01, "max": 0.5,
                "step": 0.01, "default": 0.05,
            },
        ],
        "steps": [],
    },
    "mare2dem": {
        "title": "MARE2DEM 2.5-D Setup",
        "icon": "bi-water",
        "color": "var(--sapphire)",
        "desc": (
            "Generate MARE2DEM input files: "
            "emdata (TE/TM responses), starting "
            "resistivity, and settings."
        ),
        "fields": [
            {
                "id": {"type": "am-pf",
                       "key": "output_modes"},
                "key": "output_modes",
                "label": "Data modes",
                "type": "radio",
                "options": [
                    {"label": "TE + TM",
                     "value": "all impedance"},
                    {"label": "TE only",
                     "value": "TE"},
                    {"label": "TM only",
                     "value": "TM"},
                    {"label": "TE+TM+tipper",
                     "value": "all"},
                ],
                "default": "all impedance",
                "inline": True,
            },
            {
                "id": {"type": "am-pf",
                       "key": "error_floor"},
                "key": "error_floor",
                "label": "Error floor (rel.)",
                "type": "number",
                "min": 0.01, "max": 0.5,
                "step": 0.01, "default": 0.05,
            },
            {
                "id": {"type": "am-pf",
                       "key": "initial_rho"},
                "key": "initial_rho",
                "label": "Start ρ (Ω·m)",
                "type": "number",
                "min": 0.1, "max": 10000,
                "step": 0.1, "default": 100.0,
            },
            {
                "id": {"type": "am-pf",
                       "key": "target_rms"},
                "key": "target_rms",
                "label": "Target RMS",
                "type": "number",
                "min": 0.5, "max": 5.0,
                "step": 0.1, "default": 1.0,
            },
            {
                "id": {"type": "am-pf",
                       "key": "max_iterations"},
                "key": "max_iterations",
                "label": "Max iterations",
                "type": "slider",
                "min": 10, "max": 300,
                "step": 10, "default": 150,
            },
        ],
        "steps": [],
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


# ── Plotting tasks (shared field fragments) ───────────────────────────────
_FID_PLOT_LINES = {"type": "am-pf", "key": "lines"}

# Line selector — auto-injected by the modal when >1 line is loaded (see
# _prepare_dynamic_fields). Lets the user pick the profile(s) for a
# pseudo-section, and drives the dependent station list.
_PLOT_FIELD_LINES = {
    "id": _FID_PLOT_LINES,
    "key": "lines",
    "label": "Line(s)",
    "type": "multiselect",
    "options": [],
    "default": [],
    "placeholder": "all lines",
    "help": (
        "Survey line(s) / profile(s) to use. Leave empty for all lines; "
        "picking lines also narrows the station list below."
    ),
}
_PLOT_FIELD_STATIONS = {
    "id": _FID_PLOT_STATIONS,
    "key": "stations",
    "label": "Stations",
    "type": "multiselect",
    "options": [],
    "default": [],
    "placeholder": "all stations",
    "help": (
        "Stations to include. Leave empty for every station (in the "
        "selected line(s))."
    ),
}
_PLOT_FIELD_PMIN = {
    "id": _FID_PLOT_PMIN, "key": "period_min",
    "label": "Min period (s)", "type": "number",
    "min": 1e-6, "max": 1e5, "step": None, "default": None,
    "help": "Optional lower period bound. Blank = no limit.",
}
_PLOT_FIELD_PMAX = {
    "id": _FID_PLOT_PMAX, "key": "period_max",
    "label": "Max period (s)", "type": "number",
    "min": 1e-6, "max": 1e5, "step": None, "default": None,
    "help": "Optional upper period bound. Blank = no limit.",
}
_PLOT_FIELD_PUB = {
    "id": _FID_PLOT_PUB, "key": "publication",
    "label": "Publication style", "type": "radio",
    "options": [
        {"label": "Standard", "value": "off"},
        {"label": "Publication", "value": "on"},
    ],
    "inline": True, "default": "off",
    "help": "Larger fonts, tighter layout and 300-dpi-ready styling.",
}
_PLOT_FIELD_COMPONENTS = {
    "id": _FID_PLOT_COMPONENT, "key": "components",
    "label": "Components", "type": "select",
    "options": [
        {"label": "xy + yx", "value": "xy,yx"},
        {"label": "xy", "value": "xy"},
        {"label": "yx", "value": "yx"},
    ],
    "default": "xy,yx",
    "help": "Off-diagonal impedance components to display.",
}

_SCHEMAS["rhophi"] = {
    "title": "Rho/Phi Sounding Curves",
    "icon": "bi-graph-up",
    "color": "var(--blue)",
    "desc": (
        "Apparent-resistivity and phase versus period, per station. "
        "Pick the stations and components to view."
    ),
    "fields": [
        _PLOT_FIELD_STATIONS,
        dict(_PLOT_FIELD_COMPONENTS,
             options=[
                 {"label": "xy + yx", "value": "xy,yx"},
                 {"label": "xy", "value": "xy"},
                 {"label": "yx", "value": "yx"},
                 {"label": "determinant", "value": "det"},
             ]),
        {
            "id": _FID_PLOT_ERRBAR, "key": "errorbar",
            "label": "Error bars", "type": "radio",
            "options": [
                {"label": "On", "value": "on"},
                {"label": "Off", "value": "off"},
            ],
            "inline": True, "default": "on",
            "help": "Show measurement error bars when available.",
        },
        _PLOT_FIELD_PMIN,
        _PLOT_FIELD_PMAX,
        _PLOT_FIELD_PUB,
    ],
    "steps": [],
}

_SCHEMAS["phase_psection"] = {
    "title": "Phase Pseudo-section",
    "icon": "bi-grid-3x3",
    "color": "var(--mauve)",
    "desc": (
        "Scalar phase (deg) as a station-versus-period pseudo-section. "
        "One panel per component."
    ),
    "fields": [
        _PLOT_FIELD_COMPONENTS,
        _PLOT_FIELD_STATIONS,
        _PLOT_FIELD_PMIN,
        _PLOT_FIELD_PMAX,
        _PLOT_FIELD_PUB,
    ],
    "steps": [],
}

_SCHEMAS["pt_psection"] = {
    "title": "Phase-tensor (Φ) Pseudo-section",
    "icon": "bi-circle-half",
    "color": "var(--green)",
    "desc": (
        "Phase-tensor ellipses per station and period, coloured by an "
        "invariant. Choose the colour mapping and ellipse scale."
    ),
    "fields": [
        _PLOT_FIELD_STATIONS,
        {
            "id": _FID_PLOT_COLORBY, "key": "color_by",
            "label": "Colour by", "type": "select",
            "options": [
                {"label": "Skew (β)", "value": "skew"},
                {"label": "φ max", "value": "phimax"},
                {"label": "φ min", "value": "phimin"},
                {"label": "Ellipticity", "value": "ellipticity"},
            ],
            "default": "skew",
            "help": "Phase-tensor invariant used for the fill colour.",
        },
        {
            "id": _FID_PLOT_SCALE, "key": "scale",
            "label": "Ellipse scale", "type": "number",
            "min": 0.1, "max": 100.0, "step": None, "default": None,
            "help": "Optional ellipse size factor. Blank = auto.",
        },
        _PLOT_FIELD_PMIN,
        _PLOT_FIELD_PMAX,
        _PLOT_FIELD_PUB,
    ],
    "steps": [],
}

_SCHEMAS["strike"] = {
    "title": "Strike Analyzer",
    "icon": "bi-compass",
    "color": "var(--peach)",
    "desc": (
        "Estimate geoelectric strike per station (with a rose/analysis "
        "figure). Note the inherent 90° ambiguity."
    ),
    "fields": [
        {
            "id": _FID_PLOT_METHOD, "key": "method",
            "label": "Method", "type": "select",
            "options": [
                {"label": "Consensus", "value": "consensus"},
                {"label": "Impedance sweep", "value": "sweep"},
                {"label": "Phase tensor", "value": "pt"},
            ],
            "default": "consensus",
            "help": "Strike estimator.",
        },
        _PLOT_FIELD_STATIONS,
        _PLOT_FIELD_PMIN,
        _PLOT_FIELD_PMAX,
    ],
    "steps": [],
}

_SCHEMAS["dimensionality"] = {
    "title": "Dimensionality Classifier",
    "icon": "bi-diagram-2",
    "color": "var(--mauve)",
    "desc": (
        "Classify each station × period as 1-D / 2-D / 3-D from phase-tensor "
        "skew and ellipticity, with a pseudo-section figure."
    ),
    "fields": [
        {
            "id": _FID_SKEW_TH, "key": "skew_th",
            "label": "Skew threshold (°)", "type": "number",
            "min": 0.1, "max": 30.0, "step": None, "default": 3.0,
            "help": "|β| above this flags 3-D.",
        },
        {
            "id": _FID_ELLIPT_TH, "key": "ellipt_th",
            "label": "Ellipticity threshold", "type": "number",
            "min": 0.01, "max": 1.0, "step": None, "default": 0.2,
            "help": "Below this is treated as 1-D.",
        },
        _PLOT_FIELD_STATIONS,
    ],
    "steps": [],
}

_SCHEMAS["validator"] = {
    "title": "EDI Validator",
    "icon": "bi-check2-square",
    "color": "var(--green)",
    "desc": (
        "Per-station quality checklist — flags missing impedance, missing "
        "coordinates and low QC scores."
    ),
    "fields": [
        _PLOT_FIELD_STATIONS,
    ],
    "steps": [],
}

# ── Data / IO tools (Wave C) ───────────────────────────────────────────────
_SCHEMAS["coords"] = {
    "title": "Coordinate Transformer",
    "icon": "bi-pin-map",
    "color": "var(--sapphire)",
    "desc": (
        "Transform each station's geographic latitude/longitude to UTM "
        "easting/northing. The zone is auto-detected unless overridden."
    ),
    "fields": [
        {
            "id": _FID_DATUM, "key": "datum",
            "label": "Datum", "type": "select",
            "options": [
                {"label": "WGS84", "value": "WGS84"},
                {"label": "NAD83", "value": "NAD83"},
                {"label": "GRS80", "value": "GRS80"},
            ],
            "default": "WGS84",
            "help": "Geodetic datum for the UTM projection.",
        },
        {
            "id": _FID_ZONE, "key": "zone",
            "label": "UTM zone (0 = auto)", "type": "number",
            "min": 0, "max": 60, "step": 1, "default": 0,
            "help": "Force a UTM zone, or leave 0 to detect from longitude.",
        },
        {
            "id": _FID_HEMISPHERE, "key": "hemisphere",
            "label": "Hemisphere", "type": "radio",
            "options": [
                {"label": "North", "value": "N"},
                {"label": "South", "value": "S"},
            ],
            "inline": True, "default": "N",
            "help": "Northern or southern hemisphere.",
        },
        _PLOT_FIELD_STATIONS,
    ],
    "steps": [],
}

_SCHEMAS["elevation"] = {
    "title": "Elevation Enrichment",
    "icon": "bi-graph-up",
    "color": "var(--teal)",
    "desc": (
        "Fetch ground elevation for every station with coordinates from an "
        "open web service. Requires an internet connection."
    ),
    "fields": [
        {
            "id": _FID_API, "key": "api",
            "label": "Elevation service", "type": "radio",
            "options": [
                {"label": "Open-Meteo", "value": "open_meteo"},
                {"label": "Open-Topo-Data", "value": "open_topo_data"},
            ],
            "inline": True, "default": "open_meteo",
            "help": "Public elevation API queried over the network.",
        },
        _PLOT_FIELD_STATIONS,
    ],
    "steps": [],
}

_SCHEMAS["converter"] = {
    "title": "Format Converter",
    "icon": "bi-arrow-left-right",
    "color": "var(--yellow)",
    "desc": (
        "Re-export the loaded survey. CSV/JSON write per-station metadata; "
        "EDI re-writes one .edi per station. Files are saved to the folder "
        "you choose below."
    ),
    "fields": [
        {
            "id": _FID_FORMAT, "key": "format",
            "label": "Output format", "type": "select",
            "options": [
                {"label": "CSV  (station metadata)", "value": "csv"},
                {"label": "JSON (station metadata)", "value": "json"},
                {"label": "EDI  (re-export)", "value": "edi"},
            ],
            "default": "csv",
            "help": "What to write for each station.",
        },
        {
            "id": _FID_OUTDIR, "key": "output_dir",
            "label": "Output folder", "type": "text",
            "default": "",
            "placeholder": "blank = ~/pycsamt_export",
            "help": "Destination folder (created if missing).",
        },
        _PLOT_FIELD_STATIONS,
    ],
    "steps": [],
}

_SCHEMAS["batch_export"] = {
    "title": "Batch Plot Export",
    "icon": "bi-images",
    "color": "var(--mauve)",
    "desc": (
        "Render a bundle of standard plots from the loaded data and save "
        "them to a folder at the chosen format and resolution."
    ),
    "fields": [
        {
            "id": _FID_PLOTS, "key": "plots",
            "label": "Plot bundle", "type": "select",
            "options": [
                {"label": "Overview (ρ/φ + phase + Φ section)",
                 "value": "overview"},
                {"label": "Phase tensor (Φ section + map)",
                 "value": "phase_tensor"},
                {"label": "All of the above", "value": "all"},
                {"label": "ρ/φ sounding curves", "value": "rhophi"},
                {"label": "Phase pseudo-section", "value": "phase_psection"},
                {"label": "Phase-tensor (Φ) pseudo-section",
                 "value": "pt_psection"},
                {"label": "Phase-tensor map", "value": "pt_map"},
            ],
            "default": "overview",
            "help": "Which standard plots to render and save.",
        },
        {
            "id": _FID_FORMAT, "key": "format",
            "label": "Image format", "type": "select",
            "options": [
                {"label": "PNG", "value": "png"},
                {"label": "PDF", "value": "pdf"},
                {"label": "SVG", "value": "svg"},
                {"label": "TIFF", "value": "tiff"},
            ],
            "default": "png",
            "help": "Saved figure format.",
        },
        {
            "id": _FID_DPI, "key": "dpi",
            "label": "Resolution (dpi)", "type": "number",
            "min": 72, "max": 600, "step": 1, "default": 150,
            "help": "Raster resolution (72–600).",
        },
        {
            "id": _FID_OUTDIR, "key": "output_dir",
            "label": "Output folder", "type": "text",
            "default": "",
            "placeholder": "blank = ~/pycsamt_figures",
            "help": "Destination folder (created if missing).",
        },
    ],
    "steps": [],
}

# ── Stateful tools (Wave D) ────────────────────────────────────────────────
_SCHEMAS["freq_editor"] = {
    "title": "Frequency Editor",
    "icon": "bi-sliders",
    "color": "var(--red)",
    "desc": (
        "Confidence-based frequency QC. Edits run out-of-place; afterwards "
        "you can apply the edited survey to the session or export it."
    ),
    "fields": [
        {
            "id": _FID_MODE, "key": "mode",
            "label": "Mode", "type": "radio",
            "options": [
                {"label": "Recover", "value": "recover"},
                {"label": "Drop", "value": "drop"},
                {"label": "Mask", "value": "mask"},
            ],
            "inline": True, "default": "recover",
            "help": "Recover interpolates weak periods; drop/mask remove them.",
        },
        {
            "id": _FID_PLOT_METHOD, "key": "method",
            "label": "Confidence method", "type": "select",
            "options": [
                {"label": "Composite", "value": "composite"},
                {"label": "SNR", "value": "snr"},
                {"label": "Phase slope", "value": "phase_slope"},
                {"label": "Coherence", "value": "coherence"},
            ],
            "default": "composite",
            "help": "How per-period confidence is scored.",
        },
        {
            "id": _FID_THRESHOLD, "key": "threshold",
            "label": "Threshold", "type": "number",
            "min": 0.0, "max": 1.0, "step": None, "default": 0.5,
            "help": "Periods below this confidence are acted on.",
        },
        {
            "id": _FID_CI_HI, "key": "ci_hi",
            "label": "CI high", "type": "number",
            "min": 0.0, "max": 1.0, "step": None, "default": 0.9,
            "help": "Upper confidence band.",
        },
        {
            "id": _FID_CI_LO, "key": "ci_lo",
            "label": "CI low", "type": "number",
            "min": 0.0, "max": 1.0, "step": None, "default": 0.5,
            "help": "Lower confidence band.",
        },
        {
            "id": _FID_ALSO, "key": "also",
            "label": "Apply to", "type": "select",
            "options": [
                {"label": "Z + tipper", "value": "both"},
                {"label": "Impedance (Z)", "value": "z"},
                {"label": "Tipper", "value": "tipper"},
            ],
            "default": "both",
            "help": "Which transfer functions to edit.",
        },
        {
            "id": _FID_REJECT, "key": "reject",
            "label": "Reject rows", "type": "radio",
            "options": [
                {"label": "Drop", "value": "drop"},
                {"label": "Mask", "value": "mask"},
                {"label": "Keep", "value": "keep"},
            ],
            "inline": True, "default": "drop",
            "help": "What to do with rejected periods.",
        },
        _PLOT_FIELD_STATIONS,
    ],
    "steps": [],
}

_SCHEMAS["layered_model"] = {
    "title": "Layered Model Builder",
    "icon": "bi-stack",
    "color": "var(--mauve)",
    "desc": (
        "Build and preview a synthetic 1-D layered-earth resistivity model. "
        "Use a preset, or enter resistivities and thicknesses directly. "
        "No loaded data required."
    ),
    "fields": [
        {
            "id": _FID_PRESET, "key": "preset",
            "label": "Source", "type": "select",
            "options": [
                {"label": "Custom (enter values below)", "value": "custom"},
                {"label": "Random preset", "value": "random"},
                {"label": "Blocky preset", "value": "blocky"},
                {"label": "Smooth preset", "value": "smooth"},
            ],
            "default": "custom",
            "help": "Custom uses the fields below; presets are auto-generated.",
        },
        {
            "id": _FID_N_LAYERS, "key": "n_layers",
            "label": "Number of layers (presets)", "type": "slider",
            "min": 2, "max": 20, "step": 1, "default": 3,
            "help": "Layer count for the random/blocky/smooth presets.",
        },
        {
            "id": _FID_RESISTIVITIES, "key": "resistivities",
            "label": "Resistivities (Ω·m)", "type": "text",
            "default": "",
            "placeholder": "e.g. 100, 10, 500  (top → halfspace)",
            "help": "Comma-separated, one per layer. Blank = 100, 10, 500.",
        },
        {
            "id": _FID_THICKNESSES, "key": "thicknesses",
            "label": "Thicknesses (m)", "type": "text",
            "default": "",
            "placeholder": "e.g. 300, 800  (one fewer than layers)",
            "help": "Comma-separated; one fewer than resistivities (no "
                    "halfspace thickness). Blank = 300, 800.",
        },
        {
            "id": _FID_DEPTH_MAX, "key": "depth_max",
            "label": "Max depth (m, random preset)", "type": "number",
            "min": 100, "max": 50000, "step": 100, "default": 2000.0,
            "help": "Total depth spanned by the random preset.",
        },
    ],
    "steps": [],
}

_SCHEMAS["station_response"] = {
    "title": "Station Response Inspector",
    "icon": "bi-activity",
    "color": "var(--blue)",
    "desc": (
        "Per-station impedance response (apparent resistivity & phase with "
        "the Bode-predicted phase) for distortion checks."
    ),
    "fields": [
        dict(_PLOT_FIELD_STATIONS,
             label="Station",
             help="Station to inspect. Blank = first station."),
        dict(_PLOT_FIELD_COMPONENTS,
             options=[
                 {"label": "xy + yx", "value": "xy,yx"},
                 {"label": "xy", "value": "xy"},
                 {"label": "yx", "value": "yx"},
             ]),
        _PLOT_FIELD_PMIN,
        _PLOT_FIELD_PMAX,
        _PLOT_FIELD_PUB,
    ],
    "steps": [],
}

_SCHEMAS["strike_profile"] = {
    "title": "Strike Profile Viewer",
    "icon": "bi-bar-chart-steps",
    "color": "var(--peach)",
    "desc": (
        "Geoelectric strike angle versus station position along the line, "
        "with an inter-quartile ribbon."
    ),
    "fields": [
        {
            "id": _FID_PLOT_METHOD, "key": "method",
            "label": "Strike method", "type": "select",
            "options": [
                {"label": "Consensus", "value": "consensus"},
                {"label": "Impedance sweep", "value": "sweep"},
                {"label": "Phase tensor", "value": "pt"},
            ],
            "default": "consensus",
            "help": "Strike estimator to drive the profile.",
        },
        {
            "id": _FID_PLOT_SORTBY, "key": "sort_by",
            "label": "Order stations by", "type": "select",
            "options": [
                {"label": "Auto", "value": "auto"},
                {"label": "Longitude", "value": "lon"},
                {"label": "Latitude", "value": "lat"},
                {"label": "Name", "value": "name"},
            ],
            "default": "auto",
            "help": "Horizontal ordering of stations along the profile.",
        },
        _PLOT_FIELD_STATIONS,
        _PLOT_FIELD_PMIN,
        _PLOT_FIELD_PMAX,
        _PLOT_FIELD_PUB,
    ],
    "steps": [],
}

_SCHEMAS["phase_tensor_map"] = {
    "title": "Phase-tensor Map",
    "icon": "bi-geo-alt",
    "color": "var(--green)",
    "desc": (
        "Geographic map of phase-tensor ellipses at one period, coloured "
        "by an invariant (tipper arrows overlaid when available)."
    ),
    "fields": [
        {
            "id": _FID_PLOT_PERIOD, "key": "period",
            "label": "Period (s)", "type": "number",
            "min": 1e-6, "max": 1e5, "step": None, "default": 1.0,
            "help": "Period at which to draw the ellipse map.",
        },
        {
            "id": _FID_PLOT_COLORBY, "key": "color_by",
            "label": "Colour by", "type": "select",
            "options": [
                {"label": "Skew (β)", "value": "skew"},
                {"label": "φ max", "value": "phimax"},
                {"label": "φ min", "value": "phimin"},
                {"label": "Ellipticity", "value": "ellipticity"},
            ],
            "default": "skew",
            "help": "Phase-tensor invariant used for the fill colour.",
        },
        {
            "id": _FID_PLOT_SCALE, "key": "scale",
            "label": "Ellipse scale", "type": "number",
            "min": 0.1, "max": 100.0, "step": None, "default": None,
            "help": "Optional ellipse size factor. Blank = auto.",
        },
        _PLOT_FIELD_STATIONS,
        _PLOT_FIELD_PUB,
    ],
    "steps": [],
}

_SCHEMAS["tipper_plot"] = {
    "title": "Tipper Plot",
    "icon": "bi-arrow-up-right",
    "color": "var(--peach)",
    "desc": (
        "Vertical-field tipper: component curves (Tx/Ty) or induction "
        "arrows. Only available when the data contains a tipper."
    ),
    "fields": [
        {
            "id": _FID_PLOT_VIEW, "key": "view",
            "label": "View", "type": "radio",
            "options": [
                {"label": "Components (Tx/Ty)", "value": "components"},
                {"label": "Induction arrows", "value": "arrows"},
            ],
            "inline": True, "default": "components",
            "help": "Component curves vs period, or a map of induction arrows.",
        },
        _PLOT_FIELD_STATIONS,
        {
            "id": _FID_PLOT_PARTS, "key": "parts",
            "label": "Parts (components view)", "type": "select",
            "options": [
                {"label": "real + imaginary", "value": "real,imag"},
                {"label": "real", "value": "real"},
                {"label": "imaginary", "value": "imag"},
            ],
            "default": "real,imag",
            "help": "Real / imaginary tipper parts to draw.",
        },
        {
            "id": _FID_PLOT_CONV, "key": "convention",
            "label": "Arrow convention (arrows view)", "type": "select",
            "options": [
                {"label": "Parkinson", "value": "park"},
                {"label": "Wiese", "value": "wiese"},
            ],
            "default": "park",
            "help": "Sign convention for induction arrows.",
        },
        {
            "id": _FID_PLOT_PERIOD, "key": "period",
            "label": "Period for arrows (s)", "type": "number",
            "min": 1e-6, "max": 1e5, "step": None, "default": 1.0,
            "help": "Period at which to draw induction arrows.",
        },
        _PLOT_FIELD_PUB,
    ],
    "steps": [],
}


# ── correction-method schemas (generated from the catalogue) ──────────────────
# Each correction workflow's parameter modal is built directly from its
# CorrectionController ParamSpec list, so "different correction with different
# parameter control" needs no hand-written schema per method.

def _corr_field_from_spec(ps) -> dict:
    """Map a CorrectionController ``ParamSpec`` to a param-modal field dict."""
    fid = {"type": "am-pf", "key": ps.name}
    field = {
        "id": fid, "key": ps.name, "label": ps.label,
        "help": getattr(ps, "tip", "") or "",
    }
    kind = ps.kind
    if kind == "combo":
        field.update(
            type="select",
            options=[{"label": str(o), "value": o} for o in (ps.opts or [])],
            default=ps.default,
        )
    elif kind == "check":
        field.update(
            type="radio", inline=True,
            options=[{"label": "Yes", "value": True},
                     {"label": "No", "value": False}],
            default=bool(ps.default),
        )
    else:  # "spin" (int) / "dspin" (float)
        lo, hi, st = (tuple(ps.opts) + (None, None, None))[:3] if ps.opts \
            else (None, None, None)
        field.update(
            type="number", min=lo, max=hi,
            step=(int(st) if (kind == "spin" and st) else st),
            default=ps.default,
        )
    return field


def _build_correction_schemas() -> None:
    try:
        from pycsamt.agents._corrections import (
            CORRECTION_METHODS, param_specs, method_desc,
        )
    except Exception:  # noqa: BLE001 — corrections optional
        return
    for wf, meta in CORRECTION_METHODS.items():
        try:
            specs = param_specs(wf)
        except Exception:  # noqa: BLE001
            continue
        fields = [_corr_field_from_spec(ps) for ps in specs]
        fields.append(_PLOT_FIELD_STATIONS)  # optional station targeting
        _SCHEMAS[wf] = {
            "title": meta.get("title", wf),
            "icon":  meta.get("icon", "bi-sliders"),
            "color": meta.get("color", "var(--blue)"),
            "desc":  method_desc(wf) or meta.get("title", wf),
            "fields": fields,
            "steps": [],
        }


_build_correction_schemas()


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
    if ft == "multiselect":
        # tolerate a stale scalar value from older localStorage
        mv = val if isinstance(val, list) else ([val] if val else [])
        return dcc.Dropdown(
            id=fid,
            options=f.get("options", []),
            value=mv,
            multi=True,
            placeholder=f.get("placeholder", ""),
            className="am-pf-select am-pf-multi",
        )
    if ft == "radio":
        return dbc.RadioItems(
            id=fid,
            options=f["options"],
            value=val,
            inline=f.get("inline", False),
            className="am-pf-radio",
        )
    if ft == "text":
        return dbc.Input(
            id=fid,
            type="text",
            value=val or "",
            placeholder=f.get(
                "placeholder", ""
            ),
            size="sm",
            debounce=True,
            className="am-pf-text",
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


# ── Dynamic line / station options ────────────

def _line_station_options(groups, edi_path):
    """Return ``(line_options, station_options, {line: [stations]})`` for the
    param modal's dependent dropdowns.

    Derived from the loaded line groups (or by scanning the EDI path when
    ungrouped). Station names are the EDI file stems, so no EDI parsing is
    needed."""
    from pathlib import Path
    line_to_st: dict[str, list[str]] = {}
    if groups:
        for ln, files in groups.items():
            line_to_st[str(ln)] = sorted(
                {Path(str(f)).stem for f in (files or [])}
            )
    elif edi_path:
        try:
            p = Path(str(edi_path))
            if p.is_dir():
                stems = (sorted({f.stem for f in p.rglob("*.edi")})
                         or sorted({f.stem for f in p.rglob("*.EDI")}))
                if stems:
                    line_to_st["(all)"] = stems
            elif p.is_file():
                line_to_st["(all)"] = [p.stem]
        except Exception:  # noqa: BLE001
            pass
    line_opts = [
        {"label": f"{ln}  ({len(sts)})", "value": ln}
        for ln, sts in line_to_st.items()
    ]
    all_st = sorted({s for sts in line_to_st.values() for s in sts})
    station_opts = [{"label": s, "value": s} for s in all_st]
    return line_opts, station_opts, line_to_st


def _station_options_for_lines(lines, line_to_st):
    """Station dropdown options for the selected *lines* (all when empty)."""
    line_to_st = line_to_st or {}
    if lines:
        sts = sorted({s for ln in lines for s in line_to_st.get(ln, [])})
    else:
        sts = sorted({s for vals in line_to_st.values() for s in vals})
    return [{"label": s, "value": s} for s in sts]


def _prepare_dynamic_fields(fields, groups, edi_path, preselect):
    """Inject runtime options into the station field and auto-prepend a line
    selector when the form filters stations and more than one line is loaded.

    Returns ``(fields, {line: [stations]})``."""
    line_opts, station_opts, line_to_st = _line_station_options(
        groups, edi_path
    )
    has_stations = any(f.get("key") == "stations" for f in fields)
    out: list = []
    if has_stations and len(line_opts) >= 2:
        out.append(dict(
            _PLOT_FIELD_LINES,
            options=line_opts,
            default=[ln for ln in (preselect or []) if ln in line_to_st],
        ))
    for f in fields:
        if f.get("key") == "stations" and f.get("type") == "multiselect":
            f = dict(f, options=station_opts)
        out.append(f)
    return out, line_to_st


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
        Output(IDs.STORE_LINE_STATIONS, "data"),
        Input(IDs.STORE_PENDING, "data"),
        State(IDs.STORE_INV_CONFIG, "data"),
        State(IDs.STORE_EDI, "data"),
        prevent_initial_call=True,
    )
    def _open_param_modal(pending, inv_config, edi_store):
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
        # Loaded line groups + station names drive the line/station selectors.
        _groups = (
            (edi_store or {}).get("groups", {})
            or pending.get("edi_groups", {}) or {}
        )
        _edi_path = (
            (edi_store or {}).get("path", "")
            or pending.get("edi_path", "")
        )
        _presel = pending.get("selected_lines", []) or []
        form: list = []
        line_to_st: dict = {}
        wf_fields = schema.get("fields", [])
        if wf_fields:
            prepared, line_to_st = _prepare_dynamic_fields(
                wf_fields, _groups, _edi_path, _presel
            )
            form.extend(_render_form(prepared, ic))
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
            True, title, schema["desc"], form, line_to_st
        )

    # Dependent dropdown: selected line(s) → station options.
    @app.callback(
        Output({"type": "am-pf", "key": "stations"}, "options"),
        Input({"type": "am-pf", "key": "lines"}, "value"),
        State(IDs.STORE_LINE_STATIONS, "data"),
        prevent_initial_call=True,
    )
    def _stations_for_lines(lines, line_to_st):
        return _station_options_for_lines(lines, line_to_st)

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

        # Apply line filter if user had pre-
        # selected line(s) before the param modal.
        # Fall back to edi_path snapshotted in
        # pending if STORE_EDI is now empty.
        jid = _new_job()
        _edi_use = dict(edi_store or {})
        if not _edi_use.get("path"):
            # STORE_EDI lost — recover from pending
            _fb_path = pending.get("edi_path", "")
            if _fb_path:
                _edi_use["path"] = _fb_path
            _fb_grp = pending.get("edi_groups", {})
            if _fb_grp:
                _edi_use["groups"] = _fb_grp
        # Line selection made *in the modal* wins, else any pre-selection
        # carried over from the line picker.
        _modal_lines = pf_vals.get("lines") or []
        _sel = _modal_lines or pending.get("selected_lines", [])
        if _sel:
            _edi_use["selected_lines"] = _sel
        threading.Thread(
            target=_run_agent,
            args=(
                jid,
                pending.get("text", ""),
                _edi_use,
                settings or {},
                new_ic,
            ),
            daemon=True,
        ).start()

        # Persist inversion PARAMS only — never the transient workflow
        # selection. STORE_INV_CONFIG is localStorage; persisting
        # "workflow" here would force every later request through this
        # workflow (e.g. all requests routed to ai_inversion). The
        # current run still gets its workflow via new_ic passed to the
        # thread above.
        persist_ic = {
            k: v for k, v in new_ic.items()
            if k not in ("workflow", "lines", "stations")
        }
        return (
            msgs,
            {"jid": jid},
            False,
            stored_msgs or [],
            persist_ic,
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
