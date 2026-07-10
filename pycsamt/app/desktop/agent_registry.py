# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
agent_registry.py — metadata catalogue for all pycsamt agents and
processing operations exposed through the RunAgentDialog / AgentWorker.

No Qt imports — importable in tests and the Dash web app.

Schema for each entry
---------------------
type          'llm' | 'processing'
category      (llm only) logical group label shown in RunAgentDialog
description   one-line tooltip shown in RunAgentDialog
class_name    (llm only)  name of the class in pycsamt.agents
fn_name       (processing) name of the function in pycsamt.emtools
params        ordered dict of parameter descriptors (see below)
result_plot   (optional) emtools function name to plot the result

Parameter descriptor keys
--------------------------
type      'combo' | 'float' | 'int' | 'bool' | 'str'
default   default value
options   (combo only) list of strings
range     (float/int) (min, max) inclusive
step      (float/int) step for spinbox
label     display label (defaults to the key name)
tip       optional tooltip
"""

from __future__ import annotations

from typing import Any

# ──────────────────────────────────────────────────────────────────────────────

ParamSpec  = dict[str, Any]
AgentEntry = dict[str, Any]

AGENT_REGISTRY: dict[str, AgentEntry] = {

    # ══════════════════════════════════════════════════════════════════════════
    # LLM agents (require API key)
    # ══════════════════════════════════════════════════════════════════════════

    # ── Data Loading ──────────────────────────────────────────────────────────

    "MT Loader": {
        "type": "llm",
        "category": "Data Loading",
        "class_name": "MTLoaderAgent",
        "description": "Load EDI / AVG / J files into a Sites object with a per-station QC report.",
        "params": {
            "recursive": {
                "type": "bool",
                "default": True,
                "label": "Recursive search",
                "tip": "Search subdirectories recursively for data files.",
            },
            "on_dup": {
                "type": "combo",
                "options": ["replace", "skip", "error"],
                "default": "replace",
                "label": "Duplicate handling",
                "tip": "What to do when a station name appears more than once.",
            },
        },
        "result_plot": None,
    },

    "Context Parser": {
        "type": "llm",
        "category": "Data Loading",
        "class_name": "ContextInputAgent",
        "description": "Parse a natural-language MT workflow request into a structured config.",
        "params": {},
        "result_plot": None,
    },

    # ── Quality Control ───────────────────────────────────────────────────────

    "QC Analysis": {
        "type": "llm",
        "category": "Quality Control",
        "class_name": "DataQCAgent",
        "description": "Quality-control scoring via LLM interpretation of SNR, skew and coverage metrics.",
        "params": {
            "method": {
                "type": "combo",
                "options": ["composite", "snr", "skew"],
                "default": "composite",
                "label": "QC method",
                "tip": "Composite combines SNR and skew; others score one criterion only.",
            },
            "min_frac_ok": {
                "type": "float",
                "default": 0.6,
                "range": (0.0, 1.0),
                "step": 0.05,
                "label": "Min fraction OK",
                "tip": "Minimum fraction of frequencies that must pass QC (0–1).",
            },
            "min_snr_med": {
                "type": "float",
                "default": 2.0,
                "range": (0.1, 20.0),
                "step": 0.5,
                "label": "Min median SNR",
            },
            "max_skew_med": {
                "type": "float",
                "default": 6.0,
                "range": (0.5, 30.0),
                "step": 0.5,
                "label": "Max median skew (°)",
            },
        },
        "result_plot": "plot_qc_quicklook",
    },

    "Anomaly Detection": {
        "type": "llm",
        "category": "Quality Control",
        "class_name": "AnomalyDetectionAgent",
        "description": "Unsupervised CAE anomaly flagging per (station, frequency).",
        "params": {
            "threshold_percentile": {
                "type": "float",
                "default": 95.0,
                "range": (50.0, 99.9),
                "step": 1.0,
                "label": "Anomaly threshold (%ile)",
                "tip": "Reconstruction-error percentile above which a sample is flagged.",
            },
            "latent_dim": {
                "type": "int",
                "default": 32,
                "range": (4, 256),
                "step": 4,
                "label": "CAE latent dim",
            },
            "epochs": {
                "type": "int",
                "default": 50,
                "range": (5, 500),
                "step": 5,
                "label": "Training epochs",
            },
        },
        "result_plot": None,
    },

    # ── Pre-processing ─────────────────────────────────────────────────────────

    "Static Shift": {
        "type": "llm",
        "category": "Pre-processing",
        "class_name": "StaticShiftAgent",
        "description": "Detect and correct static-shift distortion using spatial averaging with LLM review.",
        "params": {
            "method": {
                "type": "combo",
                "options": ["ama", "refmedian", "bilateral", "loess"],
                "default": "ama",
                "label": "Correction method",
                "tip": "AMA = azimuthal moving average; loess = robust local polynomial.",
            },
            "half_window": {
                "type": "int",
                "default": 3,
                "range": (1, 15),
                "step": 1,
                "label": "Half-window (stations)",
                "tip": "Number of stations on each side used for spatial averaging.",
            },
            "inplace": {
                "type": "bool",
                "default": False,
                "label": "Correct in-place",
                "tip": "If checked, modifies the loaded Sites object directly.",
            },
        },
        "result_plot": "plot_ss_station_curves",
    },

    "Tensor Rotation": {
        "type": "llm",
        "category": "Pre-processing",
        "class_name": "TensorRotationAgent",
        "description": "Rotate impedance tensors and tipper vectors by a fixed geoelectric strike angle.",
        "params": {
            "strike_deg": {
                "type": "float",
                "default": 0.0,
                "range": (-180.0, 180.0),
                "step": 1.0,
                "label": "Strike angle (°)",
                "tip": "Clockwise rotation angle from geographic north.",
            },
        },
        "result_plot": None,
    },

    "Denoising": {
        "type": "llm",
        "category": "Pre-processing",
        "class_name": "DenoisingAgent",
        "description": "RPCA / Hampel / EMAP / AI-CAE denoising of MT impedance data.",
        "params": {
            "method": {
                "type": "combo",
                "options": ["rpca", "hampel", "emap", "cae"],
                "default": "rpca",
                "label": "Denoising method",
                "tip": "RPCA = robust PCA; CAE = convolutional autoencoder.",
            },
            "rank": {
                "type": "int",
                "default": 2,
                "range": (1, 20),
                "step": 1,
                "label": "RPCA rank",
                "tip": "Low-rank approximation for RPCA (ignored for other methods).",
            },
            "half_window": {
                "type": "int",
                "default": 3,
                "range": (1, 15),
                "step": 1,
                "label": "Hampel half-window",
                "tip": "Window size for Hampel identifier (ignored for RPCA/CAE).",
            },
        },
        "result_plot": None,
    },

    "Frequency Decimation": {
        "type": "llm",
        "category": "Pre-processing",
        "class_name": "FrequencyDecimationAgent",
        "description": "Select optimal periods from MT data for inversion (logarithmic decimation + SNR filter).",
        "params": {
            "n_per_decade": {
                "type": "int",
                "default": 6,
                "range": (1, 20),
                "step": 1,
                "label": "Periods per decade",
            },
            "snr_threshold": {
                "type": "float",
                "default": 3.0,
                "range": (0.5, 20.0),
                "step": 0.5,
                "label": "Min SNR",
                "tip": "Periods with median SNR below this value are excluded.",
            },
            "component": {
                "type": "combo",
                "options": ["xy", "yx", "det"],
                "default": "xy",
                "label": "Impedance component",
            },
            "period_range": {
                "type": "str",
                "default": "",
                "label": "Period range (s)",
                "tip": "Comma-separated min,max (e.g. '1e-4,1e2'). Leave blank to use full range.",
            },
        },
        "result_plot": None,
    },

    # ── Phase & Dimensionality ─────────────────────────────────────────────────

    "Phase Analysis": {
        "type": "llm",
        "category": "Phase & Dimensionality",
        "class_name": "PhaseAnalysisAgent",
        "description": "Dimensionality analysis (skew + ellipticity) with LLM geological interpretation.",
        "params": {
            "skew_th": {
                "type": "float",
                "default": 5.0,
                "range": (0.5, 30.0),
                "step": 0.5,
                "label": "Skew threshold (°)",
            },
            "ellipt_th": {
                "type": "float",
                "default": 0.1,
                "range": (0.01, 1.0),
                "step": 0.01,
                "label": "Ellipticity threshold",
            },
        },
        "result_plot": "plot_dimensionality_psection",
    },

    "Tipper Analysis": {
        "type": "llm",
        "category": "Phase & Dimensionality",
        "class_name": "TipperAnalysisAgent",
        "description": "Analyse tipper vectors: induction arrows, strike, phase split.",
        "params": {
            "convention": {
                "type": "combo",
                "options": ["wiese", "parkinson"],
                "default": "wiese",
                "label": "Arrow convention",
                "tip": "Wiese: real arrows point away from conductors; Parkinson: toward.",
            },
            "use_imag": {
                "type": "bool",
                "default": False,
                "label": "Plot imaginary part",
            },
            "period_ref": {
                "type": "float",
                "default": 1.0,
                "range": (1e-5, 1e5),
                "step": 1.0,
                "label": "Reference period (s)",
                "tip": "Period used for single-frequency map. Leave default to use median.",
            },
        },
        "result_plot": None,
    },

    "Sensitivity": {
        "type": "llm",
        "category": "Phase & Dimensionality",
        "class_name": "SensitivityAgent",
        "description": "Bostick sensitivity kernels and vertical resolution analysis.",
        "params": {
            "component": {
                "type": "combo",
                "options": ["xy", "yx", "det"],
                "default": "xy",
                "label": "Impedance component",
            },
            "depth_max": {
                "type": "float",
                "default": 50000.0,
                "range": (100.0, 500000.0),
                "step": 1000.0,
                "label": "Max depth (m)",
                "tip": "Maximum depth for sensitivity kernel plot.",
            },
        },
        "result_plot": None,
    },

    # ── Forward Modelling ─────────────────────────────────────────────────────

    "Forward Model": {
        "type": "llm",
        "category": "Forward Modelling",
        "class_name": "ForwardModelAgent",
        "description": "Run a 1-D, 2-D, or 3-D MT forward model.",
        "params": {
            "dim": {
                "type": "combo",
                "options": ["1", "2", "3"],
                "default": "1",
                "label": "Dimension",
                "tip": "Dimensionality of the forward model.",
            },
        },
        "result_plot": None,
    },

    # ── Inversion Preparation ─────────────────────────────────────────────────

    "Inversion Prep": {
        "type": "llm",
        "category": "Inversion",
        "class_name": "InversionPrepAgent",
        "description": "Prepare MT data files for Occam2D, ModEM3D or MARE2DEM inversion codes.",
        "params": {
            "code": {
                "type": "combo",
                "options": ["occam2d", "modem", "mare2dem"],
                "default": "occam2d",
                "label": "Inversion code",
            },
            "error_floor": {
                "type": "float",
                "default": 0.05,
                "range": (0.001, 0.5),
                "step": 0.005,
                "label": "Error floor",
                "tip": "Minimum relative data error applied before writing input files.",
            },
        },
        "result_plot": None,
    },

    "Occam2D": {
        "type": "llm",
        "category": "Inversion",
        "class_name": "Occam2DAgent",
        "description": "Generate a complete Occam2D inversion input file set and optionally run the code.",
        "params": {
            "error_floor": {
                "type": "float",
                "default": 0.05,
                "range": (0.001, 0.5),
                "step": 0.005,
                "label": "Error floor",
            },
            "target_rms": {
                "type": "float",
                "default": 1.0,
                "range": (0.5, 10.0),
                "step": 0.1,
                "label": "Target RMS",
            },
            "modes": {
                "type": "str",
                "default": "TE,TM",
                "label": "Modes (comma-separated)",
                "tip": "e.g. 'TE,TM' or 'det'",
            },
        },
        "result_plot": None,
    },

    "ModEM": {
        "type": "llm",
        "category": "Inversion",
        "class_name": "ModEmAgent",
        "description": "Write a ModEM3D MT data file from EDI sources and optionally launch the code.",
        "params": {
            "error_floor": {
                "type": "float",
                "default": 0.05,
                "range": (0.001, 0.5),
                "step": 0.005,
                "label": "Error floor",
            },
            "component_types": {
                "type": "str",
                "default": "Full_Impedance",
                "label": "Component types",
                "tip": "Comma-separated ModEM component type strings.",
            },
        },
        "result_plot": None,
    },

    "MARE2DEM": {
        "type": "llm",
        "category": "Inversion",
        "class_name": "Mare2DEMAgent",
        "description": "Orchestrate a complete MARE2DEM 2.5-D EM inversion workflow.",
        "params": {
            "n_procs": {
                "type": "int",
                "default": 4,
                "range": (1, 128),
                "step": 1,
                "label": "MPI processes",
            },
            "use_mpi": {
                "type": "bool",
                "default": True,
                "label": "Use MPI",
            },
            "initial_rho": {
                "type": "float",
                "default": 1.0,
                "range": (0.001, 10000.0),
                "step": 1.0,
                "label": "Initial ρ (Ω·m)",
            },
            "target_rms": {
                "type": "float",
                "default": 1.0,
                "range": (0.5, 10.0),
                "step": 0.1,
                "label": "Target RMS",
            },
            "max_iterations": {
                "type": "int",
                "default": 150,
                "range": (10, 1000),
                "step": 10,
                "label": "Max iterations",
            },
        },
        "result_plot": None,
    },

    "Inversion Backend": {
        "type": "llm",
        "category": "Inversion",
        "class_name": "InversionBackendAgent",
        "description": "Drive pycsamt physics-based inversion backends (Occam, NLCG, L-BFGS, Marquardt).",
        "params": {
            "backend": {
                "type": "combo",
                "options": ["builtin", "occam2d", "modem", "mare2dem"],
                "default": "builtin",
                "label": "Backend",
            },
            "dimension": {
                "type": "combo",
                "options": ["1d", "2d", "3d"],
                "default": "1d",
                "label": "Dimension",
            },
            "method": {
                "type": "combo",
                "options": ["mt", "amt", "csamt", "tem"],
                "default": "mt",
                "label": "EM method",
            },
            "n_layers": {
                "type": "int",
                "default": 5,
                "range": (2, 50),
                "step": 1,
                "label": "Number of layers",
            },
            "max_iter": {
                "type": "int",
                "default": 80,
                "range": (5, 500),
                "step": 5,
                "label": "Max iterations",
            },
            "regularization": {
                "type": "combo",
                "options": ["smooth", "sharp", "tv", "sparse"],
                "default": "smooth",
                "label": "Regularization",
            },
            "error_floor": {
                "type": "float",
                "default": 0.05,
                "range": (0.001, 0.5),
                "step": 0.005,
                "label": "Error floor",
            },
        },
        "result_plot": None,
    },

    # ── AI Inversion ──────────────────────────────────────────────────────────

    "AI Inversion (1-D)": {
        "type": "llm",
        "category": "AI Inversion",
        "class_name": "AIInversionAgent",
        "description": "End-to-end 1-D AI inversion: train on synthetic data, predict on observed sites.",
        "params": {
            "arch": {
                "type": "combo",
                "options": ["resnet", "cnn", "fcn"],
                "default": "resnet",
                "label": "Architecture",
            },
            "n_layers": {
                "type": "int",
                "default": 5,
                "range": (2, 30),
                "step": 1,
                "label": "Output layers",
                "tip": "Number of 1-D earth model layers predicted.",
            },
            "n_train_samples": {
                "type": "int",
                "default": 2000,
                "range": (100, 50000),
                "step": 100,
                "label": "Training samples",
            },
            "epochs": {
                "type": "int",
                "default": 30,
                "range": (1, 500),
                "step": 5,
                "label": "Training epochs",
            },
            "pretrained": {
                "type": "str",
                "default": "",
                "label": "Pretrained model path",
                "tip": "Leave blank to train from scratch.",
            },
        },
        "result_plot": None,
    },

    "AI Inversion (2-D)": {
        "type": "llm",
        "category": "AI Inversion",
        "class_name": "Inv2DAgent",
        "description": "U-Net 2-D profile inversion with lateral continuity regularisation.",
        "params": {
            "arch": {
                "type": "combo",
                "options": ["unet", "unet++"],
                "default": "unet",
                "label": "Architecture",
            },
            "n_depth": {
                "type": "int",
                "default": 40,
                "range": (10, 200),
                "step": 5,
                "label": "Depth cells",
            },
            "n_freqs": {
                "type": "int",
                "default": 32,
                "range": (8, 128),
                "step": 4,
                "label": "Frequency channels",
            },
            "n_train_profiles": {
                "type": "int",
                "default": 200,
                "range": (10, 2000),
                "step": 10,
                "label": "Training profiles",
            },
            "epochs": {
                "type": "int",
                "default": 30,
                "range": (1, 500),
                "step": 5,
                "label": "Training epochs",
            },
        },
        "result_plot": None,
    },

    "AI Inversion (3-D)": {
        "type": "llm",
        "category": "AI Inversion",
        "class_name": "Inv3DAgent",
        "description": "GCN 3-D spatial inversion using inter-station graph message-passing.",
        "params": {
            "n_layers": {
                "type": "int",
                "default": 5,
                "range": (2, 30),
                "step": 1,
                "label": "Output depth layers",
            },
            "n_freqs": {
                "type": "int",
                "default": 32,
                "range": (8, 128),
                "step": 4,
                "label": "Frequency channels",
            },
            "n_train_profiles": {
                "type": "int",
                "default": 150,
                "range": (10, 2000),
                "step": 10,
                "label": "Training profiles",
            },
            "epochs": {
                "type": "int",
                "default": 30,
                "range": (1, 500),
                "step": 5,
                "label": "Training epochs",
            },
            "radius": {
                "type": "float",
                "default": 5000.0,
                "range": (100.0, 100000.0),
                "step": 500.0,
                "label": "Graph edge radius (m)",
                "tip": "Stations within this distance are connected in the GCN graph.",
            },
        },
        "result_plot": None,
    },

    "Ensemble Inversion": {
        "type": "llm",
        "category": "AI Inversion",
        "class_name": "EnsembleAgent",
        "description": "Ensemble 1-D MT inversion with conformal uncertainty quantification.",
        "params": {
            "n_estimators": {
                "type": "int",
                "default": 5,
                "range": (2, 50),
                "step": 1,
                "label": "Ensemble size",
            },
            "arch": {
                "type": "combo",
                "options": ["resnet", "cnn", "fcn"],
                "default": "resnet",
                "label": "Architecture",
            },
            "n_layers": {
                "type": "int",
                "default": 5,
                "range": (2, 30),
                "step": 1,
                "label": "Output layers",
            },
            "epochs": {
                "type": "int",
                "default": 30,
                "range": (1, 500),
                "step": 5,
                "label": "Training epochs",
            },
            "calibrate": {
                "type": "bool",
                "default": True,
                "label": "Conformal calibration",
                "tip": "Apply conformal prediction to calibrate uncertainty bands.",
            },
        },
        "result_plot": None,
    },

    "Joint Inversion": {
        "type": "llm",
        "category": "AI Inversion",
        "class_name": "JointInversionAgent",
        "description": "Multi-modal MT joint inversion (DRCNN): MT + TEM / CSAMT / gravity.",
        "params": {
            "modalities": {
                "type": "str",
                "default": "mt,tem",
                "label": "Modalities (comma-sep)",
                "tip": "e.g. 'mt,tem' or 'mt,csamt,gravity'.",
            },
            "n_layers": {
                "type": "int",
                "default": 5,
                "range": (2, 30),
                "step": 1,
                "label": "Output layers",
            },
            "epochs": {
                "type": "int",
                "default": 30,
                "range": (1, 500),
                "step": 5,
                "label": "Training epochs",
            },
        },
        "result_plot": None,
    },

    "Model Zoo": {
        "type": "llm",
        "category": "AI Inversion",
        "class_name": "ModelZooAgent",
        "description": "Browse, download, and deploy pre-trained EM inverter checkpoints.",
        "params": {
            "force_download": {
                "type": "bool",
                "default": False,
                "label": "Force re-download",
                "tip": "Re-download even if the checkpoint is already cached.",
            },
        },
        "result_plot": None,
    },

    # ── Post-processing & Evaluation ──────────────────────────────────────────

    "Inversion Evaluation": {
        "type": "llm",
        "category": "Post-processing",
        "class_name": "InversionEvaluationAgent",
        "description": "Compute RMS, residual phase-tensor section, and misfit pseudosection.",
        "params": {},
        "result_plot": None,
    },

    "Inversion Comparison": {
        "type": "llm",
        "category": "Post-processing",
        "class_name": "InversionComparisonAgent",
        "description": "Compare two resistivity inversion results side-by-side.",
        "params": {},
        "result_plot": None,
    },

    "Resistivity Map": {
        "type": "llm",
        "category": "Post-processing",
        "class_name": "ResistivityMapAgent",
        "description": "Build horizontal resistivity depth-slice maps from 1-D inversion results.",
        "params": {
            "interp_method": {
                "type": "combo",
                "options": ["linear", "nearest", "cubic"],
                "default": "linear",
                "label": "Interpolation method",
            },
            "grid_n": {
                "type": "int",
                "default": 50,
                "range": (10, 500),
                "step": 10,
                "label": "Grid resolution (n×n)",
            },
            "depth_indices": {
                "type": "str",
                "default": "",
                "label": "Depth slice indices",
                "tip": "Comma-separated 0-based depth indices. Leave blank for all.",
            },
        },
        "result_plot": None,
    },

    # ── Interpretation & Reporting ─────────────────────────────────────────────

    "Interpretation": {
        "type": "llm",
        "category": "Interpretation & Reporting",
        "class_name": "InterpretationAgent",
        "description": "Map resistivity ranges to lithology; correlate with borehole logs.",
        "params": {
            "context": {
                "type": "str",
                "default": "",
                "label": "Geological context",
                "tip": "Free-text hint for the LLM (e.g. 'sedimentary basin, Cretaceous').",
            },
        },
        "result_plot": None,
    },

    "Report": {
        "type": "llm",
        "category": "Interpretation & Reporting",
        "class_name": "ReportAgent",
        "description": "Assemble all figures and tables into Markdown / HTML / PDF.",
        "params": {
            "report_title": {
                "type": "str",
                "default": "MT/AMT Survey Report",
                "label": "Report title",
            },
            "formats": {
                "type": "str",
                "default": "markdown,html",
                "label": "Output formats",
                "tip": "Comma-separated: markdown, html, pdf.",
            },
        },
        "result_plot": None,
    },

    "Code Generation": {
        "type": "llm",
        "category": "Interpretation & Reporting",
        "class_name": "CodeGenerationAgent",
        "description": "Emit a standalone reproducible Python script from the workflow config.",
        "params": {
            "script_title": {
                "type": "str",
                "default": "pycsamt MT Processing Workflow",
                "label": "Script title",
            },
        },
        "result_plot": None,
    },

    "EDI Export": {
        "type": "llm",
        "category": "Interpretation & Reporting",
        "class_name": "EDIExportAgent",
        "description": "Write processed Sites to EDI files on disk.",
        "params": {
            "file_pattern": {
                "type": "str",
                "default": "{station}.edi",
                "label": "Filename pattern",
                "tip": "Use {station} as placeholder for the station name.",
            },
            "overwrite": {
                "type": "bool",
                "default": False,
                "label": "Overwrite existing",
            },
        },
        "result_plot": None,
    },

    # ── Workflow & Orchestration ───────────────────────────────────────────────

    "Workflow Orchestrator": {
        "type": "llm",
        "category": "Workflow",
        "class_name": "WorkflowOrchestratorAgent",
        "description": "NL → classify workflow → build and run the correct agent chain.",
        "params": {
            "default_workflow": {
                "type": "combo",
                "options": [
                    "qc",
                    "phase_analysis",
                    "pre_inversion",
                    "ai_inversion",
                    "inv2d",
                    "inv3d",
                    "ensemble_inversion",
                    "joint_inversion",
                    "modem",
                    "tipper",
                    "sensitivity",
                    "rotation",
                    "full",
                ],
                "default": "qc",
                "label": "Default workflow",
                "tip": "Fallback workflow when the NL request is ambiguous.",
            },
        },
        "result_plot": None,
    },

    "Pipeline": {
        "type": "llm",
        "category": "Workflow",
        "class_name": "PipelineAgent",
        "description": "LLM-assisted MT processing pipeline: select preset, run, and interpret.",
        "params": {
            "preset": {
                "type": "combo",
                "options": ["basic_qc", "full_qc", "inversion_prep", "ai_inversion", "full"],
                "default": "basic_qc",
                "label": "Pipeline preset",
                "tip": "Named preset defining the sequence of processing steps.",
            },
        },
        "result_plot": None,
    },

    "Batch Survey": {
        "type": "llm",
        "category": "Workflow",
        "class_name": "BatchSurveyAgent",
        "description": "Process multiple MT profiles through a shared agent chain.",
        "params": {
            "workflow": {
                "type": "combo",
                "options": [
                    "qc",
                    "phase_analysis",
                    "pre_inversion",
                    "ai_inversion",
                    "full",
                ],
                "default": "qc",
                "label": "Workflow",
            },
            "n_jobs": {
                "type": "int",
                "default": 1,
                "range": (1, 32),
                "step": 1,
                "label": "Parallel jobs",
                "tip": "Number of parallel profiles to process (requires joblib).",
            },
        },
        "result_plot": None,
    },

    # ══════════════════════════════════════════════════════════════════════════
    # Processing operations (no API key required)
    # ══════════════════════════════════════════════════════════════════════════

    "QC Quicklook": {
        "type": "processing",
        "fn_name": "plot_qc_quicklook",
        "description": "Fast QC overview plot: SNR, skew, and coverage across all stations.",
        "params": {},
        "result_plot": "plot_qc_quicklook",
    },

    "Dimensionality": {
        "type": "processing",
        "fn_name": "classify_dimensionality",
        "description": "Classify 1-D / 2-D / 3-D structure from Swift skew and phase tensor.",
        "params": {
            "skew_th": {
                "type": "float",
                "default": 3.0,
                "range": (0.5, 20.0),
                "step": 0.5,
                "label": "Skew threshold (°)",
            },
            "ellipt_th": {
                "type": "float",
                "default": 0.2,
                "range": (0.01, 1.0),
                "step": 0.01,
                "label": "Ellipticity threshold",
            },
        },
        "result_plot": "plot_dimensionality_psection",
    },

    "Static Shift (fast)": {
        "type": "processing",
        "fn_name": "correct_static_shift",
        "description": "Correct static shift using spatial averaging (no LLM review).",
        "params": {
            "window_m": {
                "type": "float",
                "default": 1500.0,
                "range": (100.0, 10000.0),
                "step": 100.0,
                "label": "Window radius (m)",
            },
            "spacing_m": {
                "type": "float",
                "default": 200.0,
                "range": (10.0, 2000.0),
                "step": 10.0,
                "label": "Station spacing (m)",
            },
            "comp": {
                "type": "combo",
                "options": ["det", "xy", "yx"],
                "default": "det",
                "label": "Impedance component",
            },
            "inplace": {
                "type": "bool",
                "default": False,
                "label": "Correct in-place",
            },
        },
        "result_plot": "plot_rho_phase_bode",
    },
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper functions used by RunAgentDialog and AgentWorker
# ──────────────────────────────────────────────────────────────────────────────

def agent_names() -> list[str]:
    """Return all registered agent display names."""
    return list(AGENT_REGISTRY.keys())


def llm_agents() -> list[str]:
    return [n for n, e in AGENT_REGISTRY.items() if e["type"] == "llm"]


def processing_agents() -> list[str]:
    return [n for n, e in AGENT_REGISTRY.items() if e["type"] == "processing"]


def agents_by_category() -> dict[str, list[str]]:
    """Return {category: [display_names]} for LLM agents (preserves insertion order)."""
    cats: dict[str, list[str]] = {}
    for name, entry in AGENT_REGISTRY.items():
        if entry["type"] != "llm":
            continue
        cat = entry.get("category", "Other")
        cats.setdefault(cat, []).append(name)
    return cats


def get_entry(name: str) -> AgentEntry | None:
    return AGENT_REGISTRY.get(name)


def default_params(name: str) -> dict[str, Any]:
    """Return a {param_name: default_value} dict for the named agent."""
    entry = AGENT_REGISTRY.get(name, {})
    return {
        k: v.get("default")
        for k, v in entry.get("params", {}).items()
    }
