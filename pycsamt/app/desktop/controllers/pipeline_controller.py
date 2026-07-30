# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PipelineController — step definitions, state machine and execution logic
for the pycsamt processing pipeline.

Data flow
─────────
Each step receives ``input_sites`` (a Sites object) and returns
``output_sites`` — which becomes the next step's input.  The controller
stores a snapshot after every completed step so the user can roll back
or review intermediate results.

Step catalogue (8 steps)
────────────────────────
  0  Load          — use current / load from folder
  1  QC Screening  — flag / drop low-confidence frequencies
  2  Freq Edit     — edit / drop bad frequencies
  3  Static Shift  — AMA / LOESS / bilateral / ref-median
  4  Noise Removal — remove_noise_pipeline / EMAP filter / RPCA
  5  Strike        — consensus / phase-tensor / sweep estimate
  6  Rotation      — rotate to strike / Z rotation / skip
  7  Export        — save processed EDI files
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable

# ── Step status ───────────────────────────────────────────────────────────────


class StepStatus(Enum):
    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    DONE = "done"
    ERROR = "error"
    SKIPPED = "skipped"


STATUS_ICON = {
    StepStatus.PENDING: "⬜",
    StepStatus.QUEUED: "⏳",
    StepStatus.RUNNING: "▶",
    StepStatus.DONE: "✅",
    StepStatus.ERROR: "❌",
    StepStatus.SKIPPED: "⏭",
}


# ── Method / parameter specs ──────────────────────────────────────────────────


@dataclass
class ParamSpec:
    """Describes one configurable parameter."""

    kind: str  # "combo" | "float" | "int" | "bool" | "path"
    label: str
    default: Any
    options: list = field(default_factory=list)  # for kind == "combo"
    min_val: float = 0.0
    max_val: float = 1e9
    decimals: int = 3


@dataclass
class MethodSpec:
    name: str
    label: str
    params: dict[str, ParamSpec] = field(default_factory=dict)


# ── Pipeline step definition ──────────────────────────────────────────────────


@dataclass
class PipelineStep:
    id: int
    name: str
    description: str
    methods: list[MethodSpec]
    default_method: str  # MethodSpec.name
    can_skip: bool = True
    abort_on_error: bool = False

    # Runtime state
    status: StepStatus = StepStatus.PENDING
    active_method: str = ""  # filled from default_method at init
    params: dict[str, Any] = field(default_factory=dict)
    output_sites: Any = None
    result_info: str = ""
    error_msg: str = ""
    elapsed_s: float = 0.0
    diag_fn: str | None = None  # emtools fn name for the preview plot

    def __post_init__(self):
        if not self.active_method:
            self.active_method = self.default_method
        # Seed params from the default method
        self.reset_params()

    def reset_params(self) -> None:
        m = self.get_method(self.active_method)
        if m:
            self.params = {k: p.default for k, p in m.params.items()}

    def get_method(self, name: str) -> MethodSpec | None:
        for m in self.methods:
            if m.name == name:
                return m
        return None

    @property
    def active_method_spec(self) -> MethodSpec | None:
        return self.get_method(self.active_method)


# ── Step catalogue ────────────────────────────────────────────────────────────


def _build_steps() -> list[PipelineStep]:
    return [
        # ── 0  Load ───────────────────────────────────────────────────────────
        PipelineStep(
            id=0,
            name="Load Data",
            description=(
                "Load survey EDI files as the pipeline input.\n"
                "Choose whether to use the data already loaded in the main "
                "window, or browse to a fresh EDI folder."
            ),
            can_skip=False,
            methods=[
                MethodSpec("current", "Use current data"),
                MethodSpec(
                    "folder",
                    "Load from folder",
                    params={
                        "folder": ParamSpec("path", "EDI folder", ""),
                    },
                ),
            ],
            default_method="current",
            diag_fn=None,
        ),
        # ── 1  QC Screening ───────────────────────────────────────────────────
        PipelineStep(
            id=1,
            name="QC Screening",
            description=(
                "Assess data quality and flag / remove frequencies with "
                "low confidence.  Results are visible in the Preview tab."
            ),
            methods=[
                MethodSpec(
                    "drop_low_conf",
                    "Drop low-confidence frequencies",
                    params={
                        "threshold": ParamSpec(
                            "float",
                            "Min confidence",
                            0.5,
                            min_val=0.0,
                            max_val=1.0,
                            decimals=2,
                        ),
                        "method": ParamSpec(
                            "combo",
                            "Score method",
                            "coverage",
                            options=["coverage", "snr", "coherence"],
                        ),
                    },
                ),
                MethodSpec(
                    "flag_only",
                    "Flag only (no drop)",
                    params={
                        "threshold": ParamSpec(
                            "float",
                            "Min confidence",
                            0.5,
                            min_val=0.0,
                            max_val=1.0,
                            decimals=2,
                        ),
                    },
                ),
            ],
            default_method="drop_low_conf",
            diag_fn="plot_confidence_band_summary",
        ),
        # ── 2  Frequency Edit ─────────────────────────────────────────────────
        PipelineStep(
            id=2,
            name="Frequency Edit",
            description=(
                "Edit frequency set: remove noisy bands, close skew gaps, "
                "or trim to the period range of interest."
            ),
            methods=[
                MethodSpec(
                    "edit_by_conf",
                    "Edit by confidence",
                    params={
                        "mode": ParamSpec(
                            "combo",
                            "Mode",
                            "recover",
                            options=["recover", "drop", "mask"],
                        ),
                        "ci_hi": ParamSpec(
                            "float",
                            "High CI threshold",
                            0.8,
                            min_val=0.0,
                            max_val=1.0,
                            decimals=2,
                        ),
                    },
                ),
                MethodSpec(
                    "close_gaps",
                    "Close skew gaps",
                    params={
                        "max_gap": ParamSpec(
                            "int",
                            "Max gap (decades×10)",
                            3,
                            min_val=1,
                            max_val=20,
                        ),
                    },
                ),
                MethodSpec(
                    "regrid",
                    "Regrid to log-space",
                    params={
                        "n_freq": ParamSpec(
                            "int",
                            "Target N frequencies",
                            30,
                            min_val=5,
                            max_val=200,
                        ),
                    },
                ),
            ],
            default_method="edit_by_conf",
            diag_fn="plot_frequency_confidence_psection",
        ),
        # ── 3  Static Shift ───────────────────────────────────────────────────
        PipelineStep(
            id=3,
            name="Static Shift Correction",
            description=(
                "Correct near-surface static shift distortions using "
                "spatial averaging (AMA), LOESS, bilateral or reference-"
                "median methods."
            ),
            methods=[
                MethodSpec(
                    "ama",
                    "Azimuthal Moving Average (AMA)",
                    params={
                        "half_window": ParamSpec(
                            "int",
                            "Half window (stations)",
                            3,
                            min_val=1,
                            max_val=20,
                        ),
                        "weights": ParamSpec(
                            "combo",
                            "Weights",
                            "tri",
                            options=["tri", "gauss", "uniform"],
                        ),
                        "max_skew": ParamSpec(
                            "float",
                            "Max skew (°)",
                            6.0,
                            min_val=0.0,
                            max_val=20.0,
                            decimals=1,
                        ),
                    },
                ),
                MethodSpec(
                    "loess",
                    "LOESS spatial smoothing",
                    params={
                        "half_window": ParamSpec(
                            "int",
                            "Half window (stations)",
                            3,
                            min_val=1,
                            max_val=20,
                        ),
                    },
                ),
                MethodSpec(
                    "bilateral",
                    "Bilateral filter",
                    params={
                        "sig_dist": ParamSpec(
                            "float",
                            "Distance σ (stations)",
                            2.0,
                            min_val=0.5,
                            max_val=20.0,
                            decimals=1,
                        ),
                        "sig_val": ParamSpec(
                            "float",
                            "Value σ (log-ρ)",
                            0.5,
                            min_val=0.05,
                            max_val=5.0,
                            decimals=2,
                        ),
                    },
                ),
                MethodSpec(
                    "refmedian",
                    "Reference median",
                    params={
                        "pband": ParamSpec(
                            "combo",
                            "Period band",
                            "none",
                            options=["none", "short", "long", "all"],
                        ),
                    },
                ),
            ],
            default_method="ama",
            diag_fn="ss_qc_profile",
        ),
        # ── 4  Noise Removal ──────────────────────────────────────────────────
        PipelineStep(
            id=4,
            name="Noise Removal",
            description=(
                "Remove harmonic power-line noise, smooth via moving "
                "average, or apply robust PCA / EMAP filter."
            ),
            methods=[
                MethodSpec(
                    "pipeline",
                    "Auto noise pipeline",
                    params={
                        "mains_hz": ParamSpec(
                            "float",
                            "Mains frequency (Hz)",
                            50.0,
                            min_val=0.0,
                            max_val=400.0,
                            decimals=1,
                        ),
                        "smooth_win": ParamSpec(
                            "int", "Smooth window", 3, min_val=1, max_val=15
                        ),
                    },
                ),
                MethodSpec(
                    "emap",
                    "EMAP spatial filter",
                    params={
                        "method": ParamSpec(
                            "combo",
                            "Filter method",
                            "median",
                            options=["median", "mean", "gaussian"],
                        ),
                        "window": ParamSpec(
                            "int",
                            "Window (stations)",
                            5,
                            min_val=3,
                            max_val=21,
                        ),
                    },
                ),
                MethodSpec(
                    "rpca",
                    "Robust PCA denoising",
                    params={
                        "keep_phase": ParamSpec(
                            "bool", "Preserve phase", True
                        ),
                    },
                ),
            ],
            default_method="pipeline",
            diag_fn="plot_snr_hist",
        ),
        # ── 5  Strike Analysis ────────────────────────────────────────────────
        PipelineStep(
            id=5,
            name="Strike Analysis",
            description=(
                "Estimate the regional geoelectric strike direction.  "
                "The result is stored and used by the rotation step."
            ),
            methods=[
                MethodSpec(
                    "consensus",
                    "Consensus (sweep + PT)",
                    params={
                        "w_sweep": ParamSpec(
                            "float",
                            "Sweep weight",
                            0.5,
                            min_val=0.0,
                            max_val=1.0,
                            decimals=2,
                        ),
                        "w_pt": ParamSpec(
                            "float",
                            "Phase-tensor weight",
                            0.5,
                            min_val=0.0,
                            max_val=1.0,
                            decimals=2,
                        ),
                    },
                ),
                MethodSpec("phase_tensor", "Phase-tensor only", params={}),
                MethodSpec(
                    "sweep",
                    "Swift/Bahr sweep",
                    params={
                        "metric": ParamSpec(
                            "combo",
                            "Metric",
                            "diag_ratio",
                            options=["diag_ratio", "off_diag", "det"],
                        ),
                    },
                ),
            ],
            default_method="consensus",
            diag_fn="plot_strike_rose",
        ),
        # ── 6  Rotation ───────────────────────────────────────────────────────
        PipelineStep(
            id=6,
            name="Strike Rotation",
            description=(
                "Rotate impedance tensors to the estimated strike direction "
                "from step 5, or supply a manual angle."
            ),
            methods=[
                MethodSpec(
                    "to_strike",
                    "Auto rotate to strike",
                    params={
                        "method": ParamSpec(
                            "combo",
                            "Strike method",
                            "consensus",
                            options=["consensus", "phase_tensor", "sweep"],
                        ),
                    },
                ),
                MethodSpec(
                    "z_rotation",
                    "Z tensor rotation (Swift)",
                    params={
                        "method": ParamSpec(
                            "combo",
                            "Method",
                            "swift",
                            options=["swift", "bahr"],
                        ),
                    },
                ),
                MethodSpec("none", "No rotation (skip)", params={}),
            ],
            default_method="to_strike",
            diag_fn="plot_phase_tensor_skewmap",
        ),
        # ── 7  Export ─────────────────────────────────────────────────────────
        PipelineStep(
            id=7,
            name="Export",
            description=(
                "Save the processed Sites collection as EDI files "
                "to a chosen output folder."
            ),
            can_skip=False,
            methods=[
                MethodSpec(
                    "edi",
                    "Export EDI files",
                    params={
                        "folder": ParamSpec("path", "Output folder", ""),
                    },
                ),
            ],
            default_method="edi",
            diag_fn=None,
        ),
    ]


# ── Controller ────────────────────────────────────────────────────────────────


class PipelineController:
    """
    Manages the 8-step pipeline state, data snapshots and step execution.
    """

    def __init__(self) -> None:
        self.steps: list[PipelineStep] = _build_steps()
        self._sites_input = None  # Sites passed in from main window
        self._sites_chain: list = [None] * 8  # snapshot after each step

    # ── Data binding ──────────────────────────────────────────────────────────

    def set_input_sites(self, sites) -> None:
        """Called by MainWindow when EDI data is available."""
        self._sites_input = sites
        # Pre-fill step 0 with "Use current data" if already set
        step0 = self.steps[0]
        if step0.active_method == "current" and sites is not None:
            step0.result_info = (
                f"Ready — {getattr(sites, 'n_sites', '?')} stations available"
            )

    def is_ready(self) -> bool:
        return self._sites_input is not None

    # ── Step execution ────────────────────────────────────────────────────────

    def execute_step(
        self,
        step_id: int,
        log_cb: Callable[[str], None] | None = None,
    ) -> Any:
        """
        Execute one step and return the output Sites.

        Parameters
        ----------
        step_id : int
            Index into self.steps.
        log_cb : callable(str), optional
            Called with each log message during execution.

        Returns
        -------
        Sites
            Processed survey (or original if step modifies nothing).
        """
        step = self.steps[step_id]

        def log(msg: str) -> None:
            if log_cb:
                log_cb(msg)

        # Determine input
        if step_id == 0:
            sites_in = self._sites_input
        else:
            # Walk back to find the last successful output
            sites_in = None
            for prev in range(step_id - 1, -1, -1):
                if self._sites_chain[prev] is not None:
                    sites_in = self._sites_chain[prev]
                    break
            if sites_in is None:
                sites_in = self._sites_input

        t0 = time.perf_counter()
        step.status = StepStatus.RUNNING
        step.error_msg = ""
        step.result_info = ""

        try:
            result = self._run_step(step, sites_in, log)
        except Exception as exc:
            step.status = StepStatus.ERROR
            step.error_msg = str(exc)
            step.elapsed_s = time.perf_counter() - t0
            log(f"ERROR: {exc}")
            raise

        step.elapsed_s = time.perf_counter() - t0
        step.status = StepStatus.DONE
        step.output_sites = result
        self._sites_chain[step_id] = result
        return result

    def _run_step(self, step: PipelineStep, sites_in, log) -> Any:
        import pycsamt.emtools as et

        sid = step.id
        meth = step.active_method
        prm = step.params

        # ── 0: Load ───────────────────────────────────────────────────────────
        if sid == 0:
            if meth == "current":
                log("Using currently loaded survey data.")
                step.result_info = f"{self._n(sites_in)} stations"
                return sites_in
            else:
                folder = prm.get("folder", "")
                if not folder:
                    raise ValueError("No folder selected for Load step.")
                log(f"Loading EDI files from: {folder}")
                from pycsamt.app.desktop.controllers.data_controller import (
                    DataController,
                )

                dc = DataController()
                import glob

                paths = glob.glob(str(Path(folder) / "*.edi"), recursive=False)
                if not paths:
                    raise ValueError(f"No EDI files found in {folder}")
                result = dc.load(paths)
                self._sites_input = result
                step.result_info = (
                    f"{len(paths)} files → {self._n(result)} stations"
                )
                log(f"Loaded {len(paths)} EDI files.")
                return result

        # ── 1: QC Screening ───────────────────────────────────────────────────
        if sid == 1:
            if meth == "flag_only":
                log("Flagging low-confidence frequencies (no drop).")
                step.result_info = "Frequencies flagged — data unchanged."
                return sites_in
            else:
                thresh = float(prm.get("threshold", 0.5))
                log(f"Dropping low-confidence frequencies (thresh={thresh}).")
                result = et.drop_low_confidence_frequencies(
                    sites_in, threshold=thresh, inplace=False
                )
                result = self._coerce_sites(result, sites_in)
                step.result_info = f"{self._n(result)} stations after QC"
                log(f"Done — {self._n(result)} stations retained.")
                return result

        # ── 2: Frequency Edit ─────────────────────────────────────────────────
        if sid == 2:
            if meth == "close_gaps":
                log("Closing skew gaps in frequency set.")
                result = et.close_skew_gaps(sites_in)
                result = self._coerce_sites(result, sites_in)
                step.result_info = "Gaps closed."
                return result
            elif meth == "regrid":
                n = int(prm.get("n_freq", 30))
                log(f"Regridding to {n} log-spaced frequencies.")
                result = et.regrid_logspace(sites_in, n_freq=n)
                result = self._coerce_sites(result, sites_in)
                step.result_info = f"Regridded to {n} frequencies."
                return result
            else:
                mode = prm.get("mode", "auto")
                ci_hi = float(prm.get("ci_hi", 0.8))
                log(
                    f"Editing frequencies by confidence (mode={mode}, ci_hi={ci_hi})."
                )
                result = et.edit_frequencies_by_confidence(
                    sites_in, mode=mode, ci_hi=ci_hi
                )
                result = self._coerce_sites(result, sites_in)
                step.result_info = f"Frequencies edited ({mode})."
                return result

        # ── 3: Static Shift ───────────────────────────────────────────────────
        if sid == 3:
            if meth == "ama":
                hw = int(prm.get("half_window", 3))
                wt = prm.get("weights", "tri")
                sk = float(prm.get("max_skew", 6.0))
                log(
                    f"AMA static-shift correction (window={hw}, weights={wt})."
                )
                result = et.correct_ss_ama(
                    sites_in,
                    half_window=hw,
                    weights=wt,
                    max_skew=sk,
                    inplace=False,
                )
            elif meth == "loess":
                hw = int(prm.get("half_window", 3))
                log(f"LOESS static-shift correction (half_window={hw}).")
                tbl = et.estimate_ss_loess(sites_in, half_window=hw)
                result = et.apply_ss_factors(sites_in, tbl, inplace=False)
            elif meth == "bilateral":
                sd = float(prm.get("sig_dist", 2.0))
                sv = float(prm.get("sig_val", 0.5))
                log(
                    f"Bilateral static-shift correction (sig_dist={sd}, sig_val={sv})."
                )
                tbl = et.estimate_ss_bilateral(
                    sites_in, sig_dist=sd, sig_val=sv
                )
                result = et.apply_ss_factors(sites_in, tbl, inplace=False)
            else:
                log("Reference-median static-shift correction.")
                tbl = et.estimate_ss_refmedian(sites_in)
                result = et.apply_ss_factors(sites_in, tbl, inplace=False)
            result = self._coerce_sites(result, sites_in)
            step.result_info = f"SS corrected — {self._n(result)} stations"
            log("Static shift correction complete.")
            return result

        # ── 4: Noise Removal ──────────────────────────────────────────────────
        if sid == 4:
            if meth == "emap":
                filt = prm.get("method", "median")
                win = int(prm.get("window", 5))
                log(f"EMAP filter ({filt}, window={win}).")
                result = et.apply_emap_filter(
                    sites_in, method=filt, window=win
                )
            elif meth == "rpca":
                kp = bool(prm.get("keep_phase", True))
                log(f"RPCA denoising (keep_phase={kp}).")
                try:
                    result = et.rpca_offdiag_denoise(
                        sites_in, keep_phase=kp, inplace=False
                    )
                except Exception:
                    # fallback to pipeline if RPCA fails
                    log("RPCA failed — falling back to noise pipeline.")
                    result = et.remove_noise_pipeline(sites_in)
            else:
                mhz = float(prm.get("mains_hz", 50.0))
                sw = int(prm.get("smooth_win", 3))
                log(f"Noise pipeline (mains={mhz} Hz, smooth_win={sw}).")
                result = et.remove_noise_pipeline(
                    sites_in, mains_hz=mhz, smooth_win=sw
                )
            result = self._coerce_sites(result, sites_in)
            step.result_info = "Noise removed."
            log("Noise removal complete.")
            return result

        # ── 5: Strike Analysis ────────────────────────────────────────────────
        if sid == 5:
            if meth == "phase_tensor":
                log("Estimating strike via phase-tensor method.")
                tbl = et.estimate_strike_phase_tensor(sites_in)
            elif meth == "sweep":
                metric = prm.get("metric", "diag_ratio")
                log(f"Estimating strike via sweep (metric={metric}).")
                tbl = et.estimate_strike_sweep(sites_in, metric=metric)
            else:
                ws = float(prm.get("w_sweep", 0.5))
                wp = float(prm.get("w_pt", 0.5))
                log(f"Consensus strike estimate (w_sweep={ws}, w_pt={wp}).")
                tbl = et.estimate_strike_consensus(
                    sites_in, w_sweep=ws, w_pt=wp
                )
            import pandas as pd

            if isinstance(tbl, pd.DataFrame) and not tbl.empty:
                angle_col = [
                    c
                    for c in ("strike", "angle", "theta", "azimuth")
                    if c in tbl.columns
                ]
                if angle_col:
                    import numpy as np

                    angle = float(np.nanmedian(tbl[angle_col[0]].values))
                    step.result_info = f"Strike ≈ {angle:.1f}°"
                    log(f"Estimated strike: {angle:.1f}°")
                else:
                    step.result_info = "Strike table computed."
            else:
                step.result_info = "Strike analysis complete."
            return sites_in  # strike step doesn't transform sites

        # ── 6: Rotation ───────────────────────────────────────────────────────
        if sid == 6:
            if meth == "none":
                log("Rotation skipped by user.")
                step.result_info = "No rotation applied."
                return sites_in
            elif meth == "z_rotation":
                swift_meth = prm.get("method", "swift")
                log(f"Z-tensor rotation (method={swift_meth}).")
                result = et.rotate_z_to_strike(
                    sites_in, method=swift_meth, inplace=False
                )
            else:
                rot_meth = prm.get("method", "consensus")
                log(f"Rotating to strike (method={rot_meth}).")
                result = et.rotate_to_strike(
                    sites_in, method=rot_meth, inplace=False
                )
            result = self._coerce_sites(result, sites_in)
            step.result_info = f"Rotated — {self._n(result)} stations"
            log("Rotation complete.")
            return result

        # ── 7: Export ─────────────────────────────────────────────────────────
        if sid == 7:
            folder = prm.get("folder", "")
            if not folder:
                raise ValueError("No output folder selected for Export step.")
            Path(folder).mkdir(parents=True, exist_ok=True)
            log(f"Exporting EDI files to: {folder}")
            edis = [s.edi for s in sites_in if hasattr(s, "edi")]
            if not edis:
                edis = list(sites_in.as_list())
            et.export_edis(edis, savepath=folder)
            n = len(edis)
            step.result_info = f"Exported {n} EDI files to {folder}"
            log(f"Export complete — {n} files written.")
            return sites_in

        raise ValueError(f"No handler for step {sid}")

    # ── Reset / config ────────────────────────────────────────────────────────

    def reset(self) -> None:
        """Reset all steps to PENDING, clear intermediate results."""
        for step in self.steps:
            step.status = StepStatus.PENDING
            step.output_sites = None
            step.result_info = ""
            step.error_msg = ""
            step.elapsed_s = 0.0
        self._sites_chain = [None] * 8

    def to_config(self) -> dict:
        """Serialise current method + param choices to a dict."""
        return {
            "steps": [
                {
                    "id": s.id,
                    "name": s.name,
                    "method": s.active_method,
                    "params": dict(s.params),
                }
                for s in self.steps
            ]
        }

    def from_config(self, cfg: dict) -> None:
        """Load method + param choices from a previously saved dict."""
        for entry in cfg.get("steps", []):
            sid = entry.get("id")
            if sid is None or sid >= len(self.steps):
                continue
            step = self.steps[sid]
            step.active_method = entry.get("method", step.default_method)
            for k, v in entry.get("params", {}).items():
                if k in step.params:
                    step.params[k] = v

    def get_output_sites(self) -> Any:
        """Return the final Sites after all completed steps."""
        for i in range(len(self._sites_chain) - 1, -1, -1):
            if self._sites_chain[i] is not None:
                return self._sites_chain[i]
        return self._sites_input

    # ── Helpers ───────────────────────────────────────────────────────────────

    @staticmethod
    def _n(sites) -> int:
        try:
            return len(list(sites))
        except Exception:
            return 0

    @staticmethod
    def _coerce_sites(result, fallback):
        """Return result if it looks like Sites, otherwise fallback."""
        if result is None:
            return fallback
        if hasattr(result, "as_list"):
            return result
        # APIResult wraps Sites in some functions
        inner = getattr(result, "sites", None) or getattr(result, "data", None)
        if inner is not None and hasattr(inner, "as_list"):
            return inner
        return fallback
