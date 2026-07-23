# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
CorrectionController — non-destructive correction stack for AMT/CSAMT data.

Architecture
────────────
Raw Sites (immutable) → step₁ → step₂ → … → current corrected Sites

Each ``apply()`` call computes a new Sites and pushes a CorrectionStep to the
stack.  Removing any step triggers a replay from that point so the chain
remains consistent.  The raw Sites is never modified.

``preview()`` computes a result without touching the stack — useful for
the "Before/After" canvas update before the user commits.

``commit()`` is handled by the window (emitting a Signal to MainWindow).
"""

from __future__ import annotations

from collections import namedtuple
from dataclasses import dataclass, field
from typing import Any

import numpy as np

# ── Parameter specification ────────────────────────────────────────────────────
# kind: "spin" | "dspin" | "combo" | "check"
# opts: (min, max, step) for spin/dspin; [choices] for combo; None for check

ParamSpec = namedtuple(
    "ParamSpec",
    ["name", "label", "kind", "default", "opts", "tip"],
    defaults=[""],
)

# ── Correction catalogue ───────────────────────────────────────────────────────
# Structure: {category_name: {correction_label: {fn, desc, params}}}

CATALOGUE: dict[str, dict[str, dict]] = {
    "Static Shift": {
        "AMA (spatial average)": {
            "fn": "correct_ss_ama",
            "desc": "Array Moving-Average: spatially averages impedances across "
            "neighbouring stations to estimate and remove the static-shift "
            "multiplier on each Z column.",
            "params": [
                ParamSpec(
                    "sort_by",
                    "Sort by",
                    "combo",
                    "lon",
                    ["lon", "lat", "id"],
                    "Axis used to order stations before spatial averaging.",
                ),
                ParamSpec(
                    "half_window",
                    "Half window",
                    "spin",
                    3,
                    (1, 20, 1),
                    "Number of neighbours on each side (total = 2k+1).",
                ),
                ParamSpec(
                    "weights",
                    "Kernel",
                    "combo",
                    "tri",
                    ["tri", "box", "cos"],
                    "Spatial weight kernel.",
                ),
                ParamSpec(
                    "max_skew",
                    "Max skew (°)",
                    "dspin",
                    6.0,
                    (0.5, 45.0, 0.5),
                    "Stations with |skew| > this value are excluded.",
                ),
            ],
        },
        "LOESS": {
            "fn": "_correct_ss_loess",
            "desc": "LOESS trend estimation: fits a local polynomial regression "
            "to estimate the underlying smooth ρ_a trend and derives the "
            "static-shift factor per station.",
            "params": [
                ParamSpec(
                    "half_window", "Half window", "spin", 3, (1, 20, 1)
                ),
                ParamSpec("poly", "Poly degree", "spin", 1, (1, 3, 1)),
                ParamSpec(
                    "it",
                    "Iterations",
                    "spin",
                    2,
                    (1, 10, 1),
                    "Robust re-weighting iterations.",
                ),
                ParamSpec(
                    "max_skew", "Max skew (°)", "dspin", 6.0, (0.5, 45.0, 0.5)
                ),
            ],
        },
        "Bilateral filter": {
            "fn": "_correct_ss_bilateral",
            "desc": "Bilateral spatial filter combining distance and value "
            "kernels so that sharp lateral ρ_a contrasts (geology) are "
            "preserved while smoothing the static contribution.",
            "params": [
                ParamSpec(
                    "half_window", "Half window", "spin", 4, (1, 20, 1)
                ),
                ParamSpec(
                    "max_skew", "Max skew (°)", "dspin", 6.0, (0.5, 45.0, 0.5)
                ),
            ],
        },
        "Reference median": {
            "fn": "_correct_ss_refmedian",
            "desc": "Uses a global reference (median over all stations) to "
            "derive per-station static shift. Robust for profiles without "
            "clear spatial trend.",
            "params": [
                ParamSpec(
                    "smooth_sites",
                    "Smooth sites",
                    "spin",
                    0,
                    (0, 5, 1),
                    "Optional pre-smoothing iterations (0 = none).",
                ),
                ParamSpec(
                    "max_skew", "Max skew (°)", "dspin", 6.0, (0.5, 45.0, 0.5)
                ),
            ],
        },
        "Hanning EMAP (Torres-Verdín)": {
            "fn": "correct_static_shift",
            "desc": "Hanning adaptive moving-average spatial filter (Torres-Verdín "
            "& Bostick 1992). Removes static shift by low-pass filtering "
            "the log(ρ_a) spatial profile at each frequency.",
            "params": [
                ParamSpec(
                    "window_m",
                    "Window (m)",
                    "dspin",
                    1500.0,
                    (100.0, 10000.0, 100.0),
                ),
                ParamSpec(
                    "spacing_m",
                    "Spacing (m)",
                    "dspin",
                    200.0,
                    (10.0, 2000.0, 10.0),
                ),
                ParamSpec(
                    "comp", "Component", "combo", "det", ["det", "xy", "yx"]
                ),
            ],
        },
    },
    "Noise Removal": {
        "Notch filter (power line)": {
            "fn": "notch_powerline",
            "desc": "Removes power-line harmonics (50 or 60 Hz and overtones) "
            "from the impedance tensor via masking or linear interpolation "
            "over affected frequency bins.",
            "params": [
                ParamSpec(
                    "mains_hz", "Mains (Hz)", "dspin", 50.0, (40.0, 70.0, 1.0)
                ),
                ParamSpec(
                    "n_harm",
                    "Harmonics",
                    "spin",
                    30,
                    (1, 100, 1),
                    "Number of harmonics to suppress.",
                ),
                ParamSpec(
                    "tol_hz",
                    "Tolerance (Hz)",
                    "dspin",
                    0.08,
                    (0.01, 2.0, 0.01),
                ),
                ParamSpec(
                    "mode",
                    "Mode",
                    "combo",
                    "interp",
                    ["interp", "mask"],
                    "interp: fill gaps; mask: leave as NaN.",
                ),
                ParamSpec(
                    "also",
                    "Apply to",
                    "combo",
                    "both",
                    ["both", "z", "tipper"],
                ),
            ],
        },
        "Log-freq smoothing": {
            "fn": "smooth_logfreq",
            "desc": "Triangular or boxcar kernel smoothing in log-frequency "
            "space. Reduces high-frequency noise while preserving "
            "broad spectral shape.",
            "params": [
                ParamSpec(
                    "win",
                    "Window size",
                    "spin",
                    5,
                    (2, 21, 1),
                    "Smoothing window half-width in frequency bins.",
                ),
                ParamSpec("kind", "Kernel", "combo", "tri", ["tri", "box"]),
                ParamSpec(
                    "also",
                    "Apply to",
                    "combo",
                    "both",
                    ["both", "z", "tipper"],
                ),
            ],
        },
        "ρ/φ trend smoothing": {
            "fn": "smooth_rho_phase",
            "desc": "Polynomial trend smoothing of apparent resistivity and "
            "phase versus log-frequency. ρ is fitted in log space and "
            "phase is unwrapped before fitting; the complex impedance "
            "tensor is then rebuilt so ρ and φ remain consistent.",
            "params": [
                ParamSpec(
                    "components",
                    "Components",
                    "combo",
                    "offdiag",
                    ["offdiag", "all", "xy", "yx", "xx", "yy", "diagonal"],
                    "Tensor components to smooth. offdiag = xy and yx.",
                ),
                ParamSpec(
                    "degree",
                    "Poly degree",
                    "spin",
                    3,
                    (0, 8, 1),
                    "Polynomial degree in log10(frequency).",
                ),
                ParamSpec(
                    "smooth_rho",
                    "Smooth ρ",
                    "check",
                    True,
                    None,
                    "Smooth apparent-resistivity amplitude.",
                ),
                ParamSpec(
                    "smooth_phase",
                    "Smooth φ",
                    "check",
                    True,
                    None,
                    "Smooth unwrapped impedance phase.",
                ),
                ParamSpec(
                    "blend",
                    "Blend",
                    "dspin",
                    1.0,
                    (0.0, 1.0, 0.05),
                    "0 keeps original curves; 1 fully applies the trend.",
                ),
                ParamSpec(
                    "robust",
                    "Robust fit",
                    "check",
                    True,
                    None,
                    "Down-weight isolated spikes during polynomial fitting.",
                ),
            ],
        },
    },
    "Source Effects": {
        "Near-field correction": {
            "fn": "_correct_near_field",
            "desc": "Corrects CSAMT impedance for near-field contamination. "
            "Divides each Z element by the complex transfer-function ratio "
            "F(p) = 1 − 3/p² + 3/p³ (equatorial HED).  In the far field "
            "F → 1 (no correction applied).",
            "params": [
                ParamSpec(
                    "source_offset",
                    "Source offset (m)",
                    "dspin",
                    500.0,
                    (10.0, 50000.0, 10.0),
                    "Source-receiver offset in metres.",
                ),
            ],
        },
    },
    "Tensor Rotation": {
        "Rotate by fixed angle": {
            "fn": "_wrap_rotate",
            "desc": "Rotate all Z (and tipper) components by a fixed angle θ.\n"
            "Z'(θ) = R(θ) Z R(−θ) where R is the 2×2 rotation matrix.\n"
            "Positive angles rotate counter-clockwise (geophysics convention).",
            "params": [
                ParamSpec(
                    "angle",
                    "Angle (°)",
                    "dspin",
                    0.0,
                    (-180.0, 180.0, 0.5),
                    "Rotation angle in degrees. Positive = CCW.",
                ),
            ],
        },
        "Rotate to geoelectric strike": {
            "fn": "_wrap_rotate_to_strike",
            "desc": "Automatically estimates the geoelectric strike per station and "
            "rotates each Z tensor to align the principal TE/TM axes with "
            "x → strike.  Crucial before 2-D inversion.",
            "params": [
                ParamSpec(
                    "method",
                    "Strike method",
                    "combo",
                    "swift",
                    ["swift", "bahr", "phase_diff"],
                    "Swift (default): minimises |Zxx|+|Zyy|.\n"
                    "Bahr: Bahr (1988) decomposition.\n"
                    "phase_diff: phase difference criterion.",
                ),
            ],
        },
        "Rotate to phase-tensor strike": {
            "fn": "_wrap_rotate_pt_strike",
            "desc": "Estimates strike from the phase tensor (Caldwell et al. 2004) "
            "per station and band, then rotates.  More robust than Swift for "
            "3-D structures.",
            "params": [
                ParamSpec(
                    "band_lo",
                    "Band T-min (s)",
                    "dspin",
                    0.001,
                    (1e-6, 1e4, 0.001),
                    "Lower period bound of the estimation band.",
                ),
                ParamSpec(
                    "band_hi",
                    "Band T-max (s)",
                    "dspin",
                    1000.0,
                    (0.001, 1e6, 1.0),
                    "Upper period bound of the estimation band.",
                ),
                ParamSpec("robust", "Robust (median)", "check", True, None),
            ],
        },
        "Rotate to profile azimuth": {
            "fn": "_wrap_rotate_to_profile",
            "desc": "Rotates Z so that the x-axis aligns with the survey profile "
            "azimuth.  Negative of the profile bearing is applied.\n"
            "Use azimuth=-1 to auto-infer from station positions via PCA.",
            "params": [
                ParamSpec(
                    "azimuth",
                    "Profile azimuth (°, auto=-1)",
                    "dspin",
                    -1.0,
                    (-1.0, 360.0, 0.5),
                    "Survey-line azimuth from north.  -1 = infer via PCA.",
                ),
            ],
        },
        "Antisymmetrize (2-D prep)": {
            "fn": "antisymmetrize",
            "desc": "Enforces off-diagonal antisymmetry: Zxx = −Zyy = 0,  "
            "Zxy = −Zyx.  Removes the non-1D/2D part of the tensor; "
            "required for strictly 2-D inversions.",
            "params": [
                ParamSpec(
                    "how",
                    "Average method",
                    "combo",
                    "rms",
                    ["rms", "mean"],
                    "rms: root-mean-square antisymmetric part; mean: arithmetic.",
                ),
            ],
        },
    },
    "Coordinates": {
        "Profile projection": {
            "fn": "_coord_profile_projection",
            "desc": "Projects each station onto the best-fit survey line, "
            "removing cross-profile scatter.  Azimuth can be inferred "
            "automatically via PCA or specified manually.",
            "params": [
                ParamSpec(
                    "azimuth",
                    "Azimuth (°, auto=-1)",
                    "dspin",
                    -1.0,
                    (-1.0, 360.0, 0.5),
                    "Profile azimuth from north. Use -1 for automatic PCA inference.",
                ),
                ParamSpec(
                    "keep_elevation",
                    "Keep elevation",
                    "check",
                    True,
                    None,
                    "Preserve original elevation (uncheck to interpolate).",
                ),
            ],
        },
        "Spacing regularization": {
            "fn": "_coord_spacing_regularize",
            "desc": "Redistributes stations at uniform spacing along the profile "
            "while keeping the same profile extent and station count.",
            "params": [
                ParamSpec(
                    "spacing_m",
                    "Target spacing (m)",
                    "dspin",
                    200.0,
                    (10.0, 5000.0, 10.0),
                    "Desired station spacing in metres.",
                ),
                ParamSpec(
                    "azimuth",
                    "Azimuth (°, auto=-1)",
                    "dspin",
                    -1.0,
                    (-1.0, 360.0, 0.5),
                ),
                ParamSpec(
                    "preserve_extent",
                    "Preserve extent",
                    "check",
                    True,
                    None,
                    "Keep the same start/end chainages.",
                ),
            ],
        },
        "Outlier snap-to-line": {
            "fn": "_coord_outlier_snap",
            "desc": "Detects stations that deviate more than the threshold from "
            "the profile line and snaps them to their nearest on-line point. "
            "Stations within the threshold are left untouched.",
            "params": [
                ParamSpec(
                    "threshold_m",
                    "Threshold (m)",
                    "dspin",
                    50.0,
                    (1.0, 2000.0, 1.0),
                    "Maximum allowed perpendicular offset in metres.",
                ),
                ParamSpec(
                    "azimuth",
                    "Azimuth (°, auto=-1)",
                    "dspin",
                    -1.0,
                    (-1.0, 360.0, 0.5),
                ),
            ],
        },
        "Elevation smoothing": {
            "fn": "_coord_elevation_smooth",
            "desc": "Smooths noisy GPS elevation values along the profile. "
            "GPS vertical accuracy is typically 3–5× worse than horizontal; "
            "LOESS or running-mean smoothing removes high-frequency noise.",
            "params": [
                ParamSpec(
                    "method",
                    "Method",
                    "combo",
                    "loess",
                    ["loess", "mean"],
                    "loess: local polynomial fit; mean: symmetric running average.",
                ),
                ParamSpec(
                    "window",
                    "Window",
                    "spin",
                    5,
                    (2, 21, 1),
                    "Half-window (LOESS) or full window (mean) in station count.",
                ),
            ],
        },
        "Coordinate shift": {
            "fn": "_coord_shift",
            "desc": "Applies a uniform (Δlat, Δlon, Δelev) offset to all stations. "
            "Corrects systematic datum errors or known GPS receiver biases.",
            "params": [
                ParamSpec(
                    "delta_lat",
                    "Δ Latitude (°)",
                    "dspin",
                    0.0,
                    (-1.0, 1.0, 0.0001),
                ),
                ParamSpec(
                    "delta_lon",
                    "Δ Longitude (°)",
                    "dspin",
                    0.0,
                    (-1.0, 1.0, 0.0001),
                ),
                ParamSpec(
                    "delta_elev",
                    "Δ Elevation (m)",
                    "dspin",
                    0.0,
                    (-500.0, 500.0, 0.1),
                ),
            ],
        },
        "Interpolate missing": {
            "fn": "_coord_interpolate_missing",
            "desc": "Fills NaN or zero coordinates by linear interpolation from "
            "neighbouring stations.  Useful for stations with GPS dropout.",
            "params": [
                ParamSpec(
                    "fill_nan_only",
                    "Fill NaN only",
                    "check",
                    True,
                    None,
                    "Uncheck to also replace exact-zero coordinates.",
                ),
            ],
        },
    },
    # ── Stratagem — Zonge-format AMT pipeline ────────────────────────────────
    # These corrections operate on edi_objects_ (EDIFile list) via the
    # stratagem subpackage.  The controller holds a parallel edi_objects_
    # state alongside the Sites chain so both can coexist.
    "Stratagem": {
        "QC Report": {
            "fn": "_strat_qc",
            "desc": "Station-level quality-control report.  Computes fraction of "
            "valid impedance cells, median SNR and phase-tensor skew per "
            "station.  Flags stations below the acceptance thresholds.",
            "params": [
                ParamSpec(
                    "min_frac_ok",
                    "Min valid fraction",
                    "dspin",
                    0.6,
                    (0.1, 1.0, 0.05),
                    "Flag station if valid-impedance fraction < this value.",
                ),
                ParamSpec(
                    "min_snr_med",
                    "Min SNR (median)",
                    "dspin",
                    2.0,
                    (0.5, 20.0, 0.5),
                    "Flag station if median SNR across all frequencies < this.",
                ),
                ParamSpec(
                    "max_skew_med",
                    "Max skew (°)",
                    "dspin",
                    6.0,
                    (1.0, 45.0, 0.5),
                    "Flag station if median |β| > this value.",
                ),
                ParamSpec(
                    "include_skew",
                    "Include skew",
                    "check",
                    True,
                    None,
                    "Compute and report phase-tensor skew column.",
                ),
            ],
        },
        "Static Shift (AMA)": {
            "fn": "_strat_static_shift",
            "desc": "Stratagem-native Array Moving-Average static-shift correction. "
            "Uses StaticShiftCorrector from pycsamt.stratagem; applies "
            "the AMA estimator with triangular, Gaussian or uniform spatial "
            "weighting and optional skew-based station masking.",
            "params": [
                ParamSpec(
                    "sort_by",
                    "Sort by",
                    "combo",
                    "lon",
                    ["lon", "lat", "name"],
                    "Spatial ordering axis before AMA computation.",
                ),
                ParamSpec(
                    "half_window",
                    "Half window",
                    "spin",
                    3,
                    (1, 20, 1),
                    "Neighbours on each side.  Total window = 2k + 1.",
                ),
                ParamSpec(
                    "weights",
                    "Kernel",
                    "combo",
                    "tri",
                    ["tri", "gauss", "uniform"],
                    "Spatial weight kernel applied to AMA.",
                ),
                ParamSpec(
                    "pband_lo",
                    "Period min (s)",
                    "dspin",
                    0.0,
                    (0.0, 1000.0, 0.001),
                    "Lower period bound for SS estimation band. 0 = no limit.",
                ),
                ParamSpec(
                    "pband_hi",
                    "Period max (s)",
                    "dspin",
                    0.0,
                    (0.0, 10000.0, 1.0),
                    "Upper period bound. 0 = no limit.",
                ),
                ParamSpec(
                    "max_skew",
                    "Max skew (°)",
                    "dspin",
                    6.0,
                    (0.5, 45.0, 0.5),
                    "Exclude stations with |β| > this from AMA calculation.",
                ),
            ],
        },
        "Noise Removal": {
            "fn": "_strat_noise_removal",
            "desc": "Three-stage Stratagem noise pipeline: (1) powerline notch "
            "filter suppresses 50/60 Hz harmonics via interpolation or "
            "NaN masking; (2) Hampel outlier filter removes spikes in the "
            "log-frequency domain; (3) optional log-freq triangular "
            "smoothing reduces residual high-freq scatter.",
            "params": [
                ParamSpec(
                    "mains_hz",
                    "Mains (Hz)",
                    "dspin",
                    50.0,
                    (40.0, 70.0, 1.0),
                    "Power-line fundamental (50 Hz Europe/Asia, 60 Hz Americas).",
                ),
                ParamSpec(
                    "n_harm",
                    "Harmonics",
                    "spin",
                    30,
                    (1, 100, 1),
                    "Number of harmonics to suppress.",
                ),
                ParamSpec(
                    "tol_hz",
                    "Tolerance (Hz)",
                    "dspin",
                    0.08,
                    (0.01, 2.0, 0.01),
                    "Half-width around each notch in Hz.",
                ),
                ParamSpec(
                    "notch_mode",
                    "Notch mode",
                    "combo",
                    "interp",
                    ["interp", "zero", "nan"],
                    "interp: fill by interpolation; zero: set to zero; "
                    "nan: leave as NaN.",
                ),
                ParamSpec(
                    "hampel_win",
                    "Hampel window",
                    "spin",
                    3,
                    (1, 15, 1),
                    "Half-window (number of freq bins) for Hampel filter.",
                ),
                ParamSpec(
                    "hampel_nsig",
                    "Hampel σ",
                    "dspin",
                    3.0,
                    (1.0, 10.0, 0.5),
                    "Outlier threshold in units of local MAD.",
                ),
                ParamSpec(
                    "smooth",
                    "Log-freq smooth",
                    "check",
                    False,
                    None,
                    "Apply optional log-frequency triangular smoothing.",
                ),
                ParamSpec(
                    "smooth_win",
                    "Smooth window",
                    "spin",
                    3,
                    (2, 21, 1),
                    "Half-window for log-freq smooth (only used if enabled).",
                ),
            ],
        },
        "Frequency Filter": {
            "fn": "_strat_freq_filter",
            "desc": "Three-stage band selection and incoherence masking: "
            "(1) optional hardware mask from Stratagem raw files (stack > 0); "
            "(2) hard band limits [f_min, f_max]; "
            "(3) statistical mask: drops frequency rows where SNR < "
            "threshold in > (1 − min_frac) of stations.  Zero values for "
            "limits are treated as 'no limit'.",
            "params": [
                ParamSpec(
                    "fmin",
                    "f min (Hz)",
                    "dspin",
                    0.0,
                    (0.0, 10000.0, 0.001),
                    "Low-frequency cut.  0 = no lower limit.",
                ),
                ParamSpec(
                    "fmax",
                    "f max (Hz)",
                    "dspin",
                    0.0,
                    (0.0, 10000.0, 100.0),
                    "High-frequency cut.  0 = no upper limit.",
                ),
                ParamSpec(
                    "snr_thresh",
                    "SNR threshold",
                    "dspin",
                    2.5,
                    (0.5, 20.0, 0.5),
                    "Minimum acceptable median SNR for a frequency row.",
                ),
                ParamSpec(
                    "min_frac",
                    "Min station frac",
                    "dspin",
                    0.4,
                    (0.1, 1.0, 0.05),
                    "Frequency kept only if OK in this fraction of stations.",
                ),
            ],
        },
        "Full Pipeline": {
            "fn": "_strat_full_pipeline",
            "desc": "Runs the complete Stratagem correction sequence in one step: "
            "QC → Frequency filter → Static-shift AMA → Noise removal.  "
            "All sub-step parameters use their default values shown above. "
            "Ideal for rapid first-pass processing of a new survey.",
            "params": [
                ParamSpec(
                    "mains_hz",
                    "Mains (Hz)",
                    "dspin",
                    50.0,
                    (40.0, 70.0, 1.0),
                    "Power-line fundamental for the noise-removal sub-step.",
                ),
                ParamSpec(
                    "half_window",
                    "SS half window",
                    "spin",
                    3,
                    (1, 20, 1),
                    "AMA half-window for the static-shift sub-step.",
                ),
                ParamSpec(
                    "max_skew",
                    "Max skew (°)",
                    "dspin",
                    6.0,
                    (0.5, 45.0, 0.5),
                    "Skew exclusion threshold shared by QC and SS sub-steps.",
                ),
            ],
        },
    },
}

# Category names and which ones use coordinate (map) view instead of ρ_a curves
CATEGORIES = list(CATALOGUE.keys())
COORD_CATEGORIES = {"Coordinates"}
ROTATION_CATEGORIES = {"Tensor Rotation"}
STRATAGEM_CATEGORIES = {"Stratagem"}
STATIC_SHIFT_CATEGORIES = {"Static Shift"}

# Coord correction function names (use DataFrame path, never touch EDI headers)
_COORD_FN_NAMES = {
    "_coord_profile_projection",
    "_coord_spacing_regularize",
    "_coord_outlier_snap",
    "_coord_elevation_smooth",
    "_coord_shift",
    "_coord_interpolate_missing",
}

# Stratagem function names (use edi_objects_ path, not Sites)
_STRAT_FN_NAMES = {
    "_strat_qc",
    "_strat_static_shift",
    "_strat_noise_removal",
    "_strat_freq_filter",
    "_strat_full_pipeline",
}


# ── Correction step ────────────────────────────────────────────────────────────


@dataclass
class CorrectionStep:
    label: str
    fn_name: str
    kwargs: dict[str, Any]
    sites_after: Any  # Sites — result (unchanged for coord steps)
    coords_df: Any = field(
        default=None
    )  # pd.DataFrame — only for coord corrections
    index: int = field(default=0)


# ── Controller ─────────────────────────────────────────────────────────────────


class CorrectionController:
    """
    Non-destructive correction chain.

    Coordinate corrections operate entirely on DataFrames (EDI headers are
    never touched during the workflow; they are only written at commit time).
    Impedance corrections use _apply_each which creates independent Z copies.
    """

    def __init__(self) -> None:
        self._raw_sites = None
        self._raw_coords_df = None  # snapshot taken at set_raw_sites
        self._preview_coords_df = (
            None  # result of last Preview for coord corrections
        )
        self._stack: list[CorrectionStep] = []
        self.dark: bool = True

        # ── Stratagem parallel state ──────────────────────────────────────────
        # edi_objects_ list is kept separately from the Sites chain so that
        # stratagem corrections operate natively on EDIFile objects and convert
        # back to Sites for display/commit.
        self._strat_edi_objects: list = []  # working (mutable) copy
        self._strat_raw_edi_objects: list = []  # immutable raw copy
        self._strat_edi_dir: str = ""
        self._strat_qc_report = None  # pd.DataFrame after run_qc
        self._strat_ss_factors = None  # pd.DataFrame after static shift

    # ── Data binding ──────────────────────────────────────────────────────────

    def set_raw_sites(self, sites) -> None:
        self._raw_sites = sites
        self._stack.clear()
        self._preview_coords_df = None
        try:
            from pycsamt.gis.coord_correction import (
                _get_coords_df,
            )

            self._raw_coords_df = _get_coords_df(sites).copy()
        except Exception:
            self._raw_coords_df = None

    def clear(self) -> None:
        self._raw_sites = None
        self._raw_coords_df = None
        self._preview_coords_df = None
        self._stack.clear()

    # ── Properties ────────────────────────────────────────────────────────────

    @property
    def has_data(self) -> bool:
        return self._raw_sites is not None

    @property
    def has_corrections(self) -> bool:
        return bool(self._stack)

    @property
    def has_coord_corrections(self) -> bool:
        return any(s.coords_df is not None for s in self._stack)

    @property
    def current_sites(self):
        return self._stack[-1].sites_after if self._stack else self._raw_sites

    @property
    def raw_sites(self):
        return self._raw_sites

    @property
    def raw_coords_df(self):
        return self._raw_coords_df

    def current_coords_df(self):
        """Coords after all applied coord corrections, or raw snapshot."""
        for step in reversed(self._stack):
            if step.coords_df is not None:
                return step.coords_df
        return self._raw_coords_df

    @property
    def stack(self) -> list[CorrectionStep]:
        return list(self._stack)

    @property
    def n_steps(self) -> int:
        return len(self._stack)

    # ── Correction operations ─────────────────────────────────────────────────

    def preview(self, fn_name: str, kwargs: dict) -> Any | None:
        """Compute correction without adding to the stack.
        For coord corrections returns the corrected DataFrame; otherwise Sites."""
        if fn_name in _STRAT_FN_NAMES:
            # Stratagem preview: run on a deep copy, don't mutate working state
            import copy

            saved_edis = copy.deepcopy(self._strat_edi_objects)
            saved_qc = self._strat_qc_report
            saved_factors = self._strat_ss_factors
            try:
                fn = getattr(self, fn_name)
                return fn(**kwargs)
            finally:
                self._strat_edi_objects = saved_edis
                self._strat_qc_report = saved_qc
                self._strat_ss_factors = saved_factors
        if self._raw_sites is None:
            return None
        if fn_name in _COORD_FN_NAMES:
            try:
                from pycsamt.gis.coord_correction import (
                    apply_coord_correction_to_df,
                )

                _cdf = self.current_coords_df()
                base_df = (
                    _cdf if _cdf is not None else self._raw_coords_df
                ).copy()
                corrected_df = apply_coord_correction_to_df(
                    fn_name, base_df, **kwargs
                )
                self._preview_coords_df = corrected_df
                return corrected_df  # window uses this as the "after" df
            except Exception:
                self._preview_coords_df = None
                return None
        try:
            result = self._call_fn(fn_name, self.current_sites, **kwargs)
            self._preview_coords_df = None
            return result
        except Exception:
            self._preview_coords_df = None
            raise  # re-raise so window can show the message

    def apply(
        self, fn_name: str, kwargs: dict, label: str
    ) -> CorrectionStep | None:
        """Compute correction and push it to the stack. Returns the new step."""
        if fn_name in _STRAT_FN_NAMES:
            try:
                fn = getattr(self, fn_name)
                result = fn(**kwargs)  # mutates _strat_edi_objects in-place
                step = CorrectionStep(
                    label=label,
                    fn_name=fn_name,
                    kwargs=dict(kwargs),
                    sites_after=result,
                    index=len(self._stack),
                )
                self._stack.append(step)
                self._raw_sites = (
                    self._raw_sites or result
                )  # bootstrap if needed
                return step
            except Exception:
                raise
        if self._raw_sites is None:
            return None
        if fn_name in _COORD_FN_NAMES:
            try:
                from pycsamt.gis.coord_correction import (
                    apply_coord_correction_to_df,
                )

                _cdf = self.current_coords_df()
                base_df = (
                    _cdf if _cdf is not None else self._raw_coords_df
                ).copy()
                corrected_df = apply_coord_correction_to_df(
                    fn_name, base_df, **kwargs
                )
                step = CorrectionStep(
                    label=label,
                    fn_name=fn_name,
                    kwargs=dict(kwargs),
                    sites_after=self.current_sites,  # Sites unchanged
                    coords_df=corrected_df,
                    index=len(self._stack),
                )
                self._stack.append(step)
                self._preview_coords_df = None
                return step
            except Exception:
                raise
        try:
            result = self._call_fn(fn_name, self.current_sites, **kwargs)
            step = CorrectionStep(
                label=label,
                fn_name=fn_name,
                kwargs=dict(kwargs),
                sites_after=result,
                index=len(self._stack),
            )
            self._stack.append(step)
            return step
        except Exception:
            return None

    def undo_last(self) -> None:
        if self._stack:
            self._stack.pop()
        self._preview_coords_df = None

    def remove_step(self, index: int) -> None:
        """Remove the step at *index* and replay all subsequent steps."""
        if index < 0 or index >= len(self._stack):
            return
        del self._stack[index]
        self._replay_from(index)
        self._preview_coords_df = None

    def revert_all(self) -> None:
        self._stack.clear()
        self._preview_coords_df = None

    # ── Stratagem EDI-directory loader ────────────────────────────────────────

    def load_edi_dir(self, edi_dir: str) -> int:
        """Load an EDI directory via EDIBatch; reset the correction stack.

        Returns the number of stations loaded.
        """
        import copy

        from pycsamt.stratagem.io import EDIBatch

        batch = EDIBatch().fit(edi_dir=edi_dir)
        self._strat_raw_edi_objects = copy.deepcopy(batch.edi_objects_)
        self._strat_edi_objects = batch.edi_objects_  # working copy
        self._strat_edi_dir = edi_dir
        self._strat_qc_report = None
        self._strat_ss_factors = None
        # Convert to Sites so the rest of the correction chain can work too
        sites = self._edi_objects_to_sites(self._strat_edi_objects)
        self.set_raw_sites(sites)
        return batch.n_stations_

    def reset_strat_to_raw(self) -> None:
        """Revert edi_objects_ to the raw state without touching the Sites stack."""
        import copy

        self._strat_edi_objects = copy.deepcopy(self._strat_raw_edi_objects)
        self._strat_qc_report = None
        self._strat_ss_factors = None

    @staticmethod
    def _edi_objects_to_sites(edi_objects):
        """Convert a list of EDIFile objects to a pycsamt Sites instance."""
        from pycsamt.emtools._core import ensure_sites

        return ensure_sites(
            edi_objects,
            recursive=True,
            on_dup="replace",
            strict=False,
            verbose=0,
        )

    # ── Stratagem correction methods ──────────────────────────────────────────

    def _strat_qc(
        self,
        min_frac_ok=0.6,
        min_snr_med=2.0,
        max_skew_med=6.0,
        include_skew=True,
        **_,
    ):
        """Run QualityController; stores report_; returns corrected Sites."""
        if not self._strat_edi_objects:
            raise RuntimeError(
                "No EDI directory loaded — use 'Load EDI Dir' first."
            )
        from pycsamt.stratagem.qc import QualityController

        qc = QualityController(
            min_frac_ok=float(min_frac_ok),
            min_snr_med=float(min_snr_med),
            max_skew_med=float(max_skew_med),
            include_skew=bool(include_skew),
        ).fit(self._strat_edi_objects)
        self._strat_qc_report = qc.report_
        # QC is diagnostic only — Sites are unchanged
        return self._edi_objects_to_sites(self._strat_edi_objects)

    def _strat_static_shift(
        self,
        sort_by="lon",
        half_window=3,
        weights="tri",
        pband_lo=0.0,
        pband_hi=0.0,
        max_skew=6.0,
        **_,
    ):
        """Apply StaticShiftCorrector; stores factors_; returns corrected Sites."""
        if not self._strat_edi_objects:
            raise RuntimeError(
                "No EDI directory loaded — use 'Load EDI Dir' first."
            )
        import copy

        from pycsamt.stratagem.process import (
            StaticShiftCorrector,
        )

        pband = None
        if float(pband_lo) > 0 and float(pband_hi) > float(pband_lo):
            pband = (float(pband_lo), float(pband_hi))
        corrector = StaticShiftCorrector(
            sort_by=sort_by,
            half_window=int(half_window),
            weights=weights,
            pband=pband,
            max_skew=float(max_skew),
        )
        edi_copy = copy.deepcopy(self._strat_edi_objects)
        corrector.fit(edi_copy, copy=False)
        self._strat_edi_objects = corrector.edi_objects_
        self._strat_ss_factors = corrector.factors_
        return self._edi_objects_to_sites(self._strat_edi_objects)

    def _strat_noise_removal(
        self,
        mains_hz=50.0,
        n_harm=30,
        tol_hz=0.08,
        notch_mode="interp",
        hampel_win=3,
        hampel_nsig=3.0,
        smooth=False,
        smooth_win=3,
        **_,
    ):
        """Apply NoiseRemover (notch + Hampel + optional smoothing)."""
        if not self._strat_edi_objects:
            raise RuntimeError(
                "No EDI directory loaded — use 'Load EDI Dir' first."
            )
        import copy

        from pycsamt.stratagem.process import NoiseRemover

        remover = NoiseRemover(
            mains_hz=float(mains_hz),
            n_harm=int(n_harm),
            tol_hz=float(tol_hz),
            notch_mode=notch_mode,
            hampel_win=int(hampel_win),
            hampel_nsig=float(hampel_nsig),
            smooth=bool(smooth),
            smooth_win=int(smooth_win),
        )
        edi_copy = copy.deepcopy(self._strat_edi_objects)
        remover.fit(edi_copy, copy=False)
        self._strat_edi_objects = remover.edi_objects_
        return self._edi_objects_to_sites(self._strat_edi_objects)

    def _strat_freq_filter(
        self, fmin=0.0, fmax=0.0, snr_thresh=2.5, min_frac=0.4, **_
    ):
        """Apply FrequencyFilter (band limits + incoherence mask)."""
        if not self._strat_edi_objects:
            raise RuntimeError(
                "No EDI directory loaded — use 'Load EDI Dir' first."
            )
        import copy

        from pycsamt.stratagem.qc import FrequencyFilter

        ff = FrequencyFilter(
            fmin=float(fmin) if float(fmin) > 0 else None,
            fmax=float(fmax) if float(fmax) > 0 else None,
            snr_thresh=float(snr_thresh),
            min_frac=float(min_frac),
        )
        edi_copy = copy.deepcopy(self._strat_edi_objects)
        ff.fit(edi_copy, copy=False)
        self._strat_edi_objects = ff.edi_objects_
        return self._edi_objects_to_sites(self._strat_edi_objects)

    def _strat_full_pipeline(
        self, mains_hz=50.0, half_window=3, max_skew=6.0, **_
    ):
        """Run the complete Stratagem pipeline: QC → filter → SS → noise."""
        if not self._strat_edi_objects:
            raise RuntimeError(
                "No EDI directory loaded — use 'Load EDI Dir' first."
            )
        import copy

        from pycsamt.stratagem.process import (
            NoiseRemover,
            StaticShiftCorrector,
        )
        from pycsamt.stratagem.qc import (
            FrequencyFilter,
            QualityController,
        )

        edis = copy.deepcopy(self._strat_edi_objects)

        # 1. QC (diagnostic only — does not modify edis)
        try:
            qc = QualityController(max_skew_med=float(max_skew)).fit(edis)
            self._strat_qc_report = qc.report_
        except Exception:
            pass

        # 2. Frequency filter
        try:
            ff = FrequencyFilter()
            ff.fit(edis, copy=False)
            edis = ff.edi_objects_
        except Exception:
            pass

        # 3. Static shift
        try:
            ss = StaticShiftCorrector(
                half_window=int(half_window), max_skew=float(max_skew)
            )
            ss.fit(edis, copy=False)
            edis = ss.edi_objects_
            self._strat_ss_factors = ss.factors_
        except Exception:
            pass

        # 4. Noise removal
        try:
            nr = NoiseRemover(mains_hz=float(mains_hz))
            nr.fit(edis, copy=False)
            edis = nr.edi_objects_
        except Exception:
            pass

        self._strat_edi_objects = edis
        return self._edi_objects_to_sites(edis)

    def export_stratagem(self, savepath: str, overwrite: bool = False) -> int:
        """Write the current (corrected) edi_objects_ to *savepath*."""
        if not self._strat_edi_objects:
            raise RuntimeError("No EDI objects to export.")
        from pathlib import Path

        # Build a minimal EDIBatch-like exporter directly
        out_path = Path(savepath)
        out_path.mkdir(parents=True, exist_ok=True)
        written = 0
        for edi in self._strat_edi_objects:
            try:
                name = getattr(edi, "station", None) or f"S{written:03d}"
                dest = out_path / f"{name}.edi"
                if dest.exists() and not overwrite:
                    dest = out_path / f"{name}_corr.edi"
                edi.write_edifile(
                    savepath=str(out_path), new_edifilename=dest.stem
                )
                written += 1
            except Exception:
                pass
        return written

    # ── Stratagem plot helpers ────────────────────────────────────────────────

    def plot_strat_qc(self, fig) -> None:
        """Draw QC report heatmap (frac_ok + flag marks) on *fig*."""
        s = _DARK if self.dark else _LIGHT
        fig.patch.set_facecolor(s["fig_bg"])
        fig.clf()
        if self._strat_qc_report is None or self._strat_qc_report.empty:
            ax = fig.add_subplot(111)
            ax.set_facecolor(s["bg"])
            ax.set_axis_off()
            ax.text(
                0.5,
                0.5,
                "Run QC Report first",
                transform=ax.transAxes,
                ha="center",
                va="center",
                color=s["muted"],
                fontsize=11,
            )
            return
        df = self._strat_qc_report
        stations = (
            df["station"].values
            if "station" in df.columns
            else np.arange(len(df))
        )
        n = len(stations)
        from pycsamt.compat.matplotlib import get_cmap

        cols = ["frac_ok", "snr_med"]
        if "skew_med" in df.columns:
            cols.append("skew_med")
        n_cols = len(cols)
        axes = fig.subplots(1, n_cols, sharey=True)
        if n_cols == 1:
            axes = [axes]
        labels = {
            "frac_ok": "Valid fraction",
            "snr_med": "Median SNR",
            "skew_med": "Median skew (°)",
        }
        cmaps = {
            "frac_ok": "RdYlGn",
            "snr_med": "plasma",
            "skew_med": "RdYlBu_r",
        }
        for ax, col in zip(axes, cols):
            vals = df[col].values.astype(float)
            ax.barh(
                np.arange(n),
                vals,
                color=get_cmap(cmaps[col])(
                    vals / (np.nanmax(vals) + 1e-12)
                ),
                edgecolor=s["spine"],
                lw=0.5,
            )
            ax.set_facecolor(s["bg"])
            ax.tick_params(colors=s["tick"], labelsize=6)
            ax.set_title(labels[col], fontsize=7, color=s["fg"])
            ax.set_xlabel(col, fontsize=6, color=s["fg"])
            for sp in ax.spines.values():
                sp.set_edgecolor(s["spine"])
            # Flag markers
            if "flagged" in df.columns:
                for i, (flag, v) in enumerate(
                    zip(df["flagged"].values, vals)
                ):
                    if flag:
                        ax.text(
                            v,
                            i,
                            " ✗",
                            va="center",
                            ha="left",
                            color="#f38ba8",
                            fontsize=7,
                        )
        axes[0].set_yticks(np.arange(n))
        axes[0].set_yticklabels(
            [str(s_)[:12] for s_ in stations], fontsize=5, color=s["tick"]
        )
        fig.suptitle("Stratagem QC Report", fontsize=9, color=s["fg"])
        fig.tight_layout()

    def plot_strat_ss_factors(self, ax) -> None:
        """Draw per-station static-shift factors on *ax*."""
        s = _DARK if self.dark else _LIGHT
        ax.cla()
        ax.set_facecolor(s["bg"])
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.tick_params(colors=s["tick"], labelsize=7)
        if self._strat_ss_factors is None or self._strat_ss_factors.empty:
            ax.text(
                0.5,
                0.5,
                "Run Static Shift first",
                transform=ax.transAxes,
                ha="center",
                color=s["muted"],
            )
            return
        df = self._strat_ss_factors
        x = np.arange(len(df))
        key = (
            "delta_log10_rho"
            if "delta_log10_rho" in df.columns
            else df.columns[1]
        )
        vals = df[key].values.astype(float)
        colors = ["#f38ba8" if v < 0 else "#a6e3a1" for v in vals]
        ax.bar(x, vals, color=colors, edgecolor=s["spine"], linewidth=0.5)
        ax.axhline(0, color=s["spine"], lw=0.8, ls="--")
        ax.set_xlabel("Station index", fontsize=8, color=s["fg"])
        ax.set_ylabel("Δ log₁₀ ρ", fontsize=8, color=s["fg"])
        ax.set_title("Static-Shift Factors", fontsize=9, color=s["fg"])
        ax.xaxis.label.set_color(s["fg"])
        ax.yaxis.label.set_color(s["fg"])
        ax.title.set_color(s["fg"])

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _replay_from(self, index: int) -> None:
        # Coord steps don't need replay (they use cumulative DataFrame ops)
        # Only replay impedance/rotation steps
        base = (
            self._stack[index - 1].sites_after
            if index > 0
            else self._raw_sites
        )
        for step in self._stack[index:]:
            if step.fn_name in _COORD_FN_NAMES:
                step.sites_after = base  # sites unchanged for coord steps
                continue
            try:
                step.sites_after = self._call_fn(
                    step.fn_name, base, **step.kwargs
                )
            except Exception:
                step.sites_after = base
            base = step.sites_after

    def _call_fn(self, fn_name: str, sites, **kwargs):
        """Dispatch to emtools function or internal wrapper."""
        # Stratagem wrappers (bypass normal sites path)
        if fn_name in _STRAT_FN_NAMES:
            fn = getattr(self, fn_name)
            return fn(**kwargs)
        # Internal wrappers (not in emtools public API)
        _wrappers = {
            # Static-shift wrappers
            "_correct_ss_loess": self._wrap_ss_loess,
            "_correct_ss_bilateral": self._wrap_ss_bilateral,
            "_correct_ss_refmedian": self._wrap_ss_refmedian,
            "_correct_near_field": self._wrap_near_field,
            # Tensor-rotation wrappers
            "_wrap_rotate": self._wrap_rotate,
            "_wrap_rotate_to_strike": self._wrap_rotate_to_strike,
            "_wrap_rotate_pt_strike": self._wrap_rotate_pt_strike,
            "_wrap_rotate_to_profile": self._wrap_rotate_to_profile,
            # Coordinate-correction wrappers
            "_coord_profile_projection": self._coord_profile_projection,
            "_coord_spacing_regularize": self._coord_spacing_regularize,
            "_coord_outlier_snap": self._coord_outlier_snap,
            "_coord_elevation_smooth": self._coord_elevation_smooth,
            "_coord_shift": self._coord_shift,
            "_coord_interpolate_missing": self._coord_interpolate_missing,
        }
        if fn_name in _wrappers:
            return _wrappers[fn_name](sites, **kwargs)

        import pycsamt.emtools as et

        fn = getattr(et, fn_name)
        return fn(sites, inplace=False, verbose=0, **kwargs)

    # ── Correction wrappers (for functions without one-shot correct_* entry) ──

    @staticmethod
    def _wrap_ss_loess(sites, half_window=3, poly=1, it=2, max_skew=6.0, **_):
        import pycsamt.emtools as et

        tbl = et.estimate_ss_loess(
            sites,
            half_window=half_window,
            poly=poly,
            it=it,
            max_skew=max_skew,
            verbose=0,
        )
        return et.apply_ss_factors(
            sites, tbl, key="fac_z", inplace=False, verbose=0
        )

    @staticmethod
    def _wrap_ss_bilateral(sites, half_window=4, max_skew=6.0, **_):
        import pycsamt.emtools as et

        tbl = et.estimate_ss_bilateral(
            sites,
            half_window=half_window,
            max_skew=max_skew,
            verbose=0,
        )
        return et.apply_ss_factors(
            sites, tbl, key="fac_z", inplace=False, verbose=0
        )

    @staticmethod
    def _wrap_ss_refmedian(sites, smooth_sites=0, max_skew=6.0, **_):
        import pycsamt.emtools as et

        tbl = et.estimate_ss_refmedian(
            sites,
            smooth_sites=smooth_sites,
            max_skew=max_skew,
            verbose=0,
        )
        return et.apply_ss_factors(
            sites, tbl, key="fac_z", inplace=False, verbose=0
        )

    @staticmethod
    def _wrap_near_field(sites, source_offset=500.0, **_):
        import pycsamt.emtools as et

        return et.correct_near_field(
            sites, source_offset, inplace=False, verbose=0
        )

    # ── Tensor rotation wrappers ──────────────────────────────────────────────

    @staticmethod
    def _wrap_rotate(sites, angle=0.0, **_):
        # _edit.rotate works on a single EDI; _edit.rotate_all broadcasts
        # over a Sites collection and returns a new Sites.
        from pycsamt.emtools._core import ensure_sites
        from pycsamt.site import edit as _edit

        S = ensure_sites(
            sites, recursive=True, on_dup="replace", strict=False, verbose=0
        )
        return _edit.rotate_all(S, float(angle), inplace=False)

    @staticmethod
    def _wrap_rotate_to_strike(sites, method="swift", **_):
        # Estimate per-station strike then rotate each site individually.
        from pycsamt.emtools._core import (
            _iter_items,
            ensure_sites,
        )
        from pycsamt.site import compute as _compute
        from pycsamt.site import edit as _edit
        from pycsamt.site.base import Sites

        S = ensure_sites(
            sites, recursive=True, on_dup="replace", strict=False, verbose=0
        )
        rotated_items = []
        for _i, ed in enumerate(_iter_items(S)):
            Si = ensure_sites(
                [getattr(ed, "edi", ed)], recursive=False, strict=False
            )
            ang = 0.0
            try:
                r = _compute.strike_estimate(Si, method=method)
                ang = (
                    float(r)
                    if isinstance(r, (int, float))
                    else float(
                        r.get("angle", 0.0)
                        if isinstance(r, dict)
                        else getattr(r, "angle", 0.0)
                    )
                )
            except Exception:
                ang = 0.0
            ri = _edit.rotate_all(Si, ang, inplace=False)
            rotated_items.extend(list(_iter_items(ri)))
        # Rebuild Sites from rotated EDIs
        from pycsamt.seg.collection import EDICollection

        raw_edis = [getattr(it, "edi", it) for it in rotated_items]
        try:
            return Sites(EDICollection(items=raw_edis))
        except TypeError:
            return Sites(edic=EDICollection(items=raw_edis))

    @staticmethod
    def _wrap_rotate_pt_strike(
        sites, band_lo=0.001, band_hi=1000.0, robust=True, **_
    ):
        import pycsamt.emtools as et
        from pycsamt.emtools._core import (
            _iter_items,
            _name,
            ensure_sites,
        )
        from pycsamt.seg.collection import EDICollection
        from pycsamt.site import edit as _edit
        from pycsamt.site.base import Sites

        S = ensure_sites(
            sites, recursive=True, on_dup="replace", strict=False, verbose=0
        )
        band = (float(band_lo), float(band_hi))
        tbl = et.estimate_strike_phase_tensor(
            S, band=band, robust=bool(robust), verbose=0
        )
        if tbl.empty:
            return S
        angle_map = dict(zip(tbl["station"], tbl["ang"]))

        rotated_items = []
        for i, ed in enumerate(_iter_items(S)):
            st = _name(ed, i)
            ang = float(angle_map.get(st, 0.0))
            Si = ensure_sites(
                [getattr(ed, "edi", ed)], recursive=False, strict=False
            )
            ri = _edit.rotate_all(Si, ang, inplace=False)
            rotated_items.extend(list(_iter_items(ri)))

        raw_edis = [getattr(it, "edi", it) for it in rotated_items]
        try:
            return Sites(EDICollection(items=raw_edis))
        except TypeError:
            return Sites(edic=EDICollection(items=raw_edis))

    @staticmethod
    def _wrap_rotate_to_profile(sites, azimuth=-1.0, **_):
        from pycsamt.emtools._core import ensure_sites
        from pycsamt.gis.coord_correction import (
            _get_coords_df,
            _pca_azimuth,
            _to_utm_arrays,
        )
        from pycsamt.site import edit as _edit

        S = ensure_sites(
            sites, recursive=True, on_dup="replace", strict=False, verbose=0
        )
        if float(azimuth) < 0:
            df = _get_coords_df(S)
            v = df.dropna(subset=["lat", "lon"])
            if len(v) >= 2:
                east, north, _ = _to_utm_arrays(
                    v["lat"].values, v["lon"].values
                )
                azimuth = _pca_azimuth(east, north)
            else:
                azimuth = 0.0
        return _edit.rotate_all(S, -float(azimuth), inplace=False)

    # ── Coordinate correction wrappers ────────────────────────────────────────

    @staticmethod
    def _coord_profile_projection(
        sites, azimuth=-1.0, keep_elevation=True, **_
    ):
        from pycsamt.gis.coord_correction import (
            correct_profile_projection,
        )

        az = None if azimuth < 0 else float(azimuth)
        return correct_profile_projection(
            sites,
            azimuth=az,
            keep_elevation=bool(keep_elevation),
            inplace=False,
        )

    @staticmethod
    def _coord_spacing_regularize(
        sites, spacing_m=200.0, azimuth=-1.0, preserve_extent=True, **_
    ):
        from pycsamt.gis.coord_correction import (
            correct_spacing_regularize,
        )

        az = None if azimuth < 0 else float(azimuth)
        return correct_spacing_regularize(
            sites,
            spacing_m=float(spacing_m),
            azimuth=az,
            preserve_extent=bool(preserve_extent),
            inplace=False,
        )

    @staticmethod
    def _coord_outlier_snap(sites, threshold_m=50.0, azimuth=-1.0, **_):
        from pycsamt.gis.coord_correction import (
            correct_outlier_snap,
        )

        az = None if azimuth < 0 else float(azimuth)
        return correct_outlier_snap(
            sites, threshold_m=float(threshold_m), azimuth=az, inplace=False
        )

    @staticmethod
    def _coord_elevation_smooth(sites, method="loess", window=5, **_):
        from pycsamt.gis.coord_correction import (
            correct_elevation_smooth,
        )

        return correct_elevation_smooth(
            sites, method=str(method), window=int(window), inplace=False
        )

    @staticmethod
    def _coord_shift(
        sites, delta_lat=0.0, delta_lon=0.0, delta_elev=0.0, **_
    ):
        from pycsamt.gis.coord_correction import (
            correct_coordinate_shift,
        )

        return correct_coordinate_shift(
            sites,
            delta_lat=float(delta_lat),
            delta_lon=float(delta_lon),
            delta_elev=float(delta_elev),
            inplace=False,
        )

    @staticmethod
    def _coord_interpolate_missing(sites, fill_nan_only=True, **_):
        from pycsamt.gis.coord_correction import (
            correct_interpolate_missing,
        )

        return correct_interpolate_missing(
            sites, fill_nan_only=bool(fill_nan_only), inplace=False
        )

    # ── Plotting ──────────────────────────────────────────────────────────────

    def plot_rho_curves(self, sites, ax, title: str = "") -> None:
        """Multi-line log(ρ_a xy) vs log(T) for all stations."""
        s = _DARK if self.dark else _LIGHT
        ax.set_facecolor(s["bg"])
        fig = ax.get_figure()
        if fig:
            fig.patch.set_facecolor(s["fig_bg"])
        ax.clear()
        ax.set_facecolor(s["bg"])

        if sites is None:
            _empty(ax, "No data loaded", s)
            return

        try:
            from pycsamt.emtools._core import (
                _get_z_block,
                _iter_items,
            )

            items = list(_iter_items(sites))
            n = len(items)
            if n == 0:
                _empty(ax, "Empty dataset", s)
                return
            import matplotlib.pyplot as plt

            cmap = plt.cm.turbo
            plotted = 0
            for idx, ed in enumerate(items):
                _, z, freqs = _get_z_block(ed)
                if z is None or freqs is None or freqs.size == 0:
                    continue
                Z_xy = z[:, 0, 1]
                # ρ_a = 0.2·|Z|²/f — field-unit form (Z is mV/km/nT). The SI
                # form |Z|²/(ω·μ₀) over-estimates ρ_a by ~6.3·10⁵ here.
                rho = 0.2 * np.abs(Z_xy) ** 2 / np.maximum(freqs, 1e-30)
                T = 1.0 / np.maximum(freqs, 1e-30)
                ok = np.isfinite(rho) & (rho > 0) & np.isfinite(T)
                if ok.sum() < 2:
                    continue
                ax.loglog(
                    T[ok],
                    rho[ok],
                    color=cmap(idx / max(n - 1, 1)),
                    alpha=0.55,
                    lw=0.9,
                )
                plotted += 1
            if plotted == 0:
                _empty(ax, "No valid Z data", s)
                return
        except Exception as exc:
            _empty(ax, f"Error: {exc}", s)
            return

        ax.set_xlabel("Period (s)", fontsize=8, color=s["fg"])
        ax.set_ylabel("ρ_a  [Ω·m]", fontsize=8, color=s["fg"])
        ax.set_title(title, fontsize=9, color=s["title"], pad=5)
        ax.tick_params(colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.grid(
            True, which="both", color=s["grid"], alpha=0.25, ls="--", lw=0.5
        )

    def plot_rho_pseudosection(
        self,
        sites,
        fig,
        *,
        affected_stations: list | None = None,
        title: str = "ρ_a Pseudosection  (XY)",
    ) -> None:
        """Period × station pcolormesh of log₁₀(ρ_a XY).

        Draws into *fig* (cleared first).  If ``affected_stations`` is given,
        vertical lines mark those stations.  When ``PYCSAMT_TOPO.enabled`` is
        True a terrain strip is added above the pcolormesh.
        """
        s = _DARK if self.dark else _LIGHT
        fig.clear()
        fig.patch.set_facecolor(s["fig_bg"])

        if sites is None:
            ax = fig.add_subplot(111)
            _init_ax(ax, s)
            _empty(ax, "No data loaded", s)
            return

        try:
            from pycsamt.emtools._core import (
                _get_z_block,
                _iter_items,
                _name,
            )

            items = list(_iter_items(sites))
            n_st = len(items)
            if n_st == 0:
                ax = fig.add_subplot(111)
                _init_ax(ax, s)
                _empty(ax, "Empty dataset", s)
                return

            station_names: list[str] = []
            station_data: list = []  # (T_sorted, rho_sorted) | None

            for i, ed in enumerate(items):
                station_names.append(_name(ed, i))
                _, z, freqs = _get_z_block(ed)
                if z is None or freqs is None or freqs.size == 0:
                    station_data.append(None)
                    continue
                Z_xy = z[:, 0, 1]
                # field-unit ρ_a = 0.2·|Z|²/f (see plot_rho_curves note)
                rho = 0.2 * np.abs(Z_xy) ** 2 / np.maximum(freqs, 1e-30)
                T = 1.0 / np.maximum(freqs, 1e-30)
                ok = np.isfinite(rho) & (rho > 0) & np.isfinite(T)
                if ok.sum() < 2:
                    station_data.append(None)
                    continue
                order = np.argsort(T[ok])
                station_data.append((T[ok][order], rho[ok][order]))

            valid = [d for d in station_data if d is not None]
            if not valid:
                ax = fig.add_subplot(111)
                _init_ax(ax, s)
                _empty(ax, "No valid Z data", s)
                return

            all_T = np.concatenate([d[0] for d in valid])
            T_lo = float(np.nanpercentile(all_T, 1))
            T_hi = float(np.nanpercentile(all_T, 99))
            n_T = 64
            T_grid = np.logspace(np.log10(T_lo), np.log10(T_hi), n_T)

            rho_2d = np.full((n_T, n_st), np.nan)
            for i, data in enumerate(station_data):
                if data is None:
                    continue
                T_s, rho_s = data
                rho_2d[:, i] = np.interp(
                    np.log10(T_grid),
                    np.log10(np.maximum(T_s, 1e-30)),
                    np.log10(np.maximum(rho_s, 1e-10)),
                    left=np.nan,
                    right=np.nan,
                )
        except Exception as exc:
            ax = fig.add_subplot(111)
            _init_ax(ax, s)
            _empty(ax, f"Error: {exc}", s)
            return

        # ── pcolormesh ────────────────────────────────────────────────────────
        import numpy.ma as ma

        ax = fig.add_subplot(111)
        _init_ax(ax, s)

        x_edges = np.arange(n_st + 1, dtype=float) - 0.5
        lT = np.log10(T_grid)
        dlT = np.diff(lT)
        T_edges = 10 ** np.concatenate(
            [
                [lT[0] - dlT[0] / 2],
                (lT[:-1] + lT[1:]) / 2,
                [lT[-1] + dlT[-1] / 2],
            ]
        )

        finite = rho_2d[np.isfinite(rho_2d)]
        vmin = float(np.percentile(finite, 2)) if finite.size else 0.0
        vmax = float(np.percentile(finite, 98)) if finite.size else 4.0

        im = ax.pcolormesh(
            x_edges,
            T_edges,
            ma.masked_invalid(rho_2d),
            shading="flat",
            cmap="jet",
            vmin=vmin,
            vmax=vmax,
        )

        # y-axis: shorter period (higher freq, near-surface) at TOP,
        # longer period (lower freq, deep) at BOTTOM — geophysics convention.
        ax.set_yscale("log")
        ax.set_ylim(T_edges[0], T_edges[-1])  # small T at bottom …
        ax.invert_yaxis()  # … flip: small T → visual top
        ax.set_xlim(-0.5, n_st - 0.5)

        # ── Station axis at TOP (conventional pseudosection layout) ──────────
        # Move the x-axis (stations) to the top so they sit above the pcolormesh,
        # which physically corresponds to the surface.
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        tick_step = max(1, n_st // 12)
        tick_idx = np.arange(0, n_st, tick_step)
        ax.set_xticks(tick_idx)
        ax.set_xticklabels(
            [station_names[i] for i in tick_idx],
            rotation=45,
            ha="left",
            fontsize=6,
        )
        ax.set_xlabel("Station", fontsize=8, color=s["fg"])

        # Downward-triangle markers (▼) at the top edge of each station column.
        # Uses the blended xaxis transform: x in data coords, y in axes fraction.
        _trans = ax.get_xaxis_transform()
        ax.plot(
            np.arange(n_st),
            np.ones(n_st),
            linestyle="none",
            marker="v",
            ms=4,
            color=s["tick"],
            clip_on=False,
            zorder=10,
            transform=_trans,
        )

        # ── Colorbar ──────────────────────────────────────────────────────────
        try:
            cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.025)
            cbar.set_label("log₁₀(ρ_a  [Ω·m])", fontsize=7, color=s["fg"])
            cbar.ax.tick_params(labelsize=6, colors=s["tick"])
        except Exception:
            pass

        # ── Affected-station highlights ───────────────────────────────────────
        if affected_stations:
            aff_set = set(affected_stations)
            first = True
            for i, nm in enumerate(station_names):
                if nm in aff_set:
                    ax.axvline(
                        i,
                        color="#f38ba8",
                        lw=1.5,
                        alpha=0.85,
                        zorder=5,
                        label="affected" if first else None,
                    )
                    first = False
            if not first:
                ax.legend(
                    fontsize=7,
                    framealpha=0.3,
                    labelcolor=s["fg"],
                    facecolor=s["bg"],
                    loc="lower right",
                )

        ax.set_ylabel("Period  (s)", fontsize=8, color=s["fg"])
        ax.set_title(title, fontsize=9, color=s["title"], pad=14)
        ax.tick_params(which="both", colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])

        # ── Optional topo strip ───────────────────────────────────────────────
        # When PYCSAMT_TOPO.enabled the terrain profile + station pins are drawn
        # in a strip ABOVE the pcolormesh.  The station x-axis labels of the main
        # axes are then hidden to avoid duplication with the strip's own labels.
        try:
            from pycsamt.topo.config import PYCSAMT_TOPO

            if PYCSAMT_TOPO.enabled:
                from pycsamt.topo.extract import (
                    extract_elevation,
                )
                from pycsamt.topo.extract import (
                    extract_station_names as _ext_names,
                )
                from pycsamt.topo.overlay import (
                    draw_topo_strip,
                )

                elev = extract_elevation(sites)
                names = _ext_names(sites)
                if len(elev) == n_st and np.any(elev > 0):
                    chain_idx = np.arange(n_st, dtype=float)
                    ax_topo = draw_topo_strip(
                        fig,
                        ax,
                        chain_idx,
                        elev,
                        names,
                        dark=self.dark,
                    )
                    if ax_topo is not None:
                        ax_topo.set_xlim(-0.5, n_st - 0.5)
                    # Topo strip sits above — hide the main-axes station labels
                    # but keep the ▼ pin markers for visual continuity.
                    ax.tick_params(axis="x", which="both", labeltop=False)
                    ax.set_xlabel("")
        except Exception:
            pass

    def plot_overlay(self, before, after, ax) -> None:
        """Before (thin dashed) + after (solid) on the same axes."""
        s = _DARK if self.dark else _LIGHT
        ax.set_facecolor(s["bg"])
        fig = ax.get_figure()
        if fig:
            fig.patch.set_facecolor(s["fig_bg"])
        ax.clear()
        ax.set_facecolor(s["bg"])

        try:
            from pycsamt.emtools._core import (
                _get_z_block,
                _iter_items,
            )

            def _draw(sites, ls, alpha, lw):
                items = list(_iter_items(sites))
                n = max(len(items), 1)
                import matplotlib.pyplot as plt

                cmap = plt.cm.turbo
                for idx, ed in enumerate(items):
                    _, z, freqs = _get_z_block(ed)
                    if z is None or freqs is None or freqs.size == 0:
                        continue
                    Z_xy = z[:, 0, 1]
                    # field-unit ρ_a = 0.2·|Z|²/f (see plot_rho_curves note)
                    rho = 0.2 * np.abs(Z_xy) ** 2 / np.maximum(freqs, 1e-30)
                    T = 1.0 / np.maximum(freqs, 1e-30)
                    ok = np.isfinite(rho) & (rho > 0)
                    if ok.sum() < 2:
                        continue
                    ax.loglog(
                        T[ok],
                        rho[ok],
                        color=cmap(idx / (n - 1) if n > 1 else 0),
                        alpha=alpha,
                        lw=lw,
                        ls=ls,
                    )

            _draw(before, "--", 0.35, 0.8)
            _draw(after, "-", 0.75, 1.1)

            from matplotlib.lines import Line2D

            ax.legend(
                [
                    Line2D(
                        [0], [0], color="white", ls="--", lw=0.8, alpha=0.6
                    ),
                    Line2D(
                        [0], [0], color="white", ls="-", lw=1.1, alpha=0.9
                    ),
                ],
                ["Before", "After"],
                fontsize=7,
                loc="lower left",
                facecolor=s["bg"],
                edgecolor=s["spine"],
                labelcolor=s["fg"],
            )
        except Exception as exc:
            _empty(ax, f"Overlay error: {exc}", s)
            return

        ax.set_xlabel("Period (s)", fontsize=8, color=s["fg"])
        ax.set_ylabel("ρ_a  [Ω·m]", fontsize=8, color=s["fg"])
        ax.set_title(
            "Before / After overlay", fontsize=9, color=s["title"], pad=5
        )
        ax.tick_params(colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.grid(
            True, which="both", color=s["grid"], alpha=0.25, ls="--", lw=0.5
        )

    def plot_diff(self, before, after, ax) -> None:
        """Relative change (ρ_after − ρ_before) / ρ_before × 100 % per station."""
        s = _DARK if self.dark else _LIGHT
        ax.set_facecolor(s["bg"])
        fig = ax.get_figure()
        if fig:
            fig.patch.set_facecolor(s["fig_bg"])
        ax.clear()
        ax.set_facecolor(s["bg"])

        if before is None or after is None:
            _empty(ax, "Need both before and after to show diff", s)
            return

        try:
            from pycsamt.emtools._core import (
                _get_z_block,
                _iter_items,
            )

            b_items = list(_iter_items(before))
            a_items = list(_iter_items(after))
            n = min(len(b_items), len(a_items))
            if n == 0:
                _empty(ax, "No common stations", s)
                return

            T_all, delta_all, idx_all = [], [], []
            for i in range(n):
                _, zb, fb = _get_z_block(b_items[i])
                _, za, fa = _get_z_block(a_items[i])
                if zb is None or za is None or fb is None or fa is None:
                    continue
                # field-unit ρ_a = 0.2·|Z|²/f (see plot_rho_curves note)
                rho_b = 0.2 * np.abs(zb[:, 0, 1]) ** 2 / np.maximum(fb, 1e-30)
                rho_a_arr = (
                    0.2 * np.abs(za[:, 0, 1]) ** 2 / np.maximum(fa, 1e-30)
                )
                T = 1.0 / np.maximum(fb, 1e-30)
                ok = np.isfinite(rho_b) & (rho_b > 0) & np.isfinite(rho_a_arr)
                if ok.sum() < 2:
                    continue
                delta = 100.0 * (rho_a_arr[ok] - rho_b[ok]) / rho_b[ok]
                T_all.extend(T[ok])
                delta_all.extend(delta)
                idx_all.extend([i] * ok.sum())

            if not T_all:
                _empty(ax, "No valid diff data", s)
                return

            sc = ax.scatter(
                T_all,
                idx_all,
                c=delta_all,
                s=4,
                cmap="RdBu_r",
                vmin=-50,
                vmax=50,
                alpha=0.7,
            )
            cb = ax.get_figure().colorbar(sc, ax=ax, pad=0.02)
            cb.set_label("Δρ_a / ρ_a [%]", fontsize=7, color=s["fg"])
            cb.ax.yaxis.set_tick_params(color=s["tick"], labelcolor=s["tick"])

        except Exception as exc:
            _empty(ax, f"Diff error: {exc}", s)
            return

        ax.set_xscale("log")
        ax.set_xlabel("Period (s)", fontsize=8, color=s["fg"])
        ax.set_ylabel("Station index", fontsize=8, color=s["fg"])
        ax.set_title(
            "Relative change (ρ_after − ρ_before) / ρ_before",
            fontsize=9,
            color=s["title"],
            pad=5,
        )
        ax.tick_params(colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.grid(True, axis="x", color=s["grid"], alpha=0.25, ls="--", lw=0.5)

    # ── Coordinate / map plots ────────────────────────────────────────────────

    def _df_from(self, data) -> pd.DataFrame | None:
        """Normalise *data* to a coords DataFrame."""
        import pandas as pd

        if isinstance(data, pd.DataFrame):
            return data
        if data is None:
            return None
        try:
            from pycsamt.gis.coord_correction import (
                _get_coords_df,
            )

            return _get_coords_df(data)
        except Exception:
            return None

    def plot_station_map(self, data, ax, title: str = "") -> None:
        """Scatter of station positions.  *data* may be a coords DataFrame or Sites."""
        from pycsamt.gis.coord_correction import (
            _to_utm_arrays,
        )

        s = _DARK if self.dark else _LIGHT
        _init_ax(ax, s)

        df = self._df_from(data)
        if df is None or df.empty:
            _empty(ax, "No data loaded", s)
            return
        valid = df.dropna(subset=["lat", "lon"])
        if valid.empty:
            _empty(ax, "No valid coordinates", s)
            return

        try:
            east, north, _zone = _to_utm_arrays(
                valid["lat"].values, valid["lon"].values
            )
            n = len(east)
            import matplotlib.pyplot as plt

            cmap = plt.cm.turbo
            colors = [cmap(i / max(n - 1, 1)) for i in range(n)]
            stations = (
                valid["station"].values
                if "station" in valid.columns
                else [str(i + 1) for i in range(n)]
            )

            # Draw connecting line first (profile trend)
            if n >= 2:
                from numpy.polynomial import polynomial as P

                c = P.polyfit(east, north, 1)
                xs = np.linspace(east.min(), east.max(), 100)
                ax.plot(
                    xs,
                    P.polyval(xs, c),
                    color=s["spine"],
                    lw=0.8,
                    ls="--",
                    alpha=0.45,
                    zorder=1,
                )

            # Station markers + labels
            offset_e = (east.max() - east.min()) * 0.015 + 1.0
            offset_n = (north.max() - north.min()) * 0.015 + 1.0
            for i, (e, no, st) in enumerate(zip(east, north, stations)):
                ax.scatter(
                    e,
                    no,
                    color=colors[i],
                    s=60,
                    zorder=3,
                    edgecolors=s["spine"],
                    linewidths=0.6,
                    marker="o",
                )
                # Abbreviated label — show last part of station name
                short = st.split("_")[-1] if "_" in st else st
                ax.annotate(
                    short,
                    xy=(e, no),
                    xytext=(e + offset_e, no + offset_n),
                    fontsize=6,
                    color=s["fg"],
                    arrowprops=dict(
                        arrowstyle="-", color=s["spine"], lw=0.4, alpha=0.5
                    ),
                    clip_on=True,
                    zorder=4,
                )
        except Exception as exc:
            _empty(ax, f"Map error: {exc}", s)
            return

        ax.set_xlabel("Easting (m)", fontsize=8, color=s["fg"])
        ax.set_ylabel("Northing (m)", fontsize=8, color=s["fg"])
        ax.set_title(title, fontsize=9, color=s["title"], pad=5)
        ax.tick_params(colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.grid(True, color=s["grid"], alpha=0.25, ls="--", lw=0.5)
        ax.set_aspect("equal", adjustable="datalim")

    def plot_station_map_overlay(self, before_data, after_data, ax) -> None:
        """Before (gray ▲) + after (coloured ●) on same axes, with station names."""
        from pycsamt.gis.coord_correction import (
            _to_utm_arrays,
        )

        s = _DARK if self.dark else _LIGHT
        _init_ax(ax, s)

        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        def _draw(data, marker, color_fn, alpha, ls_line):
            df = self._df_from(data)
            if df is None or df.empty:
                return [], []
            v = df.dropna(subset=["lat", "lon"])
            if v.empty:
                return [], []
            e, n, _ = _to_utm_arrays(v["lat"].values, v["lon"].values)
            stations = (
                v["station"].values
                if "station" in v.columns
                else [str(i + 1) for i in range(len(e))]
            )
            n_st = len(e)
            for i, (ei, ni, _st) in enumerate(zip(e, n, stations)):
                ax.scatter(
                    ei,
                    ni,
                    marker=marker,
                    color=color_fn(i, n_st),
                    alpha=alpha,
                    s=55,
                    edgecolors=s["spine"],
                    linewidths=0.5,
                    zorder=3,
                )
            if n_st >= 2:
                from numpy.polynomial import polynomial as P

                c = P.polyfit(e, n, 1)
                xs = np.linspace(e.min(), e.max(), 100)
                ax.plot(
                    xs,
                    P.polyval(xs, c),
                    color=color_fn(0, 1),
                    lw=0.9,
                    ls=ls_line,
                    alpha=0.4,
                    zorder=1,
                )
            return e, n

        _draw(before_data, "^", lambda i, n: "#6c7086", 0.7, "--")
        _draw(
            after_data,
            "o",
            lambda i, n: plt.cm.turbo(i / max(n - 1, 1)),
            0.9,
            "-",
        )

        ax.legend(
            [
                Line2D(
                    [0],
                    [0],
                    marker="^",
                    color="w",
                    markerfacecolor="#6c7086",
                    ms=8,
                ),
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor="#89b4fa",
                    ms=8,
                ),
            ],
            ["Before", "After"],
            fontsize=7,
            loc="best",
            facecolor=s["bg"],
            edgecolor=s["spine"],
            labelcolor=s["fg"],
        )
        ax.set_xlabel("Easting (m)", fontsize=8, color=s["fg"])
        ax.set_ylabel("Northing (m)", fontsize=8, color=s["fg"])
        ax.set_title(
            "Station positions — overlay", fontsize=9, color=s["title"], pad=5
        )
        ax.tick_params(colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.grid(True, color=s["grid"], alpha=0.25, ls="--", lw=0.5)
        ax.set_aspect("equal", adjustable="datalim")

    def plot_station_elevation(self, data, ax, title: str = "") -> None:
        """Elevation profile: chainage along profile vs elevation, coloured + contoured."""
        from pycsamt.gis.coord_correction import (
            _pca_azimuth,
            _to_utm_arrays,
        )

        s = _DARK if self.dark else _LIGHT
        _init_ax(ax, s)

        df = self._df_from(data)
        if df is None or df.empty:
            _empty(ax, "No data", s)
            return
        valid = df.dropna(subset=["lat", "lon"]).copy()
        if valid.empty:
            _empty(ax, "No valid coordinates", s)
            return

        try:
            east, north, _ = _to_utm_arrays(
                valid["lat"].values, valid["lon"].values
            )
            az_rad = np.radians(_pca_azimuth(east, north))
            ue, un = np.sin(az_rad), np.cos(az_rad)
            oe, on_ = east.mean(), north.mean()
            chain = (east - oe) * ue + (north - on_) * un
            elev = valid["elev"].values.astype(float)
            stations = (
                valid["station"].values
                if "station" in valid.columns
                else [str(i + 1) for i in range(len(chain))]
            )

            order = np.argsort(chain)
            ch_s, el_s, st_s = chain[order], elev[order], stations[order]

            import matplotlib.pyplot as plt

            from pycsamt.compat.numpy import ptp as _ptp

            cmap = plt.cm.terrain
            norm = (
                plt.Normalize(el_s.min(), el_s.max())
                if _ptp(el_s) > 0
                else plt.Normalize(el_s.min() - 1, el_s.min() + 1)
            )
            colors = [cmap(norm(e)) for e in el_s]

            # Filled area under profile
            ax.fill_between(
                ch_s / 1000.0,
                el_s,
                el_s.min() - 5,
                alpha=0.18,
                color=s["grid"],
            )
            ax.plot(ch_s / 1000.0, el_s, color=s["fg"], lw=1.2, zorder=2)
            for ci, (ch, el, st) in enumerate(zip(ch_s, el_s, st_s)):
                ax.scatter(
                    ch / 1000.0,
                    el,
                    color=colors[ci],
                    s=55,
                    edgecolors=s["spine"],
                    linewidths=0.5,
                    zorder=3,
                )
                short = st.split("_")[-1] if "_" in st else st
                ax.annotate(
                    short,
                    xy=(ch / 1000.0, el),
                    xytext=(0, 5),
                    textcoords="offset points",
                    fontsize=5.5,
                    color=s["fg"],
                    ha="center",
                    clip_on=True,
                )

            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cb = ax.get_figure().colorbar(sm, ax=ax, pad=0.02, fraction=0.03)
            cb.set_label("Elevation (m)", fontsize=7, color=s["fg"])
            cb.ax.yaxis.set_tick_params(color=s["tick"], labelcolor=s["tick"])

        except Exception as exc:
            _empty(ax, f"Elevation error: {exc}", s)
            return

        ax.set_xlabel("Profile distance (km)", fontsize=8, color=s["fg"])
        ax.set_ylabel("Elevation (m)", fontsize=8, color=s["fg"])
        ax.set_title(
            title or "Elevation profile", fontsize=9, color=s["title"], pad=5
        )
        ax.tick_params(colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.grid(True, color=s["grid"], alpha=0.25, ls="--", lw=0.5)

    # ── Elevation overlay (before + after on one axis) ───────────────────────

    def plot_station_elevation_overlay(
        self, before_data, after_data, ax
    ) -> None:
        """
        Plot both the before and after elevation profiles on a single axis.

        Both profiles share the SAME reference frame (azimuth + origin computed
        from the before dataset) so chainage values are directly comparable.

        Before → blue dashed line.   After → green solid line.
        """
        from pycsamt.gis.coord_correction import (
            _pca_azimuth,
            _to_utm_arrays,
        )

        s = _DARK if self.dark else _LIGHT
        _init_ax(ax, s)

        # ── 1. Build the shared reference frame from the before dataset ────────
        df_b = self._df_from(before_data)
        ref_az_rad, ref_oe, ref_on = None, None, None
        if df_b is not None and not df_b.empty:
            valid_b = df_b.dropna(subset=["lat", "lon"])
            if not valid_b.empty:
                try:
                    eb, nb, _ = _to_utm_arrays(
                        valid_b["lat"].values, valid_b["lon"].values
                    )
                    ref_az_rad = np.radians(_pca_azimuth(eb, nb))
                    ref_oe, ref_on = eb.mean(), nb.mean()
                except Exception:
                    pass

        if ref_az_rad is None:
            # Fallback: derive reference from after data if before is unusable
            df_a = self._df_from(after_data)
            if df_a is not None and not df_a.empty:
                valid_a = df_a.dropna(subset=["lat", "lon"])
                if not valid_a.empty:
                    try:
                        ea, na, _ = _to_utm_arrays(
                            valid_a["lat"].values, valid_a["lon"].values
                        )
                        ref_az_rad = np.radians(_pca_azimuth(ea, na))
                        ref_oe, ref_on = ea.mean(), na.mean()
                    except Exception:
                        pass
            if ref_az_rad is None:
                _empty(ax, "No coordinate data available", s)
                return

        ue, un = np.sin(ref_az_rad), np.cos(ref_az_rad)

        # ── 2. Project a dataset onto the shared axis ──────────────────────────
        def _chain_elev(data):
            df = self._df_from(data)
            if df is None or df.empty:
                return None, None, None
            valid = df.dropna(subset=["lat", "lon"]).copy()
            if valid.empty:
                return None, None, None
            try:
                east, north, _ = _to_utm_arrays(
                    valid["lat"].values, valid["lon"].values
                )
                # Use the shared origin and direction — same X-axis for both
                chain = (east - ref_oe) * ue + (north - ref_on) * un
                elev = valid["elev"].values.astype(float)
                order = np.argsort(chain)
                sts = (
                    valid["station"].values
                    if "station" in valid.columns
                    else [str(i + 1) for i in range(len(chain))]
                )
                return chain[order] / 1000.0, elev[order], sts[order]
            except Exception:
                return None, None, None

        # ── 3. Draw ────────────────────────────────────────────────────────────
        clr_before = "#5b9cf6"  # blue — Before
        clr_after = "#3fc96d"  # green — After

        ch_b, el_b, _ = _chain_elev(before_data)
        ch_a, el_a, _ = _chain_elev(after_data)

        drawn = False
        if ch_b is not None and np.isfinite(el_b).any():
            base_b = np.nanmin(el_b) - 5
            ax.fill_between(ch_b, el_b, base_b, alpha=0.10, color=clr_before)
            ax.plot(
                ch_b,
                el_b,
                color=clr_before,
                lw=1.6,
                ls="--",
                label="Before",
                zorder=2,
            )
            ax.scatter(
                ch_b,
                el_b,
                color=clr_before,
                s=28,
                edgecolors=s["spine"],
                linewidths=0.4,
                zorder=3,
            )
            drawn = True

        if ch_a is not None and np.isfinite(el_a).any():
            base_a = np.nanmin(el_a) - 5
            ax.fill_between(ch_a, el_a, base_a, alpha=0.12, color=clr_after)
            ax.plot(
                ch_a,
                el_a,
                color=clr_after,
                lw=2.0,
                ls="-",
                label="After",
                zorder=4,
            )
            ax.scatter(
                ch_a,
                el_a,
                color=clr_after,
                s=42,
                edgecolors=s["spine"],
                linewidths=0.5,
                zorder=5,
            )
            drawn = True

        if not drawn:
            _empty(ax, "No elevation data", s)
            return

        ax.legend(
            fontsize=8,
            facecolor=s["bg"],
            edgecolor=s["spine"],
            labelcolor=s["fg"],
        )
        ax.set_xlabel("Profile distance (km)", fontsize=8, color=s["fg"])
        ax.set_ylabel("Elevation (m)", fontsize=8, color=s["fg"])
        ax.set_title(
            "Elevation profile — Before / After overlay",
            fontsize=9,
            color=s["title"],
            pad=5,
        )
        ax.tick_params(colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.grid(True, color=s["grid"], alpha=0.25, ls="--", lw=0.5)

    # ── Tensor rotation rose diagram ──────────────────────────────────────────

    def plot_rotation_rose(self, before_sites, after_sites, fig) -> None:
        """
        Rose diagram showing geoelectric strike distribution before and after rotation.

        Left: before — Right: after.  Bars = count per 10° bin, coloured by |Zxy|.
        The dominant orientation of bars reveals the effective strike direction.
        """
        s = _DARK if self.dark else _LIGHT
        fig.clear()
        fig.patch.set_facecolor(s["fig_bg"])

        ax_b = fig.add_subplot(121, projection="polar")
        ax_a = fig.add_subplot(122, projection="polar")
        for ax, sites, label in (
            (ax_b, before_sites, "Before"),
            (ax_a, after_sites, "After"),
        ):
            ax.set_facecolor(s["bg"])
            ax.set_title(label, fontsize=9, color=s["title"], pad=10)
            ax.tick_params(colors=s["tick"], labelsize=6)
            ax.spines["polar"].set_edgecolor(s["spine"])
            if sites is None:
                continue
            try:
                self._draw_strike_rose(sites, ax, s)
            except Exception:
                pass

    def _draw_strike_rose(self, sites, ax, s) -> None:
        """Compute per-station |Zxy| amplitude weighted by azimuth and plot rose."""
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
        )

        angles, weights = [], []
        for ed in _iter_items(sites):
            _, z, freqs = _get_z_block(ed)
            if z is None or freqs is None or freqs.size == 0:
                continue
            # Apparent resistivity for xy and yx — ratio tells us strike proximity
            # field-unit ρ_a = 0.2·|Z|²/f (used here only as a relative weight)
            rho_xy = 0.2 * np.abs(z[:, 0, 1]) ** 2 / np.maximum(freqs, 1e-30)
            0.2 * np.abs(z[:, 1, 0]) ** 2 / np.maximum(freqs, 1e-30)
            # Angle: direction of maximum off-diagonal element
            # Sweep through 0–180° and find angle that maximises |Zxy|
            from pycsamt.emtools.strike import _rotate_tensor

            best_ang = 0.0
            best_val = -np.inf
            for deg in np.arange(-90, 91, 5):
                zr = _rotate_tensor(z, deg)
                val = float(np.nanmean(np.abs(zr[:, 0, 1])))
                if val > best_val:
                    best_val = val
                    best_ang = deg
            angles.append(best_ang % 180)
            weights.append(float(np.nanmean(rho_xy)))

        if not angles:
            ax.text(
                0.5,
                0.5,
                "No Z data",
                transform=ax.transAxes,
                ha="center",
                color=s["muted"],
                fontsize=9,
            )
            return

        angles_rad = np.radians(angles)
        bins = np.linspace(0, np.pi, 19)  # 10° bins over [0, 180°]
        counts, _ = np.histogram(angles_rad, bins=bins)
        width = bins[1] - bins[0]
        theta = (bins[:-1] + bins[1:]) / 2  # bin centres

        import matplotlib.pyplot as plt

        cmap = plt.cm.plasma
        max_c = counts.max() if counts.max() > 0 else 1
        ax.bar(
            theta,
            counts,
            width=width,
            bottom=0.0,
            color=[cmap(c / max_c) for c in counts],
            edgecolor=s["spine"],
            linewidth=0.4,
            alpha=0.85,
        )
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)  # clockwise = geographic convention
        ax.set_xticks(np.radians([0, 45, 90, 135]))
        ax.set_xticklabels(
            ["N", "NE", "E", "SE"], fontsize=6, color=s["tick"]
        )
        ax.yaxis.set_tick_params(labelsize=5, colors=s["tick"])
        ax.grid(True, color=s["grid"], alpha=0.3, lw=0.5)

    def plot_displacement_diff(self, before, after, ax) -> None:
        """Horizontal bar chart of displacement per station (metres moved)."""
        from pycsamt.gis.coord_correction import (
            _get_coords_df,
            _to_utm_arrays,
        )

        s = _DARK if self.dark else _LIGHT
        ax.set_facecolor(s["bg"])
        fig = ax.get_figure()
        if fig:
            fig.patch.set_facecolor(s["fig_bg"])
        ax.clear()
        ax.set_facecolor(s["bg"])

        if before is None or after is None:
            _empty(ax, "Need before and after to show displacement", s)
            return
        try:
            db = _get_coords_df(before)
            da = _get_coords_df(after)
            merged = db.merge(da, on="station", suffixes=("_b", "_a"))
            if merged.empty:
                _empty(ax, "No matching stations", s)
                return
            merged = merged.dropna(
                subset=["lat_b", "lon_b", "lat_a", "lon_a"]
            )
            if merged.empty:
                _empty(ax, "No valid coordinates to compare", s)
                return

            eb, nb, _ = _to_utm_arrays(
                merged["lat_b"].values, merged["lon_b"].values
            )
            ea, na, _ = _to_utm_arrays(
                merged["lat_a"].values, merged["lon_a"].values
            )
            horiz_disp = np.sqrt((ea - eb) ** 2 + (na - nb) ** 2)
            elev_disp = merged["elev_a"].values - merged["elev_b"].values

            y = np.arange(len(merged))
            ax.barh(
                y - 0.2,
                horiz_disp,
                height=0.35,
                color="#89b4fa",
                alpha=0.8,
                label="Horiz. (m)",
            )
            ax.barh(
                y + 0.2,
                np.abs(elev_disp),
                height=0.35,
                color="#a6e3a1",
                alpha=0.8,
                label="Elev. |Δ| (m)",
            )
            ax.set_yticks(y)
            ax.set_yticklabels(merged["station"].values, fontsize=6)
            ax.legend(
                fontsize=7,
                facecolor=s["bg"],
                edgecolor=s["spine"],
                labelcolor=s["fg"],
            )
        except Exception as exc:
            _empty(ax, f"Displacement error: {exc}", s)
            return

        ax.set_xlabel("Displacement (m)", fontsize=8, color=s["fg"])
        ax.set_title(
            "Position change per station", fontsize=9, color=s["title"], pad=5
        )
        ax.tick_params(colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.grid(True, axis="x", color=s["grid"], alpha=0.25, ls="--", lw=0.5)


import pandas as pd


def _init_ax(ax, s: dict) -> None:
    """Reset axes with theme background."""
    ax.set_facecolor(s["bg"])
    fig = ax.get_figure()
    if fig:
        fig.patch.set_facecolor(s["fig_bg"])
    ax.clear()
    ax.set_facecolor(s["bg"])


# ── Theme dicts ────────────────────────────────────────────────────────────────

_DARK = dict(
    bg="#1e1e2e",
    fig_bg="#181825",
    fg="#cdd6f4",
    title="#cdd6f4",
    tick="#a6adc8",
    spine="#45475a",
    grid="#313244",
    muted="#585b70",
)
_LIGHT = dict(
    bg="#eff1f5",
    fig_bg="#e6e9ef",
    fg="#4c4f69",
    title="#4c4f69",
    tick="#6c6f85",
    spine="#bcc0cc",
    grid="#ccd0da",
    muted="#9ca0b0",
)


def _empty(ax, msg: str, s: dict) -> None:
    ax.text(
        0.5,
        0.5,
        msg,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=10,
        color=s["muted"],
    )
