# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents._corrections
============================

**Single source of truth** that maps the desktop *correction section*
(:data:`pycsamt.app.desktop.controllers.correction_controller.CATALOGUE`)
onto Agent-Master workflow ids.

Each entry turns one catalogue method into a chat task with full parameter
control: the chat routes a natural-language request to the workflow id, the
parameter modal renders the method's :class:`ParamSpec` list, and a single
parameterised ``correction`` :class:`~pycsamt.agents.tooling.ToolAgent` kind
drives :class:`CorrectionController` so the algorithm is never re-implemented.
The corrected ``Sites`` flow back through ``AgentResult.data["corrected_sites"]``
into the existing post-processing modal (apply-to-session / export), the same
pathway used by the static-shift / denoise / freq_editor workflows.

This module is import-light: it holds only static metadata. The actual
catalogue (and therefore numpy / pycsamt internals) is imported lazily inside
the helpers, so importing :data:`CORRECTION_METHODS` from the routing core
(``_workflows``, ``router``) costs nothing at start-up.

Registry schema
---------------
``CORRECTION_METHODS[wf_id]`` →

``fn``            catalogue function name (``CATALOGUE[category][label]["fn"]``)
``category``      catalogue category key
``label``         catalogue method label
``title``         param-modal title
``icon``          bootstrap-icon class for the modal header
``running_label`` friendly verb phrase for the "executing…" header
``description``   one-line description for ``WORKFLOW_DESCRIPTIONS``
``keywords``      ordered natural-language phrases (specific → generic)

Waves A–E extend :data:`CORRECTION_METHODS` only; the wiring in ``tooling``,
``chat``, ``params``, ``_workflows`` and ``router`` is generic.
"""
from __future__ import annotations

from typing import Any

__all__ = [
    "CORRECTION_METHODS",
    "param_specs",
    "method_desc",
    "coerce_kwargs",
    "fn_for",
]


# ── registry ────────────────────────────────────────────────────────────────────
# Wave 0 ships one representative method end-to-end (AMA static shift). Later
# waves append their catalogue methods here; nothing else needs editing.
CORRECTION_METHODS: dict[str, dict] = {
    # ── Static Shift (Wave A) ────────────────────────────────────────────────
    "corr_ss_ama": {
        "fn":            "correct_ss_ama",
        "category":      "Static Shift",
        "label":         "AMA (spatial average)",
        "title":         "Static Shift — AMA",
        "icon":          "bi-arrows-expand-vertical",
        "running_label": "AMA static-shift correction",
        "description":   "Static-shift correction (AMA spatial average)",
        "keywords": [
            "ama static shift",
            "static shift ama",
            "ama correction",
            "array moving average",
            "spatial average static shift",
        ],
    },
    "corr_ss_loess": {
        "fn":            "_correct_ss_loess",
        "category":      "Static Shift",
        "label":         "LOESS",
        "title":         "Static Shift — LOESS",
        "icon":          "bi-arrows-expand-vertical",
        "running_label": "LOESS static-shift correction",
        "description":   "Static-shift correction (LOESS local regression)",
        "keywords": [
            "loess static shift",
            "static shift loess",
            "loess correction",
            "local regression static shift",
        ],
    },
    "corr_ss_bilateral": {
        "fn":            "_correct_ss_bilateral",
        "category":      "Static Shift",
        "label":         "Bilateral filter",
        "title":         "Static Shift — Bilateral filter",
        "icon":          "bi-arrows-expand-vertical",
        "running_label": "bilateral static-shift correction",
        "description":   "Static-shift correction (bilateral spatial filter)",
        "keywords": [
            "bilateral static shift",
            "static shift bilateral",
            "bilateral filter static",
            "bilateral correction",
        ],
    },
    "corr_ss_refmedian": {
        "fn":            "_correct_ss_refmedian",
        "category":      "Static Shift",
        "label":         "Reference median",
        "title":         "Static Shift — Reference median",
        "icon":          "bi-arrows-expand-vertical",
        "running_label": "reference-median static-shift correction",
        "description":   "Static-shift correction (global reference median)",
        "keywords": [
            "reference median static",
            "ref median static shift",
            "reference-median static shift",
            "median reference static shift",
        ],
    },
    "corr_ss_emap": {
        "fn":            "correct_static_shift",
        "category":      "Static Shift",
        "label":         "Hanning EMAP (Torres-Verdín)",
        "title":         "Static Shift — Hanning EMAP",
        "icon":          "bi-arrows-expand-vertical",
        "running_label": "Hanning EMAP static-shift correction",
        "description":   "Static-shift correction (Hanning EMAP spatial filter)",
        "keywords": [
            "hanning emap",
            "emap static shift",
            "torres-verdin",
            "torres verdin",
            "hanning static shift",
            "emap filter",
        ],
    },
    # ── Noise Removal (Wave B) ───────────────────────────────────────────────
    "corr_notch": {
        "fn":            "notch_powerline",
        "category":      "Noise Removal",
        "label":         "Notch filter (power line)",
        "title":         "Noise Removal — Notch filter",
        "icon":          "bi-soundwave",
        "running_label": "power-line notch filtering",
        "description":   "Power-line harmonic notch filter (50/60 Hz)",
        "keywords": [
            "notch filter",
            "power line filter",
            "powerline filter",
            "power-line notch",
            "remove 50 hz",
            "remove 60 hz",
            "mains noise",
            "powerline harmonics",
        ],
    },
    "corr_smooth_logfreq": {
        "fn":            "smooth_logfreq",
        "category":      "Noise Removal",
        "label":         "Log-freq smoothing",
        "title":         "Noise Removal — Log-freq smoothing",
        "icon":          "bi-soundwave",
        "running_label": "log-frequency smoothing",
        "description":   "Log-frequency kernel smoothing of the tensor",
        "keywords": [
            "log frequency smoothing",
            "log-freq smoothing",
            "logfreq smooth",
            "smooth in log frequency",
            "frequency domain smoothing",
        ],
    },
    "corr_smooth_rho_phase": {
        "fn":            "smooth_rho_phase",
        "category":      "Noise Removal",
        "label":         "ρ/φ trend smoothing",
        "title":         "Noise Removal — ρ/φ trend smoothing",
        "icon":          "bi-soundwave",
        "running_label": "ρ/φ trend smoothing",
        "description":   "Polynomial ρ/φ trend smoothing (Z rebuilt consistently)",
        "keywords": [
            "rho phi trend smoothing",
            "rho/phi trend smoothing",
            "smooth rho phase",
            "smooth rho/phi trend",
            "rho phase trend smoothing",
            "polynomial rho phi smoothing",
        ],
    },
    # ── Tensor Rotation (Wave C) ─────────────────────────────────────────────
    "corr_rotate_angle": {
        "fn":            "_wrap_rotate",
        "category":      "Tensor Rotation",
        "label":         "Rotate by fixed angle",
        "title":         "Tensor Rotation — Fixed angle",
        "icon":          "bi-arrow-clockwise",
        "running_label": "fixed-angle tensor rotation",
        "description":   "Rotate Z/tipper by a fixed angle θ",
        "keywords": [
            "rotate by angle",
            "rotate by",
            "fixed angle rotation",
            "rotate the tensor by",
            "rotate impedance by",
        ],
    },
    "corr_rotate_strike": {
        "fn":            "_wrap_rotate_to_strike",
        "category":      "Tensor Rotation",
        "label":         "Rotate to geoelectric strike",
        "title":         "Tensor Rotation — Geoelectric strike",
        "icon":          "bi-arrow-clockwise",
        "running_label": "rotation to geoelectric strike",
        "description":   "Rotate each Z to its estimated geoelectric strike",
        "keywords": [
            "rotate to strike",
            "rotate to geoelectric strike",
            "align to strike",
            "rotate to principal axes",
            "rotate into strike",
        ],
    },
    "corr_rotate_pt_strike": {
        "fn":            "_wrap_rotate_pt_strike",
        "category":      "Tensor Rotation",
        "label":         "Rotate to phase-tensor strike",
        "title":         "Tensor Rotation — Phase-tensor strike",
        "icon":          "bi-arrow-clockwise",
        "running_label": "rotation to phase-tensor strike",
        "description":   "Rotate using the phase-tensor strike (Caldwell 2004)",
        "keywords": [
            "rotate to phase tensor strike",
            "phase tensor strike rotation",
            "rotate to pt strike",
            "rotate using phase tensor",
        ],
    },
    "corr_rotate_profile": {
        "fn":            "_wrap_rotate_to_profile",
        "category":      "Tensor Rotation",
        "label":         "Rotate to profile azimuth",
        "title":         "Tensor Rotation — Profile azimuth",
        "icon":          "bi-arrow-clockwise",
        "running_label": "rotation to profile azimuth",
        "description":   "Rotate Z so x aligns with the survey profile azimuth",
        "keywords": [
            "rotate to profile",
            "align to profile azimuth",
            "rotate to profile azimuth",
            "rotate along the profile",
            "rotate to survey line",
        ],
    },
    "corr_antisymmetrize": {
        "fn":            "antisymmetrize",
        "category":      "Tensor Rotation",
        "label":         "Antisymmetrize (2-D prep)",
        "title":         "Tensor Rotation — Antisymmetrize",
        "icon":          "bi-arrow-clockwise",
        "running_label": "tensor antisymmetrization",
        "description":   "Enforce off-diagonal antisymmetry for 2-D inversion",
        "keywords": [
            "antisymmetrize",
            "antisymmetrise",
            "anti-symmetrize",
            "2d prep tensor",
            "enforce antisymmetry",
        ],
    },
    # ── Coordinates (Wave D) ─────────────────────────────────────────────────
    "corr_coord_projection": {
        "fn":            "_coord_profile_projection",
        "category":      "Coordinates",
        "label":         "Profile projection",
        "title":         "Coordinates — Profile projection",
        "icon":          "bi-geo-alt",
        "running_label": "profile projection",
        "description":   "Project stations onto the best-fit survey line",
        "keywords": [
            "project to profile",
            "project stations onto",
            "profile projection",
            "project onto the line",
            "remove cross-profile scatter",
        ],
    },
    "corr_coord_spacing": {
        "fn":            "_coord_spacing_regularize",
        "category":      "Coordinates",
        "label":         "Spacing regularization",
        "title":         "Coordinates — Spacing regularization",
        "icon":          "bi-geo-alt",
        "running_label": "station spacing regularization",
        "description":   "Redistribute stations at uniform spacing",
        "keywords": [
            "regularize spacing",
            "regularise spacing",
            "uniform station spacing",
            "even out station spacing",
            "regularize station spacing",
        ],
    },
    "corr_coord_snap": {
        "fn":            "_coord_outlier_snap",
        "category":      "Coordinates",
        "label":         "Outlier snap-to-line",
        "title":         "Coordinates — Outlier snap-to-line",
        "icon":          "bi-geo-alt",
        "running_label": "outlier snap-to-line",
        "description":   "Snap off-line outlier stations to the profile",
        "keywords": [
            "snap outliers",
            "snap outlier stations",
            "snap to line",
            "snap stations to the line",
            "outlier snap",
        ],
    },
    "corr_coord_elevation": {
        "fn":            "_coord_elevation_smooth",
        "category":      "Coordinates",
        "label":         "Elevation smoothing",
        "title":         "Coordinates — Elevation smoothing",
        "icon":          "bi-geo-alt",
        "running_label": "elevation smoothing",
        "description":   "Smooth noisy GPS elevation along the profile",
        "keywords": [
            "smooth elevation",
            "elevation smoothing",
            "smooth the elevations",
            "smooth gps elevation",
            "denoise elevation",
        ],
    },
    "corr_coord_shift": {
        "fn":            "_coord_shift",
        "category":      "Coordinates",
        "label":         "Coordinate shift",
        "title":         "Coordinates — Coordinate shift",
        "icon":          "bi-geo-alt",
        "running_label": "coordinate shift",
        "description":   "Apply a uniform datum/GPS offset to all stations",
        "keywords": [
            "shift coordinates",
            "coordinate shift",
            "offset coordinates",
            "datum shift",
            "apply gps offset",
        ],
    },
    "corr_coord_interpolate": {
        "fn":            "_coord_interpolate_missing",
        "category":      "Coordinates",
        "label":         "Interpolate missing",
        "title":         "Coordinates — Interpolate missing",
        "icon":          "bi-geo-alt",
        "running_label": "missing-coordinate interpolation",
        "description":   "Fill missing/zero coordinates by interpolation",
        "keywords": [
            "interpolate missing coordinates",
            "fill missing coordinates",
            "interpolate missing positions",
            "fill gps dropout",
            "interpolate coordinates",
        ],
    },
    # ── Source Effects (Wave E) ──────────────────────────────────────────────
    "corr_near_field": {
        "fn":            "_correct_near_field",
        "category":      "Source Effects",
        "label":         "Near-field correction",
        "title":         "Source Effects — Near-field correction",
        "icon":          "bi-broadcast-pin",
        "running_label": "near-field correction",
        "description":   "Correct CSAMT impedance for near-field contamination",
        "keywords": [
            "near field correction",
            "near-field correction",
            "correct near field",
            "csamt near field",
            "remove near field",
        ],
    },
    # ── Stratagem EDI-native pipeline (Wave E) ───────────────────────────────
    "corr_strat_qc": {
        "fn":            "_strat_qc",
        "category":      "Stratagem",
        "label":         "QC Report",
        "title":         "Stratagem — QC Report",
        "icon":          "bi-shield-check",
        "running_label": "Stratagem QC report",
        "description":   "Stratagem station-level quality-control report",
        "keywords": [
            "stratagem qc",
            "stratagem quality control",
            "stratagem qc report",
        ],
    },
    "corr_strat_static_shift": {
        "fn":            "_strat_static_shift",
        "category":      "Stratagem",
        "label":         "Static Shift (AMA)",
        "title":         "Stratagem — Static Shift (AMA)",
        "icon":          "bi-arrows-expand-vertical",
        "running_label": "Stratagem AMA static-shift",
        "description":   "Stratagem-native AMA static-shift correction",
        "keywords": [
            "stratagem static shift",
            "stratagem ama",
            "stratagem static-shift",
        ],
    },
    "corr_strat_noise": {
        "fn":            "_strat_noise_removal",
        "category":      "Stratagem",
        "label":         "Noise Removal",
        "title":         "Stratagem — Noise Removal",
        "icon":          "bi-soundwave",
        "running_label": "Stratagem noise removal",
        "description":   "Stratagem notch + Hampel + smoothing noise pipeline",
        "keywords": [
            "stratagem noise removal",
            "stratagem denoise",
            "stratagem noise",
        ],
    },
    "corr_strat_freq_filter": {
        "fn":            "_strat_freq_filter",
        "category":      "Stratagem",
        "label":         "Frequency Filter",
        "title":         "Stratagem — Frequency Filter",
        "icon":          "bi-funnel",
        "running_label": "Stratagem frequency filter",
        "description":   "Stratagem band selection + incoherence masking",
        "keywords": [
            "stratagem frequency filter",
            "stratagem freq filter",
            "stratagem band filter",
        ],
    },
    "corr_strat_full": {
        "fn":            "_strat_full_pipeline",
        "category":      "Stratagem",
        "label":         "Full Pipeline",
        "title":         "Stratagem — Full Pipeline",
        "icon":          "bi-diagram-3",
        "running_label": "Stratagem full pipeline",
        "description":   "Complete Stratagem QC → filter → SS → noise pipeline",
        "keywords": [
            "stratagem full pipeline",
            "stratagem pipeline",
            "full stratagem pipeline",
            "run stratagem",
        ],
    },
}


# ── catalogue access (lazy) ──────────────────────────────────────────────────────

def _catalogue() -> dict:
    """Return the desktop correction CATALOGUE (imported lazily)."""
    from ..app.desktop.controllers.correction_controller import (
        CATALOGUE,
    )
    return CATALOGUE


def _entry(wf_id: str) -> dict:
    meta = CORRECTION_METHODS[wf_id]
    return _catalogue()[meta["category"]][meta["label"]]


def param_specs(wf_id: str) -> list:
    """The ordered ``ParamSpec`` list for a correction workflow id."""
    return list(_entry(wf_id).get("params", []))


def method_desc(wf_id: str) -> str:
    """The catalogue method description (used as the param-modal blurb)."""
    return str(_entry(wf_id).get("desc", "") or "")


def fn_for(wf_id: str) -> str:
    """The catalogue function name a workflow id applies."""
    return CORRECTION_METHODS[wf_id]["fn"]


# ── parameter coercion ──────────────────────────────────────────────────────────
# The parameter modal stores raw values (strings / numbers / bools) in
# STORE_INV_CONFIG; coerce them back to the type each catalogue function
# expects, falling back to the ParamSpec default on any bad value.

def _coerce_check(v: Any) -> bool:
    if isinstance(v, str):
        return v.strip().lower() in ("true", "1", "yes", "on")
    return bool(v)


_COERCE = {
    "spin":  lambda v: int(round(float(v))),
    "dspin": lambda v: float(v),
    "combo": lambda v: str(v),
    "check": _coerce_check,
}


def coerce_kwargs(wf_id: str, values: dict) -> dict:
    """Build the kwargs for a correction call from raw modal *values*.

    Only the method's own ``ParamSpec`` names are kept, so stray keys carried
    in STORE_INV_CONFIG (``fn_name``, ``corr_wf``, ``stations``, ``workflow``…)
    never leak into the catalogue function call.
    """
    out: dict[str, Any] = {}
    for ps in param_specs(wf_id):
        raw = values.get(ps.name)
        if raw is None or raw == "":
            out[ps.name] = ps.default
            continue
        try:
            out[ps.name] = _COERCE.get(ps.kind, lambda v: v)(raw)
        except (TypeError, ValueError):
            out[ps.name] = ps.default
    return out
