"""Step registry for the pyCSAMT processing pipeline.

Every processing operation exposed by :mod:`pycsamt.emtools` is described by a
:class:`StepSpec` entry in :data:`STEP_REGISTRY`.  Each entry carries:

- A short **code** (e.g. ``"NR001"``) and a human **name** (e.g.
  ``"notch_powerline"``) so users can look up steps either way.
- The **module path** and **function name** needed to import the transform
  lazily (avoiding circular imports at startup).
- An optional **override_fn** for two-phase steps (SS002/SS003) and wrapper
  steps (NR007) that cannot be described by a single function reference.
- A list of **(module, fn_name)** pairs for QC / diagnostic plot functions
  that the pipeline calls automatically after each step.
- **Default parameters** merged with any user-supplied overrides.
- A ``returns_sites`` flag: ``False`` for diagnostic-only (plot) steps that
  do not transform the Sites object.
"""

from __future__ import annotations

import importlib
from dataclasses import dataclass, field
from typing import Any, Callable

# ---------------------------------------------------------------------------
# Special-case wrappers
# ---------------------------------------------------------------------------


def _emap_confidence_fn(sites: Any, **kw: Any) -> Any:
    """Wrap confidence_gated_emap_filter → return result.sites."""
    from pycsamt.emtools.remove_noise import (
        confidence_gated_emap_filter,
    )

    return confidence_gated_emap_filter(sites, **kw).sites


def _ss_loess_fn(sites: Any, **kw: Any) -> Any:
    """Two-phase SS: estimate via LOESS, then apply."""
    from pycsamt.emtools.ss import (
        apply_ss_factors,
        estimate_ss_loess,
    )

    factors = estimate_ss_loess(sites, **kw)
    return apply_ss_factors(sites, factors)


def _ss_refmedian_fn(sites: Any, **kw: Any) -> Any:
    """Two-phase SS: estimate via reference-median, then apply."""
    from pycsamt.emtools.ss import (
        apply_ss_factors,
        estimate_ss_refmedian,
    )

    factors = estimate_ss_refmedian(sites, **kw)
    return apply_ss_factors(sites, factors)


def _ss_bilateral_fn(sites: Any, **kw: Any) -> Any:
    """Two-phase SS: estimate via bilateral filter, then apply."""
    from pycsamt.emtools.ss import (
        apply_ss_factors,
        estimate_ss_bilateral,
    )

    factors = estimate_ss_bilateral(sites, **kw)
    return apply_ss_factors(sites, factors)


def _qc_snapshot_fn(sites: Any, **kw: Any) -> Any:
    """QC snapshot — diagnostic only, returns sites unchanged."""
    return sites


# ---------------------------------------------------------------------------
# StepSpec
# ---------------------------------------------------------------------------


@dataclass
class StepSpec:
    """Immutable descriptor for one pipeline step.

    Parameters
    ----------
    code:
        Short uppercase identifier, e.g. ``"NR001"``.
    name:
        Snake-case name, e.g. ``"notch_powerline"``.
    label:
        Human-readable label shown in pipeline ``__repr__``.
    category:
        Logical group (``"frequency"``, ``"noise_removal"``, ``"static_shift"``,
        ``"tensor"``, ``"dimensionality"``, ``"skew"``, ``"source_effects"``,
        ``"qc"``).
    defaults:
        Keyword arguments passed to the transform function when not overridden
        by the user.
    returns_sites:
        ``True``  – step transforms the Sites object (default).
        ``False`` – diagnostic-only step; the Sites pass through unchanged.
    mod:
        Dotted module path for the primary transform function.  ``None`` when
        *override_fn* is provided.
    fn_name:
        Name of the transform function inside *mod*.  ``None`` when
        *override_fn* is provided.
    qc_defs:
        List of ``(module_path, function_name)`` pairs.  The pipeline calls
        these after each step to generate QC figures.
    override_fn:
        Direct callable override.  When set, *mod* / *fn_name* are ignored.
    """

    code: str
    name: str
    label: str
    category: str
    defaults: dict = field(default_factory=dict)
    returns_sites: bool = True
    mod: str | None = None
    fn_name: str | None = None
    qc_defs: list[tuple[str, str]] = field(default_factory=list)
    override_fn: Callable | None = field(default=None, repr=False)

    # ------------------------------------------------------------------
    # Resolution helpers
    # ------------------------------------------------------------------

    def get_fn(self) -> Callable:
        """Return (lazily imported) transform function."""
        if self.override_fn is not None:
            return self.override_fn
        if self.mod is None or self.fn_name is None:
            raise RuntimeError(
                f"StepSpec {self.code!r} has neither override_fn "
                "nor (mod, fn_name)."
            )
        m = importlib.import_module(self.mod)
        return getattr(m, self.fn_name)

    def get_qc_fns(self) -> list[tuple[str, Callable]]:
        """Return list of ``(fn_name, callable)`` QC plot functions."""
        out: list[tuple[str, Callable]] = []
        for mod_path, fn_name in self.qc_defs:
            try:
                m = importlib.import_module(mod_path)
                fn = getattr(m, fn_name, None)
                if fn is not None:
                    out.append((fn_name, fn))
            except ImportError:
                pass
        return out


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_EMTOOLS = "pycsamt.emtools"

STEP_REGISTRY: dict[str, StepSpec] = {
    spec.code: spec
    for spec in [
        # ── Frequency management ────────────────────────────────────────
        StepSpec(
            code="FREQ001",
            name="select_band",
            label="Frequency Band Select",
            category="frequency",
            mod=f"{_EMTOOLS}.frequency",
            fn_name="select_band",
            defaults={"band_hz": (1e-3, 1e4)},
            qc_defs=[
                (f"{_EMTOOLS}.frequency", "plot_band_microstrips"),
                (f"{_EMTOOLS}.frequency", "plot_coverage_quality_heatmap"),
            ],
        ),
        StepSpec(
            code="FREQ002",
            name="drop_duplicates",
            label="Drop Duplicate Frequencies",
            category="frequency",
            mod=f"{_EMTOOLS}.frequency",
            fn_name="drop_duplicates",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.frequency", "plot_coverage_quality_heatmap"),
            ],
        ),
        StepSpec(
            code="FREQ003",
            name="drop_low_confidence",
            label="Drop Low-Confidence Frequencies",
            category="frequency",
            mod=f"{_EMTOOLS}.frequency",
            fn_name="drop_low_confidence_frequencies",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.frequency", "plot_frequency_edit_decisions"),
                (f"{_EMTOOLS}.frequency", "plot_frequency_edit_summary"),
            ],
        ),
        StepSpec(
            code="FREQ004",
            name="align_grid",
            label="Frequency Grid Alignment",
            category="frequency",
            mod=f"{_EMTOOLS}.frequency",
            fn_name="align_grid",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.frequency", "plot_coverage_quality_heatmap"),
            ],
        ),
        StepSpec(
            code="FREQ005",
            name="regrid_logspace",
            label="Log-Space Frequency Regrid",
            category="frequency",
            mod=f"{_EMTOOLS}.frequency",
            fn_name="regrid_logspace",
            defaults={"n_per_decade": 6},
            qc_defs=[
                (f"{_EMTOOLS}.frequency", "plot_apparent_depth_psection"),
            ],
        ),
        StepSpec(
            code="FREQ006",
            name="decimate",
            label="Frequency Decimation",
            category="frequency",
            mod=f"{_EMTOOLS}.frequency",
            fn_name="decimate_step",
            defaults={"step": 2},
            qc_defs=[
                (f"{_EMTOOLS}.frequency", "plot_coverage_quality_heatmap"),
            ],
        ),
        StepSpec(
            code="FREQ007",
            name="smooth_freq",
            label="Frequency Moving Average Smooth",
            category="frequency",
            mod=f"{_EMTOOLS}.frequency",
            fn_name="smooth_mavg",
            defaults={"window": 3},
            qc_defs=[
                (f"{_EMTOOLS}.frequency", "plot_coverage_quality_heatmap"),
            ],
        ),
        StepSpec(
            code="FREQ008",
            name="mask_low_confidence",
            label="Mask Low-Confidence Frequencies (NaN)",
            category="frequency",
            mod=f"{_EMTOOLS}.frequency",
            fn_name="mask_low_confidence_frequencies",
            defaults={"method": "composite", "threshold": 0.5},
            qc_defs=[
                (f"{_EMTOOLS}.frequency", "plot_frequency_edit_summary"),
                (f"{_EMTOOLS}.qc", "plot_coverage_psection"),
            ],
        ),
        StepSpec(
            code="FREQ009",
            name="recover_low_confidence",
            label="Recover Low-Confidence Frequencies (Interpolate)",
            category="frequency",
            mod=f"{_EMTOOLS}.frequency",
            fn_name="recover_low_confidence_frequencies",
            defaults={
                "method": "composite",
                "ci_hi": 0.9,
                "ci_lo": 0.5,
                "interpolation": "linear",
            },
            qc_defs=[
                (f"{_EMTOOLS}.frequency", "plot_frequency_edit_summary"),
                (f"{_EMTOOLS}.frequency", "plot_coverage_quality_heatmap"),
            ],
        ),
        # ── Noise removal ───────────────────────────────────────────────
        StepSpec(
            code="NR001",
            name="notch_powerline",
            label="Power-line Harmonic Notch",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="notch_powerline",
            defaults={"mains_hz": 50, "n_harm": 30, "tol_hz": 0.08},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "nr_qc_harmonic_waterfall"),
                (f"{_EMTOOLS}.remove_noise", "nr_qc_snr_gain_profile"),
            ],
        ),
        StepSpec(
            code="NR002",
            name="smooth_logfreq",
            label="Log-Frequency Smooth",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="smooth_logfreq",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "nr_qc_station_offdiag_curves"),
            ],
        ),
        StepSpec(
            code="NR003",
            name="shrink_group_trend",
            label="Shrink Outliers to Group Trend",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="shrink_to_group_trend",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "nr_qc_delta_offdiag_psection"),
            ],
        ),
        StepSpec(
            code="NR004",
            name="hampel_filter",
            label="Hampel Outlier Filter",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="hampel_filter_freq",
            defaults={"win": 3, "nsig": 3.0},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "nr_qc_snr_gain_profile"),
            ],
        ),
        StepSpec(
            code="NR005",
            name="spatial_median",
            label="Spatial Median Filter",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="spatial_median_filter",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "nr_qc_delta_offdiag_psection"),
            ],
        ),
        StepSpec(
            code="NR006",
            name="emap_filter",
            label="EMAP Spatial Filter",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="apply_emap_filter",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "plot_emap_filter_profile"),
                (f"{_EMTOOLS}.remove_noise", "plot_emap_filter_psection"),
            ],
        ),
        StepSpec(
            code="NR007",
            name="emap_confidence",
            label="EMAP with Confidence Gating",
            category="noise_removal",
            override_fn=_emap_confidence_fn,
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "plot_emap_filter_profile"),
                (f"{_EMTOOLS}.remove_noise", "plot_emap_filter_psection"),
            ],
        ),
        StepSpec(
            code="NR008",
            name="rpca_offdiag",
            label="RPCA Off-Diagonal Denoise",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="rpca_offdiag_denoise",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "nr_qc_delta_offdiag_psection"),
            ],
        ),
        StepSpec(
            code="NR009",
            name="enforce_offdiag",
            label="Enforce Off-Diagonal Consistency",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="enforce_offdiag_consistency",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "nr_qc_delta_offdiag_psection"),
                (f"{_EMTOOLS}.remove_noise", "nr_qc_station_offdiag_curves"),
            ],
        ),
        StepSpec(
            code="NR010",
            name="mask_incoherent",
            label="Mask Incoherent Frequencies",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="mask_incoherent_freqs",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.qc", "plot_qc_quicklook"),
            ],
        ),
        StepSpec(
            code="NR011",
            name="fixed_length_mavg",
            label="Fixed-Length Moving Average Smooth",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="fixed_length_moving_average",
            defaults={"window": 5, "component": "all"},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "nr_qc_station_offdiag_curves"),
            ],
        ),
        StepSpec(
            code="NR012",
            name="trimmed_mavg",
            label="Trimmed Moving Average Smooth (Robust)",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="trimmed_moving_average",
            defaults={"window": 5, "component": "all"},
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "nr_qc_station_offdiag_curves"),
                (f"{_EMTOOLS}.remove_noise", "nr_qc_delta_offdiag_psection"),
            ],
        ),
        StepSpec(
            code="NR013",
            name="correct_static_shift_spatial",
            label="Spatial Window Static Shift Correction",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="correct_static_shift",
            defaults={"window_m": 1500.0, "spacing_m": 200.0, "comp": "det"},
            qc_defs=[
                (f"{_EMTOOLS}.ss", "ss_qc_psection"),
                (f"{_EMTOOLS}.ss", "ss_qc_profile"),
            ],
        ),
        StepSpec(
            code="NR014",
            name="denoise_pipeline",
            label="All-in-One Denoising Pipeline",
            category="noise_removal",
            mod=f"{_EMTOOLS}.remove_noise",
            fn_name="remove_noise_pipeline",
            defaults={
                "mains_hz": 50.0,
                "n_harm": 30,
                "tol_hz": 0.08,
                "smooth_win": 5,
                "gate_snr": 2.5,
            },
            qc_defs=[
                (f"{_EMTOOLS}.remove_noise", "nr_qc_harmonic_waterfall"),
                (f"{_EMTOOLS}.remove_noise", "nr_qc_snr_gain_profile"),
                (f"{_EMTOOLS}.remove_noise", "nr_qc_delta_offdiag_psection"),
            ],
        ),
        # ── Static shift correction ──────────────────────────────────────
        StepSpec(
            code="SS001",
            name="correct_ss_ama",
            label="Static Shift (AMA)",
            category="static_shift",
            mod=f"{_EMTOOLS}.ss",
            fn_name="correct_ss_ama",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.ss", "plot_ss_delta_psection"),
                (f"{_EMTOOLS}.ss", "plot_ss_summary"),
            ],
        ),
        StepSpec(
            code="SS002",
            name="correct_ss_loess",
            label="Static Shift (LOESS)",
            category="static_shift",
            override_fn=_ss_loess_fn,
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.ss", "plot_ss_delta_psection"),
                (f"{_EMTOOLS}.ss", "plot_ss_comparison_psection"),
            ],
        ),
        StepSpec(
            code="SS003",
            name="correct_ss_refmedian",
            label="Static Shift (Reference Median)",
            category="static_shift",
            override_fn=_ss_refmedian_fn,
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.ss", "plot_ss_delta_psection"),
                (f"{_EMTOOLS}.ss", "plot_ss_comparison_psection"),
            ],
        ),
        StepSpec(
            code="SS004",
            name="correct_ss_bilateral",
            label="Static Shift (Bilateral Filter)",
            category="static_shift",
            override_fn=_ss_bilateral_fn,
            defaults={"half_window": 4, "max_skew": 6.0, "summary": "median"},
            qc_defs=[
                (f"{_EMTOOLS}.ss", "plot_ss_delta_psection"),
                (f"{_EMTOOLS}.ss", "plot_ss_comparison_psection"),
                (f"{_EMTOOLS}.ss", "plot_ss_summary"),
            ],
        ),
        # ── Tensor analysis ──────────────────────────────────────────────
        StepSpec(
            code="TZ001",
            name="rotate_strike",
            label="Strike Rotation",
            category="tensor",
            mod=f"{_EMTOOLS}.tensor",
            fn_name="rotate_to_strike",
            defaults={"method": "swift"},
            qc_defs=[
                (f"{_EMTOOLS}.strike", "plot_strike_rose"),
                (f"{_EMTOOLS}.tensor", "plot_phase_tensor_psection"),
            ],
        ),
        StepSpec(
            code="TZ002",
            name="antisymmetrize",
            label="Antisymmetrise Impedance Tensor",
            category="tensor",
            mod=f"{_EMTOOLS}.tensor",
            fn_name="antisymmetrize",
            defaults={"how": "rms"},
            qc_defs=[
                (f"{_EMTOOLS}.impedance", "plot_offdiag_antisym_residual"),
            ],
        ),
        StepSpec(
            code="TZ003",
            name="sigma_clip",
            label="Sigma-Clip Outlier Frequencies",
            category="tensor",
            mod=f"{_EMTOOLS}.tensor",
            fn_name="sigma_clip_z",
            defaults={"sigma": 3},
            qc_defs=[
                (f"{_EMTOOLS}.qc", "plot_qc_quicklook"),
            ],
        ),
        StepSpec(
            code="TZ004",
            name="balance_offdiag",
            label="Balance Off-Diagonal Components",
            category="tensor",
            mod=f"{_EMTOOLS}.tensor",
            fn_name="balance_offdiag",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.impedance", "plot_offdiag_antisym_residual"),
            ],
        ),
        StepSpec(
            code="TZ005",
            name="rotate_fixed",
            label="Fixed-Angle Rotation",
            category="tensor",
            mod=f"{_EMTOOLS}.tensor",
            fn_name="rotate",
            defaults={"angle": 0.0},
            qc_defs=[
                (f"{_EMTOOLS}.tensor", "plot_phase_tensor_psection"),
                (f"{_EMTOOLS}.strike", "plot_strike_rose"),
            ],
        ),
        StepSpec(
            code="TZ006",
            name="invert_tensor",
            label="Impedance Tensor Inversion (Z → Z⁻¹)",
            category="tensor",
            mod=f"{_EMTOOLS}.tensor",
            fn_name="invert",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.impedance", "plot_determinant_track"),
                (f"{_EMTOOLS}.impedance", "plot_phasor_wheel"),
            ],
        ),
        StepSpec(
            code="TZ007",
            name="orient_from_sensors",
            label="Sensor Orientation Correction",
            category="tensor",
            mod=f"{_EMTOOLS}.tensor",
            fn_name="orient_from_sensors",
            defaults={
                "ex": 0.0,
                "ey": 0.0,
                "bx": 0.0,
                "by": 0.0,
                "degrees": True,
            },
            qc_defs=[
                (f"{_EMTOOLS}.tensor", "plot_phase_tensor_psection"),
            ],
        ),
        # ── Dimensionality ───────────────────────────────────────────────
        StepSpec(
            code="DIM001",
            name="classify_dim",
            label="Dimensionality Classification",
            category="dimensionality",
            returns_sites=False,
            mod=f"{_EMTOOLS}.dimensionality",
            fn_name="classify_dimensionality",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.dimensionality", "plot_dim_map"),
                (f"{_EMTOOLS}.dimensionality", "plot_dim_confidence_grid"),
            ],
        ),
        StepSpec(
            code="DIM002",
            name="mask_by_dim",
            label="Mask by Dimensionality Class",
            category="dimensionality",
            mod=f"{_EMTOOLS}.dimensionality",
            fn_name="mask_by_dimensionality",
            defaults={"dim_type": "2d"},
            qc_defs=[
                (f"{_EMTOOLS}.tensor", "plot_dimensionality_psection"),
            ],
        ),
        StepSpec(
            code="DIM003",
            name="project_2d",
            label="Project to 2-D (Remove 3-D Components)",
            category="dimensionality",
            mod=f"{_EMTOOLS}.dimensionality",
            fn_name="project_to_2d",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.dimensionality", "plot_dim_map"),
            ],
        ),
        # ── Skewness ─────────────────────────────────────────────────────
        StepSpec(
            code="SK001",
            name="mask_by_skew",
            label="Mask by Bahr Skewness Threshold",
            category="skew",
            mod=f"{_EMTOOLS}.skew",
            fn_name="mask_by_skew",
            defaults={"threshold": 0.3},
            qc_defs=[
                (f"{_EMTOOLS}.skew", "plot_skew_traffic_psection"),
            ],
        ),
        StepSpec(
            code="SK002",
            name="longest_low_skew",
            label="Keep Longest Low-Skew Band",
            category="skew",
            mod=f"{_EMTOOLS}.skew",
            fn_name="keep_longest_low_skew",
            defaults={"threshold": 0.3},
            qc_defs=[
                (f"{_EMTOOLS}.skew", "plot_skew_percentile_ribbon"),
            ],
        ),
        StepSpec(
            code="SK003",
            name="select_skew_band",
            label="Select Lowest-Skew Band",
            category="skew",
            mod=f"{_EMTOOLS}.skew",
            fn_name="select_low_skew_band",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.skew", "plot_skew_vote_band"),
            ],
        ),
        StepSpec(
            code="SK004",
            name="close_skew_gaps",
            label="Close Small Gaps in Low-Skew Band",
            category="skew",
            mod=f"{_EMTOOLS}.skew",
            fn_name="close_skew_gaps",
            defaults={"thresh": 3.0, "max_gap": 1},
            qc_defs=[
                (f"{_EMTOOLS}.skew", "plot_skew_traffic_psection"),
                (f"{_EMTOOLS}.skew", "plot_skew_percentile_ribbon"),
            ],
        ),
        # ── Source effects ───────────────────────────────────────────────
        StepSpec(
            code="SRC001",
            name="correct_near_field",
            label="Near-Field Correction",
            category="source_effects",
            mod=f"{_EMTOOLS}.source_effects",
            fn_name="correct_near_field",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.source_effects", "plot_overprint_section"),
            ],
        ),
        StepSpec(
            code="SRC002",
            name="normalize_response",
            label="Normalize Source Response",
            category="source_effects",
            mod=f"{_EMTOOLS}.source_effects",
            fn_name="normalize_response",
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.source_effects", "plot_normalized_response"),
            ],
        ),
        # ── QC / diagnostics ─────────────────────────────────────────────
        StepSpec(
            code="QC001",
            name="qc_snapshot",
            label="QC Quick-Look Snapshot",
            category="qc",
            returns_sites=False,
            override_fn=_qc_snapshot_fn,
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.qc", "plot_qc_quicklook"),
                (f"{_EMTOOLS}.qc", "plot_station_confidence_dashboard"),
                (f"{_EMTOOLS}.qc", "plot_coverage_psection"),
            ],
        ),
        StepSpec(
            code="QC002",
            name="field_zone_snapshot",
            label="Field Zone Classification Snapshot",
            category="qc",
            returns_sites=False,
            override_fn=_qc_snapshot_fn,
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.fieldzone", "plot_field_zones"),
            ],
        ),
        StepSpec(
            code="QC003",
            name="source_overprint_snapshot",
            label="Source Overprint Diagnostic Snapshot",
            category="qc",
            returns_sites=False,
            override_fn=_qc_snapshot_fn,
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.source_effects", "plot_overprint_section"),
            ],
        ),
        StepSpec(
            code="QC004",
            name="depth_section_snapshot",
            label="Bostick Depth Section Snapshot",
            category="qc",
            returns_sites=False,
            override_fn=_qc_snapshot_fn,
            defaults={},
            qc_defs=[
                (f"{_EMTOOLS}.csumt", "plot_depth_section"),
            ],
        ),
    ]
}

# Build a secondary index: name → StepSpec (for lookup by human name)
_NAME_INDEX: dict[str, StepSpec] = {
    spec.name: spec for spec in STEP_REGISTRY.values()
}


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def lookup_step(code_or_name: str) -> StepSpec:
    """Return the :class:`StepSpec` for *code_or_name*.

    Parameters
    ----------
    code_or_name:
        Either the uppercase code (``"NR001"``) or the snake-case name
        (``"notch_powerline"``).

    Raises
    ------
    KeyError
        When neither a code nor a name matches.
    """
    if code_or_name in STEP_REGISTRY:
        return STEP_REGISTRY[code_or_name]
    if code_or_name in _NAME_INDEX:
        return _NAME_INDEX[code_or_name]
    raise KeyError(
        f"No pipeline step found for {code_or_name!r}.  "
        f"Call list_steps() to see all available steps."
    )


def list_steps(category: str | None = None) -> list[StepSpec]:
    """Return all registered :class:`StepSpec` objects, optionally filtered.

    Parameters
    ----------
    category:
        When supplied, only steps whose :attr:`~StepSpec.category` matches
        this string are returned.
    """
    specs = list(STEP_REGISTRY.values())
    if category is not None:
        specs = [s for s in specs if s.category == category]
    return specs


def step_codes() -> list[str]:
    """Return a sorted list of all step codes."""
    return sorted(STEP_REGISTRY)


def step_names() -> list[str]:
    """Return a sorted list of all step names."""
    return sorted(_NAME_INDEX)


def categories() -> list[str]:
    """Return a sorted list of distinct step categories."""
    return sorted({s.category for s in STEP_REGISTRY.values()})


__all__ = [
    "StepSpec",
    "STEP_REGISTRY",
    "lookup_step",
    "list_steps",
    "step_codes",
    "step_names",
    "categories",
]
