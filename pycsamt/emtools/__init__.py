"""
pycsamt.emtools — public API
============================

All user-facing functions and constants from the emtools sub-modules,
organised by workflow stage.
"""
from __future__ import annotations

# ─── Rose diagram style system ────────────────────────────────────────────────
from ._rose_style import RoseStyle, resolve_rose_style

# ─── Input helpers ────────────────────────────────────────────────────────────
from ._core import ensure_sites

# ─── Survey design — Bostick depth & frequency schedules ─────────────────────
from .csumt import (
    bostick_depth_from_rho,
    vertical_resolution_pair,
    frequency_for_depth,
    frequency_schedule,
    bostick_depth,
    vertical_resolution,
    depth_coverage_table,
    plot_depth_section,
)

# ─── Source array design ──────────────────────────────────────────────────────
from .source_array import (
    wavenumber,
    sdas_element_pattern,
    array_factor,
    pas_pattern,
    beam_steer,
    steering_angles,
    sdas_directivity,
    snr_gain_db,
    plot_radiation_pattern,
)

# ─── Field zone classification ────────────────────────────────────────────────
from .fieldzone import (
    classify_field_zones,
    near_field_factor,
    plot_field_zones,
)

# ─── Source effects and near-field correction ─────────────────────────────────
from .source_effects import (
    BETA_THRESH_PCT,
    overprint_beta,
    detect_source_overprint,
    source_overprint_table,
    plot_overprint_section,
    normalize_response,
    correct_near_field,
    plot_normalized_response,
)

# ─── Diagnostics and quality control ─────────────────────────────────────────
from .diag import (
    COVERAGE_THRESH,
    coverage_score,
    rho_coverage,
    rho_error_stats,
    coverage_table,
    plot_polar_coverage,
    plot_polar_errors,
    plot_width_drift,
)
from .qc import (
    build_qc_table,
    qc_flags,
    plot_confidence_profile,
    plot_coverage_psection,
    plot_snr_hist,
    plot_qc_quicklook,
    plot_consistency_fan,
    plot_xyyx_crossover_map,
    overlay_noise_cone,
    overlay_spectral_holes,
)

# ─── Noise removal ────────────────────────────────────────────────────────────
from .remove_noise import (
    snr_table,
    notch_powerline,
    smooth_logfreq,
    shrink_to_group_trend,
    remove_noise_pipeline,
    hampel_filter_freq,
    spatial_median_filter,
    rpca_offdiag_denoise,
    enforce_offdiag_consistency,
    mask_incoherent_freqs,
    correct_static_shift,
    nr_qc_delta_offdiag_psection,
    nr_qc_snr_gain_profile,
    nr_qc_harmonic_waterfall,
    nr_qc_station_offdiag_curves,
)

# ─── Static shift correction ──────────────────────────────────────────────────
from .ss import (
    estimate_ss_ama,
    apply_ss_factors,
    correct_ss_ama,
    estimate_ss_loess,
    estimate_ss_bilateral,
    estimate_ss_refmedian,
    plot_ss_delta_psection,
    plot_ss_station_curves,
    plot_ss_delta_profile,
    ss_qc_psection,
    ss_qc_station_curves,
    ss_qc_profile,
    plot_ss_comparison_psection,
    plot_ss_1d_curves,
    plot_ss_summary,
    ss_comparison_psection,
    plot_ss_radar,
    detect_near_surface,
    plot_ns_detection,
)

# ─── Frequency grid management ────────────────────────────────────────────────
from .frequency import (
    select_band,
    drop_duplicates,
    regrid_to,
    regrid_logspace,
    decimate_step,
    smooth_mavg,
    align_grid,
    plot_coverage_quality_heatmap,
    plot_apparent_depth_psection,
    plot_band_microstrips,
)

# ─── Skewness analysis ────────────────────────────────────────────────────────
from .skew import (
    bahr_skewness,
    skew_table,
    mask_by_skew,
    keep_longest_low_skew,
    close_skew_gaps,
    select_low_skew_band,
    plot_skew_traffic_psection,
    plot_skew_percentile_ribbon,
    plot_skew_vote_band,
    plot_skewness,
)

# ─── Anisotropy analysis ──────────────────────────────────────────────────────
from .anisotropy import (
    ANISO_RATIO_THRESH,
    SWIFT_SKEW_THRESH,
    analyze_anisotropy,
    anisotropy_table,
    plot_anisotropy,
)

# ─── Dimensionality analysis ──────────────────────────────────────────────────
from .dimensionality import (
    phase_features_table,
    classify_dimensionality,
    mask_by_dimensionality,
    project_to_2d,
    learn_dim_dictionary,
    encode_dimensionality,
    mask_by_dictionary,
    plot_atom_psection,
    plot_dim_confidence_grid,
    plot_dim_occupancy_area,
    plot_dim_map,
)

# ─── Strike estimation ────────────────────────────────────────────────────────
from .strike import (
    estimate_strike_sweep,
    estimate_strike_phase_tensor,
    estimate_strike_consensus,
    rotate_to_strike,
    strike_curve_sweep,
    plot_strike_rose,
    plot_strike_rose_by_line,
    plot_strike_ribbon,
    plot_strike_profile,
    plot_strike_mapsticks,
)

# ─── Tensor operations ────────────────────────────────────────────────────────
# rotate_to_strike conflicts with strike.rotate_to_strike; alias as rotate_z_to_strike
from .tensor import (
    rotate_to_strike as rotate_z_to_strike,
    rotate,
    rotate_by_map,
    antisymmetrize,
    invert,
    orient_from_sensors,
    sigma_clip_z,
    balance_offdiag,
    build_phase_tensor_table,
    plot_phase_tensor_psection,
    plot_phase_tensor_skewmap,
    plot_theta_vs_period,
    plot_ellipticity_psection,
    plot_dimensionality_psection,
    plot_phase_tensor_rose,
    plot_phase_tensor_map,
    plot_phase_tensor_summary,
    phase_tensor_legend,
    plot_dimensionality_grid,
    plot_theta_stability_stripe,
    plot_skew_ellipt_density,
    plot_theta_rose_grid,
)

# ─── Impedance diagnostics ────────────────────────────────────────────────────
from .impedance import (
    plot_phasor_wheel,
    plot_offdiag_antisym_residual,
    plot_determinant_track,
)

# ─── Gradient imaging (Zhang et al. 2021) ────────────────────────────────────
from .gradient_imaging import (
    rho_spatial_gradient,
    rho_frequency_gradient,
    rho_joint_gradient,
    plot_gradient_section,
)

# ─── Data inspection ──────────────────────────────────────────────────────────
from .inspect import (
    sites_summary,
    list_missing_sections,
    frequency_coverage,
    plot_coverage,
    plot_rhoa_phi,
    plot_tipper_components,
    pseudosection,
)

# ─── Transfer functions — tipper ──────────────────────────────────────────────
from .tf import (
    plot_tipper_hodograms,
    plot_induction_arrows,
)

# ─── Site-level visualisation ─────────────────────────────────────────────────
from .plot import (
    plot_sites_panels,
    plot_sites_compare,
    plot_sites_fit_grid,
)

# ─── L-curve ──────────────────────────────────────────────────────────────────
from .lcurve import (
    lcurve_table,
    plot_lcurve,
)

# ─── Legacy / compatibility shims ─────────────────────────────────────────────
from .legacy import (
    check_em_kind,
    extract_z_list,
    parse_tensor,
    compute_qc,
    full_freq,
    tensor2d,
    align_tensor,
    export_edis,
    plot_confidence,
    plot_skew_1d,
    plot_skew_2d,
    plot_strike,
    plot_tensors,
    plot_station_tensors,
    wrap_phase,
    plot_skew,
)

__all__ = [
    # rose style
    "RoseStyle",
    "resolve_rose_style",
    # input helpers
    "ensure_sites",
    # survey design
    "bostick_depth_from_rho",
    "vertical_resolution_pair",
    "frequency_for_depth",
    "frequency_schedule",
    "bostick_depth",
    "vertical_resolution",
    "depth_coverage_table",
    "plot_depth_section",
    # source array design
    "wavenumber",
    "sdas_element_pattern",
    "array_factor",
    "pas_pattern",
    "beam_steer",
    "steering_angles",
    "sdas_directivity",
    "snr_gain_db",
    "plot_radiation_pattern",
    # field zone classification
    "classify_field_zones",
    "near_field_factor",
    "plot_field_zones",
    # source effects
    "BETA_THRESH_PCT",
    "overprint_beta",
    "detect_source_overprint",
    "source_overprint_table",
    "plot_overprint_section",
    "normalize_response",
    "correct_near_field",
    "plot_normalized_response",
    # diagnostics & QC
    "COVERAGE_THRESH",
    "coverage_score",
    "rho_coverage",
    "rho_error_stats",
    "coverage_table",
    "plot_polar_coverage",
    "plot_polar_errors",
    "plot_width_drift",
    "build_qc_table",
    "qc_flags",
    "plot_confidence_profile",
    "plot_coverage_psection",
    "plot_snr_hist",
    "plot_qc_quicklook",
    "plot_consistency_fan",
    "plot_xyyx_crossover_map",
    "overlay_noise_cone",
    "overlay_spectral_holes",
    # noise removal
    "snr_table",
    "notch_powerline",
    "smooth_logfreq",
    "shrink_to_group_trend",
    "remove_noise_pipeline",
    "hampel_filter_freq",
    "spatial_median_filter",
    "rpca_offdiag_denoise",
    "enforce_offdiag_consistency",
    "mask_incoherent_freqs",
    "correct_static_shift",
    "nr_qc_delta_offdiag_psection",
    "nr_qc_snr_gain_profile",
    "nr_qc_harmonic_waterfall",
    "nr_qc_station_offdiag_curves",
    # static shift
    "estimate_ss_ama",
    "apply_ss_factors",
    "correct_ss_ama",
    "estimate_ss_loess",
    "estimate_ss_bilateral",
    "estimate_ss_refmedian",
    "plot_ss_delta_psection",
    "plot_ss_station_curves",
    "plot_ss_delta_profile",
    "ss_qc_psection",
    "ss_qc_station_curves",
    "ss_qc_profile",
    "plot_ss_comparison_psection",
    "plot_ss_1d_curves",
    "plot_ss_summary",
    "ss_comparison_psection",
    "plot_ss_radar",
    "detect_near_surface",
    "plot_ns_detection",
    # frequency grid
    "select_band",
    "drop_duplicates",
    "regrid_to",
    "regrid_logspace",
    "decimate_step",
    "smooth_mavg",
    "align_grid",
    "plot_coverage_quality_heatmap",
    "plot_apparent_depth_psection",
    "plot_band_microstrips",
    # skewness
    "bahr_skewness",
    "skew_table",
    "mask_by_skew",
    "keep_longest_low_skew",
    "close_skew_gaps",
    "select_low_skew_band",
    "plot_skew_traffic_psection",
    "plot_skew_percentile_ribbon",
    "plot_skew_vote_band",
    "plot_skewness",
    # anisotropy
    "ANISO_RATIO_THRESH",
    "SWIFT_SKEW_THRESH",
    "analyze_anisotropy",
    "anisotropy_table",
    "plot_anisotropy",
    # dimensionality
    "phase_features_table",
    "classify_dimensionality",
    "mask_by_dimensionality",
    "project_to_2d",
    "learn_dim_dictionary",
    "encode_dimensionality",
    "mask_by_dictionary",
    "plot_atom_psection",
    "plot_dim_confidence_grid",
    "plot_dim_occupancy_area",
    "plot_dim_map",
    # strike
    "estimate_strike_sweep",
    "estimate_strike_phase_tensor",
    "estimate_strike_consensus",
    "rotate_to_strike",
    "strike_curve_sweep",
    "plot_strike_rose",
    "plot_strike_rose_by_line",
    "plot_strike_ribbon",
    "plot_strike_profile",
    "plot_strike_mapsticks",
    # tensor operations
    "rotate_z_to_strike",
    "rotate",
    "rotate_by_map",
    "antisymmetrize",
    "invert",
    "orient_from_sensors",
    "sigma_clip_z",
    "balance_offdiag",
    "build_phase_tensor_table",
    "plot_phase_tensor_psection",
    "plot_phase_tensor_skewmap",
    "plot_theta_vs_period",
    "plot_ellipticity_psection",
    "plot_dimensionality_psection",
    "plot_phase_tensor_rose",
    "plot_phase_tensor_map",
    "phase_tensor_legend",
    "plot_dimensionality_grid",
    "plot_theta_stability_stripe",
    "plot_skew_ellipt_density",
    "plot_theta_rose_grid",
    # impedance diagnostics
    "plot_phasor_wheel",
    "plot_offdiag_antisym_residual",
    "plot_determinant_track",
    # gradient imaging
    "rho_spatial_gradient",
    "rho_frequency_gradient",
    "rho_joint_gradient",
    "plot_gradient_section",
    # data inspection
    "sites_summary",
    "list_missing_sections",
    "frequency_coverage",
    "plot_coverage",
    "plot_rhoa_phi",
    "plot_tipper_components",
    "pseudosection",
    # transfer functions
    "plot_tipper_hodograms",
    "plot_induction_arrows",
    # site-level visualisation
    "plot_sites_panels",
    "plot_sites_compare",
    "plot_sites_fit_grid",
    # l-curve
    "lcurve_table",
    "plot_lcurve",
    # legacy
    "check_em_kind",
    "extract_z_list",
    "parse_tensor",
    "compute_qc",
    "full_freq",
    "tensor2d",
    "align_tensor",
    "export_edis",
    "plot_confidence",
    "plot_skew_1d",
    "plot_skew_2d",
    "plot_strike",
    "plot_tensors",
    "plot_station_tensors",
    "wrap_phase",
    "plot_skew",
]
