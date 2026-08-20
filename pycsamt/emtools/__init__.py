"""
pycsamt.emtools — public API
============================

All user-facing functions and constants from the emtools sub-modules,
organised by workflow stage.
"""

from __future__ import annotations

# ─── Rose diagram style system (canonical: pycsamt.api._rose_style) ──────────
from ..api._rose_style import RoseStyle, resolve_rose_style

# ─── Input helpers ────────────────────────────────────────────────────────────
from ._core import ensure_sites
from .advanced import (
    plot_apparent_anisotropy_section,
    plot_apparent_resistivity_polar,
    plot_dimensionality_depth_profile,
    plot_dimensionality_ternary,
    plot_distortion_radar,
    plot_impedance_mohr_circles,
    plot_mt_composite_section,
    plot_pt_period_clock,
    plot_rho_phase_bode,
    plot_sensitivity_depth_section,
    plot_snr_section,
    plot_strike_stability_bands,
    plot_survey_fingerprint,
    plot_tf_coherence_network,
    plot_z_invariants_section,
    plot_zt_argand,
)

# ─── AFMAG (motion-coupling physics, tilt diagnostics, QC) ───────────────────
from .afmag import (
    afmag_tilt_angles,
    coil_normal_direction,
    correct_motion_induced_noise,
    euler_rotation_matrix,
    flag_motion_susceptible_band,
    geomagnetic_field_direction,
    motion_coupling_angle,
    motion_coupling_cosine,
    motion_susceptibility_table,
    plot_afmag_correction_comparison,
    plot_afmag_tilt_polar,
    plot_afmag_tilt_profile,
    plot_afmag_tilt_psection,
    plot_motion_coupling,
    plot_motion_susceptibility_map,
    simulate_motion_induced_voltage,
)

# ─── Anisotropy analysis ──────────────────────────────────────────────────────
from .anisotropy import (
    ANISO_RATIO_THRESH,
    SWIFT_SKEW_THRESH,
    analyze_anisotropy,
    anisotropy_table,
    plot_anisotropy,
)

# ─── Survey design — Bostick depth & frequency schedules ─────────────────────
from .csumt import (
    bostick_depth,
    bostick_depth_from_rho,
    depth_coverage_table,
    frequency_for_depth,
    frequency_schedule,
    plot_depth_coverage_ranking,
    plot_depth_section,
    plot_frequency_schedule,
    vertical_resolution,
    vertical_resolution_pair,
)

# ─── Diagnostics and quality control ─────────────────────────────────────────
from .diag import (
    COVERAGE_THRESH,
    coverage_score,
    coverage_table,
    plot_polar_coverage,
    plot_polar_errors,
    plot_width_drift,
    rho_coverage,
    rho_error_stats,
)

# ─── Dimensionality analysis ──────────────────────────────────────────────────
from .dimensionality import (
    classify_dimensionality,
    encode_dimensionality,
    learn_dim_dictionary,
    mask_by_dictionary,
    mask_by_dimensionality,
    phase_features_table,
    plot_atom_psection,
    plot_dim_confidence_grid,
    plot_dim_map,
    plot_dim_occupancy_area,
    pre2d_inversion_assessment,
    project_to_2d,
)

# ─── Field zone classification ────────────────────────────────────────────────
from .fieldzone import (
    classify_field_zones,
    near_field_factor,
    plot_field_zones,
)

# ─── Frequency grid management ────────────────────────────────────────────────
from .frequency import (
    FrequencyEditResult,
    align_grid,
    decimate_step,
    drop_duplicates,
    drop_low_confidence_frequencies,
    edit_frequencies_by_confidence,
    frequency_edit_decision_table,
    frequency_edit_report,
    mask_low_confidence_frequencies,
    plot_apparent_depth_psection,
    plot_band_microstrips,
    plot_coverage_quality_heatmap,
    plot_frequency_edit_decisions,
    plot_frequency_edit_summary,
    recover_low_confidence_frequencies,
    regrid_logspace,
    regrid_to,
    select_band,
    smooth_mavg,
)
from .gb import (
    GroomBaileyResult,
    apply_groom_bailey,
    groom_bailey_decomposition,
    groom_bailey_table,
)

# ─── Gradient imaging (Zhang et al. 2021) ────────────────────────────────────
from .gradient_imaging import (
    plot_gradient_section,
    rho_frequency_gradient,
    rho_joint_gradient,
    rho_spatial_gradient,
)

# ─── Impedance diagnostics ────────────────────────────────────────────────────
from .impedance import (
    plot_determinant_track,
    plot_offdiag_antisym_residual,
    plot_phasor_wheel,
)

# ─── Data inspection ──────────────────────────────────────────────────────────
from .inspect import (
    frequency_coverage,
    list_missing_sections,
    plot_coverage,
    plot_rhoa_phi,
    plot_station_response,
    plot_survey_inventory_overview,
    plot_tipper_components,
    pseudosection,
    sites_summary,
)

# ─── L-curve ──────────────────────────────────────────────────────────────────
from .lcurve import (
    LCurveData,
    lcurve_from_mare2dem,
    lcurve_from_modem,
    lcurve_from_occam2d,
    lcurve_table,
    plot_lcurve,
)

# ─── Legacy / compatibility shims ─────────────────────────────────────────────
from .legacy import (
    align_tensor,
    check_em_kind,
    compute_qc,
    export_edis,
    extract_z_list,
    full_freq,
    parse_tensor,
    plot_confidence,
    plot_station_tensors,
    plot_strike,
    plot_tensors,
    tensor2d,
    wrap_phase,
)

# ─── Site-level visualisation ─────────────────────────────────────────────────
from .overview import plot_response_overview
from .plot import (
    plot_raw_sites_1d,
    plot_response_tipper,
    plot_sites_compare,
    plot_sites_fit_grid,
    plot_sites_panels,
)
from .qc import (
    build_qc_table,
    confidence_ratio,
    frequency_confidence_table,
    overlay_noise_cone,
    overlay_spectral_holes,
    plot_confidence_band_summary,
    plot_confidence_profile,
    plot_consistency_fan,
    plot_coverage_psection,
    plot_frequency_confidence_psection,
    plot_qc_quicklook,
    plot_snr_hist,
    plot_station_confidence_dashboard,
    plot_station_confidence_spectrum,
    plot_xyyx_crossover_map,
    qc_flags,
    station_confidence_table,
)

# ─── Noise removal ────────────────────────────────────────────────────────────
from .remove_noise import (
    EMAPFilterResult,
    apply_emap_filter,
    confidence_gated_emap_filter,
    correct_static_shift,
    drop_freqs_manual,
    emap_filter_report,
    emi_mitigation_report,
    enforce_offdiag_consistency,
    fixed_length_moving_average,
    hampel_filter_freq,
    mask_incoherent_freqs,
    notch_powerline,
    nr_qc_delta_offdiag_psection,
    nr_qc_harmonic_waterfall,
    nr_qc_snr_gain_profile,
    nr_qc_station_offdiag_curves,
    plot_emap_filter_profile,
    plot_emap_filter_psection,
    remove_noise_pipeline,
    rpca_offdiag_denoise,
    shrink_to_group_trend,
    smooth_logfreq,
    smooth_rho_phase,
    snr_table,
    spatial_median_filter,
    trimmed_moving_average,
)

# ─── Skewness analysis ────────────────────────────────────────────────────────
from .skew import (
    bahr_skewness,
    close_skew_gaps,
    keep_longest_low_skew,
    mask_by_skew,
    plot_skew_percentile_ribbon,
    plot_skew_traffic_psection,
    plot_skew_vote_band,
    plot_skewness,
    select_low_skew_band,
    skew_table,
)

# ─── Source array design ──────────────────────────────────────────────────────
from .source_array import (
    array_factor,
    beam_steer,
    pas_pattern,
    plot_radiation_pattern,
    sdas_directivity,
    sdas_element_pattern,
    snr_gain_db,
    steering_angles,
    wavenumber,
)

# ─── Source effects and near-field correction ─────────────────────────────────
from .source_effects import (
    BETA_THRESH_PCT,
    correct_near_field,
    detect_source_overprint,
    normalize_response,
    overprint_beta,
    plot_normalized_response,
    plot_overprint_section,
    source_overprint_table,
)

# ─── MobileMT (admittance tables, determinant σ_a, skew, band QC) ────────────
from .mobilemt import (
    admittance_determinant_table,
    admittance_skew_table,
    admittance_table,
    ensure_mobilemt_dataset,
    mask_outside_mobilemt_band,
    plot_mobilemt_admittance_profile,
    plot_mobilemt_conductivity_psection,
    plot_mobilemt_skew_profile,
)

# ─── ZTEM (along-profile divergence, phase rotation, band QC) ────────────────
from .ztem import (
    mask_outside_ztem_band,
    phase_rotate_table,
    plot_ztem_band_mask_psection,
    plot_ztem_divergence_profile,
    plot_ztem_divergence_psection,
    plot_ztem_divergence_psection_grid,
    plot_ztem_flight_lines,
    plot_ztem_map,
    plot_ztem_phase_rotation_profile,
    plot_ztem_tipper_profile,
    total_divergence_table,
    ztem_crossover_diagnostics,
)

# ─── Cross-spectra analysis and visualisation ────────────────────────────────
from .spectra import (
    band_select,
    coherence_matrix,
    coherence_table,
    mask_low_coherence,
    plot_coherence,
    plot_coherence_section,
    plot_psd,
    plot_psd_section,
    plot_spectra_matrix,
    plot_tipper_from_spectra,
    plot_z_from_spectra,
    psd_table,
    snr_table,
    spectra_summary,
)

# ─── Static shift correction ──────────────────────────────────────────────────
from .ss import (
    apply_ss_factors,
    correct_ss_ama,
    detect_near_surface,
    estimate_ss_ama,
    estimate_ss_bilateral,
    estimate_ss_loess,
    estimate_ss_refmedian,
    plot_ns_detection,
    plot_ss_1d_curves,
    plot_ss_comparison_psection,
    plot_ss_delta_profile,
    plot_ss_delta_psection,
    plot_ss_radar,
    plot_ss_station_curves,
    plot_ss_summary,
    ss_comparison_psection,
    ss_qc_profile,
    ss_qc_psection,
    ss_qc_station_curves,
)

# ─── Strike estimation ────────────────────────────────────────────────────────
from .strike import (
    estimate_strike_consensus,
    estimate_strike_phase_tensor,
    estimate_strike_sweep,
    plot_strike_analysis,
    plot_strike_mapsticks,
    plot_strike_profile,
    plot_strike_ribbon,
    plot_strike_rose,
    plot_strike_rose_by_line,
    rotate_to_strike,
    strike_curve_sweep,
)
from .tensor import (
    antisymmetrize,
    balance_offdiag,
    build_phase_tensor_table,
    invert,
    orient_from_sensors,
    phase_tensor_legend,
    plot_dimensionality_grid,
    plot_dimensionality_psection,
    plot_ellipticity_psection,
    plot_phase_tensor_map,
    plot_phase_tensor_psection,
    plot_phase_tensor_rose,
    plot_phase_tensor_skewmap,
    plot_phase_tensor_strip,
    plot_phase_tensor_strip_grid,
    plot_phase_tensor_summary,
    plot_skew_ellipt_density,
    plot_strike_director_field,
    plot_theta_rose_grid,
    plot_theta_stability_stripe,
    plot_theta_vs_period,
    rotate,
    rotate_by_map,
    sigma_clip_z,
)

# ─── Tensor operations ────────────────────────────────────────────────────────
# rotate_to_strike conflicts with strike.rotate_to_strike; alias as rotate_z_to_strike
from .tensor import (
    rotate_to_strike as rotate_z_to_strike,
)

# ─── Transfer functions — tipper & induction arrows ──────────────────────────
from .tf import (
    plot_induction_arrows,
    plot_induction_convention,
    # new: richer induction arrow views
    plot_induction_map,
    # spectra-direct wrappers (no Sites needed)
    plot_induction_map_from_spectra,
    # multi-period stacked map (Boukhalfa 2020 style)
    plot_induction_multiperiod_map,
    plot_induction_rose,
    plot_induction_rose_from_spectra,
    plot_induction_section,
    plot_tipper_hodograms,
    plot_tipper_polar,
    plot_tipper_polar_from_spectra,
)

__all__ = [
    # rose style
    "RoseStyle",
    "resolve_rose_style",
    # input helpers
    "ensure_sites",
    # AFMAG
    "euler_rotation_matrix",
    "geomagnetic_field_direction",
    "coil_normal_direction",
    "motion_coupling_cosine",
    "motion_coupling_angle",
    "simulate_motion_induced_voltage",
    "correct_motion_induced_noise",
    "afmag_tilt_angles",
    "motion_susceptibility_table",
    "flag_motion_susceptible_band",
    "plot_afmag_tilt_profile",
    "plot_afmag_tilt_psection",
    "plot_afmag_tilt_polar",
    "plot_motion_coupling",
    "plot_motion_susceptibility_map",
    "plot_afmag_correction_comparison",
    # survey design
    "bostick_depth_from_rho",
    "vertical_resolution_pair",
    "frequency_for_depth",
    "frequency_schedule",
    "plot_frequency_schedule",
    "bostick_depth",
    "vertical_resolution",
    "depth_coverage_table",
    "plot_depth_section",
    "plot_depth_coverage_ranking",
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
    "confidence_ratio",
    "frequency_confidence_table",
    "qc_flags",
    "station_confidence_table",
    "plot_confidence_band_summary",
    "plot_confidence_profile",
    "plot_frequency_confidence_psection",
    "plot_station_confidence_dashboard",
    "plot_station_confidence_spectrum",
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
    "smooth_rho_phase",
    "shrink_to_group_trend",
    "remove_noise_pipeline",
    "hampel_filter_freq",
    "spatial_median_filter",
    "fixed_length_moving_average",
    "trimmed_moving_average",
    "apply_emap_filter",
    "confidence_gated_emap_filter",
    "EMAPFilterResult",
    "emap_filter_report",
    "plot_emap_filter_profile",
    "plot_emap_filter_psection",
    "rpca_offdiag_denoise",
    "enforce_offdiag_consistency",
    "mask_incoherent_freqs",
    "drop_freqs_manual",
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
    "drop_low_confidence_frequencies",
    "edit_frequencies_by_confidence",
    "FrequencyEditResult",
    "frequency_edit_decision_table",
    "frequency_edit_report",
    "mask_low_confidence_frequencies",
    "plot_frequency_edit_decisions",
    "plot_frequency_edit_summary",
    "recover_low_confidence_frequencies",
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
    "pre2d_inversion_assessment",
    "mask_by_dimensionality",
    "project_to_2d",
    "learn_dim_dictionary",
    "encode_dimensionality",
    "mask_by_dictionary",
    "plot_atom_psection",
    "plot_dim_confidence_grid",
    "plot_dim_occupancy_area",
    "plot_dim_map",
    "GroomBaileyResult",
    "apply_groom_bailey",
    "groom_bailey_decomposition",
    "groom_bailey_table",
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
    "plot_strike_analysis",
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
    "plot_strike_director_field",
    "plot_ellipticity_psection",
    "plot_dimensionality_psection",
    "plot_phase_tensor_rose",
    "plot_phase_tensor_map",
    "phase_tensor_legend",
    "plot_dimensionality_grid",
    "plot_theta_stability_stripe",
    "plot_skew_ellipt_density",
    "plot_theta_rose_grid",
    "plot_phase_tensor_strip",
    "plot_phase_tensor_strip_grid",
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
    "plot_survey_inventory_overview",
    "plot_rhoa_phi",
    "plot_tipper_components",
    "pseudosection",
    "plot_station_response",
    # transfer functions / tipper / induction arrows
    "plot_tipper_hodograms",
    "plot_induction_arrows",
    "plot_induction_map",
    "plot_induction_section",
    "plot_induction_convention",
    "plot_tipper_polar",
    "plot_induction_rose",
    "plot_induction_map_from_spectra",
    "plot_tipper_polar_from_spectra",
    "plot_induction_rose_from_spectra",
    "plot_induction_multiperiod_map",
    # site-level visualisation
    "plot_response_overview",
    "plot_response_tipper",
    "plot_raw_sites_1d",
    "plot_sites_panels",
    "plot_sites_compare",
    "plot_sites_fit_grid",
    # l-curve
    "lcurve_table",
    "plot_lcurve",
    "LCurveData",
    "lcurve_from_occam2d",
    "lcurve_from_modem",
    "lcurve_from_mare2dem",
    # cross-spectra
    "coherence_matrix",
    "psd_table",
    "coherence_table",
    "snr_table",
    "band_select",
    "mask_low_coherence",
    "spectra_summary",
    "plot_psd",
    "plot_coherence",
    "plot_spectra_matrix",
    "plot_z_from_spectra",
    "plot_tipper_from_spectra",
    "plot_psd_section",
    "plot_coherence_section",
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
    "plot_strike",
    "plot_tensors",
    "plot_station_tensors",
    "wrap_phase",
    #Advanced plots 
    "plot_apparent_anisotropy_section",
    "plot_apparent_resistivity_polar",
    "plot_dimensionality_depth_profile",
    "plot_dimensionality_ternary",
    "plot_distortion_radar",
    "plot_impedance_mohr_circles",
    "plot_mt_composite_section",
    "plot_pt_period_clock",
    "plot_rho_phase_bode",
    "plot_sensitivity_depth_section",
    "plot_snr_section",
    "plot_strike_stability_bands",
    "plot_survey_fingerprint",
    "plot_tf_coherence_network",
    "plot_z_invariants_section",
    "plot_zt_argand",
    # ZTEM
    "ztem_crossover_diagnostics",
    "total_divergence_table",
    "phase_rotate_table",
    "mask_outside_ztem_band",
    "plot_ztem_tipper_profile",
    "plot_ztem_divergence_profile",
    "plot_ztem_divergence_psection",
    "plot_ztem_divergence_psection_grid",
    "plot_ztem_phase_rotation_profile",
    "plot_ztem_band_mask_psection",
    "plot_ztem_flight_lines",
    "plot_ztem_map",
    # MobileMT
    "ensure_mobilemt_dataset",
    "admittance_table",
    "admittance_determinant_table",
    "admittance_skew_table",
    "mask_outside_mobilemt_band",
    "plot_mobilemt_admittance_profile",
    "plot_mobilemt_conductivity_psection",
    "plot_mobilemt_skew_profile",
]
