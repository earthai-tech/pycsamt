"""
pycsamt.ai.processing
=====================

ML-based EM data preparation and quality control.

Phase 3 modules
---------------
:class:`~pycsamt.ai.processing.denoise.EMDenoiser`
    Convolutional autoencoder for MT impedance tensor denoising.
    Suppresses broadband Gaussian noise and narrowband artefacts while
    preserving the smooth spectral shape of the true signal.

:class:`~pycsamt.ai.processing.qc.EMQCScorer`
    Automated per-(site, frequency) QC scoring.  Combines hard SNR /
    Swift-skew thresholds with an IsolationForest anomaly model.

:class:`~pycsamt.ai.processing.anomaly.AnomalyDetector`
    Profile-level unsupervised anomaly detection using a fully-connected
    autoencoder (PCA fallback when PyTorch is unavailable).

:class:`~pycsamt.ai.processing.classify.DimensionalityClassifier`
    MLP classifier that labels MT observations as 1-D, 2-D, or 3-D and
    predicts the geoelectric strike direction for 2-D sites.

:class:`~pycsamt.ai.processing.imputer.EMImputer`
    Masked-reconstruction gap filling for genuinely missing (not
    merely noisy) impedance cells, trained by synthetically hiding
    observed cells the way :class:`EMDenoiser` adds synthetic noise.

:class:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier`
    MLP triage classifier (clean / static-shift-only / distorted)
    routing a station toward :mod:`pycsamt.emtools.ss` or
    :mod:`pycsamt.emtools.gb` rather than re-solving the distortion
    physics itself.

:class:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator`
    Learned per-cell error-floor re-estimation: regresses a
    field-processing-consistent fractional error from the same
    signal-quality features :class:`EMQCScorer` extracts, then
    rescales ``Z.z_err`` toward that calibrated level via
    :meth:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator.apply`.

:class:`~pycsamt.ai.processing.tsdenoise.TimeSeriesDenoiser`
    MMF-SVM-K-SVD raw field-time-series denoising (Gui et al., 2024),
    one stage upstream of everything else here -- it works on a
    :class:`~pycsamt.ts.TSData` record (Ex, Ey, Hx, Hy, Hz samples)
    before :func:`~pycsamt.ts.process.ts_to_spectra` turns it into the
    impedance spectra :class:`EMDenoiser` and friends operate on.
    Combines mathematical morphological filtering, an SVM
    signal/noise triage classifier
    (:class:`~pycsamt.ai.processing.tsdenoise.SignalQualityClassifier`),
    and K-SVD dictionary-learning denoising
    (:class:`~pycsamt.ai.processing.tsdenoise.KSVDDenoiser`).

Visualisation
-------------
:mod:`pycsamt.ai.processing.plot` (all re-exported here) provides
dedicated figures for each tool above -- per-station bar charts,
station x frequency heat-maps, score distributions, denoiser
before/after spectra, a strike rose, a missing-data map with held-out
reconstruction validation for :class:`EMImputer`, a log-scaled
error-level map with a before/after comparison for
:class:`UncertaintyCalibrator`, and a per-station regime map with a
feature-space scatter for :class:`DistortionTypeClassifier` -- plus
:func:`~pycsamt.ai.processing.plot.plot_training_history`, shared by
every network-based estimator's ``history_`` attribute.
"""

from .anomaly import AnomalyDetector
from .classify import DimensionalityClassifier
from .denoise import EMDenoiser, prepare_z_features
from .distortion import (
    DistortionTypeClassifier,
    build_distortion_features_table,
)
from .imputer import EMImputer
from .plot import (
    plot_anomaly_score_distribution,
    plot_anomaly_scores,
    plot_anomaly_summary,
    plot_denoise_noise_reduction,
    plot_denoise_spectra,
    plot_denoise_summary,
    plot_dimensionality_map,
    plot_dimensionality_summary,
    plot_distortion_feature_space,
    plot_distortion_map,
    plot_distortion_summary,
    plot_imputer_gaps,
    plot_imputer_reconstruction,
    plot_imputer_summary,
    plot_imputer_validation,
    plot_predicted_strike_rose,
    plot_qc_feature_heatmap,
    plot_qc_heatmap,
    plot_qc_score_distribution,
    plot_qc_score_spread,
    plot_qc_scores,
    plot_qc_summary,
    plot_training_history,
    plot_ts_denoise_mmf_split,
    plot_ts_denoise_segments,
    plot_ts_denoise_summary,
    plot_uncertainty_map,
    plot_uncertainty_summary,
    plot_uncertainty_validation,
)
from .qc import EMQCScorer
from .tsdenoise import (
    KSVDDenoiser,
    SignalQualityClassifier,
    TimeSeriesDenoiser,
    compute_entropy_features,
    generate_synthetic_library,
    mmf_split,
    snr_db,
    time_domain_ncc,
)
from .uncertainty import (
    UncertaintyCalibrator,
    build_uncertainty_features_table,
)

__all__ = [
    "EMDenoiser",
    "prepare_z_features",
    "EMQCScorer",
    "AnomalyDetector",
    "DimensionalityClassifier",
    "EMImputer",
    "UncertaintyCalibrator",
    "DistortionTypeClassifier",
    "build_distortion_features_table",
    "build_uncertainty_features_table",
    "TimeSeriesDenoiser",
    "SignalQualityClassifier",
    "KSVDDenoiser",
    "mmf_split",
    "compute_entropy_features",
    "generate_synthetic_library",
    "snr_db",
    "time_domain_ncc",
    "plot_qc_scores",
    "plot_qc_heatmap",
    "plot_qc_feature_heatmap",
    "plot_qc_score_distribution",
    "plot_qc_score_spread",
    "plot_qc_summary",
    "plot_denoise_spectra",
    "plot_denoise_noise_reduction",
    "plot_denoise_summary",
    "plot_anomaly_scores",
    "plot_anomaly_score_distribution",
    "plot_anomaly_summary",
    "plot_dimensionality_map",
    "plot_dimensionality_summary",
    "plot_predicted_strike_rose",
    "plot_imputer_gaps",
    "plot_imputer_validation",
    "plot_imputer_reconstruction",
    "plot_imputer_summary",
    "plot_uncertainty_map",
    "plot_uncertainty_validation",
    "plot_uncertainty_summary",
    "plot_distortion_map",
    "plot_distortion_feature_space",
    "plot_distortion_summary",
    "plot_ts_denoise_mmf_split",
    "plot_ts_denoise_segments",
    "plot_ts_denoise_summary",
    "plot_training_history",
]
