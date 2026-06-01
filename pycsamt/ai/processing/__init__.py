# -*- coding: utf-8 -*-
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
"""
from .denoise import EMDenoiser, prepare_z_features
from .qc import EMQCScorer
from .anomaly import AnomalyDetector
from .classify import DimensionalityClassifier
from .plot import (
    plot_qc_scores,
    plot_qc_heatmap,
    plot_qc_feature_heatmap,
    plot_qc_summary,
)

__all__ = [
    "EMDenoiser",
    "prepare_z_features",
    "EMQCScorer",
    "AnomalyDetector",
    "DimensionalityClassifier",
    "plot_qc_scores",
    "plot_qc_heatmap",
    "plot_qc_feature_heatmap",
    "plot_qc_summary",
]
