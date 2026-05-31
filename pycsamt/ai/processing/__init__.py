# -*- coding: utf-8 -*-
"""
pycsamt.ai.processing
=====================

ML-based EM data preparation and QC.

Planned modules (Phase 3)
--------------------------
denoise.py
    :class:`EMDenoiser` — deep learning noise suppression on impedance
    tensor data.  Replaces or supplements conventional robust processing
    (Egbert 1997).

qc.py
    :class:`EMQCScorer` — automated per-frequency QC scoring that
    detects and flags noise bursts, static shifts, near-field
    contamination, and galvanic distortion.

anomaly.py
    :class:`AnomalyDetector` — unsupervised anomaly detection on
    tensor data across a profile (autoencoder-based).

classify.py
    :class:`DimensionalityClassifier` — classifies MT data as 1-D,
    2-D, or 3-D; also predicts strike direction.
    Replaces heuristic phase tensor / skew thresholds.
"""
__all__: list = []
