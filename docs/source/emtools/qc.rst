.. _emtools_qc:

pycsamt.emtools.qc — Quality-Control Confidence Scoring
=============================================================

:mod:`pycsamt.emtools.qc` turns "does this transfer function look
trustworthy?" into numbers: per-station and per-frequency composite
confidence scores built from data coverage, tensor uncertainty,
off-diagonal consistency, diagonal leakage, phase smoothness, and
spatial coherence with neighboring stations — plus a wide family of
plots built on top of those scores, from a Kouadio et al. (2024)-style
along-line profile to an uncertainty-propagated apparent-resistivity
fan chart.

Functions
---------

- :func:`~pycsamt.emtools.qc.build_qc_table`
- :func:`~pycsamt.emtools.qc.qc_flags`
- :func:`~pycsamt.emtools.qc.station_confidence_table`
- :func:`~pycsamt.emtools.qc.frequency_confidence_table`
- :func:`~pycsamt.emtools.qc.plot_confidence_profile`
- :func:`~pycsamt.emtools.qc.plot_frequency_confidence_psection`
- :func:`~pycsamt.emtools.qc.plot_station_confidence_spectrum`
- :func:`~pycsamt.emtools.qc.plot_station_confidence_dashboard`
- :func:`~pycsamt.emtools.qc.plot_confidence_band_summary`
- :func:`~pycsamt.emtools.qc.plot_qc_quicklook`
- :func:`~pycsamt.emtools.qc.plot_consistency_fan`
- :func:`~pycsamt.emtools.qc.plot_xyyx_crossover_map`
- :func:`~pycsamt.emtools.qc.overlay_noise_cone`
- :func:`~pycsamt.emtools.qc.overlay_spectral_holes`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``), whose real ``z_err``
tensors make the uncertainty-aware scores and the Monte Carlo fan
chart genuine rather than illustrative. It moves from simple summary
tables (where a surprising finding shows up immediately: every station
gets flagged ``high_skew`` under the default threshold, since this
line's real Bibby skew runs 30-50 degrees — a structural signal, not a
data defect) through the ``presence``-vs-``composite`` confidence
contrast, the along-line profile and period-vs-station pseudo-section,
a worst-station/best-station dashboard comparison that separates *data
quality* from the *structural complexity* already explored in the
``anisotropy``/``impedance`` examples, a survey-wide period-band
summary, a coverage/SNR quicklook, and finishes with genuine
Monte Carlo error propagation and crossover/hole-detection diagnostics
on the same strongly anisotropic station featured throughout the
gallery.

.. include:: auto_examples/plot_qc.rst
