.. _emtools_ss:

pycsamt.emtools.ss — Static-Shift Estimation, Correction, and QC
=====================================================================

:mod:`pycsamt.emtools.ss` is the largest ``emtools`` module: four
independent static-shift estimators (adaptive moving-average, LOESS,
bilateral filtering, reference-median), factor application, a full
suite of before/after QC plots from simple bar charts to
publication-quality multi-panel pseudo-sections, a polar per-station
radar view, and a completely separate diagnostic — the Lei et al.
(2017) near-surface-vs-static-shift classifier, which distinguishes
distortion conventional static-shift correction *can* fix from
distortion it cannot.

Functions
---------

- :func:`~pycsamt.emtools.ss.estimate_ss_ama`
- :func:`~pycsamt.emtools.ss.estimate_ss_loess`
- :func:`~pycsamt.emtools.ss.estimate_ss_bilateral`
- :func:`~pycsamt.emtools.ss.estimate_ss_refmedian`
- :func:`~pycsamt.emtools.ss.apply_ss_factors`
- :func:`~pycsamt.emtools.ss.correct_ss_ama`
- :func:`~pycsamt.emtools.ss.plot_ss_delta_psection`
- :func:`~pycsamt.emtools.ss.plot_ss_station_curves`
- :func:`~pycsamt.emtools.ss.plot_ss_delta_profile`
- :func:`~pycsamt.emtools.ss.plot_ss_comparison_psection`
- :func:`~pycsamt.emtools.ss.plot_ss_1d_curves`
- :func:`~pycsamt.emtools.ss.plot_ss_summary`
- :func:`~pycsamt.emtools.ss.ss_qc_psection`
- :func:`~pycsamt.emtools.ss.ss_qc_station_curves`
- :func:`~pycsamt.emtools.ss.ss_qc_profile`
- :func:`~pycsamt.emtools.ss.ss_comparison_psection`
- :func:`~pycsamt.emtools.ss.plot_ss_radar`
- :func:`~pycsamt.emtools.ss.detect_near_surface`
- :func:`~pycsamt.emtools.ss.plot_ns_detection`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``), with a real bug fix
along the way: station ordering by ``sort_by="lon"``/``"lat"``
silently fell back to alphabetical-by-name for every real ``Site``
object (it only checked flat ``.lon``/``.lat`` attributes that do not
exist; real coordinates live in ``.coords``), now fixed — and, once
fixed, revealing that this line runs almost due north-south, so the
module's own default ``sort_by="lon"`` would still be the wrong axis
choice even with the bug gone. It moves from the raw station-to-station
spread motivating static-shift correction in the first place, through
the flagship AMA estimator (including the real, severe effect of this
survey's high phase-tensor skew starving the default threshold of
data), exact verification that factor application does what it
claims, a four-estimator cross-check that finds three methods agreeing
closely and a fourth disagreeing sharply — but only because of one
single station with a sign-flipped estimate, one-shot QC wrappers and
publication-quality multi-panel comparison figures, a polar radar view
shown first calm and then honestly noisy once per-frequency
phase-tensor rotation is turned on, and finishes with the Lei (2017)
near-surface-vs-static classifier and a check of how its own internal
consistency (rather than independent agreement) explains why it
correlates almost perfectly with the AMA table.

.. include:: auto_examples/plot_ss.rst
