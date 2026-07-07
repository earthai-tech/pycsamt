.. _emtools_frequency:

pycsamt.emtools.frequency — Frequency Editing, Resampling, and QC
======================================================================

:mod:`pycsamt.emtools.frequency` is the toolbox for everything that
touches the frequency axis rather than a single station's values: band
selection, confidence-based dropping/masking/recovery of bad rows,
regridding/resampling, and coverage/quality/depth visualizations. It is
one of the largest ``emtools`` modules and central to any real
processing workflow, so this page is the most thorough one in the
gallery so far.

.. warning::

   Building this example found **two** independent, previously-untested
   bugs, both now fixed:

   1. :func:`~pycsamt.emtools.frequency.mask_low_confidence_frequencies`
      was a complete, silent no-op on the real ``Z`` class used
      throughout this library (0 rows masked at three different
      thresholds, no error raised) — a helper refused to write any
      NaN-containing array back to a "strict" Z-like object. The same
      guard also silently defeated
      :func:`~pycsamt.emtools.frequency.recover_low_confidence_frequencies`
      whenever a column couldn't be fully recovered. The existing tests
      missed this because their fake fixture lacked the attribute that
      triggers the guard.
   2. :func:`~pycsamt.emtools.frequency.select_band` called a
      single-site function from :mod:`pycsamt.site.edit` on a
      multi-site ``Sites`` container, silently leaving every station's
      frequency grid unchanged regardless of ``fmin``/``fmax`` — the
      same "single vs. broadcast" mix-up as the ``dimensionality``
      example's ``project_to_2d`` bug. The broadcast variant already
      existed and is now used instead.

Functions
---------

- :func:`~pycsamt.emtools.frequency.select_band`
- :func:`~pycsamt.emtools.frequency.drop_duplicates`
- :func:`~pycsamt.emtools.frequency.drop_low_confidence_frequencies`
- :func:`~pycsamt.emtools.frequency.mask_low_confidence_frequencies`
- :func:`~pycsamt.emtools.frequency.recover_low_confidence_frequencies`
- :func:`~pycsamt.emtools.frequency.edit_frequencies_by_confidence`
- :func:`~pycsamt.emtools.frequency.frequency_edit_report`
- :func:`~pycsamt.emtools.frequency.frequency_edit_decision_table`
- :func:`~pycsamt.emtools.frequency.plot_frequency_edit_summary`
- :func:`~pycsamt.emtools.frequency.plot_frequency_edit_decisions`
- :func:`~pycsamt.emtools.frequency.regrid_to`
- :func:`~pycsamt.emtools.frequency.regrid_logspace`
- :func:`~pycsamt.emtools.frequency.decimate_step`
- :func:`~pycsamt.emtools.frequency.smooth_mavg`
- :func:`~pycsamt.emtools.frequency.align_grid`
- :func:`~pycsamt.emtools.frequency.plot_coverage_quality_heatmap`
- :func:`~pycsamt.emtools.frequency.plot_apparent_depth_psection`
- :func:`~pycsamt.emtools.frequency.plot_band_microstrips`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``), with **KAP03**
(``data/MT/``) brought in for section 8 to show a real
heterogeneous-grid case. It moves from band selection and the
underlying confidence table (whose ceiling — 0.863, never reaching the
module's own default ``ci_hi=0.90`` — drives several of the findings
below), through drop/mask/recover (now working correctly; with the
default ``ci_hi`` every one of 1484 rows now correctly comes back NaN
rather than silently unchanged), the high-level
``edit_frequencies_by_confidence`` entry point and its report/decision
views, and finishes with regridding, decimation, smoothing, grid
alignment (where ``align_grid(mode="intersection")`` turns out to
silently no-op on real, independently-processed multi-station data —
by contract, not by bug, since no single frequency is bit-for-bit
identical across all of KAP03's 26 stations), and three richer
visualizations: a coverage/quality heatmap, an apparent-depth
pseudo-section, and band-collapsed microstrips.

.. include:: ../examples/emtools/plot_frequency.rst
   :start-after: .. _sphx_glr_examples_emtools_plot_frequency.py:
   :end-before: .. _sphx_glr_download_examples_emtools_plot_frequency.py:
