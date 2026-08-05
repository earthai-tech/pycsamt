.. _emtools_qc:

Quality-Control Confidence Scoring
==================================

``pycsamt.emtools.qc`` turns transfer-function :term:`quality control`
into tables and figures. It does two related jobs:

* summarize station-level data quality, coverage, :term:`SNR`, tipper
  presence, and phase-tensor :term:`skew`;
* compute confidence ratios from several finite, bounded scores so that
  stations and individual station-frequency samples can be ranked,
  plotted, masked, or down-weighted before inversion.

Full callable signatures live in the :doc:`API reference <../../api/emtools>`.
This page explains the confidence-ratio formulation, the table
workflow, the plotting workflow, and how confidence values should be
read in practice. Every example below uses pyCSAMT's bundled
:term:`AMT` line, ``data/AMT/WILLY_DATA/L18PLT``, since all 28 stations
there carry real impedance-error tensors -- the ingredient several of
the component scores need.

Why QC Is More Than Coverage
----------------------------

A station can be complete and still be questionable. ``frac_ok=1.0``
only says that impedance rows are finite. It does not say that the
off-diagonal modes agree, that diagonal leakage is small, that phase is
smooth, that uncertainty is low, or that the station is coherent with
its neighbours.

The QC module therefore separates three ideas:

``coverage``
    Are the required transfer-function rows finite?

``confidence``
    How trustworthy is a station or a station-frequency cell after
    combining coverage, uncertainty, tensor-shape, phase, and spatial
    criteria?

``flags``
    Which simple thresholds does a station or frequency cell fail?

Load A Survey
-------------

All QC functions accept the usual pyCSAMT inputs, but a reproducible
script should normalize once with ``ensure_sites``.

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.emtools import ensure_sites
   >>> survey = ensure_sites(
   ...     Path("data/AMT/WILLY_DATA/L18PLT"),
   ...     recursive=True,
   ...     on_dup="replace",
   ...     strict=True,
   ...     verbose=1,
   ... )

Use ``strict=True`` for reports and automated processing. Use
``strict=False`` in exploratory notebooks if you want empty plots to
render as "no data" messages.

The Confidence Ratio
--------------------

The composite confidence ratio is a weighted mean of available
component scores:

.. math::

   \mathrm{CR}_{i,f} =
   { \sum_k w_k s_{k,i,f}
     \mathbf{1}_{s_{k,i,f}\ \mathrm{finite}} \over
     \sum_k w_k
     \mathbf{1}_{s_{k,i,f}\ \mathrm{finite}} },
   \qquad 0 \le s_k \le 1.

Missing scores are ignored. Finite scores are clipped to ``[0, 1]``.
The default weights are:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Score
     - Weight
     - Meaning
   * - ``coverage``
     - ``0.35``
     - Finite impedance rows or components.
   * - ``uncertainty``
     - ``0.20``
     - Median relative impedance error.
   * - ``offdiag``
     - ``0.15``
     - Similarity of ``Zxy`` and ``Zyx`` amplitudes.
   * - ``diagonal``
     - ``0.10``
     - Penalty for diagonal leakage into off-diagonal response.
   * - ``phase``
     - ``0.10``
     - Penalty for abrupt off-diagonal phase jumps.
   * - ``spatial``
     - ``0.10``
     - Coherence with neighbouring stations.

Every component follows the same shape. Each one reduces the tensor to a
single non-negative "badness" statistic :math:`m_k`, compares it against
a tolerance :math:`\tau_k`, and folds the excess back into a
:math:`[0, 1]` trust score:

.. math::

   s_k = \mathrm{clip}_{[0,1]}\!\left(1 - \frac{m_k}{\tau_k}\right).

A statistic of ``0`` scores a perfect ``1``; as :math:`m_k` grows toward
:math:`\tau_k` the score decays linearly to ``0`` and stays there beyond
it, so nothing is ever penalized twice or rewarded for being
implausibly clean. Concretely, the six statistics are:

- ``coverage`` needs no tolerance -- it is already the
  :term:`Finite coverage` fraction of the tensor, in ``[0, 1]``.
- ``uncertainty`` uses :math:`m = \mathrm{median}(|Z_{\mathrm{err}}| /
  |Z|)` against ``relerr_threshold`` (default ``0.20``): how large the
  reported error is relative to the signal.
- ``offdiag`` uses :math:`m = \left|\log_{10}(|Z_{xy}| /
  |Z_{yx}|)\right|` against ``offdiag_tolerance_log10`` (default
  ``0.35``): how far the two :term:`off-diagonal component` modes
  disagree in log amplitude.
- ``diagonal`` uses :math:`m = \mathrm{med}(|Z_{xx}|, |Z_{yy}|) \big/
  \left[\mathrm{med}(|Z_{xx}|, |Z_{yy}|) + \mathrm{med}(|Z_{xy}|,
  |Z_{yx}|)\right]` against ``diagonal_leakage_max`` (default ``0.35``):
  the fraction of tensor amplitude sitting on the diagonal, where a 1-D
  or 2-D earth should keep most of it off-diagonal.
- ``phase`` uses :math:`m = \mathrm{med}|\Delta\varphi_{xy,yx}|`, the
  median absolute jump in degrees between the unwrapped phase of
  neighbouring frequencies, against ``phase_jump_tolerance_deg``
  (default ``90``).
- ``spatial`` uses :math:`m = \left|\log_{10}\rho -
  \mathrm{median}(\log_{10}\rho_{\mathrm{neighbours}})\right|` against
  ``spatial_tolerance_log10`` (default ``0.60``), where :math:`\rho` is a
  determinant-style :term:`apparent resistivity` proxy compared against
  the immediately adjacent stations (or, at the frequency level, the
  same frequency at the two nearest stations by distance).

At the station level each :math:`m_k` is one median taken over every
frequency row, which is why ``station_confidence_table`` returns one
score per station. At the frequency level, further down, the same six
formulas are evaluated one row at a time instead, so a station that
looks coherent overall can still carry a handful of untrustworthy
individual frequencies underneath.

The default confidence bands are:

* ``CR >= 0.95``: safe / retained;
* ``0.85 <= CR < 0.95``: recoverable or marginal;
* ``CR < 0.85``: reject, down-weight, or manually review.

Compute A Confidence Ratio Directly
-----------------------------------

Use ``confidence_ratio`` when you already have scores and want to apply
the same weighted formula used by the tables.

.. code-block:: pycon

   >>> from pycsamt.emtools.qc import confidence_ratio
   >>> scores = {
   ...     "coverage": 1.00,
   ...     "uncertainty": 0.82,
   ...     "offdiag": 0.76,
   ...     "diagonal": 0.55,
   ...     "phase": 0.90,
   ...     "spatial": 0.88,
   ... }
   >>> cr, cr_err = confidence_ratio(scores, n_freq=53, return_error=True)
   >>> print(f"CR={cr:.3f} +/- {cr_err:.3f}")
   CR=0.861 +/- 0.141

``confidence_err`` is not a formal statistical uncertainty on
:math:`\mathrm{CR}` -- it is a cheap, honest stand-in for one:

.. math::

   \sigma_{\mathrm{CR}} =
   \begin{cases}
   \mathrm{std}\bigl(\{s_k\}\bigr), & \text{two or more finite } s_k, \\[4pt]
   \sqrt{\mathrm{CR}\,(1 - \mathrm{CR}) / n_{\mathrm{freq}}}, & \text{otherwise.}
   \end{cases}

With six components on hand, as in the example above,
:math:`\sigma_{\mathrm{CR}}` is simply the population standard deviation
of ``{1.00, 0.82, 0.76, 0.55, 0.90, 0.88}``, which is exactly the
``0.141`` printed above -- the components disagree with each other, and
that disagreement *is* the reported uncertainty. When only one
component survives (most often ``coverage`` alone, when no error tensor
or neighbouring station exists to score against), there is nothing to
disagree with, so the formula falls back to the binomial-style standard
error :math:`\sqrt{\mathrm{CR}(1-\mathrm{CR})/n_{\mathrm{freq}}}` instead.

Station QC Summary
------------------

``build_qc_table`` is the first station-level table. It reports the
number of frequencies, finite-row coverage, tipper availability, median
row :term:`SNR` when ``z_err`` exists, period range, and optional
phase-tensor :term:`skew`.

.. code-block:: pycon

   >>> from pycsamt.emtools import build_qc_table, ensure_sites
   >>> survey = ensure_sites("data/AMT/WILLY_DATA/L18PLT", strict=True)
   >>> qc = build_qc_table(survey, include_skew=True, api=False)
   >>> qc[["station", "n_freq", "n_ok", "frac_ok", "n_tip", "snr_med", "pmin", "pmax", "skew_med"]].head()
      station  n_freq  n_ok  frac_ok  ...    snr_med      pmin      pmax   skew_med
   0  18-001A      53    53      1.0  ...  17.658396  0.000096  0.992063   4.809977
   1  18-002U      53    53      1.0  ...  16.687366  0.000096  0.992063   6.596544
   2  18-003A      53    53      1.0  ...  12.031672  0.000096  0.992063  11.018588
   3  18-004A      53    53      1.0  ...  10.430580  0.000096  0.992063  14.509282
   4  18-005U      53    53      1.0  ...  14.360341  0.000096  0.992063   9.171151
   <BLANKLINE>
   [5 rows x 9 columns]

Read ``frac_ok`` as completeness, not as full confidence. Read
``snr_med`` as a row-level signal-to-error ratio when impedance error
tensors are available. Read ``skew_med`` as structural/tensor
complexity, not automatically as bad acquisition -- it is the median
absolute :term:`phase tensor` asymmetry angle :math:`|\beta|` across a
station's frequencies, the same quantity the :term:`skew` glossary
entry defines.

Station Flags
-------------

``qc_flags`` adds simple threshold labels to the station QC table.

.. code-block:: pycon

   >>> from pycsamt.emtools import qc_flags
   >>> flags = qc_flags(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     min_frac_ok=0.60,
   ...     min_snr_med=2.0,
   ...     max_skew_med=6.0,
   ... )
   >>> flagged = flags[flags["flags"] != ""]
   >>> print(len(flagged), "of", len(flags), "stations flagged")
   27 of 28 stations flagged
   >>> flagged[["station", "frac_ok", "snr_med", "skew_med", "flags"]].head()
      station  frac_ok    snr_med   skew_med      flags
   1  18-002U      1.0  16.687366   6.596544  high_skew
   2  18-003A      1.0  12.031672  11.018588  high_skew
   3  18-004A      1.0  10.430580  14.509282  high_skew
   4  18-005U      1.0  14.360341   9.171151  high_skew
   5  18-006A      1.0  13.272516  12.375357  high_skew

Every station is fully covered (``frac_ok=1.0`` throughout) and every
station's median row SNR clears ``min_snr_med=2.0`` by a wide margin, so
``high_skew`` is the only flag this survey ever raises, and it raises it
for all but one station: only ``18-001A`` (``skew_med=4.81``) falls
under the default ``max_skew_med=6.0`` threshold, while the rest run
from ``6.60`` up past ``50`` degrees. That is not a data defect -- it is
a genuinely 2-D/3-D structural line, expressed through a threshold tuned
for near-1-D settings. Possible station flags include ``low_coverage``,
``low_snr``, and ``high_skew``. A high-skew flag can be a real
structural signal. It does not mean the station should automatically be
deleted.

Presence Confidence Versus Composite Confidence
-----------------------------------------------

``station_confidence_table`` has two modes. ``method="presence"`` is
coverage-only. ``method="composite"`` combines all available component
scores.

.. code-block:: pycon

   >>> from pycsamt.emtools import station_confidence_table
   >>> presence = station_confidence_table(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     method="presence",
   ...     api=False,
   ... )
   >>> composite = station_confidence_table(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     method="composite",
   ...     api=False,
   ... )
   >>> print("presence range", presence["confidence"].min(), presence["confidence"].max())
   presence range 1.0 1.0
   >>> print("composite range", composite["confidence"].min(), composite["confidence"].max())
   composite range 0.5440199944767435 0.8119223880398303
   >>> ranked = composite.sort_values("confidence")
   >>> ranked[["station", "confidence", "coverage", "uncertainty", "offdiag", "diagonal", "phase", "spatial"]].head()
       station  confidence  coverage  ...  diagonal     phase   spatial
   22  18-022U    0.544020       1.0  ...  0.000000  0.941638  0.000000
   21  18-021B    0.574114       1.0  ...  0.000000  0.937891  0.232771
   17  18-018A    0.578410       1.0  ...  0.086979  0.957546  0.000000
   20  18-021U    0.594344       1.0  ...  0.000000  0.938165  0.233576
   16  18-017U    0.595479       1.0  ...  0.261699  0.969783  0.038160
   <BLANKLINE>
   [5 rows x 8 columns]

Presence confidence is a flat ``1.0`` everywhere -- every row is
finite, so a coverage-only view has nothing left to say about this
particular survey. Composite confidence spreads from
:math:`\approx 0.54` to :math:`\approx 0.81` instead: the coverage-only
view was hiding real, measurable quality variation. If presence
confidence is high everywhere but composite confidence varies, the
survey is complete but not equally trustworthy everywhere. That is
common in real EM data.

Customize Confidence Weights
----------------------------

Use custom weights when a project has a clear processing policy. For
example, inversion preparation may emphasize uncertainty and coverage,
while structural interpretation may care more about off-diagonal and
spatial coherence.

.. code-block:: pycon

   >>> weights = {
   ...     "coverage": 0.40,
   ...     "uncertainty": 0.30,
   ...     "offdiag": 0.10,
   ...     "diagonal": 0.05,
   ...     "phase": 0.05,
   ...     "spatial": 0.10,
   ... }
   >>> table = station_confidence_table(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     method="composite",
   ...     weights=weights,
   ...     relerr_threshold=0.25,
   ...     offdiag_tolerance_log10=0.40,
   ...     diagonal_leakage_max=0.40,
   ...     phase_jump_tolerance_deg=90.0,
   ...     spatial_tolerance_log10=0.60,
   ...     api=False,
   ... )
   >>> table.sort_values("confidence")[["station", "distance_m", "confidence", "coverage", "uncertainty", "offdiag", "diagonal", "phase", "spatial"]].head()
       station  distance_m  confidence  ...  diagonal     phase   spatial
   22  18-022U      4400.0    0.630495  ...  0.000000  0.941638  0.000000
   21  18-021B      4200.0    0.658344  ...  0.000000  0.937891  0.232771
   17  18-018A      3400.0    0.666682  ...  0.201107  0.957546  0.000000
   16  18-017U      3200.0    0.672222  ...  0.353986  0.969783  0.038160
   20  18-021U      4000.0    0.682870  ...  0.000000  0.938165  0.233576
   <BLANKLINE>
   [5 rows x 9 columns]

The worst station is still ``18-022U`` under either weighting -- its
penalty is concentrated enough (zero ``diagonal`` and near-zero
``spatial``) that it stays worst regardless of emphasis. Changing
thresholds changes the meaning of the scores. Record custom weights and
thresholds in reports so another user can reproduce your confidence
classes.

Frequency-Level Confidence
--------------------------

``frequency_confidence_table`` returns one row per station-frequency
sample, scored with the same six formulas and the same
:math:`\mathrm{CR}_{i,f}` weighting as the station table above -- the
only change is that each :math:`m_k` is now evaluated at a single
frequency rather than medianed across the whole station. This is the
table to use for period-band decisions, masks, and inversion
down-weighting.

.. code-block:: pycon

   >>> from pycsamt.emtools import frequency_confidence_table
   >>> freq_qc = frequency_confidence_table(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     method="composite",
   ...     ci_hi=0.95,
   ...     ci_lo=0.85,
   ...     api=False,
   ... )
   >>> freq_qc.columns.tolist()
   ['station', 'station_index', 'distance_m', 'frequency_hz', 'period_s', 'log10_period', 'confidence', 'confidence_err', 'method', 'n_components', 'coverage', 'uncertainty', 'offdiag', 'diagonal', 'phase', 'spatial', 'logrho_proxy', 'flags']
   >>> freq_qc[["station", "frequency_hz", "period_s", "confidence", "flags"]].head()
      station  ...                                              flags
   0  18-001A  ...  recoverable,high_error,offdiag_mismatch,diagon...
   1  18-001A  ...  reject,high_error,offdiag_mismatch,diagonal_le...
   2  18-001A  ...  reject,high_error,offdiag_mismatch,diagonal_le...
   3  18-001A  ...  reject,high_error,offdiag_mismatch,diagonal_le...
   4  18-001A  ...  reject,high_error,offdiag_mismatch,diagonal_le...
   <BLANKLINE>
   [5 rows x 5 columns]
   >>> rejected = freq_qc[freq_qc["flags"].str.contains("reject", na=False)]
   >>> print("rejected cells:", len(rejected), "of", len(freq_qc))
   rejected cells: 1479 of 1484

Frequency flags can include ``reject``, ``recoverable``, ``missing``,
``high_error``, ``offdiag_mismatch``, ``diagonal_leakage``,
``phase_jump``, and ``spatial_outlier``. Nearly every one of this
survey's 1484 station-frequency cells reads ``reject`` at the frequency
level, a much harsher picture than the station-level composite range of
0.54-0.81 above -- station scores are medians over frequency, so a
station can stay "recoverable" overall while most of its individual
frequencies fail one or more component thresholds underneath.

Build A Mask From Confidence
----------------------------

The QC module does not force one masking policy. You can build one from
the confidence table.

.. code-block:: pycon

   >>> table = frequency_confidence_table(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     method="composite",
   ...     ci_lo=0.85,
   ...     api=False,
   ... )
   >>> keep = table["confidence"] >= 0.85
   >>> review = (table["confidence"] >= 0.70) & (table["confidence"] < 0.85)
   >>> drop = table["confidence"] < 0.70
   >>> print("keep:", int(keep.sum()))
   keep: 5
   >>> print("review:", int(review.sum()))
   review: 533
   >>> print("drop:", int(drop.sum()))
   drop: 946

Use a review band instead of a hard delete when the flagged frequencies
line up with known structural complexity. Low confidence is a prompt
for inspection, not always proof of bad data.

Station Confidence Profile
--------------------------

``plot_confidence_profile`` plots station confidence along the line. It
uses green, pink, and red markers for safe, recoverable, and rejected
stations.

.. code-block:: pycon

   >>> import matplotlib.pyplot as plt
   >>> from pycsamt.emtools import plot_confidence_profile
   >>> fig, ax = plt.subplots(figsize=(9.5, 4.2))
   >>> _ = plot_confidence_profile(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     method="composite",
   ...     ci_hi=0.95,
   ...     ci_lo=0.85,
   ...     shade_mode="score",
   ...     station_label_step=2,
   ...     show_errorbars=True,
   ...     ax=ax,
   ... )
   >>> fig.tight_layout()
   >>> fig.savefig("confidence_profile_l18plt.png", dpi=200)
   >>> plt.close(fig)

.. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-09.png
   :width: 100%

Every station in this survey lands in the pink "recoverable" band
(:math:`0.50 \le \mathrm{CR} < 0.95`), consistent with the 0.54-0.81
range printed earlier -- none rejected, none fully safe. If no station
coordinate metadata are available, distance falls back to regular
spacing controlled by ``spacing_m``, which is exactly what happens here
since these EDI objects expose latitude/longitude rather than
projected easting/northing, so the x-axis reads as station order along
the line rather than a surveyed distance.

Frequency Confidence Pseudo-Section
-----------------------------------

``plot_frequency_confidence_psection`` shows confidence by station and
period. It can plot any metric column from the frequency table, not only
``confidence``.

.. code-block:: pycon

   >>> from pycsamt.emtools import plot_frequency_confidence_psection
   >>> fig, ax = plt.subplots(figsize=(10.0, 4.8))
   >>> _ = plot_frequency_confidence_psection(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     method="composite",
   ...     metric="confidence",
   ...     station_label_step=2,
   ...     ax=ax,
   ... )
   >>> fig.savefig("frequency_confidence_psection.png", dpi=200)
   >>> plt.close(fig)

.. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-10.png
   :width: 100%

Change ``metric`` to ``"uncertainty"``, ``"offdiag"``, ``"diagonal"``,
``"phase"``, or ``"spatial"`` to see which component is driving low
confidence at a given station and period, rather than only the combined
score.

Single-Station Spectrum
-----------------------

``plot_station_confidence_spectrum`` overlays the overall confidence
curve and component scores for one station.

.. code-block:: pycon

   >>> from pycsamt.emtools import plot_station_confidence_spectrum
   >>> fig, ax = plt.subplots(figsize=(7.5, 4.2))
   >>> _ = plot_station_confidence_spectrum(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     station="18-022U",
   ...     method="composite",
   ...     ax=ax,
   ... )
   >>> fig.savefig("station_confidence_spectrum_18-022U.png", dpi=200)
   >>> plt.close(fig)

.. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-11.png
   :width: 100%

``18-022U`` is this survey's lowest-confidence station from the ranked
table above. Use this plot when you know a station is weak and want to
see whether the problem comes from uncertainty, diagonal leakage,
off-diagonal mismatch, phase jumps, or spatial incoherence, before
reaching for the dashboard's split-panel view.

Single-Station Dashboard
------------------------

``plot_station_confidence_dashboard`` breaks the same information into
separate panels. It is easier to read than the overlaid spectrum when
many score components overlap.

.. code-block:: pycon

   >>> from pycsamt.emtools import plot_station_confidence_dashboard
   >>> fig = plot_station_confidence_dashboard(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     station="18-022U",
   ...     method="composite",
   ...     ci_hi=0.95,
   ...     ci_lo=0.85,
   ...     figsize=(11.0, 7.0),
   ... )
   >>> fig.savefig("station_confidence_dashboard_18-022U.png", dpi=200)
   >>> plt.close(fig)

.. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-12.png
   :width: 100%

Split into six panels, the story reads clearly: "Data coverage" stays a
flat, uninformative 1.0 throughout (top-middle), while "Offdiag
consistency", "Diagonal leakage", and "Phase + spatial coherence"
(bottom row) all repeatedly collapse toward zero -- the composite
penalty is concentrated in tensor-shape and spatial diagnostics, not
missing data. Use dashboards for station-by-station review before
deciding whether a low-confidence station should be edited,
down-weighted, or retained.

Period-Band Summary
-------------------

``plot_confidence_band_summary`` collapses the frequency table by
period. It plots median and mean confidence and shades the fraction of
stations in rejected or recoverable bands.

.. code-block:: pycon

   >>> from pycsamt.emtools import plot_confidence_band_summary
   >>> fig, ax = plt.subplots(figsize=(8.5, 4.2))
   >>> _ = plot_confidence_band_summary(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     method="composite",
   ...     ci_hi=0.95,
   ...     ci_lo=0.85,
   ...     ax=ax,
   ... )
   >>> fig.savefig("confidence_band_summary.png", dpi=200)
   >>> plt.close(fig)

.. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-13.png
   :width: 100%

This view is useful when deciding whether a whole period band should be
edited, down-weighted, or treated cautiously -- it is also the figure
that makes the "1479 of 1484 rejected cells" number from the frequency
table above legible instead of abstract: the rejected-fraction shading
runs high across nearly the entire period range shown here.

Coverage And SNR Quicklook
--------------------------

``plot_qc_quicklook`` combines three first-pass plots: a presence
pseudo-section, an SNR pseudo-section, and an SNR histogram.

.. code-block:: pycon

   >>> from pycsamt.emtools.qc import plot_qc_quicklook
   >>> fig = plot_qc_quicklook(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     figsize=(10.0, 8.0),
   ... )
   >>> fig.savefig("qc_quicklook_l18plt.png", dpi=200)
   >>> plt.close(fig)

.. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-14.png
   :width: 100%

The top panel is solid green: 100% row presence everywhere, the same
``frac_ok=1.0`` finding from the station table. The bottom-left SNR
pseudo-section is where the real texture is -- brighter (higher row
:term:`SNR`, :math:`|Z|/\sigma`) around the shorter periods and near a
few stations, fading elsewhere -- and the histogram on the right shows
the underlying distribution peaking in the low tens with a long tail. If
the SNR histogram says error tensors are not available, the survey can
still be inspected, but uncertainty-based scores will be missing from
the composite confidence ratio.

Coverage Pseudo-Section And SNR Histogram
-----------------------------------------

For more control, call the helper plots separately.

.. code-block:: pycon

   >>> from pycsamt.emtools.qc import plot_coverage_psection, plot_snr_hist
   >>> fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5))
   >>> _ = plot_coverage_psection(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     metric="presence",
   ...     ax=axes[0],
   ... )
   >>> _ = axes[0].set_title("Finite-row presence")
   >>> _ = plot_snr_hist(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     bins=40,
   ...     ax=axes[1],
   ... )
   >>> _ = axes[1].set_title("Row SNR distribution")
   >>> fig.tight_layout()
   >>> fig.savefig("coverage_and_snr.png", dpi=200)
   >>> plt.close(fig)

.. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-15.png
   :width: 100%

``metric="presence"`` shows finite rows. ``metric="snr"`` colours by
row SNR when ``z_err`` exists. ``metric="offdiag"`` shows an
off-diagonal amplitude proxy.

Consistency Fan
---------------

Apparent resistivity is a nonlinear function of impedance, so a
symmetric error on :math:`Z` does not become a symmetric error on
:math:`\rho_a` -- the ``uncertainty`` score above summarizes that error
as one number, but it cannot show its shape. ``plot_consistency_fan``
shows the shape directly by Monte Carlo sampling: it draws
:math:`n_{\mathrm{draws}}` complex Gaussian perturbations with standard
deviation equal to the reported error,

.. math::

   Z^{(d)} = Z + E^{(d)}, \qquad
   E^{(d)} \sim \mathcal{CN}(0,\, |Z_{\mathrm{err}}|^2), \qquad
   d = 1, \dots, n_{\mathrm{draws}},

recomputes :math:`\rho_a^{(d)} = 0.2\,|Z^{(d)}|^2 / f` for each draw, and
plots the requested percentiles of the resulting distribution -- the
``10``/``50``/``90`` default draws a band around the median rather than
a single curve. It can compare ``xy`` and ``yx``
:term:`apparent resistivity` bands for one station, which is exactly
where the ``offdiag`` score above would flag a problem, but here you
see the two bands and can judge for yourself whether they merely
disagree or actually fail to overlap.

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.emtools.qc import overlay_noise_cone, plot_consistency_fan
   >>> fig, ax = plt.subplots(figsize=(8.8, 4.5))
   >>> _ = plot_consistency_fan(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     station="18-016A",
   ...     comps=("xy", "yx"),
   ...     pcts=(10.0, 50.0, 90.0),
   ...     n_draws=300,
   ...     ax=ax,
   ... )
   >>> ax.set_yscale("log")
   >>> period = np.logspace(-4, 0, 30)
   >>> overlay_noise_cone(
   ...     ax,
   ...     period,
   ...     lo=np.full(period.size, 10.0),
   ...     hi=np.full(period.size, 100.0),
   ...     color="0.5",
   ...     alpha=0.20,
   ... )
   >>> fig.savefig("consistency_fan_18-016A.png", dpi=200)
   >>> plt.close(fig)

.. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-16.png
   :width: 100%

``18-016A`` is the same station flagged elsewhere in this survey for
strong ratio anisotropy: on this log axis, :math:`\rho_{a,xy}` climbs
into the tens of thousands of :math:`\Omega\,\mathrm{m}` while
:math:`\rho_{a,yx}` stays two to three decades lower throughout, and the
shaded Monte Carlo bands around each curve are genuinely propagated
from the EDI's own error tensor rather than a linearized approximation.
The grey noise cone overlay is a visual reference band, not estimated
by the QC module -- it happens to bracket most of this station's real
:math:`\rho_{a,yx}` values here while sitting far below
:math:`\rho_{a,xy}`. Supply project-specific lower and upper bounds when
you use it in a report.

XY/YX Crossover Map
-------------------

Away from a 1-D earth, ``rho_xy`` and ``rho_yx`` need not agree, but
which one is larger should not flip back and forth erratically with
period. ``plot_xyyx_crossover_map`` tracks the sign of
:math:`d(f) = \rho_{a,xy}(f) - \rho_{a,yx}(f)` along each station's
sounding and marks every period where it changes,

.. math::

   \operatorname{sign}\bigl(d(f_j)\bigr) \ne
   \operatorname{sign}\bigl(d(f_{j+1})\bigr),

placing the marker at the log-period linearly interpolated between
:math:`f_j` and :math:`f_{j+1}`, weighted by how close each side came to
zero. A station with one or two crossovers across a smooth sounding is
unremarkable; a station peppered with them is a cheap, effective
anisotropy and mode-consistency diagnostic worth following up with the
consistency fan above.

.. code-block:: pycon

   >>> from pycsamt.emtools.qc import overlay_spectral_holes, plot_xyyx_crossover_map
   >>> fig, ax = plt.subplots(figsize=(9.5, 4.8))
   >>> _ = plot_xyyx_crossover_map(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     ax=ax,
   ... )
   >>> overlay_spectral_holes(
   ...     ax,
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     thresh_dec=0.30,
   ... )
   >>> fig.savefig("xy_yx_crossover_map.png", dpi=200)
   >>> plt.close(fig)

.. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-17.png
   :width: 100%

``overlay_spectral_holes`` shades large gaps in log-period sampling on
top of a pseudo-section-style axis. L18PLT's real frequency grid is
dense enough that the default 0.30-decade threshold finds nothing to
shade here -- an honest negative result rather than a broken overlay.
Lower ``thresh_dec`` only when you intentionally want to reveal small
grid-spacing differences instead.

Propagation To Inversion
------------------------

For :term:`MARE2DEM` exports created from :term:`EDI` data, CR-derived
uncertainty propagation can be enabled with:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.models.mare2dem.edi import make_mt_data_from_edi
   >>> out_path = Path("mare2dem_data_with_confidence.emdata")
   >>> emd = make_mt_data_from_edi(
   ...     survey,
   ...     out_path,
   ...     confidence_weighting=True,
   ... )
   >>> print(emd.n_mt_receivers, "receivers,", emd.n_mt_frequencies, "frequencies,", emd.n_data, "data")
   28 receivers, 53 frequencies, 5936 data

The effective relative impedance error is inflated as confidence
decreases:

.. math::

   \epsilon_{Z,\mathrm{eff}} =
   \epsilon_Z
   \left[{1 \over \max(\mathrm{CR}, \mathrm{CR}_{\min})}\right]^p .

The defaults are ``CR_min=0.05`` and ``p=1``. The usual propagation is
then:

.. math::

   \sigma_{\rho_a,\mathrm{eff}} =
   2\rho_a\,\epsilon_{Z,\mathrm{eff}},
   \qquad
   \sigma_{\phi,\mathrm{eff}} =
   {180 \over \pi}\epsilon_{Z,\mathrm{eff}} .

Confidence weighting should increase uncertainty for low-confidence
data. It should not make any datum artificially more precise.

Build A QC Report Bundle
------------------------

The following script writes station tables, frequency tables, and the
main QC figures for one line, reusing the survey already loaded above.

.. code-block:: pycon

   >>> from pathlib import Path
   >>> import matplotlib.pyplot as plt
   >>> from pycsamt.emtools import (
   ...     build_qc_table,
   ...     frequency_confidence_table,
   ...     plot_confidence_band_summary,
   ...     plot_confidence_profile,
   ...     plot_frequency_confidence_psection,
   ...     qc_flags,
   ...     station_confidence_table,
   ... )
   >>> from pycsamt.emtools.qc import plot_qc_quicklook
   >>> out = Path("qc_report_l18plt")
   >>> out.mkdir(parents=True, exist_ok=True)
   >>> build_qc_table(survey, api=False).to_csv(out / "station_qc_summary.csv", index=False)
   >>> qc_flags(survey).to_csv(out / "station_qc_flags.csv", index=False)
   >>> station_confidence_table(survey, method="composite", api=False).to_csv(
   ...     out / "station_confidence.csv", index=False,
   ... )
   >>> frequency_confidence_table(survey, method="composite", api=False).to_csv(
   ...     out / "frequency_confidence.csv", index=False,
   ... )
   >>> fig, ax = plt.subplots(figsize=(9.5, 4.2))
   >>> _ = plot_confidence_profile(survey, method="composite", ax=ax)
   >>> fig.tight_layout()
   >>> fig.savefig(out / "confidence_profile.png", dpi=200)
   >>> plt.close(fig)
   >>> fig, ax = plt.subplots(figsize=(10.0, 4.8))
   >>> _ = plot_frequency_confidence_psection(survey, method="composite", ax=ax)
   >>> fig.savefig(out / "frequency_confidence_psection.png", dpi=200)
   >>> plt.close(fig)
   >>> fig, ax = plt.subplots(figsize=(8.5, 4.2))
   >>> _ = plot_confidence_band_summary(survey, method="composite", ax=ax)
   >>> fig.savefig(out / "confidence_band_summary.png", dpi=200)
   >>> plt.close(fig)
   >>> fig = plot_qc_quicklook(survey)
   >>> fig.savefig(out / "qc_quicklook.png", dpi=200)
   >>> plt.close(fig)
   >>> sorted(p.name for p in out.iterdir())
   ['confidence_band_summary.png', 'confidence_profile.png', 'frequency_confidence.csv', 'frequency_confidence_psection.png', 'qc_quicklook.png', 'station_confidence.csv', 'station_qc_flags.csv', 'station_qc_summary.csv']

.. grid:: 1 1 2 2
   :gutter: 2

   .. grid-item::

      .. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-19-01.png
         :width: 100%

   .. grid-item::

      .. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-19-02.png
         :width: 100%

   .. grid-item::

      .. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-19-03.png
         :width: 100%

   .. grid-item::

      .. image:: ../../images/user_guide/emtools/user-guide-emtools-qc-19-04.png
         :width: 100%

Eight files land in the output directory -- four tables, four figures --
in the same order the calls were written. A directory listing like the
one above is a cheap sanity check that a batch script actually produced
everything it claims to before a report gets built from it.

Reading QC Results
------------------

Use QC scores as evidence, not as an automatic delete button.

``coverage`` is low
    The station or frequency is genuinely incomplete. Investigate
    loading, editing, or acquisition gaps.

``uncertainty`` is low
    Error tensors are large relative to impedance magnitude. This is a
    strong reason to down-weight or review.

``offdiag`` is low
    ``Zxy`` and ``Zyx`` amplitudes disagree beyond the selected
    tolerance. Compare with anisotropy, impedance, and tensor tools.

``diagonal`` is low
    Diagonal leakage is high. Check dimensionality and coordinate
    orientation before calling it acquisition noise.

``phase`` is low
    Off-diagonal phase changes abruptly with frequency. Check frequency
    editing and processing history.

``spatial`` is low
    The station differs from neighbouring stations at the same frequency
    or in median response. Check station metadata and local geology.

Worked Example
--------------

The gallery example applies the station, frequency, profile,
dashboard, quicklook, fan, crossover, and hole-overlay workflows to the
bundled L18PLT survey.

Open the rendered gallery page here:
:ref:`sphx_glr_examples_emtools_plot_qc.py`.
