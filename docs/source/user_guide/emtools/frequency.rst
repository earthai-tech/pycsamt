.. _emtools_frequency:

Frequency Editing, Resampling, And QC
=====================================

``pycsamt.emtools.frequency`` is the user-guide hub for operations that
change, evaluate, or visualize the frequency axis. It is used before
inversion, before section plotting, after :term:`quality control`, and
whenever multiple stations need to share a consistent frequency grid.

The module covers four practical jobs:

- selecting a usable frequency band;
- dropping, masking, or recovering low-confidence rows;
- regridding, decimating, smoothing, and aligning station grids;
- plotting frequency coverage, apparent depth, and band-level summary
  diagnostics.

Full callable signatures live in the :doc:`API reference <../../api/emtools>`.
This page focuses on how to use the workflows safely, how to read their
tables, and how to avoid common frequency-axis traps.

Why Frequency Editing Matters
-----------------------------

CSAMT/AMT/MT processing is rarely limited by one bad station value. More
often, a workflow succeeds or fails because the frequency axis is
inconsistent, too noisy at the band edges, or incompatible across
stations. A 2-D inversion, pseudo-section, or model-comparison workflow
needs clear decisions about:

- which band is physically meaningful;
- which frequency rows are trusted;
- whether bad rows are removed or kept as missing values;
- whether stations should be interpolated to a shared grid;
- whether smoothing or decimation is acceptable for the next step.

The frequency module makes those decisions explicit and auditable.

Data Contract
-------------

All editing functions accept the usual ``emtools`` input:

- a directory containing EDI files,
- one EDI-like object,
- a ``Sites`` container,
- an iterable of site-like objects.

Internally the functions call ``ensure_sites``. The familiar
``recursive``, ``on_dup``, ``strict``, and ``verbose`` arguments behave
the same way as in the rest of the user guide.

Most editing functions default to ``inplace=False``. Keep that default
while exploring. Switch to ``inplace=True`` only when you intentionally
want to mutate the input object.

Selecting A Band
----------------

Use ``select_band`` when you already know the frequency limits.

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.emtools.frequency import select_band
   >>> edi_dir = Path("data/AMT/WILLY_DATA/L18PLT")
   >>> mid_band = select_band(
   ...     edi_dir,
   ...     fmin=10.0,
   ...     fmax=1000.0,
   ...     inplace=False,
   ... )

``fmin`` and ``fmax`` are in hertz. The convenience argument
``band_hz=(10.0, 1000.0)`` is also accepted. If both forms are provided,
explicit ``fmin`` and ``fmax`` take precedence.

Use ``keep`` when you want specific frequencies rather than a continuous
band.

.. code-block:: pycon

   >>> selected = select_band(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     keep=[10.0, 100.0, 1000.0],
   ... )

Band selection removes rows outside the chosen band. It changes the
frequency grid length.

Removing Duplicate Frequencies
------------------------------

Use ``drop_duplicates`` when a station has repeated or unsorted
frequency rows.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import drop_duplicates
   >>> cleaned = drop_duplicates(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     tol=1e-10,
   ...     inplace=False,
   ... )

The function keeps the first occurrence of each frequency within the
tolerance and applies the same row selection to impedance and tipper
blocks when present.

Frequency Confidence
--------------------

Confidence-based editing uses ``frequency_confidence_table`` from the
QC workflow (:ref:`emtools_qc`). Every station-frequency sample is
scored with the same composite confidence ratio defined there, which
folds six bounded diagnostics into one number:

.. math::

   \mathrm{CR}_{i,f} =
   { \sum_k w_k\, s_{k,i,f}\, \mathbf{1}_{s_{k,i,f}\ \mathrm{finite}}
     \over
     \sum_k w_k\, \mathbf{1}_{s_{k,i,f}\ \mathrm{finite}} },
   \qquad 0 \le s_k \le 1,

where :math:`i` indexes the station, :math:`f` the frequency, and
:math:`k` runs over ``coverage``, ``uncertainty``, ``offdiag``,
``diagonal``, ``phase``, and ``spatial``, weighted respectively
``0.35, 0.20, 0.15, 0.10, 0.10, 0.10``. A missing component drops out of
both the numerator and the denominator, so a station with no
:term:`tipper` is not penalized twice for the same gap. ``ci_hi`` and ``ci_lo`` then cut
the ``[0, 1]`` range of :math:`\mathrm{CR}` into three working classes:
at or above ``ci_hi`` a row is trusted and kept as-is, between ``ci_lo``
and ``ci_hi`` it is a candidate for interpolation, and below ``ci_lo``
it is a rejection candidate. The frequency module computes this table
internally, but you should look at its distribution before editing
anything.

.. code-block:: pycon

   >>> import matplotlib.pyplot as plt
   >>> from pycsamt.emtools.qc import frequency_confidence_table
   >>> survey = "data/AMT/WILLY_DATA/L18PLT"
   >>> confidence = frequency_confidence_table(survey, api=True).to_pandas()
   >>> confidence[["station", "frequency_hz", "confidence", "flags"]].head()
     station  ...                                              flags
   0  18-001A  ...  recoverable,high_error,offdiag_mismatch,diagon...
   1  18-001A  ...  reject,high_error,offdiag_mismatch,diagonal_le...
   2  18-001A  ...  reject,high_error,offdiag_mismatch,diagonal_le...
   3  18-001A  ...  reject,high_error,offdiag_mismatch,diagonal_le...
   4  18-001A  ...  reject,high_error,offdiag_mismatch,diagonal_le...
   [5 rows x 4 columns]
   >>> confidence["confidence"].describe()
   count    1484.000000
   mean        0.659129
   std         0.096056
   min         0.417464
   25%         0.587712
   50%         0.656305
   75%         0.737060
   max         0.863041
   Name: confidence, dtype: float64
   >>> station = "18-001A"
   >>> one = confidence.loc[confidence["station"] == station].sort_values("period_s")
   >>> fig, ax = plt.subplots(figsize=(7, 4.5))
   >>> _ = ax.semilogx(one["period_s"], one["confidence"], "o-")
   >>> _ = ax.axhline(0.90, color="tab:red", linestyle="--", label="ci_hi")
   >>> _ = ax.axhline(0.50, color="tab:orange", linestyle="--", label="ci_lo")
   >>> _ = ax.set_xlabel("Period (s)")
   >>> _ = ax.set_ylabel("Confidence")
   >>> _ = ax.legend()
   >>> fig.tight_layout()

.. image:: ../../images/user_guide/emtools/user-guide-emtools-frequency-04.png
   :width: 100%

The key lesson is simple: do not choose ``ci_hi``, ``ci_lo``, or
``threshold`` before looking at the confidence distribution above. In
this survey the median :math:`\mathrm{CR}` sits near ``0.66`` and the
maximum barely reaches ``0.86``. Leaving ``ci_hi`` at a textbook value
like ``0.90`` would mean no row is ever trusted enough to act as an
interpolation donor, and recovery would silently have nothing to work
from. Read the ``describe()`` output first, then set ``ci_hi`` and
``ci_lo`` (or ``threshold``) to values this particular survey can
actually reach.

Drop Or Mask Low-Confidence Rows
--------------------------------

Once a threshold is chosen, both functions apply the same keep rule to
row :math:`(i, f)`:

.. math::

   \text{keep}_{i,f} =
   \bigl(\mathrm{CR}_{i,f}\ \text{undefined}\bigr)
   \ \lor\
   \bigl(\mathrm{CR}_{i,f} \ge \text{threshold}\bigr).

Rows with no confidence score at all (for example, a frequency present
in the impedance block but absent from the confidence table) are kept
rather than silently discarded — undefined is not the same as bad. What
differs between the two functions is what happens to the rows that fail
the rule. Dropping removes them from the frequency grid entirely.
Masking keeps the grid the same length but replaces the corresponding
tensor rows with ``NaN``.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import (
   ...     drop_low_confidence_frequencies,
   ...     mask_low_confidence_frequencies,
   ... )
   >>> dropped = drop_low_confidence_frequencies(
   ...     survey,
   ...     threshold=0.50,
   ...     also="both",
   ...     inplace=False,
   ... )
   >>> masked = mask_low_confidence_frequencies(
   ...     survey,
   ...     threshold=0.50,
   ...     also="both",
   ...     inplace=False,
   ... )

Use ``drop`` when downstream code cannot handle missing rows. Use
``mask`` when preserving the original frequency grid matters, for
example when comparing before and after processing. Dropping shrinks
the impedance array before the frequency array, in that order, to
avoid a hard failure in the underlying setter; you may see one harmless
``ERROR``-level log line per edited station as a result ("Failed to
compute rho/phi after setting Z") -- it self-corrects one line later
once the frequency array is updated to match, and does not indicate
that anything actually failed.

The ``also`` argument controls which blocks are edited:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Value
     - Meaning
   * - ``"z"``
     - Edit impedance only.
   * - ``"tipper"``
     - Edit tipper only.
   * - ``"both"``
     - Edit impedance and tipper when present.

Recover Low-Confidence Rows
---------------------------

``recover_low_confidence_frequencies`` splits a station's rows into
three sets from the same :math:`\mathrm{CR}_{i,f}` used above:

.. math::

   G = \{f : \mathrm{CR}_{i,f} \ge \text{ci\_hi}\}, \qquad
   R = \{f : \text{ci\_lo} \le \mathrm{CR}_{i,f} < \text{ci\_hi}\}, \qquad
   B = \{f : \mathrm{CR}_{i,f} < \text{ci\_lo}\},

the trusted donors :math:`G`, the recoverable rows :math:`R`, and the
rejected rows :math:`B`. Every row in :math:`R` is filled by
interpolating the complex tensor value at that frequency from the rows
in :math:`G`, in log-frequency space and independently for the real and
imaginary parts:

.. math::

   Z(f_i) = \mathrm{interp}\bigl(\log_{10} f_i;\
   \{\log_{10} f_g\}_{f_g \in G},\ \{Z(f_g)\}_{f_g \in G}\bigr),
   \qquad f_i \in R,

with ``interpolation="linear"`` performing ordinary linear interpolation
in :math:`\log_{10} f` and ``"nearest"`` copying the closest donor
instead. Working in :math:`\log_{10} f` rather than raw hertz matters
because CSAMT/AMT bands are sampled logarithmically — linear
interpolation in hertz would be dominated by the widely spaced
high-frequency end and would barely constrain the low-frequency end at
all. Rows in :math:`B` are not interpolated; they are handled entirely
by ``reject``.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import recover_low_confidence_frequencies
   >>> recovered = recover_low_confidence_frequencies(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     ci_hi=0.72,
   ...     ci_lo=0.50,
   ...     interpolation="linear",
   ...     reject="drop",
   ...     also="both",
   ...     inplace=False,
   ... )

Valid interpolation modes are ``"linear"`` and ``"nearest"``. Valid
``reject`` modes are:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Reject mode
     - Meaning
   * - ``"drop"``
     - Remove rows below ``ci_lo``.
   * - ``"mask"``
     - Keep rows below ``ci_lo`` but set them to ``NaN``.
   * - ``"keep"``
     - Leave rejected rows unchanged.

Recovery is interpolation, not magic: it can only be as good as the
donors in :math:`G`. When a recoverable row sits outside the span of
:math:`\log_{10} f_g` covered by the donors, ``interp`` clamps to the
nearest edge, so "linear" quietly degrades to "nearest" exactly at the
band edges where recovery is riskiest. Always inspect the result when
many rows are recovered, and treat edge-of-band recoveries with more
suspicion than interior ones.

High-Level Editing Workflow
---------------------------

Use ``edit_frequencies_by_confidence`` when you want the edited sites,
a station-level report, and a row-level decision table in one object.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import edit_frequencies_by_confidence
   >>> result = edit_frequencies_by_confidence(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     mode="recover",
   ...     ci_hi=0.72,
   ...     ci_lo=0.50,
   ...     reject="drop",
   ...     interpolation="linear",
   ...     api=False,
   ... )
   >>> result.summary()
   "FrequencyEditResult(mode='recover', method='composite', dropped=68, masked=102, recovered=867)"
   >>> result.report.head()
      station  n_freq_before  ...  n_masked_or_unfinite  confidence_delta
   0  18-001A             53  ...                     0          0.027413
   1  18-002U             53  ...                     0         -0.038885
   2  18-003A             53  ...                     3          0.061644
   3  18-004A             53  ...                     3         -0.026442
   4  18-005U             53  ...                     1          0.026781
   [5 rows x 18 columns]
   >>> result.decisions.head()
      station  frequency_hz  period_s  ...  present_after  finite_after action
   0  18-001A       10400.0  0.000096  ...           True          True   kept
   1  18-001A        8707.0  0.000115  ...           True          True   kept
   2  18-001A        7289.0  0.000137  ...           True          True   kept
   3  18-001A        6102.0  0.000164  ...           True          True   kept
   4  18-001A        5108.0  0.000196  ...           True          True   kept
   [5 rows x 10 columns]

``result`` is a ``FrequencyEditResult`` with:

- ``sites``: edited sites.
- ``report``: station-level before/after table.
- ``decisions``: one row per original station-frequency row.
- ``n_dropped``, ``n_masked``, ``n_recovered``: convenience counts.
- ``mode``, ``method``, ``ci_hi``, ``ci_lo``, ``reject``,
  ``interpolation``: edit metadata.

For in-memory objects, pass an independent ``before_sites`` when you
need an untouched baseline for reporting. Some lower-level editors may
mutate wrapped impedance objects while constructing edited outputs.

Station-Level Report
--------------------

``frequency_edit_report`` compares before and after sites.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import (
   ...     edit_frequencies_by_confidence,
   ...     frequency_edit_report,
   ... )
   >>> survey = "data/AMT/WILLY_DATA/L18PLT"
   >>> result = edit_frequencies_by_confidence(
   ...     survey,
   ...     mode="drop",
   ...     threshold=0.50,
   ...     api=False,
   ... )
   >>> report = frequency_edit_report(
   ...     survey,
   ...     result.sites,
   ...     ci_hi=0.90,
   ...     ci_lo=0.50,
   ...     api=True,
   ... ).to_pandas()
   >>> report[[
   ...     "station", "n_freq_before", "n_freq_after", "n_dropped",
   ...     "n_finite_before", "n_finite_after",
   ...     "confidence_median_before", "confidence_median_after",
   ... ]].head()
      station  n_freq_before  ...  confidence_median_before  confidence_median_after
   0  18-001A             53  ...                  0.711753                 0.711753
   1  18-002U             53  ...                  0.749480                 0.749480
   2  18-003A             53  ...                  0.666613                 0.671639
   3  18-004A             53  ...                  0.735994                 0.743510
   4  18-005U             53  ...                  0.728841                 0.733250
   [5 rows x 8 columns]

Important report columns include:

- ``n_freq_before`` and ``n_freq_after``: grid size change.
- ``n_finite_before`` and ``n_finite_after``: finite tensor row counts.
- ``frac_finite_before`` and ``frac_finite_after``: finite row fraction.
- ``safe_fraction_*``: fraction at or above ``ci_hi``.
- ``recoverable_fraction_*``: fraction in ``[ci_lo, ci_hi)``.
- ``reject_fraction_*``: fraction below ``ci_lo``.
- ``n_dropped``: number of removed frequency rows.
- ``n_masked_or_unfinite``: finite rows lost to masking or missing data.
- ``confidence_delta``: median confidence change.

Decision Table
--------------

``frequency_edit_decision_table`` records one row per original
station-frequency sample.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import (
   ...     edit_frequencies_by_confidence,
   ...     frequency_edit_decision_table,
   ... )
   >>> survey = "data/AMT/WILLY_DATA/L18PLT"
   >>> result = edit_frequencies_by_confidence(
   ...     survey,
   ...     mode="recover",
   ...     ci_hi=0.72,
   ...     ci_lo=0.50,
   ...     reject="drop",
   ...     api=False,
   ... )
   >>> decisions = frequency_edit_decision_table(
   ...     survey,
   ...     result.sites,
   ...     ci_hi=0.72,
   ...     ci_lo=0.50,
   ...     api=True,
   ... ).to_pandas()
   >>> decisions["action"].value_counts()
   action
   recovered    867
   kept         447
   masked       102
   dropped       68
   Name: count, dtype: int64

Actions are:

- ``kept``: row survived unchanged.
- ``dropped``: row is absent after editing.
- ``masked``: row remains but is no longer finite.
- ``recovered``: row was filled or changed by recovery.

Plot Edit Results
-----------------

Use the summary plot for station-level effects and the decision plot for
row-level actions.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import (
   ...     plot_frequency_edit_decisions,
   ...     plot_frequency_edit_summary,
   ... )
   >>> fig, (ax_summary, ax_decisions) = plt.subplots(2, 1, figsize=(10, 8))
   >>> _ = plot_frequency_edit_summary(
   ...     survey, result.sites, ci_hi=0.72, ci_lo=0.50, ax=ax_summary,
   ... )
   >>> _ = plot_frequency_edit_decisions(
   ...     survey, result.sites, ci_hi=0.72, ci_lo=0.50, ax=ax_decisions,
   ... )
   >>> fig.tight_layout()

.. image:: ../../images/user_guide/emtools/user-guide-emtools-frequency-10.png
   :width: 100%

If the decision plot is dominated by ``recovered`` or ``masked`` colors,
the edit is not conservative and should be justified.

Regrid To A Target Grid
-----------------------

Use ``regrid_to`` when you already have the target frequencies.

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.emtools.frequency import regrid_to
   >>> target_freq = np.array([1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0])
   >>> regridded = regrid_to(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     target_freq,
   ...     method="linear",
   ...     inplace=False,
   ... )

``regrid_to`` reuses the same log-frequency interpolator introduced
above for recovery, except every target frequency is filled from the
*entire* native grid rather than from a confidence-selected subset of
it — there is no trusted/recoverable split here, only the target grid
you hand it. ``method="linear"`` interpolates the real and imaginary
parts in :math:`\log_{10} f`; ``method="nearest"`` copies the closest
native row instead of interpolating. Choose a target grid that is
meaningful for the inversion or comparison you are preparing for, not
merely one that is convenient to type.

Build A Log-Spaced Grid
-----------------------

Typing out a target grid by hand gets tedious once a survey needs a
regular per-decade sampling. ``regrid_logspace`` builds one automatically
by geometric progression:

.. math::

   f_k = f_{\min} \left(\frac{f_{\max}}{f_{\min}}\right)^{k / (n - 1)},
   \qquad k = 0, \dots, n - 1,
   \qquad n = \max\!\left(2,\ \left\lceil
   \log_{10}(f_{\max} / f_{\min}) \times \text{per\_decade}
   \right\rceil\right).

Equal steps in :math:`k` are equal steps in :math:`\log_{10} f`, so
``per_decade`` directly controls resolution: ``per_decade=6`` places six
points in every factor-of-ten span, matching the density a typical
CSAMT/AMT band is actually sampled at. When ``fmin``/``fmax`` are
omitted, they default to the min/max of the union grid across all
stations.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import regrid_logspace
   >>> regular = regrid_logspace(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     fmin=1.0,
   ...     fmax=10000.0,
   ...     per_decade=6,
   ...     method="linear",
   ...     inplace=False,
   ... )

``band_hz=(lo, hi)`` is accepted as an alias for ``fmin`` and ``fmax``.
``n_per_decade`` is accepted as an alias for ``per_decade`` when
``per_decade`` is left at its default.

Decimation And Moving Average Smoothing
---------------------------------------

``decimate_step`` keeps rows at indices :math:`0, \text{step}, 2\cdot
\text{step}, \dots` along the native (not log-spaced) frequency array —
it thins the grid but does not touch the tensor values that survive.
``smooth_mavg`` does the opposite: it keeps every row but replaces each
tensor value with a centered moving average of width :math:`k`,
computed independently on the real and imaginary parts,

.. math::

   Z'(f_i) = \frac{1}{k} \sum_{j=-\lfloor k/2 \rfloor}^{\lceil k/2
   \rceil - 1} Z(f_{i+j}),

so it is a smoothing filter over row order, not over :math:`\log_{10}
f`, and it changes every value it touches rather than removing any.
Near the two ends of the grid the window runs past the available rows;
the implementation still divides by the full :math:`k` rather than
shrinking the window, which pulls the first and last :math:`\lfloor
k/2\rfloor` smoothed values slightly toward zero. Keep ``k`` small
relative to the number of rows in a station's band if the edge rows
still need to be trusted afterward.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import decimate_step, smooth_mavg
   >>> decimated = decimate_step(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     step=2,
   ...     inplace=False,
   ... )
   >>> smoothed = smooth_mavg(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     k=3,
   ...     on="z",
   ...     inplace=False,
   ... )

Smoothing changes the measured tensor values. Use it as a processing
choice, not as a plotting convenience, and report the window size.

Align Station Grids
-------------------

Use ``align_grid`` when stations need a common grid.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import align_grid
   >>> union_aligned = align_grid(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     mode="union",
   ...     method="nearest",
   ... )
   >>> intersection_aligned = align_grid(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     mode="intersection",
   ...     method="nearest",
   ... )
   >>> ref_aligned = align_grid(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     mode="ref",
   ...     ref_station="18-001A",
   ...     method="nearest",
   ... )

The modes are:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Mode
     - Meaning
   * - ``union``
     - Use all frequencies present in any station.
   * - ``intersection``
     - Use only frequencies present in every station.
   * - ``ref``
     - Use the grid from ``ref_station``.

``intersection`` can be empty or nearly empty for independently
processed data because frequencies may not match bit-for-bit. If the
intersection is empty, the function returns the input sites unchanged.

Coverage And Quality Heatmap
----------------------------

``plot_coverage_quality_heatmap`` shows which station-frequency cells
exist and how reliable they are according to relative impedance error.
This is a separate, lighter-weight score than the composite
:math:`\mathrm{CR}` used above — it looks only at the off-diagonal
impedance error, so it stays cheap enough to compute for a full survey
heatmap. For each cell it averages the relative error of the two
off-diagonal modes and maps that to a quality in ``(0, 1]``:

.. math::

   r_{i,f} = \frac{1}{2}\left(
   \frac{|\sigma_{Z_{xy}}|}{|Z_{xy}|} +
   \frac{|\sigma_{Z_{yx}}|}{|Z_{yx}|}
   \right),
   \qquad
   q_{i,f} = \frac{1}{1 + r_{i,f}},

where :math:`\sigma_Z` is the reported impedance error at that cell. A
noise-free measurement has :math:`r = 0` and :math:`q = 1`; as the error
grows relative to the signal, :math:`q` decays smoothly toward ``0``
instead of blowing up, which keeps the colour scale usable even for
badly noisy cells.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import plot_coverage_quality_heatmap
   >>> fig, ax = plt.subplots(figsize=(9, 4.5))
   >>> _ = plot_coverage_quality_heatmap(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     axis="period",
   ...     ax=ax,
   ... )
   >>> fig.tight_layout()

.. image:: ../../images/user_guide/emtools/user-guide-emtools-frequency-15.png
   :width: 100%

Values closer to ``1`` are higher quality; missing cells simply mean
that station has no row at that union-grid frequency, which is a
coverage gap rather than a quality problem.

Apparent Depth Pseudo-Section
-----------------------------

``plot_apparent_depth_psection`` turns each station-frequency cell into
a rough sense of investigation depth. It first needs an :term:`apparent
resistivity`, which it gets from the determinant of the two off-diagonal
:term:`impedance tensor` modes — the same geometric-mean convention used
elsewhere in ``emtools`` (see :ref:`emtools_csumt`):

.. math::

   \rho_{a,\det} =
   \sqrt{
   \left(0.2\,\frac{|Z_{xy}|^2}{f}\right)
   \left(0.2\,\frac{|Z_{yx}|^2}{f}\right)
   },

and then converts that resistivity and frequency into a
:term:`skin depth`-style apparent depth,

.. math::

   \delta \approx 503 \sqrt{\rho_{a,\det} / f},

with :math:`\delta` in metres when :math:`\rho_{a,\det}` is in
:math:`\Omega\,\mathrm{m}` and :math:`f` is in hertz. Lower frequencies
and higher apparent resistivity both push :math:`\delta` deeper, which
is exactly the trade-off a CSAMT/AMT sounding is designed to exploit.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import plot_apparent_depth_psection
   >>> fig, ax = plt.subplots(figsize=(9, 4.5))
   >>> _ = plot_apparent_depth_psection(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     axis_y="period",
   ...     agg="median",
   ...     log_color=True,
   ...     ax=ax,
   ... )
   >>> fig.tight_layout()

.. image:: ../../images/user_guide/emtools/user-guide-emtools-frequency-16.png
   :width: 100%

This is a diagnostic visualization, not an inversion result — :math:`\delta`
is a proxy, not a layer boundary. Use it to see how apparent depth
varies across stations and periods before committing to a model
parameterization.

Band Microstrips
----------------

A full pseudo-section is the right tool for inspecting one survey
closely, but it stops being readable once you want to compare many
lines side by side. ``plot_band_microstrips`` collapses the period axis
into a handful of bands and reduces each station-band cell to three
median summary metrics:

.. math::

   \log_{10}\rho_{a,\det}, \qquad
   \varphi_{\det} = \arg\bigl[\det(Z)\bigr]
   = \arg\bigl(Z_{xx}Z_{yy} - Z_{xy}Z_{yx}\bigr), \qquad
   |T| = \sqrt{|T_x|^2 + |T_y|^2},

the log-apparent-resistivity defined above, the determinant phase (in
degrees, using the same :math:`\det(Z)` convention as the tensor pages),
and the tipper magnitude when a station has one. Each metric gets its
own marker — circle, square, triangle — and its own colour
normalization, so a single compact grid of dots can stand in for what
would otherwise be three separate pseudo-sections per line.

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import plot_band_microstrips
   >>> bands = [
   ...     (1e-4, 1e-3),
   ...     (1e-3, 1e-2),
   ...     (1e-2, 1e-1),
   ...     (1e-1, 1.0),
   ... ]
   >>> fig, ax = plt.subplots(figsize=(9, 6))
   >>> _ = plot_band_microstrips(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     bands=bands,
   ...     ax=ax,
   ... )
   >>> fig.tight_layout()

.. image:: ../../images/user_guide/emtools/user-guide-emtools-frequency-17.png
   :width: 100%

Use microstrips for compact survey summaries. They are not a substitute
for full curves when a band contains strong internal variation.

Recommended Processing Pattern
------------------------------

A conservative frequency workflow looks like this:

1. Inspect the native frequency range and confidence table.
2. Select a physically meaningful band.
3. Decide whether bad rows should be dropped, masked, or recovered.
4. Save the edit report and decision table.
5. Regrid or align only when the next workflow requires it.
6. Plot coverage and apparent depth after editing.
7. Keep the original data available for audit.

.. code-block:: pycon

   >>> from pathlib import Path
   >>> out = Path("outputs/frequency_l18plt")
   >>> out.mkdir(parents=True, exist_ok=True)
   >>> banded = select_band(survey, fmin=10.0, fmax=1000.0)
   >>> result = edit_frequencies_by_confidence(
   ...     banded,
   ...     mode="recover",
   ...     before_sites=banded,
   ...     ci_hi=0.72,
   ...     ci_lo=0.50,
   ...     reject="drop",
   ...     api=False,
   ... )
   >>> regular = regrid_logspace(result.sites, fmin=10.0, fmax=1000.0, per_decade=6)
   >>> result.report.to_csv(out / "frequency_edit_report.csv", index=False)
   >>> result.decisions.to_csv(out / "frequency_edit_decisions.csv", index=False)
   >>> fig1, ax1 = plt.subplots(figsize=(9, 4.5))
   >>> _ = plot_coverage_quality_heatmap(regular, ax=ax1)
   >>> fig1.savefig(out / "coverage_quality_heatmap.png", dpi=200)
   >>> fig2, ax2 = plt.subplots(figsize=(9, 4.5))
   >>> _ = plot_apparent_depth_psection(regular, ax=ax2)
   >>> fig2.savefig(out / "apparent_depth_psection.png", dpi=200)

.. grid:: 1 1 2 2
   :gutter: 2

   .. grid-item::

      .. image:: ../../images/user_guide/emtools/user-guide-emtools-frequency-18-01.png
         :width: 100%

   .. grid-item::

      .. image:: ../../images/user_guide/emtools/user-guide-emtools-frequency-18-02.png
         :width: 100%

Common Failure Modes
--------------------

Using default ``ci_hi`` without checking confidence
   Recovery needs trusted rows at or above ``ci_hi``. If the survey
   confidence never reaches that value, recovery cannot behave as
   expected.

Dropping rows before checking downstream needs
   Dropping shortens the grid. If the next step expects aligned grids,
   masking may be safer.

Masking rows before plotting finite-only workflows
   Some workflows ignore or fail on ``NaN`` rows. Save a decision table
   so missing values are explainable.

Regridding too early
   Interpolation can hide native acquisition problems. Inspect native
   confidence first.

Empty intersection alignment
   Independently processed stations may have no exact shared frequency.
   Use ``union`` or a reference station when appropriate.

Smoothing as hidden preprocessing
   Moving averages change tensor values. Report the window and apply it
   intentionally.

Worked Example
--------------

The gallery example uses **L18PLT** from ``data/AMT/WILLY_DATA/`` and
brings in **KAP03** for heterogeneous-grid alignment. It demonstrates
band selection, confidence editing, recovery, reports, decision plots,
regridding, decimation, smoothing, grid alignment, coverage heatmaps,
apparent-depth pseudo-sections, and band microstrips.

Open the rendered example here:
:ref:`sphx_glr_examples_emtools_plot_frequency.py`.
