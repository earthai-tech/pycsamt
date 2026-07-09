.. _emtools_ss:

Static-Shift Correction
=======================

Static shift is a station-dependent, frequency-independent offset in
apparent resistivity.  In CSAMT, AMT, and MT data it is usually caused
by shallow galvanic distortion: the curve keeps nearly the same shape,
but it is lifted or lowered by a multiplier.  Because apparent
resistivity is proportional to ``|Z|**2``, pyCSAMT estimates the shift in
``log10(rho)`` and applies the square-root correction to the impedance
tensor ``Z``.

The static-shift tools in ``pycsamt.emtools`` cover four related jobs:

.. list-table::
   :header-rows: 1
   :widths: 28 32 40

   * - Job
     - Main tools
     - Use when
   * - Estimate correction factors
     - ``estimate_ss_ama``, ``estimate_ss_loess``,
       ``estimate_ss_bilateral``, ``estimate_ss_refmedian``
     - You want a per-station table of static-shift amplitudes before
       modifying data.
   * - Apply factors
     - ``apply_ss_factors``, ``correct_ss_ama``
     - You have accepted a factor table and want corrected ``Sites``.
   * - Check the correction visually
     - ``plot_ss_delta_psection``, ``plot_ss_station_curves``,
       ``plot_ss_delta_profile``, ``ss_qc_*``,
       ``ss_comparison_psection``
     - You need before/after diagnostics rather than only a table.
   * - Separate static shift from near-surface effects
     - ``detect_near_surface``, ``plot_ns_detection``
     - The distortion may be frequency-dependent and therefore not
       removable by a conventional static-shift multiplier.

The examples below use two-level imports from ``pycsamt.emtools`` because
the static-shift functions are exported as public user-facing symbols.

Load The Survey
---------------

Start from the canonical loader.  It accepts a directory, a glob-like
collection, or an existing ``Sites`` object and skips files without valid
impedance data.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.emtools import ensure_sites

   edi_dir = Path("data/AMT/WILLY_DATA/L18PLT")
   sites = ensure_sites(edi_dir, recursive=True, verbose=0)

For static-shift work, station order matters.  The AMA method estimates a
local spatial trend from neighbouring stations, so choose ``sort_by`` to
match the survey line direction:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - ``sort_by``
     - Meaning
   * - ``"lon"``
     - Order stations by longitude.  This is useful for an east-west line.
   * - ``"lat"``
     - Order stations by latitude.  This is useful for a north-south line.
   * - ``"name"``
     - Order stations lexically by station name.  Use this only when the
       station names encode the real line order.

Estimate AMA Factors
--------------------

``estimate_ss_ama`` is the main estimator.  AMA means adaptive
moving-average: for each station, pyCSAMT builds a weighted local
reference from neighbouring stations, compares the target station with
that reference in ``log10(rho_det)``, and reduces the per-frequency
residuals to one static-shift value.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import estimate_ss_ama

   factors = estimate_ss_ama(
       sites,
       sort_by="lat",
       half_window=3,
       weights="tri",
       pband=(1e-4, 10.0),
       max_skew=45.0,
       robust_freq="median",
       robust_overall="median",
   )

   print(factors[["station", "delta_log10_rho", "fac_rho", "fac_z", "n_used"]])

The returned table has one row per station with these columns:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Column
     - Interpretation
   * - ``station``
     - Station identifier used to match the factor back to the survey.
   * - ``delta_log10_rho``
     - Estimated vertical offset of the station relative to the local
       spatial trend.  Positive means the station's apparent resistivity
       is above the trend.
   * - ``fac_rho``
     - Resistivity correction factor, computed as
       ``10 ** (-delta_log10_rho)``.
   * - ``fac_z``
     - Impedance correction factor, computed as
       ``10 ** (-0.5 * delta_log10_rho)``.  This is the column normally
       used by ``apply_ss_factors``.
   * - ``n_used``
     - Number of finite frequency samples that contributed to the station
       estimate after period-band and skew filtering.

The phase-tensor skew filter is important.  The default threshold is
conservative, so a structurally complex survey may return a sparse table
unless you choose a threshold appropriate for the line.  A quick audit is
usually better than accepting the first number silently:

.. code-block:: python
   :linenos:

   strict_factors = estimate_ss_ama(
       sites,
       sort_by="lat",
       half_window=3,
       max_skew=6.0,
   )

   survey_factors = estimate_ss_ama(
       sites,
       sort_by="lat",
       half_window=3,
       max_skew=45.0,
   )

   print("strict stations:", len(strict_factors))
   print("survey stations:", len(survey_factors))
   print("median samples:", survey_factors["n_used"].median())

If ``n_used`` is low for many stations, the result may be technically
valid but weakly constrained.  In that case, inspect the skew page,
restrict ``pband`` to the useful period range, or compare several
estimators before applying the correction.

Read Factor Signs Correctly
---------------------------

A positive ``delta_log10_rho`` means the raw apparent resistivity is high
relative to the local trend, so the correction factor is less than one.
A negative ``delta_log10_rho`` means the raw apparent resistivity is low,
so the correction factor is greater than one.

.. code-block:: python
   :linenos:

   view = factors.assign(
       direction=factors["delta_log10_rho"].map(
           lambda value: "lower Z" if value > 0 else "raise Z"
       )
   )

   print(
       view[
           ["station", "delta_log10_rho", "fac_z", "n_used", "direction"]
       ].sort_values("delta_log10_rho")
   )

``fac_z`` must be finite and positive.  ``apply_ss_factors`` guards the
impedance tensor against non-finite, zero, and negative factors by
leaving those stations unchanged, but your processing report should still
flag such rows because they indicate a failed or unreliable estimate.

.. code-block:: python
   :linenos:

   import numpy as np

   good = np.isfinite(factors["fac_z"]) & (factors["fac_z"] > 0.0)
   rejected = factors.loc[~good, ["station", "fac_z", "n_used"]]

   if not rejected.empty:
       print("Rows that will not be applied:")
       print(rejected)

Apply Factors
-------------

Use ``apply_ss_factors`` when you want an explicit estimate-then-apply
workflow.  This is the preferred pattern for production processing
because the table can be saved, reviewed, plotted, and reproduced.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import apply_ss_factors

   corrected = apply_ss_factors(
       sites,
       factors,
       key="fac_z",
       inplace=False,
   )

Keep ``inplace=False`` while developing a processing flow.  It returns a
corrected copy and leaves the original ``Sites`` object available for
before/after plots.  Use ``inplace=True`` only when you deliberately want
to mutate the object already in memory.

For exploratory notebooks, ``correct_ss_ama`` combines estimation and
application in one call:

.. code-block:: python
   :linenos:

   from pycsamt.emtools import correct_ss_ama

   corrected = correct_ss_ama(
       sites,
       sort_by="lat",
       half_window=3,
       max_skew=45.0,
       inplace=False,
   )

This convenience call is useful, but it does not show the factor table by
itself.  When documenting a survey, keep the explicit ``factors`` table
as part of the processing record.

Compare Estimators
------------------

Static-shift estimation is not a single universal recipe.  The four
estimators use different references:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Estimator
     - Reference model
     - Practical reading
   * - ``estimate_ss_ama``
     - Weighted local median/mean from neighbouring stations.
     - Good default for line data when station order is meaningful.
   * - ``estimate_ss_loess``
     - Local polynomial regression across neighbouring stations.
     - Useful when the background level changes smoothly along the line.
   * - ``estimate_ss_bilateral``
     - Neighbour trend weighted by both distance and value similarity.
     - Can protect sharp local contrasts, but may disagree at outliers.
   * - ``estimate_ss_refmedian``
     - Global frequency-wise median reference.
     - Useful as a broad check, less local than AMA or LOESS.

Run several estimators before correcting a valuable line.  Agreement is a
stronger signal than a single factor table.

.. code-block:: python
   :linenos:

   from functools import reduce

   from pycsamt.emtools import (
       estimate_ss_bilateral,
       estimate_ss_loess,
       estimate_ss_refmedian,
   )

   tables = {
       "ama": estimate_ss_ama(
           sites,
           sort_by="lat",
           half_window=3,
           max_skew=45.0,
       ),
       "loess": estimate_ss_loess(
           sites,
           half_window=3,
           max_skew=45.0,
       ),
       "bilateral": estimate_ss_bilateral(
           sites,
           half_window=4,
           max_skew=45.0,
       ),
       "refmedian": estimate_ss_refmedian(
           sites,
           max_skew=45.0,
       ),
   }

   compact = []
   for name, table in tables.items():
       compact.append(
           table[["station", "delta_log10_rho"]].rename(
               columns={"delta_log10_rho": name}
           )
       )

   comparison = reduce(
       lambda left, right: left.merge(right, on="station", how="inner"),
       compact,
   )

   print(comparison)
   print(comparison.drop(columns="station").corr().round(2))

Look for sign disagreements, not only large absolute differences.  If
AMA says a station should be lowered while bilateral filtering says it
should be raised, that station deserves a manual plot before correction.

QC Plots In One Call
--------------------

The ``ss_qc_*`` wrappers estimate, apply, and plot in one call.  They are
good first-look tools when you do not need to reuse the corrected object.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       ss_qc_profile,
       ss_qc_psection,
       ss_qc_station_curves,
   )

   ss_qc_psection(
       sites,
       method="ama",
       sort_by="lat",
       half_window=3,
       max_skew=45.0,
       pband=(1e-4, 10.0),
   )

   ss_qc_profile(
       sites,
       method="ama",
       sort_by="lat",
       half_window=3,
       max_skew=45.0,
       robust="median",
   )

   ss_qc_station_curves(
       sites,
       method="ama",
       station="18-016A",
       sort_by="lat",
       half_window=3,
       max_skew=45.0,
   )

``ss_qc_psection`` shows the pointwise
``log10(rho_after) - log10(rho_before)`` correction as a pseudosection.
For a true static-shift correction, each station tends to show a nearly
vertical band because the multiplier is constant with frequency.

``ss_qc_profile`` compresses that difference to one value per station.
Use it to find stations with unusually large corrections.

``ss_qc_station_curves`` overlays the before and after curves for one
station.  It is the quickest way to confirm that the correction shifted
the level without changing the curve shape.

Before/After Plots From Existing Sites
--------------------------------------

When you already have a corrected ``Sites`` object, call the lower-level
plotters directly.  These functions do not re-estimate the correction;
they compare the two data sets you pass in.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       plot_ss_delta_profile,
       plot_ss_delta_psection,
       plot_ss_station_curves,
   )

   plot_ss_delta_psection(
       sites,
       corrected,
       axis_y="logperiod",
       pband=(1e-4, 10.0),
   )

   plot_ss_delta_profile(
       sites,
       corrected,
       robust="median",
       pband=(1e-4, 10.0),
   )

   plot_ss_station_curves(
       sites,
       corrected,
       station="18-016A",
       pband=(1e-4, 10.0),
   )

These calls are useful in reports because they make the provenance clear:
``sites`` is the uncorrected input and ``corrected`` is the accepted
corrected output.

Publication Comparison Figures
------------------------------

``ss_comparison_psection`` is the high-level publication helper.  It
estimates the correction, applies it, builds aligned before/after
``log10(rho_det)`` arrays internally, and renders a shared-scale
pseudosection comparison.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import ss_comparison_psection

   fig = ss_comparison_psection(
       sites,
       method="ama",
       sort_by="lat",
       half_window=3,
       max_skew=45.0,
       show_delta=True,
       suptitle="Static-shift correction: AMA",
   )

Use ``show_delta=True`` when the figure must show both the corrected
resistivity field and the actual correction amplitude.  The before and
after panels share colour limits, so station-level offsets remain
visible instead of being hidden by automatic rescaling.

The lower-level ``plot_ss_comparison_psection``, ``plot_ss_1d_curves``,
and ``plot_ss_summary`` functions accept precomputed arrays:
``logRho_before`` with shape ``(n_stations, n_frequencies)``,
``logRho_after`` with the same shape, and ``freqs`` in hertz.  Use those
array-level functions when you have already built a custom resistivity
matrix outside the normal ``Sites`` workflow.

Radar View
----------

``plot_ss_radar`` shows the ``xy`` and ``yx`` apparent-resistivity
components for one station on a polar grid.  The angle around the circle
represents period or frequency, and the radius represents resistivity or
``log10(resistivity)``.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_ss_radar

   plot_ss_radar(
       sites,
       station="18-016A",
       rotate="none",
       radial="log10rho",
       theta_axis="logperiod",
       fill_between=True,
   )

   plot_ss_radar(
       sites,
       station="18-016A",
       rotate="pt",
       rotate_stat="median",
       radial="log10rho",
   )

Use ``rotate="none"`` for a raw component view, ``rotate="deg"`` with
``rotate_deg=...`` for a fixed rotation, and ``rotate="pt"`` for
phase-tensor based rotation.  The radar plot is diagnostic: it helps you
see directional imbalance and frequency structure, but it is not itself a
static-shift estimator.

Near-Surface Versus Static Shift
--------------------------------

Conventional static-shift correction assumes the distortion is a constant
multiplier.  Some shallow effects are frequency-dependent: the high
frequency part of the curve may become much noisier or steeper than the
low frequency part.  A single static multiplier cannot fix that behavior.

``detect_near_surface`` classifies each station using three diagnostics:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Diagnostic
     - Meaning
   * - ``ns_index``
     - Ratio of high-frequency residual spread to low-frequency residual
       spread.  Values above ``ns_threshold`` indicate near-surface
       frequency-dependent distortion.
   * - ``gradient_delta``
     - Difference between high-frequency and low-frequency log-log slopes.
       Larger values indicate a change in curve shape.
   * - ``ss_delta_log10``
     - AMA-like static-shift residual.  Values above ``ss_threshold`` are
       consistent with a station-level static offset.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import detect_near_surface, plot_ns_detection

   ns = detect_near_surface(
       sites,
       f_split=1.0,
       ns_threshold=2.0,
       ss_threshold=0.1,
       sort_by="lat",
       half_window=3,
       max_skew=45.0,
   )

   print(ns[[
       "station",
       "ns_index",
       "gradient_delta",
       "ss_delta_log10",
       "distortion_type",
   ]])
   print(ns["distortion_type"].value_counts())

   plot_ns_detection(
       sites,
       f_split=1.0,
       ns_threshold=2.0,
       ss_threshold=0.1,
       sort_by="lat",
       half_window=3,
       max_skew=45.0,
       show_ss=True,
   )

The ``distortion_type`` column has four possible values:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Meaning
   * - ``"clean"``
     - No strong near-surface index and no strong static-shift residual.
   * - ``"static"``
     - Static-shift residual is large, but the high-frequency residual
       spread is not.  A conventional static-shift correction may help.
   * - ``"near_surface"``
     - High-frequency residual spread is large, but static residual is not.
       Treat this as frequency-dependent distortion, not as a simple
       multiplier.
   * - ``"mixed"``
     - Both diagnostics are large.  Static correction may remove the
       constant part, but the remaining frequency-dependent effect needs
       separate geological and inversion judgment.

Recommended Processing Pattern
------------------------------

For a survey report, keep the static-shift workflow explicit:

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.emtools import (
       apply_ss_factors,
       detect_near_surface,
       ensure_sites,
       estimate_ss_ama,
       ss_comparison_psection,
   )

   edi_dir = Path("data/AMT/WILLY_DATA/L18PLT")
   sites = ensure_sites(edi_dir, recursive=True)

   factors = estimate_ss_ama(
       sites,
       sort_by="lat",
       half_window=3,
       pband=(1e-4, 10.0),
       max_skew=45.0,
   )

   factors.to_csv("static_shift_factors.csv", index=False)

   ns = detect_near_surface(
       sites,
       f_split=1.0,
       sort_by="lat",
       half_window=3,
       max_skew=45.0,
   )
   ns.to_csv("near_surface_diagnostics.csv", index=False)

   accepted = factors.merge(
       ns[["station", "distortion_type"]],
       on="station",
       how="left",
   )

   print(
       accepted.groupby("distortion_type", dropna=False)
       ["station"]
       .count()
       .rename("n_stations")
   )

   corrected = apply_ss_factors(
       sites,
       factors,
       key="fac_z",
       inplace=False,
   )

   fig = ss_comparison_psection(
       sites,
       method="ama",
       sort_by="lat",
       half_window=3,
       pband=(1e-4, 10.0),
       max_skew=45.0,
       show_delta=True,
   )

The important decisions are visible in that script: loader, station
ordering, period band, skew threshold, accepted estimator, saved factor
table, and a near-surface check before interpreting the correction.

Common Pitfalls
---------------

``sort_by`` does not describe the coordinate you like best; it describes
the along-line order used by the neighbourhood estimator.  A north-south
line usually wants ``sort_by="lat"``.  An east-west line usually wants
``sort_by="lon"``.

``max_skew`` is a filter, not a quality score.  A strict value can leave
too few samples for a complex survey.  Check ``n_used`` before trusting
the factors.

Static shift changes level, not shape.  If the before/after station plot
would need a frequency-dependent correction to align the curves, use the
near-surface diagnostics and interpret with caution.

Do not apply an empty or nearly empty factor table just because the code
returns without error.  Empty estimates are a valid no-op outcome for
single-station inputs or surveys where no station has usable neighbours.

Prefer ``fac_z`` for impedance correction.  ``fac_rho`` is the
resistivity-domain factor and is useful for interpretation, but
``apply_ss_factors`` scales ``Z``.

Worked Example
--------------

The gallery example below uses the L18PLT survey and walks through raw
curve spread, AMA estimation, exact factor application, estimator
comparison, QC wrappers, radar plots, and near-surface classification.

.. literalinclude:: ../../../examples/emtools/plot_ss.py
   :language: python
   :linenos:
