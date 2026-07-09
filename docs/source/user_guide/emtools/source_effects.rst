.. _emtools_source_effects:

Source Effects And Near-Field Correction
========================================

``pycsamt.emtools.source_effects`` helps diagnose when the artificial
CSAMT transmitter is influencing the measured response. Natural-source
MT interpretation assumes a plane-wave source. CSAMT does not have that
luxury: the transmitter has a finite offset from each receiver, and the
offset can control whether a station-frequency row behaves like near
field, transition field, or far field.

The module contains two related but independent families of tools:

* Yan and Fu / Da et al. source-overprint diagnostics, based on the
  ground-wave to surface-wave amplitude ratio :math:`\beta_{Ey}`.
* Wang and Lin normalized-response and near-field correction tools,
  based on skin-depth field zones and an equatorial horizontal electric
  dipole correction factor.

Full function signatures and parameter defaults are maintained in the
:doc:`API reference <../../api/emtools>`. This guide uses the public
two-level imports from ``pycsamt.emtools``.

Why Offset Matters
------------------

Every source-effect calculation needs a source-receiver offset ``r``.
Standard EDI files usually do not store the CSAMT transmitter geometry,
so you must provide the offset explicitly unless your station objects
already carry an attribute such as ``source_offset``, ``offset``, or
``dist``.

The offset can be supplied as a scalar or as a station dictionary:

.. code-block:: python
   :linenos:

   source_offset = 2000.0

   source_offset_by_station = {
       "18-001A": 1800.0,
       "18-002A": 1950.0,
       "18-003A": 2100.0,
   }

Use a scalar only when the same representative offset is justified for
all stations. For field processing, prefer a station dictionary derived
from transmitter and receiver coordinates.

Workflow Map
------------

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - Goal
     - Use this
     - Output
   * - Evaluate pure overprint beta
     - ``overprint_beta``
     - :math:`\beta_{Ey}` in percent for arrays or scalars.
   * - Build per-frequency overprint table
     - ``detect_source_overprint``
     - Long-form table with ``beta_pct``, ``kr``, and flags.
   * - Summarize overprint by station
     - ``source_overprint_table``
     - Station table with max/mean beta, fraction flagged, and slopes.
   * - Plot beta pseudo-section
     - ``plot_overprint_section``
     - Station-period map of :math:`\beta_{Ey}`.
   * - Normalize response
     - ``normalize_response``
     - Apparent-resistivity ratio, phase residual, field zone, and ``kr``.
   * - Correct near-field response
     - ``correct_near_field``
     - ``Sites`` object with impedance divided by the near-field factor.
   * - Plot normalized response
     - ``plot_normalized_response``
     - Two-panel pseudo-section of normalized resistivity and phase.

Loading A Survey
----------------

Load the survey once with ``ensure_sites``. Keep the raw object unchanged
while you inspect source effects and test correction settings.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import ensure_sites

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True)
   source_offset = 2000.0

The examples below use ``2000.0`` meters only as a concrete processing
example. In a real CSAMT project, replace it with measured geometry.

Overprint Beta
--------------

``overprint_beta`` is the pure mathematical interface. It does not need
EDI files. It evaluates the Yan and Fu ground-wave to surface-wave ratio
and returns :math:`\beta_{Ey}` in percent.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.emtools import BETA_THRESH_PCT, overprint_beta

   freq = np.logspace(-1, 3, 60)
   rho = 300.0

   for offset in (500.0, 2000.0, 8000.0):
       beta_pct = overprint_beta(rho=rho, freq=freq, offset=offset)
       contaminated = freq[beta_pct > BETA_THRESH_PCT]
       if contaminated.size:
           print(
               f"offset={offset:g} m: beta>{BETA_THRESH_PCT:g}% "
               f"up to {contaminated.max():.3g} Hz"
           )

``BETA_THRESH_PCT`` is ``3.0``. Values above that threshold indicate
potential source overprint under the Yan and Fu criterion. The threshold
is useful, but the exact result depends strongly on ``rho``, frequency,
and offset.

Per-Frequency Overprint Detection
---------------------------------

``detect_source_overprint`` applies ``overprint_beta`` to every
station-frequency row using apparent resistivity computed from the
observed impedance tensor.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import detect_source_overprint, ensure_sites

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True)

   detail = detect_source_overprint(
       sites,
       source_offset=2000.0,
       beta_threshold=3.0,
   )

   print(detail.head())
   print(detail["beta_pct"].describe())
   print(detail["overprint_flag"].value_counts())

The returned table has one row per station and frequency:

.. code-block:: text

   station, freq_hz, period_s, offset_m, rho_a_ohmm,
   kr, beta_pct, overprint_flag

Rows with unknown offset keep the station and frequency information, but
``kr`` and ``beta_pct`` are ``NaN``. That is intentional: source-effect
diagnostics cannot be inferred honestly without geometry.

Station-Level Summary
---------------------

``source_overprint_table`` summarizes the long-form table by station. It
adds maximum and mean :math:`\beta`, the number and fraction of flagged
rows, and a low-/high-frequency slope comparison inspired by Da et al.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import ensure_sites, source_overprint_table

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True)

   summary = source_overprint_table(
       sites,
       source_offset=2000.0,
       beta_threshold=3.0,
       f_split=50.0,
   )

   cols = [
       "station",
       "beta_max_pct",
       "beta_mean_pct",
       "n_overprint",
       "overprint_frac",
       "lf_slope",
       "hf_slope",
       "slope_delta",
       "overprint_flag",
   ]
   print(summary[cols].sort_values("overprint_frac", ascending=False).head())

``f_split`` separates low- and high-frequency bands for slope analysis.
Choose it from the actual survey frequency range. If the split falls
outside the sampled range, one of the slope columns will be ``NaN``.

Overprint Pseudo-Section
------------------------

``plot_overprint_section`` maps :math:`\beta_{Ey}` across station and
period. It can contour key beta levels, including the 3 percent
threshold.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import ensure_sites, plot_overprint_section

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True)

   fig, ax = plt.subplots(figsize=(10, 5))
   plot_overprint_section(
       sites,
       source_offset=2000.0,
       beta_threshold=3.0,
       beta_levels=(1.0, 3.0, 10.0, 30.0),
       period_axis=True,
       log_y=True,
       ax=ax,
   )

Use this plot to see whether source contamination is localized to
specific stations, periods, or broad regions of the line. If most of the
plot sits above the threshold, the assumed offset may place much of the
survey outside a clean far-field regime.

Normalized Response
-------------------

``normalize_response`` implements the Wang and Lin view of source
effects. It computes:

.. math::

   \rho_n = \rho_\mathrm{obs} / \rho_\mathrm{ref}

.. math::

   \phi_\mathrm{diff} = \phi_\mathrm{obs} - \phi_\mathrm{ref}

It also classifies each row using the skin-depth relation
:math:`\delta = 503\sqrt{\rho_a/f}` and the source offset:

* ``near`` when ``r / delta < 0.5``;
* ``transition`` when ``0.5 <= r / delta < 4``;
* ``far`` when ``r / delta >= 4``.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import ensure_sites, normalize_response

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True)

   norm = normalize_response(
       sites,
       rho_ref=300.0,
       source_offset=2000.0,
       comp="det",
       phi_ref_deg=45.0,
   )

   print(norm.head())
   print(norm["zone"].value_counts(dropna=False))
   print(norm[["station", "freq_hz", "rho_n", "phi_diff_deg", "zone", "kr"]].head())

Use ``comp="det"`` for a determinant-style response, or ``"xy"`` /
``"yx"`` when a specific off-diagonal component is the interpretation
target.

Normalized-Response Plot
------------------------

``plot_normalized_response`` draws the normalized resistivity and
subtracted phase as two side-by-side pseudo-sections.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import ensure_sites, plot_normalized_response

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True)

   fig, axes = plt.subplots(1, 2, figsize=(13, 5))
   plot_normalized_response(
       sites,
       rho_ref=300.0,
       source_offset=2000.0,
       comp="det",
       phi_ref_deg=45.0,
       axes=axes,
   )
   fig.tight_layout()

The left panel answers whether apparent resistivity is high or low
relative to the reference half-space. The right panel answers whether
phase is above or below the reference phase. Read both panels alongside
the field-zone column from ``normalize_response``.

Near-Field Correction
---------------------

``correct_near_field`` divides each impedance tensor row by a complex
near-field factor:

.. math::

   Z_\mathrm{corrected} = Z_\mathrm{observed} / F(p)

.. math::

   F(p) = 1 - 3/p^2 + 3/p^3

The factor tends toward ``1`` in the far field. In the near field it can
be very large, so the correction can strongly change apparent
resistivity.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import correct_near_field, ensure_sites

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True)

   corrected = correct_near_field(
       sites,
       source_offset=2000.0,
       inplace=False,
   )

Use ``inplace=False`` while testing. If a correction changes a station by
orders of magnitude, treat that as a diagnostic result, not just a
processed output. It means the raw response was far from the plane-wave
assumption under the supplied offset.

Comparing Before And After
--------------------------

You can compare source-effect diagnostics before and after correction
without reaching into private helpers. For example, compare normalized
response tables:

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       correct_near_field,
       ensure_sites,
       normalize_response,
   )

   raw = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True)
   corrected = correct_near_field(raw, source_offset=2000.0, inplace=False)

   before = normalize_response(raw, rho_ref=300.0, source_offset=2000.0)
   after = normalize_response(corrected, rho_ref=300.0, source_offset=2000.0)

   joined = before.merge(
       after,
       on=["station", "freq_hz"],
       suffixes=("_raw", "_corrected"),
   )
   joined["rho_n_ratio"] = joined["rho_n_corrected"] / joined["rho_n_raw"]

   print(
       joined[
           ["station", "freq_hz", "zone_raw", "rho_n_raw", "rho_n_corrected", "rho_n_ratio"]
       ].head()
   )

This style keeps the comparison in public tables and is easier to
document than extracting impedance arrays directly.

Combining The Two Diagnostics
-----------------------------

The Yan/Fu beta flag and Wang/Lin field-zone label come from different
physical arguments. Agreement between them is a strong warning that
source geometry is controlling part of the response.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       detect_source_overprint,
       ensure_sites,
       normalize_response,
   )

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True)

   beta = detect_source_overprint(sites, source_offset=2000.0)
   zones = normalize_response(sites, rho_ref=300.0, source_offset=2000.0)

   merged = beta.merge(
       zones[["station", "freq_hz", "zone", "kr"]],
       on=["station", "freq_hz"],
       how="left",
   )

   print(merged.groupby("zone")["overprint_flag"].mean())
   print(merged.groupby("zone")["beta_pct"].describe())

If ``near`` and ``transition`` rows are usually overprint-flagged while
``far`` rows are mostly unflagged, the two methods are telling a
consistent story. If they disagree, inspect the assumed offset,
reference resistivity, and frequency range.

Choosing Offsets
----------------

For real processing, offsets should come from survey geometry. A useful
pattern is to build a station dictionary and pass it to every function:

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       detect_source_overprint,
       normalize_response,
       source_overprint_table,
   )

   offset_by_station = {
       "18-001A": 1800.0,
       "18-002A": 1900.0,
       "18-003A": 2050.0,
       # Continue for the full line.
   }

   detail = detect_source_overprint(sites, source_offset=offset_by_station)
   summary = source_overprint_table(sites, source_offset=offset_by_station)
   norm = normalize_response(sites, source_offset=offset_by_station)

Keep the same offset dictionary across all source-effect diagnostics so
the beta table, normalized-response table, field zones, and correction
are comparable.

Suggested Review Sequence
-------------------------

Use this sequence before applying a correction:

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       detect_source_overprint,
       normalize_response,
       source_overprint_table,
   )

   detail = detect_source_overprint(sites, source_offset=2000.0)
   summary = source_overprint_table(sites, source_offset=2000.0, f_split=50.0)
   norm = normalize_response(sites, rho_ref=300.0, source_offset=2000.0)

   print(detail["overprint_flag"].mean())
   print(summary.sort_values("overprint_frac", ascending=False).head())
   print(norm["zone"].value_counts(dropna=False))

Then plot:

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import plot_normalized_response, plot_overprint_section

   fig, ax = plt.subplots(figsize=(10, 5))
   plot_overprint_section(sites, source_offset=2000.0, ax=ax)

   fig, axes = plt.subplots(1, 2, figsize=(13, 5))
   plot_normalized_response(
       sites,
       rho_ref=300.0,
       source_offset=2000.0,
       axes=axes,
   )

Correct only after the diagnostics show that correction is scientifically
justified and after the offset geometry has been checked.

Pitfalls
--------

Do not invent the source offset from the impedance. The offset is survey
geometry and should come from field records or transmitter/receiver
coordinates.

Do not interpret a representative scalar offset as a final result for a
line with varying transmitter distance. A scalar is fine for examples or
sensitivity tests; station-specific offsets are better for processing.

Do not treat near-field correction as harmless smoothing. It changes the
impedance tensor and can alter apparent resistivity by large factors in
near-field rows.

Do not use the Da et al. slope columns without checking ``f_split``. A
split outside the sampled frequency range produces undefined low- or
high-frequency slopes.

Worked Example
--------------

The example below uses the L18PLT survey with an explicitly stated
representative offset. It demonstrates the pure beta formula, per-row
overprint detection, station summaries, overprint pseudo-sections,
normalized response, near-field correction, and comparison between the
two independent source-effect diagnostics.

.. literalinclude:: ../../../examples/emtools/plot_source_effects.py
   :language: python
   :linenos:
