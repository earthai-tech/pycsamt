.. _emtools_strike:

Geoelectric Strike
==================

Geoelectric strike is the preferred 2-D structural direction inferred
from electromagnetic data.  It is one of the checks you make before
rotating impedances, preparing 2-D inversion inputs, or interpreting
transverse electric and transverse magnetic modes.  In pyCSAMT, the
strike tools estimate the direction from impedance rotation, phase
tensor azimuth, or a consensus of both, then visualize the result as
tables, profiles, ribbons, map-sticks, and rose diagrams.

Strike angles are axial.  A strike of ``0`` degrees and a strike of
``180`` degrees are the same geological direction.  For that reason,
pyCSAMT reports station-level strike in the compact ``[-90, 90]`` range
for estimator tables, while rose diagrams fold angles into ``0`` to
``180`` and mirror the histogram around the full polar circle.

Use this page when you need to answer concrete processing questions:

.. list-table::
   :header-rows: 1
   :widths: 30 32 38

   * - Question
     - Main tools
     - Output
   * - What strike does each station prefer?
     - ``estimate_strike_sweep``,
       ``estimate_strike_phase_tensor``,
       ``estimate_strike_consensus``
     - One row per station with angle, stability, period band, and sample
       count.
   * - Does strike change with period?
     - ``strike_curve_sweep``, ``plot_strike_ribbon``
     - One angle per station and frequency, plus a ribbon image.
   * - What strike should I rotate to?
     - ``rotate_to_strike``
     - A corrected ``Sites`` object, with impedance tensors rotated by the
       selected station-level angles.
   * - How do lines or bands compare?
     - ``plot_strike_rose``, ``plot_strike_rose_by_line``
     - Axial rose diagrams with weighted mean directions.
   * - Is the strike spatially coherent?
     - ``plot_strike_profile``, ``plot_strike_mapsticks``
     - Along-line and geographic views.
   * - How does Z strike compare with phase tensor and tipper direction?
     - ``plot_strike_analysis``
     - Three-panel strike, PT azimuth, and tipper rose diagram.

The examples use the public two-level import style:
``from pycsamt.emtools import ...``.

Load Data
---------

Load the survey with ``ensure_sites`` first.  This gives every strike
function the same clean input and avoids repeating EDI parsing options in
every call.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.emtools import ensure_sites

   edi_dir = Path("data/AMT/WILLY_DATA/L18PLT")
   sites = ensure_sites(edi_dir, recursive=True, verbose=0)

For map and profile plots, coordinates matter.  ``plot_strike_profile``
can order stations by ``"lon"``, ``"lat"``, ``"name"``, or ``"auto"``.
Choose the ordering that matches the survey line, not merely the default.
For a north-south line, ``sort_by="lat"`` is often the clearer choice.
For an east-west line, ``sort_by="lon"`` is usually better.

Station-Level Estimators
------------------------

pyCSAMT provides three station-level strike estimators.  They all return
a ``pandas.DataFrame`` with the same practical columns:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Meaning
   * - ``station``
     - Station identifier.
   * - ``ang``
     - Estimated strike angle in degrees, wrapped into ``[-90, 90]``.
   * - ``iqr``
     - Interquartile range of the frequency-level estimates used to
       summarize the station.  Smaller values mean a more stable strike.
   * - ``lo`` and ``hi``
     - Period-band limits in seconds used by the estimate.
   * - ``n``
     - Number of frequency samples used.

The impedance sweep rotates each tensor through a grid of trial angles
and chooses the angle that optimizes a metric.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.emtools import estimate_strike_sweep

   sweep = estimate_strike_sweep(
       sites,
       angles=np.arange(-90.0, 91.0, 1.0),
       metric="diag_ratio",
       band=(0.001, 10.0),
   )

   print(sweep[["station", "ang", "iqr", "n"]])

The default ``metric="diag_ratio"`` searches for the rotation that
minimizes diagonal energy relative to off-diagonal energy.  This is a
useful impedance-based strike diagnostic, but it can be sensitive to
noise, 3-D effects, and weak diagonal/off-diagonal contrast.

The phase-tensor estimator summarizes the phase tensor ``theta`` angle.
It is often more stable than a raw impedance sweep because the phase
tensor is less affected by static-shift amplitude distortion.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import estimate_strike_phase_tensor

   pt = estimate_strike_phase_tensor(
       sites,
       band=(0.001, 10.0),
       robust=True,
   )

   print(pt[["station", "ang", "iqr", "n"]])

The consensus estimator blends the sweep and phase-tensor estimates.
Use it when neither method should dominate the processing decision.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import estimate_strike_consensus

   consensus = estimate_strike_consensus(
       sites,
       band=(0.001, 10.0),
       w_sweep=0.4,
       w_pt=0.6,
       metric="diag_ratio",
   )

   print(consensus[["station", "ang", "iqr", "n"]])

For all three tables, treat ``iqr`` as a stability warning.  A station
with an angle near ``30`` degrees and an ``iqr`` near ``80`` degrees does
not have a reliable single strike; it has a broad or frequency-dependent
strike population.

Compare Axial Angles Correctly
------------------------------

Do not compare strike estimates with ordinary subtraction unless you
first account for the ``180`` degree ambiguity.  The axial difference
between ``89`` and ``-89`` degrees is ``2`` degrees, not ``178`` degrees.

.. code-block:: python
   :linenos:

   merged = sweep.merge(
       pt,
       on="station",
       suffixes=("_sweep", "_pt"),
   )

   axial_diff = (
       (merged["ang_sweep"] - merged["ang_pt"] + 90.0) % 180.0
   ) - 90.0

   merged["abs_axial_diff"] = axial_diff.abs()

   print(
       merged[
           ["station", "ang_sweep", "ang_pt", "abs_axial_diff"]
       ].sort_values("abs_axial_diff", ascending=False)
   )

Use this pattern when comparing sweep, phase tensor, consensus, tipper
azimuth, or externally interpreted structural trends.  A naive Pearson
correlation of raw angles can be misleading because it treats the wrap
boundary as a real discontinuity.

Choose A Period Band
--------------------

The ``band`` argument is a period band in seconds.  It is available on
the station-level estimators and on the high-level plots that summarize
station-level strike.  Use it to separate shallow, high-frequency
behavior from deeper, long-period behavior.

.. code-block:: python
   :linenos:

   short_period = estimate_strike_consensus(
       sites,
       band=(0.001, 0.1),
   )

   long_period = estimate_strike_consensus(
       sites,
       band=(0.1, 10.0),
   )

   band_compare = short_period[["station", "ang", "iqr"]].merge(
       long_period[["station", "ang", "iqr"]],
       on="station",
       suffixes=("_short", "_long"),
   )

   band_compare["band_axial_diff"] = (
       (band_compare["ang_short"] - band_compare["ang_long"] + 90.0)
       % 180.0
   ) - 90.0

   print(band_compare)

If short- and long-period strikes disagree strongly, do not force a
single rotation across the entire band.  Review dimensionality, static
shift, near-surface diagnostics, and the inversion band before choosing a
processing strike.

Rotate Data Onto Strike
-----------------------

``rotate_to_strike`` estimates one strike angle per station and rotates
that station's impedance tensor.  Keep the original and rotated data
separate until you have checked the result.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import rotate_to_strike

   rotated = rotate_to_strike(
       sites,
       method="consensus",
       band=(0.001, 10.0),
       metric="diag_ratio",
       inplace=False,
   )

   before = estimate_strike_consensus(
       sites,
       band=(0.001, 10.0),
   )
   after = estimate_strike_consensus(
       rotated,
       band=(0.001, 10.0),
   )

   print("before mean abs strike:", before["ang"].abs().mean())
   print("after mean abs strike:", after["ang"].abs().mean())

Valid method names are ``"consensus"``, ``"sweep"``, and ``"pt"``.
Use ``inplace=False`` while building a workflow.  It returns a rotated
copy and keeps the unrotated survey available for before/after checks.

Rotation does not make a survey 2-D by itself.  If the selected band has
high skew, unstable strike, or strong station-to-station disagreement,
the rotated tensors may still be poor 2-D inversion input.

Per-Frequency Strike Curve
--------------------------

``strike_curve_sweep`` keeps the frequency dimension instead of reducing
each station to one number.  It is useful for finding unstable bands or
stations whose strike flips with period.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import strike_curve_sweep

   curve = strike_curve_sweep(
       sites,
       angles=np.arange(-90.0, 91.0, 1.0),
       metric="diag_ratio",
       smooth=5,
   )

   print(curve.head())
   print(curve.groupby("station")["ang"].agg(["median", "std", "count"]))

The table columns are ``station``, ``freq``, and ``ang``.  ``smooth``
applies a moving average to the frequency-level sweep angles before they
are wrapped back into the axial range.  Increase it only when you want a
smoother visual trend; do not use smoothing to hide genuine strike
changes.

Ribbon Plot
-----------

``plot_strike_ribbon`` converts the per-frequency strike curve to a
station-by-period image.  Hue encodes strike angle.  Saturation encodes
local stability: desaturated colors indicate high local variance.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_strike_ribbon

   plot_strike_ribbon(
       sites,
       method="sweep",
       win=5,
       show_colorbar=True,
   )

Use the ribbon before selecting a single strike for a broad period band.
If the ribbon changes color systematically from short period to long
period, the survey may need band-specific interpretation.

Rose Diagrams
-------------

``plot_strike_rose`` draws axial strike histograms.  It mirrors the
``0`` to ``180`` degree histogram around the full circle, so both halves
of the polar plot represent the same set of axes.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_strike_rose

   fig = plot_strike_rose(
       sites,
       method="consensus",
       band=(0.001, 10.0),
       bins=36,
       weight="inv_iqr",
       suptitle="Consensus geoelectric strike",
   )

``weight="inv_iqr"`` down-weights stations whose strike varies strongly
with frequency.  Use ``weight="uniform"`` when every station should
contribute equally.

When ``groups`` is omitted, pyCSAMT attempts to group stations by a
profile-like station-name prefix.  This works for names such as
``E1S01`` because the inferred group is ``E1``.  For station names that
do not encode the line this way, pass an explicit mapping.

.. code-block:: python
   :linenos:

   from pathlib import Path

   line18 = sorted(Path("data/AMT/WILLY_DATA/L18PLT").glob("*.edi"))
   line22 = sorted(Path("data/AMT/WILLY_DATA/L22PLT").glob("*.edi"))

   groups = {
       "L18PLT": [path.stem for path in line18],
       "L22PLT": [path.stem for path in line22],
   }

   fig = plot_strike_rose(
       line18 + line22,
       groups=groups,
       method="consensus",
       bins=36,
       n_cols=2,
       suptitle="Strike by profile line",
   )

``plot_strike_rose_by_line`` is a simpler line-comparison helper.  It
requires at least two stations per group; if automatic grouping produces
only singleton groups, pass the explicit ``groups`` dictionary yourself.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_strike_rose_by_line

   fig = plot_strike_rose_by_line(
       line18 + line22,
       groups=groups,
       method="consensus",
       band=(0.001, 10.0),
       weight="inv_iqr",
   )

Frequency-Band Roses
--------------------

Use ``bar_style="bands"`` when you want one rose diagram to show several
period bands.  Each band contributes its own stacked histogram.

.. code-block:: python
   :linenos:

   fig = plot_strike_rose(
       sites,
       method="consensus",
       bar_style="bands",
       freq_bands=[
           (0.001, 0.01),
           (0.01, 0.1),
           (0.1, 10.0),
       ],
       band_labels=[
           "very short period",
           "short period",
           "long period",
       ],
       suptitle="Strike by period band",
   )

If the bands stack around different mean directions, report the
band-specific behavior instead of collapsing it to one survey-wide
number.

Profile And Map-Stick Views
---------------------------

``plot_strike_profile`` shows strike angle along station order with an
IQR ribbon.  It is the best quick check for station-to-station
coherence.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_strike_profile

   plot_strike_profile(
       sites,
       method="consensus",
       band=(0.001, 10.0),
       sort_by="lat",
   )

The profile plot uses ``ang`` as the line and ``iqr`` as the uncertainty
ribbon.  A coherent 2-D line should not show random jumps from station to
station unless there is a real geological or data-quality reason.

``plot_strike_mapsticks`` draws a short line segment at each station
coordinate, oriented along the estimated strike.  It is useful for
checking whether nearby stations point in a consistent direction.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_strike_mapsticks

   plot_strike_mapsticks(
       sites,
       method="consensus",
       band=(0.001, 10.0),
       len_deg=0.02,
   )

The ``len_deg`` value is a display length in coordinate degrees, not a
geological length.  Adjust it for readability when the survey extent is
very small or very large.

Combined Strike Analysis
------------------------

``plot_strike_analysis`` creates a three-panel rose figure: impedance
strike, phase-tensor azimuth, and tipper strike.  The tipper panel is
empty when the data do not contain vertical magnetic transfer functions.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_strike_analysis

   fig = plot_strike_analysis(
       sites,
       method="consensus",
       band=(0.001, 10.0),
       bins=36,
       suptitle="Strike, phase tensor, and tipper azimuth",
   )

Use this figure as a consistency check.  If Z strike and phase-tensor
azimuth cluster around the same axial direction, confidence increases.
If tipper strike is available and points somewhere else, investigate
regional 3-D structure, coast effects, cultural noise, or sign
convention before choosing a rotation angle.

Recommended Workflow
--------------------

A robust strike workflow keeps estimation, comparison, visualization, and
rotation separate:

.. code-block:: python
   :linenos:

   from pathlib import Path

   import numpy as np

   from pycsamt.emtools import (
       ensure_sites,
       estimate_strike_consensus,
       estimate_strike_phase_tensor,
       estimate_strike_sweep,
       plot_strike_profile,
       plot_strike_rose,
       rotate_to_strike,
   )

   sites = ensure_sites(
       Path("data/AMT/WILLY_DATA/L18PLT"),
       recursive=True,
   )

   band = (0.001, 10.0)

   sweep = estimate_strike_sweep(
       sites,
       band=band,
       angles=np.arange(-90.0, 91.0, 1.0),
   )
   pt = estimate_strike_phase_tensor(sites, band=band)
   consensus = estimate_strike_consensus(
       sites,
       band=band,
       w_sweep=0.4,
       w_pt=0.6,
   )

   comparison = consensus[["station", "ang", "iqr"]].merge(
       pt[["station", "ang", "iqr"]],
       on="station",
       suffixes=("_consensus", "_pt"),
   )
   comparison["consensus_pt_diff"] = (
       (comparison["ang_consensus"] - comparison["ang_pt"] + 90.0)
       % 180.0
   ) - 90.0

   comparison.to_csv("strike_comparison.csv", index=False)
   consensus.to_csv("strike_consensus.csv", index=False)

   plot_strike_profile(
       sites,
       method="consensus",
       band=band,
       sort_by="lat",
   )
   plot_strike_rose(
       sites,
       method="consensus",
       band=band,
       weight="inv_iqr",
   )

   stable = consensus["iqr"].median() < 45.0
   if stable:
       rotated = rotate_to_strike(
           sites,
           method="consensus",
           band=band,
           inplace=False,
       )

The ``stable`` threshold in this example is only a processing rule of
thumb.  Choose the final threshold based on survey purpose, dimensionality
diagnostics, period band, and inversion assumptions.

Common Pitfalls
---------------

Strike has a ``180`` degree ambiguity.  Always compare angles with an
axial difference formula.

``iqr`` is not decoration.  A high-IQR strike is unstable across
frequency and should not be used blindly as a rotation angle.

The period band is part of the result.  A strike estimated over
``(0.001, 0.1)`` seconds may not match a strike estimated over
``(0.1, 10.0)`` seconds.

Automatic rose grouping depends on station names.  If your station names
do not encode line membership, pass ``groups`` explicitly.

``rotate_to_strike`` rotates by station-level estimates.  It does not
guarantee a single regional strike, and it does not remove 3-D structure.

Map-stick plots require usable station coordinates.  If no coordinates
are available, use profile and rose diagrams instead.

Worked Example
--------------

The gallery example below uses L18PLT, adds L22PLT for a multi-line rose
comparison, demonstrates estimator agreement, rotates data onto strike,
and builds ribbon, rose, map-stick, profile, and combined strike-analysis
figures.

.. literalinclude:: ../../../examples/emtools/plot_strike.py
   :language: python
   :linenos:
