.. _emtools_tensor:

Phase Tensor And Impedance Tensor Tools
=======================================

The tensor page is one of the central pages in the pyCSAMT user guide
because it connects interpretation and preprocessing.  The same
impedance tensor ``Z`` supports apparent resistivity, phase, skew,
geoelectric strike, static-shift checks, dimensionality checks, and 2-D
rotation decisions.  The ``pycsamt.emtools`` tensor tools expose two
main workflows:

.. list-table::
   :header-rows: 1
   :widths: 30 32 38

   * - Workflow
     - Main tools
     - Purpose
   * - Phase-tensor interpretation
     - ``build_phase_tensor_table``,
       ``plot_phase_tensor_psection``,
       ``plot_phase_tensor_map``,
       ``plot_phase_tensor_summary``
     - Diagnose dimensionality, strike, skew, ellipticity, and spatial
       coherence without being dominated by static-shift amplitude.
   * - Impedance-tensor editing
     - ``rotate``, ``rotate_by_map``, ``rotate_z_to_strike``,
       ``antisymmetrize``, ``orient_from_sensors``, ``sigma_clip_z``,
       ``balance_offdiag``, ``invert``
     - Apply deliberate preprocessing operations to ``Z`` before
       inversion, plotting, or compatibility with a 2-D assumption.

The phase tensor follows the Caldwell et al. style decomposition.  For
each frequency, pyCSAMT splits the impedance tensor into real and
imaginary parts, then computes the phase tensor
``Phi = real(Z)^-1 imag(Z)`` and its invariants.

The examples in this guide use public two-level imports from
``pycsamt.emtools``.  One name needs special attention:
``pycsamt.emtools.rotate_to_strike`` belongs to the strike module, while
the tensor module's own rotation-to-strike helper is exported as
``rotate_z_to_strike``.

Load Data
---------

Start with the canonical loader.  It returns a ``Sites`` object and skips
stations without valid impedance data.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.emtools import ensure_sites

   edi_dir = Path("data/AMT/WILLY_DATA/L18PLT")
   sites = ensure_sites(
       edi_dir,
       recursive=True,
       on_dup="replace",
       strict=False,
       verbose=0,
   )

Keep ``sites`` as the unmodified reference object.  Tensor editing
functions can return corrected copies when ``inplace=False``, so it is
easy to compare original and edited data.

Build The Phase-Tensor Table
----------------------------

``build_phase_tensor_table`` is the foundation for the plotting tools.
It returns one row per station and frequency.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import build_phase_tensor_table

   pt = build_phase_tensor_table(
       sites,
       recursive=False,
   )

   print(pt.head())
   print(pt[["station", "freq", "period", "theta", "skew", "ellipt"]])

The table contains these core columns:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Column
     - Meaning
   * - ``station``
     - Station identifier.
   * - ``freq``
     - Frequency in hertz.
   * - ``period``
     - Period in seconds, computed as ``1 / freq``.
   * - ``s1`` and ``s2``
     - Phase-tensor principal values.  They control ellipse major and
       minor axes.
   * - ``theta``
     - Principal-axis angle, interpreted as phase-tensor strike.  It is
       axial, so directions separated by 180 degrees are equivalent.
   * - ``alpha``
     - Phase-tensor coordinate angle.
   * - ``beta`` and ``skew``
     - Skew angle.  ``skew`` is an alias for ``beta``.
   * - ``ellipt``
     - Ellipticity, computed from the principal values.  Values near zero
       are closer to circular; larger values indicate stronger 2-D
       anisotropy of the phase tensor.

The table is also the best place to audit a survey numerically before
plotting:

.. code-block:: python
   :linenos:

   summary = pt.groupby("station").agg(
       n=("freq", "count"),
       median_abs_skew=("skew", lambda values: values.abs().median()),
       median_ellipt=("ellipt", "median"),
       theta_iqr=("theta", lambda values: values.quantile(0.75) - values.quantile(0.25)),
   )

   print(summary.sort_values("median_abs_skew", ascending=False))

Large ``median_abs_skew`` means the station is not behaving like a clean
1-D or 2-D response in that period range.  Large ``theta_iqr`` means
phase-tensor strike changes strongly with frequency, so one rotation
angle may be a poor summary.

Filter By Period
----------------

Most tensor interpretation should be tied to a period band.  A shallow
band and a deeper band can show different strike, skew, or ellipticity.

.. code-block:: python
   :linenos:

   period_band = (0.001, 10.0)
   band_pt = pt[
       (pt["period"] >= period_band[0])
       & (pt["period"] <= period_band[1])
   ]

   print("rows in band:", len(band_pt))
   print("stations in band:", band_pt["station"].nunique())
   print("median |skew|:", band_pt["skew"].abs().median())

When you report a tensor result, always report the period band.  A map at
``period=1.0`` second and a summary over ``0.001`` to ``10.0`` seconds
answer different questions.

Read Dimensionality From Skew And Ellipticity
---------------------------------------------

A simple phase-tensor dimensionality rule uses skew and ellipticity:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Class
     - Rule
     - Interpretation
   * - 1-D
     - ``abs(skew) <= skew_threshold`` and
       ``abs(ellipt) <= ellipt_threshold``
     - Low skew and nearly circular phase tensor.
   * - 2-D
     - ``abs(skew) <= skew_threshold`` and
       ``abs(ellipt) > ellipt_threshold``
     - Low skew but elongated phase tensor.
   * - 3-D
     - ``abs(skew) > skew_threshold``
     - High skew; strike and 2-D rotation should be treated with caution.

.. code-block:: python
   :linenos:

   import numpy as np

   skew_threshold = 3.0
   ellipt_threshold = 0.2

   work = band_pt.copy()
   abs_skew = work["skew"].abs()
   abs_ellipt = work["ellipt"].abs()

   work["dimensionality"] = np.select(
       [
           (abs_skew <= skew_threshold) & (abs_ellipt <= ellipt_threshold),
           (abs_skew <= skew_threshold) & (abs_ellipt > ellipt_threshold),
       ],
       ["1D", "2D"],
       default="3D",
   )

   print(work["dimensionality"].value_counts(normalize=True))

The default ``3`` degree skew threshold is strict.  It is useful as a
textbook 1-D/2-D screen, but it can classify many real field samples as
3-D.  That is not a failure of the function; it is a warning about the
data and the 2-D assumption.

Simple Phase-Tensor Views
-------------------------

Use the simpler plots before the full ellipse plot.  They make it easy
to identify which invariant is causing the interpretation.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       plot_dimensionality_grid,
       plot_dimensionality_psection,
       plot_ellipticity_psection,
       plot_phase_tensor_skewmap,
       plot_theta_vs_period,
   )

   plot_theta_vs_period(sites, recursive=False)
   plot_phase_tensor_skewmap(sites, recursive=False, axis_y="logperiod")
   plot_ellipticity_psection(sites, recursive=False)
   plot_dimensionality_psection(
       sites,
       skew_th=3.0,
       ellipt_th=0.2,
       recursive=False,
   )
   plot_dimensionality_grid(
       sites,
       skew_th=3.0,
       ellipt_th=0.2,
       recursive=False,
   )

``plot_theta_vs_period`` is a scatter view of strike angle by period.
It is quick, but it puts an axial angle on a linear y-axis.  Treat jumps
near the wrap boundary carefully.

``plot_phase_tensor_skewmap`` and ``plot_ellipticity_psection`` show
station-by-period heatmaps.  They are good for finding period bands that
are consistently low-skew or strongly elongated.

``plot_dimensionality_psection`` and ``plot_dimensionality_grid`` apply
the simple skew/ellipticity classification to every station-period cell.

Phase-Tensor Ellipse Pseudosection
----------------------------------

``plot_phase_tensor_psection`` is the main phase-tensor figure.  Each
cell is an ellipse:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Visual element
     - Meaning
   * - Major axis
     - ``s1``, the larger phase-tensor principal value.
   * - Minor axis
     - ``s2``, the smaller phase-tensor principal value.
   * - Orientation
     - ``theta``, the phase-tensor strike direction.
   * - Fill color
     - Controlled by ``c_by``.  The default is usually skew.
   * - Thick border
     - Optional marker for cells where ``abs(skew)`` exceeds
       ``skew_threshold``.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_phase_tensor_psection

   plot_phase_tensor_psection(
       sites,
       stations=None,
       period_range=(0.001, 10.0),
       axis_y="logperiod",
       period_up=True,
       c_by="skew",
       skew_threshold=3.0,
       mark_3d=True,
       normalise_by="cell",
       recursive=False,
   )

Useful ``c_by`` values include ``"skew"``, ``"beta"``, ``"theta"``,
``"ellipt"``, ``"s1"``, ``"s2"``, ``"|skew|"``, ``"phi_mean"``,
``"phi_max"``, and ``"phi_min"``.

Use ``normalise_by="cell"`` for most survey plots because it scales the
ellipses to the local plotting grid.  Use ``normalise_by="unity"`` when
you want the 45-degree, 1-D reference to have an explicit visual meaning.
Use ``normalise_by="abs"`` only when absolute ellipse sizes in data units
are intentional.

Strike As A Director Field
--------------------------

``theta`` is axial.  A director field is often easier to interpret than
a scatter plot because the glyph has no arrow head and therefore
respects the 180-degree ambiguity.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_strike_director_field

   plot_strike_director_field(
       sites,
       color_by="skew",
       length_by="ellipt",
       skew_max=6.0,
       streamlines=True,
       period_subsample=40,
       recursive=False,
   )

Interpret long, aligned, low-skew directors as a more coherent 2-D
strike signal.  Short directors mean the phase tensor is close to
circular and strike is poorly defined.  High-skew directors mean the
response is more 3-D or distorted, even if the direction looks visually
organized.

Rose And Stability Plots
------------------------

Rose plots are useful for summarizing axial direction.  Stability plots
are useful for seeing whether that summary hides period dependence.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       plot_phase_tensor_rose,
       plot_theta_rose_grid,
       plot_theta_stability_stripe,
   )

   plot_phase_tensor_rose(
       sites,
       band=(0.001, 10.0),
       bins=36,
       recursive=False,
   )

   plot_theta_rose_grid(
       sites,
       n_bands=6,
       bins=24,
       recursive=False,
   )

   plot_theta_stability_stripe(
       sites,
       win=5,
       recursive=False,
   )

``plot_phase_tensor_rose`` folds all selected ``theta`` values into one
axial histogram.  ``plot_theta_rose_grid`` splits the period range into
equal log-width bands.  ``plot_theta_stability_stripe`` uses hue for
``theta`` and saturation for local stability.

If the rose grid changes direction from one period band to another,
avoid selecting a single strike for the entire data set.

Skew-Ellipticity Density
------------------------

``plot_skew_ellipt_density`` shows the joint distribution of
``abs(beta)`` and ``abs(ellipt)``.  It is a compact way to see whether
the data cluster in a 1-D, 2-D, or 3-D region.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_skew_ellipt_density

   plot_skew_ellipt_density(
       sites,
       band=(0.001, 10.0),
       gridsize=40,
       recursive=False,
   )

Use this plot with the dimensionality grid.  The grid tells you where
problem cells occur; the density plot tells you how the full population
is distributed.

Summary Figure
--------------

``plot_phase_tensor_summary`` combines the main ellipse pseudosection,
the dimensionality grid, and the skew-ellipticity density into one
figure.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_phase_tensor_summary

   fig = plot_phase_tensor_summary(
       sites,
       stations=None,
       period_range=(0.001, 10.0),
       c_by="skew",
       skew_threshold=3.0,
       ellipt_threshold=0.2,
       recursive=False,
   )

This is a good report figure when the audience needs the whole
phase-tensor story in one place: what the ellipses look like, how much
of the band is 1-D/2-D/3-D, and where the skew/ellipticity population
sits.

Geographic Phase-Tensor Map
---------------------------

``plot_phase_tensor_map`` draws one phase-tensor ellipse per station at
the period nearest a requested target period.  It can also overlay tipper
arrows when vertical magnetic transfer functions are present.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_phase_tensor_map

   plot_phase_tensor_map(
       sites,
       period=1.0,
       c_by="skew",
       show_tipper=True,
       tipper_convention="parkinson",
       station_labels=True,
       recursive=False,
   )

Use ``period`` as a target; each station uses its nearest available
period.  If the EDI headers do not provide usable coordinates, pass an
explicit ``coords`` dictionary:

.. code-block:: python
   :linenos:

   coords = {
       "S001": (7.312, -5.218),
       "S002": (7.318, -5.211),
       "S003": (7.324, -5.204),
   }

   plot_phase_tensor_map(
       sites,
       period=1.0,
       coords=coords,
       show_tipper=False,
       recursive=False,
   )

The coordinate tuple is ``(lat, lon)``.  A map with no coordinates is not
a tensor failure; it is a metadata problem.  Use pseudosections and
profiles until coordinates are supplied.

Per-Station Ellipse Strips
--------------------------

``plot_phase_tensor_strip`` draws one station as an ellipse sequence
through period.  ``plot_phase_tensor_strip_grid`` tiles several such
strips, optionally grouped by profile.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       plot_phase_tensor_strip,
       plot_phase_tensor_strip_grid,
   )

   plot_phase_tensor_strip(
       sites,
       station="18-016A",
       period_range=(0.001, 10.0),
       c_by="skew",
       recursive=False,
   )

   groups = {
       "L18PLT": ["18-001A", "18-002A", "18-003A", "18-004A"],
   }

   plot_phase_tensor_strip_grid(
       sites,
       groups,
       period_range=(0.001, 10.0),
       c_by="skew",
       recursive=False,
   )

Use strips when a single station deserves close inspection.  Use the
grid when comparing stations or profile groups without the compression
of a full pseudosection.

Standalone Legend
-----------------

``phase_tensor_legend`` draws a reference ellipse that can be composed
into custom figures.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import phase_tensor_legend

   phase_tensor_legend(size=1.0)

It is not a diagnostic by itself; it is a small plotting component for
figures where the phase-tensor ellipse convention needs to be explained.

Impedance-Tensor Editing
------------------------

Tensor editing functions change ``Z``.  They should be treated as
processing operations, not as harmless plots.  Keep ``inplace=False``
unless you deliberately want to mutate the object in memory.

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Function
     - What it does
     - Typical use
   * - ``rotate``
     - Rotates all station impedance tensors by one fixed angle.
     - Apply a known regional strike or sensor-frame correction.
   * - ``rotate_by_map``
     - Rotates each station by a value from a station-angle dictionary.
     - Apply station-specific rotations from an external interpretation.
   * - ``rotate_z_to_strike``
     - Estimates a tensor strike for each station and rotates by it.
     - Exploratory tensor-level strike rotation.
   * - ``antisymmetrize``
     - Enforces off-diagonal antisymmetry.
     - Prepare data for methods that assume ``Zxy = -Zyx``.
   * - ``orient_from_sensors``
     - Corrects electric and magnetic sensor orientations.
     - Fix known sensor azimuth errors.
   * - ``sigma_clip_z``
     - Sets outlying ``Z`` entries to ``NaN``.
     - Remove isolated spikes before plotting or inversion.
   * - ``balance_offdiag``
     - Balances off-diagonal magnitudes.
     - Exploratory symmetry conditioning.
   * - ``invert``
     - Replaces each 2-by-2 impedance tensor by its matrix inverse.
     - Advanced workflows that explicitly need admittance-like tensors.

Fixed-Angle Rotation
--------------------

Use ``rotate`` when every station should receive the same angle.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import rotate

   rotated_30 = rotate(
       sites,
       30.0,
       inplace=False,
       recursive=False,
   )

Use this only when the angle is justified by strike analysis, field
geometry, or a documented processing convention.

Station-Specific Rotation
-------------------------

Use ``rotate_by_map`` when each station needs its own angle.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import rotate_by_map

   angle_by_station = {
       "18-001A": 25.0,
       "18-002A": 27.5,
       "18-003A": 24.0,
   }

   rotated_stationwise = rotate_by_map(
       sites,
       angle_by_station,
       inplace=False,
       recursive=False,
   )

Stations missing from the dictionary receive a ``0`` degree rotation.
For production work, check that every intended station is present in the
map before applying it.

Tensor Rotation To Strike
-------------------------

The tensor module's rotation-to-strike function is exported as
``rotate_z_to_strike`` at the ``pycsamt.emtools`` level.  This avoids a
name collision with the strike module's ``rotate_to_strike``.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import rotate_z_to_strike

   rotated_to_tensor_strike = rotate_z_to_strike(
       sites,
       method="swift",
       inplace=False,
       recursive=False,
   )

Use this as a tensor editing step, not as the main strike-analysis
interface.  For detailed strike estimation, use the dedicated
``strike.rst`` workflow and its station-level tables.

Antisymmetrize And Balance
--------------------------

``antisymmetrize`` enforces off-diagonal antisymmetry.  ``balance_offdiag``
balances ``|Zxy|`` and ``|Zyx|`` while preserving phase.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import antisymmetrize, balance_offdiag

   antisymmetric = antisymmetrize(
       sites,
       how="rms",
       inplace=False,
       recursive=False,
   )

   balanced = balance_offdiag(
       sites,
       mode="avgabs",
       inplace=False,
       recursive=False,
   )

These operations can be useful for controlled experiments or for
preparing data for algorithms with strong 2-D assumptions.  They can
also hide real 3-D information if applied blindly.  Always keep an
uncorrected copy.

Sensor Orientation
------------------

Use ``orient_from_sensors`` when field notes show that the electric and
magnetic sensors were not aligned with the assumed axes.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import orient_from_sensors

   oriented = orient_from_sensors(
       sites,
       ex=5.0,
       ey=95.0,
       bx=5.0,
       by=95.0,
       degrees=True,
       inplace=False,
       recursive=False,
   )

Angles are degrees by default.  Pass ``degrees=False`` only when your
inputs are radians.  This operation is meaningful only when the sensor
orientation metadata or field notes are trustworthy.

Sigma Clip Outliers
-------------------

``sigma_clip_z`` flags outlying entries in the complex impedance tensor
and sets them to ``NaN``.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import sigma_clip_z

   clipped = sigma_clip_z(
       sites,
       sigma=3.0,
       inplace=False,
       recursive=False,
   )

After clipping, rerun coverage and QC checks.  A small number of clipped
entries may remove isolated spikes; a large number means the survey or
threshold needs review.

Invert Tensor
-------------

``invert`` applies the 2-by-2 matrix inverse frequency by frequency.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import invert

   admittance_like = invert(
       sites,
       inplace=False,
       recursive=False,
   )

This is an advanced operation.  Do not use it as a generic noise
correction.  It changes the physical quantity represented by the tensor.

Audit Edits With Phase-Tensor Tables
------------------------------------

A safe editing workflow compares before and after diagnostics.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       build_phase_tensor_table,
       plot_phase_tensor_summary,
       rotate,
   )

   before = build_phase_tensor_table(sites, recursive=False)

   rotated = rotate(
       sites,
       30.0,
       inplace=False,
       recursive=False,
   )

   after = build_phase_tensor_table(rotated, recursive=False)

   audit = before[["station", "period", "theta", "skew", "ellipt"]].merge(
       after[["station", "period", "theta", "skew", "ellipt"]],
       on=["station", "period"],
       suffixes=("_before", "_after"),
   )

   audit["theta_change"] = (
       (audit["theta_after"] - audit["theta_before"] + 90.0) % 180.0
   ) - 90.0

   print(audit[["station", "period", "theta_change", "skew_before", "skew_after"]])

   plot_phase_tensor_summary(
       rotated,
       period_range=(0.001, 10.0),
       recursive=False,
   )

The axial difference formula keeps ``theta`` comparisons honest across
the 180-degree wrap boundary.

Recommended Interpretation Workflow
-----------------------------------

For a survey report, keep the phase-tensor interpretation explicit:

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.emtools import (
       build_phase_tensor_table,
       ensure_sites,
       plot_phase_tensor_psection,
       plot_phase_tensor_summary,
       plot_skew_ellipt_density,
       plot_theta_rose_grid,
   )

   sites = ensure_sites(
       Path("data/AMT/WILLY_DATA/L18PLT"),
       recursive=True,
   )

   period_range = (0.001, 10.0)
   skew_threshold = 3.0
   ellipt_threshold = 0.2

   pt = build_phase_tensor_table(sites, recursive=False)
   pt_band = pt[
       (pt["period"] >= period_range[0])
       & (pt["period"] <= period_range[1])
   ]
   pt_band.to_csv("phase_tensor_table.csv", index=False)

   plot_phase_tensor_psection(
       sites,
       period_range=period_range,
       c_by="skew",
       skew_threshold=skew_threshold,
       mark_3d=True,
       recursive=False,
   )

   plot_theta_rose_grid(
       sites,
       n_bands=6,
       recursive=False,
   )

   plot_skew_ellipt_density(
       sites,
       band=period_range,
       recursive=False,
   )

   plot_phase_tensor_summary(
       sites,
       period_range=period_range,
       skew_threshold=skew_threshold,
       ellipt_threshold=ellipt_threshold,
       recursive=False,
   )

This workflow saves the table, plots the core ellipse section, checks
strike by band, checks skew/ellipticity distribution, and creates a
summary figure with the same thresholds.

Common Pitfalls
---------------

Do not treat ``theta`` as an ordinary linear angle.  It is axial.  Use an
axial difference formula when comparing directions.

Do not rotate data only because a strike value exists.  High skew,
unstable ``theta``, or strong period dependence can make the rotation
misleading.

Do not hide 3-D behavior by editing ``Z`` too early.  Plot the raw phase
tensor before applying antisymmetrization, balancing, clipping, or
rotation.

Do not interpret a missing map as missing tensor data.  A map also needs
coordinates.  If coordinates are missing, use pseudosections or pass an
explicit ``coords`` mapping.

Do not mix the two rotation-to-strike functions.  Use
``rotate_z_to_strike`` for the tensor editing helper exported from
``pycsamt.emtools``.  Use the dedicated strike page for the strike module
workflow.

Worked Example
--------------

The gallery example computes the phase-tensor table, builds the
main pseudosections and rose diagrams, demonstrates the map and strip
views, and audits impedance-tensor editing operations on real data.

Open the rendered gallery page here:
:ref:`sphx_glr_examples_emtools_plot_tensor.py`.
