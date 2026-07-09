.. _emtools_tf:

Transfer Functions And Tipper Diagnostics
=========================================

The transfer-function tools in ``pycsamt.emtools`` focus on the
vertical-field response, usually called the tipper.  The tipper relates
the vertical magnetic field to the horizontal magnetic field:

.. code-block:: text

   Hz = Tx Hx + Ty Hy

where ``Tx`` and ``Ty`` are complex, frequency-dependent transfer
functions.  When the subsurface is laterally uniform, the vertical field
is weak.  When current is channelled by lateral conductors or sharp
resistivity contrasts, the tipper grows and induction arrows become one
of the fastest qualitative diagnostics for conductor position, strike,
and period-dependent structure.

This page covers two related workflows:

.. list-table::
   :header-rows: 1
   :widths: 30 32 38

   * - Workflow
     - Main tools
     - Purpose
   * - EDI or ``Sites`` tipper diagnostics
     - ``plot_tipper_hodograms``, ``plot_tipper_polar``,
       ``plot_induction_map``, ``plot_induction_section``,
       ``plot_induction_rose``
     - Work from assembled transfer functions, usually EDI files.
   * - Spectra-direct tipper diagnostics
     - ``plot_induction_map_from_spectra``,
       ``plot_tipper_polar_from_spectra``,
       ``plot_induction_rose_from_spectra``
     - Work from spectral estimates before a final EDI has been written.

All public functions used below are exported from ``pycsamt.emtools``,
so the examples use the two-level import style.

Use A Dataset With Tipper
-------------------------

Many AMT/CSAMT data sets contain only horizontal electric and magnetic
channels.  Those files can have valid impedance but no vertical magnetic
transfer function.  The tipper functions will then return graceful
``"no tipper"`` messages.

For induction-vector work, first verify that the selected survey really
has tipper data.  The bundled KAP03 long-period MT profile is useful for
examples because it includes vertical-field measurements.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.emtools import ensure_sites

   edi_dir = Path("data/MT/kap03lmt_edis")
   sites = ensure_sites(
       edi_dir,
       recursive=True,
       on_dup="replace",
       strict=False,
       verbose=0,
   )

If a plot says ``"no tipper"``, check the data before changing plotting
options.  Missing tipper is a data-content issue, not necessarily a
failed plot.

What The Tipper Stores
----------------------

For each station and frequency, pyCSAMT expects a two-component complex
tipper:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Quantity
     - Meaning
   * - ``Tx``
     - Complex coefficient relating ``Hx`` to ``Hz``.
   * - ``Ty``
     - Complex coefficient relating ``Hy`` to ``Hz``.
   * - ``real(T)``
     - In-phase part.  Commonly used for Parkinson induction arrows.
   * - ``imag(T)``
     - Quadrature part.  Useful for checking frequency-dependent or
       inductive behavior that is out of phase with the horizontal field.
   * - ``abs(T)``
     - Magnitude, often summarized as
       ``sqrt(abs(Tx)**2 + abs(Ty)**2)``.

You usually do not need to extract these arrays manually.  The plotting
functions read them from the site objects.  Still, understanding the
components helps interpret the figures.

Choose Periods And Bands
------------------------

Tipper diagnostics are strongly period-dependent.  A station may be weak
at short period, strong at mid-period, and weak again at long period.
Choose periods and period bands deliberately.

.. code-block:: python
   :linenos:

   import numpy as np

   # Example period choices for a broad-band MT profile.
   periods = [25.0, 650.0, 2000.0, 17000.0]
   short_band = (25.0, 200.0)
   long_band = (2000.0, 17000.0)

   print("periods:", periods)
   print("short band:", short_band)
   print("long band:", long_band)

Use the same period choices across maps, roses, and sections when you
want the figures to support one interpretation.

Single-Station Hodograms
------------------------

Start with ``plot_tipper_hodograms`` when inspecting one station.  It
plots ``Tx`` and ``Ty`` in the complex plane, with colors split by period
band.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_tipper_hodograms

   fig = plot_tipper_hodograms(
       sites,
       station="kap151",
       bands=[
           (25.0, 200.0),
           (200.0, 2000.0),
           (2000.0, 17000.0),
       ],
       unit_circle=True,
       normalize=False,
   )

Read a hodogram before reading arrows.  It shows whether a station has a
smooth, coherent complex response or a scattered cloud.  A large loop
outside the unit circle can be real for strong 3-D/lateral induction; it
is not automatically an error.

Set ``normalize=True`` only when comparing shape rather than amplitude.
For conductor-strength interpretation, keep the raw amplitude.

Single-Station Polar View
-------------------------

``plot_tipper_polar`` converts one station's tipper into azimuth and
magnitude versus period.  The polar angle is the tipper direction, radius
is magnitude, and color is log-period.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_tipper_polar

   plot_tipper_polar(
       sites,
       station="kap151",
       component="real",
   )

Valid components are ``"real"``, ``"imag"``, and ``"abs"``.  Use
``"real"`` for a Parkinson-style conductor-direction reading, use
``"imag"`` to inspect the quadrature response, and use ``"abs"`` when
you mainly care about magnitude.

Induction Map At One Period
---------------------------

``plot_induction_map`` draws real and imaginary induction arrows at a
single target period.  The function picks the nearest available period
for each station.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_induction_map

   plot_induction_map(
       sites,
       period=2000.0,
       convention="park",
       show_real=True,
       show_imag=True,
       scale=4.0,
       station_labels=True,
       reference_arrow=0.1,
   )

The station coordinates come from easting/northing, x/y, or lon/lat when
available.  If none are present, pyCSAMT falls back to an index along a
line.  That fallback is still useful for a profile, but do not interpret
the x-axis as real distance unless the source data contain real
coordinates.

``scale`` controls arrow length in plot coordinates.  If arrows are too
small or overlap badly, adjust ``scale`` rather than changing the tipper
data.

Compare Several Periods On One Axis
-----------------------------------

``plot_induction_arrows`` overlays arrows from several requested periods
on one profile axis.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_induction_arrows

   plot_induction_arrows(
       sites,
       periods=[25.0, 650.0, 2000.0, 17000.0],
       convention="park",
       scale=1.0,
       normalize=True,
       strike_ticks=False,
   )

Use this for a fast period comparison, not as the final publication
figure.  Many periods on one axis can become visually crowded.  If the
period behavior matters, follow with a period section or a multi-period
map.

Sign Conventions
----------------

Induction-vector interpretation depends on convention.  The two common
views are Parkinson and Wiese.  They are rotated relative to each other,
so a figure can be misread if the convention is not stated.

``plot_induction_convention`` puts Parkinson/Wiese and real/imaginary
components in one 2-by-2 figure.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_induction_convention

   plot_induction_convention(
       sites,
       period=650.0,
       station_labels=False,
   )

Use this plot when communicating with collaborators or comparing to a
paper.  It makes sign and component choices visible instead of leaving
them implicit.

Period Pseudosection
--------------------

``plot_induction_section`` shows tipper magnitude or component strength
over stations and period.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_induction_section

   plot_induction_section(
       sites,
       component="abs",
       n_periods=30,
       cmap="RdBu_r",
       section="pseudosection",
   )

Use ``component="abs"`` for anomaly strength, ``"real"`` for in-phase
strength, and ``"imag"`` for quadrature strength.  A section is the best
single view for answering: where along the line is the tipper strong,
and at what periods?

Induction Rose
--------------

``plot_induction_rose`` summarizes arrow azimuths over all stations and
selected periods.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_induction_rose

   plot_induction_rose(
       sites,
       component="real",
       pband=(25.0, 200.0),
       nbins=36,
       title="Short-period induction azimuths",
   )

   plot_induction_rose(
       sites,
       component="real",
       pband=(2000.0, 17000.0),
       nbins=36,
       title="Long-period induction azimuths",
   )

Compare short and long period roses before claiming a regional conductor.
A short-period rose may be scattered because shallow heterogeneity points
in many directions.  A long-period rose that tightens into one sector can
support a deeper, more coherent conductive structure.

Multi-Period Map
----------------

``plot_induction_multiperiod_map`` stacks one map panel per period and
is the most report-ready induction-vector figure.  It can use real EDI
tipper, or an explicit ``tipper_data`` override.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_induction_multiperiod_map

   fig, axes = plot_induction_multiperiod_map(
       sites,
       periods=[25.0, 650.0, 2000.0, 17000.0],
       convention="park",
       arrow_scale=6.0,
       reference_arrow=0.1,
       show_background_cbar=False,
       station_labels=False,
       title="Induction vectors across period",
   )

When ``background`` is not supplied, the function draws a synthetic
terrain-like background.  That background is a visual placeholder, not a
real DEM.  For a report, pass your own ``background`` and
``background_extent``.

The fallback EDI read path in this function can only use a single tipper
component in some situations.  When you need full two-component vectors,
pass ``tipper_data`` explicitly as a dictionary keyed by period:

.. code-block:: python
   :linenos:

   import numpy as np

   # Each value is an array with shape (n_stations, 2):
   # column 0 is Tx, column 1 is Ty. Replace these synthetic
   # rows with processed Tx/Ty values from your own workflow.
   n_stations = 26
   tipper_data = {}
   for period, tx, ty in [
       (25.0, 0.08 + 0.02j, 0.02 + 0.01j),
       (650.0, 0.22 + 0.05j, 0.10 + 0.03j),
       (2000.0, 0.16 + 0.03j, 0.18 + 0.04j),
       (17000.0, 0.05 + 0.01j, 0.20 + 0.04j),
   ]:
       tipper_data[period] = np.tile(
           np.array([[tx, ty]], dtype=complex),
           (n_stations, 1),
       )

   fig, axes = plot_induction_multiperiod_map(
       sites,
       periods=list(tipper_data),
       tipper_data=tipper_data,
       arrow_scale=6.0,
       show_background_cbar=False,
   )

The station order in each ``tipper_data`` array must match the station
order returned by ``ensure_sites`` for the input survey.

Spectra-Direct Workflows
------------------------

The spectra-direct helpers work before final EDI assembly.  They expect
``Spectra`` objects or dictionaries of spectra objects and recover the
tipper from spectral estimates.

Use these functions when your workflow is still at the spectra stage:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Function
     - Use
   * - ``plot_induction_map_from_spectra``
     - Draw real and imaginary induction arrows from one or more spectra
       objects.
   * - ``plot_tipper_polar_from_spectra``
     - Inspect one spectra object's tipper azimuth and magnitude.
   * - ``plot_induction_rose_from_spectra``
     - Summarize spectra-derived induction azimuths over a period band.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       plot_induction_map_from_spectra,
       plot_induction_rose_from_spectra,
       plot_tipper_polar_from_spectra,
   )

   spectra_by_station = {
       "S001": spectra_001,
       "S002": spectra_002,
       "S003": spectra_003,
   }

   coords = {
       "S001": (0.0, 0.0),
       "S002": (500.0, 0.0),
       "S003": (1000.0, 0.0),
   }

   plot_induction_map_from_spectra(
       spectra_by_station,
       coords=coords,
       period=100.0,
   )

   plot_tipper_polar_from_spectra(
       spectra_by_station["S001"],
       component="real",
   )

   plot_induction_rose_from_spectra(
       spectra_by_station,
       component="real",
       pband=(10.0, 1000.0),
   )

For spectra maps, ``coords`` are plot coordinates ``(x, y)``.  A bare
``Spectra`` object does not carry reliable map geometry.

Recommended Workflow
--------------------

A robust tipper interpretation keeps the raw station behavior, the
period behavior, and the sign convention visible:

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.emtools import (
       ensure_sites,
       plot_induction_convention,
       plot_induction_map,
       plot_induction_multiperiod_map,
       plot_induction_rose,
       plot_induction_section,
       plot_tipper_hodograms,
       plot_tipper_polar,
   )

   sites = ensure_sites(
       Path("data/MT/kap03lmt_edis"),
       recursive=True,
   )

   strongest_station = "kap151"
   periods = [25.0, 650.0, 2000.0, 17000.0]

   plot_tipper_hodograms(
       sites,
       station=strongest_station,
       bands=[(25.0, 200.0), (200.0, 2000.0), (2000.0, 17000.0)],
   )

   plot_tipper_polar(
       sites,
       station=strongest_station,
       component="real",
   )

   plot_induction_map(
       sites,
       period=2000.0,
       convention="park",
       show_real=True,
       show_imag=True,
       scale=4.0,
   )

   plot_induction_convention(
       sites,
       period=650.0,
       station_labels=False,
   )

   plot_induction_rose(
       sites,
       component="real",
       pband=(25.0, 200.0),
       title="Short-period induction azimuths",
   )

   plot_induction_rose(
       sites,
       component="real",
       pband=(2000.0, 17000.0),
       title="Long-period induction azimuths",
   )

   plot_induction_section(
       sites,
       component="abs",
       n_periods=30,
   )

   plot_induction_multiperiod_map(
       sites,
       periods=periods,
       convention="park",
       arrow_scale=6.0,
       show_background_cbar=False,
       station_labels=False,
   )

This sequence answers the practical questions in order: which station is
strong, whether its response is coherent, where the profile responds,
which convention is being used, whether azimuths tighten with period,
and how the anomaly migrates across period.

Common Pitfalls
---------------

Do not use tipper tools on surveys without vertical-field data and then
interpret ``"no tipper"`` as a geological result.

Always state the sign convention.  Parkinson and Wiese views are rotated
relative to each other.

Do not interpret index-based map axes as geographic distance.  If EDI
coordinates are missing, the plots may fall back to station index.

Do not collapse all periods too early.  A strong whole-band station may
be strong only over a narrow period window.

Do not treat synthetic or placeholder backgrounds as real topography in
multi-period maps.  Pass a real background raster for publication.

Worked Example
--------------

The gallery example below uses the KAP03 MT profile with real tipper
data.  It moves from station-level hodograms and polar plots to maps,
convention comparisons, roses, period sections, and a multi-period
publication-style map.

.. literalinclude:: ../../../examples/emtools/plot_tf.py
   :language: python
   :linenos:
