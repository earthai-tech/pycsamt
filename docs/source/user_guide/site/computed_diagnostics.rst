.. _site_computed_diagnostics:

Computed Diagnostics
====================

.. currentmodule:: pycsamt.site.compute

The :mod:`pycsamt.site.compute` module provides lightweight station
diagnostics for EDI-like sites and site collections. These functions do not
modify the input data. They read impedance, frequency, and tipper arrays and
return small numerical summaries that are useful during quality control,
survey preparation, and interpretation.

Use computed diagnostics before heavier processing when you need to answer
questions such as:

* does each station have a usable impedance tensor?
* what strike angle is suggested by the off-diagonal tensor structure?
* what apparent resistivity is observed near a frequency of interest?
* does phase vary smoothly across the band used for inversion?
* is the tipper response small, large, missing, or spatially variable?

Input Contract
--------------

The diagnostics accept either one EDI-like object or an iterable of EDI-like
objects. In practice this means:

* a :class:`pycsamt.seg.edi.EDIFile`;
* a :class:`pycsamt.site.base.Site`;
* a :class:`pycsamt.site.base.Sites` collection;
* a list, tuple, or other iterable of EDI-like objects.

The object must expose the data needed by the selected diagnostic:

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Diagnostic
     - Required arrays
     - Missing-data behaviour
   * - :func:`strike_estimate`
     - Frequency vector and impedance tensor shaped ``(n_freq, 2, 2)``.
     - Returns ``NaN`` when tensor or frequency data are missing.
   * - :func:`res_at_freq`
     - Frequency vector and impedance tensor shaped ``(n_freq, 2, 2)``.
     - Returns ``NaN`` resistivity values when required data are missing.
   * - :func:`phase_slope`
     - Frequency vector and impedance tensor shaped ``(n_freq, 2, 2)``.
     - Returns ``NaN`` slopes when the requested band is empty or unusable.
   * - :func:`tipper_magnitude`
     - Frequency vector and tipper array shaped ``(n_freq, 2)``.
     - Returns ``NaN`` summary values, or an empty per-frequency result.

For deterministic diagnostics on incomplete data, prepare the site first with
:func:`pycsamt.site.edit.fill_missing` or filter the survey with
:func:`pycsamt.site.selection.keep_finite_z`.

Return Types
------------

The module follows one consistent rule:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Input
     - Return value
   * - Single site
     - A scalar or dictionary.
   * - Iterable of sites
     - A :class:`pandas.DataFrame`.
   * - Iterable of sites with ``api=True``
     - A pyCSAMT ``APIFrame`` wrapper around the DataFrame, when API view
       support is enabled.

This makes the functions convenient both in notebooks and in automated
pipeline steps.

.. code-block:: python
   :linenos:

   from pycsamt.seg.edi import EDIFile
   from pycsamt.site.compute import res_at_freq, strike_estimate

   site = EDIFile("data/edi/S01.edi")

   theta = strike_estimate(site)
   rho = res_at_freq(site, 100.0)

   print(theta)
   print(rho["res_xy"], rho["res_yx"], rho["f_used"])

For a collection:

.. code-block:: python
   :linenos:

   from pycsamt.site import Sites
   from pycsamt.site.compute import phase_slope

   sites = Sites.from_path("data/edi")
   table = phase_slope(sites, band=(1.0, 1000.0))

   print(table.head())

Diagnostic Map
--------------

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Function
     - Output
     - Main use
   * - :func:`strike_estimate`
     - ``theta_deg`` in degrees.
     - Quick geoelectric strike estimate for tensor rotation checks and 2-D
       assumption screening.
   * - :func:`res_at_freq`
     - ``res_xy``, ``res_yx``, and ``f_used``.
     - Compare stations at one frequency, select target frequency slices, or
       build compact QC tables.
   * - :func:`phase_slope`
     - ``slope_xy`` and ``slope_yx`` in degrees per frequency decade.
     - Detect abrupt phase behaviour, band-edge problems, or unstable curves.
   * - :func:`tipper_magnitude`
     - Mean, median, max, or per-frequency tipper magnitude.
     - Evaluate vertical magnetic transfer response and possible 3-D
       structure indicators.

Strike Estimate
---------------

:func:`strike_estimate` estimates a 2-D geoelectric strike angle from the
impedance tensor. It supports three method names:

``"swift"``
   Grid search from 0 to 179 degrees. The selected angle minimizes diagonal
   tensor power after rotation.

``"groom"``
   Currently an alias for ``"swift"`` in this lightweight diagnostic layer.

``"phase_diff"``
   A coarse fallback returning either 0 or 90 degrees based on the relative
   median magnitude of :math:`Z_{xy}` and :math:`Z_{yx}`.

For the Swift-style criterion, the tensor is rotated by a trial angle
:math:`\theta`, then the diagonal power is summarized:

.. math::

   J(\theta) =
   \operatorname{median}
   \left(
      |Z'_{xx}|^2 + |Z'_{yy}|^2
   \right).

The reported strike is the angle that minimizes :math:`J(\theta)` over the
1-degree search grid.

.. code-block:: python
   :linenos:

   from pycsamt.site import Sites
   from pycsamt.site.compute import strike_estimate
   from pycsamt.site.edit import fill_missing

   sites = Sites.from_path("data/edi")

   # Optional: make missing tensor entries explicit before diagnostics.
   prepared = [
       fill_missing(site, how="zero", components=("Z",), inplace=False)
       for site in sites
   ]

   strike_table = strike_estimate(prepared, method="swift")
   print(strike_table[["station", "theta_deg"]])

Interpretation notes:

* strike estimates are screening values, not final structural interpretation;
* unstable or sparse tensors can produce misleading angles;
* the result is in the interval :math:`[0, 180)`;
* compare strike values across neighbouring stations before rotating an entire
  survey line;
* use :func:`pycsamt.site.edit.rotate` or
  :func:`pycsamt.site.edit.rotate_all` after deciding on a rotation angle.

Apparent Resistivity At One Frequency
-------------------------------------

:func:`res_at_freq` evaluates apparent resistivity for the off-diagonal
components :math:`Z_{xy}` and :math:`Z_{yx}` at a requested frequency.

The apparent resistivity is computed as:

.. math::

   \rho_a =
   \frac{|Z|^2}{\mu_0\,2\pi f},

where :math:`Z` is the selected complex impedance component, :math:`\mu_0` is
magnetic permeability, and :math:`f` is frequency in Hz.

Two frequency-selection modes are available:

``how="nearest"``
   Select the nearest available native frequency and report that value in
   ``f_used``.

``how="interp"``
   Compute resistivity at all native frequencies, then interpolate to the
   requested frequency using linear interpolation in frequency.

.. code-block:: python
   :linenos:

   from pycsamt.site import Sites
   from pycsamt.site.compute import res_at_freq

   sites = Sites.from_path("data/edi")

   near = res_at_freq(sites, 150.0, how="nearest")
   interp = res_at_freq(sites, 150.0, how="interp")

   print(near[["station", "res_xy", "res_yx", "f_used"]])
   print(interp[["station", "res_xy", "res_yx", "f_used"]])

Use ``nearest`` when you want to preserve the exact sampled frequency axis.
Use ``interp`` when you need one common comparison frequency across stations
whose frequency grids are slightly different.

Phase Slope
-----------

:func:`phase_slope` summarizes how phase changes across a frequency band. For
each site, it computes phase for the off-diagonal components:

.. math::

   \phi(f) = \arg(Z(f)) \times 180 / \pi,

then fits a straight line against log-frequency:

.. math::

   \phi(f) \approx a \log_{10}(f) + b.

The reported slope :math:`a` is measured in degrees per decade. The function
returns one slope for :math:`Z_{xy}` and one for :math:`Z_{yx}`.

.. code-block:: python
   :linenos:

   from pycsamt.site import Sites
   from pycsamt.site.compute import phase_slope

   sites = Sites.from_path("data/edi")

   slopes = phase_slope(sites, band=(1.0, 1000.0))
   steep = slopes[
       slopes["slope_xy"].abs().gt(30.0)
       | slopes["slope_yx"].abs().gt(30.0)
   ]

   print(steep)

Phase-slope diagnostics are useful for finding:

* phase curves that change abruptly within the inversion band;
* stations with unstable tensor estimates;
* band selections that cross noisy frequency ranges;
* components that may require visual inspection before inversion.

The function does not unwrap phase. If phase wraps are important for your
dataset, inspect the curves directly before treating the slope as a physical
trend.

Tipper Magnitude
----------------

:func:`tipper_magnitude` computes the magnitude of the tipper vector:

.. math::

   \|\mathbf{T}\| =
   \sqrt{|T_x|^2 + |T_y|^2}.

By default it returns summary statistics. Set ``per_freq=True`` to return one
row per frequency.

.. code-block:: python
   :linenos:

   from pycsamt.site import Sites
   from pycsamt.site.compute import tipper_magnitude

   sites = Sites.from_path("data/edi")

   summary = tipper_magnitude(sites)
   long_table = tipper_magnitude(sites, per_freq=True)

   print(summary[["station", "mean", "median", "max"]])
   print(long_table.head())

Summary mode returns:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Column
     - Meaning
   * - ``mean``
     - Average tipper magnitude over available frequencies.
   * - ``median``
     - Median tipper magnitude, less sensitive to isolated spikes.
   * - ``max``
     - Maximum tipper magnitude, useful for finding extreme response bands.

Per-frequency mode returns ``station``, ``freq``, and ``mag`` columns for
collection inputs. For one site, it returns a dictionary containing ``freq``
and ``mag`` arrays.

APIFrame Output
---------------

For collection inputs, pass ``api=True`` when the result should carry pyCSAMT
API metadata in addition to tabular values.

.. code-block:: python
   :linenos:

   from pycsamt.site import Sites
   from pycsamt.site.compute import (
       phase_slope,
       res_at_freq,
       strike_estimate,
       tipper_magnitude,
   )

   sites = Sites.from_path("data/edi")

   strike = strike_estimate(sites, api=True)
   rho = res_at_freq(sites, 100.0, api=True)
   slopes = phase_slope(sites, (1.0, 1000.0), api=True)
   tipper = tipper_magnitude(sites, api=True)

   print(strike.kind)
   print(rho.kind)
   print(slopes.kind)
   print(tipper.kind)

This is useful when diagnostics are emitted by CLI commands, agents, or
pipelines that preserve result provenance.

Quality-Control Workflow
------------------------

The following example combines selection, editing, and computed diagnostics.

.. code-block:: python
   :linenos:

   from pycsamt.site import Sites
   from pycsamt.site.compute import (
       phase_slope,
       res_at_freq,
       strike_estimate,
       tipper_magnitude,
   )
   from pycsamt.site.edit import fill_missing
   from pycsamt.site.selection import by_freq, keep_finite_z

   sites = Sites.from_path("data/edi")
   sites = keep_finite_z(sites)
   sites = by_freq(sites, fmin=1.0, fmax=1000.0)

   prepared = [
       fill_missing(site, how="zero", components=("Z",), inplace=False)
       for site in sites
   ]

   strike = strike_estimate(prepared)
   rho100 = res_at_freq(prepared, 100.0, how="nearest")
   slopes = phase_slope(prepared, band=(1.0, 1000.0))
   tipper = tipper_magnitude(prepared)

   qc = (
       strike
       .merge(rho100, on="station", how="outer")
       .merge(slopes, on="station", how="outer")
       .merge(tipper, on="station", how="outer")
   )

   print(qc)

Common Mistakes
---------------

Using diagnostics as final interpretation
   These functions are quick screening tools. Confirm important decisions with
   maps, pseudo-sections, tensor plots, and inversion sensitivity tests.

Mixing frequency axes without checking ``f_used``
   In ``nearest`` mode, different stations may use different native
   frequencies. Always inspect ``f_used`` before comparing values.

Ignoring missing arrays
   Missing Z or tipper data are reported as ``NaN`` rather than raising hard
   failures. Filter or fill intentionally before trusting the table.

Over-reading phase slope
   A large phase slope may indicate structure, noise, phase wrapping, or a bad
   frequency band. Use it as a prompt for inspection, not as a standalone
   classifier.

Next Pages
----------

Continue with:

* :doc:`selection` for selecting stations before diagnostics;
* :doc:`editing` for preparing tensors, frequencies, names, and coordinates;
* :doc:`export_reporting` for writing diagnostic summaries into survey
  deliverables;
* :doc:`../theory/impedance_tensor` for the physical meaning of impedance,
  apparent resistivity, phase, and tensor components.

