.. _site_selection:

Site Selection
==============

.. currentmodule:: pycsamt.site.selection

The selection helpers choose which stations remain in a survey before
editing, diagnostics, export, modelling, or inversion preparation. They are
small, stable filters: each helper accepts a site-like input, walks through
the EDI-like items in their existing order, and returns a new
:class:`pycsamt.site.base.Sites` container with the retained stations.

Use this page when you need to:

* keep stations by name, index, chainage, frequency coverage, or location;
* remove empty or unusable stations before expensive processing;
* apply a quick quality gate based on finite impedance or phase error;
* combine several simple selectors into a reproducible preparation workflow;
* decide when to use functional selectors versus
  :meth:`pycsamt.site.base.Sites.select`.

Selector Map
------------

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Function
     - Keeps stations when
     - Typical use
   * - :func:`by_names`
     - The station name matches one or more patterns.
     - Keep a named line, block, or manually reviewed station list.
   * - :func:`by_index`
     - The station position appears in an index list.
     - Keep the first, last, or manually inspected positions.
   * - :func:`by_chainage`
     - Stored chainage lies inside a closed interval.
     - Trim a profile to a modelled distance window.
   * - :func:`by_freq`
     - At least one finite frequency lies inside a closed interval.
     - Keep stations that cover a target period or frequency band.
   * - :func:`by_bbox`
     - Latitude and longitude fall inside a geographic bounding box.
     - Keep stations inside a field area or map tile.
   * - :func:`by_predicate`
     - A custom Python function returns ``True``.
     - Express project-specific filters.
   * - :func:`keep_finite_z`
     - Impedance or resistivity arrays contain finite values.
     - Drop invalid placeholders before diagnostics or inversion.
   * - :func:`mask_large_phase_err`
     - The maximum finite phase error is below a threshold.
     - Remove stations with very uncertain phase products.
   * - :func:`drop_empty`
     - Frequency and impedance-like data are structurally present.
     - Drop stations with no usable data arrays.

Selection Contract
------------------

All functional selectors follow the same broad contract.

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.selection import by_names, keep_finite_z

   collection = EDICollection.from_sources("data/edi")

   selected = by_names(collection, "L01_*")
   selected = keep_finite_z(selected)

   print(type(selected).__name__)
   print([site.name for site in selected])

The important guarantees are:

* input can be a :class:`pycsamt.site.base.Sites`, an EDI collection, a list of
  EDI files, or any iterable accepted by the site iterator;
* output is a new :class:`pycsamt.site.base.Sites` wrapper;
* relative order of retained stations follows the input order;
* duplicate index requests are collapsed because selection is membership
  based;
* the original container is not modified by the selection itself.

Selectors keep the raw EDI objects inside the returned ``Sites`` wrapper. This
means later edits can still act on the original EDI structures if an editing
function is called with ``inplace=True``.

Name Selection
--------------

:func:`by_names` keeps stations whose resolved station name matches one or
more patterns. Station names are resolved through
:func:`pycsamt.site.utils.station_name`, so names are taken from the same
header/object fields used by the rest of the site tools.

Supported pattern types are:

``str``
   Literal names, glob-like patterns such as ``"S*"`` and ``"A??"``, and
   regex-looking strings.

``re.Pattern``
   A compiled regular expression. Use this when you need strict regular
   expression control.

``callable``
   A function receiving the station name and returning ``True`` or ``False``.

``list`` or ``tuple``
   Any mix of the above. A station is kept if any pattern matches.

.. code-block:: python
   :linenos:

   import re

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.selection import by_names

   sites = EDICollection.from_sources("data/edi")

   line_a = by_names(sites, "A*")
   stations_01_to_05 = by_names(sites, re.compile(r"^S0[1-5]$"))
   reviewed = by_names(sites, ["S01", "S05", "S09"])
   ending_with_a = by_names(sites, lambda name: name.endswith("A"))

String matching uses the shared pyCSAMT name matcher, which is designed for
case-insensitive station matching. For strict case-sensitive regex behavior,
pass a compiled regular expression with the flags you want.

Index Selection
---------------

:func:`by_index` keeps stations by zero-based position. Negative indices follow
normal Python sequence rules, so ``-1`` means the last station and ``-2`` means
the station before it.

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.selection import by_index

   sites = EDICollection.from_sources("data/edi")

   first_and_last = by_index(sites, [0, -1])
   second = by_index(sites, 1)

   print([site.name for site in first_and_last])
   print([site.name for site in second])

Invalid indices are ignored. If all requested indices are invalid, the result
is an empty ``Sites`` container.

The output order follows the original survey order, not the order of the
``indices`` argument:

.. code-block:: python
   :linenos:

   subset = by_index(sites, [-1, 0])

   # The first station still appears before the last station in the result.
   print([site.name for site in subset])

Chainage Selection
------------------

:func:`by_chainage` keeps stations whose chainage lies in a closed interval
:math:`s_{min} \le s \le s_{max}`. The selector first looks for
``HEAD.chainage`` and then falls back to ``edi.chainage``.

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.selection import by_chainage

   sites = EDICollection.from_sources("data/edi")

   model_window = by_chainage(
       sites,
       smin=1_000.0,
       smax=8_000.0,
   )

   print([site.name for site in model_window])

Stations without numeric chainage are skipped. This is useful after a profile
has already been computed and chainage values have been written onto the EDI
headers or attached to the EDI objects.

If you have coordinates but no stored chainage yet, build a profile first.

.. code-block:: python
   :linenos:

   from pycsamt.site.profile import Profile
   from pycsamt.site.selection import by_chainage

   profile = Profile.from_sites(sites)

   for ed in sites:
       # Store chainage on the raw EDI object for later selectors.
       station = getattr(ed, "station", "")
       ed.chainage = profile.chainages.get(station)

   middle = by_chainage(sites, 2_000.0, 5_000.0)

Frequency Coverage Selection
----------------------------

:func:`by_freq` keeps a station if at least one finite frequency value lies
inside a closed interval.

.. math::

   f_{min} \le f_i \le f_{max}

This is a station-level filter. It does not remove frequency rows from each
station.

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.selection import by_freq

   sites = EDICollection.from_sources("data/edi")

   broadband = by_freq(sites, fmin=0.1, fmax=1_000.0)
   low_frequency = by_freq(sites, fmin=0.1, fmax=10.0)

   print(len(broadband), len(low_frequency))

To actually subset the rows of frequency-indexed arrays, use
:func:`pycsamt.site.edit.select_freq` for one site or
:func:`pycsamt.site.edit.select_freq_all` for a collection.

.. code-block:: python
   :linenos:

   from pycsamt.site.edit import select_freq_all
   from pycsamt.site.selection import by_freq

   sites = by_freq(sites, fmin=0.1, fmax=100.0)
   sliced = select_freq_all(
       sites,
       fmin=0.1,
       fmax=100.0,
       inplace=False,
   )

Bounding Box Selection
----------------------

:func:`by_bbox` keeps stations inside an inclusive latitude/longitude box.

.. math::

   minlat \le lat \le maxlat

.. math::

   minlon \le lon \le maxlon

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.selection import by_bbox

   sites = EDICollection.from_sources("data/edi")

   field_area = by_bbox(
       sites,
       minlat=34.5,
       minlon=12.0,
       maxlat=35.5,
       maxlon=13.0,
   )

Coordinates are read through :func:`pycsamt.site.utils.get_coords`. Stations
with missing or non-finite coordinates are skipped.

This selector is intentionally simple: it performs an axis-aligned test in
latitude and longitude degrees. It does not project coordinates and it does
not handle antimeridian wrapping. For projected distance or profile filters,
use the profile and location tools described in :doc:`location_profile`.

Custom Predicate Selection
--------------------------

:func:`by_predicate` keeps stations for which a user-provided function returns
``True``. The predicate receives the raw EDI-like object, not a
:class:`pycsamt.site.base.Site` wrapper.

.. code-block:: python
   :linenos:

   from pycsamt.site.base import Site
   from pycsamt.site.selection import by_predicate

   stations_with_tipper = by_predicate(
       sites,
       lambda ed: Site(ed).has_component("tipper"),
   )

   enough_frequencies = by_predicate(
       sites,
       lambda ed: len(Site(ed).freq) >= 8,
   )

Predicate exceptions are caught and treated as ``False``. This keeps survey
filters robust when one station has an unusual structure.

Use ``by_predicate`` for project logic that is too specific for the standard
selectors:

.. code-block:: python
   :linenos:

   from pycsamt.site.utils import station_name

   def keep_production_station(ed):
       name = station_name(ed)
       return (
           name.startswith("P")
           and not name.endswith("_TEST")
       )

   production = by_predicate(sites, keep_production_station)

Finite Impedance Selection
--------------------------

:func:`keep_finite_z` keeps stations that contain at least one finite
impedance value. If impedance values are unavailable, it falls back to common
resistivity array names such as ``_resistivity``, ``resistivity``, or ``rho``.

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.selection import keep_finite_z

   sites = EDICollection.from_sources("data/edi")
   valid = keep_finite_z(sites)

   print(f"{len(valid)} stations contain finite impedance data")

This is the selector to run before diagnostics that require real impedance
content. It does not inspect phase errors, frequency coverage, or coordinate
quality.

Phase Error Selection
---------------------

:func:`mask_large_phase_err` filters whole stations using their maximum finite
phase-error value. Despite the function name, the current implementation does
not mask individual data rows. A station is kept when:

* no phase-error array is found;
* the phase-error array is empty or entirely non-finite;
* the maximum finite phase error is less than or equal to ``thresh``.

.. code-block:: python
   :linenos:

   from pycsamt.site.selection import mask_large_phase_err

   conservative = mask_large_phase_err(sites, thresh=10.0)

   print([site.name for site in conservative])

Common phase-error attributes are checked on the ``Z`` container, including
``_phase_err`` and ``phase_err``.

Because stations without stored phase-error arrays are kept, this selector is
best used as a quality gate after a processing step that actually produced
phase-error products. If absence of error estimates should be considered a
failure in your workflow, combine this selector with a stricter predicate.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.site.selection import by_predicate, mask_large_phase_err

   def has_phase_error(ed):
       z = getattr(ed, "Z", None)
       arr = getattr(z, "_phase_err", None) if z is not None else None
       return arr is not None and np.asarray(arr).size > 0

   with_errors = by_predicate(sites, has_phase_error)
   clean = mask_large_phase_err(with_errors, thresh=10.0)

Empty Site Selection
--------------------

:func:`drop_empty` removes stations that are structurally empty. A station is
considered empty when:

* the frequency vector is missing or has zero length;
* the ``Z`` container is missing;
* the ``Z`` container has no impedance array and no resistivity surrogate;
* available arrays are empty, entirely non-finite, or only huge sentinel
  values.

.. code-block:: python
   :linenos:

   from pycsamt.site.selection import drop_empty, keep_finite_z

   non_empty = drop_empty(sites)
   finite = keep_finite_z(non_empty)

This is usually the first cleanup selector in a workflow, because it removes
stations that cannot participate in most downstream computations.

Functional Selectors Versus Sites.select
----------------------------------------

The :class:`pycsamt.site.base.Sites` container also exposes a compact
:meth:`pycsamt.site.base.Sites.select` method. Use it for simple name lists or
wrapper-based predicates.

.. code-block:: python
   :linenos:

   from pycsamt.site.base import Sites

   sites = Sites.from_any("data/edi")

   reviewed = sites.select(names=["S01", "S03"])
   with_zxy = sites.select(
       predicate=lambda site: site.has_component("Zxy")
   )

Use the functions in :mod:`pycsamt.site.selection` when you need richer
matching rules, index handling, chainage filtering, frequency coverage,
coordinate boxes, or data-quality filters.

Combining Selectors
-------------------

Selectors compose naturally because each one returns a ``Sites`` object.
Order the pipeline from cheap structural filters to more specific project
filters.

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.edit import select_freq_all
   from pycsamt.site.selection import (
       by_bbox,
       by_freq,
       by_names,
       drop_empty,
       keep_finite_z,
       mask_large_phase_err,
   )

   sites = EDICollection.from_sources("data/edi")

   sites = drop_empty(sites)
   sites = keep_finite_z(sites)
   sites = by_names(sites, ["L01_*", "L02_*"])
   sites = by_bbox(
       sites,
       minlat=34.5,
       minlon=12.0,
       maxlat=35.5,
       maxlon=13.0,
   )
   sites = by_freq(sites, fmin=0.1, fmax=100.0)
   sites = mask_large_phase_err(sites, thresh=10.0)

   # Now modify the frequency rows, after selecting stations that overlap
   # the requested band.
   prepared = select_freq_all(
       sites,
       fmin=0.1,
       fmax=100.0,
       inplace=False,
   )

This pattern separates station selection from data editing. The selectors
decide which stations remain; editing functions decide how station contents
change.

Selection Before Inversion
--------------------------

Before preparing files for Occam2D, ModEM, or MARE2DEM, a selection workflow
should normally answer four questions:

1. Does every retained station contain usable impedance data?
2. Does every retained station overlap the target frequency band?
3. Are the retained stations inside the intended profile or map area?
4. Are station errors acceptable for the planned inversion weighting?

For a 2-D profile workflow, combine selectors with the profile tools:

.. code-block:: python
   :linenos:

   from pycsamt.site.profile import Profile
   from pycsamt.site.selection import (
       by_chainage,
       by_freq,
       drop_empty,
       keep_finite_z,
   )

   sites = drop_empty(sites)
   sites = keep_finite_z(sites)
   sites = by_freq(sites, 0.1, 100.0)

   profile = Profile.from_sites(sites)

   for ed in sites:
       name = getattr(ed, "station", "")
       ed.chainage = profile.chainages.get(name)

   model_sites = by_chainage(sites, 500.0, 9_500.0)
   ordered = profile.sort_sites(model_sites)

For a 3-D or map-based workflow, prefer geographic selection and then export
or modelling preparation:

.. code-block:: python
   :linenos:

   from pycsamt.site.selection import by_bbox, by_freq, keep_finite_z

   cube_sites = keep_finite_z(sites)
   cube_sites = by_bbox(cube_sites, 34.5, 12.0, 35.5, 13.0)
   cube_sites = by_freq(cube_sites, 0.01, 1_000.0)

Common Mistakes
---------------

Using :func:`by_freq` as if it sliced rows
   :func:`by_freq` keeps or drops stations. Use
   :func:`pycsamt.site.edit.select_freq_all` to change frequency rows.

Expecting :func:`mask_large_phase_err` to mask cells
   The current function filters whole stations. Use editing tools or a custom
   predicate when you need row-level masking.

Running :func:`by_chainage` before chainage exists
   Chainage must already be stored on ``HEAD.chainage`` or ``edi.chainage``.
   Build a profile first when you only have coordinates.

Assuming :func:`by_bbox` projects coordinates
   It works directly in latitude/longitude degrees. Use projection helpers
   from :doc:`location_profile` for projected distances.

Using strict name case expectations with string patterns
   String matching is intentionally tolerant. Use compiled regular
   expressions when strict case behavior matters.

Next Pages
----------

* :doc:`containers` for the ``Site`` and ``Sites`` wrappers;
* :doc:`editing` for changing station names, coordinates, frequency rows, and
  data arrays after selection;
* :doc:`location_profile` for chainage, profile geometry, distance, and
  coordinate projection;
* :doc:`computed_diagnostics` for station-level diagnostics that often run
  after selection.
