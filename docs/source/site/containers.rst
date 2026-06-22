.. _site_containers:

Site Containers
===============

.. currentmodule:: pycsamt.site.base

The :mod:`pycsamt.site.base` module provides the station-centric containers
used throughout the site layer. These containers wrap SEG-EDI objects without
replacing the lower-level EDI parser. Their purpose is to make common survey
operations easier: station lookup, coordinate access, tensor inspection,
DataFrame export, bulk selection, profile preparation, and writing cleaned
sites back to disk.

The main objects are:

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Object
     - Scope
     - Main use
   * - :class:`SiteMixin`
     - One EDI-like site.
     - Shared accessors for name, coordinates, frequency, impedance,
       resistivity, phase, tipper, metadata, quality flags, and DataFrame
       export.
   * - :class:`Site`
     - One concrete station.
     - High-level wrapper around one :class:`pycsamt.seg.edi.EDIFile` with a
       stable station identity and convenience edit helpers.
   * - :class:`Sites`
     - Many stations.
     - Ordered collection wrapper with integer indexing, name lookup,
       selection, batch edits, topography alignment, profile conversion, and
       export.
   * - :func:`to_sites`
     - API boundary helper.
     - Coerce compatible EDI-like inputs into a :class:`Sites` container.
   * - :func:`to_edis`
     - API boundary helper.
     - Unwrap :class:`Site` or :class:`Sites` objects back to raw EDI objects
       or an :class:`pycsamt.seg.collection.EDICollection`.

Where Containers Fit
--------------------

The site containers sit above :mod:`pycsamt.seg` and below the higher-level
workflow layers.

.. list-table::
   :header-rows: 1
   :widths: 26 37 37

   * - Layer
     - Main object
     - Responsibility
   * - EDI parsing
     - :class:`pycsamt.seg.edi.EDIFile`
     - Read, store, and write one SEG-EDI file.
   * - EDI collection
     - :class:`pycsamt.seg.collection.EDICollection`
     - Discover and parse many EDI files from paths, folders, or glob
       patterns.
   * - Site layer
     - :class:`Site`, :class:`Sites`
     - Expose station-oriented accessors, editing, filtering, summaries, and
       survey preparation helpers.
   * - Downstream workflows
     - Selection, diagnostics, export, profiles, inversion, agents.
     - Use prepared site containers as consistent inputs.

Use :class:`EDIFile` or :class:`EDICollection` when you are primarily parsing
SEG-EDI data. Use :class:`Site` or :class:`Sites` when you are preparing,
inspecting, selecting, or editing station data.

Creating One Site
-----------------

Create a :class:`Site` from a parsed :class:`pycsamt.seg.edi.EDIFile`.

.. code-block:: python
   :linenos:

   from pycsamt.seg.edi import EDIFile
   from pycsamt.site.base import Site

   edi = EDIFile("data/edi/S01.edi")
   site = Site(edi)

   print(site.name)
   print(site.coords)
   print(site.summary())

The wrapper keeps the original EDI object available as ``site.edi`` and
through :meth:`Site.to_edi`. Prefer the :class:`Site` accessors for routine
work, and unwrap to EDI only when you need lower-level EDI functionality.

.. code-block:: python
   :linenos:

   raw = site.to_edi()
   detached = site.to_edi(copy=True)

Station Identity
----------------

Station identity is normalized from EDI ``HEAD`` fields. The name resolution
order is:

#. ``dataid``;
#. ``station``;
#. ``sitename``;
#. ``name``;
#. ``STATION``;
#. the file stem, when header fields are missing.

When a :class:`Site` is created, the constructor attempts to stabilize the EDI
header so that ``dataid`` and, when absent, ``station`` match the resolved
site stem. This makes name lookup and downstream joins more predictable.

.. code-block:: python
   :linenos:

   site = Site(EDIFile("data/edi/LINE01_012.edi"))

   print(site.name)
   print(site.meta)

Rename a site with :meth:`Site.rename`:

.. code-block:: python
   :linenos:

   renamed = site.rename("L01_012")

   print(site.name)      # original unchanged
   print(renamed.name)   # new wrapper

Use ``inplace=True`` when the existing wrapper should be modified:

.. code-block:: python
   :linenos:

   site.rename("L01_012", inplace=True)

Coordinates
-----------

Coordinates are exposed as ``(lat, lon, elev)`` in decimal degrees and meters.

.. code-block:: python
   :linenos:

   lat, lon, elev = site.coords
   print(lat, lon, elev)

Update coordinates with :meth:`Site.set_coords`:

.. code-block:: python
   :linenos:

   moved = site.set_coords(10.25, 20.75, 640.0)
   print(moved.coords)

Like :meth:`Site.rename`, coordinate updates are copy-returning by default and
in-place when requested:

.. code-block:: python
   :linenos:

   site.set_coords(10.25, 20.75, 640.0, inplace=True)

Array Accessors
---------------

The container exposes the most common station arrays as properties.

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Property
     - Meaning
     - Shape convention
   * - ``freq``
     - Frequency vector in Hz.
     - ``(n_freq,)``.
   * - ``z``
     - Complex impedance tensor.
     - Usually ``(n_freq, 2, 2)``; flattened internally as
       ``Zxx, Zxy, Zyx, Zyy`` when exported.
   * - ``z_err``
     - Impedance uncertainty.
     - Aligned with ``z`` when present.
   * - ``rho``
     - Apparent resistivity.
     - Aligned with impedance components.
   * - ``phase``
     - Impedance phase in degrees.
     - Aligned with impedance components.
   * - ``tipper``
     - Vertical magnetic transfer function.
     - ``(n_freq, 2)`` with ``Tx`` and ``Ty``.
   * - ``meta``
     - Header and optional ``INFO`` metadata.
     - Python dictionary.

Example:

.. code-block:: python
   :linenos:

   print(site.freq)
   print(site.z.shape if site.z is not None else "no Z")
   print(site.quality_flags())

Quality And Component Checks
----------------------------

Use :meth:`Site.quality_flags` for coarse array availability checks.

.. code-block:: python
   :linenos:

   flags = site.quality_flags()

   if not flags["has_z"]:
       print(f"{site.name} has no finite impedance tensor")

Use :meth:`Site.has_component` when a specific impedance or tipper component
matters.

.. code-block:: python
   :linenos:

   if site.has_component("Zxy") and site.has_component("Zyx"):
       print("off-diagonal impedance components are available")

   if site.has_component("tipper"):
       print("tipper is available")

Supported impedance component names are ``Zxx``, ``Zxy``, ``Zyx``, and
``Zyy``. Tipper checks accept ``tip``, ``tx``, ``ty``, and ``tipper``.

DataFrame Export
----------------

:meth:`Site.to_dataframe` exports core arrays into frequency-indexed
:class:`pandas.DataFrame` objects. This is the quickest way to inspect one
station in a notebook or build masks for editing.

.. code-block:: python
   :linenos:

   z = site.to_dataframe("z")
   rp = site.to_dataframe("resphase")
   tip = site.to_dataframe("tipper")

   print(z.head())
   print(rp.head())
   print(tip.head())

The accepted ``kind`` values are:

.. list-table::
   :header-rows: 1
   :widths: 24 28 48

   * - Kind
     - Aliases
     - Columns
   * - ``"z"``
     - ``"imp"``, ``"impedance"``
     - ``Zxx``, ``Zxy``, ``Zyx``, ``Zyy``.
   * - ``"resphase"``
     - ``"rp"``, ``"rho_phase"``
     - ``rho_zxx``, ``phi_zxx``, ``rho_zxy``, ``phi_zxy``, and the
       equivalent ``Zyx``/``Zyy`` columns.
   * - ``"tipper"``
     - ``"tip"``, ``"t"``
     - ``Tx``, ``Ty``.

Pass ``api=True`` when a pyCSAMT ``APIFrame`` is required:

.. code-block:: python
   :linenos:

   frame = site.to_dataframe("z", api=True)
   print(frame.kind)
   print(frame.df.head())

Creating A Collection
---------------------

Use :class:`Sites` when working with more than one station.

.. code-block:: python
   :linenos:

   from pycsamt.seg.edi import EDIFile
   from pycsamt.site.base import Sites

   edis = [
       EDIFile("data/edi/S01.edi"),
       EDIFile("data/edi/S02.edi"),
       EDIFile("data/edi/S03.edi"),
   ]

   sites = Sites(edis)

   print(len(sites))
   print([site.name for site in sites])

If your input is a directory or glob pattern, parse it first with
:class:`pycsamt.seg.collection.EDICollection`, then wrap it:

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.base import Sites

   collection = EDICollection.from_sources("data/edi")
   sites = Sites(collection)

For API boundaries that may receive several input types, use :func:`to_sites`:

.. code-block:: python
   :linenos:

   from pycsamt.site.base import to_sites

   sites = to_sites(collection)

Unwrapping Back To EDI
----------------------

The site layer is intentionally reversible. Use :func:`to_edis` when a
workflow has prepared, filtered, renamed, or edited :class:`Site` objects but
the next function expects raw EDI objects.

.. code-block:: python
   :linenos:

   from pycsamt.site.base import to_edis

   raw = to_edis(site)
   edis = to_edis(sites)
   collection = to_edis(sites, as_collection=True)

By default, unwrapping is shallow: the returned EDI objects are the same
objects stored by the wrappers. Pass ``copy=True`` when the caller should be
able to mutate the returned objects independently.

.. code-block:: python
   :linenos:

   detached_edis = to_edis(sites, copy=True)

The helper also accepts mixed inputs, path-like sources, raw EDI objects, and
existing :class:`~pycsamt.seg.collection.EDICollection` instances. Invalid
items are skipped unless ``strict=True`` is requested.

.. code-block:: python
   :linenos:

   edis = to_edis([site, raw_edi, "data/edi"], strict=False)
   checked = to_edis(sites, strict=True, progress=True)

Use the collection methods when you are already holding a :class:`Sites`
instance:

.. code-block:: python
   :linenos:

   edis = sites.to_edis()
   collection = sites.to_edicollection(copy=True)

``Sites.from_any`` is another entry point when you want pyCSAMT's session
normalization to interpret heterogeneous inputs:

.. code-block:: python
   :linenos:

   sites = Sites.from_any("data/edi")

Collection Lookup
-----------------

:class:`Sites` preserves input order and supports integer indexing,
case-insensitive name lookup, and safe lookup.

.. code-block:: python
   :linenos:

   first = sites[0]
   same = sites.by_index(0)
   station = sites["s02"]
   maybe = sites.get("missing")

   print(first.name)
   print(same.name)
   print(station.name)
   print(maybe is None)

``sites["missing"]`` raises ``KeyError``. Use :meth:`Sites.get` when missing
stations are acceptable.

The :meth:`Sites.as_list` method returns the underlying EDI objects:

.. code-block:: python
   :linenos:

   edi_objects = sites.as_list()

This is useful when passing a prepared collection into utilities that still
operate directly on EDI objects. For new API boundaries, prefer
:meth:`Sites.to_edis` or :func:`to_edis` because they also support copying,
progress display, strict validation, and collection output.

Mapping And Selection
---------------------

Use :meth:`Sites.map` to apply a lightweight operation to every station:

.. code-block:: python
   :linenos:

   names = sites.map(lambda site: site.name)
   coverage = sites.map(lambda site: site.summary()["nfreq"])

For simple collection filtering, use :meth:`Sites.select`:

.. code-block:: python
   :linenos:

   subset = sites.select(names=["S01", "S03"])
   with_zxy = sites.select(predicate=lambda site: site.has_component("Zxy"))

For richer selection by glob names, frequency coverage, coordinates, chainage,
or data quality, use the functions documented in :doc:`selection`.

Bulk Edits
----------

:meth:`Sites.edit_all` applies common operations across the collection.

.. code-block:: python
   :linenos:

   renamed = sites.edit_all(
       rename=lambda name: f"LINE01_{name}",
       inplace=False,
   )

   sliced = sites.edit_all(
       freq_slice=slice(1, None),
       inplace=False,
   )

The optional ``mask`` callback receives the output of
``site.to_dataframe("z")`` and returns a boolean mask. Rows where the mask is
``False`` are set to ``NaN`` in the impedance tensor.

.. code-block:: python
   :linenos:

   def keep_lower_frequencies(frame):
       cutoff = frame.index.to_series().median()
       return frame.index <= cutoff

   masked = sites.edit_all(mask=keep_lower_frequencies)

Use :mod:`pycsamt.site.edit` when you need more explicit editing functions
such as tensor rotation, coordinate table assignment, frequency range
selection, or resistivity/phase recomputation.

Topography Alignment
--------------------

:meth:`Sites.with_topography` updates site coordinates and elevation from a
table-like object. Station identifiers are matched by normalized names.

.. code-block:: python
   :linenos:

   import pandas as pd

   topo = pd.DataFrame(
       {
           "station": ["S01", "S02"],
           "latitude": [10.125, 10.250],
           "longitude": [20.500, 20.625],
           "elevation": [640.0, 655.0],
       }
   )

   updated = sites.with_topography(topo, inplace=False)
   print(updated["S01"].coords)

Unmatched stations are left unchanged. Use this after loading raw EDI files
when field GPS or corrected elevation tables are authoritative.

Closest Site
------------

:meth:`Sites.closest` finds the nearest station to a target coordinate.

.. code-block:: python
   :linenos:

   nearest = sites.closest(lat=10.20, lon=20.55)

   if nearest is not None:
       print(nearest.name)

Add ``tol`` in meters when a maximum distance is required:

.. code-block:: python
   :linenos:

   nearest = sites.closest(lat=10.20, lon=20.55, tol=500.0)

If the closest station is farther than ``tol``, the method returns ``None``.

Profiles
--------

:meth:`Sites.to_profile` converts station coordinates into a survey-line
profile for line-ordered work.

.. code-block:: python
   :linenos:

   profile = sites.to_profile(origin=(10.0, 20.0), azimuth=90.0)

   if hasattr(profile, "chainages"):
       print(profile.chainages)
   else:
       print([site.name for site in profile["sites"]])

When :class:`pycsamt.site.profile.Profile` is available, the method returns a
rich profile object. Otherwise it returns a fallback dictionary containing
``origin``, ``azimuth``, and line-ordered ``sites``.

For more control over profile inference, chainage, and spacing statistics,
continue to :doc:`location_profile`.

Writing Sites
-------------

:meth:`Sites.write` writes one EDI file per station into a directory.

.. code-block:: python
   :linenos:

   paths = sites.write(
       "outputs/clean_edi",
       template="{station}.edi",
       exist_ok=True,
   )

   for path in paths:
       print(path)

For richer export workflows, including manifest CSV files and zip packages,
use the functions documented in :doc:`export_reporting`.

Common Patterns
---------------

Inspect one station:

.. code-block:: python
   :linenos:

   site = sites["S01"]
   print(site.summary())
   print(site.to_dataframe("resphase").head())

Build a QC table:

.. code-block:: python
   :linenos:

   import pandas as pd

   qc = pd.DataFrame([site.summary() for site in sites])
   print(qc[["name", "nfreq", "components", "tipper"]])

Select, align, and export:

.. code-block:: python
   :linenos:

   selected = sites.select(predicate=lambda site: site.has_component("Zxy"))
   updated = selected.with_topography(topo)
   updated.write("outputs/selected_edi", exist_ok=True)

Common Mistakes
---------------

Assuming every object has ``Sites.from_path``
   Use :class:`pycsamt.seg.collection.EDICollection.from_sources`,
   :meth:`Sites.from_any`, or :func:`to_sites` depending on the input form.

Editing ``site.edi`` directly during routine workflows
   Direct EDI edits are sometimes necessary, but the site helpers keep name,
   coordinate, and array conventions more consistent.

Forgetting copy versus in-place behaviour
   :meth:`Site.rename`, :meth:`Site.set_coords`, :meth:`Site.set_empty`,
   :meth:`Sites.edit_all`, and :meth:`Sites.with_topography` return new
   objects by default unless ``inplace=True`` is used.

Comparing station names before normalization
   Wrap EDI files as :class:`Site` or :class:`Sites` before doing name-based
   joins or lookups.

Using ``summary`` as a substitute for data inspection
   ``summary`` is a quick overview. Use :meth:`Site.to_dataframe`,
   :meth:`Site.quality_flags`, and the diagnostics in
   :doc:`computed_diagnostics` for deeper checks.

Next Pages
----------

Continue with:

* :doc:`selection` for station filtering;
* :doc:`editing` for explicit single-site and collection edits;
* :doc:`location_profile` for coordinates, topography, distances, and
  chainage;
* :doc:`computed_diagnostics` for strike, resistivity, phase, and tipper
  diagnostics.
