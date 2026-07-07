Loading Map Data
================

The mapping package uses one normalized object,
:class:`pycsamt.map.MapData`, for every station map, profile
pseudosection, and 3-D view.  You rarely need to build
``MapData`` by hand.  Most public plotting functions accept the
same inputs directly and call the normalizer for you.

Use this page when you want explicit control over loading, line
grouping, or pre-flight checks before drawing a map.

Loading Philosophy
------------------

The loader follows the science API rather than the web app internals.
It delegates EDI parsing to
:func:`pycsamt.emtools._core.ensure_sites`, then extracts the small
mapping contract needed by the renderers:

* station identifier;
* latitude, longitude, and elevation when available;
* profile or survey-line name;
* the original EDI-like object, including its impedance ``Z`` object.

This means a code workflow and the web workflow can start from the
same survey files, but code users can inspect, group, transform, and
export each intermediate object.

Input Types
-----------

All loading helpers accept flexible EDI sources:

``str`` or ``pathlib.Path``
   A single EDI file or a directory containing EDI files.

sequence of paths
   A list or tuple of EDI files.  A sequence of directories is treated
   as one line per directory by :func:`pycsamt.map.load_lines`.

EDI-like iterable
   Objects with station metadata and a ``Z`` attribute, including
   objects already loaded elsewhere in pyCSAMT.

``Sites``-like container
   Any object accepted by ``ensure_sites``.  Containers exposing
   ``as_list()`` are unpacked into EDI-like objects.

``MapData``
   Returned unchanged, so it is safe to pass normalized data through
   higher-level builders.

Single-Line Loading
-------------------

For one profile line, call :func:`pycsamt.map.ensure_map_data`
directly or pass the source to a plotting helper.

.. code-block:: python

   from pycsamt.map import ensure_map_data, plot_station_map

   source = "data/AMT/WILLY_DATA/L18PLT"

   data = ensure_map_data(source, recursive=True)
   print(data.station_ids)
   print(data.lines)
   print(data.has_geo)

   fig = plot_station_map(data, frequency=100.0)

``ensure_map_data`` does not force every station to have geographic
coordinates.  ``data.has_geo`` tells you whether all normalized
stations have finite latitude and longitude.  If coordinates are
missing, profile views can still use station order or distance, while
geographic station maps may need CRS or coordinate preprocessing.

Multiple-Line Loading
---------------------

Use :func:`pycsamt.map.load_lines` when a survey is split across
several folders or when you want one combined object for 3-D fence,
block, or depth-slice views.

.. code-block:: python

   from pycsamt.map import load_lines

   data = load_lines(
       "data/AMT/WILLY_DATA",
       detect="folder",
       recursive=True,
   )

   print(data.lines)
   print(data.metadata["n_lines"])
   print(data.station_ids[:5])

Each station is re-indexed across the combined survey and tagged with
its resolved line name.  The renderers then use ``data.profiles`` to
draw separate profiles or offset lines in 3-D.

Grouping Modes
--------------

When ``load_lines`` receives one directory, the ``detect`` option
decides how files become lines.

``detect="folder"``
   Group by immediate parent folder.  This is the best option for
   layouts such as ``WILLY_DATA/L18PLT/*.edi`` and
   ``WILLY_DATA/L22PLT/*.edi``.  A flat folder becomes one line named
   after the folder.

``detect="flat"``
   Treat every discovered EDI file under the directory as one line.
   Use this for a single profile stored in a nested directory tree.

``detect="auto"``
   Group by station-ID prefix.  Numeric prefixes are normalized with an
   ``L`` prefix.  For example, station files beginning with ``22`` are
   grouped under ``L22``.  Use this when line folders are not reliable
   but station naming is consistent.

You can inspect the grouping without loading impedance data by calling
:func:`pycsamt.map.resolve_line_groups`.

.. code-block:: python

   from pycsamt.map import resolve_line_groups

   groups = resolve_line_groups(
       "data/AMT/WILLY_DATA",
       detect="folder",
   )

   for line, source in groups.items():
       print(line, len(source))

Explicit Line Mapping
---------------------

The most reproducible multi-line workflow is an explicit mapping:

.. code-block:: python

   from pycsamt.map import load_lines

   data = load_lines({
       "L18": "data/AMT/WILLY_DATA/L18PLT",
       "L22": "data/AMT/WILLY_DATA/L22PLT",
       "L26": "data/AMT/WILLY_DATA/L26PLT",
   })

Mapping values may be directories, file lists, or EDI-like iterables.
This is useful in notebooks and scripts where line names should be
stable even if folder names later change.

If your EDI-like objects do not carry line metadata, but you already
know station membership, pass ``line_map`` to
:func:`pycsamt.map.ensure_map_data`.

.. code-block:: python

   from pycsamt.map import ensure_map_data

   data = ensure_map_data(
       edis,
       line_map={
           "L18": ["S001", "S002", "S003"],
           "L22": ["S101", "S102", "S103"],
       },
   )

Line metadata embedded on the EDI object takes priority.  The fallback
order is:

* ``line``;
* ``profile``;
* ``survey_line``;
* ``Line``;
* ``Profile``;
* station lookup in ``line_map``;
* default line name ``"line"``.

Normalized Contract
-------------------

``MapData`` exposes the surface that all map renderers share:

.. code-block:: python

   data.station_ids      # tuple[str, ...]
   data.lines            # tuple[str, ...]
   data.has_geo          # bool
   data.iter_edis()      # tuple of original EDI-like objects
   data.stations         # tuple[StationRecord, ...]
   data.profiles         # tuple[ProfileLine, ...]
   data.metadata         # counts and loader metadata

Each :class:`pycsamt.map.StationRecord` contains:

.. code-block:: python

   station.id
   station.latitude
   station.longitude
   station.elevation
   station.line
   station.index
   station.source

``station.source`` is the original EDI-like object.  Keep it when you
need to inspect impedance arrays, frequencies, or file-level metadata
after loading.

Pre-Flight Checks
-----------------

A good loading workflow checks the normalized survey before plotting:

.. code-block:: python

   from pycsamt.map import (
       ensure_map_data,
       select_frequency,
       value_at_frequency_details,
   )

   data = ensure_map_data("data/AMT/WILLY_DATA/L18PLT")
   first_edi = data.iter_edis()[0]
   freq = first_edi.Z.freq
   selection = select_frequency(freq, requested=100.0)

   if selection is None:
       raise RuntimeError("No finite positive frequency was found.")

   print("requested:", selection.requested)
   print("actual:", selection.actual)
   print("relative delta:", selection.relative_delta)

   values = value_at_frequency_details(
       data,
       frequency=selection.actual,
       quantity="rho",
       component="xy",
   )
   print("stations with values:", len(values))

This is especially helpful when different stations have slightly
different frequency grids.  Map extraction chooses the nearest finite
positive frequency per station and can enforce a tolerance when you
need a strict match.

Handling Missing or Partial Data
--------------------------------

The loader is designed for field surveys, where not every file is
perfect.

* Files without valid impedance are skipped by value-extraction helpers
  instead of breaking the whole map.
* Non-finite coordinates become ``None`` on ``StationRecord``.
* If no line name is available, stations are grouped into ``"line"``.
* Empty groups are ignored by ``load_lines``.
* A ``MapData`` with no stations is still a valid object; renderers can
  return empty figures rather than crashing.

For strict workflows, validate the counts yourself:

.. code-block:: python

   data = load_lines("data/AMT/WILLY_DATA", detect="folder")

   if not data.stations:
       raise RuntimeError("No stations were loaded.")
   if not data.has_geo:
       raise RuntimeError("Station map requires geographic coordinates.")
   if len(data.lines) < 2:
       raise RuntimeError("Expected a multi-line survey.")

Using Loaded Data
-----------------

Once loaded, pass the same ``MapData`` object to station, profile, and
3-D builders.  This avoids re-reading files and keeps grouping
consistent across views.

.. code-block:: python

   from pycsamt.map import (
       load_lines,
       plot_station_map,
       plot_pseudosection,
       plot_volume_map,
   )

   data = load_lines("data/AMT/WILLY_DATA", detect="folder")

   station_fig = plot_station_map(data, frequency=100.0)
   section_fig = plot_pseudosection(data, component="xy")
   volume_fig = plot_volume_map(data, mode="fence")

Use the export helpers to save any resulting figure:

.. code-block:: python

   from pycsamt.map import write_html

   write_html(station_fig, "outputs/stations.html")

Troubleshooting
---------------

``FileNotFoundError``
   The path passed to ``load_lines`` or ``ensure_map_data`` does not
   exist.  Resolve project-relative paths before loading.

``ValueError: Unknown detect mode``
   ``detect`` must be one of ``"folder"``, ``"flat"``, or ``"auto"``.

No stations in the output
   Check that EDI files were discovered and that ``ensure_sites`` can
   read them.  Start with ``resolve_line_groups`` to confirm grouping.

Station map is empty
   Inspect ``data.has_geo``.  Geographic maps require finite latitude
   and longitude unless you first reproject or provide coordinates.

Pseudosection has gaps
   Some stations may be missing impedance values, or the requested
   frequency/period may not be present within tolerance.  Inspect
   ``value_at_frequency_details`` for selected-frequency metadata.
