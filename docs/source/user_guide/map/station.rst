Station Maps
============

Station maps show station positions, line traces, labels, and scalar
overlays.  Use them for survey QC, coverage inspection, and quick
comparison of per-station quantities.

The station-map API is intentionally code-first.  It accepts the same
sources described in :doc:`loading`: EDI folders, EDI files, file
lists, EDI-like objects, ``Sites`` containers, or an already-normalized
:class:`pycsamt.map.MapData`.

What A Station Map Shows
------------------------

A station map combines three kinds of information:

* station locations or station order;
* optional profile-line traces and station labels;
* one scalar value per station, called the ``overlay``.

If at least one finite coordinate pair is available, the Plotly backend
builds a geographic map.  If no coordinates are available, it falls
back to a profile-style station-index plot rather than failing.

Function API
------------

Use :func:`pycsamt.map.plot_station_map` for one-shot scripts and
notebooks.  Rendering options are passed with
:class:`pycsamt.map.StationMapOptions`.

.. code-block:: python

   from pycsamt.map import StationMapOptions, plot_station_map

   fig = plot_station_map(
       "data/AMT/WILLY_DATA/L18PLT",
       options=StationMapOptions(
           overlay="rho",
           frequency=10.0,
           frequency_tolerance=2.0,
           component="xy",
           show_profiles=True,
           show_labels=True,
       ),
   )

The function returns a Plotly figure by default.  You can display it in
a notebook with ``fig.show()`` or export it with the helpers described
in :doc:`export`.

Reusing Loaded Data
-------------------

For repeated maps, load once and reuse the same
:class:`pycsamt.map.MapData`.

.. code-block:: python

   from pycsamt.map import (
       StationMapOptions,
       ensure_map_data,
       plot_station_map,
   )

   data = ensure_map_data("data/AMT/WILLY_DATA/L18PLT")

   rho = plot_station_map(
       data,
       options=StationMapOptions(
           overlay="rho",
           frequency=100.0,
           component="xy",
       ),
   )

   phase = plot_station_map(
       data,
       options=StationMapOptions(
           overlay="phase",
           frequency=100.0,
           component="xy",
           cmap="Cividis",
       ),
   )

Builder API
-----------

Use :class:`pycsamt.map.StationMap` when you want an immutable builder
style.  ``with_overlay`` and ``with_options`` return new builders that
reuse the same normalized data.

.. code-block:: python

   from pycsamt.map import StationMap

   fig = (
       StationMap("data/AMT/WILLY_DATA/L18PLT")
       .with_overlay("skin_depth", frequency=10.0)
       .with_options(show_labels=False)
       .figure()
   )

This pattern is convenient when a notebook or application lets a user
switch overlays without re-reading the EDI files.

Overlays
--------

Supported overlay names include:

``index`` or ``station``
   Station index along the line.

``elevation``
   Station elevation, when present.

``rho`` or ``resistivity``
   Apparent resistivity at the nearest frequency.

``phase``
   Phase at the nearest frequency.

``skin_depth`` or ``depth``
   Skin-depth scale computed from apparent resistivity and the actual
   selected frequency.

custom column name
   If the overlay name matches a column in the station table, that
   column is used.  Otherwise the station index is used as a safe
   fallback and the colorbar is labelled with the requested overlay
   name.

Frequency-Based Overlays
------------------------

``rho``, ``phase``, and ``skin_depth`` need a frequency.  The map uses
the nearest finite positive frequency for each station.  Use
``frequency_tolerance`` when the match must be strict.

.. code-block:: python

   from pycsamt.map import StationMapOptions, plot_station_map

   fig = plot_station_map(
       "data/AMT/WILLY_DATA/L18PLT",
       options=StationMapOptions(
           overlay="rho",
           frequency=100.0,
           frequency_tolerance=5.0,
           component="xy",
           log_color=True,
           value_range=(10.0, 10000.0),
       ),
   )

The ``component`` option accepts the same component names used by the
map core helpers:

``"xy"``, ``"yx"``, ``"xx"``, ``"yy"``
   Individual impedance tensor components.

``"avg"``
   Average of ``xy`` and ``yx``.

``"det"``
   Determinant-style derived value for resistivity, or an average for
   phase.

When ``log_color=True``, only positive values are transformed with
``log10``.  Non-positive values become gaps in the color scale.

Labels, Lines, And Selection
----------------------------

Station maps can label stations, draw line traces, filter lines, and
highlight one station.

.. code-block:: python

   from pycsamt.map import StationMapOptions, load_lines, plot_station_map

   data = load_lines("data/AMT/WILLY_DATA", detect="folder")

   fig = plot_station_map(
       data,
       options=StationMapOptions(
           overlay="elevation",
           line_filter=("L18PLT", "L22PLT"),
           selected_id="S012",
           show_profiles=True,
           show_labels=False,
           marker_size=9,
       ),
   )

``line_filter`` compares against normalized line names in
``data.lines``.  ``selected_id`` increases the selected station marker
size; it does not remove other stations.

Basemaps
--------

Geographic station maps use public Plotly basemap styles by default.
The theme controls the default:

``theme="light"``
   Uses ``"open-street-map"`` unless ``basemap`` is set.

``theme="dark"``
   Uses ``"carto-darkmatter"`` unless ``basemap`` is set.

Pass ``basemap`` to choose a specific style.

.. code-block:: python

   fig = plot_station_map(
       data,
       options=StationMapOptions(
           overlay="index",
           basemap="carto-positron",
           bearing=15.0,
       ),
   )

The map package also supports token-free ESRI raster basemaps through
the shared basemap helpers.  Common style names include
``"esri-satellite"``, ``"esri-topo"``, and ``"esri-natgeo"``.

Density And Contour Layers
--------------------------

For quick spatial trends, enable a density layer:

.. code-block:: python

   fig = plot_station_map(
       data,
       options=StationMapOptions(
           overlay="rho",
           frequency=100.0,
           show_contours=True,
           contour_opacity=0.45,
       ),
   )

``show_contours`` adds a density-style Plotly map layer when at least
three finite coordinate/value triples are available.

For a filled image layer similar to a Surfer-style contour map, enable
``contour_image``:

.. code-block:: python

   fig = plot_station_map(
       data,
       options=StationMapOptions(
           overlay="rho",
           frequency=100.0,
           contour_image=True,
           contour_mode="filled+lines",
           contour_levels=16,
           contour_opacity=0.55,
       ),
   )

``contour_image`` rasterizes contours to a transparent PNG and inserts
it below station markers as a map image layer.  It returns no layer
when there are fewer than three finite points or when log scaling would
leave fewer than three positive values.

Themes And Color Scales
-----------------------

Station maps share the map theme and color utilities.

.. code-block:: python

   fig = plot_station_map(
       data,
       options=StationMapOptions(
           overlay="phase",
           frequency=100.0,
           theme="publication",
           cmap="Turbo",
           opacity=0.85,
           title="L18 phase at 100 Hz",
       ),
   )

Use ``value_range=(min, max)`` for stable comparisons across multiple
maps.  This is especially important when exporting a sequence of maps
for a report.

Backends
--------

``plotly`` is the interactive default.  ``matplotlib`` produces static
figures for reports and batch processing.

.. code-block:: python

   from pycsamt.map import StationMapOptions, plot_station_map

   fig = plot_station_map(
       "data/AMT/WILLY_DATA/L18PLT",
       options=StationMapOptions(backend="matplotlib"),
   )

Backend behavior differs slightly:

``plotly``
   Produces interactive geographic maps when coordinates exist, and
   station-index plots when they do not.

``matplotlib``
   Produces a static Matplotlib figure.  It uses longitude/latitude
   when available and station index/value otherwise.

Unknown backends raise ``ValueError``.

Coordinate Fallback
-------------------

Station maps do not require coordinates to produce a useful diagnostic
view.  If at least one finite latitude/longitude pair exists, the
Plotly backend attempts a geographic map.  If no geographic coordinates
exist, it creates a Cartesian station-index plot.

.. code-block:: python

   from pycsamt.map import ensure_map_data

   data = ensure_map_data("data/AMT/WILLY_DATA/L18PLT")

   if not data.has_geo:
       print("Using profile fallback because coordinates are incomplete.")

For production geographic maps, validate ``data.has_geo`` before
plotting when every station must appear at a real location, and fix
coordinates through your survey metadata or CRS preprocessing.

Exporting Station Maps
----------------------

Use the export helpers to save figures consistently:

.. code-block:: python

   from pycsamt.map import write_html, save_png

   write_html(fig, "outputs/stations.html")
   save_png(fig, "outputs/stations.png", width=1400, height=900)

PNG export requires a static image backend such as Kaleido for Plotly
figures.  HTML export works without Kaleido.

Troubleshooting
---------------

Empty figure
   No stations were loaded, or filtering removed every line.  Inspect
   ``data.station_ids`` and ``data.lines``.

Map appears as a station-index plot
   No finite geographic coordinates were available.  Check
   ``data.has_geo`` and your station metadata.

``rho`` or ``phase`` colors are missing
   The requested frequency may not exist within
   ``frequency_tolerance``, or the station may not have a valid ``Z``
   object.  Try removing the tolerance to inspect nearest-frequency
   behavior.

Contours do not appear
   The density and contour-image layers need at least three finite
   coordinate/value triples.  They are skipped silently when there is
   not enough data.

Colorbars differ between maps
   Set the same ``value_range`` and ``log_color`` options on each map.
