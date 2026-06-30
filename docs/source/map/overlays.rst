Overlays
========

The overlay helpers are reusable building blocks for custom maps and
application views.  They are useful when you want to add CRS
transforms, contours, labels, profile lines, or topography to figures
that you control yourself.

Use the high-level map builders first when you can:
:class:`pycsamt.map.StationMap`, :class:`pycsamt.map.ProfileMap`, and
:class:`pycsamt.map.VolumeMap` already add the common survey overlays.
Use this page when you need to assemble traces manually, add your own
computed values, or keep the same overlay logic in scripts, notebooks,
and application callbacks.

Overlay Workflow
----------------

Most overlay workflows follow the same pattern:

1. Load and normalize survey data.
2. Extract coordinates and values.
3. Reproject coordinates if needed.
4. Build one or more Plotly traces or layout settings.
5. Add the traces to a figure you own.

.. code-block:: python

   import plotly.graph_objects as go

   from pycsamt.map import (
       build_contour_overlay,
       build_profile_line_overlay,
       build_station_label_overlay,
       ensure_map_data,
   )

   data = ensure_map_data("data/AMT/WILLY_DATA/L18PLT")
   lon = [station.longitude for station in data.stations]
   lat = [station.latitude for station in data.stations]
   labels = data.station_ids
   values = [station.index for station in data.stations]

   fig = go.Figure()
   fig.add_trace(build_contour_overlay(lon, lat, values))
   fig.add_trace(build_profile_line_overlay(lon, lat, geo=False))
   fig.add_trace(build_station_label_overlay(lon, lat, labels))

Coordinate Reference Systems
----------------------------

Many CSAMT surveys store station positions in projected coordinates,
while web maps usually expect WGS84 longitude and latitude.  The CRS
helpers keep that conversion explicit.

``CRSConfig``
   Describes a source CRS, target CRS, and coordinate-axis convention.

``normalize_epsg``
   Converts integers or strings such as ``4326`` into ``"EPSG:4326"``.

``transform_xy``
   Transforms coordinate arrays between arbitrary CRS definitions.

``reproject_xy_to_lonlat``
   Convenience wrapper for projected coordinates to WGS84 lon/lat.

CRS Helpers
-----------

.. code-block:: python

   import numpy as np

   from pycsamt.map import CRSConfig, transform_xy

   lon, lat = transform_xy(
       np.array([500000.0, 500250.0]),
       np.array([850000.0, 850200.0]),
       crs=CRSConfig(source=32630, target=4326),
   )

The returned arrays are NumPy arrays.  With the default
``always_xy=True``, the input order is always ``x, y`` and the output
for EPSG:4326 is ``lon, lat``.  This avoids the common latitude/longitude
axis-order surprise in modern PROJ/GDAL stacks.

For display text or logs, use :func:`pycsamt.map.resolve_crs_info`.

.. code-block:: python

   from pycsamt.map import resolve_crs_info

   print(resolve_crs_info("geo"))
   print(resolve_crs_info("utm", zone=30, hemisphere="N"))
   print(resolve_crs_info("custom", epsg=32630))

Basemaps
--------

``build_basemap_layout`` creates the geographic map settings shared by
station maps and application views.  It does not return a full Plotly
layout dictionary; it returns a :class:`pycsamt.map.BasemapConfig`
dataclass so callers can apply it to either new Plotly map APIs or
legacy Mapbox-backed figures.

.. code-block:: python

   from pycsamt.map import build_basemap_layout

   basemap = build_basemap_layout(
       lon,
       lat,
       dark=False,
       bearing=0.0,
   )

   print(basemap.style)
   print(basemap.center)
   print(basemap.zoom)

Default styles are public and token-free:

``dark=False``
   Uses ``"open-street-map"``.

``dark=True``
   Uses ``"carto-darkmatter"``.

Pass ``style=...`` when your application has its own tile policy:

.. code-block:: python

   basemap = build_basemap_layout(
       lon,
       lat,
       style="carto-positron",
   )

If all coordinates are missing or non-finite, the helper returns a
world-scale center at ``{"lat": 0.0, "lon": 0.0}`` and zoom ``1``.

Contours
--------

Contours interpolate scattered values to a regular grid.  SciPy is
used when available; otherwise the helper falls back to a nearest
neighbour grid.  This makes contour overlays useful in lightweight
environments where SciPy is not installed, while still giving smoother
linear interpolation when SciPy is available.

.. code-block:: python

   import numpy as np

   from pycsamt.map import build_contour_overlay

   contour = build_contour_overlay(
       np.array([2.0, 2.1, 2.2]),
       np.array([1.0, 1.1, 1.0]),
       np.array([100.0, 120.0, 80.0]),
       levels=8,
       mode="both",
   )

The coordinate arrays and value array must contain at least three
finite points.  Non-finite points are ignored before interpolation.

Important options:

``levels``
   Number of contour bands.  Values below two are promoted to two.

``cmap``
   Plotly colorscale name, such as ``"Viridis"``, ``"Cividis"``, or
   ``"Turbo"``.

``opacity``
   Trace opacity from ``0`` to ``1``.

``grid_size``
   Number of grid nodes along each axis.  Larger grids look smoother
   but cost more memory and rendering time.

``mode``
   One of ``"lines"``, ``"fill"``, ``"filled"``, ``"heatmap"``, or
   ``"both"``.  ``"lines"`` draws contour lines only; ``"both"`` uses
   heatmap coloring with contour lines.

When you need the interpolated grid rather than a Plotly trace, call
:func:`pycsamt.map.interpolate_overlay_grid`.

.. code-block:: python

   from pycsamt.map import interpolate_overlay_grid

   xi, yi, grid = interpolate_overlay_grid(
       lon,
       lat,
       values,
       grid_size=120,
   )

Labels And Profile Lines
------------------------

Station labels and profile lines are simple Plotly traces.  They can be
added on top of a station map, a custom contour figure, or a geographic
scatter map.

.. code-block:: python

   from pycsamt.map import (
       build_profile_line_overlay,
       build_station_label_overlay,
   )

   line = build_profile_line_overlay(lon, lat, geo=True)
   labels = build_station_label_overlay(
       lon,
       lat,
       ["S00", "S01"],
       geo=True,
   )

Use ``geo=True`` when ``x`` and ``y`` are longitude and latitude for a
Plotly map trace.  Use ``geo=False`` for ordinary Cartesian figures,
including profile plots and local projected-coordinate maps.

.. code-block:: python

   labels = build_station_label_overlay(
       lon,
       lat,
       labels=data.station_ids,
       geo=True,
       color="#0f172a",
   )

   profile = build_profile_line_overlay(
       lon,
       lat,
       geo=True,
       name="L18",
       color="#2563eb",
       width=2.5,
   )

The helpers use Plotly's modern ``Scattermap`` trace when available and
fall back to ``Scattermapbox`` for older Plotly versions.

Topography
----------

``build_topography_overlay`` returns ``Mesh3d`` for scattered
elevations and ``Surface`` for 2-D elevation grids.

Scattered station elevations:

.. code-block:: python

   from pycsamt.map import build_topography_overlay

   topo = build_topography_overlay(
       x,
       y,
       elevation,
       opacity=0.55,
       colorscale="Earth",
   )

Regular elevation grid:

.. code-block:: python

   topo = build_topography_overlay(
       grid_x,
       grid_y,
       elevation_grid,
       opacity=0.40,
   )

Use topography overlays in 3-D figures where the vertical dimension is
meaningful.  For station maps, use elevation as a station overlay or a
contour value instead of adding a 3-D terrain trace.

Adding Overlays To Figures
--------------------------

The helpers return traces, so they can be added directly to a Plotly
figure:

.. code-block:: python

   import plotly.graph_objects as go

   fig = go.Figure()
   fig.add_trace(contour)
   fig.add_trace(line)
   fig.add_trace(labels)

For geographic station maps, apply the basemap settings to the figure
layout:

.. code-block:: python

   fig.update_layout(
       map=dict(
           style=basemap.style,
           center=basemap.center,
           zoom=basemap.zoom,
           bearing=basemap.bearing,
       )
   )

For older Plotly versions that use Mapbox layout keys, applications may
need to map the same dataclass to ``mapbox`` instead:

.. code-block:: python

   fig.update_layout(
       mapbox=dict(
           style=basemap.style,
           center=basemap.center,
           zoom=basemap.zoom,
           bearing=basemap.bearing,
       )
   )

Practical Examples
------------------

Contour station resistivity at one frequency:

.. code-block:: python

   from pycsamt.map import (
       build_contour_overlay,
       ensure_map_data,
       value_at_frequency_details,
   )

   data = ensure_map_data("data/AMT/WILLY_DATA/L18PLT")
   values_by_station = value_at_frequency_details(
       data,
       frequency=100.0,
       quantity="rho",
       component="xy",
   )

   lon = []
   lat = []
   values = []
   for station in data.stations:
       if station.id not in values_by_station:
           continue
       lon.append(station.longitude)
       lat.append(station.latitude)
       values.append(values_by_station[station.id].value)

   contour = build_contour_overlay(
       lon,
       lat,
       values,
       levels=12,
       cmap="Viridis",
       mode="both",
   )

Build labels for every profile in a multi-line survey:

.. code-block:: python

   from pycsamt.map import (
       build_profile_line_overlay,
       build_station_label_overlay,
       load_lines,
   )

   data = load_lines("data/AMT/WILLY_DATA", detect="folder")
   traces = []

   for profile in data.profiles:
       lon = [station.longitude for station in profile.stations]
       lat = [station.latitude for station in profile.stations]
       labels = [station.id for station in profile.stations]

       traces.append(
           build_profile_line_overlay(
               lon,
               lat,
               geo=True,
               name=profile.name,
           )
       )
       traces.append(
           build_station_label_overlay(
               lon,
               lat,
               labels,
               geo=True,
           )
       )

Troubleshooting
---------------

``ImportError`` from CRS transforms
   ``transform_xy`` requires ``pyproj``.  Install the ``geo`` or
   ``full`` extra before using CRS conversion helpers.

``ValueError: At least three finite points are required.``
   Contour interpolation needs at least three finite coordinate/value
   triples.  Drop incomplete stations or use a scatter trace instead.

Empty or world-scale basemap
   The input longitude/latitude arrays have no finite coordinate pairs.
   Inspect ``data.has_geo`` before building geographic overlays.

Contour is blocky
   SciPy may be unavailable, so interpolation falls back to nearest
   neighbour.  Install SciPy or increase ``grid_size`` if appropriate.

Labels do not appear on a map
   Check the ``geo`` flag.  Geographic Plotly maps need ``geo=True``;
   Cartesian figures need ``geo=False``.
