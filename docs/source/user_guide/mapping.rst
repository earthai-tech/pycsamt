Code-First Maps
===============

The :mod:`pycsamt.map` package provides the open-source mapping API for
users who want direct control from Python code instead of using only
the web application.  The package can render station maps,
pseudosections, 3-D quick-look views, overlays, and exported figures
from scripts, notebooks, or applications.

Loading One Or More Lines
-------------------------

For one line, pass the EDI directory directly to a plotting function.
For a survey made of several line folders, use
:func:`pycsamt.map.load_lines` or :class:`pycsamt.map.MapView`.

.. code-block:: python

   from pycsamt.map import MapView, load_lines

   data = load_lines(
       "data/AMT/WILLY_DATA",
       detect="folder",
   )

   mv = MapView(data)
   print(mv.lines)
   print(mv.stations)

Station Maps
------------

Use :func:`pycsamt.map.plot_station_map` to inspect station positions,
survey lines, labels, and scalar overlays.  Overlays include station
index, elevation, apparent resistivity, phase, and skin depth.

.. code-block:: python

   from pycsamt.map import StationMapOptions, plot_station_map

   fig = plot_station_map(
       "data/AMT/WILLY_DATA/L18PLT",
       options=StationMapOptions(
           overlay="rho",
           frequency=10.0,
           frequency_tolerance=2.0,
           component="xy",
           basemap="open-street-map",
           show_labels=True,
           show_profiles=True,
       ),
   )
   fig.write_html("station_map.html")

Pseudosections
--------------

Use :func:`pycsamt.map.plot_pseudosection` for interactive
apparent-resistivity or phase sections.  Set ``x_axis="distance"``
when station coordinates are available and you want profile distance
rather than station names on the horizontal axis.

.. code-block:: python

   from pycsamt.map import ProfileMapOptions, plot_pseudosection

   fig = plot_pseudosection(
       "data/AMT/WILLY_DATA/L18PLT",
       options=ProfileMapOptions(
           quantity="rho",
           components=("xy", "yx"),
           period_range=(0.001, 10.0),
           x_axis="distance",
       ),
   )
   fig.show()

3-D Quick-Look Maps
-------------------

Use :func:`pycsamt.map.plot_volume_map` for fence, block, depth-slice,
and isosurface views.

.. code-block:: python

   from pycsamt.map import VolumeMapOptions, plot_volume_map

   fig = plot_volume_map(
       "data/AMT/WILLY_DATA/L18PLT",
       options=VolumeMapOptions(
           mode="fence",
           quantity="rho",
           component="xy",
           depth_range=(0.0, 2000.0),
           period_range=(0.001, 10.0),
           opacity=0.85,
       ),
   )
   fig.write_html("survey_3d.html")

The available 3-D modes are:

``fence``
    One pseudo-depth surface per survey line.
``block``
    A sparse volume rendering from all finite samples.
``depth``
    Horizontal pseudo-depth slices.
``surface``
    Isosurfaces across the pseudo-depth point cloud.

These are pseudo-depth products derived from apparent resistivity and
period.  They are useful for QC, target screening, and communication,
but they are not a replacement for inversion.

Builder API
-----------

Builder objects keep normalized survey data in memory while allowing
small option changes.

.. code-block:: python

   from pycsamt.map import ProfileMap, StationMap, VolumeMap

   station_fig = (
       StationMap("data/AMT/WILLY_DATA/L18PLT")
       .with_overlay("skin_depth", frequency=10.0)
       .with_options(show_labels=False)
       .figure()
   )

   pseudo_fig = (
       ProfileMap("data/AMT/WILLY_DATA/L18PLT")
       .with_quantity("phase")
       .with_component("xy")
       .pseudosection()
   )

   volume_fig = (
       VolumeMap("data/AMT/WILLY_DATA/L18PLT")
       .with_mode("surface")
       .with_quantity("phase")
       .figure()
   )

Overlays
--------

The overlay helpers are reusable in custom figures.

.. code-block:: python

   import numpy as np

   from pycsamt.map import (
       build_contour_overlay,
       build_profile_line_overlay,
       build_station_label_overlay,
       reproject_xy_to_lonlat,
   )

   lon, lat = reproject_xy_to_lonlat(
       np.array([500000.0, 500250.0, 500100.0]),
       np.array([850000.0, 850200.0, 850350.0]),
       epsg=32630,
   )

   line = build_profile_line_overlay(lon, lat, geo=True)
   labels = build_station_label_overlay(
       lon,
       lat,
       ["S00", "S01", "S02"],
       geo=True,
   )

   contour = build_contour_overlay(
       lon,
       lat,
       np.array([100.0, 120.0, 80.0]),
       grid_size=20,
   )

Export
------

Use :func:`pycsamt.map.export_figure` when the output format should be
inferred from the path or supplied explicitly.

.. code-block:: python

   from pycsamt.map import ExportOptions, export_figure

   export_figure(
       station_fig,
       ExportOptions(path="station.html"),
   )

   export_figure(
       volume_fig,
       ExportOptions(path="volume", format="html"),
   )

Static Plotly image export requires Kaleido.  If Kaleido is not
installed, export to HTML or install the map/image export dependencies.

Relationship to the Web App
---------------------------

The web application keeps UI controls, callbacks, uploads, and session
state in :mod:`pycsamt.app.web`.  Reusable rendering belongs in
:mod:`pycsamt.map`, so users can generate the same scientific map
objects from scripts and notebooks.
