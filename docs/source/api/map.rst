pycsamt.map
===========

Code-first station, profile, and 3-D mapping tools for pyCSAMT.

The :mod:`pycsamt.map` package is the reusable mapping layer for
scripts, notebooks, and applications.  It accepts the same survey
inputs as :func:`pycsamt.emtools._core.ensure_sites`: an EDI
directory, one or more EDI files, EDI-like objects, or a
:class:`pycsamt.site.Sites` container.  For multi-line surveys, use
:func:`pycsamt.map.load_lines` or :class:`pycsamt.map.MapView`.

Quick Use
---------

Function-style use is compact and works well in notebooks.

.. code-block:: python

   from pycsamt.map import (
       ProfileMapOptions,
       StationMapOptions,
       VolumeMapOptions,
       plot_pseudosection,
       plot_station_map,
       plot_volume_map,
   )

   line = "data/AMT/WILLY_DATA/L18PLT"

   station_fig = plot_station_map(
       line,
       options=StationMapOptions(
           overlay="rho",
           frequency=10.0,
           component="xy",
           frequency_tolerance=2.0,
       ),
   )

   pseudo_fig = plot_pseudosection(
       line,
       options=ProfileMapOptions(
           quantity="phase",
           components=("xy", "yx"),
           x_axis="distance",
       ),
   )

   volume_fig = plot_volume_map(
       line,
       options=VolumeMapOptions(
           mode="fence",
           component="xy",
           depth_range=(0.0, 2000.0),
       ),
   )

Builder API
-----------

Builder classes keep one normalized survey in memory and let you
create variants without reloading the EDI data.

.. code-block:: python

   from pycsamt.map import ProfileMap, StationMap, VolumeMap

   station = (
       StationMap("data/AMT/WILLY_DATA/L18PLT")
       .with_overlay("skin_depth", frequency=10.0)
       .with_options(show_labels=False)
       .figure()
   )

   profile = (
       ProfileMap("data/AMT/WILLY_DATA/L18PLT")
       .with_quantity("rho")
       .with_components("xy", "yx")
       .pseudosection()
   )

   volume = (
       VolumeMap("data/AMT/WILLY_DATA/L18PLT")
       .with_mode("surface")
       .with_quantity("phase")
       .figure()
   )

MapView
-------

:class:`pycsamt.map.MapView` is a session facade for code users who
want the same survey to drive station, profile, and 3-D views.

.. code-block:: python

   from pycsamt.map import MapView

   mv = MapView.from_folder(
       "data/AMT/WILLY_DATA",
       detect="folder",
   )

   station = mv.station(overlay="elevation")
   pseudo = mv.pseudosection(component="xy")
   fence = mv.map3d(mode="fence", depth_range=(0.0, 2000.0))

Export
------

Use :func:`pycsamt.map.export_figure` for format-aware export, or
the direct helpers for explicit output types.

.. code-block:: python

   from pycsamt.map import ExportOptions, export_figure

   export_figure(
       station,
       ExportOptions(path="station.html"),
   )

   export_figure(
       fence,
       ExportOptions(path="fence", format="html"),
   )

Scientific Note
---------------

The 3-D map modes in :mod:`pycsamt.map` are quick-look
visualizations.  Fence, block, depth-slice, and isosurface views use
pseudo-depth grids derived from apparent resistivity and period.  They
are useful for inspection and communication, but they are not a
replacement for a constrained 2-D or 3-D inversion model.

Public API
----------

.. automodule:: pycsamt.map
   :members:
   :show-inheritance:

Modules
-------

.. autosummary::
   :toctree: generated

   pycsamt.map.config
   pycsamt.map._core
   pycsamt.map.loader
   pycsamt.map.view
   pycsamt.map.station
   pycsamt.map.profile
   pycsamt.map.volume
   pycsamt.map.overlays
   pycsamt.map.styles
   pycsamt.map.export
   pycsamt.map.inversion
   pycsamt.map.topo

Detailed Module API
-------------------

.. automodule:: pycsamt.map._core
   :members:

.. automodule:: pycsamt.map.loader
   :members:

.. automodule:: pycsamt.map.view
   :members:

.. automodule:: pycsamt.map.station
   :members:

.. automodule:: pycsamt.map.profile
   :members:

.. automodule:: pycsamt.map.volume
   :members:

.. automodule:: pycsamt.map.config
   :members:

.. automodule:: pycsamt.map.overlays
   :members:

.. automodule:: pycsamt.map.export
   :members:
