MapView Session
===============

:class:`pycsamt.map.MapView` keeps one normalized survey in memory and
renders multiple map views from it.  It is the most convenient API for
multi-line surveys and repeated exports.

Create A Session
----------------

.. code-block:: python

   from pycsamt.map import MapView

   mv = MapView.from_folder(
       "data/AMT/WILLY_DATA",
       detect="folder",
   )

   print(mv.lines)
   print(mv.stations)
   print(mv.n_stations)

Render Views
------------

.. code-block:: python

   station = mv.station(overlay="elevation")
   pseudo = mv.pseudosection(component="xy")
   volume = mv.map3d(mode="fence", depth_range=(0.0, 2000.0))

Generic Figure Dispatch
-----------------------

.. code-block:: python

   fig = mv.figure("station", overlay="rho", frequency=10.0)

Supported view names are ``station``, ``profile``, ``pseudosection``,
and ``map3d``.

Export
------

.. code-block:: python

   mv.export(
       "station.html",
       view="station",
       overlay="rho",
       frequency=10.0,
   )

   mv.export_all("out/maps", fmt="html")

Launch The Platform
-------------------

``MapView.launch`` and :func:`pycsamt.map.open_app` hand the session to
the optional map-view platform.  The GUI dependencies are optional and
are not required for code-first rendering.
