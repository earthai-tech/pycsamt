Map Tools
=========

The :mod:`pycsamt.map` package is the code-first mapping layer for
pyCSAMT surveys.  It gives script and notebook users direct control
over station maps, profile pseudosections, 3-D quick-look maps,
overlays, multi-line survey loading, and figure export.

Use this section when you want to work with maps from Python code.  The
web and desktop applications can build on the same rendering layer, but
UI state and callbacks do not live here.

.. warning::

   3-D maps in :mod:`pycsamt.map` are pseudo-depth visualizations
   derived from apparent resistivity and period.  They are useful for
   QC, screening, and communication, but they are not constrained
   inversion models.

.. toctree::
   :numbered: 4
   :maxdepth: 1
   :class: pycsamt-guide-toc

   loading
   station
   profile
   volume
   overlays
   export
   mapview
   api

Typical workflow
------------------

.. code-block:: python

   from pycsamt.map import MapView

   mv = MapView.from_folder(
       "data/AMT/WILLY_DATA",
       detect="folder",
   )

   station = mv.station(overlay="rho", frequency=10.0)
   section = mv.pseudosection(component="xy")
   fence = mv.map3d(mode="fence", depth_range=(0.0, 2000.0))

   mv.export("station.html", view="station", overlay="rho")
