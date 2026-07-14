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

Typical workflow
----------------

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

Guide Sections
--------------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Loading Surveys
      :link: loading
      :link-type: doc
      :img-top: ../../_static/icons/loading.svg
      :class-card: pycsamt-card sd-text-center

      Build a ``MapView`` from folders, globs, or site collections —
      single lines or whole multi-line surveys.

   .. grid-item-card:: Station Maps
      :link: station
      :link-type: doc
      :img-top: ../../_static/icons/location.svg
      :class-card: pycsamt-card sd-text-center

      Plan-view station maps with resistivity, phase, or elevation
      overlays and basemap styles.

   .. grid-item-card:: Profiles & Pseudosections
      :link: profile
      :link-type: doc
      :img-top: ../../_static/icons/map-and-profile.svg
      :class-card: pycsamt-card sd-text-center

      Distance-frequency pseudosections and profile views along any
      survey line.

   .. grid-item-card:: 3-D Quick-Look
      :link: volume
      :link-type: doc
      :img-top: ../../_static/icons/model.svg
      :class-card: pycsamt-card sd-text-center

      Pseudo-depth volumes, fence diagrams, and topography-draped 3-D
      survey scenes.

   .. grid-item-card:: Overlays
      :link: overlays
      :link-type: doc
      :img-top: ../../_static/icons/controls.svg
      :class-card: pycsamt-card sd-text-center

      Contours, annotations, auxiliary layers, and styling controls on
      top of any map view.

   .. grid-item-card:: Export
      :link: export
      :link-type: doc
      :img-top: ../../_static/icons/export.svg
      :class-card: pycsamt-card sd-text-center

      Save interactive HTML or static PNG figures for reports and
      publications.

   .. grid-item-card:: MapView Workbench
      :link: mapview
      :link-type: doc
      :img-top: ../../_static/icons/mapview-app.svg
      :class-card: pycsamt-card sd-text-center

      The interactive Dash workbench built on the same rendering
      layer.

   .. grid-item-card:: API Reference
      :link: api
      :link-type: doc
      :img-top: ../../_static/icons/api-reference-icon-braces.svg
      :class-card: pycsamt-card sd-text-center

      Classes, options, and helpers of :mod:`pycsamt.map` in one
      reference page.

.. toctree::
   :maxdepth: 1
   :hidden:

   loading
   station
   profile
   volume
   overlays
   export
   mapview
   api
