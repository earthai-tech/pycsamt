.. _applications_mapview_overview:

Overview
========

When To Use MapView
-----------------------

MapView is the right surface when the **spatial picture** is the point:

* place a multi-line survey on a satellite, street, topographic, or light/dark
  basemap and confirm its geometry;
* colour stations by index, elevation, apparent resistivity, phase, or skin
  depth, and drape a Surfer-style contour over the basemap;
* build an interactive 3-D scene — fence panels along each line, a filled
  block, depth slices, or iso-surfaces — with topography and a resistivity
  visibility band to isolate conductive or resistive structure;
* drop a ModEM 3-D **inversion result** into the same scene as the survey
  geometry and inspect the model in place.

For per-station response curves, pseudosections, QC, corrections, and
modelling, use the :doc:`web app <../web/index>` or the :doc:`desktop GUI
<../desktop/index>`. MapView deliberately does one thing — space — and does it
well.

The Two Views
----------------

Everything in MapView lives under **two views**, switched from the left rail:

* **Map** — the 2-D basemap view: stations, profile lines, contour overlay,
  basemap, and coordinate system.
* **3-D** — the interactive resistivity scene: fence, block, depth-slice, and
  iso-surface renderings with topography and depth/resistivity filters.

Both views read the same loaded survey and the same line selection, so
switching between them never changes *what* is selected — only how it is drawn.

Run MapView
--------------

.. code-block:: bash

   pycsamt-mapview

The launcher starts a local Dash server on ``http://127.0.0.1:8770`` and opens
the app in the default browser. The module form works everywhere the package
is importable:

.. code-block:: bash

   python -m pycsamt.app.mapview

Useful launch options:

.. code-block:: bash

   pycsamt-mapview --no-browser
   pycsamt-mapview --port 8780
   pycsamt-mapview --host 0.0.0.0
   pycsamt-mapview --data path/to/edi_folder
   pycsamt-mapview --debug

``--data`` preloads an EDI folder (one or more survey lines) so the app opens
with the survey already on the map. See :doc:`installation` for the full
option list.

Main Concepts
----------------

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Concept
     - Meaning In MapView
   * - Survey lines
     - Named groups of stations (one folder or filename prefix per line),
       toggled from the line picker; the selection applies to both views.
   * - View
     - **Map** (2-D) or **3-D** — the two ways MapView draws the same survey.
   * - Colour field / overlay
     - The quantity drawn on the map: station index, elevation, apparent
       resistivity, phase, or skin depth.
   * - Contour overlay
     - A Surfer-style filled/line contour draped over the basemap.
   * - 3-D mode
     - Fence, Block, Depth slices, or Iso-surface — how the 3-D scene renders
       the resistivity volume.
   * - Resistivity band
     - A ρ visibility window that isolates the conductive or resistive part of
       the 3-D volume.
   * - Inversion result
     - A ModEM 3-D model folder imported and shown in the survey's own scene.
   * - Session cache
     - Loaded data held per browser tab on the server; the app never writes
       into your data folders.

.. seealso::

   :doc:`/user_guide/map/index`
       The Python mapping layer behind the app — build the same figures
       from code.

   :ref:`3-D maps gallery <sphx_glr_examples_map_3d>`
       Executed examples of every station-map and 3-D view, with code.
