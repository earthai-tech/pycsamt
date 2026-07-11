.. _applications-mapview-exports:

Exports
=======

MapView produces two kinds of durable output: an **interactive HTML** file of
the current view for sharing, and — through the Python layer it shares — fully
scripted **static figures** for publication.

.. note::

   **Export** and **Session** are different.  **Export** (below) writes a
   figure of the *current view* to share with others.  **Session** (see
   :doc:`loading_and_sessions`) saves and restores your *working state* inside
   MapView.  Use Export to hand a scene to someone; use Session to resume your
   own work.

Interactive HTML
----------------

The **Export** button on the top bar saves the **current view** as a
standalone interactive HTML file (``mapview_<view>.html``) — the full Plotly
figure with pan, zoom, rotate, and hover intact.  The file loads plotly.js from
a CDN, so it stays small and can be mailed or dropped into a shared drive;
opening it needs only a browser, not a pyCSAMT installation.

This works for every view, including the 3-D scene — a fence-section scene
exported this way remains fully rotatable for the recipient.  For a quick
raster instead, the Plotly modebar's **camera** button on the Map canvas saves
a PNG of the current framing.

Static Figures And Beyond
-------------------------

For publication-quality static exports (PNG/SVG at a chosen size and
camera), drive the same figures from code with the
:doc:`Python mapping layer </user_guide/map/index>`:

.. code-block:: python

   from pycsamt.map import MapView

   view = MapView.from_folder("path/to/edi_folder", recursive=True)
   fig = view.map3d(mode="fence", quantity="rho")
   fig.write_image("survey_fence.png", width=1200, height=800)

Everything MapView draws is reproducible this way — the app and the API
share one figure factory.
