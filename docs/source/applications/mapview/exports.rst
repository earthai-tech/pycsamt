Exports
=======

Interactive HTML
----------------

The export button saves the **current view** as a standalone
interactive HTML file (``mapview_<view>.html``) — the full Plotly
figure with pan, zoom, rotate, and hover intact. The file loads
plotly.js from a CDN, so it stays small and can be mailed or dropped
into a shared drive; opening it needs a browser but no pyCSAMT
installation.

This works for every view, including the 3-D scene — a fence-section
scene exported this way remains fully rotatable for the recipient.

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
