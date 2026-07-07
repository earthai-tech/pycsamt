Exporting Map Figures
=====================

Use :func:`pycsamt.map.export_figure` when you want the export format
to be inferred from the path or supplied through
:class:`pycsamt.map.ExportOptions`.

The export helpers work with the figure objects returned by station,
profile, 3-D, overlay, and :class:`pycsamt.map.MapView` workflows.
They intentionally stay small: they create parent directories, route to
the right backend method, and return the output path.

Choosing An Export Helper
-------------------------

Use the helper that matches your goal:

``write_html``
   Save an interactive Plotly figure as HTML.

``save_png``
   Save a Matplotlib figure or a Plotly figure as PNG.

``write_image``
   Save a Plotly static image in any format supported by Plotly image
   export, such as PNG, SVG, PDF, or WebP.

``write_json``
   Save a Plotly figure JSON specification.

``write_dict``
   Save a serialized figure dictionary as JSON text.

``figure_to_dict``
   Return a Python dictionary without writing a file.

``export_figure``
   Dispatch to the correct helper from an
   :class:`pycsamt.map.ExportOptions` object.

Output Directories
------------------

All write helpers create parent directories automatically.

.. code-block:: python

   from pycsamt.map import write_html

   write_html(fig, "outputs/maps/station.html")

If the directory does not exist, it is created before the figure is
written.

HTML
----

HTML is the safest default for Plotly station, profile, and 3-D maps.
It preserves hover text, legends, zooming, panning, map tiles, and 3-D
rotation.

.. code-block:: python

   from pycsamt.map import ExportOptions, export_figure

   export_figure(
       fig,
       ExportOptions(path="station.html"),
   )

You can also call the direct helper:

.. code-block:: python

   from pycsamt.map import write_html

   write_html(
       fig,
       "station.html",
       include_plotlyjs="cdn",
   )

``include_plotlyjs`` is passed to Plotly's ``write_html`` method.
Common values are:

``"cdn"``
   Keep files smaller by loading Plotly JavaScript from a CDN.

``True``
   Embed Plotly JavaScript inside the HTML file for offline sharing.

``False``
   Do not include Plotly JavaScript.  Use only when the surrounding
   page already loads Plotly.

Explicit Format
---------------

If the path has no suffix, set ``format``.  ``format`` may be provided
with or without a leading dot.

.. code-block:: python

   export_figure(
       fig,
       ExportOptions(path="station", format="html"),
   )

Without a suffix or explicit format, ``export_figure`` raises
``ValueError`` rather than guessing.

PNG
---

Matplotlib figures are written with ``savefig``.  Plotly figures use
``write_image`` and require Kaleido for static image export.

.. code-block:: python

   export_figure(
       fig,
       ExportOptions(path="station.png", scale=2.0),
   )

Direct helper:

.. code-block:: python

   from pycsamt.map import save_png

   save_png(
       fig,
       "station.png",
       width=1400,
       height=900,
       scale=2.0,
   )

For Plotly figures, ``width``, ``height``, and ``scale`` are forwarded
to Plotly image export.  For Matplotlib figures, Plotly-specific keys
are filtered out before calling ``savefig``.

.. code-block:: python

   save_png(
       matplotlib_fig,
       "station_static.png",
       dpi=200,
       bbox_inches="tight",
   )

Static Plotly images require Kaleido.  If Kaleido is missing, pyCSAMT
raises an ``ImportError`` with a clearer message.

Other Static Image Formats
--------------------------

Use :func:`pycsamt.map.write_image` for non-PNG Plotly image formats:

.. code-block:: python

   from pycsamt.map import write_image

   write_image(fig, "station.svg", width=1200, height=800)
   write_image(fig, "profile.pdf", width=1400, height=900)

``export_figure`` also routes unknown suffixes to ``write_image``.
This lets Plotly decide whether a format is supported.

JSON And Dict Export
--------------------

``write_json`` uses Plotly's JSON methods when available and falls
back to dictionary export.  ``write_dict`` writes the figure dictionary
as JSON text.

.. code-block:: python

   from pycsamt.map import write_dict, write_json

   write_json(fig, "station.json")
   write_dict(fig, "station.dict")

Use JSON or dict exports when you want to archive the figure
specification, compare figure output in tests, or post-process traces
in another tool.

.. code-block:: python

   from pycsamt.map import figure_to_dict

   spec = figure_to_dict(fig)
   print(len(spec.get("data", [])))

``figure_to_dict`` requires a Plotly-like object with ``to_dict`` or
``to_plotly_json``.  Matplotlib figures do not expose this contract.

ExportOptions Reference
-----------------------

:class:`pycsamt.map.ExportOptions` stores the common export controls:

``path``
   Output path.

``format``
   Optional format override.  Use this when ``path`` has no suffix or
   when you want to force a specific export route.

``width`` and ``height``
   Static image dimensions for Plotly image export.

``scale``
   Static Plotly image scale factor.  Defaults to ``2.0``.

``include_plotlyjs``
   Passed to Plotly HTML export.  Defaults to ``"cdn"``.

Batch Export
------------

Because the helpers return :class:`pathlib.Path`, batch workflows can
record or log output paths easily.

.. code-block:: python

   from pycsamt.map import (
       ExportOptions,
       StationMapOptions,
       export_figure,
       plot_station_map,
   )

   overlays = ["index", "elevation", "rho", "phase"]
   outputs = []

   for overlay in overlays:
       fig = plot_station_map(
           data,
           options=StationMapOptions(
               overlay=overlay,
               frequency=100.0,
           ),
       )
       outputs.append(
           export_figure(
               fig,
               ExportOptions(
                   path=f"outputs/station_{overlay}.html",
               ),
           )
       )

   print(outputs)

MapView Export
--------------

:class:`pycsamt.map.MapView` uses the same export machinery.  It builds
the requested view, creates ``ExportOptions``, and calls
``export_figure`` internally.

.. code-block:: python

   from pycsamt.map import MapView

   mv = MapView(data)
   mv.export(
       "outputs/fence.html",
       view="volume",
       mode="fence",
       depth_range=(0.0, 2500.0),
   )

This keeps application-style workflows and direct function workflows on
the same export path.

Round-Trip And Testing Workflows
--------------------------------

JSON and dictionary exports are useful when you want to inspect or test
figure content without rendering images.

.. code-block:: python

   from pycsamt.map import figure_to_dict

   spec = figure_to_dict(fig)
   trace_types = [trace.get("type") for trace in spec["data"]]
   assert "scatter" in trace_types or "scattermap" in trace_types

For persistent artifacts, prefer ``write_json``:

.. code-block:: python

   from pycsamt.map import write_json

   write_json(fig, "outputs/station_spec.json")

Recommended Formats
-------------------

Station maps
   Use HTML for inspection and PNG/SVG for reports.

Profile pseudosections
   Use HTML while tuning options, then PNG/SVG/PDF for publication
   drafts.

3-D maps
   Prefer HTML.  Static images capture only one camera angle.

Automated tests
   Use ``figure_to_dict`` or ``write_json`` rather than pixel images.

Troubleshooting
---------------

``TypeError: HTML export requires a Plotly-like figure``
   ``write_html`` only works with figures exposing ``write_html``.
   Matplotlib figures should be saved with ``save_png`` or
   ``fig.savefig``.

``ImportError`` mentioning Kaleido
   Plotly static image export needs Kaleido.  Install the full/docs
   environment or export HTML instead.

``ValueError`` about extension or format
   ``export_figure`` needs a file suffix such as ``.html`` or
   ``.png``, or an explicit ``ExportOptions.format``.

Unexpected Matplotlib keyword behavior
   ``save_png`` removes Plotly-only keys such as ``scale``, ``width``,
   and ``height`` before calling Matplotlib ``savefig``.  Use
   Matplotlib options like ``dpi`` and ``bbox_inches`` for static
   Matplotlib figures.

Dictionary export fails
   The object must expose ``to_dict`` or ``to_plotly_json``.  Use JSON
   or dict export for Plotly-like figures, not plain Matplotlib figures.
