.. _applications-mapview:

MapView
=======

MapView is the dedicated map workbench of pyCSAMT: a Dash/Plotly browser
application focused entirely on **seeing a survey in space** — basemap
station maps, profile overlays, pseudosections, and interactive 3-D
survey scenes with fence sections and draped topography. It is the
point-and-click layer over the :doc:`Python mapping API
</user_guide/map/index>`: everything the app draws comes from the same
:class:`~pycsamt.map.MapView` façade you can drive from code.

Run MapView
-----------

.. code-block:: bash

   pycsamt-mapview

The launcher starts a local Dash server on ``http://127.0.0.1:8770`` and
opens the app in the default browser. The module form works everywhere
the package is importable:

.. code-block:: bash

   python -m pycsamt.app.mapview

Useful launch options:

.. code-block:: bash

   pycsamt-mapview --no-browser
   pycsamt-mapview --port 8780
   pycsamt-mapview --host 0.0.0.0
   pycsamt-mapview --data path/to/edi_folder
   pycsamt-mapview --debug

``--data`` preloads an EDI folder (one or more survey lines) so the app
opens with the survey already on the map.

Documentation Pages
-------------------

.. toctree::
   :maxdepth: 1

   installation
   loading_and_sessions
   views
   exports
   troubleshooting

.. seealso::

   :doc:`/user_guide/map/index`
       The Python mapping layer behind the app — build the same figures
       from code.

   :ref:`3-D maps gallery <sphx_glr_examples_map_3d>`
       Executed examples of every station-map and 3-D view, with code.
