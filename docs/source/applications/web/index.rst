.. _applications-web:

Web Application
===============

The pyCSAMT web application is a Dash/Plotly browser interface for survey
loading, map and profile visualization, processing pages, pipeline workflows,
and browser-friendly exports.

Run The Web App
---------------

.. code-block:: bash

   pycsamt-web

The launcher starts a local Dash server and opens the app in the default
browser. It prefers ``http://127.0.0.1:8050`` and chooses another free port if
needed.

Useful launch options:

.. code-block:: bash

   pycsamt-web --no-browser
   pycsamt-web --port 8060
   pycsamt-web --port 0
   pycsamt-web --host 0.0.0.0
   pycsamt-web --debug

Documentation Pages
-------------------

.. toctree::
   :maxdepth: 1

   installation
   loading_and_sessions
   navigation
   maps_and_profiles
   processing_pages
   deployment
   exports
   troubleshooting
