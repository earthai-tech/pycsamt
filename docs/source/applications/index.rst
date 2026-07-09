.. _applications:

Applications
============

pyCSAMT includes four user-facing applications built on the same
scientific core — the same readers, processing tools, plotting
conventions, and configuration system you use from Python. They differ
only in working style: a native desktop GUI for local interactive
review, a browser app for shared servers and team demonstrations, a
conversational surface that delegates workflows to the pyCSAMT agents,
and a dedicated map workbench for seeing a survey in space. Results are
interchangeable — a survey processed in one surface can be picked up in
any other, or in plain Python.

Pick a card below to open an application's guide; each one starts with
installation and launch and ends with exports and troubleshooting.

Application Surfaces
--------------------

.. grid:: 1 1 2 2
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Desktop app
      :link: desktop/index
      :link-type: doc
      :img-top: ../_static/icons/desktop-app.svg
      :class-card: pycsamt-card sd-text-center

      Local interactive survey review, processing, and plotting in a
      native GUI — sessions, workspaces, and exports on your machine.

   .. grid-item-card:: Web app
      :link: web/index
      :link-type: doc
      :img-top: ../_static/icons/web-app.svg
      :class-card: pycsamt-card sd-text-center

      Browser-based workflows, dashboards, maps, and processing pages —
      ideal for team demonstrations and shared servers.

   .. grid-item-card:: Agent Master
      :link: agent_master/index
      :link-type: doc
      :img-top: ../_static/icons/agent-master.svg
      :class-card: pycsamt-card sd-text-center

      Conversational workflow delegation, guided automation, and report
      help — the chat surface over the pyCSAMT agents.

   .. grid-item-card:: MapView
      :link: mapview/index
      :link-type: doc
      :img-top: ../_static/icons/mapview-app.svg
      :class-card: pycsamt-card sd-text-center

      The dedicated map workbench — station maps, pseudosections, and
      interactive 3-D survey scenes with fence sections and topography.

Install The App Extra
---------------------

The application surfaces use optional GUI and web dependencies:

.. code-block:: bash

   pip install "pycsamt[app]"

For development from a source checkout:

.. code-block:: bash

   pip install -e ".[app,dev]"

.. toctree::
   :maxdepth: 1
   :hidden:

   desktop/index
   web/index
   agent_master/index
   mapview/index

.. toctree::
   :hidden:

   desktop
   web
