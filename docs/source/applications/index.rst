.. _applications:

Applications
============

pyCSAMT includes three user-facing applications. They share the same
scientific package, loaders, processing tools, and plotting conventions, but
serve different working styles.

Application Surfaces
--------------------

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Application
     - Best for
     - Documentation area
   * - Desktop GUI
     - Local interactive survey review, processing, and plotting.
     - :doc:`desktop/index`
   * - Web application
     - Browser-based workflows, dashboards, maps, and team demonstrations.
     - :doc:`web/index`
   * - Agent Master
     - Conversational workflow delegation, guided automation, and report help.
     - :doc:`agent_master/index`

Install The App Extra
---------------------

The application surfaces use optional GUI and web dependencies:

.. code-block:: bash

   pip install "pycsamt[app]"

For development from a source checkout:

.. code-block:: bash

   pip install -e ".[app,dev]"

Page Plan
---------

Each application folder starts with the same documentation rhythm:

* installation and launch;
* data loading and session behavior;
* workspace or navigation guide;
* core workflows;
* export workflows and reproducibility notes;
* screenshots gallery placeholders;
* troubleshooting.

The screenshots will be added later under
``docs/source/_static/applications/``.

.. toctree::
   :maxdepth: 2

   desktop/index
   web/index
   agent_master/index

.. toctree::
   :hidden:

   desktop
   web
