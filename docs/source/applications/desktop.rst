.. _desktop_application:

Desktop Application
===================

The pyCSAMT desktop application is a PySide6 interface for interactive
survey work. It is intended for users who need to inspect stations,
visualize profiles, run QC, launch processing tools, prepare inversion
inputs, review results, and coordinate agent-assisted workflows without
writing code for every operation.

.. note::

   This page is a stub. It records the current entry points and major
   application areas so the documentation navigation is complete. A full
   user guide will be written later.

Run The Desktop App
-------------------

Install the application extra:

.. code-block:: bash
   :linenos:

   pip install "pycsamt[app]"

Start the application from an installed environment:

.. code-block:: bash
   :linenos:

   pycsamt-gui

From a source checkout, the module entry point is also available:

.. code-block:: bash
   :linenos:

   python -m pycsamt.app.desktop

Current Architecture
--------------------

The desktop application lives under :mod:`pycsamt.app.desktop`.

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Package area
     - Role
   * - ``main_window.py``
     - Main application shell and top-level window composition.
   * - ``controllers/``
     - Application logic for data, plotting, QC, correction, forward
       modelling, interpretation, pipeline, settings, and TDEM workflows.
   * - ``panels/``
     - Reusable UI panels for stations, maps, profiles, sections, QC,
       agents, and logs.
   * - ``windows/``
     - Dedicated workflow windows for advanced tools, corrections, forward
       modelling, interpretation, inversion, maps, profiles, pipelines, QC,
       TDEM, and agents.
   * - ``tools/``
     - Utility tools such as converters, coordinate transforms, validators,
       frequency editing, elevation, dimensionality, phase tensor maps, and
       strike/profile analysis.
   * - ``workers/``
     - Background execution wrappers for loading, agents, forward modelling,
       inversion, AI inversion, and pipelines.
   * - ``models/``
     - Session and station data models.
   * - ``resources/``
     - Themes, icons, and application assets.

Feature Areas To Document
-------------------------

The detailed desktop manual should later cover:

* opening EDI/AVG/J/TDEM survey data;
* survey session management and persistence;
* station table filtering and station-detail inspection;
* map viewer, profile viewer, and section viewer behavior;
* QC summaries, coverage plots, signal-to-noise views, and diagnostics;
* correction tools including static shift and source-effect workflows;
* forward modelling windows and synthetic response inspection;
* inversion preparation, execution workers, and result viewers;
* interpretation and export tools;
* pipeline window, preset selection, and run monitoring;
* AI agent panel, agent registry, and background execution logs;
* desktop preferences, themes, and plotting settings.

Screenshot Assets
-----------------

Desktop screenshot assets are already available under
``docs/source/_static/app/desktop/``. They should be used when this page is
expanded into a full visual guide.

.. image:: ../_static/app/desktop/home.png
   :alt: pyCSAMT desktop application home screen
   :width: 100%

Related API
-----------

The application API stub is available at :doc:`../api/app`.
