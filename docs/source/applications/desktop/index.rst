.. _applications-desktop:

Desktop GUI
===========

The pyCSAMT desktop GUI is the local application for interactive CSAMT, AMT,
and MT survey work.  It is designed for the day-to-day workflow of loading EDI
data, checking station geometry, inspecting profiles and maps, running QC,
testing corrections, preparing forward and inversion models, exporting figures
and processed EDIs, and documenting the processing chain.

The desktop app uses the same scientific package as the Python API.  The GUI
is not a separate science layer; it is an interactive surface over pyCSAMT
loaders, ``Sites`` objects, processing functions, plotting tools, inversion
builders, and export helpers.

When To Use The Desktop
-----------------------

Use the desktop GUI when you need to:

* inspect a survey line station by station;
* compare maps, profiles, QC plots, and correction previews side by side;
* decide whether static shift, frequency editing, source-effect handling, or
  tensor rotation is justified;
* build visual evidence for reports and reproducible processing notes;
* prepare forward models and inversion input folders;
* export corrected EDIs, figures, pipeline settings, and interpretation files.

For automated scripts, batch jobs, or notebooks, use the Python API directly.
For browser-based demonstrations or shared dashboards, use the web
application.  For guided workflow reasoning, use Agent Master.

Run The Desktop GUI
-------------------

.. code-block:: bash

   pycsamt-desktop

The shorter historical alias still launches the same desktop app:

.. code-block:: bash

   pycsamt-gui

From a source checkout, the module entry point is also available:

.. code-block:: bash

   python -m pycsamt.app.desktop

Start Here
----------

Read the desktop documentation in this order if you are setting up the app for
the first time:

1. :doc:`installation` -- install dependencies, launch the desktop, and confirm
   the console commands.
2. :doc:`loading_and_sessions` -- load EDI/AVG/J files, understand the active
   ``Sites`` object, and save session state.
3. :doc:`workspace` -- learn the main window, station table, toolbar, floating
   panels, log dock, status bar, theme, and preferences.
4. :doc:`maps_and_profiles` -- inspect station geometry, contours, response
   curves, pseudosections, topography, and phase tensors.
5. :doc:`processing_workflows` -- run evidence-based QC, corrections, advanced
   diagnostics, forward modelling, inversion preparation, pipeline processing,
   and agent-assisted review.
6. :doc:`exports` -- save figures, metadata, EDIs, pipeline JSON, inversion
   outputs, interpretation products, and reproducibility notes.
7. :doc:`troubleshooting` -- diagnose launch, loading, plotting, export,
   correction, pipeline, and inversion problems.

Workflow Map
------------

The desktop workflow is intentionally conservative:

.. code-block:: text

   load -> inspect -> QC -> preview correction -> commit -> recheck -> model
        -> prepare inversion -> export -> document

The map/profile viewers should be used before and after any state-changing
operation.  **Apply**, **Commit to Main**, **Run All**, and inversion **Run**
are not routine navigation buttons; they change data state or generate new
products.  Use them only when the current diagnostics explain what should
happen next.

Main Concepts
-------------

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Concept
     - Meaning In The Desktop
   * - Active survey
     - The currently loaded ``Sites`` object shared by the station table,
       map/profile viewers, QC, corrections, modelling, inversion, and agents.
   * - Panel windows
     - Independent floating windows such as Profile, Map, QC, Corrections,
       Forward, Inversion, Interpretation, Pipeline, Advanced Tools, and
       Agent Master.
   * - Correction stack
     - Non-destructive correction steps previewed and applied locally before
       **Commit to Main** replaces the active survey.
   * - Pipeline configuration
     - JSON record of processing methods and parameters for repeatable
       load/QC/edit/correct/rotate/export workflows.
   * - Session state
     - Convenience state stored under ``~/.pycsamt/session.json``: theme,
       recent files, selected station, geometry, solver paths, and preferences.
   * - Exported product
     - A file or folder meant to be reused outside the current GUI state:
       figures, corrected EDIs, metadata tables, solver inputs, inversion
       results, interpretation exports, or reports.

Screenshots
-----------

Screenshots are embedded inside the workflow pages where they explain a real
decision or user action.  There is no separate screenshot gallery for the
desktop app.  For example:

* loader and launch screenshots are in :doc:`installation` and
  :doc:`loading_and_sessions`;
* map/profile screenshots are in :doc:`maps_and_profiles`;
* QC, correction, forward, inversion, pipeline, and agent screenshots are in
  :doc:`processing_workflows`.

Documentation Pages
-------------------

.. toctree::
   :maxdepth: 1

   installation
   loading_and_sessions
   workspace
   maps_and_profiles
   processing_workflows
   exports
   troubleshooting
