.. _applications_web_overview:

Overview
========

When To Use The Web App
---------------------------

Use the web application when you need to:

* run the full CSAMT/AMT/MT workflow from a browser, without a local Python
  environment on every workstation;
* share a running instance on an internal server so several people can review
  the same survey;
* demonstrate loading, QC, correction, modelling, and interpretation to a
  team or class;
* combine interactive Plotly maps and sections with matplotlib publication
  figures on one screen;
* drive the pyCSAMT agents (LLM, processing, and inversion) from a chat or a
  runnable agent list.

For scripted or batch work, use the Python API. For local, single-user
interactive review with native windows, use the :doc:`desktop GUI
<../desktop/index>`. For a dedicated conversational surface, use
:doc:`Agent Master <../agent_master/index>`.

Run The Web App
-------------------

.. code-block:: bash

   pycsamt-web

The launcher starts a local Dash server and opens the app in your default
browser. It prefers ``http://127.0.0.1:8050`` and automatically chooses
another free port if that one is busy.

From a source checkout, the module entry point launches the same app:

.. code-block:: bash

   python -m pycsamt.app.web

Common launch options:

.. code-block:: bash

   pycsamt-web --no-browser        # start the server without opening a browser
   pycsamt-web --port 8060         # prefer a specific port
   pycsamt-web --port 0            # always ask the OS for a free port
   pycsamt-web --host 0.0.0.0      # allow access from other machines on the LAN
   pycsamt-web --debug             # enable Dash dev tools and callback debugging

See :doc:`installation` for the full option list and :doc:`deployment` for
host, port, and network guidance.

Navigation Map
------------------

The left rail groups every page under a small number of headings. The order
mirrors a conservative survey workflow: look at the survey, check the data,
analyse it, process it, model it, then review results.

.. list-table::
   :header-rows: 1
   :widths: 18 22 60

   * - Group
     - Page
     - Purpose
   * - (top)
     - Home
     - Survey dashboard, KPIs, station list, and quick workflows.
   * - Survey
     - Map View
     - Interactive station basemap, colouring, contours, and coordinates.
   * - Survey
     - Profile View
     - Per-station responses, pseudosections, tipper, and phase tensor.
   * - Data
     - Quality Control
     - Coverage, noise/SNR, skew/dimensionality, static-shift, distortion QC.
   * - Data
     - Correction
     - The 25-method non-destructive correction chain with preview and undo.
   * - Analysis
     - Advanced Plots
     - Strike, phase tensor, induction, impedance, and depth diagnostics.
   * - Analysis
     - TDEM
     - Time-domain decay, section, map, and dashboard plots.
   * - Processing
     - Pipeline
     - The eight-step load → QC → edit → correct → strike → export workflow.
   * - Modelling
     - Forward Model
     - 1-D, 2-D, and 3-D forward responses.
   * - Modelling
     - Inversion
     - Traditional, AI neural, PINN, and hybrid inversion in 1-D/2-D/3-D.
   * - Results
     - Interpretation
     - Geological, hydrological, and diagnostic interpretation workflows.
   * - Results
     - 3D Map
     - Fence sections, volume, and depth slices in an interactive 3-D scene.
   * - Results
     - Results View
     - Browse ModEM, Occam2D, and MARE2DEM inversion result folders.
   * - Results
     - AI Agents
     - Runnable agents and a conversational assistant over the survey.

Main Concepts
----------------

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Concept
     - Meaning In The Web App
   * - Active survey
     - The loaded survey shared by every page — station list, maps, QC,
       corrections, modelling, and agents all read the same data.
   * - Survey lines
     - Named subgroups of stations (usually one folder per line). Lines can
       be toggled **Active** so a page acts only on the lines you select.
   * - Command bar
     - The persistent top bar: page indicator, **Tools**, **Settings**,
       **Help**, and the survey actions **Load Data**, **+Lines**,
       **Recompute**, **Lines**, **Session**, and **Theme**.
   * - Navigation rail
     - The collapsible left sidebar that routes between pages; it can be
       collapsed to icons to widen the plotting area.
   * - Correction chain
     - A non-destructive stack of correction steps, previewed before
       **Apply**, with per-step **Undo** and full **Reset All**.
   * - Browser session
     - Auto-saved state kept in the browser (theme, API key, last view) plus
       a downloadable session JSON for sharing across machines.
   * - Exported product
     - A figure, corrected EDI, inversion product, or session file meant to be
       reused outside the current browser state.

Screenshots
--------------

Screenshots are embedded inside the pages where they explain a real action or
decision, not collected in a separate gallery. Launch and welcome screens are
in :doc:`installation`; the command bar, rail, and drawers are in
:doc:`navigation`; loading and line management are in
:doc:`loading_and_sessions`; maps and profiles are in :doc:`maps_and_profiles`;
and the processing, modelling, results, and agent surfaces are in
:doc:`processing_pages`.
