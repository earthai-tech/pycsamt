.. _applications_desktop_overview:

Overview
========

Who This Guide Is For
-------------------------

This guide is written for users who need to make decisions from survey data,
not only click through the interface. The desktop app is useful to several
roles, and each role usually cares about a different question:

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - User
     - Main Question
     - Desktop Pages To Read First
   * - Field geophysicist
     - Did the line load correctly, and which stations need attention?
     - :doc:`loading_and_sessions`, :doc:`maps_and_profiles`,
       :doc:`processing_workflows`
   * - Processing specialist
     - Which corrections are justified by the evidence?
     - :doc:`workspace`, :doc:`processing_workflows`, :doc:`exports`
   * - Inversion modeller
     - Is this data set ready to become solver input?
     - :doc:`maps_and_profiles`, :doc:`processing_workflows`,
       :doc:`exports`
   * - Interpreter or reviewer
     - Can I trace every figure and conclusion back to a reproducible state?
     - :doc:`processing_workflows`, :doc:`exports`,
       :doc:`troubleshooting`
   * - New pyCSAMT user
     - What should I do first, and what should I avoid changing too early?
     - :doc:`installation`, :doc:`loading_and_sessions`,
       :doc:`workspace`

The pages below therefore explain both **how** to use each section and **why**
that section exists in the scientific workflow. When a page describes a
button such as **Commit to Main**, **Run All**, **Export**, or inversion
**Run**, read it as a data-state decision, not just an interface action.

Guided First Session
------------------------

If this is your first time using the desktop app, do not start with every
tool. Start with one small survey line and work through the minimum path that
teaches the interface and protects the data state:

.. list-table::
   :header-rows: 1
   :widths: 12 24 36 28

   * - Step
     - Do This
     - Why It Matters
     - Read Next
   * - 1
     - Launch the desktop and load a known EDI folder.
     - Confirms the environment, file loader, and active survey state before
       any interpretation begins.
     - :doc:`installation`, :doc:`loading_and_sessions`
   * - 2
     - Check the station table, overview card, and station detail card.
     - Catches wrong folders, missing coordinates, and suspicious frequency
       coverage early.
     - :doc:`workspace`
   * - 3
     - Open the map and profile viewers.
     - Verifies spatial geometry, station order, response curves, and
       pseudosection continuity before processing.
     - :doc:`maps_and_profiles`
   * - 4
     - Run first-pass QC diagnostics.
     - Separates data-health problems from features that might be geology.
     - :doc:`processing_workflows`
   * - 5
     - Preview one justified correction, if needed.
     - Shows the effect before the active survey is changed.
     - :doc:`processing_workflows`
   * - 6
     - Export figures, corrected data, or notes.
     - Turns the interactive session into a reproducible record.
     - :doc:`exports`

This first pass is successful when you can explain what data were loaded,
which stations looked reliable, which figures support that judgement, and
whether the active survey is still raw or has been corrected. It is better to
finish one small line carefully than to open every panel and lose track of the
state.

When To Use The Desktop
---------------------------

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
application. For guided workflow reasoning, use Agent Master.

What A Complete Desktop Session Produces
--------------------------------------------

A useful desktop session should leave more behind than an open window. At the
end of a review, the user should be able to point to:

* the input files or folder that defined the active survey;
* the station count, coordinate sanity checks, and frequency coverage reviewed
  before processing;
* the QC evidence used to keep, reject, or mark stations;
* the correction previews used before any correction was committed;
* the exported corrected EDIs or derived products, when the data changed;
* the figures, maps, profiles, or tables used in a report;
* the pipeline JSON or notes needed to repeat the same work later;
* the inversion or modelling folder, if the session prepared solver input.

If one of those items is missing, the session may still be useful for
exploration, but it is not yet a reproducible processing record. The export
page explains how to turn an interactive review into files that another user,
script, or future you can inspect.

Run The Desktop GUI
-----------------------

.. code-block:: bash

   pycsamt-desktop

The shorter historical alias still launches the same desktop app:

.. code-block:: bash

   pycsamt-gui

From a source checkout, the module entry point is also available:

.. code-block:: bash

   python -m pycsamt.app.desktop

Workflow Map
---------------

The desktop workflow is intentionally conservative:

.. code-block:: text

   load -> inspect -> QC -> preview correction -> commit -> recheck -> model
        -> prepare inversion -> export -> document

The map/profile viewers should be used before and after any state-changing
operation. **Apply**, **Commit to Main**, **Run All**, and inversion **Run**
are not routine navigation buttons; they change data state or generate new
products. Use them only when the current diagnostics explain what should
happen next.

Decision Points In The Workflow
------------------------------------

The desktop app is organized around a few repeated decisions:

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Decision
     - Why It Matters
     - Evidence To Check
   * - Load one line or many
     - A focused line is easier to QC and prepare for inversion; a broader
       folder is useful for reconnaissance.
     - Station count, line labels, coordinate spread, and map geometry.
   * - Trust station metadata
     - Bad coordinates, elevation, or frequency metadata can make later plots
       look persuasive but wrong.
     - Station table, detail card, map view, and profile ordering.
   * - Apply a correction
     - Corrections should solve a visible problem, not simply improve a curve
       aesthetically.
     - Before/after response plots, QC metrics, neighbouring stations, and
       physical plausibility.
   * - Commit corrected data
     - Committing replaces the active survey state used by all downstream
       panels.
     - Correction stack, preview plots, log messages, and export plan.
   * - Prepare inversion input
     - Solver inputs should be built only after station selection, errors,
       modes, and frequency bands are defensible.
     - QC summary, selected stations, profile geometry, errors, and exported
       configuration.
   * - Export results
     - Files are the handoff from an interactive session to a report,
       notebook, solver, or reviewer.
     - Output folder, metadata, figures, corrected EDIs, and pipeline JSON.

Main Concepts
----------------

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
--------------

Screenshots are embedded inside the workflow pages where they explain a real
decision or user action. There is no separate screenshot gallery for the
desktop app. For example:

* loader and launch screenshots are in :doc:`installation` and
  :doc:`loading_and_sessions`;
* map/profile screenshots are in :doc:`maps_and_profiles`;
* QC, correction, forward, inversion, pipeline, and agent screenshots are in
  :doc:`processing_workflows`.
