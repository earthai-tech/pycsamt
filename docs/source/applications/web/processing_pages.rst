.. _applications-web-processing:

Processing Pages
================

The processing pages are where the web app does scientific work: quality
control, corrections, advanced diagnostics, TDEM, the processing pipeline,
forward and inversion modelling, interpretation, the results viewer, and the
agents.  Every page is an interactive surface over the same
``pycsamt.app.desktop.controllers`` used by the desktop app — the callbacks
delegate to the package, so results are identical to the Python API and
interchangeable with the desktop GUI.

Most pages share a common shape: a left-hand control column with a primary
action button (**Generate**, **Run**, **Plot**, **Apply**…), tabs or a group
selector across the plotting area, and per-page **Lines** and **Stations**
pickers so an operation can be scoped to a subset of the survey.

.. _web-processing-shared:

Shared Controls
---------------

Before the individual pages, three patterns recur:

* **Lines / Stations pickers** — most pages let you restrict the operation to
  chosen lines and stations.  These are independent of, but seeded by, the
  global :doc:`Lines drawer <loading_and_sessions>`.
* **The primary action** — pages do not recompute continuously.  A single
  button (**Plot**, **Generate**, **Run …**, **Apply**) triggers the work,
  which keeps large surveys responsive.
* **Result tabs** — modelling and results pages separate the plot, the
  numerical summary, and the run log into tabs so you can inspect each without
  losing the others.

Quality Control
---------------

The **Quality Control** page groups every QC diagnostic behind one control
column.  Use it to decide whether data are healthy before correcting or
modelling.

.. figure:: ../../_static/applications/web/quality-control.png
   :alt: The Quality Control page with group buttons, plot list, and parameters
   :class: pycsamt-screenshot

   Quality Control.  Choose a **Group**, pick a **Plot**, set **Parameters**,
   and press **Plot** — or use **Overview** for a multi-panel quicklook.

* **Filters** — restrict to **Lines**, **Stations**, and a **Component**
  (``XY``, ``YX``, or ``Both``).
* **Groups** — *Overview*, *Coverage*, *Noise / SNR*, *Skew / Dim*,
  *Static Shift*, and *Distortion*.  Each group exposes a set of plots.
* **Plot list** — for example *QC Quicklook (multi-panel)*, *Coverage quality
  heatmap*, *Frequency confidence section*, *Confidence band summary*, and
  *Confidence profile*.
* **Parameters** — plot-specific inputs such as the **metric** (for example
  coherence), the outlier **method** (for example IQR), and a **threshold**.
* **Refresh / Auto** — redraw on demand, or auto-redraw as settings change.

Read QC before deciding on corrections: sparse frequency coverage, erratic
error bars, or isolated low-SNR stations are the kind of findings that justify
frequency editing or static-shift work on the next page.

Correction
----------

The **Correction** page applies pyCSAMT's catalogue of corrections as a
**non-destructive chain**: 25 methods across 6 categories, previewed before
they are committed, with per-step undo.

.. figure:: ../../_static/applications/web/correction.png
   :alt: The Correction page previewing an AMA static-shift correction
   :class: pycsamt-screenshot

   Correction.  Pick a **Category** and **Method**, set **Parameters**,
   **Preview** the before/after, then **Apply** to add the step to the chain.

* **Category** and **Correction method** — choose from the catalogue (for
  example *Static Shift → AMA (spatial average)*).  A short description
  explains what the method does.
* **Parameters** — method-specific inputs (for example sort key, half-window,
  and kernel for a spatial-average static-shift correction).
* **View modes** — *Before / After*, *ρ/φ View*, *Overlay*, and *Difference*
  compare the raw and corrected responses.
* **Show** and **Comp** — display *Raw*, *Corrected*, or *Both*, for the
  ``XY``, ``YX``, or both components, across up to *Max stn* stations.
* **Actions** — **Preview** shows the effect without committing; **Apply** adds
  the step to the chain; **Undo** removes the last step; **Reset All** returns
  to raw data.

The chain is the important idea: corrections are staged and reversible, so you
can preview a static-shift method, apply it, preview a noise method on top, and
still undo back to raw at any point.  Corrected EDIs are exported from this
page; see :doc:`exports`.

Advanced Plots
--------------

The **Advanced Plots** page produces the survey-scale diagnostics used for
strike, dimensionality, and depth reasoning.

.. figure:: ../../_static/applications/web/advanced-plots.png
   :alt: The Advanced Plots page showing a phase-tensor pseudosection
   :class: pycsamt-screenshot

   Advanced Plots on the Phase Tensor tab: a station × period ellipse
   pseudosection coloured by skew β.

* **Tabs** — *Strike*, *Phase Tensor*, *Induction*, *Impedance/Z*, *Depth*,
  and *Survey Tools*.
* **Active Lines** and **Stations** — scope the diagnostic to chosen lines and
  stations.
* **Parameters** — plot-specific inputs and **Figure size** presets.
* **Generate** — build the selected diagnostic; a short **Description**
  explains what the current plot shows.

These are whole-survey versions of the per-line diagnostics in Profile View —
use them for strike consistency and dimensionality across all lines at once.

TDEM
----

The **TDEM** page handles time-domain electromagnetic data with its own folder
browser and a fixed tab bar of plot categories.

* **Tabs** — *Decay / Rho*, *Survey Section*, *Map & Overview*, and
  *Dashboard*.
* Each tab exposes a set of plots, a figure-size preset, and a colour map.
* Point the folder browser at a TDEM dataset, choose a plot, and generate it.

TDEM is grouped under **Analysis** in the rail alongside Advanced Plots.

.. _web-processing-pipeline:

Processing Pipeline
-------------------

The **Pipeline** page runs the standard survey workflow as an ordered,
inspectable sequence.  It is the eight-step ``load → QC → edit → correct →
strike → export`` chain, run step by step or all at once.

The steps are:

.. list-table::
   :header-rows: 1
   :widths: 8 30 62

   * - #
     - Step
     - What it does
   * - 1
     - Load Data
     - Use the loaded survey, or browse to a fresh EDI folder.
   * - 2
     - QC Screening
     - Drop low-confidence frequencies by threshold.
   * - 3
     - Frequency Edit
     - Edit frequencies by confidence (recover or trim).
   * - 4
     - Static Shift Correction
     - Apply a static-shift correction (for example AMA).
   * - 5
     - Noise Removal
     - Denoise the impedance response.
   * - 6
     - Strike Analysis
     - Estimate geoelectric strike.
   * - 7
     - Strike Rotation
     - Rotate to the strike frame.
   * - 8
     - Export
     - Write processed EDIs and products to an export folder.

**Controls**

* **Run Step** runs the current step; **Run All** runs the remaining steps in
  order; **Skip** advances without running; **Reset** returns to the start.
* For each step you choose a **Method** and see a **Description**; the
  **Steps** list and the numbered progress track show where you are.
* The output pane has **Log**, **Preview**, and **Status** tabs.  The log
  reports each step (``— Step 2 done —``) and surfaces failures (for example
  ``ERROR: SVD did not converge`` if a static-shift solve fails) so you can
  adjust parameters and re-run.
* The **Export** step writes to a folder you choose with **Select Export
  Folder**.

The pipeline is the repeatable counterpart to the interactive pages: the same
QC, correction, and strike operations, run as a recorded sequence you can
re-run on new data.

Forward Modelling
-----------------

The **Forward Model** page computes model responses in 1-D, 2-D, and 3-D.

.. figure:: ../../_static/applications/web/forward-model.png
   :alt: The Forward Model page showing 3-D MT model slices for a halfspace
   :class: pycsamt-screenshot

   Forward Model on the 3-D MT tab: model slices for a halfspace, with grid,
   background, and frequency-range controls on the left.

* **Tabs** — *1-D*, *2-D MT*, and *3-D MT (quasi-3D)*.
* **Model Type** — for example *Halfspace* or a layered/preset model.
* **Background & Grid** — background resistivity, station layout, grid
  dimensions (``Nx × Ny × Nz``), and extents.
* **Frequency Range** — ``log10(Hz)`` min/max and the number of points.
* **Run Forward** computes the response; the plotting area shows the model and
  its response for the selected type, component, and frequency.

Use forward modelling to build intuition for what a target looks like in the
data, and to generate training or test responses for the inversion page.

Inversion
---------

The **Inversion** page runs traditional, AI-neural, PINN, and hybrid inversions
in 1-D, 2-D, and 3-D.

.. figure:: ../../_static/applications/web/inversion-ai-neural.png
   :alt: The Inversion page showing AI training convergence
   :class: pycsamt-screenshot

   Inversion configured for a 2-D AI-neural (U-Net) run, with the training
   convergence curve on the Convergence tab.

* **Problem dimension** — for example *2-D Profile (U-Net)*.
* **Architecture** — the network or solver (for example *UNet2D
  (encoder–decoder)*).
* **Forward solver (training data)** — the physics used to generate training
  responses (for example *MT 1-D*).
* **Network Config** and **Training** — collapsible panels for network and
  optimisation settings.
* **Frequency Range** — the band used for the inversion.
* **Result tabs** — *Result*, *Convergence*, *Statistics*, *Log*, and
  *Data Fit* separate the model, the loss/RMS curve, the numbers, the run log,
  and the observed-vs-predicted fit.
* **Run Inversion** starts the run.

The header names the available families — *Traditional · AI Neural · PINN ·
Hybrid* — and the supported dimensionalities.  To browse the output of an
external solver run instead, use the Results View page below.

Results View
------------

The **Results View** page (*Inversion Results Viewer*) browses, inspects, and
exports inversion outputs from external solvers — **ModEM**, **Occam2D**, and
**MARE2DEM**.

.. figure:: ../../_static/applications/web/results-section.png
   :alt: The Results View page showing a ModEM resistivity section
   :class: pycsamt-screenshot

   Results View on the Section tab: a ModEM N–S resistivity section, with the
   solver auto-detected and the run metadata (iterations, RMS, grid) shown in
   the header.

* **Results Folder** — point at a solver output folder and **Load** it.
* **Solver** — *Auto-detect*, or force *ModEM (3-D)*, *Occam2D*, or
  *MARE2DEM*.  The panel confirms what was loaded (for example
  ``MODEM loaded · 74 iter · RMS=3.057``).
* **Tabs** — *Convergence*, *Section*, *Depth Map*, *All Profiles*,
  *Covariance*, *Response*, and *Pseudo*.
* **Display** — resistivity range and maximum depth.
* **Generate Plot** redraws with the current settings; **Export PNG** saves the
  figure.

.. figure:: ../../_static/applications/web/results-convergence.png
   :alt: The Results View Convergence tab showing RMS misfit versus iteration
   :class: pycsamt-screenshot

   The Convergence tab: RMS misfit versus iteration, with the best RMS and the
   RMS=1 target marked.

.. figure:: ../../_static/applications/web/results-response.png
   :alt: The Results View Response tab showing observed versus predicted data
   :class: pycsamt-screenshot

   The Response tab: observed-versus-predicted apparent resistivity and phase
   per station, with per-station RMS — the direct check of how well the model
   fits the data.

Interpretation
--------------

The **Interpretation** page turns processed data and models into geological and
hydrological interpretation products — 42 workflows across 9 categories.

.. figure:: ../../_static/applications/web/interpretation.png
   :alt: The Interpretation page with category buttons and a plot selector
   :class: pycsamt-screenshot

   Interpretation.  Pick a **category**, choose a **Plot / Analysis**, select a
   **data source**, set **display options**, and **Generate**.

* **Categories** — *Setup*, *Geology*, *Hydrology*, *Field*, *EM Diag*,
  *Uncert.*, *Advanced*, *Fusion*, and *Export*.
* **Plot / Analysis** — the specific workflow (for example *Constraint
  misfit — observed vs modelled field data*).
* **Data source** — *Raw (loaded)* or *Corrected*.
* **Display options** — figure size, colour map, and workflow-specific
  controls.
* **Generate** builds the product; **Export** (the Export category) saves
  interpretation outputs.

AI Agents
---------

The **AI Agents** page exposes pyCSAMT's agents — 36 agents across 10
categories, spanning LLM, processing, and inversion — in two modes.

**Agent Runner** lists agents you can configure and run directly.

.. figure:: ../../_static/applications/web/agents-runner.png
   :alt: The AI Agents page in Agent Runner mode running the Dimensionality agent
   :class: pycsamt-screenshot

   Agent Runner.  Search and pick an agent, scope it by line and station, set
   parameters, and **Run Agent**; the run history and last result (log, figure,
   summary) appear on the right.

* The agent list is grouped — **LLM** agents (Interpretation, Report, Code
  Generation, EDI Export), **Workflow** agents (Workflow Orchestrator,
  Pipeline, Batch Survey), and **Processing** agents (QC Quicklook,
  Dimensionality, Static Shift, and more).
* Selecting an agent shows its description, **Filter by line** and **Station
  filter** pickers, and its **Parameters**.
* **Run Agent** runs it; **Run History** logs each run, and **Last Result**
  shows the **Log**, **Figure**, and **Summary**.

**Chat** is a conversational assistant over the loaded survey.

.. figure:: ../../_static/applications/web/agents-chat.png
   :alt: The AI Agents Chat tab with a proposed plan and quick actions
   :class: pycsamt-screenshot

   The Chat tab.  Describe what you want in plain language; the assistant
   proposes a plan it can run, and **Quick Actions** cover common tasks.

Type a request such as *"Run quality control on all stations and identify bad
SNR"* and the assistant proposes a runnable plan.  **Quick Actions** —
*Run full QC*, *Correct static shift*, *Dimensionality*, *Prep inversion*,
*Interpret results* — trigger common workflows directly.

The agents need an LLM provider and API key for the LLM-backed features;
configure these in the **Settings** drawer (see :doc:`navigation`).  For a
dedicated, full-screen conversational surface, use the standalone
:doc:`Agent Master <../agent_master/index>` application.

Tools
-----

The **Tools** menu on the command bar collects standalone utilities that open
in their own drawer, independent of the current page.

.. figure:: ../../_static/applications/web/tools-menu.png
   :alt: The Tools dropdown grouped into analysis, conversion, geospatial, and more
   :class: pycsamt-screenshot

   The Tools menu, grouped by purpose.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Group
     - Tools
   * - Analysis
     - Strike Analyzer, EDI Validator
   * - Conversion & Export
     - Format Converter, Batch Export Plots
   * - Geospatial
     - Coordinate Transformer, Elevation Enrichment
   * - Visualisation
     - Station Response Inspector, Strike Profile Viewer, Phase Tensor Map
   * - Classification & Editing
     - Dimensionality Classifier, Frequency Editor
   * - Modelling
     - (modelling utilities)

Selecting a tool opens it in a drawer with an icon rail of the available tools
and the selected tool's controls and output.

.. figure:: ../../_static/applications/web/tools-strike-analyzer.png
   :alt: The Strike Analyzer tool drawer showing a consensus strike rose diagram
   :class: pycsamt-screenshot

   The **Strike Analyzer** tool: a consensus-strike rose diagram and a
   per-line strike distribution, computed from the currently loaded survey.

The **Batch Export Plots** tool is the fastest way to export many figures at
once; see :doc:`exports`.

How Pages Delegate To The Package
---------------------------------

None of these pages implement science of their own.  A control change updates a
Dash store, and the primary-action callback calls the matching
``pycsamt.app.desktop.controllers`` function — the same one the desktop GUI
calls and the same computation the Python API performs.  This is why a survey
processed in the web app can be picked up unchanged in the desktop app or in a
script, and why the results are reproducible outside the browser.

Next Steps
----------

* :doc:`exports` -- save figures, corrected EDIs, and inversion products.
* :doc:`maps_and_profiles` -- the spatial views that precede processing.
* :doc:`troubleshooting` -- when a page is blank, a run fails, or an agent
  errors.
