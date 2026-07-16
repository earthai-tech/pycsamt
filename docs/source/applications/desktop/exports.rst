.. _applications-desktop-exports:

Exports And Reproducibility
===========================

The desktop app produces several kinds of outputs: figures, re-exported EDI
files, station metadata tables, pipeline configuration files, solver input
folders, interpretation products, logs, and the desktop session state.  Treat
these as different artifacts.  A figure explains a decision; an EDI folder
stores a processed data product; a JSON configuration makes the workflow
repeatable.

The guiding rule is simple: every exported result should carry enough context
for another user to understand what data state it came from and which settings
produced it.

Think In Data States
--------------------

Before exporting, name the data state you are working from.  Most confusion in
desktop projects comes from mixing raw, corrected, recomputed, and pipeline
outputs in the same folder or figure caption.

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Data State
     - Meaning
     - Export Habit
   * - ``raw``
     - Original files loaded from the field or source archive.
     - Keep immutable; do not overwrite or mix with processed EDIs.
   * - ``qc_reviewed``
     - Raw data with QC figures and notes, but no data values changed.
     - Save figures, station inventory, and review notes.
   * - ``corrected``
     - Active survey after a correction was committed to the desktop session.
     - Export corrected EDIs plus before/after figures and parameters.
   * - ``recomputed``
     - EDIs rewritten by the recompute tool after rotation, trimming, filling,
       or impedance-derived value updates.
     - Save the recompute manifest and reload the output before using it.
   * - ``pipeline``
     - Output from an automated processing chain.
     - Save pipeline JSON, logs, preview figures, and exported EDIs together.
   * - ``inversion_input``
     - Files prepared for an external or internal inversion run.
     - Save input files with the QC/correction evidence that defines them.

Use these state names in folder names and captions.  A file named
``L30_corrected_static_shift_ama_rhoa_phase.png`` is easier to audit than
``profile_final.png``.

Export Types
------------

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Output
     - Desktop Surface
     - Typical Use
   * - Figure files
     - **More > Export Figure**, panel export buttons
     - Save the current map, profile, QC, correction, forward, or inversion
       plot for reports and processing notes.
   * - Batch figures
     - **Tools > Batch Export Plots...**
     - Save every open canvas figure into one folder.
   * - Survey metadata
     - **Tools > Format Converter...**
     - Export station metadata to CSV or JSON.
   * - Re-exported EDI files
     - **Tools > Format Converter...**, correction tools, pipeline export
     - Save a processed or converted survey as station EDI files.
   * - Recomputed EDIs
     - **Tools > Recompute EDIs...**
     - Rotate tensors, trim frequency range, fill missing values, recompute
       resistivity/phase, and write an optional manifest CSV.
   * - Pipeline configuration
     - **Processing Pipeline > Save config**
     - Preserve method and parameter choices for a repeatable processing chain.
   * - Inversion input/output
     - **Inversion Wizard**
     - Store solver input files, logs, and results in a run-specific working
       directory.
   * - Interpretation products
     - **Interpretation Studio**
     - Export Oasis Montaj XYZ, LAS, CSV tables, and VTK grids.
   * - Desktop session
     - **File > Save Session**
     - Restore theme, recent files, selected station, window geometry, and
       preferences.

Figure Export
-------------

Most desktop plots are matplotlib figures.  Use **More > Export Figure** or a
panel's own **Export** button to save the active figure.  The export dialog
delegates to ``Figure.savefig`` and supports:

.. list-table::
   :header-rows: 1
   :widths: 18 28 54

   * - Format
     - Type
     - Best Use
   * - ``PNG``
     - Raster, lossless
     - Screenshots, documentation, quick sharing, and web pages.
   * - ``TIFF``
     - Raster, high resolution
     - Journal or report workflows that require high-DPI raster figures.
   * - ``PDF``
     - Vector
     - Reports, print, and figures that need sharp text and lines.
   * - ``SVG``
     - Vector
     - Web documentation or editing in vector graphics tools.
   * - ``EPS``
     - Vector, legacy
     - Older publishing pipelines that still request EPS.

The dialog defaults to 300 DPI, which is a good report setting for raster
exports.  Use 150 DPI for quick review images, 300 DPI for most reports, and
600 DPI only when a publisher or print workflow needs it.  Vector formats are
usually better for line plots, rose diagrams, and profile curves because labels
remain sharp at any zoom level.

When exporting a figure for a processing decision, put the decision in the
filename:

.. code-block:: text

   qc/L30_coverage_before_frequency_edit.png
   qc/L30_snr_after_noise_filter.png
   corrections/L30_static_shift_before_after.pdf
   strike/L30_strike_stability_bands.svg

Avoid generic names such as ``figure1.png``.  They become untraceable as soon
as several windows are open.

Figure Naming Pattern
---------------------

A useful figure name usually contains:

.. code-block:: text

   <line>_<data-state>_<view>_<component-or-band>_<purpose>.<ext>

Examples:

.. code-block:: text

   L30_raw_station_map_geometry_check.png
   L30_qc_reviewed_snr_before_filtering.png
   L30_corrected_static_shift_ama_before_after.pdf
   L30_recomputed_rot30_profile_xy_1s.svg
   L30_pipeline_qc_ss_rotate_coverage_after.png

Use short names, but include the part that would otherwise be forgotten:
line, data state, view, component or band, and purpose.

Batch Plot Export
-----------------

Use **Tools > Batch Export Plots...** when several panel windows are open and
you want to capture the current desktop evidence in one folder.  The batch
exporter collects visible matplotlib canvases from open panel windows and saves
them as ``<label>.<format>``.

Batch export supports ``PNG``, ``PDF``, ``SVG``, ``EPS``, and ``TIFF``.  It
also lets you choose DPI and whether the background should be transparent.
Transparent backgrounds are useful for slide decks, but they can make dark
theme labels hard to read on a light page.  For documentation and reports,
prefer a normal opaque background unless you know the final page color.

Use batch export after a review pass:

1. Open the map, profile, QC, and correction windows you want to document.
2. Set each window to the exact plot and station/band selection you want.
3. Run batch export to a dated folder.
4. Rename or annotate key figures if the automatic label is not specific
   enough.

Batch export saves visual evidence.  It does not save the processing settings
that produced corrected data.  Pair it with pipeline JSON, correction notes,
or exported EDIs when reproducibility matters.

Minimum Evidence Set
--------------------

When a processing decision changes the data state, save at least:

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - Decision
     - Minimum Figures
     - Minimum Machine-Readable Files
   * - Frequency or band selection
     - Coverage and SNR before/after.
     - Pipeline JSON or notes listing selected bands.
   * - Static-shift correction
     - Before/after station curves and pseudosection.
     - Corrected EDIs or correction parameters.
   * - Rotation
     - Strike/stability evidence and profile comparison.
     - Recomputed EDIs plus manifest or pipeline JSON.
   * - Source-effect correction
     - Diagnostic figure and before/after response comparison.
     - Corrected EDIs in a clearly named output folder.
   * - Inversion setup
     - QC summary, profile/map figure, and starting model figure.
     - Solver input files, run config, and log.

This minimum set is what lets another user understand both what changed and
why it changed.

Survey Metadata Export
----------------------

Use **Tools > Format Converter...** to export station metadata from the loaded
survey.  The converter can write:

* ``survey_stations.csv`` for spreadsheet review and quick station inventory;
* ``survey_stations.json`` for structured downstream tools;
* EDI re-exports when the loaded station objects expose an EDI writer.

The CSV and JSON metadata include station name, latitude, longitude, frequency
count, minimum and maximum period, and whether impedance error data are
present.  Use these files as an inventory, not as the full MT/CSAMT data
product.

Metadata export is useful before sharing a survey because it lets a reviewer
check station count, coordinate coverage, and frequency availability without
opening every EDI file.

Metadata Review Uses
--------------------

Station metadata exports are especially useful for:

* checking whether all expected stations are present;
* comparing field sheets with loaded station IDs;
* reviewing coordinate or elevation problems in a spreadsheet;
* documenting frequency coverage before processing;
* sharing a lightweight inventory when the EDI files are too large or private;
* giving inversion collaborators a quick summary of the profile.

Metadata are not a replacement for EDIs or QC figures.  Treat them as the
project index.

EDI Data Products
-----------------

EDI outputs are the data products that other tools can reload.  The desktop
can write EDI files from several places:

* **Format Converter** can re-export the currently loaded survey when station
  objects support EDI writing.
* **Data Corrections** can export Stratagem-corrected EDIs from the Stratagem
  correction path.
* **Processing Pipeline** writes processed EDIs during its final **Export**
  step using ``pycsamt.emtools.export_edis``.
* **Recompute EDIs** writes a transformed EDI set with optional tensor
  rotation, frequency trimming, missing-value fill, and resistivity/phase
  recomputation.

Always keep raw EDIs separate from corrected or recomputed EDIs.  A practical
layout is:

.. code-block:: text

   raw_edi/
   corrected_edi/static_shift_ama/
   corrected_edi/source_effect_test/
   recomputed_edi/rotated_to_strike/
   pipeline/L30_qc_ss_rotate/exported_edi/

After exporting an EDI product, reload it in a fresh desktop session or clear
mental state by reopening the folder.  Check the station count, map geometry,
profile curves, and QC coverage before using the exported product for
inversion.

EDI Product Safety Rules
------------------------

Use these rules whenever the desktop writes EDI files:

* never write corrected EDIs into the raw input folder;
* keep one output folder per correction, recompute, or pipeline attempt;
* include method names in folder names, not only dates;
* keep failed or experimental outputs separate from accepted outputs;
* reload exported EDIs before using them as inversion input;
* keep the figures and parameter files beside the exported EDIs;
* do not rename station files in a way that breaks station identity.

An EDI folder is a data product.  Once another program consumes it, the folder
name and nearby evidence become the processing history.

Recomputed EDIs
---------------

The **Recompute EDIs...** tool is for controlled data-product generation.  It
can use the loaded stations or a selected EDI file/folder as input, then write
a new EDI set with chosen transformations.

Important options:

* **Rotation angle** -- rotates impedance and/or tipper components when a
  strike decision has already been justified.
* **Frequency limits** -- trims the output to a trusted frequency band.
* **Fill missing values** -- fills missing data according to the selected
  strategy; document this choice carefully.
* **Recompute resistivity/phase** -- updates derived values from impedance.
* **Manifest CSV** -- writes a station-level record of the operation.
* **Filename template** -- controls output names, for example
  ``{source_stem}.edi``.

Use recomputed EDIs when you need a clean, explicit handoff to external
software.  Do not use this tool as a hidden correction step; save the manifest
and the diagnostic figures that justify every transformation.

Recompute Audit Checklist
-------------------------

Before accepting a recomputed EDI folder:

* confirm the output folder name states the operation, for example
  ``rotated_30deg_trim_1e-3_1e2hz``;
* save the manifest CSV when available;
* reload the recomputed folder in the desktop;
* compare station count with the source folder;
* inspect one early, one middle, and one late station profile;
* verify frequency limits, rotation, or filled values are visible as expected;
* keep the source EDIs unchanged and nearby.

If the recomputed folder cannot pass this audit, keep it as an experiment and
do not use it as the accepted inversion input.

Pipeline Configurations
-----------------------

The processing pipeline can save and load its configuration as JSON.  The JSON
stores each step ID, step name, selected method, and method parameters.  It
does not replace the exported EDIs; it explains how those EDIs were generated.

Save the pipeline JSON beside the pipeline output:

.. code-block:: text

   pipeline/
     L30_qc_ss_rotate.json
     L30_qc_ss_rotate.log
     exported_edi/
     figures/

When you rerun a saved configuration, confirm that step 1 is using the
intended input source.  A valid JSON file can still produce the wrong output if
it is run against a different active survey.

Pipeline Export Package
-----------------------

A pipeline export is complete when it contains:

* the pipeline JSON;
* the run log;
* exported EDIs, when the final step writes data;
* preview or summary figures for important steps;
* a short note naming the input data state;
* a clear output folder name that includes line and major methods.

For example:

.. code-block:: text

   pipeline/
     L30_raw_to_qc_ss_rotate/
       L30_raw_to_qc_ss_rotate.json
       L30_raw_to_qc_ss_rotate.log
       figures/
       exported_edi/
       README.md

The ``README.md`` can be short.  One or two sentences about input state,
major choices, and accepted output are enough to make the folder easier to
review later.

Inversion Exports
-----------------

The inversion wizard writes into the selected working directory.  For external
engines such as Occam2D, ModEM, and MARE2DEM, input-file generation is already
an export worth preserving.  Keep one folder per run attempt:

.. code-block:: text

   inversion/
     L30_occam2d_corrected_ss_v01/
       inputs/
       logs/
       results/
       figures/

Use run names that identify line, data state, engine, and attempt number.  A
folder such as ``modem_run`` is not enough once several corrected data states
exist.

Before archiving an inversion folder, save:

* the input files generated by the builder;
* the solver log or desktop run log;
* the starting model or forward-model parameters;
* result plots;
* the QC and correction figures that define the input data state.

Inversion Run Naming
--------------------

Use run names that can survive months of project history:

.. code-block:: text

   <line>_<engine>_<data-state>_<model-or-band>_<attempt>

Examples:

.. code-block:: text

   L30_occam2d_corrected_ss_v01
   L30_mare2dem_recomputed_rot30_v02
   L18_modem_pipeline_qc_ss_rotate_v01

Avoid names such as ``test``, ``final``, ``new_final``, or ``run2``.  They
only make sense on the day they are created.

Interpretation Exports
----------------------

The Interpretation Studio exports derived geological products.  Available
exports include:

* **XYZ** for Oasis Montaj style workflows;
* **LAS** for station or borehole-style logs;
* **CSV** for flat interpretation tables;
* **VTK** for gridded model visualization.

These files are downstream interpretation products.  Save them under a folder
that points back to the inversion run and corrected EDI set that produced
them:

.. code-block:: text

   interpretation/
     from_L30_occam2d_corrected_ss_v01/
       L30_interpretation.csv
       L30_model.vtk
       station_logs/

Interpretation outputs should not be separated from their parent inversion
context.  A VTK grid or CSV classification is only meaningful if the reviewer
can find the model, inversion run, corrected data state, and assumptions that
produced it.

Session And Settings
--------------------

Use **File > Save Session** to persist the desktop session.  The session is
stored at ``~/.pycsamt/session.json`` and includes theme, recent files, last
data directory, selected station, frequency preferences, window geometry,
solver paths, default inversion work directory, logging level, map tile
provider, and API key setting.

Session state is useful for convenience, but it is not a project archive.  It
does not guarantee that another user can reconstruct a processing product.
For project reproducibility, save exported data, figures, pipeline JSON, and
inversion/interpretation folders in the project directory.

Settings profiles are separate JSON snapshots of application/API settings.
Use them when you want to move a configuration between machines or preserve a
known set of plotting and processing defaults.

What Not To Rely On
-------------------

Do not rely on these as the only record of a project:

* screenshots copied from the screen;
* the desktop session file alone;
* recent-file history;
* a solver result folder without input files;
* corrected EDIs without before/after figures;
* pipeline EDIs without the JSON configuration;
* an interpretation export without its parent inversion run.

These items are useful, but alone they do not explain the processing chain.
Pair visual, data, and configuration outputs whenever a result matters.

Recommended Project Layout
--------------------------

A complete desktop project folder can look like this:

.. code-block:: text

   L30_project/
     raw_edi/
     metadata/
       survey_stations.csv
     qc/
       L30_coverage_before_processing.png
       L30_snr_before_processing.png
     corrections/
       static_shift_ama/
         figures/
         corrected_edi/
     strike/
       L30_strike_analysis.pdf
       L30_strike_stability_bands.svg
     forward/
       geothermal_starting_model.json
       synthetic_response.png
     pipeline/
       L30_qc_ss_rotate.json
       L30_qc_ss_rotate.log
       exported_edi/
     inversion/
       L30_occam2d_corrected_ss_v01/
     interpretation/
       from_L30_occam2d_corrected_ss_v01/
     session/
       notes.md

The exact folder names can change, but the structure should separate raw data,
processed data, diagnostic figures, configurations, solver runs, and final
interpretation products.

Export Checklist
----------------

Before sharing or archiving a desktop result, check that you have:

* raw input EDIs or a link to their immutable location;
* exported corrected/recomputed EDIs if the active data were changed;
* QC figures before and after major corrections;
* strike, dimensionality, or topography figures when they affected modelling;
* pipeline JSON and logs for automated workflows;
* inversion input files, logs, and result figures;
* interpretation exports with a clear parent inversion run;
* a short note explaining the data state and major assumptions.

This is the difference between "I saved the plot" and "I can reproduce the
result."  The second one is the one future you will thank present you for.
