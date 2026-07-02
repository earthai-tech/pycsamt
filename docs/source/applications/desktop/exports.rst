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
