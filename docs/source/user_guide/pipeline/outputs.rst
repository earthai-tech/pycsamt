.. _pipeline-outputs:

Pipeline Outputs
================

Pipeline outputs are the files and programmatic objects produced by a
:class:`pycsamt.pipeline.Pipeline` run.  A normal run writes a reproducible
output directory containing processed EDI files, QC figures, reports, and the
exact pipeline configuration used for that run.  The same run also returns a
:class:`pycsamt.pipeline.PipelineResult` object for Python workflows.

Use this page when you need to know:

* where processed EDI files are written;
* how plot directories are named;
* what ``summary.txt`` and ``report.html`` contain;
* how ``pipeline.yaml`` supports reproducibility;
* how to disable selected output families;
* how to run completely in memory;
* how to inspect output paths from Python.

Output Philosophy
-----------------

The pipeline output system is built around three rules.

Outputs are explicit
    A run writes files only when an output directory is resolved.  Passing
    ``outdir=None`` to :meth:`pycsamt.pipeline.Pipeline.run` is an explicit
    in-memory run and produces no filesystem side effects.

Outputs are reproducible
    When output is enabled, pyCSAMT writes ``pipeline.yaml`` before processing
    starts.  This file captures the resolved step list and parameters.

Outputs are organized by responsibility
    Final processed data go under ``processed/``.  Per-step QC figures go
    under ``plots/``.  Human-readable summaries live at the output root.

Canonical Directory Tree
------------------------

A typical run writes this directory tree:

.. code-block:: text
   :linenos:

   results/basic_qc/
   |-- processed/
   |   |-- station001.edi
   |   |-- station002.edi
   |   `-- station003.edi
   |-- plots/
   |   |-- 01_notch/
   |   |   |-- nr_qc_harmonic_waterfall.png
   |   |   `-- nr_qc_snr_gain_profile.png
   |   |-- 02_drop_duplicates/
   |   |   `-- plot_coverage_quality_heatmap.png
   |   |-- 03_select_band/
   |   |   |-- plot_band_microstrips.png
   |   |   `-- plot_coverage_quality_heatmap.png
   |   `-- 05_qc_snapshot/
   |       |-- plot_qc_quicklook.png
   |       |-- plot_station_confidence_dashboard.png
   |       `-- plot_coverage_psection.png
   |-- pipeline.yaml
   |-- report.html
   `-- summary.txt

The exact plot names depend on the QC functions attached to each registered
step.  A step may produce zero, one, or several figures.

Output Directory Resolution
---------------------------

When using Python, the output root is resolved in this order:

1. explicit ``outdir`` passed to ``Pipeline.run``;
2. ``output_dir`` stored on a pipeline loaded from a config file;
3. global ``PYCSAMT_PIPE.output_root`` default, which is ``pipe_results``.

Examples:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_preset("basic_qc")

   # Explicit output directory.
   result = pipe.run(sites, outdir="results/basic_qc")

   # Use the pipeline or global default.
   result = pipe.run(sites)

   # Run in memory and write nothing.
   result = pipe.run(sites, outdir=None)

From the CLI, ``--out`` controls the output root:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/basic_qc

If ``--out`` is omitted, the CLI uses the pipeline config value when
available, otherwise the global default ``pipe_results``.

OutputDir
---------

The on-disk tree is managed internally by ``OutputDir``.  It creates the root,
``processed/``, and ``plots/`` directories when a run starts.

The default subdirectory names come from
:class:`pycsamt.api.pipe.PipelineAPIConfig`:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Config field
     - Default
     - Meaning
   * - ``output_root``
     - ``pipe_results``
     - Default output root when no explicit ``outdir`` or config
       ``output_dir`` is available.
   * - ``processed_subdir``
     - ``processed``
     - Subdirectory for final processed EDI files.
   * - ``plots_subdir``
     - ``plots``
     - Subdirectory for per-step QC figures.
   * - ``plot_dpi``
     - ``150``
     - DPI used when saving Matplotlib figures.
   * - ``plot_fmt``
     - ``png``
     - Plot format.  CLI choices are ``png``, ``pdf``, and ``svg``.
   * - ``report_formats``
     - ``("html", "txt")``
     - Report files to write when reports are enabled.
   * - ``save_intermediate``
     - ``False``
     - When true, save EDI snapshots after successful intermediate steps.

Configure global defaults:

.. code-block:: python
   :linenos:

   from pycsamt.api.pipe import configure_pipe

   configure_pipe(
       output_root="results/default_pipe",
       processed_subdir="edis_processed",
       plots_subdir="figures",
       plot_dpi=300,
       plot_fmt="pdf",
   )

Or override them temporarily:

.. code-block:: python
   :linenos:

   from pycsamt.api.pipe import PYCSAMT_PIPE

   with PYCSAMT_PIPE.context(plot_dpi=300, plot_fmt="svg"):
       result = pipe.run(sites, outdir="results/svg_qc")

Processed EDI Files
-------------------

Final processed EDI files are written under ``<outdir>/processed/`` when
``save_edis=True``.  The pipeline writes the site collection produced after
the last step.

CLI:

.. code-block:: console
   :linenos:

   # Write processed EDIs.
   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/basic_qc

   # Skip processed EDI export.
   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/no_edi \
       --no-edi

Python:

.. code-block:: python
   :linenos:

   result = pipe.run(
       sites,
       outdir="results/basic_qc",
       save_edis=True,
   )

   print(result.processed_paths)

If EDI export fails, pyCSAMT warns and returns an empty or partial
``processed_paths`` list.  The run may still contain valid step results and
QC plots, so inspect ``summary.txt`` and ``report.html`` before discarding the
whole run.

QC Figures
----------

Each registered step may declare one or more QC plot functions.  When
``save_plots=True`` and the step succeeds, pyCSAMT calls those functions and
saves the resulting Matplotlib figures under:

.. code-block:: text
   :linenos:

   <outdir>/plots/<step_index>_<step_label>/<qc_function_name>.<plot_fmt>

Example:

.. code-block:: text
   :linenos:

   results/basic_qc/plots/01_notch/nr_qc_harmonic_waterfall.png
   results/basic_qc/plots/01_notch/nr_qc_snr_gain_profile.png
   results/basic_qc/plots/05_qc_snapshot/plot_qc_quicklook.png

Control plotting from the CLI:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset publication_ready \
       --out results/publication_ready_pdf \
       --dpi 300 \
       --plot-fmt pdf

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/no_plots \
       --no-plots

Control plotting from Python:

.. code-block:: python
   :linenos:

   result = pipe.run(
       sites,
       outdir="results/no_plots",
       save_plots=False,
   )

   print(result.plots)

QC plot failures are skipped so that a successful transform does not become a
failed processing step only because a diagnostic figure could not be drawn.
If an expected figure is missing, check whether the step succeeded, whether
plots were disabled, and whether the current data contain the quantities
required by that QC function.

Reports
-------

When ``save_report=True``, the pipeline writes reports according to
``PYCSAMT_PIPE.report_formats``.  By default it writes both:

``summary.txt``
    Plain-text report for terminals, CI logs, quick review, and diffable
    processing notes.

``report.html``
    Self-contained HTML report with summary metadata, per-step cards, linked
    plot thumbnails, errors, parameters, and the pipeline configuration.

The text report contains:

* pipeline name;
* run timestamp;
* input and output site counts;
* total runtime;
* one row per step;
* status, step code, elapsed time, site counts, and plot count;
* error messages for failed steps;
* output directory pointers.

The HTML report contains the same run summary plus richer per-step cards:

* step label and registry code;
* step status;
* human-readable registry label;
* input and output site counts;
* elapsed time;
* parameter values;
* linked plot thumbnails when figures were saved;
* embedded pipeline YAML.

Disable report writing:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/no_report \
       --no-report

Python equivalent:

.. code-block:: python
   :linenos:

   result = pipe.run(
       sites,
       outdir="results/no_report",
       save_report=False,
   )

Write only one report format:

.. code-block:: python
   :linenos:

   from pycsamt.api.pipe import PYCSAMT_PIPE

   with PYCSAMT_PIPE.context(report_formats=("txt",)):
       result = pipe.run(sites, outdir="results/text_only")

Pipeline Snapshot
-----------------

Every output-enabled run saves:

.. code-block:: text
   :linenos:

   <outdir>/pipeline.yaml

This file is written before the main processing loop starts.  It is the
resolved pipeline configuration for the run and should be treated as the
source of truth for reproducing the processing sequence.

Reload a saved pipeline:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_yaml("results/basic_qc/pipeline.yaml")
   rerun = pipe.run(sites, outdir="results/basic_qc_rerun")

Use ``pipeline.yaml`` to:

* rerun a workflow after changing data loading code;
* review exactly which step parameters were active;
* compare two output directories;
* archive a processing recipe with reports and figures;
* debug a CLI run from Python.

In-Memory Runs
--------------

Pass ``outdir=None`` when you want a pure Python result without writing files.

.. code-block:: python
   :linenos:

   result = pipe.run(
       sites,
       outdir=None,
       save_plots=False,
       save_edis=False,
       save_report=False,
   )

   assert result.outdir is None
   assert result.processed_paths == []

This is useful in tests, notebooks, and exploratory workflows where the
processed ``Sites`` object is enough.

Intermediate EDI Snapshots
--------------------------

The global ``save_intermediate`` option can save EDI snapshots after
successful intermediate steps.  These snapshots are written inside the step's
plot directory:

.. code-block:: text
   :linenos:

   <outdir>/plots/03_select_band/edi_snapshot/

Enable snapshots temporarily:

.. code-block:: python
   :linenos:

   from pycsamt.api.pipe import PYCSAMT_PIPE

   with PYCSAMT_PIPE.context(save_intermediate=True):
       result = pipe.run(sites, outdir="results/debug_snapshots")

Use this option for debugging only.  It can create many files, especially for
large surveys or long pipelines.

PipelineResult
--------------

The return value of :meth:`pycsamt.pipeline.Pipeline.run` is a
``PipelineResult``.  It is the programmatic companion to the files written on
disk.

Important fields:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Field
     - Meaning
   * - ``sites_in``
     - Original input site collection.
   * - ``sites_out``
     - Site collection after the final step.
   * - ``step_results``
     - One ``StepResult`` per step.
   * - ``outdir``
     - Output root path, or ``None`` for in-memory runs.
   * - ``elapsed_sec``
     - Total wall-clock runtime.
   * - ``processed_paths``
     - Paths returned by final EDI export.
   * - ``pipeline_name``
     - Pipeline label.
   * - ``plots``
     - Derived list of all saved plot paths.
   * - ``ok``
     - ``True`` when every step completed without error.
   * - ``n_errors``
     - Number of failed steps.

Example:

.. code-block:: python
   :linenos:

   result = pipe.run(sites, outdir="results/basic_qc")

   print(result.summary())
   print(result.ok)
   print(result.outdir)
   print(result.processed_paths)

   for step_result in result.step_results:
       print(
           step_result.step_idx,
           step_result.step_name,
           step_result.step_code,
           step_result.ok,
           step_result.n_sites_in,
           step_result.n_sites_out,
           len(step_result.plots),
       )

Output Control Matrix
---------------------

The output flags are independent, but they only have an effect when an output
directory exists.

.. list-table::
   :header-rows: 1
   :widths: 30 26 44

   * - Control
     - Applies to
     - Result
   * - ``outdir=None``
     - Python API
     - No output directory and no files.
   * - ``--out DIR``
     - CLI
     - Sets the output root for the run.
   * - ``save_plots=False``
     - Python API
     - Do not generate or save QC figures.
   * - ``--no-plots``
     - CLI
     - Do not generate or save QC figures.
   * - ``save_edis=False``
     - Python API
     - Do not write final processed EDI files.
   * - ``--no-edi``
     - CLI
     - Do not write final processed EDI files.
   * - ``save_report=False``
     - Python API
     - Do not write ``summary.txt`` or ``report.html``.
   * - ``--no-report``
     - CLI
     - Do not write ``summary.txt`` or ``report.html``.
   * - ``plot_dpi`` / ``--dpi``
     - Figures
     - Controls saved figure resolution.
   * - ``plot_fmt`` / ``--plot-fmt``
     - Figures
     - Controls saved figure extension and Matplotlib output format.
   * - ``report_formats``
     - Reports
     - Selects ``html`` and/or ``txt`` when reports are enabled.

Recommended Output Layout
-------------------------

Keep raw data, configs, and outputs separate:

.. code-block:: text
   :linenos:

   project/
   |-- data/
   |   `-- raw_edis/
   |-- config/
   |   |-- basic_qc.yaml
   |   `-- publication_ready.yaml
   `-- results/
       |-- basic_qc/
       |-- noise_reduction/
       `-- publication_ready/

Run:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/raw_edis \
       --config config/basic_qc.yaml \
       --out results/basic_qc

This layout prevents accidental overwrites of raw files and makes it easy to
compare multiple processing strategies.

Comparing Two Runs
------------------

When comparing output directories, inspect the same files in each run:

.. code-block:: text
   :linenos:

   results/basic_qc/summary.txt
   results/noise_reduction/summary.txt
   results/basic_qc/report.html
   results/noise_reduction/report.html
   results/basic_qc/pipeline.yaml
   results/noise_reduction/pipeline.yaml

Useful comparisons:

* step status and error count in ``summary.txt``;
* site counts before and after each step;
* plot counts per step;
* parameter differences in ``pipeline.yaml``;
* visual differences in matching ``plots/<step>/`` folders;
* EDI differences under ``processed/``.

Stratagem Output Note
---------------------

:class:`pycsamt.pipeline.stratagem.StratagemPipeline` follows the normal
pipeline output tree.  When ``rename_basename`` is configured, it can also
copy or rename processed EDI files from ``processed/`` into a ``renamed/``
directory, or into a custom ``rename_dir``.

For raw Stratagem convenience workflows using ``run_stratagem_preset``, the
function writes a Stratagem-oriented output layout, including ``corrected``
and ``renamed`` directories under the requested output root.

Troubleshooting
---------------

No output directory was created
    In Python, check whether ``outdir=None`` was passed.  That is an explicit
    no-files run.  If using the CLI, check that the command reached the run
    phase and was not a ``--dry-run``.

``pipeline.yaml`` exists but reports are missing
    The pipeline writes ``pipeline.yaml`` before the main step loop.  Reports
    are written after the run only when ``save_report=True`` and the selected
    ``report_formats`` include ``html`` or ``txt``.

Plots are missing
    Check that ``--no-plots`` or ``save_plots=False`` was not used.  Also
    confirm that the step succeeded and that the step has registered QC plot
    functions.

Processed EDIs are missing
    Check that ``--no-edi`` or ``save_edis=False`` was not used.  If export
    failed, pyCSAMT warns and ``result.processed_paths`` may be empty.

Only some plots are present
    QC plot functions are skipped individually when they cannot produce a
    figure for the current site collection.  Inspect the report and run with
    verbose CLI output if needed.

The CLI dry-run output directory looks generic
    ``--dry-run`` reports the explicit ``--out`` value when supplied.  For a
    real run, output resolution still follows CLI ``--out``, config
    ``output_dir``, then the global default.

Output files were overwritten
    The output manager creates directories with ``exist_ok=True`` and writes
    standard filenames such as ``pipeline.yaml``, ``summary.txt``, and
    ``report.html``.  Use a new output root for each experimental run.

Related Pages
-------------

* :doc:`concepts` explains the run lifecycle and ``PipelineResult`` object.
* :doc:`configuration_files` explains ``output_dir`` and config-driven runs.
* :doc:`cli_pipe` explains CLI output flags such as ``--out``,
  ``--no-plots``, ``--no-edi``, ``--no-report``, ``--dpi``, and
  ``--plot-fmt``.
* :doc:`steps` explains which steps generate QC figures.
* :doc:`presets` explains how preset workflows produce comparable output
  directories.
