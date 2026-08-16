.. _pipeline-outputs:

Pipeline Outputs
================

Pipeline outputs are the on-disk files and in-memory objects produced by
:meth:`pycsamt.pipeline.Pipeline.run`.  A normal run creates one
:term:`pipeline output directory` containing a reproduced
:term:`canonical pipeline snapshot`, optional :term:`processed EDI` files,
optional :term:`QC figure` files, and text or HTML reports.  The same run also
returns a :term:`PipelineResult` for Python workflows.

The useful way to think about outputs is not simply "files were written".
A pipeline run transforms an input :term:`site collection`
:math:`S_0` through an ordered set of processing steps
:math:`T_1,\ldots,T_n`:

.. math::

   S_j = T_j(S_{j-1}; \theta_j), \qquad j=1,\ldots,n.

The returned ``PipelineResult`` carries :math:`S_0`, :math:`S_n`, and one
:term:`StepResult` for each transform.  The output directory records the same
run in a form that a reviewer can open without Python: ``pipeline.yaml`` for
the recipe, ``processed/`` for the final data state, ``plots/`` for per-step
diagnostics, and reports for a human-readable audit trail.

Use this page when you need to know where files are written, which flags
control each :term:`output artifact`, how to run in memory, and how to inspect
the result object from Python.

Output Lifecycle
----------------

An output-enabled run follows a fixed order:

1. Resolve the output root.
2. Create the root, ``processed/``, and ``plots/`` directories.
3. Write ``pipeline.yaml`` before the first processing step starts.
4. Run each step, recording one ``StepResult``.
5. Save QC figures for successful steps when plotting is enabled.
6. Optionally save intermediate EDI snapshots after successful steps.
7. Write final processed EDIs after the last step when EDI export is enabled.
8. Write ``summary.txt`` and/or ``report.html`` when reports are enabled.
9. Return the ``PipelineResult`` object.

This order matters for reproducibility.  If a later step fails and the run is
configured to continue, the output directory can still contain the original
pipeline snapshot, the completed step records, any figures produced before the
failure, and the final in-memory state used by the returned result.

Canonical Directory Tree
------------------------

A typical run writes this tree:

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

The exact plot names depend on the QC functions registered for each step.  A
step may produce zero, one, or several figures.  The step folder name is
deterministic:

.. math::

   f_j = \mathrm{plots}/\mathrm{format}("\%02d\_\%s", j, \ell_j),

where :math:`j` is the 1-based step index and :math:`\ell_j` is the configured
step label.  Stable labels therefore make output folders easier to compare
between processing runs.

Output Directory Resolution
---------------------------

Python resolves the output root with three cases:

.. math::

   d_{\mathrm{root}} =
   \begin{cases}
   \varnothing, & \text{if } \mathrm{outdir}=\texttt{None},\\
   \mathrm{outdir}, & \text{if an explicit path is passed},\\
   d_{\mathrm{pipeline}} \text{ or } d_{\mathrm{global}},
      & \text{if outdir is omitted}.
   \end{cases}

The first case is an explicit in-memory run.  The third case uses the
``output_dir`` stored on a pipeline loaded from a configuration file, or falls
back to the global ``PYCSAMT_PIPE.output_root`` default, ``pipe_results``.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> pipe = Pipeline.from_preset("basic_qc")
   >>> result = pipe.run(sites, outdir="results/basic_qc")
   >>> result.outdir
   WindowsPath('results/basic_qc')

   >>> result = pipe.run(sites)
   >>> result.outdir
   WindowsPath('pipe_results')

   >>> result = pipe.run(sites, outdir=None)
   >>> result.outdir is None
   True

From the CLI, ``--out`` controls the output root:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
       --preset basic_qc \
       --out results/basic_qc

If ``--out`` is omitted, the CLI uses the config file's ``output_dir`` when
available, otherwise the global default.

OutputDir And Defaults
----------------------

The on-disk tree is managed internally by ``OutputDir``.  It creates the root,
``processed/``, and ``plots/`` directories when a run starts.  The defaults
come from :class:`pycsamt.api.pipe.PipelineAPIConfig`:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Config field
     - Default
     - Meaning
   * - ``output_root``
     - ``pipe_results``
     - Default output root when neither an explicit ``outdir`` nor config
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
     - Report files to write when reports are enabled.  ``"dashboard"`` is
       a third, opt-in value — see `Dashboard Report`_.
   * - ``save_intermediate``
     - ``False``
     - Save EDI snapshots after successful intermediate steps.

Configure global defaults:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.api.pipe import configure_pipe, reset_pipe, PYCSAMT_PIPE
   >>> configure_pipe(
   ...     output_root="results/default_pipe",
   ...     processed_subdir="edis_processed",
   ...     plots_subdir="figures",
   ...     plot_dpi=300,
   ...     plot_fmt="pdf",
   ... )
   >>> PYCSAMT_PIPE.plot_fmt
   'pdf'
   >>> reset_pipe()

For temporary changes, use the context manager so the previous settings are
restored automatically:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.api.pipe import PYCSAMT_PIPE
   >>> with PYCSAMT_PIPE.context(plot_dpi=300, plot_fmt="svg"):
   ...     result = pipe.run(sites, outdir="results/svg_qc")
   >>> PYCSAMT_PIPE.plot_fmt
   'png'

Captured Minimal Run
--------------------

The following transcript was run against ``data/3edis`` with plots and EDI
export disabled so the output is small:

.. code-block:: pycon
   :linenos:

   >>> from pathlib import Path
   >>> from pycsamt.api import read_edis
   >>> from pycsamt.api.pipe import PYCSAMT_PIPE
   >>> from pycsamt.pipeline import Pipeline
   >>> sites = read_edis("data/3edis", strict=False).to_collection()
   >>> pipe = Pipeline.from_preset("basic_qc", pipeline_name="basic_qc")
   >>> with PYCSAMT_PIPE.context(show_progress=False, plot_dpi=72):
   ...     result = pipe.run(
   ...         sites,
   ...         outdir=".tmp/docs_outputs/basic_qc",
   ...         save_plots=False,
   ...         save_edis=False,
   ...         save_report=True,
   ...     )
   >>> print(result.summary())
   PipelineResult  'basic_qc'
     Sites   : 3 in -> 3 out
     Steps   : 5 (5 ok, 0 err)
     Time    : 0.95 s
     Plots   : 0
     Output  : .tmp\docs_outputs\basic_qc
   >>> sorted(p.name for p in Path(".tmp/docs_outputs/basic_qc").iterdir())
   ['pipeline.yaml', 'plots', 'processed', 'report.html', 'summary.txt']

Even with ``save_plots=False`` and ``save_edis=False``, the directories are
created because the run is output-enabled.  The disabled families simply do
not add files beneath them.

Processed EDI Files
-------------------

Final processed EDIs are written under ``<outdir>/processed/`` when
``save_edis=True``.  They represent :math:`S_n`, the site collection after the
last step, and should not be mixed with raw field EDIs.

CLI:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
       --preset basic_qc \
       --out results/basic_qc

   pycsamt pipe run data/3edis \
       --preset basic_qc \
       --out results/no_edi \
       --no-edi

Python:

.. code-block:: pycon
   :linenos:

   >>> result = pipe.run(
   ...     sites,
   ...     outdir="results/basic_qc",
   ...     save_edis=True,
   ... )
   >>> len(result.processed_paths) >= 0
   True

If EDI export fails, pyCSAMT warns and returns an empty or partial
``processed_paths`` list.  The transform may still have succeeded, so inspect
``result.ok``, the step results, and the reports before discarding the run.

QC Figures
----------

Each registered step may declare QC plotting functions.  When
``save_plots=True`` and the step succeeds, pyCSAMT calls those functions and
saves the returned :term:`Matplotlib figure` objects under:

.. code-block:: text
   :linenos:

   <outdir>/plots/<step_index>_<step_label>/<qc_function_name>.<plot_fmt>

Examples:

.. code-block:: text
   :linenos:

   results/basic_qc/plots/01_notch/nr_qc_harmonic_waterfall.png
   results/basic_qc/plots/01_notch/nr_qc_snr_gain_profile.png
   results/basic_qc/plots/05_qc_snapshot/plot_qc_quicklook.png

Control plotting from the CLI:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
       --preset publication_ready \
       --out results/publication_ready_pdf \
       --dpi 300 \
       --plot-fmt pdf

   pycsamt pipe run data/3edis \
       --preset basic_qc \
       --out results/no_plots \
       --no-plots

Control plotting from Python:

.. code-block:: pycon
   :linenos:

   >>> result = pipe.run(
   ...     sites,
   ...     outdir="results/no_plots",
   ...     save_plots=False,
   ... )
   >>> result.plots
   []

QC plot failures are skipped individually so that a successful transform does
not become a failed processing step only because a diagnostic figure could not
be drawn.  Missing figures usually mean one of four things: plots were
disabled, the step failed, the step has no registered QC functions, or the
current data do not contain the quantities required by that QC function.

Reports
-------

When ``save_report=True``, the pipeline writes reports according to
``PYCSAMT_PIPE.report_formats``.  By default it writes both ``summary.txt`` and
``report.html``.  A third, opt-in format, ``"dashboard"``, writes a richer
branded report with KPI stat tiles and charts — see `Dashboard Report`_
below.

``summary.txt``
    Plain-text report for terminals, CI logs, quick review, and diffable
    processing notes.

``report.html``
    Self-contained HTML report with run metadata, per-step cards, linked plot
    thumbnails, errors, parameter values, and embedded pipeline YAML.

A text report starts like this:

.. code-block:: text
   :linenos:

   pyCSAMT Pipeline Report
   Pipeline : basic_qc
   Run at   : 2026-07-18 20:43:38
   Sites    : 3 in -> 3 out
   Total    : 0.95s

   Step results
     #  Name                   Code     Status Time(s)       Sites  Plots

Disable report writing:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
       --preset basic_qc \
       --out results/no_report \
       --no-report

Python equivalent:

.. code-block:: pycon
   :linenos:

   >>> result = pipe.run(
   ...     sites,
   ...     outdir="results/no_report",
   ...     save_report=False,
   ... )

Write only one report format:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.api.pipe import PYCSAMT_PIPE
   >>> with PYCSAMT_PIPE.context(report_formats=("txt",)):
   ...     result = pipe.run(sites, outdir="results/text_only")

Dashboard Report
----------------

``report_formats`` also accepts ``"dashboard"``, a third, richer report
tier written to ``<outdir>/dashboard.html`` alongside — not instead of —
``summary.txt`` and ``report.html``.  Where ``report.html`` stays a plain,
cheap-to-render step-card list, the dashboard adds pyCSAMT's own logo and
brand colors, KPI stat tiles, and three native inline-SVG charts built
from the same ``step_results`` the other two reports already use.  There
is no external JavaScript or CDN dependency, so the file stays
self-contained: it opens directly from disk in a browser, or travels as a
single email attachment.

Enable it from Python by adding ``"dashboard"`` to ``report_formats``:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.api.pipe import PYCSAMT_PIPE
   >>> with PYCSAMT_PIPE.context(report_formats=("html", "txt", "dashboard")):
   ...     result = pipe.run(sites, outdir="results/with_dashboard")

or from the CLI with ``--dashboard``:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
       --preset basic_qc \
       --out results/with_dashboard \
       --dashboard

The dashboard adds four things beyond ``report.html``:

Stat tiles
    Steps ok/total, error count, total elapsed time, sites in → out, cache
    hit rate, and total figures generated — one glance at whether the run
    needs attention.

Step status timeline
    One colored block per step (green = OK, red = error), with a small
    gold-ringed dot marking any step whose result was replayed from the
    :doc:`step cache <caching>` instead of recomputed.  Hovering a block
    shows the step name, code, elapsed time, and cache status via a native
    SVG ``<title>`` — no script required.

Step duration bars
    One bar per step, in the brand blue by default.  A step at or above
    the run's own 80th-percentile elapsed time is drawn in gold instead,
    so an unusually slow step stands out without a fixed, run-independent
    threshold.

Site-count flow
    A two-series line chart of sites in vs. sites out per step, so a step
    that silently drops stations is visible at a glance rather than buried
    in a table column.

Every chart sits above a plain ``<table>`` restating the same per-step
numbers — the table is not an afterthought; it is the accessibility twin
required for anyone who cannot read the charts, and it is what a
:kbd:`Ctrl-F` search or a text diff actually matches against.

A run against three real EDIs, captured the same way as
`Captured Minimal Run`_ above:

.. code-block:: pycon
   :linenos:

   >>> with PYCSAMT_PIPE.context(
   ...     show_progress=False,
   ...     plot_dpi=72,
   ...     report_formats=("html", "txt", "dashboard"),
   ... ):
   ...     result = pipe.run(
   ...         sites,
   ...         outdir=".tmp/docs_outputs/basic_qc_dashboard",
   ...         save_plots=False,
   ...         save_edis=False,
   ...         save_report=True,
   ...     )
   >>> sorted(p.name for p in Path(".tmp/docs_outputs/basic_qc_dashboard").iterdir())
   ['dashboard.html', 'pipeline.yaml', 'plots', 'processed', 'report.html', 'summary.txt']

The stat-tiles block from that same run (SVG icon path data elided for
brevity — each tile embeds one small ``currentColor`` icon so it inherits
the surrounding text color in both light and dark mode):

.. code-block:: html
   :linenos:

   <div class="tiles">
     <div class="tile"><svg class="icon" ...></svg>
       <div class="label">Steps</div><div class="value">5/5 ok</div>
     </div>
     <div class="tile"><svg class="icon" ...></svg>
       <div class="label">Errors</div><div class="value">0</div>
     </div>
     <div class="tile"><svg class="icon" ...></svg>
       <div class="label">Total time</div><div class="value">28.50s</div>
     </div>
     ...
   </div>

The dashboard's palette is not an arbitrary restyle: the brand blue/orange
pair used for the site-count-flow chart, and the good/warning/critical
status colors used throughout, were both checked against pyCSAMT's own
light and dark surfaces with the project's color-blindness and contrast
validator before being wired in, and both passed without substitution.
Where a status color's contrast is intentionally low against a light
surface (the gold "slow step" / "cached" marker), the mitigation is a
visible caption or tooltip beside it, never color alone.

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

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> pipe = Pipeline.from_yaml("results/basic_qc/pipeline.yaml")
   >>> rerun = pipe.run(sites, outdir="results/basic_qc_rerun")
   >>> rerun.pipeline_name
   'basic_qc'

Use ``pipeline.yaml`` to rerun a workflow, review active parameters, compare
two output directories, archive a processing recipe with reports and figures,
or debug a CLI run from Python.

In-Memory Runs
--------------

Pass ``outdir=None`` when you want a pure Python result without writing files:

.. code-block:: pycon
   :linenos:

   >>> result = pipe.run(
   ...     sites,
   ...     outdir=None,
   ...     save_plots=False,
   ...     save_edis=False,
   ...     save_report=False,
   ... )
   >>> result.outdir is None
   True
   >>> result.processed_paths
   []

This is useful in tests, notebooks, and exploratory workflows where the
processed ``Sites`` object is enough.  The mathematical run state still exists
in memory as :math:`S_n`; only the filesystem projection is disabled.

Intermediate EDI Snapshots
--------------------------

The global ``save_intermediate`` option can save EDI snapshots after
successful intermediate steps.  These snapshots are written inside the step's
plot directory:

.. code-block:: text
   :linenos:

   <outdir>/plots/03_select_band/edi_snapshot/

Enable snapshots temporarily:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.api.pipe import PYCSAMT_PIPE
   >>> with PYCSAMT_PIPE.context(save_intermediate=True):
   ...     result = pipe.run(sites, outdir="results/debug_snapshots")

Use this option for debugging only.  It can create many files, especially for
large surveys or long pipelines.

PipelineResult And StepResult
-----------------------------

The return value of ``Pipeline.run`` is the programmatic companion to the
files on disk.  Important ``PipelineResult`` fields are:

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

Inspect the run:

.. code-block:: pycon
   :linenos:

   >>> print(result.summary())  # doctest: +ELLIPSIS
   PipelineResult  'basic_qc'
     Sites   : 3 in -> 3 out
     Steps   : 5 (5 ok, 0 err)
     ...
   >>> result.ok
   True
   >>> result.n_errors
   0
   >>> [(sr.step_idx, sr.step_name, sr.step_code, sr.ok)
   ...  for sr in result.step_results]
   [(1, 'notch', 'NR001', True), ..., (5, 'qc_snapshot', 'QC001', True)]

Each ``StepResult`` also stores the parameters passed to the step, elapsed
time, input and output site counts, saved plot paths, and any captured error.
For review, the most important relation is

.. math::

   \mathrm{ok}_{\mathrm{run}}
   = \bigwedge_{j=1}^{n} \mathrm{ok}_j,

so ``result.ok`` is true only when every step result is true.

Output Control Matrix
---------------------

The output flags are independent, but they only write files when an output
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
     - Selects ``html`` and/or ``txt`` when reports are enabled; add
       ``dashboard`` for the richer branded report.
   * - ``--dashboard``
     - CLI
     - Adds ``dashboard.html`` for this run without disabling the default
       ``html``/``txt`` reports.

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
compare multiple processing strategies.  Use a new output root for each
experimental branch, especially when changing step order, frequency bands,
notch parameters, static-shift corrections, or inversion-preparation logic.

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

Useful comparisons include step status and error count in ``summary.txt``,
site counts before and after each step, plot counts per step, parameter
differences in ``pipeline.yaml``, visual differences in matching
``plots/<step>/`` folders, and EDI differences under ``processed/``.

Stratagem Output Note
---------------------

:class:`pycsamt.pipeline.stratagem.StratagemPipeline` follows the normal
pipeline output tree.  When ``rename_basename`` is configured, it can also
copy or rename processed EDI files from ``processed/`` into a ``renamed/``
directory, or into a custom ``rename_dir``.

For raw Stratagem convenience workflows using ``run_stratagem_preset``, the
function writes a Stratagem-oriented output layout, including ``corrected`` and
``renamed`` directories under the requested output root.

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
  ``--no-plots``, ``--no-edi``, ``--no-report``, ``--dpi``, ``--plot-fmt``,
  and ``--dashboard``.
* :doc:`steps` explains which steps generate QC figures.
* :doc:`presets` explains how preset workflows produce comparable output
  directories.
