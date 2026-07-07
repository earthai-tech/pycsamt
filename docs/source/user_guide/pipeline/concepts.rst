.. _pipeline-concepts:

Pipeline Concepts
=================

The pyCSAMT pipeline system is a reproducible processing engine for MT, AMT,
and CSAMT site collections.  A pipeline is an ordered list of registered
processing steps.  Each step receives the current site collection, transforms
it or inspects it, and passes the result to the next step.

At the end of a run, pyCSAMT returns a :class:`pycsamt.pipeline.PipelineResult`
and, when an output directory is enabled, writes a reproducible directory
containing processed EDI files, QC figures, reports, and the exact pipeline
YAML used for the run.

Why Pipelines Exist
-------------------

Field data processing often grows from a notebook into a sequence of repeated
operations: remove power-line harmonics, drop duplicate frequencies, trim to a
band, align stations onto a common frequency grid, correct static shift, make
quality-control plots, and prepare the data for inversion.

Pipelines make that sequence explicit.  They help you:

* repeat the same workflow across multiple survey lines;
* save the processing recipe alongside the outputs;
* inspect which step changed the data;
* run the same workflow from Python and the CLI;
* generate reports and QC figures in a predictable directory tree;
* share processing decisions with collaborators.

Core Objects
------------

The pipeline package is exposed through :mod:`pycsamt.pipeline`.  Users should
normally import from this public namespace:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline, Step

The main objects are:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Object
     - Role
   * - :class:`pycsamt.pipeline.Pipeline`
     - Ordered sequence of ``(label, Step)`` entries.  It can be built from
       code, a preset, or a YAML/JSON/Python config file.
   * - :class:`pycsamt.pipeline.Step`
     - Configured wrapper around one registered processing operation.  It
       binds a registry code such as ``"NR001"`` to parameter overrides.
   * - ``StepSpec``
     - Registry descriptor for a step: code, registry name, category,
       function path, default parameters, QC plot functions, and whether the
       step transforms the site collection.
   * - ``Preset``
     - Named, ordered collection of steps for common workflows such as
       ``"basic_qc"``, ``"full_processing"``, or ``"publication_ready"``.
   * - :class:`pycsamt.pipeline.StepResult`
     - One per executed step.  Records timing, input/output station count,
       saved plot paths, parameters, and any stored exception.
   * - :class:`pycsamt.pipeline.PipelineResult`
     - Returned by ``Pipeline.run``.  Contains the original sites, final
       sites, all step results, output directory, processed file paths, total
       runtime, and status helpers.

The Mental Model
----------------

A pipeline run is a left-to-right data flow:

.. code-block:: text
   :linenos:

   input Sites
      |
      v
   Step 1: transform or inspect
      |
      v
   Step 2: transform or inspect
      |
      v
   ...
      |
      v
   final Sites + PipelineResult + optional files

Most steps transform the site collection and return a new or modified
collection.  Diagnostic-only steps run a function for side effects or checks
and pass the input collection through unchanged.

Step Registry
-------------

Pipeline steps are not arbitrary strings.  They are registered in the pipeline
step registry.  Each registered step has:

* a short code, for example ``NR001``;
* a registry name, for example ``notch_powerline``;
* a category, for example ``noise_removal`` or ``frequency``;
* default parameters;
* a transform function;
* optional QC plot functions;
* a ``returns_sites`` flag.

The code and registry name both identify the same step:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Step

   notch_by_code = Step("NR001", mains_hz=50.0)
   notch_by_name = Step("notch_powerline", mains_hz=50.0)

Use the code form in pipeline configuration files and reports.  Use the
registry name when it improves readability in exploratory Python code.

Discover available steps from Python:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   print(Pipeline.catalogue())
   print(Pipeline.catalogue("frequency"))
   print(Pipeline.step_info("NR001"))

Discover the same information from the CLI:

.. code-block:: console
   :linenos:

   pycsamt pipe steps
   pycsamt pipe steps --category frequency
   pycsamt pipe steps --info NR001

Configured Steps
----------------

A :class:`pycsamt.pipeline.Step` combines a registry entry with user parameter
overrides.  The registry defaults are merged with your overrides.

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Step

   # Uses registry defaults except mains_hz.
   step = Step("NR001", mains_hz=60.0)
   print(step)

If the registry default for ``NR001`` includes ``n_harm`` and ``tol_hz``, the
configured step still carries those defaults.  You only need to provide the
values that should change for the workflow.

Pipeline Structure
------------------

A pipeline stores steps as ``(label, Step)`` tuples.  The label is the name of
this occurrence in this workflow.  It appears in printed summaries, output
subdirectories, reports, and CLI slicing options.

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline, Step

   pipe = Pipeline(
       [
           ("notch", Step("NR001", mains_hz=50.0)),
           ("select_band", Step("FREQ001", band_hz=(0.001, 10000.0))),
           ("align_grid", Step("FREQ004")),
           ("qc_snapshot", Step("QC001")),
       ],
       name="first_qc",
   )

   print(pipe)

Labels should be short, stable, and meaningful.  Prefer ``select_amt_band`` or
``correct_ss`` over vague labels such as ``step1``.

Building A Pipeline
-------------------

There are four common ways to build a pipeline.

Build directly in Python:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline, Step

   pipe = Pipeline([
       ("notch", Step("NR001")),
       ("drop_duplicates", Step("FREQ002")),
       ("select_band", Step("FREQ001")),
       ("qc_snapshot", Step("QC001")),
   ])

Build from a preset:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_preset("basic_qc")

Build from a config file:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_yaml("config/basic_qc.yaml")
   pipe = Pipeline.from_json("config/basic_qc.json")
   pipe = Pipeline.from_py("config/basic_qc.py")

Build from the CLI:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc
   pycsamt pipe run data/edis --config config/basic_qc.yaml
   pycsamt pipe run data/edis --steps NR001,FREQ002,FREQ001,QC001

Configuration files are documented in :doc:`configuration_files`.

Presets
-------

Presets are named pipelines for common processing intentions.  They are useful
when you want a known baseline without writing every step manually.

Examples include:

``basic_qc``
    Minimal denoising and frequency cleanup.  Good for first-pass inspection.

``noise_reduction``
    Stacked noise-removal chain for high-EMI environments.

``full_processing``
    Standard chain for noise removal, frequency cleanup, skew gate,
    static-shift correction, and strike rotation.

``publication_ready``
    A longer chain for publication-quality outputs.

Use a preset directly:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_preset("publication_ready")

Or generate an editable config from a preset:

.. code-block:: console
   :linenos:

   pycsamt pipe init --preset publication_ready \
       --name line22_publication \
       --output config/line22_publication.yaml

See :doc:`presets` for the dedicated preset guide.

Mutable Until Run
-----------------

A pipeline can be edited before it starts running:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline, Step

   pipe = Pipeline.from_preset("full_processing")
   pipe.remove("mask_skew")
   pipe.append("final_qc", Step("QC001"))
   pipe.replace("notch", Step("NR001", mains_hz=60.0))

During ``Pipeline.run``, the step list is protected from mutation.  This
prevents accidental changes while step results and reports are being produced.

Run Lifecycle
-------------

Calling ``Pipeline.run`` performs these operations in order:

1. Resolve the runtime configuration.
2. Resolve the output directory.
3. Save a canonical ``pipeline.yaml`` snapshot when output is enabled.
4. For each configured step:

   * count input sites;
   * run the step transform;
   * handle errors according to ``on_step_error``;
   * generate and save QC plots when enabled;
   * optionally save intermediate EDI snapshots;
   * create a :class:`pycsamt.pipeline.StepResult`.

5. Write final processed EDI files when ``save_edis=True``.
6. Write HTML and/or text reports when ``save_report=True``.
7. Return a :class:`pycsamt.pipeline.PipelineResult`.

Example:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis
   from pycsamt.pipeline import Pipeline

   survey = read_edis("data/edis", strict=False)
   pipe = Pipeline.from_preset("basic_qc")

   result = pipe.run(
       survey.to_collection(),
       outdir="results/basic_qc",
       save_plots=True,
       save_edis=True,
       save_report=True,
   )

   print(result.summary())

Output Resolution
-----------------

The output directory is resolved in this order:

1. explicit ``outdir`` passed to ``Pipeline.run``;
2. ``output_dir`` stored on a pipeline loaded from a config file;
3. global ``PYCSAMT_PIPE.output_root``.

Passing ``outdir=None`` is an explicit opt-out: the pipeline runs in memory
and writes no output files.

.. code-block:: python
   :linenos:

   # Write to the config/default output directory.
   result = pipe.run(sites)

   # Override output directory for this run.
   result = pipe.run(sites, outdir="results/experiment_01")

   # In-memory run: no files are written.
   result = pipe.run(sites, outdir=None)

Output Directory Contract
-------------------------

When output is enabled, pyCSAMT writes a predictable run directory:

.. code-block:: text
   :linenos:

   results/basic_qc/
   |-- processed/
   |   `-- *.edi
   |-- plots/
   |   |-- 01_notch/
   |   |-- 02_drop_duplicates/
   |   `-- ...
   |-- pipeline.yaml
   |-- report.html
   `-- summary.txt

``pipeline.yaml``
    Reproducible snapshot of the exact pipeline that was run.

``processed/``
    Final processed EDI files when ``save_edis=True``.

``plots/``
    QC figures generated after individual steps when ``save_plots=True``.

``report.html`` and ``summary.txt``
    Run reports when ``save_report=True`` and the corresponding report formats
    are enabled.

The output-directory details are documented in :doc:`outputs`.

Error Handling
--------------

Pipeline error behavior is controlled by ``PYCSAMT_PIPE.on_step_error`` or by
the CLI ``--on-error`` option.

``"raise"``
    Re-raise the step exception immediately and stop the run.

``"warn"``
    Store the exception in the step result, warn, continue with the previous
    site collection, and mark the final ``PipelineResult`` as not OK.

``"skip"``
    Store the exception and continue silently with the previous site
    collection.

Use ``"raise"`` during debugging and strict production validation.  Use
``"warn"`` for exploratory processing when you want a full report showing
which steps failed.

Runtime Configuration
---------------------

Pipeline runtime defaults live in :data:`pycsamt.pipeline.PYCSAMT_PIPE`.
Configure them globally:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import configure_pipe

   configure_pipe(
       output_root="results",
       on_step_error="warn",
       plot_dpi=200,
       plot_fmt="png",
       show_progress=True,
   )

Or temporarily with a context manager:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import PYCSAMT_PIPE

   with PYCSAMT_PIPE.context(plot_dpi=300, plot_fmt="pdf"):
       result = pipe.run(sites, outdir="results/high_resolution")

Important runtime settings include:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Setting
     - Meaning
   * - ``output_root``
     - Default output root when no explicit run output is provided.
   * - ``processed_subdir``
     - Name of the subdirectory for processed EDI files.
   * - ``plots_subdir``
     - Name of the subdirectory for QC figures.
   * - ``on_step_error``
     - ``"raise"``, ``"warn"``, or ``"skip"``.
   * - ``save_intermediate``
     - Whether to write EDI snapshots after each successful step.
   * - ``show_progress``
     - Whether to print progress while running.
   * - ``plot_dpi`` and ``plot_fmt``
     - Saved figure resolution and format.
   * - ``report_formats``
     - Report types to write, usually ``("html", "txt")``.

PipelineResult
--------------

``Pipeline.run`` returns a :class:`pycsamt.pipeline.PipelineResult`.  Use it as
the programmatic summary of the run:

.. code-block:: python
   :linenos:

   result = pipe.run(sites, outdir="results/basic_qc")

   print(result.ok)
   print(result.n_errors)
   print(result.plots)
   print(result.processed_paths)
   print(result.summary())

``result.sites_in``
    Original site collection passed to ``Pipeline.run``.

``result.sites_out``
    Final site collection after all steps.

``result.step_results``
    Ordered list of step records.

``result.plots``
    All saved plot paths across every step.

``result.processed_paths``
    Written processed EDI files.

``result.ok``
    ``True`` when every step completed without error.

StepResult
----------

Each :class:`pycsamt.pipeline.StepResult` records what happened during one
step:

.. code-block:: python
   :linenos:

   for step_result in result.step_results:
       print(step_result.summary_line())
       if not step_result.ok:
           print(step_result.error)

Useful fields include ``step_idx``, ``step_name``, ``step_code``,
``step_label``, ``params``, ``elapsed_sec``, ``plots``, ``n_sites_in``,
``n_sites_out``, and ``error``.

CLI And Python Equivalence
--------------------------

The CLI and Python API use the same pipeline engine.

This Python call:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis
   from pycsamt.pipeline import Pipeline

   survey = read_edis("data/edis")
   pipe = Pipeline.from_yaml("config/basic_qc.yaml")
   result = pipe.run(survey.to_collection(), outdir="results/basic_qc")

is conceptually equivalent to:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/basic_qc.yaml \
       --out results/basic_qc

Use Python when the pipeline is part of a larger analysis script.  Use the CLI
when the workflow should be easy to repeat from a terminal, automation script,
or processing log.

How Concepts Connect
--------------------

The pipeline documentation is organized around these ideas:

* :doc:`configuration_files` explains how to store pipelines as YAML, JSON, or
  Python files.
* :doc:`presets` explains built-in workflows and when to start from each one.
* :doc:`steps` explains registered operations and extension patterns.
* :doc:`outputs` explains the generated files and reports.
* :doc:`cli_pipe` explains the ``pycsamt pipe`` command group.

In Short
--------

A pyCSAMT pipeline is a reproducible chain of registered steps:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline, Step

   pipe = Pipeline([
       ("notch", Step("NR001")),
       ("band", Step("FREQ001")),
       ("qc", Step("QC001")),
   ])

   result = pipe.run(sites, outdir="results/basic_qc")

The key ideas are simple: registered step codes define what can run, labels
define how a workflow is reported, configs define reproducibility, runtime
settings define output/error behavior, and ``PipelineResult`` records what
happened.
