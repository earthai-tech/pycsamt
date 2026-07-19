.. _pipeline-concepts:

Pipeline Concepts
=================

The pyCSAMT :term:`processing pipeline` is a reproducible engine for
:term:`MT`, :term:`AMT`, and :term:`CSAMT` site collections.  A pipeline is an
ordered sequence of registered operations.  Each operation receives the
current :term:`site collection`, either transforms it or inspects it, records
what happened, and passes the appropriate collection to the next operation.

At the end of a run, pyCSAMT returns a
:class:`pycsamt.pipeline.PipelineResult`.  When an output directory is enabled,
the same run also writes processed :term:`EDI` files, quality-control figures,
reports, and the exact ``pipeline.yaml`` snapshot used for the run.  The
important point is that the data product and the processing recipe travel
together.

Why Pipelines Exist
-------------------

Field processing often begins as a notebook: remove power-line harmonics, drop
duplicate frequencies, trim to an interpretation band, align stations onto a
shared frequency grid, correct :term:`static shift`, inspect tensor rotation,
and make quality-control plots.  That is a natural way to explore, but it can
become fragile when the same decisions must be repeated on another line or
reviewed months later.

Pipelines make the sequence explicit.  They help you:

* repeat the same workflow across multiple survey lines;
* save the processing recipe alongside the outputs;
* inspect which step changed the data;
* run the same workflow from Python and the CLI;
* generate reports and QC figures in a predictable directory tree;
* share processing decisions with collaborators.

Mathematically, a processing pipeline is a composition.  If the loaded survey
is :math:`S_0` and the pipeline has step transforms
:math:`f_1,\ldots,f_n`, then the final collection is

.. math::

   S_n = f_n\left(f_{n-1}\left(\cdots f_1(S_0)\right)\right).

Diagnostic steps fit the same mental model by returning the input collection
unchanged while recording plots, tables, or warnings.  A failed step under a
non-strict error policy also passes forward the last valid collection; that is
useful for diagnosis, but the report must be read as a failed or partial run.

Core Objects
------------

The public pipeline namespace is :mod:`pycsamt.pipeline`:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step
   >>> Pipeline.__name__, Step.__name__
   ('Pipeline', 'Step')

The main objects are:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Object
     - Role
   * - :class:`pycsamt.pipeline.Pipeline`
     - Ordered sequence of ``(label, Step)`` entries.  It can be built from
       Python code, a preset, or a YAML/JSON/Python config file.
   * - :class:`pycsamt.pipeline.Step`
     - Configured wrapper around one registered processing operation.  It
       binds a registry code such as ``"NR001"`` to parameter overrides.
   * - ``StepSpec``
     - :term:`step registry` descriptor: code, registry name, category,
       function path, default parameters, QC plot functions, and whether the
       step returns a modified site collection.
   * - ``Preset``
     - Named, ordered collection of steps for common workflows such as
       ``"basic_qc"``, ``"full_processing"``, or ``"publication_ready"``.
   * - :class:`pycsamt.pipeline.StepResult`
     - One per executed step.  Records timing, input/output station counts,
       saved plot paths, parameters, and any stored exception.
   * - :class:`pycsamt.pipeline.PipelineResult`
     - Returned by ``Pipeline.run``.  Contains the original sites, final
       sites, all step results, output directory, processed file paths, total
       runtime, and status helpers.

The Mental Model
----------------

A pipeline run is left-to-right data flow:

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
collection.  Diagnostic-only steps run checks or plots and pass the input
collection through unchanged.  The pipeline report is the companion record:
it says which operation ran, with which parameters, how long it took, what it
saved, and whether it failed.

Step Registry
-------------

Pipeline steps are not arbitrary strings.  They are registered in the
:term:`step registry`, and each registered step has:

* a short :term:`pipeline step code`, for example ``NR001``;
* a registry name, for example ``notch_powerline``;
* a category, for example ``noise_removal`` or ``frequency``;
* default parameters;
* a transform function;
* optional QC plot functions;
* a ``returns_sites`` flag.

The code and registry name both identify the same operation:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Step
   >>> notch_by_code = Step("NR001", mains_hz=50.0)
   >>> notch_by_name = Step("notch_powerline", mains_hz=50.0)
   >>> notch_by_code.code, notch_by_name.code
   ('NR001', 'NR001')

Use the code form in configuration files and reports because codes are compact
and stable.  Use registry names in exploratory Python when they make intent
easier to read.

Discover available steps from Python:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> print(Pipeline.catalogue("frequency").splitlines()[0])
   Available pipeline steps
   >>> print(Pipeline.step_info("NR001").splitlines()[0])
   NR001  notch_powerline

Discover the same information from the CLI:

.. code-block:: console
   :linenos:

   pycsamt pipe steps
   pycsamt pipe steps --category frequency
   pycsamt pipe steps --info NR001

Configured Steps
----------------

A :class:`pycsamt.pipeline.Step` combines a registry entry with user parameter
overrides.  The registry defaults are merged with your overrides:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Step
   >>> step = Step("NR001", mains_hz=60.0)
   >>> step.code
   'NR001'
   >>> step.params["mains_hz"]
   60.0

If the registry default for ``NR001`` includes ``n_harm`` and ``tol_hz``, the
configured step still carries those defaults at run time.  You only need to
provide the values that should change for the workflow.  This merge is why a
step should be inspected before editing its ``params``: a short configuration
may still imply several default behaviours.

Pipeline Structure
------------------

A pipeline stores steps as ``(label, Step)`` tuples.  The label names this
occurrence in this workflow.  It appears in printed summaries, output
subdirectories, reports, and CLI slicing options.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step
   >>> pipe = Pipeline(
   ...     [
   ...         ("notch", Step("NR001", mains_hz=50.0)),
   ...         ("select_band", Step("FREQ001", band_hz=(0.001, 10000.0))),
   ...         ("align_grid", Step("FREQ004")),
   ...         ("qc_snapshot", Step("QC001")),
   ...     ],
   ...     name="first_qc",
   ... )
   >>> print(pipe)
   Pipeline  'first_qc'  -  4 steps
     ( 1) notch        [NR001]    Power-line Harmonic Notch  mains_hz=50.0  n_harm=30  tol_hz=0.08
     ( 2) select_band  [FREQ001]  Frequency Band Select      band_hz=(0.001, 10000.0)
     ( 3) align_grid   [FREQ004]  Frequency Grid Alignment
     ( 4) qc_snapshot  [QC001]    QC Quick-Look Snapshot

Labels should be short, stable, and meaningful.  Prefer
``select_amt_band`` or ``correct_ss`` over vague labels such as ``step1``.
Changing a label changes report names and slicing handles, even when the
underlying step code is unchanged.

Building A Pipeline
-------------------

There are four common ways to build a pipeline.

Build directly in Python:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step
   >>> pipe = Pipeline([
   ...     ("notch", Step("NR001")),
   ...     ("drop_duplicates", Step("FREQ002")),
   ...     ("select_band", Step("FREQ001")),
   ...     ("qc_snapshot", Step("QC001")),
   ... ])
   >>> type(pipe).__name__
   'Pipeline'

Build from a preset:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> pipe = Pipeline.from_preset("basic_qc")
   >>> pipe.name
   'basic_qc'

Build from a config file:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> # pipe = Pipeline.from_yaml("config/basic_qc.yaml")
   >>> # pipe = Pipeline.from_json("config/basic_qc.json")
   >>> # pipe = Pipeline.from_py("config/basic_qc.py")

Build from the CLI:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc
   pycsamt pipe run data/edis --config config/basic_qc.yaml
   pycsamt pipe run data/edis --steps NR001,FREQ002,FREQ001,QC001

Configuration files are documented in :doc:`configuration_files`.

Presets
-------

:term:`Pipeline preset`\ s are named pipelines for common processing
intentions.  They are useful when you want a known baseline without writing
every step manually.

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

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> pipe = Pipeline.from_preset("publication_ready")
   >>> pipe.name
   'publication_ready'

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

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step
   >>> pipe = Pipeline.from_preset("full_processing")
   >>> pipe.remove("mask_skew")
   >>> pipe.append("final_qc", Step("QC001"))
   >>> pipe.replace("notch", Step("NR001", mains_hz=60.0))

During ``Pipeline.run``, the step list is protected from mutation.  This
prevents accidental changes while step results and reports are being produced.
Think of the run as an immutable transaction: once execution begins, the
pipeline must be the same object that later appears in ``pipeline.yaml``.

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

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.emtools import ensure_sites
   >>> from pycsamt.pipeline import Pipeline
   >>> sites = ensure_sites("data/3edis", recursive=True, verbose=0)
   >>> pipe = Pipeline.from_preset("basic_qc")
   >>> # result = pipe.run(
   >>> #     sites,
   >>> #     outdir="results/basic_qc",
   >>> #     save_plots=True,
   >>> #     save_edis=True,
   >>> #     save_report=True,
   >>> # )
   >>> len(sites)
   3

The final run is commented here because it writes a full output directory.
Use the same call in a project script when you want processed files, plots,
and reports to be created.

Output Resolution
-----------------

The output directory is resolved in this order:

1. explicit ``outdir`` passed to ``Pipeline.run``;
2. ``output_dir`` stored on a pipeline loaded from a config file;
3. global ``PYCSAMT_PIPE.output_root``.

Passing ``outdir=None`` is an explicit opt-out: the pipeline runs in memory
and writes no output files.

.. code-block:: pycon
   :linenos:

   >>> # Write to the config/default output directory.
   >>> # result = pipe.run(sites)
   >>> # Override output directory for this run.
   >>> # result = pipe.run(sites, outdir="results/experiment_01")
   >>> # In-memory run: no files are written.
   >>> # result = pipe.run(sites, outdir=None)

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
which steps failed.  A run that continued after a failed step is diagnostic
evidence, not a final processing product.

Runtime Configuration
---------------------

Pipeline runtime defaults live in :data:`pycsamt.pipeline.PYCSAMT_PIPE`.
Configure them globally:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import configure_pipe
   >>> configure_pipe(
   ...     output_root="results",
   ...     on_step_error="warn",
   ...     plot_dpi=200,
   ...     plot_fmt="png",
   ...     show_progress=True,
   ... )

Or temporarily with a context manager:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import PYCSAMT_PIPE
   >>> # with PYCSAMT_PIPE.context(plot_dpi=300, plot_fmt="pdf"):
   >>> #     result = pipe.run(sites, outdir="results/high_resolution")

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

.. code-block:: pycon
   :linenos:

   >>> # result = pipe.run(sites, outdir="results/basic_qc")
   >>> # result.ok
   >>> # result.n_errors
   >>> # result.plots
   >>> # result.processed_paths
   >>> # print(result.summary())

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

.. code-block:: pycon
   :linenos:

   >>> # for step_result in result.step_results:
   >>> #     print(step_result.summary_line())
   >>> #     if not step_result.ok:
   >>> #         print(step_result.error)

Useful fields include ``step_idx``, ``step_name``, ``step_code``,
``step_label``, ``params``, ``elapsed_sec``, ``plots``, ``n_sites_in``,
``n_sites_out``, and ``error``.  These fields are the fine-grained audit trail
behind a pipeline run: they explain why the final result is OK, warning-only,
or failed.

CLI And Python Equivalence
--------------------------

The CLI and Python API use the same pipeline engine.

This Python call:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.emtools import ensure_sites
   >>> from pycsamt.pipeline import Pipeline
   >>> sites = ensure_sites("data/3edis", recursive=True, verbose=0)
   >>> pipe = Pipeline.from_yaml("config/basic_qc.yaml")  # doctest: +SKIP
   >>> result = pipe.run(sites, outdir="results/basic_qc")  # doctest: +SKIP

is conceptually equivalent to:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
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

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step
   >>> pipe = Pipeline([
   ...     ("notch", Step("NR001")),
   ...     ("band", Step("FREQ001")),
   ...     ("qc", Step("QC001")),
   ... ])
   >>> # result = pipe.run(sites, outdir="results/basic_qc")
   >>> type(pipe).__name__
   'Pipeline'

The key ideas are simple: registered step codes define what can run, labels
define how a workflow is reported, configs define reproducibility, runtime
settings define output and error behavior, and ``PipelineResult`` records what
happened.
