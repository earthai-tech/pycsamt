.. _pipeline-cli:

Pipeline CLI
============

The ``pycsamt pipe`` command group is the command-line interface to the
pyCSAMT :term:`processing pipeline`.  It uses the same pipeline engine as the
Python API, so a workflow can move between a terminal, a notebook, and a
project configuration file without changing meaning.  The CLI is therefore not
a lighter or separate implementation: it is a reproducible shell around the
same :class:`pycsamt.pipeline.Pipeline`, :class:`pycsamt.pipeline.Step`, and
``PipelineResult`` objects used from Python.

Use the CLI when the work needs to leave a clear command trail.  A notebook is
excellent for exploration, but a terminal command is easier to rerun in
continuous integration, paste into a field-processing log, or compare across
survey revisions.  The important idea is that the command, resolved step list,
input survey, output directory, generated ``pipeline.yaml``, and saved reports
together form the audit record for the processing decision.

Use ``pycsamt pipe`` when you want to:

* list available :term:`pipeline preset`\ s;
* inspect registered :term:`pipeline step code`\ s;
* scaffold a YAML, JSON, or Python pipeline config;
* preview a config before running it;
* run a preset, config file, or ad-hoc step list against an EDI directory;
* save processed :term:`EDI` files, plots, reports, and a reproducible
  ``pipeline.yaml`` snapshot.

Command Overview
----------------

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Command
     - Purpose
   * - ``pycsamt pipe presets``
     - List built-in presets and optionally expand one preset into its step
       sequence.
   * - ``pycsamt pipe steps``
     - List registered step codes, filter by category, or show detailed
       information for one step.
   * - ``pycsamt pipe init``
     - Generate a ready-to-edit YAML, JSON, or Python config file.
   * - ``pycsamt pipe show``
     - Pretty-print a config file or preset before running it.
   * - ``pycsamt pipe run``
     - Run a pipeline against MT/AMT EDI data.
   * - ``pycsamt pipe plugins``
     - Discover and list third-party :term:`pipeline plugin` steps.
   * - ``pycsamt pipe history``
     - List runs previously logged with ``pycsamt pipe run --history``.

The ``pipe`` group itself also accepts ``--with-plugins``, which discovers
plugin steps before any subcommand runs, so a plugin step code can be passed
to ``--steps`` in the same command.  It is off by default because the
underlying scan can take several seconds in a large environment; see
:doc:`extending` for the full plugin mechanism and when this flag is needed.

``pycsamt pipe run`` also accepts ``--cache``/``--cache-dir``, which caches
each step's output so a rerun of the identical command replays already-
completed steps instead of recomputing them -- also how an interrupted run
resumes.  Off by default; see :doc:`caching` for the full mechanism and its
caveats.

``--live`` renders a live-updating per-step status table instead of a
static progress bar; ``--history``/``--history-file`` logs the run for
``pycsamt pipe history`` to list later.  Both off by default; see
:doc:`observability` for the full mechanism, including the ``on_step``
Python hook these CLI flags build on.

``--dashboard`` adds a richer, branded ``dashboard.html`` report alongside
the default ``report.html``/``summary.txt`` -- KPI stat tiles and inline-SVG
charts built from the same per-step data.  Off by default; see
:doc:`outputs` for what it contains.

Start With Help
---------------

The CLI help is the quickest reference for exact options in your installed
version:

.. code-block:: console
   :linenos:

   pycsamt pipe --help
   pycsamt pipe run --help
   pycsamt pipe init --help
   pycsamt pipe steps --help
   pycsamt pipe presets --help
   pycsamt pipe show --help
   pycsamt pipe plugins --help
   pycsamt pipe history --help

The examples below assume the survey data live in ``data/edis`` and outputs
should be written under ``results/``.

If you are running directly from a source checkout where the console script is
not installed on ``PATH``, the equivalent entry point is:

.. code-block:: console
   :linenos:

   python -m pycsamt.cli pipe --help

The documentation uses ``pycsamt pipe`` because that is the installed command
users normally run.

Typical First Run
-----------------

For a first survey, run the ``basic_qc`` preset:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/basic_qc \
       --on-error warn \
       --dpi 200 \
       --plot-fmt png \
       -v

This command:

* reads the :term:`EDI` directory;
* resolves the ``basic_qc`` preset;
* runs each step in order;
* writes processed EDI files to ``results/basic_qc/processed``;
* writes QC figures under ``results/basic_qc/plots``;
* writes ``report.html`` and ``summary.txt``;
* saves ``results/basic_qc/pipeline.yaml`` for reproducibility.

Under the hood the run can be thought of as a composition of site-collection
transforms.  If :math:`S_0` is the loaded survey and the resolved pipeline has
step functions :math:`f_1,\ldots,f_n`, then a successful run produces

.. math::

   S_n = f_n(\cdots f_2(f_1(S_0))\cdots).

Each step also produces a ``StepResult`` containing status, runtime, messages,
errors, and output artefacts.  The final ``PipelineResult.ok`` flag is true
only when the executed step results are acceptable under the selected error
policy.

Use ``--dry-run`` first when you want to confirm the pipeline and site count
without writing files.  A dry run resolves the survey, config or preset,
slicing flags, and output target, then stops before mutating data or creating
artefacts:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/basic_qc \
       --dry-run

Input Survey Resolution
-----------------------

``pycsamt pipe run`` can receive survey data in three ways.  The priority is:

1. positional ``EDI_DIR`` argument;
2. ``--survey EDI_DIR`` option;
3. active survey context, if one was configured through the survey CLI.

This priority prevents accidental reuse of stale state.  A positional path is
the most explicit statement of intent, so it wins over ``--survey`` and over a
stored :term:`active survey context`.  The active context is convenient for
interactive command sequences, but it should be treated as session state, not
as a substitute for a documented input path in a reproducible script.

Examples:

.. code-block:: console
   :linenos:

   # Positional source.
   pycsamt pipe run data/edis --preset basic_qc

   # Equivalent explicit survey option.
   pycsamt pipe run --survey data/edis --preset basic_qc

   # Use the active survey context.
   pycsamt pipe run --preset basic_qc

Use ``--fresh`` when the survey has been edited on disk and should be
re-parsed rather than resolved from a cached survey context:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc --fresh

Pipeline Definition Priority
----------------------------

``pycsamt pipe run`` can build the pipeline from a config file, a preset, or
an ad-hoc comma-separated step list.  The priority is:

1. ``--config FILE``;
2. ``--preset NAME``;
3. ``--steps CODE,CODE,...``.

If ``--config`` is present, the config file is the source of truth and
``--preset`` or ``--steps`` are ignored.

The priority is deliberately conservative.  A configuration file can carry a
pipeline name, output directory, labels, and per-step ``params``; allowing an
extra ``--preset`` or ``--steps`` value to partially override it would make the
run harder to reconstruct.  Use CLI overrides for execution controls such as
``--out``, ``--format``, ``--dry-run``, and slicing, then keep the processing
definition itself in one place.

Run from a config file:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/basic_qc.yaml \
       --out results/basic_qc

Run from a preset:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset full_processing \
       --out results/full_processing

Run an ad-hoc step list:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --steps NR001,FREQ002,FREQ001,FREQ004,QC001 \
       --out results/ad_hoc_qc

Use ad-hoc step lists for quick experiments.  Use config files for workflows
that should be repeated or reviewed.

List Presets
------------

List all built-in presets:

.. code-block:: console
   :linenos:

   pycsamt pipe presets

Expand one preset to see its steps:

.. code-block:: console
   :linenos:

   pycsamt pipe presets --expand full_processing

Get machine-readable output:

.. code-block:: console
   :linenos:

   pycsamt pipe presets --format json
   pycsamt pipe presets --expand basic_qc --format json
   pycsamt pipe presets --format csv

A shortened JSON listing looks like this:

.. code-block:: json
   :linenos:

   [
     {
       "name": "basic_qc",
       "description": "Minimal denoising + frequency cleanup.  Good for quick-look inspection.",
       "n_steps": 5,
       "codes": ["NR001", "FREQ002", "FREQ001", "FREQ004", "QC001"]
     },
     {
       "name": "full_processing",
       "description": "Standard end-to-end workflow: noise -> frequency -> skew gate -> static-shift -> strike rotation.",
       "n_steps": 8,
       "codes": ["NR001", "FREQ002", "FREQ001", "FREQ004", "SK001", "TZ001", "SS001", "QC001"]
     }
   ]

The full command also lists ``noise_reduction``, ``tensor_analysis``,
``dimensionality_filter``, ``publication_ready``, and ``stratagem_mt`` in the
current package.

Common preset choices:

``basic_qc``
    First-pass quality control and frequency cleanup.

``noise_reduction``
    Noise-removal chain for high-EMI environments.

``full_processing``
    Standard processing sequence: noise, frequency cleanup, skew gate,
    static-shift correction, and strike rotation.

``publication_ready``
    Longer chain intended for polished processing outputs.

See :doc:`presets` for the full preset guide.

List And Inspect Steps
----------------------

List all registered step codes grouped by category:

.. code-block:: console
   :linenos:

   pycsamt pipe steps

Filter by category:

.. code-block:: console
   :linenos:

   pycsamt pipe steps --category frequency
   pycsamt pipe steps --category noise_removal
   pycsamt pipe steps --category static_shift

Show detailed information for one step:

.. code-block:: console
   :linenos:

   pycsamt pipe steps --info NR001
   pycsamt pipe steps --info notch_powerline

Print only step codes:

.. code-block:: console
   :linenos:

   pycsamt pipe steps --codes-only
   pycsamt pipe steps --category static_shift --codes-only

Captured output begins with the frequency-editing steps and then continues
through noise removal, static shift, tensor/strike tools, source checks, and
QC:

.. code-block:: text
   :linenos:

   FREQ001
   FREQ002
   FREQ003
   FREQ004
   FREQ005
   FREQ006
   FREQ007
   FREQ008
   FREQ009
   NR001
   NR002
   ...
   QC001
   QC002
   QC003
   QC004

Machine-readable formats:

.. code-block:: console
   :linenos:

   pycsamt pipe steps --format json
   pycsamt pipe steps --category frequency --format csv
   pycsamt pipe steps --info FREQ001 --format json

Use ``pycsamt pipe steps --info CODE`` before editing ``params`` in a config
file.  It shows the registry defaults, transform function, QC plots, and
whether the step returns a modified site collection.  This matters because a
pipeline is not only a list of labels.  Each step has a declared transform,
parameter defaults, category, and optional review plots; changing ``params``
changes the scientific operation applied to the survey.

Generate Config Files
---------------------

Use ``pycsamt pipe init`` to create a valid starter config.

Default YAML scaffold:

.. code-block:: console
   :linenos:

   pycsamt pipe init

Scaffold from a preset:

.. code-block:: console
   :linenos:

   pycsamt pipe init \
       --preset basic_qc \
       --name line22_basic_qc \
       --outdir results/line22_basic_qc \
       --output config/line22_basic_qc.yaml

Generate Python or JSON:

.. code-block:: console
   :linenos:

   pycsamt pipe init --format py \
       --preset basic_qc \
       --output config/line22_basic_qc.py

   pycsamt pipe init --format json \
       --preset basic_qc \
       --output config/line22_basic_qc.json

Print the generated config instead of writing it:

.. code-block:: console
   :linenos:

   pycsamt pipe init --preset full_processing --print

When ``--output`` points to a directory, the filename is derived from
``--name`` and ``--format``.  When ``--output`` is omitted, the file is written
as ``<name>.<format>`` in the current directory.

See :doc:`configuration_files` for the config schema and editing guidance.

Prefer generated starter files over hand-written configs for production
workflows.  The scaffold captures the schema expected by the installed
version, while still leaving every step parameter visible for review.  A
typical review question is not "can the command run?", but "can someone tell
which operation happened at each station and frequency?"

Show A Pipeline Before Running
------------------------------

Use ``pycsamt pipe show`` to inspect a config or preset without running data
through it.

Show a config file:

.. code-block:: console
   :linenos:

   pycsamt pipe show config/line22_basic_qc.yaml

Show a preset:

.. code-block:: console
   :linenos:

   pycsamt pipe show --preset publication_ready

Preview a sliced pipeline:

.. code-block:: console
   :linenos:

   pycsamt pipe show config/line22_basic_qc.yaml --n-steps 3
   pycsamt pipe show config/line22_basic_qc.yaml --from-step select_band
   pycsamt pipe show config/line22_basic_qc.yaml --until-step QC001

Use JSON or CSV when another tool needs the resolved step list:

.. code-block:: console
   :linenos:

   pycsamt pipe show --preset full_processing --format json
   pycsamt pipe show config/line22_basic_qc.yaml --format csv

``show`` is also the safest place to verify pipeline slicing.  Slicing changes
which transforms are actually executed, so it should be visible before any
processed EDI files or plots are written.

Run Controls
------------

``pycsamt pipe run`` supports slicing options that make debugging easier.

Slicing is evaluated after the pipeline has been resolved from ``--config``,
``--preset``, or ``--steps``.  If the full ordered sequence is
:math:`(s_1,\ldots,s_n)`, then ``--from-step`` and ``--until-step`` choose a
contiguous subsequence and ``--n-steps`` truncates that subsequence.  The
resulting run is still a valid pipeline run, but its ``pipeline.yaml`` and
summary should be read as a partial execution, not as evidence that the full
workflow passed.

Run only the first N steps:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/line22_basic_qc.yaml \
       --n-steps 3 \
       --out results/debug_first3

Start from a label, registry name, or code:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/line22_basic_qc.yaml \
       --from-step select_band \
       --out results/debug_from_band

Stop after a label, registry name, or code:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/line22_basic_qc.yaml \
       --until-step QC001 \
       --out results/debug_until_qc

Combine slicing options:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/full_processing.yaml \
       --from-step FREQ001 \
       --until-step SS001 \
       --out results/frequency_to_static_shift

Use ``pycsamt pipe show`` with the same slicing flags before running if you
want to verify the resolved step list.

Error Policy
------------

Use ``--on-error`` to control what happens when a step raises an exception.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Value
     - Behavior
   * - ``raise``
     - Stop immediately and report the exception.
   * - ``warn``
     - Warn, store the error in the step result, continue with the previous
       site collection, and exit nonzero at the end if any step failed.
   * - ``skip``
     - Store the error and continue silently with the previous site
       collection.  The final command still exits nonzero when the result is
       not OK.

Examples:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc --on-error raise
   pycsamt pipe run data/edis --preset basic_qc --on-error warn
   pycsamt pipe run data/edis --preset basic_qc --on-error skip

Use ``raise`` while debugging.  Use ``warn`` for exploratory runs where the
report should show every step that failed.

The error policy controls continuation, not scientific acceptance.  With
``warn`` or ``skip``, a later step may run on the last valid site collection
instead of on the failed step's intended output.  That is useful for diagnosis
because the report can show several failures in one pass, but the resulting
processed files should not be treated as final deliverables until the failed
steps have been resolved.

Output Controls
---------------

Choose the run directory:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/line22_basic_qc.yaml \
       --out results/line22_basic_qc

Skip specific output families:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc --no-plots
   pycsamt pipe run data/edis --preset basic_qc --no-edi
   pycsamt pipe run data/edis --preset basic_qc --no-report

Control saved figure format:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset publication_ready \
       --out results/publication_ready_pdf \
       --dpi 300 \
       --plot-fmt pdf

Supported plot formats are ``png``, ``pdf``, and ``svg``.

For reproducible reports, choose ``png`` for quick review and lightweight
HTML, ``pdf`` for publication-oriented vector handoff, and ``svg`` when plots
must remain inspectable in text-based review.  The numeric plot content comes
from the step result; ``--dpi`` and ``--plot-fmt`` only control how the figure
is written.

Output Summary Formats
----------------------

The ``--format`` option controls the terminal output summary, not the saved
report files.

Text summary:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc --format text

JSON summary:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc --format json

CSV summary:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc --format csv

Use JSON or CSV in automation when another script needs step status,
runtime, output directory, and error information.

The terminal summary format is intentionally separate from saved artefacts.
Changing ``--format`` from ``text`` to ``json`` changes what is printed to
stdout, but it does not remove ``summary.txt``, ``report.html``, figures, or
the reproduced ``pipeline.yaml`` when those output families are enabled.

Verbose And Color Options
-------------------------

Use ``-v`` or ``--verbose`` to show more progress information:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc -v
   pycsamt pipe run data/edis --preset basic_qc -vv

Use ``--no-color`` for plain logs in terminals or CI environments that do not
handle ANSI color:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc --no-color

The ``--jobs`` option is accepted as a shared CLI option.  Current pipeline
steps run through the pipeline engine in order; treat ``--jobs`` as a
forward-compatible option unless a specific step documents parallel behavior.

Sequential execution is important for correctness.  Many processing steps
consume the edited site collection produced by the previous step, so the
pipeline is ordered dataflow rather than a bag of independent tasks.  Parallel
work is only safe inside a step that explicitly documents how it partitions
stations, frequencies, or plots.

Exit Status
-----------

``pycsamt pipe run`` exits with status ``0`` when the final
``PipelineResult.ok`` is true.  If one or more steps failed and the command
continued under ``--on-error warn`` or ``--on-error skip``, the command prints
the summary and exits nonzero.

In automation, treat the exit status as the hard gate and the JSON/CSV summary
as the explanation.  A successful process exit means the pipeline accepted the
run under the selected policy; a nonzero exit means a script should stop,
archive the summary, and surface the failing step names.

This makes the command useful in automation:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/line22_basic_qc.yaml \
       --out results/line22_basic_qc \
       --format json

Troubleshooting
---------------

No pipeline specified
    Provide one of ``--config``, ``--preset``, or ``--steps``.  Run
    ``pycsamt pipe presets`` and ``pycsamt pipe steps`` to discover valid
    choices.

Unknown preset
    Run ``pycsamt pipe presets`` to list names.  Preset names are exact, for
    example ``basic_qc`` or ``publication_ready``.

Unknown step code
    Run ``pycsamt pipe steps`` or ``pycsamt pipe steps --info CODE``.  Step
    identifiers can be codes such as ``NR001`` or names such as
    ``notch_powerline``.

Unsupported config format
    ``--config`` accepts ``.yaml``, ``.yml``, ``.json``, and ``.py`` files.

Dry-run output directory looks generic
    Pass ``--out`` explicitly when the dry-run display must show the exact
    target directory.  During a real run, the pipeline resolves the output
    directory from the CLI override, then the config ``output_dir``, then the
    global pipeline default.

Command exits nonzero after warnings
    At least one step failed.  Inspect the text summary, JSON output, or
    saved report.  Re-run with ``--on-error raise`` to stop at the first
    failing step.

Generated plots are missing
    Check that ``--no-plots`` was not used, that an output directory is
    enabled, and that the selected steps have QC plot functions.

Processed EDI files are missing
    Check that ``--no-edi`` was not used and inspect warnings from the export
    step.  Some malformed or partial site collections may still run through
    diagnostic steps but fail during final EDI export.

Recommended Workflow
--------------------

For a real survey, a robust command-line workflow looks like this:

.. code-block:: console
   :linenos:

   pycsamt pipe presets
   pycsamt pipe presets --expand basic_qc
   pycsamt pipe init --preset basic_qc \
       --name line22_basic_qc \
       --outdir results/line22_basic_qc \
       --output config/line22_basic_qc.yaml
   pycsamt pipe show config/line22_basic_qc.yaml
   pycsamt pipe run data/edis \
       --config config/line22_basic_qc.yaml \
       --dry-run
   pycsamt pipe run data/edis \
       --config config/line22_basic_qc.yaml \
       --out results/line22_basic_qc \
       --on-error warn \
       -v

After the run, inspect ``results/line22_basic_qc/report.html``,
``results/line22_basic_qc/summary.txt``, ``results/line22_basic_qc/plots``,
and ``results/line22_basic_qc/pipeline.yaml``.

Related Pages
-------------

* :doc:`concepts` explains the pipeline object model and run lifecycle.
* :doc:`configuration_files` explains YAML, JSON, and Python config files.
* :doc:`presets` explains built-in workflows and how to choose one.
* :doc:`steps` explains registered operations and extension patterns.
* :doc:`extending` explains ``pycsamt pipe plugins`` and ``--with-plugins``
  in depth.
* :doc:`caching` explains ``--cache``/``--cache-dir`` and how a rerun
  resumes an interrupted pipeline.
* :doc:`observability` explains ``--live``, ``--history``, and the
  ``on_step`` Python hook.
* :doc:`outputs` explains the generated output directory.
