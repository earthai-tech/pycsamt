.. _pipeline-cli:

Pipeline CLI
============

The ``pycsamt pipe`` command group is the command-line interface to the
pyCSAMT processing pipeline.  It uses the same pipeline engine as the Python
API, so a workflow can move between a terminal, a notebook, and a project
configuration file without changing meaning.

Use ``pycsamt pipe`` when you want to:

* list available pipeline presets;
* inspect registered processing step codes;
* scaffold a YAML, JSON, or Python pipeline config;
* preview a config before running it;
* run a preset, config file, or ad-hoc step list against an EDI directory;
* save processed EDI files, plots, reports, and a reproducible
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

The examples below assume the survey data live in ``data/edis`` and outputs
should be written under ``results/``.

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

* reads the EDI directory;
* resolves the ``basic_qc`` preset;
* runs each step in order;
* writes processed EDI files to ``results/basic_qc/processed``;
* writes QC figures under ``results/basic_qc/plots``;
* writes ``report.html`` and ``summary.txt``;
* saves ``results/basic_qc/pipeline.yaml`` for reproducibility.

Use ``--dry-run`` first when you want to confirm the pipeline and site count
without writing files:

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

Machine-readable formats:

.. code-block:: console
   :linenos:

   pycsamt pipe steps --format json
   pycsamt pipe steps --category frequency --format csv
   pycsamt pipe steps --info FREQ001 --format json

Use ``pycsamt pipe steps --info CODE`` before editing ``params`` in a config
file.  It shows the registry defaults, transform function, QC plots, and
whether the step returns a modified site collection.

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

Run Controls
------------

``pycsamt pipe run`` supports slicing options that make debugging easier.

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

Exit Status
-----------

``pycsamt pipe run`` exits with status ``0`` when the final
``PipelineResult.ok`` is true.  If one or more steps failed and the command
continued under ``--on-error warn`` or ``--on-error skip``, the command prints
the summary and exits nonzero.

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
* :doc:`outputs` explains the generated output directory.
