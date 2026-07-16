CLI Overview
============

The ``pycsamt`` command is the terminal interface to pyCSAMT's survey
loading, inspection, conversion, processing, modelling, inversion,
mapping, and interpretation tools. It is built on Click and delegates to
the same package APIs used by Python workflows, so a command-line run can
usually be translated into Python code when you need finer control.

Use this page to understand the root command, shared options, path
resolution rules, output conventions, and the way command families fit
together. The command-specific pages document every family in detail.

Top-Level Command
-----------------

Start with the root command:

.. code-block:: console

   pycsamt
   pycsamt --help
   pycsamt --version

When invoked without a subcommand, ``pycsamt`` prints grouped help. If
``rich`` is installed, the help is rendered as a colorized panel;
otherwise Click's plain text help is used.

The command form is:

.. code-block:: console

   pycsamt [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS] [ARGUMENTS]

Examples:

.. code-block:: console

   pycsamt edi validate data/edis/
   pycsamt site info data/edis/
   pycsamt pipe run --preset basic_qc --survey data/edis/ --out results/qc/
   pycsamt invert build data/edis/ --solver occam2d --workdir runs/occam/

Global Options
--------------

Global options must appear before the command family.

.. list-table::
   :header-rows: 1
   :widths: 24 28 48

   * - Option
     - Example
     - Meaning
   * - ``--version``, ``-V``
     - ``pycsamt -V``
     - Print the installed pyCSAMT version and exit.
   * - ``--verbose``, ``-v``
     - ``pycsamt -v site info data/edis/``
     - Increase verbosity. Use ``-v`` for informational output and
       ``-vv`` for debug-level output where commands support it.
   * - ``--no-color``
     - ``pycsamt --no-color pipe presets``
     - Disable ANSI color output. Useful for logs, CI, and terminals
       that do not handle color escape sequences.

The root command reads ``~/.pycsamt.toml`` before dispatching to a
subcommand. Command-line flags win over values loaded from the
configuration file.

Configuration Loading
---------------------

At startup, the CLI looks for:

.. code-block:: text

   ~/.pycsamt.toml

If the file exists, pyCSAMT loads the configured sections and applies
them to the runtime CLI settings. The ``config`` command group is the
preferred way to inspect and edit this file:

.. code-block:: console

   pycsamt config show
   pycsamt config list
   pycsamt config set plot.dpi 300
   pycsamt config reset plot --yes

Use persistent configuration for preferences that should follow you
across sessions, such as plotting style, output format, default figure
DPI, pipeline behavior, or agent provider preferences. Use command
options for one run only.

Command Families
----------------

The CLI is organized by workflow domain.

.. list-table::
   :header-rows: 1
   :widths: 18 28 54

   * - Command
     - Typical first command
     - Purpose
   * - ``config``
     - ``pycsamt config show``
     - Inspect and persist CLI preferences in ``~/.pycsamt.toml``.
   * - ``survey``
     - ``pycsamt survey set data/edis/``
     - Register, show, rebuild, and clear the active survey context.
   * - ``edi``
     - ``pycsamt edi validate data/edis/``
     - Validate, inspect, rotate, select, and profile EDI files.
   * - ``avg``
     - ``pycsamt avg info raw/L18.avg``
     - Inspect, validate, correct, and export Zonge AVG files.
   * - ``jones``
     - ``pycsamt jones validate raw_j/``
     - Inspect, validate, select, and examine Jones ``.j`` files.
   * - ``transform``
     - ``pycsamt transform avg raw_avg/ --output-dir edis/``
     - Convert spectra, AVG, and Jones products into impedance EDI-style
       outputs.
   * - ``site``
     - ``pycsamt site info data/edis/``
     - Load station collections, summarize usable data, compute
       attributes, select stations, edit metadata, and export products.
   * - ``map``
     - ``pycsamt map stations data/edis/``
     - List station coordinates and create station map figures.
   * - ``pipe``
     - ``pycsamt pipe presets``
     - Discover processing steps, initialize workflows, preview
       pipelines, and run reproducible processing chains.
   * - ``forward``
     - ``pycsamt forward model geology``
     - Run synthetic forward responses and generate modelling datasets.
   * - ``invert``
     - ``pycsamt invert build data/edis/``
     - Build inversion workdirs, launch supported solvers, inspect
       status/results, and plot inversion products.
   * - ``interp``
     - ``pycsamt interp rocks --rho 250``
     - Classify inversion models, query resistivity labels, and export
       interpreted products.
   * - ``tdem``
     - ``pycsamt tdem info raw_tdem/``
     - Inspect, convert, and plot transient electromagnetic survey data.
   * - ``info``
     - ``pycsamt info data/edis/``
     - Legacy/general quick EDI information command.
   * - ``convert``
     - ``pycsamt convert raw_j/ --out converted/``
     - Legacy/general conversion entry point. Prefer ``transform`` for
       the newer grouped conversion workflows.

The detailed reference pages focus on the main workflow groups:
``config``, ``edi``, ``avg``, ``jones``, ``site``, ``transform``,
``forward``, ``invert``, ``map``, ``interp``, ``pipe``, ``tdem``, and
``survey``.

Active Survey Context
---------------------

Many commands can work from an explicit path:

.. code-block:: console

   pycsamt site info data/AMT/WILLY_DATA/L18PLT
   pycsamt map plot data/AMT/WILLY_DATA/L18PLT --output figures/l18.png
   pycsamt invert build data/AMT/WILLY_DATA/L18PLT --workdir runs/l18_occam/

For exploratory work, register an active survey once:

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L18PLT
   pycsamt survey show
   pycsamt site info
   pycsamt map stations
   pycsamt pipe run --preset basic_qc --out results/l18_qc/

The active survey stores a path and a cached ``Sites`` object so repeated
commands do not need to reload the same EDI directory every time.

Use explicit paths in scripts, notebooks converted to shell commands,
CI jobs, and published methods sections. Use the active survey for
interactive terminal sessions where convenience matters more than
self-contained command history.

Survey Cache
------------

``pycsamt survey set`` pre-builds a cache for the active survey. This is
useful because loading and validating a full station collection can be
more expensive than reading a simple configuration value.

Common cache commands:

.. code-block:: console

   pycsamt survey show
   pycsamt survey rebuild
   pycsamt survey rebuild --force
   pycsamt survey cache list
   pycsamt survey cache purge
   pycsamt survey clear

Use ``rebuild`` after editing EDI files, adding stations, replacing
metadata, or when ``survey show`` reports that the cache is stale.

Input Paths
-----------

Commands generally accept one of these path types:

* a single station file, such as ``S01.edi`` or ``S01.j``;
* a directory of station files;
* a work directory created by an inversion command;
* an output directory for converted, corrected, or exported products.

Path meaning depends on the command. For example:

.. code-block:: console

   pycsamt edi info data/edis/
   pycsamt avg info raw/L18.avg
   pycsamt invert status runs/line01_occam/
   pycsamt interp export runs/line01_occam/ --output-dir exports/

As a rule, keep raw field data separate from generated products:

.. code-block:: text

   raw_edis/
   clean_edis/
   converted_edis/
   runs/
   results/
   figures/
   exports/

This makes it clear which files came from the field and which files were
created by pyCSAMT.

Output Formats
--------------

Many inspection commands support ``--format``:

.. code-block:: console

   pycsamt edi validate data/edis/ --format json
   pycsamt edi stations data/edis/ --format csv
   pycsamt pipe show --preset basic_qc --format json
   pycsamt invert status runs/line01_occam/ --format json

Use the default text output for terminal reading. Use JSON when another
program will parse the result. Use CSV when the output is tabular and is
headed toward a spreadsheet, GIS tool, or report table.

When a command supports ``--output`` or ``--output-dir``, prefer writing
directly to a named file or folder instead of redirecting terminal text:

.. code-block:: console

   pycsamt map stations data/edis/ --format csv --output station_coords.csv
   pycsamt map plot data/edis/ --output figures/stations.png
   pycsamt interp export runs/line01_occam/ --format csv --output-dir exports/

Dry Runs
--------

Commands that can modify files or launch larger workflows often provide
``--dry-run``. Use it to confirm inputs, selections, and planned outputs
before writing anything:

.. code-block:: console

   pycsamt edi select raw_edis/ --min-nfreq 20 --dry-run
   pycsamt edi rotate raw_edis/ --angle 30 --dry-run
   pycsamt avg correct raw/L18.avg --method static-shift --dry-run
   pycsamt pipe run --config workflow.yaml --dry-run
   pycsamt tdem convert raw_tdem/ --output-dir converted_tdem/ --dry-run

Dry-run output should be treated as a planning report. Once the plan
looks right, rerun the same command with the output option enabled and
without ``--dry-run``.

Typical Workflows
-----------------

New EDI survey inventory:

.. code-block:: console

   pycsamt edi validate raw_edis/
   pycsamt edi info raw_edis/ --format csv
   pycsamt edi stations raw_edis/ --format csv
   pycsamt edi profile raw_edis/ --distances

Clean station selection:

.. code-block:: console

   pycsamt edi select raw_edis/ --min-nfreq 20 --keep-finite --dry-run
   pycsamt edi select raw_edis/ --min-nfreq 20 --keep-finite --output-dir clean_edis/
   pycsamt site info clean_edis/

AVG or Jones conversion:

.. code-block:: console

   pycsamt avg validate raw/L18.avg
   pycsamt avg correct raw/L18.avg --dry-run
   pycsamt avg correct raw/L18.avg --output-dir corrected_avg/
   pycsamt transform avg corrected_avg/ --output-dir converted_edis/
   pycsamt site info converted_edis/

Pipeline processing:

.. code-block:: console

   pycsamt survey set clean_edis/
   pycsamt pipe steps
   pycsamt pipe presets --expand basic_qc
   pycsamt pipe init --preset basic_qc --name line01_qc -o line01_qc.yaml
   pycsamt pipe run --config line01_qc.yaml --dry-run
   pycsamt pipe run --config line01_qc.yaml --out results/line01_qc/

Inversion preparation and review:

.. code-block:: console

   pycsamt invert build clean_edis/ --solver occam2d --workdir runs/line01_occam/
   pycsamt invert run runs/line01_occam/ --max-iter 100 --target-misfit 1.05
   pycsamt invert status runs/line01_occam/
   pycsamt invert results runs/line01_occam/ --format json
   pycsamt invert plot model runs/line01_occam/ --save figures/model.png
   pycsamt invert plot response runs/line01_occam/ --save figures/response.png

Interpretation export:

.. code-block:: console

   pycsamt interp rocks --rho 250
   pycsamt interp classify runs/line01_occam/ --output interpreted_layers.csv
   pycsamt interp export runs/line01_occam/ --format vtk --output-dir vtk/

Choosing Commands
-----------------

Start with the narrowest command that answers your question.

.. list-table::
   :header-rows: 1
   :widths: 32 30 38

   * - Question
     - First command
     - Next command
   * - Are these EDI files usable?
     - ``pycsamt edi validate``
     - ``pycsamt site info``
   * - What stations and coordinates are present?
     - ``pycsamt edi stations``
     - ``pycsamt map stations``
   * - What is the profile order?
     - ``pycsamt edi profile``
     - ``pycsamt site info``
   * - Should I keep or reject stations?
     - ``pycsamt edi select --dry-run``
     - ``pycsamt site select``
   * - What processing steps exist?
     - ``pycsamt pipe steps``
     - ``pycsamt pipe presets``
   * - How do I run a standard workflow?
     - ``pycsamt pipe presets --expand NAME``
     - ``pycsamt pipe run --preset NAME``
   * - How do I build inversion files?
     - ``pycsamt invert build``
     - ``pycsamt invert run``
   * - How do I inspect inversion output?
     - ``pycsamt invert status``
     - ``pycsamt invert results`` and ``invert plot``
   * - How do I make quick station figures?
     - ``pycsamt map stations``
     - ``pycsamt map plot``
   * - How do I classify model resistivity?
     - ``pycsamt interp rocks``
     - ``pycsamt interp classify``

Scripting Guidelines
--------------------

For reproducible shell scripts:

* pass explicit input paths instead of relying on active survey state;
* pass explicit output paths or output directories;
* use ``--format json`` for programmatic inspection;
* use ``--no-color`` when capturing logs;
* start destructive or expensive workflows with ``--dry-run``;
* keep raw inputs and generated products in separate directories;
* save the exact command sequence with important outputs.

Example script-friendly commands:

.. code-block:: console

   pycsamt --no-color edi validate raw_edis/ --format json
   pycsamt --no-color pipe run --survey clean_edis/ --config workflow.yaml --out results/qc/
   pycsamt --no-color invert status runs/line01_occam/ --format json

Troubleshooting
---------------

``pycsamt`` is not found
   Confirm that the environment containing pyCSAMT is active. During
   development, install from the repository root with
   ``python -m pip install -e ".[full]"``.

The command prints no colored tables
   ``rich`` may not be installed, or ``--no-color`` may be set. The
   plain text output is intentional fallback behavior.

A command cannot find survey data
   Pass the survey path explicitly or run ``pycsamt survey show`` to
   check the active survey. If no active survey is set, use
   ``pycsamt survey set PATH``.

The active survey looks stale
   Run ``pycsamt survey rebuild`` after changing station files. Use
   ``pycsamt survey rebuild --force`` when you want an unconditional
   reload.

JSON or CSV is needed but the command prints a table
   Check the command-specific help for ``--format``. If the output you
   need is not exposed as JSON or CSV, use the Python API instead of
   parsing terminal tables.

An output directory is empty
   Re-run with ``--dry-run`` if supported, check the input selector
   options, and confirm that the command was given an explicit
   ``--output``, ``--output-dir``, ``--out``, or ``--save`` argument as
   required by that command.

Getting Help
------------

Every command group and subcommand has local help:

.. code-block:: console

   pycsamt COMMAND --help
   pycsamt edi validate --help
   pycsamt pipe run --help
   pycsamt invert plot model --help

When in doubt, start with the group help, then inspect the subcommand
you intend to run. The command-specific documentation pages in this
section provide examples and workflow notes beyond the terminal help.

Relationship To The Python API
-------------------------------

Each command family is a convenience layer over package modules. The CLI
should be used when command-line composition is enough; the Python API
should be used when you need custom loops, plotting control, in-memory
data access, or integration with notebooks and applications.

Typical transitions are:

* from ``pycsamt site info`` to ``pycsamt.emtools._core.ensure_sites``;
* from ``pycsamt pipe run`` to the pipeline Python API;
* from ``pycsamt invert results`` to inversion result loader objects;
* from ``pycsamt map plot`` to the map plotting builders and view
  classes;
* from ``pycsamt interp classify`` to the interpretation and export
  functions.

The command pages include API notes where that transition is especially
useful.

