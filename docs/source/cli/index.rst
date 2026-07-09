.. _cli:

CLI
===

The pyCSAMT command-line interface is the shell entry point for the same
survey workflows exposed by the Python API: data inspection, format
conversion, site inventory, processing pipelines, forward modelling,
inversion preparation, interpretation, mapping, and TDEM utilities.

Use the CLI when you need repeatable commands for a field dataset,
batch-friendly outputs for automation, or a quick way to validate that a
survey directory can be loaded before writing Python code. The commands
are organized by workflow domain, but they share the same conventions:
explicit paths, structured output formats, dry-run modes where writes
are involved, and a project survey registry for remembered defaults.

Quick Start
-----------

The top-level command is ``pycsamt``:

.. code-block:: console

   pycsamt --help

Most workflows start by inspecting a directory of EDI files:

.. code-block:: console

   pycsamt edi validate data/edis/
   pycsamt site info data/edis/
   pycsamt map stations data/edis/

For repeated work on one survey line, store the active survey once:

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L18PLT
   pycsamt site info
   pycsamt map plot --output figures/l18_station_map.png
   pycsamt pipe run --preset basic_qc --out results/l18_qc/

The active survey is a convenience, not a hidden requirement. Commands
that accept a survey path can still be called with an explicit path, and
explicit paths should be preferred in scripts, reports, and CI jobs.

Start Here
----------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: CLI overview
      :link: overview
      :link-type: doc
      :img-top: ../_static/icons/cli-icon-terminal.svg
      :class-card: pycsamt-card sd-text-center

      High-level orientation to the CLI structure and the command
      patterns every family shares.

   .. grid-item-card:: Configuration
      :link: config
      :link-type: doc
      :img-top: ../_static/icons/developer-notes-icon-tools.svg
      :class-card: pycsamt-card sd-text-center

      Inspect runtime configuration, paths, defaults, and environment
      information.

   .. grid-item-card:: Survey registry
      :link: survey
      :link-type: doc
      :img-top: ../_static/icons/quickstart-cli.svg
      :class-card: pycsamt-card sd-text-center

      Manage active survey paths and registry-backed line defaults for
      repeated work on one dataset.

Data And Formats
----------------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: EDI commands
      :link: edi
      :link-type: doc
      :img-top: ../_static/icons/impedance.svg
      :class-card: pycsamt-card sd-text-center

      Validate, inventory, rotate, select, and profile EDI survey
      files.

   .. grid-item-card:: AVG commands
      :link: avg
      :link-type: doc
      :img-top: ../_static/icons/user-guide-data-loading.svg
      :class-card: pycsamt-card sd-text-center

      Inspect Zonge AVG files, validate components, apply corrections,
      and export cleaned products.

   .. grid-item-card:: Jones commands
      :link: jones
      :link-type: doc
      :img-top: ../_static/icons/glossary-icon-book.svg
      :class-card: pycsamt-card sd-text-center

      Inspect, validate, select, and read Jones ``.j`` files before
      conversion or analysis.

   .. grid-item-card:: Site commands
      :link: site
      :link-type: doc
      :img-top: ../_static/icons/user-guide-site.svg
      :class-card: pycsamt-card sd-text-center

      Load station collections, check usable impedance, compute
      summaries, and prepare site-level products.

   .. grid-item-card:: Transform commands
      :link: transform
      :link-type: doc
      :img-top: ../_static/icons/external-validation.svg
      :class-card: pycsamt-card sd-text-center

      Convert spectra, AVG, and Jones data into impedance EDI-style
      products.

   .. grid-item-card:: TDEM commands
      :link: tdem
      :link-type: doc
      :img-top: ../_static/icons/tdem.svg
      :class-card: pycsamt-card sd-text-center

      Work with transient electromagnetic survey files and conversion
      workflows.

Processing And Modelling
------------------------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Pipeline commands
      :link: pipe
      :link-type: doc
      :img-top: ../_static/icons/user-guide-pipeline.svg
      :class-card: pycsamt-card sd-text-center

      List steps, expand presets, dry-run pipelines, and execute
      reproducible processing chains.

   .. grid-item-card:: Forward commands
      :link: forward
      :link-type: doc
      :img-top: ../_static/icons/user-guide-forward.svg
      :class-card: pycsamt-card sd-text-center

      Run synthetic forward responses and generate training or
      benchmark datasets.

   .. grid-item-card:: Inversion commands
      :link: invert
      :link-type: doc
      :img-top: ../_static/icons/user-guide-inversion.svg
      :class-card: pycsamt-card sd-text-center

      Build inversion workdirs, launch solver runs, check status, read
      results, and plot products.

   .. grid-item-card:: Interpretation commands
      :link: interp
      :link-type: doc
      :img-top: ../_static/icons/user-guide-interpretation.svg
      :class-card: pycsamt-card sd-text-center

      Classify inversion results, export interpreted products, and map
      resistivity to rock or hydrogeologic labels.

   .. grid-item-card:: Map commands
      :link: map
      :link-type: doc
      :img-top: ../_static/icons/user-guide-mapping.svg
      :class-card: pycsamt-card sd-text-center

      List station coordinates and create station map figures from
      survey directories or the active survey.

Core Concepts
-------------

The CLI is designed around a few simple ideas.

Survey inputs
   Most science workflows operate on a directory of station files. For
   EDI-backed MT, AMT, and CSAMT workflows, the canonical loader is the
   same package loader used by the Python API. Stations without usable
   impedance are skipped or reported according to the command.

Active survey
   ``pycsamt survey set`` stores the current survey path for convenience.
   Commands such as ``site info``, ``map stations``, and ``pipe run`` can
   use that stored path when no path is supplied. Use ``--fresh`` on
   commands that support it when you want to reload from disk instead of
   reusing cached survey context.

Structured output
   Inspection commands commonly support human-readable tables plus
   ``--format json`` or ``--format csv``. Use JSON for automation, CSV
   for spreadsheets and GIS handoff, and the default table format for
   interactive terminal work.

Dry runs
   Commands that copy, convert, rotate, correct, or run a workflow often
   provide ``--dry-run``. Use it before modifying files or launching a
   longer operation. A good CLI workflow proves the selected inputs first
   and writes outputs second.

Output directories
   Prefer explicit output directories such as ``results/l18_qc/``,
   ``converted_edis/``, or ``figures/``. This keeps generated products
   separate from raw field data and makes reruns easier to compare.

Local APIs
   The CLI wraps package APIs; it is not a separate scientific engine.
   When you need more control, move from the CLI command to the Python
   API shown in the matching page's "Python API" section.

Common Workflows
----------------

Inventory a new EDI survey:

.. code-block:: console

   pycsamt edi validate raw_edis/
   pycsamt edi info raw_edis/ --format csv
   pycsamt edi stations raw_edis/ --format csv
   pycsamt edi profile raw_edis/ --distances

Select usable files and verify the cleaned directory:

.. code-block:: console

   pycsamt edi select raw_edis/ --min-nfreq 20 --keep-finite --dry-run
   pycsamt edi select raw_edis/ --min-nfreq 20 --keep-finite --output-dir clean_edis/
   pycsamt site info clean_edis/

Convert non-EDI products into impedance EDI outputs:

.. code-block:: console

   pycsamt transform avg raw_avg/ --output-dir converted_edis/
   pycsamt transform j raw_j/ --output-dir converted_edis/
   pycsamt site info converted_edis/

Register a survey and run a processing preset:

.. code-block:: console

   pycsamt survey set clean_edis/
   pycsamt pipe presets
   pycsamt pipe presets --expand basic_qc
   pycsamt pipe run --preset basic_qc --out results/basic_qc/

Build and inspect an inversion workspace:

.. code-block:: console

   pycsamt invert build clean_edis/ --solver occam2d --workdir runs/line01_occam/
   pycsamt invert status runs/line01_occam/
   pycsamt invert results runs/line01_occam/ --format json

Plot station and inversion outputs:

.. code-block:: console

   pycsamt map plot clean_edis/ --output figures/line01_stations.png
   pycsamt invert plot model runs/line01_occam/ --save figures/line01_model.png
   pycsamt invert plot response runs/line01_occam/ --save figures/line01_response.png

Interpret and export model products:

.. code-block:: console

   pycsamt interp rocks --rho 250
   pycsamt interp classify runs/line01_occam/ --output interpreted_layers.csv
   pycsamt interp export runs/line01_occam/ --format csv --output-dir exports/

Command Selection Guide
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 32 28 40

   * - I need to...
     - Start with
     - Then usually use
   * - Check whether field EDI files are loadable.
     - ``pycsamt edi validate``
     - ``pycsamt edi info`` and ``pycsamt site info``
   * - See station positions and profile order.
     - ``pycsamt edi stations``
     - ``pycsamt edi profile`` or ``pycsamt map stations``
   * - Convert AVG or Jones files before processing.
     - ``pycsamt transform avg`` or ``pycsamt transform j``
     - ``pycsamt site info`` on the output directory
   * - Keep only usable stations.
     - ``pycsamt edi select`` or ``pycsamt jones select``
     - ``--dry-run`` first, then write to a clean output directory
   * - Run standard processing steps.
     - ``pycsamt pipe presets``
     - ``pycsamt pipe run --preset ...``
   * - Design a custom processing chain.
     - ``pycsamt pipe steps``
     - ``pycsamt pipe init`` and ``pycsamt pipe show``
   * - Prepare a classical inversion.
     - ``pycsamt invert build``
     - ``pycsamt invert run``, ``status``, ``results``, and ``plot``
   * - Create synthetic data.
     - ``pycsamt forward run``
     - ``pycsamt forward generate``
   * - Turn inversion models into interpreted layers.
     - ``pycsamt interp classify``
     - ``pycsamt interp export``
   * - Make quick map figures.
     - ``pycsamt map stations``
     - ``pycsamt map plot``
   * - Work repeatedly on the same survey.
     - ``pycsamt survey set``
     - Any command that accepts the active survey

Output Conventions
------------------

Interactive table output is meant for a person at the terminal. It is
compact, sorted for readability, and may omit fields that are available
in structured formats.

JSON output is the safest choice for downstream scripts:

.. code-block:: console

   pycsamt edi validate raw_edis/ --format json
   pycsamt pipe show --preset basic_qc --format json
   pycsamt invert status runs/line01_occam/ --format json

CSV output is useful when the result is naturally tabular:

.. code-block:: console

   pycsamt edi stations raw_edis/ --format csv
   pycsamt map stations raw_edis/ --format csv --output station_coords.csv
   pycsamt interp classify runs/line01_occam/ --format csv

For commands that write files, use explicit output paths:

.. code-block:: console

   pycsamt map plot raw_edis/ --output figures/stations.png
   pycsamt transform avg raw_avg/ --output-dir converted_edis/
   pycsamt pipe run --preset basic_qc --out results/basic_qc/

Avoid mixing raw files and generated files in the same directory. Keep
``raw_*``, ``clean_*``, ``converted_*``, ``runs/``, ``results/``, and
``figures/`` separate unless a command page explicitly documents a
different convention.

Safety Practices
----------------

Use dry-run modes before writing:

.. code-block:: console

   pycsamt edi select raw_edis/ --min-nfreq 20 --dry-run
   pycsamt edi rotate raw_edis/ --angle 30 --dry-run
   pycsamt pipe run --config workflow.yaml --dry-run

Write corrected or transformed data to a new folder:

.. code-block:: console

   pycsamt avg correct raw/L18.avg --output-dir corrected_avg/
   pycsamt transform avg corrected_avg/ --output-dir clean_edis/

Keep command history with the output. For important processing or
inversion work, save the command line in a README, report, or pipeline
configuration next to the generated products.

For scripts and CI, prefer explicit paths over the active survey:

.. code-block:: console

   pycsamt pipe run --survey clean_edis/ --preset basic_qc --out results/basic_qc/

The active survey is excellent for exploratory terminal work, but
explicit paths make automated runs easier to reproduce.

Relationship to the Python API
------------------------------

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

Troubleshooting
---------------

Command not found
   Confirm that pyCSAMT is installed in the active environment. During
   development, install from the repository root with
   ``python -m pip install -e ".[full]"``.

No stations loaded
   Check the input path, file extension, and whether the station files
   contain usable impedance data. Start with ``edi validate`` and
   ``site info``.

Unexpected active survey
   Run ``pycsamt survey show`` or set the path again with
   ``pycsamt survey set PATH``. In scripts, pass the survey path
   explicitly.

Empty output directory
   Re-run with a smaller command and more explicit paths. For selectors
   and conversions, use ``--dry-run`` to confirm what would be written.

Long-running inversion or pipeline
   Use ``invert status`` or pipeline manifests where available. Keep run
   workdirs under ``runs/`` or ``results/`` so partial outputs are easy
   to inspect and remove intentionally.

Need machine-readable output
   Check whether the command supports ``--format json`` or
   ``--format csv``. If neither format is available for the exact output
   you need, use the Python API rather than parsing terminal tables.

.. toctree::
   :maxdepth: 1
   :hidden:

   overview
   config
   edi
   avg
   jones
   site
   transform
   forward
   invert
   map
   interp
   pipe
   tdem
   survey

