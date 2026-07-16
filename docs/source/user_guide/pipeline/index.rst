.. _pipeline:

Pipeline System
===============

The pyCSAMT pipeline system is the reproducible processing layer for MT/AMT
survey workflows.  It connects registered processing steps, named presets,
configuration files, command-line execution, output directories, reports, and
Python result objects into one coherent workflow.

Use the pipeline when you want to:

* run the same processing chain on multiple surveys;
* inspect and compare quality-control figures after each step;
* save processed EDI files without modifying raw inputs;
* keep an exact ``pipeline.yaml`` snapshot for reproducibility;
* move between the Python API and the ``pycsamt pipe`` CLI;
* standardize processing recipes for a project, publication, or team.

The pipeline API follows a scikit-learn-style shape: a pipeline is an ordered
sequence of labelled steps, and each step is a configured operation from the
registry.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   concepts
   configuration_files
   cli_pipe
   presets
   steps
   outputs

Minimal Example
---------------

Run a built-in preset from the command line:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/basic_qc \
       --on-error warn \
       -v

The same idea from Python:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_preset("basic_qc")
   result = pipe.run(sites, outdir="results/basic_qc")

   print(result.summary())

A typical output directory contains:

.. code-block:: text
   :linenos:

   results/basic_qc/
   |-- processed/
   |-- plots/
   |-- pipeline.yaml
   |-- report.html
   `-- summary.txt

Documentation Map
-----------------

.. list-table::
   :header-rows: 1
   :widths: 28 42 30

   * - Page
     - What it explains
     - Read this when
   * - :doc:`concepts`
     - The pipeline object model, run lifecycle, runtime configuration,
       ``PipelineResult``, and CLI/Python equivalence.
     - You want the mental model before writing workflows.
   * - :doc:`configuration_files`
     - YAML, JSON, and Python config files; ``preset`` expansion; explicit
       step lists; output directories; project organization.
     - You want reproducible processing recipes.
   * - :doc:`cli_pipe`
     - ``pycsamt pipe presets``, ``steps``, ``init``, ``show``, and ``run``.
     - You process surveys from the terminal or automation scripts.
   * - :doc:`presets`
     - Built-in recipes such as ``basic_qc``, ``noise_reduction``,
       ``full_processing``, ``publication_ready``, and ``stratagem_mt``.
     - You want a safe starting workflow.
   * - :doc:`steps`
     - All registered step codes, categories, defaults, ordering guidance,
       QC plots, and extension policy.
     - You need to choose, tune, or add processing operations.
   * - :doc:`outputs`
     - Processed EDIs, plots, reports, ``pipeline.yaml``, in-memory runs,
       ``PipelineResult``, and output troubleshooting.
     - You need to understand what a run writes and how to inspect it.

Recommended Reading Paths
-------------------------

New users
    Start with :doc:`concepts`, then run ``basic_qc`` from :doc:`cli_pipe`,
    then inspect the files described in :doc:`outputs`.

Users preparing a real survey
    Read :doc:`presets`, generate a config with ``pycsamt pipe init``, edit it
    using :doc:`configuration_files`, and validate the selected operations
    with :doc:`steps`.

Users building automated workflows
    Read :doc:`cli_pipe` for terminal execution and :doc:`outputs` for
    machine-readable result summaries and directory contracts.

Python users
    Read :doc:`concepts` first, then :doc:`steps` for ``Step`` construction,
    and :doc:`outputs` for ``PipelineResult`` inspection.

Developers adding operations
    Read :doc:`steps` for the registry and extension policy, then update
    :doc:`presets` if the new operation belongs in a built-in recipe.

Core Concepts At A Glance
-------------------------

``Pipeline``
    Ordered sequence of labelled ``Step`` objects.  It can be built directly
    in Python, from a preset, from a config file, or from the CLI.

``Step``
    Configured processing operation.  It references a registry code such as
    ``NR001`` or a snake-case name such as ``notch_powerline`` and carries
    parameter overrides.

``StepSpec``
    Registry metadata for a step: code, name, label, category, defaults,
    transform function, QC functions, and whether the step returns modified
    sites.

``Preset``
    Named ordered step list for common workflows.  Presets are convenient
    starting points, not hidden black boxes.

``PipelineResult``
    Programmatic record returned by ``Pipeline.run``.  It contains input and
    output sites, per-step results, plot paths, processed EDI paths, status,
    errors, runtime, and output directory.

``pipeline.yaml``
    Saved snapshot of the resolved pipeline.  It is written into every
    output-enabled run and can be reloaded for reproduction.

Common Workflows
----------------

First-pass QC:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/basic_qc

Generate an editable config from a preset:

.. code-block:: console
   :linenos:

   pycsamt pipe init \
       --preset publication_ready \
       --name line22_publication_ready \
       --outdir results/line22_publication_ready \
       --output config/line22_publication_ready.yaml

Preview before processing:

.. code-block:: console
   :linenos:

   pycsamt pipe show config/line22_publication_ready.yaml
   pycsamt pipe run data/edis \
       --config config/line22_publication_ready.yaml \
       --dry-run

Run a reviewed config:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/line22_publication_ready.yaml \
       --out results/line22_publication_ready \
       --on-error warn \
       -v

Run in memory from Python:

.. code-block:: python
   :linenos:

   result = pipe.run(
       sites,
       outdir=None,
       save_plots=False,
       save_edis=False,
       save_report=False,
   )

Pipeline Design Rules
---------------------

Keep raw and processed data separate
    Raw EDI directories should remain unchanged.  Write processed EDIs,
    plots, reports, and snapshots to a dedicated ``results/`` tree.

Prefer configs for repeatable work
    Direct presets are useful for exploration.  Reviewed survey processing
    should usually live in a YAML, JSON, or Python config file.

Use presets as starting points
    If a preset needs changed parameters, export it to an explicit config
    rather than appending duplicate steps that look similar but do not modify
    the preset step in place.

Inspect before running
    Use ``pycsamt pipe show`` and ``--dry-run`` before expensive processing
    or before writing large output directories.

Name step labels carefully
    Labels appear in reports, plot folders, and slicing options such as
    ``--from-step`` and ``--until-step``.

Archive ``pipeline.yaml``
    The saved pipeline snapshot is the processing recipe that connects raw
    inputs, outputs, figures, and reports.

Typical Project Layout
----------------------

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

This layout makes it easy to compare runs without mixing raw observations,
configuration files, and generated outputs.

Next Steps
----------

If you are learning the system, continue with :doc:`concepts`.  If you
already want to run data, jump to :doc:`cli_pipe`.  If you are ready to write
project recipes, start with :doc:`configuration_files`.
