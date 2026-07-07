.. _pipeline-presets:

Pipeline Presets
================

Pipeline presets are named, opinionated processing recipes.  They are useful
when you want a tested starting point without writing a full step list by
hand.  A preset is not a hidden mode: internally it is just an ordered list of
``(label, Step)`` tuples, so it can be inspected, exported to a config file,
customized, and reviewed like any other pipeline.

Use presets when you want to:

* run a first-pass QC workflow quickly;
* compare a few standard processing strategies;
* scaffold a reproducible YAML, JSON, or Python config;
* teach new users a safe default order for common MT/AMT processing tasks;
* keep command-line workflows concise while preserving a saved
  ``pipeline.yaml`` in the output directory.

Preset Model
------------

Each built-in preset is represented by a ``Preset`` object with three fields:

``name``
    Stable identifier used by the CLI and Python API, for example
    ``basic_qc``.

``description``
    Short explanation shown in preset catalogues.

``steps``
    Ordered list of ``(label, Step)`` tuples.  The labels become report names
    and output subdirectory names; the ``Step`` objects hold registry codes
    and parameter defaults.

The preset API is intentionally small:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import get_preset, list_presets, preset_catalogue

   print(preset_catalogue())

   preset = get_preset("basic_qc")
   for label, step in preset.steps:
       print(label, step.spec.code, step.spec.name, step.params)

   for preset in list_presets():
       print(preset.name, len(preset.steps))

Run A Preset
------------

From the command line:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/basic_qc \
       -v

From Python:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_preset("basic_qc")
   result = pipe.run(sites, outdir="results/basic_qc")

Inspect Presets Before Running
------------------------------

List all presets:

.. code-block:: console
   :linenos:

   pycsamt pipe presets

Expand one preset into its step sequence:

.. code-block:: console
   :linenos:

   pycsamt pipe presets --expand full_processing
   pycsamt pipe show --preset full_processing

Use JSON or CSV when another tool needs the preset list:

.. code-block:: console
   :linenos:

   pycsamt pipe presets --format json
   pycsamt pipe presets --expand basic_qc --format json
   pycsamt pipe presets --format csv

Built-In Preset Summary
-----------------------

The normal pipeline registry currently provides seven presets.

.. list-table::
   :header-rows: 1
   :widths: 24 16 36 24

   * - Preset
     - Steps
     - Best for
     - Sequence
   * - ``basic_qc``
     - 5
     - First-pass inspection and quick survey sanity checks.
     - ``NR001``, ``FREQ002``, ``FREQ001``, ``FREQ004``, ``QC001``
   * - ``noise_reduction``
     - 6
     - High-EMI data where denoising is the main question.
     - ``NR001``, ``NR004``, ``NR005``, ``NR003``, ``NR010``, ``QC001``
   * - ``full_processing``
     - 8
     - Standard end-to-end processing before interpretation.
     - ``NR001``, ``FREQ002``, ``FREQ001``, ``FREQ004``, ``SK001``,
       ``TZ001``, ``SS001``, ``QC001``
   * - ``tensor_analysis``
     - 5
     - Tensor-focused cleanup after data already have an acceptable
       frequency grid.
     - ``TZ001``, ``TZ002``, ``TZ003``, ``TZ004``, ``QC001``
   * - ``dimensionality_filter``
     - 4
     - Classifying dimensionality and keeping/projecting 2-D-compatible
       intervals.
     - ``DIM001``, ``DIM002``, ``DIM003``, ``QC001``
   * - ``publication_ready``
     - 9
     - Longer reviewed workflow for polished reports and figures.
     - ``NR001``, ``FREQ002``, ``FREQ001``, ``FREQ004``, ``SS001``,
       ``TZ001``, ``TZ002``, ``SK001``, ``QC001``
   * - ``stratagem_mt``
     - 7
     - Stratagem AMT data already loaded as a ``Sites`` object.
     - ``SS001``, ``FREQ001``, ``FREQ002``, ``NR001``, ``NR004``,
       ``NR010``, ``QC001``

Choosing A Preset
-----------------

Start with the narrowest preset that answers your current question.

.. list-table::
   :header-rows: 1
   :widths: 34 33 33

   * - Situation
     - Start with
     - Why
   * - You just received a survey and need a quick check.
     - ``basic_qc``
     - It does only basic notch, frequency cleanup, alignment, and QC.
   * - Harmonic and spatial noise dominate the data.
     - ``noise_reduction``
     - It stacks targeted noise-removal steps before a QC snapshot.
   * - You want a general processing run before inversion preparation.
     - ``full_processing``
     - It combines denoising, frequency cleanup, skew gating, rotation,
       static-shift correction, and QC.
   * - You already trust the frequency grid and want tensor diagnostics.
     - ``tensor_analysis``
     - It avoids frequency and noise steps and focuses on tensor operations.
   * - You are deciding which intervals are compatible with 2-D assumptions.
     - ``dimensionality_filter``
     - It classifies dimensionality, masks by class, projects to 2-D, and
       generates QC.
   * - You need a polished, repeatable processing chain.
     - ``publication_ready``
     - It is longer and more opinionated, with static-shift correction,
       tensor cleanup, skew gating, and final QC.
   * - You are working with Stratagem AMT data already represented as
       ``Sites``.
     - ``stratagem_mt``
     - It applies Stratagem-oriented AMT band selection, static shift,
       denoising, and QC at the emtools pipeline level.

Preset Details
--------------

basic_qc
~~~~~~~~

``basic_qc`` is the safest first preset for most surveys.  It does not try to
solve every processing problem; it prepares a clean enough view to understand
what the data need next.

Sequence:

.. code-block:: text
   :linenos:

   notch            NR001   notch_powerline
   drop_duplicates  FREQ002 drop_duplicates
   select_band      FREQ001 select_band
   align_grid       FREQ004 align_grid
   qc_snapshot      QC001   qc_snapshot

Run:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/basic_qc \
       --on-error warn

Use ``basic_qc`` when:

* you need quick figures before committing to a processing plan;
* you want to verify that EDI loading and output generation work;
* you are comparing surveys and want the same minimal cleanup everywhere.

Move beyond ``basic_qc`` when:

* static shift is obvious;
* skew or dimensionality gates are needed;
* power-line removal is not enough for the noise environment;
* you need a processing chain suitable for inversion preparation.

noise_reduction
~~~~~~~~~~~~~~~

``noise_reduction`` concentrates on denoising.  It is useful when the first
inspection shows power-line harmonics, local spikes, spatially coherent
outliers, or incoherent frequency bins.

Sequence:

.. code-block:: text
   :linenos:

   notch         NR001 notch_powerline
   hampel        NR004 hampel_filter
   spatial_med   NR005 spatial_median
   shrink_trend  NR003 shrink_group_trend
   mask_incoher  NR010 mask_incoherent
   qc_snapshot   QC001 qc_snapshot

Run:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset noise_reduction \
       --out results/noise_reduction

Use this preset to compare denoising impact against ``basic_qc``:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc \
       --out results/compare/basic_qc
   pycsamt pipe run data/edis --preset noise_reduction \
       --out results/compare/noise_reduction

full_processing
~~~~~~~~~~~~~~~

``full_processing`` is the standard end-to-end workflow.  It starts with
noise and frequency cleanup, then applies a skew gate, strike rotation,
static-shift correction, and QC.

Sequence:

.. code-block:: text
   :linenos:

   notch          NR001   notch_powerline
   drop_dup       FREQ002 drop_duplicates
   select_band    FREQ001 select_band
   align_grid     FREQ004 align_grid
   mask_skew      SK001   mask_by_skew
   rotate_strike  TZ001   rotate_strike
   correct_ss     SS001   correct_ss_ama
   qc_snapshot    QC001   qc_snapshot

Run:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset full_processing \
       --out results/full_processing \
       -v

Use this preset when:

* the survey needs a broad processing pass;
* you want an auditable default before building an inversion-specific config;
* you need one chain that exercises the main processing families.

tensor_analysis
~~~~~~~~~~~~~~~

``tensor_analysis`` assumes the data are already in reasonable condition and
focuses on tensor operations.

Sequence:

.. code-block:: text
   :linenos:

   rotate_strike  TZ001 rotate_strike
   antisymm       TZ002 antisymmetrize
   sigma_clip     TZ003 sigma_clip
   balance        TZ004 balance_offdiag
   qc_snapshot    QC001 qc_snapshot

Use it when you want to inspect tensor behavior without changing the
frequency selection or applying the broader denoising chain.

dimensionality_filter
~~~~~~~~~~~~~~~~~~~~~

``dimensionality_filter`` is for 1-D / 2-D / 3-D screening and 2-D projection
workflows.

Sequence:

.. code-block:: text
   :linenos:

   classify_dim  DIM001 classify_dim
   mask_dim      DIM002 mask_by_dim
   project_2d    DIM003 project_2d
   qc_snapshot   QC001  qc_snapshot

Run:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset dimensionality_filter \
       --out results/dimensionality

Use it after basic cleanup when the main question is whether the remaining
intervals are compatible with a 2-D interpretation or inversion assumption.

publication_ready
~~~~~~~~~~~~~~~~~

``publication_ready`` is the longest built-in general-purpose preset.  It is
designed for polished processing output rather than quick exploration.

Sequence:

.. code-block:: text
   :linenos:

   notch          NR001   notch_powerline
   drop_dup       FREQ002 drop_duplicates
   select_band    FREQ001 select_band
   align_grid     FREQ004 align_grid
   correct_ss     SS001   correct_ss_ama
   rotate_strike  TZ001   rotate_strike
   antisymm       TZ002   antisymmetrize
   mask_skew      SK001   mask_by_skew
   qc_snapshot    QC001   qc_snapshot

Run:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset publication_ready \
       --out results/publication_ready \
       --dpi 300 \
       --plot-fmt pdf \
       -v

Use this preset when:

* you already inspected the data with a lighter preset;
* the default step order is scientifically acceptable for your survey;
* you need high-quality saved figures and a complete run report.

stratagem_mt
~~~~~~~~~~~~

``stratagem_mt`` is a normal emtools pipeline preset specialized for
Stratagem AMT data that are already loaded as a ``Sites`` object.  It does
not perform raw-coordinate injection, raw hardware-file parsing, or station
renaming by itself.

Sequence:

.. code-block:: text
   :linenos:

   correct_ss    SS001   correct_ss_ama
   select_band   FREQ001 select_band   band_hz=(10.0, 100000.0)
   drop_dup      FREQ002 drop_duplicates
   notch         NR001   notch_powerline
   hampel        NR004   hampel_filter
   mask_incoher  NR010   mask_incoherent
   qc_snapshot   QC001   qc_snapshot

Use this preset with the normal pipeline API when your input is already a
site collection:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_preset("stratagem_mt")
   result = pipe.run(sites, outdir="results/stratagem_mt")

Use :class:`pycsamt.pipeline.stratagem.StratagemPipeline` or
``run_stratagem_preset`` when you also need the full raw EDI plus GPS CSV
workflow.

Export A Preset To A Config
---------------------------

For serious work, use a preset to generate a config and then commit the
expanded recipe to your project.  This makes the workflow auditable and easy
to rerun.

Generate YAML:

.. code-block:: console
   :linenos:

   pycsamt pipe init \
       --preset publication_ready \
       --name line22_publication_ready \
       --outdir results/line22_publication_ready \
       --output config/line22_publication_ready.yaml

Generate Python:

.. code-block:: console
   :linenos:

   pycsamt pipe init \
       --preset basic_qc \
       --format py \
       --output config/basic_qc.py

Preview before running:

.. code-block:: console
   :linenos:

   pycsamt pipe show config/line22_publication_ready.yaml
   pycsamt pipe run data/edis \
       --config config/line22_publication_ready.yaml \
       --dry-run

Customizing Presets Safely
--------------------------

There are three ways to customize a preset.

Use the preset directly, then edit the pipeline in Python:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline, Step

   pipe = Pipeline.from_preset("basic_qc")
   pipe.replace("notch", Step("NR001", mains_hz=60, n_harm=25))
   pipe.append("static_shift", Step("SS001"))

Use ``pycsamt pipe init`` to expand a preset into an explicit config, then
edit the generated step parameters:

.. code-block:: yaml
   :linenos:

   name: basic_qc_60hz
   output_dir: results/basic_qc_60hz

   steps:
     - name: notch
       code: NR001
       params:
         mains_hz: 60
         n_harm: 25
         tol_hz: 0.08
     - name: drop_duplicates
       code: FREQ002
     - name: select_band
       code: FREQ001
     - name: align_grid
       code: FREQ004
     - name: qc_snapshot
       code: QC001

Append extra steps after a preset in a config:

.. code-block:: yaml
   :linenos:

   name: basic_qc_plus_static_shift
   output_dir: results/basic_qc_plus_static_shift
   preset: basic_qc

   steps:
     - name: static_shift
       code: SS001
     - name: final_qc
       code: QC001

Important: in a config file, ``preset: basic_qc`` loads all preset steps
first, then appends the explicit ``steps`` list.  It does not modify an
existing preset step.  If you need to change ``NR001`` from 50 Hz to 60 Hz,
use an explicit expanded step list instead of ``preset: basic_qc`` plus a
second ``NR001`` step.

CLI Priority
------------

``pycsamt pipe run`` resolves the pipeline definition in this order:

1. ``--config FILE``;
2. ``--preset NAME``;
3. ``--steps CODE,CODE,...``.

If ``--config`` is supplied, ``--preset`` and ``--steps`` are ignored because
the config file is the source of truth.

Examples:

.. code-block:: console
   :linenos:

   # Uses the config.  The preset argument is ignored.
   pycsamt pipe run data/edis \
       --config config/basic_qc.yaml \
       --preset publication_ready

   # Uses the preset.
   pycsamt pipe run data/edis \
       --preset publication_ready

   # Uses the ad-hoc steps.
   pycsamt pipe run data/edis \
       --steps FREQ002,FREQ001,FREQ004,NR001,QC001

Compare Presets
---------------

A useful way to choose a recipe is to run several presets into separate
output directories and compare the reports:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/compare/basic_qc

   pycsamt pipe run data/edis \
       --preset noise_reduction \
       --out results/compare/noise_reduction

   pycsamt pipe run data/edis \
       --preset full_processing \
       --out results/compare/full_processing

Compare:

* ``summary.txt`` for step failures, runtime, and site counts;
* ``report.html`` for per-step status and embedded pipeline YAML;
* ``plots/`` for visual differences between cleanup strategies;
* ``processed/`` for exported EDI differences.

Stratagem Presets
-----------------

There are two related but different Stratagem preset systems.

``stratagem_mt``
    A normal :class:`pycsamt.pipeline.Pipeline` preset.  It expects data that
    can already be processed as ``Sites`` and runs normal registered pipeline
    steps.

``StratagemPreset``
    A convenience workflow in :mod:`pycsamt.pipeline.stratagem` for raw
    Stratagem EDI directories plus coordinate CSV files.  It calls
    ``StratagemSurvey`` methods such as ``remove_static_shift``,
    ``drop_frequencies``, ``remove_noises``, ``export``, and ``rename``.

The Stratagem convenience presets are:

.. list-table::
   :header-rows: 1
   :widths: 24 28 48

   * - Preset
     - Main workflow
     - Best for
   * - ``basic``
     - coordinate injection, AMA static shift, frequency trim, noise removal,
       export, rename
     - Direct replacement for the legacy Stratagem processing script.
   * - ``full_processing``
     - QC, AMA static shift, hardware-aware frequency filtering, smoothed
       noise removal, export, rename
     - Raw Stratagem workflows with hardware files and a full QC pass.
   * - ``publication_ready``
     - stricter QC, hardware SNR masking, AMT band trimming, stronger
       smoothing, export, rename
     - Polished Stratagem outputs after the basic workflow has been reviewed.

Run a raw Stratagem convenience preset:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline.stratagem import run_stratagem_preset

   survey = run_stratagem_preset(
       "full_processing",
       edi_dir="2/2EDI",
       coord_file="2.csv",
       raw_dir="raw/2HX",
       outdir="results/stratagem",
       epsg=32649,
       utm_zone="49N",
       rename_basename="T2.",
       overwrite=True,
       verbose=1,
   )

Build a Stratagem pipeline object from a normal emtools preset:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline.stratagem import StratagemPipeline

   pipe = StratagemPipeline.from_preset(
       "stratagem_mt",
       coord_file="2.csv",
       raw_dir="raw/2HX",
       epsg=32649,
       utm_zone="49N",
       rename_basename="T2.",
   )

   result = pipe.run("2/2EDI", outdir="results/stratagem_mt")

Troubleshooting
---------------

Unknown preset
    Run ``pycsamt pipe presets``.  Preset names are exact and lowercase, for
    example ``basic_qc`` or ``publication_ready``.

I changed ``preset: basic_qc`` but the notch is still 50 Hz
    A config ``preset`` expands the preset first.  Explicit ``steps`` are
    appended; they do not edit existing preset steps.  Generate an expanded
    config with ``pycsamt pipe init --preset basic_qc`` and edit the
    ``NR001`` parameters directly.

The preset is too aggressive
    Move to a narrower preset such as ``basic_qc`` or export the preset to a
    config and remove the steps that are not justified by the data.

The preset does not include a step I need
    Append the step in Python, or add it to an explicit config after
    generating the preset scaffold.

I need raw Stratagem coordinate injection and renaming
    Use :mod:`pycsamt.pipeline.stratagem` rather than the normal
    ``stratagem_mt`` preset alone.

Related Pages
-------------

* :doc:`steps` explains the registered step codes used by each preset.
* :doc:`configuration_files` explains how presets expand inside YAML, JSON,
  and Python configs.
* :doc:`cli_pipe` explains ``pycsamt pipe presets``, ``pycsamt pipe show``,
  and ``pycsamt pipe run``.
* :doc:`concepts` explains pipeline execution and result objects.
* :doc:`outputs` explains generated reports, plots, processed EDIs, and saved
  pipeline snapshots.
