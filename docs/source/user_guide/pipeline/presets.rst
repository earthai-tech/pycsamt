.. _pipeline-presets:

Pipeline Presets
================

A :term:`pipeline preset` is a named, built-in processing recipe.  It gives a
survey a tested starting order without forcing the user to write every
:term:`pipeline step code` by hand.  A preset is not a hidden processing mode:
inside pyCSAMT it is an ordered list of ``(label, Step)`` tuples.  That list
can be inspected, expanded, exported to a :term:`pipeline configuration file`,
edited, reviewed, and archived like any other :term:`processing pipeline`.

Use presets for quick first-pass QC, standard comparisons between processing
strategies, teaching safe default step order, and generating starter configs.
For serious survey processing, treat the preset as the beginning of the
conversation, not the final scientific decision.  The reproducible record is
the expanded step sequence written to ``pipeline.yaml`` or to a committed
configuration file.

Preset Mental Model
-------------------

Each built-in preset is represented by a small object with three fields:

``name``
    Stable identifier used by the CLI and Python API, for example
    ``basic_qc``.

``description``
    Short explanation shown in preset catalogues.

``steps``
    Ordered list of ``(label, Step)`` tuples.  The label becomes the
    :term:`step label` used in reports and plot directories.  The ``Step``
    object holds the registry code and the parameters passed to the transform.

Mathematically, a preset defines a fixed composition of transforms:

.. math::

   P(S_0) =
   T_n(\cdots T_2(T_1(S_0;\theta_1);\theta_2)\cdots;\theta_n),

where :math:`S_0` is the input :term:`site collection`, :math:`T_j` is the
registered operation at position :math:`j`, and :math:`\theta_j` is the
resolved parameter dictionary for that step.  The order is part of the
meaning.  Running ``FREQ001`` before ``NR001`` is not the same workflow as
running ``NR001`` before ``FREQ001`` because each step receives a different
intermediate survey state.

The preset API is intentionally small:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import get_preset, list_presets, preset_catalogue
   >>> preset = get_preset("basic_qc")
   >>> preset.name
   'basic_qc'
   >>> [(label, step.spec.code) for label, step in preset.steps]
   [('notch', 'NR001'), ('drop_duplicates', 'FREQ002'), ('select_band', 'FREQ001'), ('align_grid', 'FREQ004'), ('qc_snapshot', 'QC001')]
   >>> [(p.name, len(p.steps)) for p in list_presets()]
   [('basic_qc', 5), ('noise_reduction', 6), ('full_processing', 8), ('tensor_analysis', 5), ('dimensionality_filter', 4), ('publication_ready', 9), ('stratagem_mt', 7)]

The CLI exposes the same information:

.. code-block:: console
   :linenos:

   pycsamt pipe presets
   pycsamt pipe presets --format json
   pycsamt pipe presets --expand basic_qc --format json

Captured expansion excerpt:

.. code-block:: json
   :linenos:

   {
     "name": "basic_qc",
     "description": "Minimal denoising + frequency cleanup.  Good for quick-look inspection.",
     "n_steps": 5,
     "steps": [
       {
         "label": "notch",
         "code": "NR001",
         "name": "notch_powerline",
         "category": "noise_removal",
         "params": {"mains_hz": 50, "n_harm": 30, "tol_hz": 0.08}
       },
       {
         "label": "drop_duplicates",
         "code": "FREQ002",
         "name": "drop_duplicates",
         "category": "frequency",
         "params": {}
       }
     ]
   }

The excerpt is deliberately partial.  In practice, use the full expansion when
reviewing or exporting a preset.

Run A Preset
------------

From the command line:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
       --preset basic_qc \
       --out results/basic_qc \
       -v

From Python:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> pipe = Pipeline.from_preset("basic_qc")
   >>> [(label, step.spec.code) for label, step in pipe]
   [('notch', 'NR001'), ('drop_duplicates', 'FREQ002'), ('select_band', 'FREQ001'), ('align_grid', 'FREQ004'), ('qc_snapshot', 'QC001')]
   >>> result = pipe.run(sites, outdir="results/basic_qc")

When output is enabled, the run writes ``results/basic_qc/pipeline.yaml``.
That file is the run-specific :term:`canonical pipeline snapshot`; it records
the step sequence that actually ran.

Built-In Preset Summary
-----------------------

The normal pipeline registry currently provides seven presets.

.. list-table::
   :header-rows: 1
   :widths: 24 14 38 24

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
     - Tensor-focused cleanup after frequency and noise handling are already
       acceptable.
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

Start with the narrowest preset that answers the current question.  A narrow
preset is easier to diagnose because fewer transforms stand between
:math:`S_0` and the output state :math:`S_n`.

.. list-table::
   :header-rows: 1
   :widths: 34 30 36

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

.. code-block:: text
   :linenos:

   notch            NR001   notch_powerline
   drop_duplicates  FREQ002 drop_duplicates
   select_band      FREQ001 select_band
   align_grid       FREQ004 align_grid
   qc_snapshot      QC001   qc_snapshot

The default transform is

.. math::

   S_5 =
   T_{\mathrm{QC001}}(
   T_{\mathrm{FREQ004}}(
   T_{\mathrm{FREQ001}}(
   T_{\mathrm{FREQ002}}(
   T_{\mathrm{NR001}}(S_0))))).

Use ``basic_qc`` when you need quick figures before committing to a processing
plan, want to verify that EDI loading and output generation work, or need the
same minimal cleanup across several surveys.  Move beyond it when static
shift, skew gating, dimensionality filtering, or stronger denoising is clearly
justified by the data.

noise_reduction
~~~~~~~~~~~~~~~

``noise_reduction`` concentrates on denoising.  It is useful when the first
inspection shows power-line harmonics, local spikes, spatially coherent
outliers, or incoherent frequency bins.

.. code-block:: text
   :linenos:

   notch         NR001 notch_powerline
   hampel        NR004 hampel_filter
   spatial_med   NR005 spatial_median
   shrink_trend  NR003 shrink_group_trend
   mask_incoher  NR010 mask_incoherent
   qc_snapshot   QC001 qc_snapshot

Run it next to ``basic_qc`` rather than replacing the first pass silently:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis --preset basic_qc \
       --out results/compare/basic_qc
   pycsamt pipe run data/3edis --preset noise_reduction \
       --out results/compare/noise_reduction

full_processing
~~~~~~~~~~~~~~~

``full_processing`` is the standard end-to-end workflow.  It starts with noise
and frequency cleanup, then applies a skew gate, strike rotation,
static-shift correction, and QC.

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

Use this preset when the survey needs a broad processing pass, when you want
an auditable default before building an inversion-specific config, or when you
need one chain that exercises the main processing families.

tensor_analysis
~~~~~~~~~~~~~~~

``tensor_analysis`` assumes the data are already in reasonable condition and
focuses on tensor operations.

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

.. code-block:: text
   :linenos:

   classify_dim  DIM001 classify_dim
   mask_dim      DIM002 mask_by_dim
   project_2d    DIM003 project_2d
   qc_snapshot   QC001  qc_snapshot

Use it after basic cleanup when the main question is whether the remaining
intervals are compatible with a :term:`2-D` interpretation or inversion
assumption.

publication_ready
~~~~~~~~~~~~~~~~~

``publication_ready`` is the longest built-in general-purpose preset.  It is
designed for polished processing output rather than quick exploration.

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

   pycsamt pipe run data/3edis \
       --preset publication_ready \
       --out results/publication_ready \
       --dpi 300 \
       --plot-fmt pdf \
       -v

Use this preset after a lighter pass has shown that its assumptions are
reasonable for the survey.  The name does not make the output publication
ready by itself; the review comes from inspecting the report, figures,
processed EDIs, and saved pipeline snapshot.

stratagem_mt
~~~~~~~~~~~~

``stratagem_mt`` is a normal emtools pipeline preset specialized for Stratagem
AMT data that are already loaded as a ``Sites`` object.  It does not perform
raw-coordinate injection, raw hardware-file parsing, or station renaming by
itself.

.. code-block:: text
   :linenos:

   correct_ss    SS001   correct_ss_ama
   select_band   FREQ001 select_band   band_hz=(10.0, 100000.0)
   drop_dup      FREQ002 drop_duplicates
   notch         NR001   notch_powerline
   hampel        NR004   hampel_filter
   mask_incoher  NR010   mask_incoherent
   qc_snapshot   QC001   qc_snapshot

Use it with the normal pipeline API when your input is already a site
collection:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> pipe = Pipeline.from_preset("stratagem_mt")
   >>> result = pipe.run(sites, outdir="results/stratagem_mt")

Use :class:`pycsamt.pipeline.stratagem.StratagemPipeline` or
``run_stratagem_preset`` when you also need the full raw EDI plus GPS CSV
workflow.

Export A Preset To A Config
---------------------------

:term:`Preset expansion` is the safest transition from exploration to
reproducible work.  Generate the preset once, review the explicit list, then
commit the config with the survey project.

.. code-block:: console
   :linenos:

   pycsamt pipe init \
       --preset publication_ready \
       --name line22_publication_ready \
       --outdir results/line22_publication_ready \
       --output config/line22_publication_ready.yaml

Preview before running:

.. code-block:: console
   :linenos:

   pycsamt pipe show config/line22_publication_ready.yaml
   pycsamt pipe run data/3edis \
       --config config/line22_publication_ready.yaml \
       --dry-run

Generate Python or JSON instead:

.. code-block:: console
   :linenos:

   pycsamt pipe init --preset basic_qc --format py \
       --output config/basic_qc.py
   pycsamt pipe init --preset basic_qc --format json \
       --output config/basic_qc.json

Customizing Presets Safely
--------------------------

There are three supported customization patterns.

Edit a pipeline object in Python:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step
   >>> pipe = Pipeline.from_preset("basic_qc")
   >>> pipe.replace("notch", Step("NR001", mains_hz=60, n_harm=25)) is pipe
   True
   >>> pipe.append("static_shift", Step("SS001")) is pipe
   True

Expand a preset into an explicit config and edit the generated parameters:

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

The third pattern appends.  It does not edit the preset.  In symbols, the
loaded step list is

.. math::

   S_{\mathrm{final}} = S_{\mathrm{preset}} \Vert S_{\mathrm{explicit}},

where ``\Vert`` means append.  If you need to change ``NR001`` from 50 Hz to
60 Hz, use an explicit expanded step list instead of ``preset: basic_qc`` plus
a second ``NR001`` step.

CLI Priority
------------

``pycsamt pipe run`` resolves the pipeline definition in this order:

1. ``--config FILE``;
2. ``--preset NAME``;
3. ``--steps CODE,CODE,...``.

If ``--config`` is supplied, ``--preset`` and ``--steps`` are ignored because
the config file is the source of truth.

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
       --config config/basic_qc.yaml \
       --preset publication_ready

   pycsamt pipe run data/3edis \
       --preset publication_ready

   pycsamt pipe run data/3edis \
       --steps FREQ002,FREQ001,FREQ004,NR001,QC001

Compare Presets
---------------

A :term:`preset comparison run` keeps the input fixed and changes only the
recipe.  Let :math:`P_a` and :math:`P_b` be two presets.  The comparison is
meaningful only when both are applied to the same initial state:

.. math::

   S_n^{(a)} = P_a(S_0), \qquad S_m^{(b)} = P_b(S_0).

Run the branches into separate output roots:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
       --preset basic_qc \
       --out results/compare/basic_qc

   pycsamt pipe run data/3edis \
       --preset noise_reduction \
       --out results/compare/noise_reduction

   pycsamt pipe run data/3edis \
       --preset full_processing \
       --out results/compare/full_processing

Compare ``summary.txt`` for step failures and site counts, ``report.html`` for
per-step status and embedded YAML, ``plots/`` for visual differences,
``processed/`` for exported EDI differences, and ``pipeline.yaml`` for the
exact recipe behind each branch.

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
     - Coordinate injection, AMA static shift, frequency trim, noise removal,
       export, rename.
     - Direct replacement for the legacy Stratagem processing script.
   * - ``full_processing``
     - QC, AMA static shift, hardware-aware frequency filtering, smoothed
       noise removal, export, rename.
     - Raw Stratagem workflows with hardware files and a full QC pass.
   * - ``publication_ready``
     - Stricter QC, hardware SNR masking, AMT band trimming, stronger
       smoothing, export, rename.
     - Polished Stratagem outputs after the basic workflow has been reviewed.

Run a raw Stratagem convenience preset:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline.stratagem import run_stratagem_preset
   >>> survey = run_stratagem_preset(
   ...     "full_processing",
   ...     edi_dir="2/2EDI",
   ...     coord_file="2.csv",
   ...     raw_dir="raw/2HX",
   ...     outdir="results/stratagem",
   ...     epsg=32649,
   ...     utm_zone="49N",
   ...     rename_basename="T2.",
   ...     overwrite=True,
   ...     verbose=1,
   ... )

Build a Stratagem pipeline object from a normal emtools preset:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline.stratagem import StratagemPipeline
   >>> pipe = StratagemPipeline.from_preset(
   ...     "stratagem_mt",
   ...     coord_file="2.csv",
   ...     raw_dir="raw/2HX",
   ...     epsg=32649,
   ...     utm_zone="49N",
   ...     rename_basename="T2.",
   ... )
   >>> result = pipe.run("2/2EDI", outdir="results/stratagem_mt")

Troubleshooting
---------------

Unknown preset
    Run ``pycsamt pipe presets``.  Preset names are exact and lowercase, for
    example ``basic_qc`` or ``publication_ready``.

I changed ``preset: basic_qc`` but the notch is still 50 Hz
    A config ``preset`` expands the preset first.  Explicit ``steps`` are
    appended; they do not edit existing preset steps.  Generate an expanded
    config and edit the ``NR001`` parameters directly.

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
