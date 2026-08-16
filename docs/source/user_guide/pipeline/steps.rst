.. _pipeline-steps:

Pipeline Steps
==============

Pipeline steps are the atomic processing operations used by
:class:`pycsamt.pipeline.Pipeline`.  Each operation is defined once in the
:term:`step registry` and can then be used from the Python API, YAML/JSON/Python
configuration files, :term:`pipeline preset` recipes, and the ``pycsamt pipe``
command line.

The registry gives every operation a stable :term:`pipeline step code`, such
as ``FREQ001`` or ``NR001``.  The code is intentionally short because it
appears in terminal output, saved reports, plot folders, and reproducible
``pipeline.yaml`` files.  Each code also has a snake-case registry name such
as ``select_band`` or ``notch_powerline``; both forms are accepted when
constructing a step.

Use this page when you need to:

* choose the right operation for a processing chain;
* understand the order in which steps should usually run;
* discover default parameters and diagnostic plots;
* write a custom pipeline configuration;
* extend the registry in a controlled way for new processing operations.

Step Contract
-------------

A registered step has two layers.

``StepSpec``
    :term:`StepSpec` registry metadata.  It stores the code, name, human
    label, category, default parameters, transform function,
    :term:`QC plot function` entries, and whether the operation returns a
    modified :term:`site collection`.

``Step``
    A configured instance used inside a pipeline.  It binds a ``StepSpec`` to
    user :term:`parameter override` values.

The relationship is deterministic.  If a registry entry has defaults
:math:`\theta_{\mathrm{default}}` and the user supplies overrides
:math:`\theta_{\mathrm{user}}`, the configured parameters are

.. math::

   \theta_{\mathrm{step}}
   = \theta_{\mathrm{default}} \cup \theta_{\mathrm{user}},

where the right-hand values win for keys that appear in both dictionaries.
This is why a short configuration can override only ``mains_hz`` for
``NR001`` while keeping the default ``n_harm`` and ``tol_hz``.

In Python, a step is created with :class:`pycsamt.pipeline.Step`:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Step
   >>> notch = Step("NR001", mains_hz=50, n_harm=30)
   >>> notch.params
   {'mains_hz': 50, 'n_harm': 30, 'tol_hz': 0.08}
   >>> band = Step("select_band", band_hz=(0.01, 10000.0))
   >>> band.spec.code
   'FREQ001'

The two lookup styles are equivalent:

.. code-block:: pycon
   :linenos:

   >>> Step("NR001").spec is Step("notch_powerline").spec
   True

The first form is better for compact configs and reports.  The second form is
often easier to read in Python code.

How Steps Execute
-----------------

During a run, the pipeline applies each step in order.  If the input state is
:math:`S_{j-1}` and the configured step is :math:`T_j`, then a normal
:term:`transform step` computes

.. math::

   S_j = T_j(S_{j-1}; \theta_j).

A :term:`diagnostic step` instead records figures or summaries while passing
the same site collection forward:

.. math::

   S_j = S_{j-1}.

Operationally, each step follows this sequence:

1. count the input sites;
2. call the step transform function;
3. generate QC plots, if enabled;
4. save intermediate EDI snapshots, if configured;
5. record a :class:`pycsamt.pipeline.StepResult`;
6. pass the output sites to the next step.

For normal processing steps, ``returns_sites=True`` and the transform output
becomes the input to the next step.  Diagnostic-only steps use
``returns_sites=False``.  They keep the site collection unchanged and are used
to generate plots or summaries at useful checkpoints.

Example result fields:

.. code-block:: pycon
   :linenos:

   >>> result = pipe.run(sites, outdir="results/basic_qc")
   >>> [(sr.step_code, sr.step_name, sr.ok, sr.n_sites_in, sr.n_sites_out)
   ...  for sr in result.step_results]
   [('NR001', 'notch', True, 3, 3), ..., ('QC001', 'qc_snapshot', True, 3, 3)]

Discover Steps From The CLI
---------------------------

The command-line catalogue is the fastest way to inspect your installed
version:

.. code-block:: console
   :linenos:

   pycsamt pipe steps
   pycsamt pipe steps --category frequency
   pycsamt pipe steps --category noise_removal --codes-only
   pycsamt pipe steps --info NR001
   pycsamt pipe steps --info notch_powerline

Machine-readable output is useful for notebooks, dashboards, or scripts:

.. code-block:: console
   :linenos:

   pycsamt pipe steps --format json
   pycsamt pipe steps --category static_shift --format csv

Captured frequency catalogue excerpt:

.. code-block:: json
   :linenos:

   {
     "frequency": [
       {
         "code": "FREQ001",
         "name": "select_band",
         "label": "Frequency Band Select",
         "defaults": {"band_hz": [0.001, 10000.0]},
         "returns_sites": true
       },
       {
         "code": "FREQ002",
         "name": "drop_duplicates",
         "label": "Drop Duplicate Frequencies",
         "defaults": {},
         "returns_sites": true
       }
     ]
   }

Captured ``NR001`` information, simplified to ASCII:

.. code-block:: text
   :linenos:

   NR001  Power-line Harmonic Notch  [noise_removal]
   name     : notch_powerline
   function : pycsamt.emtools.remove_noise.notch_powerline
   defaults : mains_hz=50  n_harm=30  tol_hz=0.08
   qc plots : nr_qc_harmonic_waterfall, nr_qc_snr_gain_profile
   returns  : Sites (transform)
   Usage:
     Step("NR001")
     Step("NR001", mains_hz=50)
     Step("notch_powerline")

Discover Steps From Python
--------------------------

The same registry is available from Python:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import categories, list_steps, lookup_step
   >>> categories()
   ['dimensionality', 'export', 'frequency', 'noise_removal', 'preview', 'qc', 'skew', 'source_effects', 'static_shift', 'tensor']
   >>> {c: len(list_steps(c)) for c in categories()}
   {'dimensionality': 3, 'export': 3, 'frequency': 9, 'noise_removal': 14, 'preview': 2, 'qc': 7, 'skew': 4, 'source_effects': 2, 'static_shift': 4, 'tensor': 7}
   >>> [(spec.code, spec.name, spec.defaults)
   ...  for spec in list_steps("frequency")[:2]]
   [('FREQ001', 'select_band', {'band_hz': (0.001, 10000.0)}), ('FREQ002', 'drop_duplicates', {})]
   >>> lookup_step("NR001") is lookup_step("notch_powerline")
   True

Use :meth:`pycsamt.pipeline.Pipeline.step_info` when you want a formatted
human-readable explanation:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> "Power-line Harmonic Notch" in Pipeline.step_info("NR001")
   True
   >>> "SS001" in Pipeline.catalogue("static_shift")
   True

Registered Categories
---------------------

The current registry contains 55 built-in processing and diagnostic steps.
The opt-in AI audit step (``AI001``, see :doc:`extending`) is not counted
here since it is not registered by default.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Category
     - Codes
     - Main purpose
   * - ``frequency``
     - ``FREQ001`` - ``FREQ009``
     - Select, clean, align, regrid, decimate, smooth, mask, and recover
       frequency samples.
   * - ``noise_removal``
     - ``NR001`` - ``NR014``
     - Remove power-line harmonics, smooth noisy responses, suppress
       outliers, apply EMAP-style filters, enforce off-diagonal consistency,
       and run bundled denoising chains.
   * - ``static_shift``
     - ``SS001`` - ``SS004``
     - Estimate and apply static-shift corrections using AMA, LOESS,
       reference-median, or bilateral-filter strategies.
   * - ``tensor``
     - ``TZ001`` - ``TZ007``
     - Rotate, antisymmetrise, clip, balance, invert, and reorient impedance
       tensors.
   * - ``dimensionality``
     - ``DIM001`` - ``DIM003``
     - Classify dimensionality, mask by dimensionality class, and project to a
       2-D representation.
   * - ``skew``
     - ``SK001`` - ``SK004``
     - Mask or select bands using Bahr skewness and related low-skew criteria.
   * - ``source_effects``
     - ``SRC001`` - ``SRC002``
     - Correct near-field source effects and normalize source response.
   * - ``qc``
     - ``QC001`` - ``QC007``
     - Generate diagnostic snapshots without changing the site collection.
       ``QC005`` - ``QC007`` are data-driven ("smart") diagnostics used by
       the method-aware presets -- see :doc:`presets`.
   * - ``preview``
     - ``PRE001`` - ``PRE002``
     - Raw-vs-processed 1-D preview of a deterministic random station
       subset, meant to open (``PRE001``) and close (``PRE002``) a preset.
   * - ``export``
     - ``EXPORT001`` - ``EXPORT003``
     - Write processed sites out to ModEM, Occam2D, or MARE2DEM solver input
       files.

Recommended Ordering
--------------------

The registry does not force a single scientific workflow, because survey
physics and data quality differ.  Still, most survey pipelines are easier to
reason about when :term:`step ordering` follows the data dependency:
operations that define the sampling axis should come before operations that
estimate physical structure from that axis.  A useful high-level order is:

1. frequency cleanup;
2. noise removal;
3. dimensionality and skew gates;
4. static-shift correction;
5. tensor rotation or tensor balancing;
6. source-effect correction, when needed;
7. QC snapshots and final output.

In symbolic form, a reviewed processing chain often looks like

.. math::

   S_{\mathrm{final}} =
   Q\!\left(
   R_{\mathrm{tensor}}\!\left(
   C_{\mathrm{ss}}\!\left(
   G_{\mathrm{dim/skew}}\!\left(
   N_{\mathrm{noise}}\!\left(
   F_{\mathrm{freq}}(S_0)\right)\right)\right)\right)\right),

where :math:`F_{\mathrm{freq}}` defines the usable frequency axis,
:math:`N_{\mathrm{noise}}` suppresses artifacts, :math:`G_{\mathrm{dim/skew}}`
masks intervals that fail the interpretation assumptions, :math:`C_{\mathrm{ss}}`
corrects static shift, :math:`R_{\mathrm{tensor}}` rotates or balances tensor
components, and :math:`Q` records final diagnostics.  This is a guide, not a
law.  For example, ``stratagem_mt`` runs static-shift correction before AMT
band trimming because that workflow expects the complete Stratagem frequency
range for the static-shift estimate.

A compact first-pass chain looks like this:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --steps FREQ002,FREQ001,FREQ004,NR001,QC001 \
       --out results/basic_manual

The equivalent Python pipeline is:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step
   >>> pipe = Pipeline(
   ...     [
   ...         ("drop_duplicates", Step("FREQ002")),
   ...         ("select_band", Step("FREQ001", band_hz=(0.001, 10000.0))),
   ...         ("align_grid", Step("FREQ004")),
   ...         ("notch", Step("NR001", mains_hz=50)),
   ...         ("qc", Step("QC001")),
   ...     ],
   ...     name="basic_manual",
   ... )
   >>> len(pipe)
   5

Ordering matters because each step receives the site collection produced by
the previous step.  For example, it is usually better to drop duplicate
frequencies before aligning the grid, and it is usually better to perform
frequency and noise cleanup before estimating dimensionality or skew masks.

Frequency Steps
---------------

Frequency steps define the usable sampling axis before more interpretive
processing begins.

.. list-table::
   :header-rows: 1
   :widths: 13 23 32 32

   * - Code
     - Name
     - Use when
     - Defaults
   * - ``FREQ001``
     - ``select_band``
     - You need to keep a scientific frequency interval and remove very low
       or very high frequencies from the chain.
     - ``band_hz=(0.001, 10000.0)``
   * - ``FREQ002``
     - ``drop_duplicates``
     - Multiple samples share the same nominal frequency and should be
       collapsed before alignment or interpolation.
     - none
   * - ``FREQ003``
     - ``drop_low_confidence``
     - Low-confidence frequencies should be removed entirely.
     - none
   * - ``FREQ004``
     - ``align_grid``
     - Sites need a common frequency grid before profile operations or
       cross-site comparisons.
     - none
   * - ``FREQ005``
     - ``regrid_logspace``
     - You want a regular log-frequency sampling density.
     - ``n_per_decade=6``
   * - ``FREQ006``
     - ``decimate``
     - You want a lighter frequency set for quick inspection or faster
       downstream runs.
     - ``step=2``
   * - ``FREQ007``
     - ``smooth_freq``
     - Small frequency-to-frequency oscillations should be smoothed.
     - ``window=3``
   * - ``FREQ008``
     - ``mask_low_confidence``
     - Low-confidence samples should remain in the grid as ``NaN`` values
       rather than being removed.
     - ``method='composite'``, ``threshold=0.5``
   * - ``FREQ009``
     - ``recover_low_confidence``
     - Low-confidence samples should be recovered by interpolation after
       confidence gating.
     - ``method='composite'``, ``ci_hi=0.9``, ``ci_lo=0.5``,
       ``interpolation='linear'``

Common frequency cleanup example:

.. code-block:: yaml
   :linenos:

   steps:
     - name: drop_duplicates
       code: FREQ002
     - name: select_band
       code: FREQ001
       params:
         band_hz: [0.01, 10000.0]
     - name: align_grid
       code: FREQ004

Noise Removal Steps
-------------------

Noise-removal steps suppress survey artifacts while preserving the
interpretable impedance structure.  Prefer targeted operations first; reserve
the bundled denoising step for quick workflows or carefully reviewed
production recipes.

.. list-table::
   :header-rows: 1
   :widths: 13 25 38 24

   * - Code
     - Name
     - Use when
     - Defaults
   * - ``NR001``
     - ``notch_powerline``
     - Power-line harmonics contaminate the response.
     - ``mains_hz=50``, ``n_harm=30``, ``tol_hz=0.08``.  ``mains_hz="auto"``
       is also accepted -- see below.
   * - ``NR002``
     - ``smooth_logfreq``
     - Responses are noisy along log-frequency.
     - none
   * - ``NR003``
     - ``shrink_group_trend``
     - Isolated outliers should be pulled toward the group trend rather than
       removed aggressively.
     - none
   * - ``NR004``
     - ``hampel_filter``
     - Spike-like outliers need robust local filtering.
     - ``win=3``, ``nsig=3.0``
   * - ``NR005``
     - ``spatial_median``
     - Neighboring stations define a stable spatial trend.
     - none
   * - ``NR006``
     - ``emap_filter``
     - You want EMAP-style spatial filtering.
     - none
   * - ``NR007``
     - ``emap_confidence``
     - EMAP filtering should be gated by confidence and return the filtered
       site collection.
     - none
   * - ``NR008``
     - ``rpca_offdiag``
     - Off-diagonal components need robust low-rank denoising.
     - none
   * - ``NR009``
     - ``enforce_offdiag``
     - Off-diagonal impedance components should be made more internally
       consistent.
     - none
   * - ``NR010``
     - ``mask_incoherent``
     - Incoherent frequency samples should be masked before interpretation.
     - none
   * - ``NR011``
     - ``fixed_length_mavg``
     - A fixed moving-average smoother is appropriate for all or selected
       components.
     - ``window=5``, ``component='all'``
   * - ``NR012``
     - ``trimmed_mavg``
     - You want a moving average that is less sensitive to extremes.
     - ``window=5``, ``component='all'``
   * - ``NR013``
     - ``correct_static_shift_spatial``
     - A spatial-window static-shift correction is needed from the noise
       module.
     - ``window_m=1500.0``, ``spacing_m=200.0``, ``comp='det'``
   * - ``NR014``
     - ``denoise_pipeline``
     - A bundled denoising recipe is acceptable for a first pass or a
       reviewed production preset.
     - ``mains_hz=50.0``, ``n_harm=30``, ``tol_hz=0.08``,
       ``smooth_win=5``, ``gate_snr=2.5``

Example:

.. code-block:: yaml
   :linenos:

   steps:
     - name: notch
       code: NR001
       params:
         mains_hz: 60
         n_harm: 20
     - name: robust_spikes
       code: NR004
       params:
         win: 5
         nsig: 3.5
     - name: incoherence_gate
       code: NR010

``mains_hz="auto"`` (``NR001`` only)
    Real EDI frequency grids are usually log-spaced rather than sampled
    exactly on multiples of 50/60 Hz, so a fixed numeric ``mains_hz`` can
    silently miss harmonics that fall outside the tight ``+-tol_hz``
    window.  Passing the literal string ``"auto"`` instead makes the
    match data-driven: pyCSAMT scores 50 Hz and 60 Hz -- the only two
    real-world AC grid frequencies -- against the survey's actual pooled
    frequency array, resolves to whichever explains more harmonics, then
    snaps each harmonic to its single nearest real sample within a
    relative ``snap_frac`` tolerance (default 1%, deliberately tight so a
    coarse, non-mains-aware grid does not get treated as evidence).  If
    neither candidate clears a minimum-evidence bar, ``"auto"`` leaves the
    data untouched and warns rather than guessing:

    .. code-block:: pycon
       :linenos:

       >>> from pycsamt.api import read_edis
       >>> from pycsamt.emtools.remove_noise import notch_powerline
       >>> survey = read_edis("data/AMT/WILLY_DATA/L18PLT", recursive=False, strict=False, progress=False)
       >>> out = notch_powerline(survey.collection, mains_hz="auto", verbose=1)

    A ``UserWarning`` reports the resolved fundamental (captured against
    the real 28-station ``L18PLT`` line):

    .. code-block:: text
       :linenos:

       notch_powerline: mains_hz="auto" resolved to 60 Hz (4/30 harmonics matched within 1%) for this survey.

    Passing a plain number is completely unaffected -- output is
    identical to every previous release. ``verbose=2`` additionally logs
    every individual snapped harmonic (station, target frequency, the
    real frequency substituted for it, and the offset). In a config file:

    .. code-block:: yaml
       :linenos:

       - name: notch
         code: NR001
         params:
           mains_hz: auto

Static-Shift Steps
------------------

Static-shift correction should normally run after basic frequency and noise
cleanup.  Static-shift steps modify the site collection and generate
comparison plots when output plotting is enabled.

.. list-table::
   :header-rows: 1
   :widths: 13 26 42 19

   * - Code
     - Name
     - Use when
     - Defaults
   * - ``SS001``
     - ``correct_ss_ama``
     - You want the standard AMA static-shift correction.
     - none
   * - ``SS002``
     - ``correct_ss_loess``
     - You want a two-phase LOESS estimate followed by factor application.
     - none
   * - ``SS003``
     - ``correct_ss_refmedian``
     - A reference median across stations is the preferred shift estimator.
     - none
   * - ``SS004``
     - ``correct_ss_bilateral``
     - Shift factors should be estimated with a bilateral-style spatial
       filter.
     - ``half_window=4``, ``max_skew=6.0``, ``summary='median'``

Example:

.. code-block:: yaml
   :linenos:

   steps:
     - name: select_band
       code: FREQ001
       params:
         band_hz: [0.01, 10000.0]
     - name: notch
       code: NR001
     - name: static_shift
       code: SS004
       params:
         half_window: 5
         max_skew: 5.0

Tensor Steps
------------

Tensor steps change tensor orientation, symmetry, or component consistency.
Run them only after the frequency axis and gross noise problems are under
control.

.. list-table::
   :header-rows: 1
   :widths: 13 24 42 21

   * - Code
     - Name
     - Use when
     - Defaults
   * - ``TZ001``
     - ``rotate_strike``
     - The impedance tensor should be rotated to estimated strike.
     - ``method='swift'``
   * - ``TZ002``
     - ``antisymmetrize``
     - Off-diagonal tensor symmetry needs to be enforced or inspected.
     - ``how='rms'``
   * - ``TZ003``
     - ``sigma_clip``
     - Tensor outlier frequencies should be sigma-clipped.
     - ``sigma=3``
   * - ``TZ004``
     - ``balance_offdiag``
     - Off-diagonal components need balancing.
     - none
   * - ``TZ005``
     - ``rotate_fixed``
     - A known fixed rotation angle should be applied.
     - ``angle=0.0``
   * - ``TZ006``
     - ``invert_tensor``
     - The impedance tensor should be inverted.
     - none
   * - ``TZ007``
     - ``orient_from_sensors``
     - Sensor orientations are known and should be corrected explicitly.
     - ``ex=0.0``, ``ey=0.0``, ``bx=0.0``, ``by=0.0``,
       ``degrees=True``

Example:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step
   >>> pipe = Pipeline(
   ...     [
   ...         ("drop_duplicates", Step("FREQ002")),
   ...         ("select_band", Step("FREQ001", band_hz=(1.0, 10000.0))),
   ...         ("static_shift", Step("SS001")),
   ...         ("rotate", Step("TZ001", method="swift")),
   ...         ("qc", Step("QC001")),
   ...     ],
   ...     name="tensor_ready",
   ... )
   >>> [label for label, step in pipe]
   ['drop_duplicates', 'select_band', 'static_shift', 'rotate', 'qc']

Dimensionality Steps
--------------------

Dimensionality steps help decide whether a frequency band or station interval
is suitable for 1-D, 2-D, or 3-D assumptions.

.. list-table::
   :header-rows: 1
   :widths: 13 24 43 20

   * - Code
     - Name
     - Use when
     - Defaults
   * - ``DIM001``
     - ``classify_dim``
     - You need a diagnostic classification of dimensionality.  This is a
       diagnostic step and does not modify sites.
     - none
   * - ``DIM002``
     - ``mask_by_dim``
     - You want to keep or mask samples by dimensionality class before
       further processing.
     - ``dim_type='2d'``
   * - ``DIM003``
     - ``project_2d``
     - The response should be projected toward a 2-D representation.
     - none

Example:

.. code-block:: yaml
   :linenos:

   steps:
     - name: classify_dimensionality
       code: DIM001
     - name: keep_2d
       code: DIM002
       params:
         dim_type: 2d
     - name: project_2d
       code: DIM003

Skew Steps
----------

Skew steps are useful when selecting lower-skew intervals for 2-D inversion
or for rejecting intervals where 3-D effects dominate.

.. list-table::
   :header-rows: 1
   :widths: 13 25 42 20

   * - Code
     - Name
     - Use when
     - Defaults
   * - ``SK001``
     - ``mask_by_skew``
     - Frequencies above a Bahr skewness threshold should be masked.
     - ``thresh=0.3``
   * - ``SK002``
     - ``longest_low_skew``
     - You want the longest continuous low-skew band.
     - ``thresh=0.3``
   * - ``SK003``
     - ``select_skew_band``
     - You want the lowest-skew band selected automatically.
     - none
   * - ``SK004``
     - ``close_skew_gaps``
     - Small gaps in a low-skew interval should be closed.
     - ``thresh=3.0``, ``max_gap=1``

Example:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --steps FREQ002,FREQ001,FREQ004,NR001,SK001,TZ001,QC001 \
       --out results/low_skew_rotation

Source-Effect Steps
-------------------

Source-effect steps are more specialized than the basic QC and denoising
chain.  Use them when survey geometry, transmitter effects, or diagnostic
plots indicate source overprint or near-field behavior.

.. list-table::
   :header-rows: 1
   :widths: 13 28 44 15

   * - Code
     - Name
     - Use when
     - Defaults
   * - ``SRC001``
     - ``correct_near_field``
     - Near-field source effects should be corrected.
     - none
   * - ``SRC002``
     - ``normalize_response``
     - Source response normalization is needed before interpretation or
       comparison.
     - none

Example:

.. code-block:: yaml
   :linenos:

   steps:
     - name: source_diagnostic
       code: QC003
     - name: near_field
       code: SRC001
     - name: normalized_source
       code: SRC002

QC And Diagnostic Steps
-----------------------

QC steps are diagnostic-only checkpoints.  They generate figures but pass the
site collection through unchanged.

.. list-table::
   :header-rows: 1
   :widths: 13 28 44 15

   * - Code
     - Name
     - Use when
     - Defaults
   * - ``QC001``
     - ``qc_snapshot``
     - You want a quick-look dashboard, station confidence view, and coverage
       pseudo-section.
     - none
   * - ``QC002``
     - ``field_zone_snapshot``
     - You want a field-zone classification plot.
     - none
   * - ``QC003``
     - ``source_overprint_snapshot``
     - You want a source-overprint diagnostic before or after source-effect
       corrections.
     - none
   * - ``QC004``
     - ``depth_section_snapshot``
     - You want a Bostick depth-section diagnostic snapshot.
     - none

Use QC snapshots at major checkpoints:

.. code-block:: yaml
   :linenos:

   steps:
     - name: raw_qc
       code: QC001
     - name: notch
       code: NR001
     - name: cleaned_qc
       code: QC001
     - name: static_shift
       code: SS001
     - name: final_qc
       code: QC001

Export Steps
------------

Export steps write processed sites out to solver-ready input files.  Unlike
QC steps they are not diagnostic-only: ``returns_sites=True``, so a real
write failure propagates through the pipeline's normal ``on_step_error``
policy rather than being silently turned into a warning.

.. list-table::
   :header-rows: 1
   :widths: 13 22 45 20

   * - Code
     - Name
     - Use when
     - Defaults
   * - ``EXPORT001``
     - ``export_modem``
     - You want ModEM data, model, and control files (plus covariance for a
       3-D run).
     - ``workdir='pipeline_exports/modem'``
   * - ``EXPORT002``
     - ``export_occam2d``
     - You want a full native Occam2D data/mesh/model/startup working
       directory.
     - ``workdir='pipeline_exports/occam2d'``
   * - ``EXPORT003``
     - ``export_mare2dem``
     - You want the MARE2DEM ``.emdata`` data file specifically.
     - ``workdir='pipeline_exports/mare2dem'``,
       ``data_filename='line.emdata'``

``EXPORT003`` deliberately writes only the data file.
``pycsamt.models.mare2dem.builder.InputBuilder.write_resistivity``/``.write_settings``
are separate calls with their own configuration, and full mesh generation
(the ``.poly`` file) is not wired into any pipeline step yet.

A step function only ever receives ``sites`` and its own parameters, with no
hook into the run's own output directory, so each export step's ``workdir``
resolves relative to the process's current working directory rather than
following ``--out``.  Point it explicitly at wherever the run's output
should live:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step
   >>> pipe = Pipeline(
   ...     [
   ...         ("notch", Step("NR001")),
   ...         (
   ...             "export_modem",
   ...             Step("EXPORT001", workdir="results/line22_export/modem"),
   ...         ),
   ...     ],
   ...     name="export_ready",
   ... )
   >>> [label for label, step in pipe]
   ['notch', 'export_modem']

Captured against the real 25-station WILLY L22 line:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/AMT/WILLY_DATA/L22PLT \
       --steps NR001,EXPORT001,EXPORT002,EXPORT003 \
       --out results/line22_export

.. code-block:: text
   :linenos:

   PipelineResult  'cli_pipeline'
     Sites   : 25 in → 25 out
     Steps   : 4 (4 ok, 0 err)
     Time    : 12.15 s
     Plots   : 2
     Output  : results\line22_export

Files actually written under each step's default ``workdir``:

.. code-block:: text
   :linenos:

   pipeline_exports/modem/     control.inv  covariance.cov  data.dat  m0.ws
   pipeline_exports/occam2d/   Occam2DMesh  Occam2DModel  OccamDataFile.dat  Startup
   pipeline_exports/mare2dem/  line.emdata

Parameters And Defaults
-----------------------

User parameters are merged on top of registry defaults.  If a default is not
overridden, the registry value is used.  The merge is shallow and
right-biased:

.. math::

   \theta_{\mathrm{active}}[k] =
   \begin{cases}
   \theta_{\mathrm{user}}[k], & k \in \theta_{\mathrm{user}},\\
   \theta_{\mathrm{default}}[k], & k \notin \theta_{\mathrm{user}}.
   \end{cases}

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Step
   >>> step = Step("NR001", mains_hz=60)
   >>> step.params["mains_hz"]
   60
   >>> step.params["n_harm"]
   30
   >>> step.params["tol_hz"]
   0.08

In YAML, place overrides under ``params``:

.. code-block:: yaml
   :linenos:

   steps:
     - name: notch_60hz
       code: NR001
       params:
         mains_hz: 60
         n_harm: 25
         tol_hz: 0.1

Keep parameter names aligned with the underlying emtools function.  Use
``pycsamt pipe steps --info CODE`` or ``Pipeline.step_info("CODE")`` before
editing an unfamiliar step.  For reviewed work, write every scientifically
important parameter explicitly even when it equals the default; that makes the
processing choice visible in the configuration file.

Step Labels
-----------

The registry code identifies the operation.  The pipeline label identifies
the position or intent of that operation inside one workflow.

.. code-block:: yaml
   :linenos:

   steps:
     - name: raw_notch
       code: NR001
     - name: post_static_shift_notch
       code: NR001
       params:
         n_harm: 10

The two entries use the same operation, but reports and plot directories can
distinguish them because their labels differ.

Diagnostic Plots
----------------

Each :term:`StepSpec` may declare :term:`QC plot function` entries.  When
plotting is enabled, the pipeline calls those functions after a successful
transform and saves the figures under the run output directory.  The files are
diagnostic artifacts; they do not feed back into the next transform unless a
separate step explicitly uses the underlying data.

The plot call order for step :math:`j` is

.. math::

   S_j = T_j(S_{j-1};\theta_j)
   \quad\Longrightarrow\quad
   \mathrm{plots}_j =
   \{q(S_j): q \in Q_j\},

where :math:`Q_j` is the set of registered QC plot functions for that step.
Diagnostic-only QC steps have :math:`S_j=S_{j-1}`, but their plots still
describe the current state of the survey at that checkpoint.

Examples:

* ``NR001`` can generate a harmonic waterfall and SNR gain profile.
* ``FREQ001`` can generate band microstrips and a coverage-quality heatmap.
* ``SS004`` can generate static-shift delta, comparison, and summary plots.
* ``QC001`` can generate the quick-look dashboard, station confidence
  dashboard, and coverage pseudo-section.

QC plot failures are skipped so that a successful processing step is not
turned into a failed run only because a diagnostic figure could not be drawn.
If a plot is missing, inspect the saved report and rerun with verbose CLI
output.

Error Handling
--------------

Step failures are controlled by the pipeline error policy.  From the CLI, use
``--on-error``:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc --on-error raise
   pycsamt pipe run data/edis --preset basic_qc --on-error warn
   pycsamt pipe run data/edis --preset basic_qc --on-error skip

``raise``
    Stop immediately at the failing step.

``warn``
    Record the error, warn, continue with the previous site collection, and
    return a nonzero CLI exit status if any step failed.

``skip``
    Record the error and continue with the previous site collection without a
    warning.  The final result is still not OK if any step failed.

Diagnostic-only steps are special: if their diagnostic transform raises, the
step warns and passes the original sites through unchanged.

Choosing Between Similar Steps
------------------------------

Use the smallest operation that matches the problem you can actually see in
the data.

For harmonic contamination
    Start with ``NR001`` and tune ``mains_hz``, ``n_harm``, and ``tol_hz``.

For frequency-grid problems
    Start with ``FREQ002`` followed by ``FREQ001`` and ``FREQ004``.

For isolated spikes
    Start with ``NR004``.  Move to ``NR012`` when the smoothing window should
    be robust to extremes across a broader interval.

For spatially coherent noise
    Compare ``NR005``, ``NR006``, and ``NR007`` with QC plots enabled.

For static shift
    Start with ``SS001`` for a standard workflow.  Use ``SS002``, ``SS003``,
    or ``SS004`` when the survey geometry and QC evidence justify a different
    estimator.

For 2-D inversion preparation
    Combine frequency cleanup, noise removal, skew or dimensionality gates,
    static-shift correction, strike rotation, and final QC.

Full Example Configuration
--------------------------

This example is explicit enough for a reviewed processing run while remaining
short enough to audit:

.. code-block:: yaml
   :linenos:

   name: line22_processing
   output_dir: results/line22_processing

   steps:
     - name: drop_duplicates
       code: FREQ002

     - name: select_band
       code: FREQ001
       params:
         band_hz: [0.01, 10000.0]

     - name: align_grid
       code: FREQ004

     - name: notch_powerline
       code: NR001
       params:
         mains_hz: 50
         n_harm: 30
         tol_hz: 0.08

     - name: robust_spike_filter
       code: NR004
       params:
         win: 3
         nsig: 3.0

     - name: low_skew_gate
       code: SK001
       params:
         thresh: 0.3

     - name: static_shift
       code: SS001

     - name: rotate_to_strike
       code: TZ001
       params:
         method: swift

     - name: final_qc
       code: QC001

Run it with:

.. code-block:: console
   :linenos:

   pycsamt pipe show config/line22_processing.yaml
   pycsamt pipe run data/edis \
       --config config/line22_processing.yaml \
       --out results/line22_processing \
       --on-error warn \
       -v

Extension Policy
----------------

The current pipeline registry is source-controlled and static.  That is
intentional: processing steps are scientific operations, so adding one should
be reviewed with tests, defaults, documentation, and diagnostic behavior.

A second, runtime path exists for steps that should not be reviewed into
pyCSAMT itself, such as a site-specific correction or a vendor-format
adapter: ``register_step`` adds a ``StepSpec`` to the same registry without
touching pyCSAMT's source tree, and the resulting entry is usable everywhere
a built-in step is.  See :doc:`extending` for the full mechanism, including
how a third-party package can announce its steps through an entry point so
they are found automatically rather than registered by hand in every script.

To add a new built-in step:

1. implement or expose the transform function in the appropriate
   ``pycsamt.emtools`` module;
2. add a :class:`pycsamt.pipeline.StepSpec` entry in
   ``pycsamt/pipeline/_registry.py``;
3. choose a code prefix that matches the category, such as ``FREQ``, ``NR``,
   ``SS``, ``TZ``, ``DIM``, ``SK``, ``SRC``, or ``QC``;
4. define conservative defaults;
5. attach QC plot functions when useful;
6. decide whether the step returns a modified site collection;
7. add tests that prove lookup, parameter merging, execution, and config
   round-tripping;
8. update this page and any presets that should use the new operation.

Minimal registry entry pattern:

.. code-block:: pycon
   :linenos:

   >>> StepSpec(
   ...     code="NR015",
   ...     name="my_new_filter",
   ...     label="My New Noise Filter",
   ...     category="noise_removal",
   ...     mod="pycsamt.emtools.remove_noise",
   ...     fn_name="my_new_filter",
   ...     defaults={"window": 5},
   ...     qc_defs=[
   ...         ("pycsamt.emtools.remove_noise", "nr_qc_station_offdiag_curves"),
   ...     ],
   ... )

For two-phase operations, use an ``override_fn`` wrapper that receives
``sites`` and keyword parameters, then returns the modified site collection.
Static-shift steps such as ``SS002`` and ``SS003`` follow this pattern.

Testing New Steps
-----------------

At minimum, a new step should be covered by tests like these:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step, lookup_step
   >>>
   >>> def test_registry_lookup_for_new_step():
   ...     spec = lookup_step("NR015")
   ...     assert spec.name == "my_new_filter"
   ...     assert lookup_step("my_new_filter") is spec
   ...
   >>> def test_step_defaults_and_overrides():
   ...     step = Step("NR015", window=7)
   ...     assert step.params["window"] == 7
   ...
   >>> def test_pipeline_accepts_new_step(sites):
   ...     pipe = Pipeline([("new_filter", Step("NR015"))])
   ...     result = pipe.run(sites, outdir=None, save_plots=False)
   ...     assert result.ok

Also test YAML loading if the step is expected to be used from configuration
files.

Troubleshooting
---------------

Unknown step code
    Run ``pycsamt pipe steps`` or ``pycsamt pipe steps --info CODE``.  Codes
    are exact and uppercase, for example ``NR001``.

Unknown step name
    Use the snake-case registry name, for example ``notch_powerline`` rather
    than ``notch``.  Pipeline labels may be arbitrary, but registry names are
    fixed.

A parameter is ignored or raises ``TypeError``
    Check the underlying emtools function signature and the output of
    ``pycsamt pipe steps --info CODE``.  Parameter names must match the
    transform function.

The site count changes unexpectedly
    Frequency-selection, masking, or dimensionality operations can change the
    effective data content.  Inspect ``StepResult.n_sites_in`` and
    ``StepResult.n_sites_out`` in the run summary.

Expected plots are missing
    Confirm that plots were not disabled with ``--no-plots``.  Some QC plots
    are skipped when the required data are absent or when the plotting
    function cannot build a figure for the current site collection.

A diagnostic step appears to do nothing
    QC steps intentionally pass sites through unchanged.  Inspect the plot
    directory and report rather than expecting data changes.

Related Pages
-------------

* :doc:`concepts` explains how pipelines execute and produce
  ``PipelineResult`` objects.
* :doc:`configuration_files` explains how to write step lists in YAML, JSON,
  and Python configs.
* :doc:`cli_pipe` explains ``pycsamt pipe steps``, ``pycsamt pipe show``, and
  ``pycsamt pipe run``.
* :doc:`presets` explains built-in step sequences.
* :doc:`extending` explains ``register_step`` and third-party plugin steps.
* :doc:`caching` explains caching a step's output and how that gives an
  interrupted run resume.
* :doc:`outputs` explains reports, plots, processed EDIs, and saved pipeline
  snapshots.
