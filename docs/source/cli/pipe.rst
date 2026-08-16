Pipeline Commands
=================

``pycsamt pipe`` runs, inspects, and scaffolds pyCSAMT processing
pipelines. A pipeline is an ordered list of processing steps applied to a
``Sites`` collection loaded from EDI data.

Use this command group when you want a reproducible survey workflow:
noise removal, frequency filtering, static-shift correction, tensor
analysis, QC plots, processed EDI outputs, and reports.

The CLI uses the same ``pycsamt.pipeline`` engine as the Python API. A
workflow can therefore move between the terminal, a notebook, and a saved
YAML/JSON/Python config file without changing meaning.

Command map
-----------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Command
     - Purpose
   * - ``pycsamt pipe steps``
     - Browse available processing step codes.
   * - ``pycsamt pipe presets``
     - List or expand named workflow presets.
   * - ``pycsamt pipe init``
     - Generate a ready-to-edit pipeline config.
   * - ``pycsamt pipe show``
     - Pretty-print a config file or preset.
   * - ``pycsamt pipe run``
     - Execute a config, preset, or ad-hoc step list.
   * - ``pycsamt pipe history``
     - List runs previously logged with ``pipe run --history``.
   * - ``pycsamt pipe plugins``
     - Discover and list third-party pipeline-step plugins.

``pycsamt pipe --with-plugins`` and ``pycsamt pipe --with-ai-steps`` are
group-level flags (placed before the subcommand) rather than subcommands of
their own -- see `Plugins and the opt-in AI step`_.

Core model
----------

The pipeline command group revolves around four public concepts:

``Pipeline``
    The ordered workflow. It may be built from a preset, a config file, or a
    short ad-hoc step list.

``Step``
    A configured processing operation. Registry codes such as ``NR001`` and
    registry names such as ``notch_powerline`` identify the operation.

``Preset``
    A named baseline workflow for common processing intentions.

``PipelineResult``
    The run record. It stores per-step status, timing, input/output site
    counts, plot paths, processed EDI paths, output directory, errors, and
    (when caching is enabled) whether each step was replayed from cache.

During ``pipe run``, pyCSAMT loads EDI data into ``Sites``, applies each
configured step in order, writes optional outputs, and prints a terminal
summary. Four opt-in flags extend this without changing default behaviour:
``--cache`` (resume from the step cache), ``--live`` (live per-step status
table), ``--history`` (log the run), and ``--dashboard`` (a richer branded
HTML report) -- see `Cache, live progress, history, and the dashboard`_.

Input resolution
----------------

Pipeline runs resolve the EDI source in this priority order:

1. Positional ``EDI_DIR``.
2. ``--survey DIR``.
3. Active survey context set by ``pycsamt survey set``.

Examples:

.. code-block:: console

   pycsamt pipe run ./edis/ --preset basic_qc
   pycsamt pipe run --survey ./edis/ --preset basic_qc
   pycsamt survey set ./edis/
   pycsamt pipe run --preset basic_qc

Use ``--fresh`` to bypass cached survey context when resolving sites.

Explore steps
-------------

.. code-block:: console

   pycsamt pipe steps
   pycsamt pipe steps --category noise_removal
   pycsamt pipe steps --category static_shift --codes-only
   pycsamt pipe steps --info NR001
   pycsamt pipe steps --info notch_powerline
   pycsamt pipe steps --format json

The registry currently holds 55 built-in steps, grouped into the
categories ``frequency``, ``noise_removal``, ``static_shift``, ``tensor``,
``dimensionality``, ``skew``, ``source_effects``, ``qc``, ``preview``, and
``export``. ``preview`` (``PRE001``/``PRE002``) and three of the ``qc``
codes (``QC005``-``QC007``) are the data-driven "smart" QC and raw/
processed-preview steps behind the method-aware presets -- see
`Explore presets`_.

Useful details:

``--info CODE_OR_NAME``
    Shows the step code, registry name, display label, category, transform
    function, default parameters, attached QC plot functions, and whether the
    step returns a modified ``Sites`` collection.

``--codes-only``
    Prints one code per line. With ``--format json`` it prints a JSON list of
    codes.

``--format csv``
    For the catalogue, emits
    ``code,name,label,category,returns_sites``.

Step codes are best for config files and reports. Registry names are useful
while exploring interactively, and are also accepted by ``--info`` and
``--steps``.

Explore presets
---------------

.. code-block:: console

   pycsamt pipe presets
   pycsamt pipe presets --expand full_processing
   pycsamt pipe presets --expand basic_qc --format json

Presets are named, opinionated workflows. Eleven are registered today:

- **Chosen by processing intent**: ``basic_qc``, ``noise_reduction``,
  ``full_processing``, ``tensor_analysis``, ``dimensionality_filter``,
  ``publication_ready``, and ``stratagem_mt`` (for Stratagem AMT data
  already loaded as ``Sites``).
- **Chosen by EM survey method**: ``mt_qc``, ``amt_qc``, ``csamt_qc``, and
  ``csumt_qc``. Each adds real near-field/transition-zone correction
  (CSAMT/CSUMT only) and data-driven QC -- phase-tensor and tipper plots
  that check the actual data before drawing themselves, instead of always
  producing the same figures regardless of what the survey contains. In
  Python, :func:`pycsamt.pipeline.get_preset_for_method` maps an explicit
  method string (``"MT"``, ``"AMT"``, ``"CSAMT"``, ``"CSUMT"``, ...) to the
  matching preset name rather than hand-building it. Full mechanism in
  :doc:`../user_guide/pipeline/presets`'s "Method-Aware Presets" section.

Use ``--expand`` to see the exact step sequence.

Machine-readable output:

``pycsamt pipe presets --format json``
    Emits a list of preset objects with ``name``, ``description``,
    ``n_steps``, and ``codes``.

``pycsamt pipe presets --expand NAME --format json``
    Emits one preset object with the ordered ``steps`` list. Each step entry
    includes the label, code, registry name, category, and parameters.

``pycsamt pipe presets --format csv``
    Emits ``name,description,n_steps,codes``.

Plugins and the opt-in AI step
------------------------------

Two families of steps are never registered by default, because discovering
them has a real cost that most invocations shouldn't pay:

.. code-block:: console

   pycsamt pipe plugins
   pycsamt pipe plugins --strict
   pycsamt pipe plugins --format json
   pycsamt pipe --with-plugins run --steps MY_PLUGIN001,NR001 --survey ./edis/
   pycsamt pipe --with-ai-steps run --steps AI001 --survey ./edis/

``pycsamt pipe plugins``
    Scans the ``pycsamt.pipeline.steps`` entry-point group (declared by a
    third-party package in its own ``pyproject.toml``), calls every
    registration callable found there, and lists every step now registered
    with ``origin=plugin``. This is the one command that always discovers
    plugins -- that is its job. ``--strict`` exits non-zero if any plugin
    failed to load (useful in CI).

``--with-plugins`` (a ``pycsamt pipe`` group flag, placed before the
subcommand)
    Discovers plugins before any other subcommand runs, so a plugin step
    code can be used with ``--steps`` in the same invocation. Off by
    default: scanning installed packages for entry points can take several
    seconds in a large environment.

``--with-ai-steps`` (same group-flag placement)
    Registers the opt-in AI domain-gap survey-audit step (``AI001``) before
    the subcommand runs. Off by default: resolving it imports ``torch``, a
    real one-time cost this CLI does not force on users who never touch AI
    steps.

Both flags are equivalent to setting an environment variable
(``PYCSAMT_PIPELINE_LOAD_PLUGINS=1`` / ``PYCSAMT_PIPELINE_LOAD_AI_STEPS=1``),
useful for a CI job or a shell profile that always wants them on.

Create a config
---------------

.. code-block:: console

   pycsamt pipe init
   pycsamt pipe init --preset full_processing --name willy_survey
   pycsamt pipe init --format py --preset basic_qc -o config/willy.py
   pycsamt pipe init --format json --print
   pycsamt pipe init --name l22_profile --outdir results/ -o config/

Supported config formats are ``yaml``, ``json``, and ``py``. A generated
config contains active steps from the selected preset plus comments or
structure for other available steps.

Path and naming rules:

``--name TEXT``
    Sets the pipeline name written into the config. The default is
    ``my_workflow``.

``--outdir DIR``
    Writes the default pipeline output directory into the config. A later
    ``pycsamt pipe run --out DIR`` overrides this value.

``-o, --output FILE_OR_DIR``
    Writes the generated config to a path. If the path is a directory, the
    filename is derived from ``--name`` and ``--format``. If omitted, the file
    is written as ``<name>.<format>`` in the current directory.

``--print``
    Prints the generated config to stdout and writes nothing.

Config schema
-------------

YAML and JSON configs use the same schema:

.. code-block:: yaml

   name: line22_basic_qc
   output_dir: results/line22_basic_qc
   steps:
     - name: notch
       code: NR001
       params:
         mains_hz: 50
         n_harm: 30
     - name: select_band
       code: FREQ001
       params:
         band_hz: [0.001, 10000.0]
     - name: qc_snapshot
       code: QC001

``name`` and ``output_dir`` are optional. Each step entry must provide a
``code`` or registry ``name``. The pipeline label is taken from the entry
``name`` when present; otherwise pyCSAMT uses the registry name.

``NR001``'s ``mains_hz`` also accepts the literal string ``auto`` --
pyCSAMT then detects 50 vs 60 Hz from the survey's own frequency grid and
snaps each harmonic to the nearest real sample instead of requiring an
exact match. See :doc:`../user_guide/pipeline/steps`.

A config may also include ``preset: basic_qc``. In that case, preset steps are
loaded first and the explicit ``steps`` entries are appended.

Python configs expose the same structure as a module-level dictionary:

.. code-block:: python

   pipeline_config = dict(
       name="line22_basic_qc",
       output_dir="results/line22_basic_qc",
       steps=[
           dict(name="notch", code="NR001", params=dict(mains_hz=50)),
           dict(
               name="select_band",
               code="FREQ001",
               params=dict(band_hz=[0.001, 10000.0]),
           ),
           dict(name="qc_snapshot", code="QC001"),
       ],
   )

Show a config or preset
-----------------------

.. code-block:: console

   pycsamt pipe show workflow.yaml
   pycsamt pipe show --preset publication_ready
   pycsamt pipe show workflow.yaml --n-steps 3
   pycsamt pipe show --preset full_processing --format json

``show`` is the safest way to verify step order and parameters before a
run. It accepts the same slicing controls as ``run``:
``--from-step``, ``--until-step``, and ``--n-steps``.

``show`` requires either a ``CONFIG_FILE`` argument or ``--preset NAME``.
Slicing matches pipeline labels, registry codes, or registry names. For
``--format json``, the output includes the pipeline ``name``, ``n_steps``, and
an ordered ``steps`` list with labels, codes, registry names, categories, long
labels, and parameters. For ``--format csv``, the header is
``idx,label,code,name,category,params``.

Run a pipeline
--------------

A pipeline definition is resolved in this priority order:

1. ``--config FILE``.
2. ``--preset NAME``.
3. ``--steps CODE,CODE,...``.

``--config`` accepts ``.yaml``, ``.yml``, ``.json``, and ``.py`` files.
If ``--config`` is present, the config file is the source of truth.
``--preset`` and ``--steps`` are ignored by resolution.

Examples:

.. code-block:: console

   pycsamt pipe run --preset full_processing
   pycsamt pipe run --config workflow.yaml --survey ./data/AMT/ --out results/
   pycsamt pipe run --steps NR001,FREQ002,FREQ001,FREQ004,SS001 --survey ./edis/
   pycsamt pipe run --config workflow.yaml --dry-run
   pycsamt pipe run --config workflow.yaml --n-steps 3 --dry-run
   pycsamt pipe run --config workflow.yaml --from-step align --until-step correct_ss
   pycsamt pipe run --config workflow.yaml --survey ./edis/ --cache
   pycsamt pipe run --preset basic_qc --live
   pycsamt pipe run --preset basic_qc --history
   pycsamt pipe run --preset basic_qc --dashboard

Execution controls:

``--dry-run``
    Resolve the sites and pipeline but do not process or write files.

``--from-step LABEL_OR_CODE``
    Start at a step label, step code, or step name.

``--until-step LABEL_OR_CODE``
    Stop after a step label, code, or name.

``--n-steps INT``
    Keep only the first ``INT`` steps after other slicing.

``--on-error raise|warn|skip``
    Choose how step failures are handled.

``--jobs INT``
    Accepted as a shared CLI option. Pipeline steps run in order; treat this
    as forward-compatible unless a specific step documents parallel behavior.

Error policy:

``raise``
    Stop immediately and report the exception.

``warn``
    Warn, store the error in the step result, continue with the previous site
    collection, and exit nonzero at the end if any step failed.

``skip``
    Store the error and continue silently with the previous site collection.
    The final command still exits nonzero when the result is not OK.

Output controls:

``--out DIR``
    Root directory for processed EDIs, plots, and reports.

``--no-plots``
    Skip QC figure generation.

``--no-edi``
    Skip processed EDI output.

``--no-report``
    Skip HTML/text run reports (also skips ``dashboard.html`` if
    ``--dashboard`` was passed).

``--dpi INT`` and ``--plot-fmt png|pdf|svg``
    Control saved QC figure resolution and format.

If ``--out`` is omitted, output directory resolution follows this order:

1. the config file ``output_dir`` value, when available;
2. the global pipeline default, ``pipe_results``.

``--dry-run`` prints the resolved pipeline and site count but writes nothing.
Pass ``--out`` explicitly when the dry-run display should show the intended
target directory.

Cache, live progress, history, and the dashboard
------------------------------------------------

Four more flags on ``pipe run`` are opt-in and off by default -- a run that
doesn't pass them behaves exactly as before:

``--cache`` / ``--cache-dir DIR``
    Keys each step's output by the exact upstream data, step code, and
    parameters, in a content-addressed store (default
    ``~/.pycsamt/pipeline_cache``, override with ``--cache-dir``, which also
    implies ``--cache``). Rerunning the identical command replays every
    already-completed step from cache instead of recomputing it -- this is
    also how an interrupted or crashed run resumes, with no separate
    checkpoint mechanism. Not every step is safe to cache -- see
    :doc:`../user_guide/pipeline/caching`.

``--live``
    Replaces the static progress bar with a live-updating status table
    (pending/running/OK/ERR/cached, one row per step, rewritten in place).
    Implies visible progress even without ``-v``.

``--history`` / ``--history-file FILE``
    Appends a one-line JSON summary of the run -- pipeline name, status,
    timing, site counts, per-step summary -- to a run history log (default
    ``~/.pycsamt/pipeline_history.jsonl``, override with ``--history-file``,
    which also implies ``--history``). Query it back with:

    .. code-block:: console

       pycsamt pipe history
       pycsamt pipe history --last 5
       pycsamt pipe history --format json

    ``--file FILE`` points ``pipe history`` at a non-default log;
    ``--last N`` trims to the N most recent runs.

``--dashboard``
    Also writes ``dashboard.html`` -- a richer, branded report with KPI stat
    tiles and inline-SVG charts (step status, per-step duration, site-count
    flow) -- alongside the default ``report.html`` and ``summary.txt``. Full
    contents in :doc:`../user_guide/pipeline/outputs`'s "Dashboard Report"
    section.

QC-figure generation is not itself cache-aware -- a cache hit skips the step
*transform* only, so the same figures regenerate on every run regardless of
``--cache``. :doc:`../user_guide/pipeline/caching` isolates the transform-only
speed-up with ``--no-plots`` if you want a clean before/after timing
comparison; the practical value of ``--cache`` day to day is resuming a long
or interrupted run without recomputing what already succeeded.

Full mechanism for all four, including the Python API equivalents
(``cache=``, ``on_step=``, ``history=``, ``report_formats=``) in
:doc:`../user_guide/pipeline/observability` and
:doc:`../user_guide/pipeline/caching`.

Output directory
----------------

A normal run with output enabled writes a directory like this:

.. code-block:: text

   results/basic_qc/
   |-- processed/
   |   |-- station001.edi
   |   `-- station002.edi
   |-- plots/
   |   |-- 01_notch/
   |   |-- 02_drop_duplicates/
   |   `-- ...
   |-- pipeline.yaml
   |-- report.html
   |-- dashboard.html        # only with --dashboard
   `-- summary.txt

``pipeline.yaml``
    Canonical snapshot of the exact resolved pipeline used for the run.

``processed/``
    Final processed EDI files, unless ``--no-edi`` was used.

``plots/``
    Per-step QC figures, unless ``--no-plots`` was used. Exact figure names
    depend on the QC functions attached to each step.

``report.html`` and ``summary.txt``
    Human-readable run reports, unless ``--no-report`` was used.

``dashboard.html``
    Richer branded report, only written when ``--dashboard`` was passed.

The step :term:`cache <step cache>` and the run history log live outside
this tree, at ``~/.pycsamt/pipeline_cache`` and
``~/.pycsamt/pipeline_history.jsonl`` by default (``--cache-dir``/
``--history-file`` to relocate either).

Output formats
--------------

Most pipeline commands support ``--format text|json|csv``. Use text for
interactive runs, JSON for automation, and CSV for quick spreadsheet or
shell pipelines.

For ``pipe run --format json``, the terminal summary includes ``pipeline``,
``ok``, ``n_errors``, ``elapsed_sec``, input/output site counts, number of
plots, ``outdir``, and a ``steps`` list -- each step entry includes a
``cached`` boolean (``true`` when that step was replayed from cache).

For ``pipe run --format csv``, the header is:

.. code-block:: text

   idx,label,code,ok,elapsed_sec,n_sites_in,n_sites_out,n_plots,cached,error

The ``--format`` option controls terminal output only. It does not change
saved ``report.html`` or ``summary.txt`` files.

Exit status and troubleshooting
-------------------------------

``pipe run`` exits with status ``0`` only when ``PipelineResult.ok`` is true.
If one or more steps failed and the command continued under ``--on-error
warn`` or ``--on-error skip``, pyCSAMT prints the summary and exits nonzero.

Common failures:

No pipeline specified
    Provide one of ``--config``, ``--preset``, or ``--steps``.

Unknown preset
    Run ``pycsamt pipe presets`` and use an exact preset name.

Unknown step
    Run ``pycsamt pipe steps`` or ``pycsamt pipe steps --info CODE``. Step
    identifiers can be codes such as ``NR001`` or registry names such as
    ``notch_powerline``.

Unsupported config format
    Use ``.yaml``, ``.yml``, ``.json``, or ``.py``.

Step not found during slicing
    ``--from-step`` and ``--until-step`` match pipeline labels, registry
    codes, or registry names. Use ``pipe show`` to inspect the exact labels.

Near-field correction warns for every station
    Expected, not a failure, when a survey's EDI headers carry no resolvable
    ``source_offset``/``offset``/``dist``. ``csamt_qc``/``csumt_qc`` run
    near-field correction with ``source_offset=None`` precisely so this
    warns and passes each station through uncorrected instead of failing
    the run.

Recommended workflow
--------------------

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L18PLT
   pycsamt pipe presets
   pycsamt pipe presets --expand basic_qc
   pycsamt pipe init --preset basic_qc --name l18_qc -o l18_qc.yaml
   pycsamt pipe show l18_qc.yaml
   pycsamt pipe run --config l18_qc.yaml --dry-run
   pycsamt pipe run --config l18_qc.yaml --out results/l18_qc/ --cache --history

Once a config stabilises, add ``--cache`` and ``--history`` as shown above so
later iterations resume instead of recomputing and every run is logged --
see `Cache, live progress, history, and the dashboard`_.

Python equivalent
-----------------

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.emtools._core import ensure_sites
   >>> from pycsamt.pipeline import Pipeline
   >>> sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True, verbose=0)
   >>> pipe = Pipeline.from_yaml("l18_qc.yaml")
   >>> result = pipe.run(sites, outdir="results/l18_qc")
   >>> print(result.summary())
   PipelineResult  'l18_qc'
     Sites   : 28 in → 28 out
     Steps   : 3 (3 ok, 0 err)
     Time    : 23.97 s
     Plots   : 7
     Output  : results/l18_qc

Related pages
-------------

For deeper background, see the pipeline pages under
``docs/source/user_guide/pipeline/``: :doc:`../user_guide/pipeline/concepts`
for the object model, :doc:`../user_guide/pipeline/configuration_files` for
config schemas, :doc:`../user_guide/pipeline/presets` for built-in workflows
including the method-aware ``mt_qc``/``amt_qc``/``csamt_qc``/``csumt_qc``
group, :doc:`../user_guide/pipeline/steps` for the step registry,
:doc:`../user_guide/pipeline/caching` for the step cache,
:doc:`../user_guide/pipeline/observability` for live progress and the run
history log, :doc:`../user_guide/pipeline/extending` for plugins and the
opt-in AI step, and :doc:`../user_guide/pipeline/outputs` for the output
directory contract including the dashboard report. The worked
:doc:`../tutorials/run_pipeline_from_config` tutorial walks all of this
end to end on real data.

