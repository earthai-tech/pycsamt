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
    counts, plot paths, processed EDI paths, output directory, and errors.

During ``pipe run``, pyCSAMT loads EDI data into ``Sites``, applies each
configured step in order, writes optional outputs, and prints a terminal
summary.

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

The step catalogue is grouped by categories such as ``frequency``,
``noise_removal``, ``static_shift``, ``tensor``, ``dimensionality``,
``skew``, ``source_effects``, and ``qc``.

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

Presets are named, opinionated workflows. Common names include
``basic_qc``, ``noise_reduction``, ``full_processing``,
``tensor_analysis``, ``dimensionality_filter``, and
``publication_ready``. ``stratagem_mt`` is also available for Stratagem AMT
data that is already loaded as ``Sites``. Use ``--expand`` to see the exact
step sequence.

Machine-readable output:

``pycsamt pipe presets --format json``
    Emits a list of preset objects with ``name``, ``description``,
    ``n_steps``, and ``codes``.

``pycsamt pipe presets --expand NAME --format json``
    Emits one preset object with the ordered ``steps`` list. Each step entry
    includes the label, code, registry name, category, and parameters.

``pycsamt pipe presets --format csv``
    Emits ``name,description,n_steps,codes``.

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
    Skip HTML/text run reports.

``--dpi INT`` and ``--plot-fmt png|pdf|svg``
    Control saved QC figure resolution and format.

If ``--out`` is omitted, output directory resolution follows this order:

1. the config file ``output_dir`` value, when available;
2. the global pipeline default, ``pipe_results``.

``--dry-run`` prints the resolved pipeline and site count but writes nothing.
Pass ``--out`` explicitly when the dry-run display should show the intended
target directory.

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

Output formats
--------------

Most pipeline commands support ``--format text|json|csv``. Use text for
interactive runs, JSON for automation, and CSV for quick spreadsheet or
shell pipelines.

For ``pipe run --format json``, the terminal summary includes ``pipeline``,
``ok``, ``n_errors``, ``elapsed_sec``, input/output site counts, number of
plots, ``outdir``, and a ``steps`` list.

For ``pipe run --format csv``, the header is:

.. code-block:: text

   idx,label,code,ok,elapsed_sec,n_sites_in,n_sites_out,n_plots,error

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

Recommended workflow
--------------------

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L18PLT
   pycsamt pipe presets
   pycsamt pipe presets --expand basic_qc
   pycsamt pipe init --preset basic_qc --name l18_qc -o l18_qc.yaml
   pycsamt pipe show l18_qc.yaml
   pycsamt pipe run --config l18_qc.yaml --dry-run
   pycsamt pipe run --config l18_qc.yaml --out results/l18_qc/

Python equivalent
-----------------

.. code-block:: python

   from pycsamt.emtools._core import ensure_sites
   from pycsamt.pipeline import Pipeline

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True, verbose=0)
   pipe = Pipeline.from_yaml("l18_qc.yaml")
   result = pipe.run(sites, outdir="results/l18_qc")

Related pages
-------------

For deeper background, see the pipeline pages under ``docs/source/pipeline``:
``concepts.rst`` for the object model, ``configuration_files.rst`` for config
schemas, ``presets.rst`` for built-in workflows, ``steps.rst`` for the step
registry, and ``outputs.rst`` for the output directory contract.

