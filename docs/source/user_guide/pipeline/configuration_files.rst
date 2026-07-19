.. _pipeline-configuration-files:

Pipeline Configuration Files
============================

A :term:`pipeline configuration file` is the portable description of a
:term:`processing pipeline`.  It says which registered operations should run,
in which order, with which :term:`parameter override` values, and where the
default outputs should go.  The same file can be reviewed in version control,
run from the CLI, loaded from Python, and archived beside a
:term:`provenance manifest` before inversion.

pyCSAMT supports three equivalent configuration formats:

* YAML, loaded with :meth:`pycsamt.pipeline.Pipeline.from_yaml`;
* JSON, loaded with :meth:`pycsamt.pipeline.Pipeline.from_json`;
* Python, loaded with :meth:`pycsamt.pipeline.Pipeline.from_py`.

YAML is the recommended project format.  It is readable, compact, easy to
diff, and directly produced by ``pycsamt pipe init``.  JSON is useful when an
external program generates the recipe.  Python is useful for trusted local
configuration files that need constants or small bits of logic, but it is
imported and executed as Python code, so do not use it for untrusted files.

When To Use A Configuration File
--------------------------------

Use a configuration file when a run should be reproducible rather than merely
convenient.  A notebook cell such as
``Pipeline.from_preset("basic_qc").run(sites)`` is fine while exploring a
small :term:`site collection`, but it leaves too much context in memory: which
output directory was chosen, which extra step was added, and which frequency
limits were changed after the first attempt.

A configuration file makes those choices explicit.  It is especially helpful
when several surveys should share the same chain, when another geophysicist
must review the processing before inversion, or when the command line and
Python API need to execute the same workflow.  Treat the file as the recipe
for the data product; raw EDI files remain unchanged, and processed products
are written under a results directory.

Basic Schema
------------

The top-level object is a mapping:

.. list-table::
   :header-rows: 1
   :widths: 24 18 58

   * - Key
     - Required
     - Meaning
   * - ``name``
     - No
     - Human-readable pipeline name used in reports and printed summaries.
       It defaults to ``"unnamed"``.
   * - ``output_dir``
     - No
     - Default output directory used when ``Pipeline.run`` is called without
       an explicit ``outdir``.
   * - ``preset``
     - No
     - Built-in :term:`pipeline preset` used to seed the pipeline before the
       explicit ``steps`` list is appended.
   * - ``steps``
     - No
     - Ordered list of step entries.  Each entry identifies one registered
       operation and optional parameter overrides.

Each item in ``steps`` is also a mapping:

.. list-table::
   :header-rows: 1
   :widths: 24 18 58

   * - Key
     - Required
     - Meaning
   * - ``code``
     - Recommended
     - :term:`Pipeline step code` such as ``"NR001"`` or registry name such
       as ``"notch_powerline"``.  The loader can fall back to ``name`` as an
       identifier, but that makes labels ambiguous.
   * - ``name``
     - No
     - :term:`Step label` for this occurrence of the operation.  Labels appear
       in reports and can be used by CLI slicing options such as
       ``--from-step``.
   * - ``params``
     - No
     - Keyword arguments passed to the step.  These override defaults from the
       :term:`step registry`.

The loader can be understood as a small deterministic transformation.  Let
``S_p`` be the ordered step list supplied by a preset, and let ``S_c`` be the
ordered step list written explicitly in the configuration file.  The final
pipeline order is

.. math::

   S_{\mathrm{final}} =
   \begin{cases}
   S_p \Vert S_c, & \text{if a preset is given},\\
   S_c, & \text{otherwise},
   \end{cases}

where ``\Vert`` means append.  A preset is therefore a seed, not a template
that can be edited in place.  If you need to change a parameter inside a
preset step, write the full step list explicitly.

For one explicit step ``j``, the identifier, label, and parameters are resolved
as

.. math::

   i_j &= \mathrm{code}_j \;\text{if present, otherwise}\; \mathrm{name}_j,\\
   \ell_j &= \mathrm{name}_j \;\text{if present, otherwise}\;
             \mathrm{registry\_name}(i_j),\\
   \theta_j &= \theta_{\mathrm{default}}(i_j)
              \cup \theta_{\mathrm{params},j}.

The last line means that values in ``params`` replace registry defaults for
the same keys, while omitted keys keep their defaults.  This makes short
configs possible, but it also means reviewers should know which defaults were
in effect for the installed pyCSAMT version.

Minimal YAML Example
--------------------

This is a complete YAML pipeline:

.. code-block:: yaml
   :linenos:

   name: first_qc
   output_dir: results/first_qc

   steps:
     - name: notch
       code: NR001
       params:
         mains_hz: 50.0
         n_harm: 30
         tol_hz: 0.08

     - name: drop_duplicates
       code: FREQ002

     - name: select_band
       code: FREQ001
       params:
         band_hz: [0.001, 10000.0]

     - name: align_grid
       code: FREQ004

     - name: qc_snapshot
       code: QC001

Load it from Python and inspect the resolved pipeline before running it:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> pipe = Pipeline.from_yaml("config/first_qc.yaml")
   >>> pipe.name
   'first_qc'
   >>> [label for label, step in pipe]
   ['notch', 'drop_duplicates', 'select_band', 'align_grid', 'qc_snapshot']
   >>> "output_dir: results/first_qc" in pipe.to_yaml_string()
   True

Then run the same file against a survey:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.api import read_edis
   >>> survey = read_edis("data/3edis", strict=False)
   >>> result = pipe.run(survey.to_collection())
   >>> print(result.summary())  # doctest: +ELLIPSIS
   Pipeline run: first_qc
   Steps: 5
   Output directory: results/first_qc
   ...

If ``outdir`` is passed to ``Pipeline.run``, it takes precedence over the
file's ``output_dir``:

.. code-block:: pycon
   :linenos:

   >>> result = pipe.run(
   ...     survey.to_collection(),
   ...     outdir="results/first_qc_experiment",
   ... )
   >>> result.outdir
   'results/first_qc_experiment'

In words, the effective output directory is chosen in this order:

.. math::

   d_{\mathrm{effective}} =
   d_{\mathrm{run}} \rightarrow d_{\mathrm{config}}
   \rightarrow d_{\mathrm{default}},

where the first available value wins.  This is useful for controlled
experiments: the configuration still describes the processing chain, while the
runtime call can place trial outputs in a separate directory.

Generate A Starter Config
-------------------------

The easiest way to create a valid file is the CLI scaffold command:

.. code-block:: console
   :linenos:

   pycsamt pipe init \
       --preset basic_qc \
       --name first_qc \
       --outdir results/first_qc \
       --output config/first_qc.yaml

Print the scaffold without writing a file:

.. code-block:: console
   :linenos:

   pycsamt pipe init --preset basic_qc --name first_qc \
       --outdir results/first_qc --print

Captured output excerpt:

.. code-block:: yaml
   :linenos:

   # pyCSAMT Pipeline Configuration
   # Generated by: Pipeline.scaffold("first_qc.yaml")

   name: first_qc
   output_dir: results/first_qc

   steps:
     - {name: select_band, code: FREQ001,
        params: {band_hz: [0.001, 10000.0]}}
     - {name: drop_duplicates, code: FREQ002}
     - {name: align_grid, code: FREQ004}
     - {name: notch_powerline, code: NR001,
        params: {mains_hz: 50, n_harm: 30, tol_hz: 0.08}}
     - {name: qc_snapshot, code: QC001}

The full scaffold also includes commented inactive steps.  Keep the active
steps you want, remove those you do not want, rename labels to match the
survey, and adjust ``params`` after checking the step information.

Generate Python or JSON instead of YAML:

.. code-block:: console
   :linenos:

   pycsamt pipe init --format py --preset basic_qc -o config/first_qc.py
   pycsamt pipe init --format json --preset basic_qc -o config/first_qc.json

YAML, JSON, And Python Formats
------------------------------

All formats express the same logical object.  YAML is concise:

.. code-block:: yaml
   :linenos:

   name: amt_line_22
   output_dir: results/line_22
   steps:
     - {name: notch, code: NR001, params: {mains_hz: 50.0}}
     - {name: select_band, code: FREQ001, params: {band_hz: [10.0, 100000.0]}}
     - {name: qc, code: QC001}

JSON is better when another program writes the file:

.. code-block:: json
   :linenos:

   {
     "name": "amt_line_22",
     "output_dir": "results/line_22",
     "steps": [
       {"name": "notch", "code": "NR001", "params": {"mains_hz": 50.0}},
       {
         "name": "select_band",
         "code": "FREQ001",
         "params": {"band_hz": [10.0, 100000.0]}
       },
       {"name": "qc", "code": "QC001"}
     ]
   }

Python config files must define a module-level ``pipeline_config``
dictionary:

.. code-block:: pycon
   :linenos:

   >>> AMT_BAND_HZ = (10.0, 100000.0)
   >>> pipeline_config = dict(
   ...     name="amt_line_22",
   ...     output_dir="results/line_22",
   ...     steps=[
   ...         dict(name="notch", code="NR001", params=dict(mains_hz=50.0)),
   ...         dict(
   ...             name="select_band",
   ...             code="FREQ001",
   ...             params=dict(band_hz=AMT_BAND_HZ),
   ...         ),
   ...         dict(name="qc", code="QC001"),
   ...     ],
   ... )
   >>> pipeline_config["steps"][1]["params"]["band_hz"]
   (10.0, 100000.0)

Save those lines in ``config/line_22.py`` and load them with:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> pipe = Pipeline.from_py("config/line_22.py")
   >>> pipe.name
   'amt_line_22'

Step Codes And Labels
---------------------

Every operation is resolved through the :term:`step registry`.  A config step
can use the short :term:`pipeline step code`:

.. code-block:: yaml
   :linenos:

   - name: notch
     code: NR001

or the registry name:

.. code-block:: yaml
   :linenos:

   - name: notch
     code: notch_powerline

The code form is compact and stable in reports.  The ``name`` field is not the
registry name; it is the :term:`step label` for this occurrence of the step.
Use labels that describe the role of the operation in this workflow:

.. code-block:: yaml
   :linenos:

   - name: remove_powerline
     code: NR001
   - name: trim_to_amt_band
     code: FREQ001

Labels are useful when slicing a run from the CLI:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis --config config/line_22.yaml \
       --from-step trim_to_amt_band

Discover Valid Steps
--------------------

Use the CLI when working in a terminal:

.. code-block:: console
   :linenos:

   pycsamt pipe steps
   pycsamt pipe steps --category frequency
   pycsamt pipe steps --info NR001
   pycsamt pipe steps --codes-only

Captured output for the frequency codes:

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

Or inspect the same registry from Python:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> "FREQ001" in Pipeline.catalogue("frequency")
   True
   >>> info = Pipeline.step_info("NR001")
   >>> "notch" in info.lower()
   True

Step defaults are merged with ``params``.  For example, if ``NR001`` has
defaults for ``mains_hz``, ``n_harm``, and ``tol_hz``, this config overrides
only ``mains_hz`` and keeps the remaining defaults:

.. code-block:: yaml
   :linenos:

   - name: notch_60hz
     code: NR001
     params:
       mains_hz: 60.0

That compact form is helpful, but for publication or hand-off work it is often
better to record every important parameter explicitly.  A future reader should
not have to guess whether a value came from the file or from the installed
registry defaults.

Preset Plus Extra Steps
-----------------------

A config may combine a preset with additional explicit steps:

.. code-block:: yaml
   :linenos:

   name: publication_with_extra_qc
   output_dir: results/publication_with_extra_qc
   preset: publication_ready

   steps:
     - name: final_frequency_confidence
       code: QC001

pyCSAMT loads ``publication_ready`` first and then appends
``final_frequency_confidence``.  It does not replace, remove, or mutate a step
inside the preset.  This pattern is good for a known baseline plus extra
diagnostics.  It is not good for changing a preset parameter, because the
result would contain the original preset step and your additional step.

Full Explicit Config From A Preset
----------------------------------

If you want ``basic_qc`` with one changed parameter, prefer an explicit file:

.. code-block:: yaml
   :linenos:

   name: basic_qc_60hz
   output_dir: results/basic_qc_60hz

   steps:
     - name: select_band
       code: FREQ001
       params:
         band_hz: [0.001, 10000.0]
     - name: drop_duplicates
       code: FREQ002
     - name: align_grid
       code: FREQ004
     - name: notch_powerline
       code: NR001
       params:
         mains_hz: 60.0
         n_harm: 30
         tol_hz: 0.08
     - name: qc_snapshot
       code: QC001

This is longer than ``preset: basic_qc``, but it is unambiguous.  The reviewer
can see the exact sequence:

.. math::

   \mathrm{select\_band}
   \rightarrow \mathrm{drop\_duplicates}
   \rightarrow \mathrm{align\_grid}
   \rightarrow \mathrm{notch\_powerline}
   \rightarrow \mathrm{qc\_snapshot}.

The same notation is useful when comparing two processing branches.  If the
only difference is ``mains_hz = 50`` versus ``mains_hz = 60``, the config file
shows that the experiment changed the notch target rather than the frequency
grid, QC logic, or output path.

Run A Config From The CLI
-------------------------

Run a config against an explicit EDI directory:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
       --config config/first_qc.yaml \
       --out results/first_qc_run \
       --on-error warn \
       --dpi 200 \
       --plot-fmt png

Dry-run before a long processing job:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/3edis \
       --config config/first_qc.yaml \
       --dry-run

The CLI chooses the pipeline definition in this priority order:

1. ``--config``;
2. ``--preset``;
3. ``--steps``.

If ``--config`` is provided, ``--preset`` and ``--steps`` are ignored because
the file is the source of truth for the processing chain.  Runtime options
such as ``--out``, ``--dpi``, ``--plot-fmt``, and ``--on-error`` still control
how that chain is executed and written.

Export An Existing Pipeline
---------------------------

You can build or modify a pipeline in Python and export a
:term:`canonical pipeline snapshot`:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline
   >>> pipe = Pipeline.from_preset("basic_qc", pipeline_name="first_qc")
   >>> pipe.to_yaml("config/first_qc.yaml")
   >>> pipe.to_json("config/first_qc.json")
   >>> pipe.to_py("config/first_qc.py")

The YAML and JSON exports are useful for reproducibility because they are data
files.  The Python export is useful when you want an editable script-style
config with comments.  Archive at least one exported snapshot beside processed
EDI files, plots, and reports so a later user can reconstruct the chain that
produced them.

Validation And Failure Modes
----------------------------

pyCSAMT validates configuration files when they are loaded.  Common failures
include:

Top level is not a mapping
    YAML must load to a mapping and JSON must load to an object.  A top-level
    list is invalid because it cannot carry the pipeline name, output
    directory, or preset.

Python file has no ``pipeline_config``
    ``Pipeline.from_py`` imports the file and looks for a module-level
    variable named ``pipeline_config``.

Step entry has no usable identifier
    Every explicit step entry must identify a registry step.  Write ``code``
    in normal configs.  Falling back to ``name`` is supported, but it prevents
    ``name`` from being a clear label.

Unknown step code
    The code or registry name does not exist.  Run ``pycsamt pipe steps`` or
    ``Pipeline.step_info(...)`` to confirm the identifier.

Unknown preset
    The value under ``preset`` is not registered.  Run
    ``pycsamt pipe presets`` to list available presets.

Parameter name is wrong
    The config may load, but the step can fail at runtime if a parameter is
    not accepted by the underlying function.  Check
    ``pycsamt pipe steps --info CODE`` and run with ``--dry-run`` before
    processing the full survey.

Output directory is ambiguous
    Prefer survey-specific directories such as ``results/line_22/basic_qc``.
    Reusing a shared directory makes it harder to know which config produced
    which plots, processed EDI files, and reports.

Recommended Project Layout
--------------------------

Keep pipeline configs near the survey project, not mixed with raw data:

.. code-block:: text
   :linenos:

   survey_line_22/
   |-- config/
   |   |-- basic_qc.yaml
   |   |-- publication_ready.yaml
   |   `-- occam2d_prep.yaml
   |-- data/
   |   `-- edis/
   |-- results/
   |   |-- basic_qc/
   |   `-- publication_ready/
   `-- notes/

Keep the raw EDI directory unchanged.  Write processed data, plots, reports,
exported snapshots, and any copied configuration file under ``results/``.
That separation matters because a processing pipeline is not the source data;
it is a deterministic transformation applied to a survey state.

Best Practices
--------------

* Commit YAML configuration files with the project when possible.
* Use Python configs only for trusted local logic.
* Give every explicit step a meaningful ``name`` label.
* Prefer explicit step lists when changing preset parameters.
* Record important parameter values instead of relying silently on defaults.
* Keep output directories survey-specific.
* Run ``pycsamt pipe run ... --dry-run`` before long jobs.
* Store raw data and processed outputs in separate directories.
* Archive the config file used to prepare inversion inputs.

In Short
--------

A pyCSAMT pipeline config is an ordered, reproducible processing recipe:

.. code-block:: yaml
   :linenos:

   name: first_qc
   output_dir: results/first_qc
   steps:
     - {name: notch, code: NR001}
     - {name: drop_duplicates, code: FREQ002}
     - {name: select_band, code: FREQ001,
        params: {band_hz: [0.001, 10000.0]}}
     - {name: qc_snapshot, code: QC001}

Load it with ``Pipeline.from_yaml`` or run it with
``pycsamt pipe run --config``.  Use ``pycsamt pipe init`` when you want a
valid starter file, and export a :term:`canonical pipeline snapshot` when a
constructed pipeline should become part of the permanent processing record.
