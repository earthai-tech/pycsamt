.. _pipeline-configuration-files:

Pipeline Configuration Files
============================

Pipeline configuration files make processing workflows reproducible.  Instead
of rebuilding a :class:`pycsamt.pipeline.Pipeline` in a notebook every time,
you can store the pipeline name, default output directory, preset seed, step
order, and step parameters in a small file.

pyCSAMT supports three configuration formats:

* YAML, loaded with :meth:`pycsamt.pipeline.Pipeline.from_yaml`;
* JSON, loaded with :meth:`pycsamt.pipeline.Pipeline.from_json`;
* Python, loaded with :meth:`pycsamt.pipeline.Pipeline.from_py`.

All three formats use the same logical schema.  YAML is the recommended
default for most survey projects because it is readable, easy to review in
version control, and directly supported by the CLI scaffold command.

When To Use A Configuration File
--------------------------------

Use a pipeline configuration file when:

* a processing run must be repeated later;
* several surveys should use the same processing chain;
* a workflow needs review by another developer or geophysicist;
* you want the command line and Python API to run the same steps;
* you need a permanent record of parameters used before inversion;
* you are preparing examples, tutorials, tests, or reports.

For quick exploration, ``Pipeline.from_preset("basic_qc")`` is fine.  For
survey processing, inversion preparation, or publication output, write the
workflow to a configuration file.

Basic Schema
------------

The top-level configuration is a mapping with these keys:

.. list-table::
   :header-rows: 1
   :widths: 24 18 58

   * - Key
     - Required
     - Meaning
   * - ``name``
     - No
     - Human-readable pipeline name used in reports and printed summaries.
       Defaults to ``"unnamed"`` when omitted.
   * - ``output_dir``
     - No
     - Default output directory used when ``Pipeline.run`` is called without
       an explicit ``outdir``.
   * - ``preset``
     - No
     - Built-in preset name used to seed the pipeline before explicit
       ``steps`` are appended.
   * - ``steps``
     - No
     - Ordered list of step entries.  Each entry identifies a registered
       pipeline step and optional parameter overrides.

Each item in ``steps`` is a mapping:

.. list-table::
   :header-rows: 1
   :widths: 24 18 58

   * - Key
     - Required
     - Meaning
   * - ``code``
     - Recommended
     - Step registry code such as ``"NR001"`` or registry name such as
       ``"notch_powerline"``.  The loader can fall back to ``name`` as the
       identifier, but the code form is clearer and is recommended for
       configs.
   * - ``name``
     - No
     - User label for this occurrence of the step.  Labels appear in reports
       and can be used by CLI slicing options such as ``--from-step``.
   * - ``params``
     - No
     - Keyword arguments passed to the step.  These override registry
       defaults.

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

Load and run it from Python:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis
   from pycsamt.pipeline import Pipeline

   survey = read_edis("data/edis", strict=False)
   pipe = Pipeline.from_yaml("config/first_qc.yaml")

   print(pipe)
   result = pipe.run(survey.to_collection())
   print(result.summary())

Because the file defines ``output_dir``, the run writes to
``results/first_qc`` unless you override it:

.. code-block:: python
   :linenos:

   result = pipe.run(
       survey.to_collection(),
       outdir="results/first_qc_experiment",
   )

The explicit ``outdir`` passed to ``Pipeline.run`` wins over the file's
``output_dir``.

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

   pycsamt pipe init --preset full_processing --print

Generate Python or JSON instead of YAML:

.. code-block:: console
   :linenos:

   pycsamt pipe init --format py --preset basic_qc -o config/first_qc.py
   pycsamt pipe init --format json --preset basic_qc -o config/first_qc.json

The scaffold includes active steps from the chosen preset and comments for
other registered steps.  Treat it as a starting point: remove steps you do
not want, rename labels, and adjust ``params`` for the survey.

YAML, JSON, And Python Formats
------------------------------

YAML is the most convenient hand-edited format:

.. code-block:: yaml
   :linenos:

   name: amt_line_22
   output_dir: results/line_22
   steps:
     - {name: notch, code: NR001, params: {mains_hz: 50.0}}
     - {name: select_band, code: FREQ001, params: {band_hz: [10.0, 100000.0]}}
     - {name: qc, code: QC001}

JSON is useful for generated configs or external tooling:

.. code-block:: json
   :linenos:

   {
     "name": "amt_line_22",
     "output_dir": "results/line_22",
     "steps": [
       {
         "name": "notch",
         "code": "NR001",
         "params": {
           "mains_hz": 50.0
         }
       },
       {
         "name": "select_band",
         "code": "FREQ001",
         "params": {
           "band_hz": [10.0, 100000.0]
         }
       },
       {
         "name": "qc",
         "code": "QC001"
       }
     ]
   }

Python config files are useful when you want comments, constants, or a small
amount of local logic.  The file must define a module-level
``pipeline_config`` dictionary:

.. code-block:: python
   :linenos:

   AMT_BAND_HZ = (10.0, 100000.0)

   pipeline_config = dict(
       name="amt_line_22",
       output_dir="results/line_22",
       steps=[
           dict(name="notch", code="NR001", params=dict(mains_hz=50.0)),
           dict(name="select_band", code="FREQ001",
                params=dict(band_hz=AMT_BAND_HZ)),
           dict(name="qc", code="QC001"),
       ],
   )

Load it with:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_py("config/line_22.py")

Use Python configs carefully.  They are imported and executed as Python code,
so they are best kept inside trusted project repositories.

Step Codes And Labels
---------------------

Every step is resolved through the pipeline registry.  A config step can use
the short registry code:

.. code-block:: yaml
   :linenos:

   - name: notch
     code: NR001

or the registry name:

.. code-block:: yaml
   :linenos:

   - name: notch
     code: notch_powerline

The code form is more compact and stable in reports.  The ``name`` field is
not the registry name; it is the label for this occurrence of the step.  Use
short labels that describe the role of the step in this workflow:

.. code-block:: yaml
   :linenos:

   - name: remove_powerline
     code: NR001
   - name: trim_to_amt_band
     code: FREQ001

Labels are useful when slicing a run from the CLI:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --config config/line_22.yaml \
       --from-step trim_to_amt_band

Discover Valid Steps
--------------------

Use the CLI:

.. code-block:: console
   :linenos:

   pycsamt pipe steps
   pycsamt pipe steps --category frequency
   pycsamt pipe steps --info NR001
   pycsamt pipe steps --codes-only

Or use Python:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   print(Pipeline.catalogue())
   print(Pipeline.catalogue("frequency"))
   print(Pipeline.step_info("NR001"))

Step defaults are merged with your ``params``.  For example, if ``NR001`` has
defaults for ``mains_hz``, ``n_harm``, and ``tol_hz``, this config overrides
only ``mains_hz`` and keeps the other defaults:

.. code-block:: yaml
   :linenos:

   - name: notch_60hz
     code: NR001
     params:
       mains_hz: 60.0

Preset Plus Extra Steps
-----------------------

A config may include ``preset`` and explicit ``steps``:

.. code-block:: yaml
   :linenos:

   name: publication_with_extra_qc
   output_dir: results/publication_with_extra_qc
   preset: publication_ready

   steps:
     - name: final_frequency_confidence
       code: QC001

When ``preset`` is present, pyCSAMT loads the preset first, then appends the
explicit ``steps`` list.  It does not replace or edit steps inside the preset.

Use this pattern when you want a known baseline plus extra diagnostics.  Do
not use it when you need to change a preset step parameter; in that case,
write the full step list explicitly so the final order and parameters are
obvious.

Full Explicit Config From A Preset
----------------------------------

If you want ``basic_qc`` with one changed parameter, prefer an explicit file:

.. code-block:: yaml
   :linenos:

   name: basic_qc_60hz
   output_dir: results/basic_qc_60hz

   steps:
     - name: notch
       code: NR001
       params:
         mains_hz: 60.0
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

This is longer than ``preset: basic_qc``, but it is unambiguous and easy to
review before processing field data.

Run A Config From The CLI
-------------------------

Run a config against an explicit EDI directory:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/first_qc.yaml \
       --out results/first_qc_run \
       --on-error warn \
       --dpi 200 \
       --plot-fmt png

Dry-run before a long processing job:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --config config/first_qc.yaml \
       --dry-run

The pipeline definition priority in the CLI is:

1. ``--config``;
2. ``--preset``;
3. ``--steps``.

If ``--config`` is provided, ``--preset`` and ``--steps`` are ignored because
the file is the source of truth for the pipeline.

Export An Existing Pipeline
---------------------------

You can build or modify a pipeline in Python and export it:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_preset("basic_qc", pipeline_name="first_qc")
   pipe.to_yaml("config/first_qc.yaml")
   pipe.to_json("config/first_qc.json")
   pipe.to_py("config/first_qc.py")

The YAML and JSON exports are useful for reproducibility.  The Python export
is useful when you want an editable script-style config with comments.

Validation And Failure Modes
----------------------------

pyCSAMT validates configuration files when they are loaded.  Common failures
include:

Top level is not a mapping
    YAML must load to a mapping and JSON must load to an object.  A top-level
    list is invalid.

Python file has no ``pipeline_config``
    ``Pipeline.from_py`` imports the file and looks for a module-level
    variable named ``pipeline_config``.

Step entry has no ``code``
    Every explicit step entry must identify a registry step.  If ``code`` is
    missing, the loader falls back to ``name`` as the step identifier, but
    this makes labels ambiguous.  Prefer always writing ``code``.

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
and exported pipeline YAML to ``results/``.

Best Practices
--------------

* Commit pipeline config files with the project when possible.
* Use YAML for shared survey workflows.
* Use Python configs only for trusted local logic.
* Give every step a meaningful ``name`` label.
* Prefer explicit step lists when changing preset parameters.
* Keep output directories survey-specific.
* Run ``pycsamt pipe run ... --dry-run`` before long jobs.
* Store raw data and processed outputs in separate directories.
* Record the config file used to prepare inversion inputs.

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
valid starter file instead of writing the schema by hand.
