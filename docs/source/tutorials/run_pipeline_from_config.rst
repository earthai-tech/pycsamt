Run a Pipeline From Config
==========================

This tutorial shows the reproducible pipeline workflow: define processing steps
in a configuration file, load that file, run the pipeline, and inspect outputs.

Why use a pipeline?
-------------------

Use a pipeline when a workflow should be repeatable across surveys or shared
with collaborators. A pipeline stores ordered processing steps and saves run
metadata with the outputs.

Minimal Python workflow
-----------------------

Load a survey and construct a pipeline from a preset:

.. code-block:: python

   from pycsamt.api import read_edis
   from pycsamt.pipeline import Pipeline

   survey = read_edis("data/edis")
   sites = survey.collection

   pipe = Pipeline.from_preset("publication_ready")
   result = pipe.run(sites, outdir="pipeline_results")

   print(result.summary())

Load from YAML
--------------

When a workflow is stored in YAML:

.. code-block:: python

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_yaml("workflow.yaml")
   result = pipe.run(sites, outdir="pipeline_results")

A typical YAML file records a name and an ordered list of steps. The exact
schema is documented in :doc:`../pipeline/configuration_files`.

Discover presets and steps
--------------------------

.. code-block:: python

   from pycsamt.pipeline import list_steps, preset_catalogue

   print(preset_catalogue())
   print(list_steps())

Use step discovery when writing a new configuration file or checking whether a
workflow uses the intended processing operations.

CLI equivalent
--------------

.. code-block:: bash

   pycsamt pipe presets
   pycsamt pipe steps
   pycsamt pipe init workflow.yaml
   pycsamt pipe run workflow.yaml --input data/edis --output pipeline_results

Expected outputs
----------------

A pipeline run can write:

- processed EDI files
- figures
- tables
- run configuration
- processing report

The returned ``PipelineResult`` also exposes step-level status and generated
paths:

.. code-block:: python

   print(result.ok)
   print(result.plots)
   print(result.processed_paths)

See Also
--------
:doc:`../pipeline/index`
    Pipeline system documentation.
:doc:`../api/pipeline`
    Pipeline API reference.
:doc:`../cli/pipe`
    Pipeline CLI reference.
