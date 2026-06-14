Orchestration, Pipeline, And Output Agents
==========================================

These agents organize complete workflows and convert processing results into
interpretable, exportable, or reproducible products.

.. _agent-workflow-orchestrator:

WorkflowOrchestratorAgent
-------------------------

``WorkflowOrchestratorAgent`` classifies a natural-language request, selects a
workflow type, builds the corresponding ``AgentCoordinator`` chain, and runs
or previews the workflow.

Supported routed workflows include QC, phase analysis, pre-inversion,
AI inversion, Inv2D, Inv3D, ensemble inversion, joint inversion, ModEM
preparation, and full workflows.

.. code-block:: python
   :linenos:

   from pycsamt.agents import WorkflowOrchestratorAgent

   result = WorkflowOrchestratorAgent().execute({
       "request": "Run phase tensor analysis and prepare a report",
       "data_path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_phase",
       "dry_run": True,
   })

   print(result["workflow_type"])
   print(result["steps"])

.. _agent-pipeline:

PipelineAgent
-------------

``PipelineAgent`` bridges :mod:`pycsamt.pipeline` with the agent framework.  It
can recommend a pipeline preset from a natural-language request or run a
specific preset or list of step codes directly.

.. code-block:: python
   :linenos:

   from pycsamt.agents import PipelineAgent

   result = PipelineAgent(preset="full_processing").execute({
       "path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_pipeline",
   })

   print(result.get("preset_used"))
   print(result.get("steps_run"))

.. _agent-batch-survey:

BatchSurveyAgent
----------------

``BatchSurveyAgent`` runs agent workflows over multiple surveys, profiles, or
directories.  Use it when the same processing recipe should be applied to many
survey lines with consistent outputs.

.. code-block:: python
   :linenos:

   result = BatchSurveyAgent().execute({
       "paths": ["/data/line01", "/data/line02", "/data/line03"],
       "workflow": "qc",
       "output_dir": "/out/batch_qc",
   })

.. _agent-interpretation:

InterpretationAgent
-------------------

``InterpretationAgent`` converts resistivity and inversion products into
geological interpretation.  It can map resistivity ranges to lithology and
summarize likely hydrogeological or structural meaning.

.. code-block:: python
   :linenos:

   result = InterpretationAgent().execute({
       "model": inversion_model,
       "context": "basement aquifer survey",
   })

   print(result.llm_interpretation)

.. _agent-resistivity-map:

ResistivityMapAgent
-------------------

``ResistivityMapAgent`` creates horizontal depth-slice maps from inversion or
resistivity products.  Use it when a model needs map-view products for
interpretation or reporting.

.. code-block:: python
   :linenos:

   maps = ResistivityMapAgent().execute({
       "models": inversion_models,
       "depths_m": [50, 100, 250],
       "output_dir": "/out/resistivity_maps",
   })

.. _agent-sensitivity:

SensitivityAgent
----------------

``SensitivityAgent`` estimates depth of investigation, vertical resolution,
and sensitivity products.  Use it to qualify which depths or model regions are
well constrained.

.. code-block:: python
   :linenos:

   sensitivity = SensitivityAgent().execute({
       "sites": processed_sites,
       "model": inversion_model,
       "output_dir": "/out/sensitivity",
   })

.. _agent-edi-export:

EDIExportAgent
--------------

``EDIExportAgent`` writes processed survey data back to EDI files.  Use it
after correction, rotation, filtering, or denoising when downstream tools need
standard EDI products.

.. code-block:: python
   :linenos:

   exported = EDIExportAgent().execute({
       "sites": processed_sites,
       "output_dir": "/out/processed_edis",
   })

.. _agent-report:

ReportAgent
-----------

``ReportAgent`` assembles workflow outputs into human-readable reports.  It is
typically the final step in QC, phase-analysis, inversion, or AI workflows.

.. code-block:: python
   :linenos:

   report = ReportAgent().execute({
       "workflow_results": results,
       "output_dir": "/out/report",
       "formats": ["markdown", "html"],
   })

.. _agent-code-generation:

CodeGenerationAgent
-------------------

``CodeGenerationAgent`` generates a standalone Python script from workflow
configuration and results.  Use it when an interactive or LLM-assisted
workflow should be made reproducible.

.. code-block:: python
   :linenos:

   script = CodeGenerationAgent().execute({
       "workflow_config": config,
       "results": results,
       "output_dir": "/out/reproducible_script",
   })

Output Workflow Pattern
-----------------------

.. code-block:: text

   WorkflowOrchestratorAgent or AgentCoordinator
   -> PipelineAgent or processing agents
   -> InterpretationAgent
   -> ReportAgent
   -> CodeGenerationAgent

Use ``EDIExportAgent`` when corrected survey data must leave pyCSAMT as EDI
files.  Use ``BatchSurveyAgent`` when the same workflow should be applied to
several profiles.
