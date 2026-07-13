Agent Catalogue
===============

The pyCSAMT agent catalogue is the navigation map for the AI-assisted
workflow layer.  It lists the public classes exported by :mod:`pycsamt.agents`,
groups them by the role they play in a survey workflow, and points to the
detail page where each agent is documented.

Use this page when you need to answer one of these questions:

* Which agent should I start with?
* Which agents can operate before inversion?
* Which agents write files for Occam2D, ModEM, or MARE2DEM?
* Which agents belong to the AI/model-zoo layer?
* Which agents produce reports, maps, scripts, or EDI exports?

The catalogue intentionally separates **agent classes** from **support
interfaces**.  Modules such as ``web.py``, ``__main__.py``, and ``_pricing.py``
support agent usage, but they are not listed as workflow agents because they
do not define an executable ``BaseAgent`` subclass.


How to read the catalogue
-------------------------

Most users should think in terms of a workflow lifecycle:

.. code-block:: text

   request or config
       -> load survey
       -> QC and correction
       -> tensor / tipper / strike diagnostics
       -> denoise / decimate / rotate if needed
       -> forward model or inversion preparation
       -> evaluate, interpret, export, report

The agent groups follow that lifecycle.  Each group page contains deeper
examples and agent-specific input/output notes.

.. toctree::
   :maxdepth: 1
   :hidden:

   foundation_agents
   processing_agents
   inversion_agents
   ai_model_zoo_agents
   orchestration_output_agents


Catalogue groups
----------------

.. list-table::
   :header-rows: 1
   :widths: 24 46 30

   * - Group
     - Use it for
     - Detail page
   * - Foundation and survey intake
     - The execution contract, LLM-aware base class, result object, natural
       language request parsing, data loading, and explicit workflow chaining.
     - :doc:`foundation_agents`
   * - Processing and diagnostics
     - Data QC, static-shift correction, tensor diagnostics, strike analysis,
       tipper products, tensor rotation, frequency decimation, and denoising.
     - :doc:`processing_agents`
   * - Forward and inversion workflows
     - Synthetic forward modelling, inversion preparation, Occam2D, ModEM,
       MARE2DEM, backend execution, result evaluation, and model comparison.
     - :doc:`inversion_agents`
   * - AI and model-zoo agents
     - Neural 1-D, 2-D, and 3-D inversion, ensemble uncertainty, joint
       inversion, anomaly detection, and checkpoint discovery.
     - :doc:`ai_model_zoo_agents`
   * - Orchestration, pipeline, and outputs
     - Natural-language workflow routing, pyCSAMT pipeline execution, batch
       processing, interpretation, maps, sensitivity, EDI export, reports, and
       reproducible scripts.
     - :doc:`orchestration_output_agents`


Choosing the right entry point
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 34 38

   * - Need
     - Start with
     - Continue with
   * - A user writes a natural-language request.
     - :ref:`ContextInputAgent <agent-context-input>`
     - :ref:`WorkflowOrchestratorAgent <agent-workflow-orchestrator>` or
       :ref:`AgentCoordinator <agent-coordinator-class>`
   * - You already know the data path.
     - :ref:`MTLoaderAgent <agent-mt-loader>`
     - :ref:`DataQCAgent <agent-data-qc>`,
       :ref:`StaticShiftAgent <agent-static-shift>`, or
       :ref:`PhaseAnalysisAgent <agent-phase-analysis>`
   * - You want a reproducible fixed chain.
     - :ref:`AgentCoordinator <agent-coordinator-class>`
     - Any ordered set of processing, inversion, and output agents
   * - You want the request to choose the workflow.
     - :ref:`WorkflowOrchestratorAgent <agent-workflow-orchestrator>`
     - The orchestrator builds and runs the matching chain
   * - You want pyCSAMT pipeline presets.
     - :ref:`PipelineAgent <agent-pipeline>`
     - :doc:`../pipeline/index`
   * - You need inversion-ready files.
     - :ref:`InversionPrepAgent <agent-inversion-prep>`
     - :ref:`Occam2DAgent <agent-occam2d>`,
       :ref:`ModEmAgent <agent-modem>`, or
       :ref:`Mare2DEMAgent <agent-mare2dem>`
   * - You need AI inversion.
     - :ref:`AIInversionAgent <agent-ai-inversion>`
     - :ref:`Inv2DAgent <agent-inv2d>`,
       :ref:`Inv3DAgent <agent-inv3d>`, or
       :ref:`ModelZooAgent <agent-model-zoo>`
   * - You need deliverables.
     - :ref:`ReportAgent <agent-report>`
     - :ref:`CodeGenerationAgent <agent-code-generation>` and
       :ref:`EDIExportAgent <agent-edi-export>`


Foundation and survey intake
----------------------------

These components are the first layer of the agent system.  They define how
agents are built, how results are returned, and how a survey enters the
workflow.

.. list-table::
   :header-rows: 1
   :widths: 26 40 20 14

   * - Component
     - Main responsibility
     - Typical input
     - Detail
   * - :ref:`BaseAgent <agent-base-agent>`
     - Shared execution base for LLM access, provider/model resolution, cost
       tracking, JSON extraction, plotting helpers, and validation helpers.
     - Subclass code
     - :doc:`foundation_agents`
   * - :ref:`AgentResult <agent-result>`
     - Standard return object with ``status``, ``summary``, ``data``,
       ``warnings``, optional LLM interpretation, elapsed time, cost, and
       failure hints.
     - Agent output
     - :doc:`foundation_agents`
   * - :ref:`ContextInputAgent <agent-context-input>`
     - Convert natural-language requests into structured workflow
       configuration.  Falls back to deterministic parsing when no LLM is
       configured.
     - ``request`` text
     - :doc:`foundation_agents`
   * - :ref:`MTLoaderAgent <agent-mt-loader>`
     - Load EDI, AVG, J, path lists, existing ``Sites``, or compatible EDI
       collections into a validated pyCSAMT survey object.
     - ``path`` or ``sites``
     - :doc:`foundation_agents`
   * - :ref:`AgentCoordinator <agent-coordinator-class>`
     - Chain explicit agent steps with input mapping, dry-run previews,
       checkpoints, and cost aggregation.
     - Step graph
     - :doc:`foundation_agents`


Processing and diagnostics
--------------------------

Processing agents operate after loading and before inversion, interpretation,
or reporting.  They are useful in notebooks, desktop/web tools, batch
pipelines, and orchestrated workflows.

.. list-table::
   :header-rows: 1
   :widths: 25 43 18 14

   * - Agent
     - Main responsibility
     - Typical output
     - Detail
   * - :ref:`DataQCAgent <agent-data-qc>`
     - Assess station coverage, dead bands, outliers, frequency gaps, and
       survey quality indicators.
     - QC table, figures
     - :doc:`processing_agents`
   * - :ref:`StaticShiftAgent <agent-static-shift>`
     - Detect and correct near-surface static-shift effects using supported
       correction strategies.
     - Corrected sites
     - :doc:`processing_agents`
   * - :ref:`PhaseAnalysisAgent <agent-phase-analysis>`
     - Produce phase-tensor, skew, strike, dimensionality, Mohr, and Argand
       diagnostics.
     - Diagnostic figures
     - :doc:`processing_agents`
   * - :ref:`TensorRotationAgent <agent-tensor-rotation>`
     - Rotate impedance tensors to a target angle or strike reference while
       preserving survey metadata.
     - Rotated sites
     - :doc:`processing_agents`
   * - :ref:`TipperAnalysisAgent <agent-tipper-analysis>`
     - Analyze tipper amplitude, phase, real/imaginary induction arrows, and
       spatial induction-vector products.
     - Tipper summaries
     - :doc:`processing_agents`
   * - :ref:`FrequencyDecimationAgent <agent-frequency-decimation>`
     - Select stable, inversion-ready periods from dense, irregular, or noisy
       frequency sampling.
     - Selected periods
     - :doc:`processing_agents`
   * - :ref:`DenoisingAgent <agent-denoising>`
     - Apply robust and AI-assisted denoising before diagnostics, inversion,
       or AI training.
     - Denoised sites
     - :doc:`processing_agents`


Forward and inversion workflows
-------------------------------

These agents connect processed survey data to modelling and inversion
workflows.  Use the general agent when you want backend-agnostic behavior, and
the specialized agents when you already know the inversion code.

.. list-table::
   :header-rows: 1
   :widths: 25 43 18 14

   * - Agent
     - Main responsibility
     - Typical output
     - Detail
   * - :ref:`ForwardModelAgent <agent-forward-model>`
     - Run 1-D, 2-D, or 3-D forward modelling from resistivity models for
       synthetic checks and sensitivity experiments.
     - Forward response
     - :doc:`inversion_agents`
   * - :ref:`InversionPrepAgent <agent-inversion-prep>`
     - Prepare inversion-ready files through a backend-agnostic interface.
     - Input directory
     - :doc:`inversion_agents`
   * - :ref:`Occam2DAgent <agent-occam2d>`
     - Write Occam2D data, mesh, model, startup, and run-control files.
     - Occam2D project
     - :doc:`inversion_agents`
   * - :ref:`ModEmAgent <agent-modem>`
     - Prepare ModEM 3-D impedance data files and related model inputs.
     - ModEM project
     - :doc:`inversion_agents`
   * - :ref:`Mare2DEMAgent <agent-mare2dem>`
     - Prepare, run, or inspect MARE2DEM 2.5-D EM inversion projects,
       including data, resistivity, settings, and optional MPI execution.
     - MARE2DEM project
     - :doc:`inversion_agents`
   * - :ref:`InversionBackendAgent <agent-inversion-backend>`
     - Drive pyCSAMT inversion backends from an agent workflow rather than
       only writing external-code input files.
     - Backend result
     - :doc:`inversion_agents`
   * - :ref:`InversionEvaluationAgent <agent-inversion-evaluation>`
     - Load inversion outputs and compute RMS, residuals, misfit sections, and
       quality summaries.
     - Evaluation report
     - :doc:`inversion_agents`
   * - :ref:`InversionComparisonAgent <agent-inversion-comparison>`
     - Compare inversion sections, parameter sweeps, before/after corrections,
       or outputs from different backends.
     - Comparison figures
     - :doc:`inversion_agents`


AI and model-zoo agents
-----------------------

AI agents are useful when a trained model or model-zoo checkpoint is available,
when rapid approximate inversion is acceptable, or when anomaly/uncertainty
screening should complement deterministic processing.

.. list-table::
   :header-rows: 1
   :widths: 25 43 18 14

   * - Agent
     - Main responsibility
     - Typical output
     - Detail
   * - :ref:`AIInversionAgent <agent-ai-inversion>`
     - Run end-to-end 1-D neural inversion from observed survey data.
     - 1-D models
     - :doc:`ai_model_zoo_agents`
   * - :ref:`Inv2DAgent <agent-inv2d>`
     - Run 2-D profile inversion with U-Net style models and lateral
       continuity.
     - 2-D section
     - :doc:`ai_model_zoo_agents`
   * - :ref:`Inv3DAgent <agent-inv3d>`
     - Run 3-D spatial inversion using graph-based neural models and
       inter-station message passing.
     - 3-D volume
     - :doc:`ai_model_zoo_agents`
   * - :ref:`EnsembleAgent <agent-ensemble>`
     - Estimate uncertainty by combining multiple model predictions or
       ensemble members.
     - Uncertainty bands
     - :doc:`ai_model_zoo_agents`
   * - :ref:`JointInversionAgent <agent-joint-inversion>`
     - Run multi-modal inversion across MT, TEM, CSAMT, gravity, or paired
       datasets when available.
     - Joint model
     - :doc:`ai_model_zoo_agents`
   * - :ref:`AnomalyDetectionAgent <agent-anomaly-detection>`
     - Flag anomalous station-frequency samples, unusual profiles, or survey
       regions needing manual review.
     - Anomaly table
     - :doc:`ai_model_zoo_agents`
   * - :ref:`ModelZooAgent <agent-model-zoo>`
     - List, download, inspect, and use pre-trained EM inversion checkpoints.
     - Checkpoint path
     - :doc:`ai_model_zoo_agents`


Orchestration, pipeline, and outputs
------------------------------------

These agents sit above individual processing steps.  They route workflows,
execute batches, bridge the pyCSAMT pipeline system, and create deliverables.

.. list-table::
   :header-rows: 1
   :widths: 25 43 18 14

   * - Agent
     - Main responsibility
     - Typical output
     - Detail
   * - :ref:`WorkflowOrchestratorAgent <agent-workflow-orchestrator>`
     - Classify a natural-language request, choose a workflow type, assemble
       a chain, and run or preview it.
     - Workflow result
     - :doc:`orchestration_output_agents`
   * - :ref:`PipelineAgent <agent-pipeline>`
     - Recommend or run pyCSAMT pipeline presets and step lists from an agent
       workflow.
     - Pipeline result
     - :doc:`orchestration_output_agents`
   * - :ref:`BatchSurveyAgent <agent-batch-survey>`
     - Apply the same workflow to multiple lines, folders, or profiles with a
       consistent output structure.
     - Batch manifest
     - :doc:`orchestration_output_agents`
   * - :ref:`InterpretationAgent <agent-interpretation>`
     - Convert resistivity, inversion, and diagnostic products into
       geological or hydrogeological interpretation.
     - Interpretation text
     - :doc:`orchestration_output_agents`
   * - :ref:`ResistivityMapAgent <agent-resistivity-map>`
     - Generate horizontal depth-slice maps from inversion or resistivity
       products.
     - Map figures
     - :doc:`orchestration_output_agents`
   * - :ref:`SensitivityAgent <agent-sensitivity>`
     - Estimate depth of investigation, vertical resolution, and constrained
       model regions.
     - Sensitivity section
     - :doc:`orchestration_output_agents`
   * - :ref:`EDIExportAgent <agent-edi-export>`
     - Export corrected, rotated, filtered, or recomputed survey objects back
       to standard EDI files.
     - EDI files
     - :doc:`orchestration_output_agents`
   * - :ref:`ReportAgent <agent-report>`
     - Assemble figures, tables, warnings, and interpretation into Markdown,
       HTML, or PDF reports.
     - Report files
     - :doc:`orchestration_output_agents`
   * - :ref:`CodeGenerationAgent <agent-code-generation>`
     - Generate a standalone Python script from workflow configuration and
       outputs so interactive work can be reproduced.
     - Python script
     - :doc:`orchestration_output_agents`


Typical chains
--------------

Quality-control chain
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   MTLoaderAgent
   -> DataQCAgent
   -> PhaseAnalysisAgent
   -> ReportAgent

Correction and inversion-preparation chain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   MTLoaderAgent
   -> DataQCAgent
   -> StaticShiftAgent
   -> TensorRotationAgent
   -> FrequencyDecimationAgent
   -> Occam2DAgent or ModEmAgent or Mare2DEMAgent

AI inversion chain
~~~~~~~~~~~~~~~~~~

.. code-block:: text

   MTLoaderAgent
   -> DataQCAgent
   -> DenoisingAgent
   -> AIInversionAgent or Inv2DAgent or Inv3DAgent
   -> EnsembleAgent
   -> InterpretationAgent
   -> ReportAgent

Batch workflow chain
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   ContextInputAgent
   -> WorkflowOrchestratorAgent
   -> BatchSurveyAgent
   -> ReportAgent
   -> CodeGenerationAgent


Minimal examples
----------------

Load and QC a survey:

.. code-block:: python
   :linenos:

   from pycsamt.agents import MTLoaderAgent, DataQCAgent

   loaded = MTLoaderAgent().execute({"path": "data/edis"})
   qc = DataQCAgent().execute({"sites": loaded["sites"]})

   print(qc.status)
   print(qc.summary)

Prepare an Occam2D project from processed sites:

.. code-block:: python
   :linenos:

   from pycsamt.agents import Occam2DAgent

   occam = Occam2DAgent().execute({
       "sites": processed_sites,
       "output_dir": "outputs/occam2d",
       "run_external": False,
   })

   print(occam["output_dir"])

Run a dry-run orchestrated workflow:

.. code-block:: python
   :linenos:

   from pycsamt.agents import WorkflowOrchestratorAgent

   plan = WorkflowOrchestratorAgent().execute({
       "request": "Load data/edis, run QC, correct static shift, and report",
       "dry_run": True,
   })

   print(plan["workflow_type"])
   print(plan["steps"])


Support interfaces
------------------

The following modules support the agent system but are not workflow agents:

``pycsamt.agents._pricing``
    Cost estimation helpers used by :data:`pycsamt.agents.AGENT_CONFIG` and
    ``BaseAgent``.

``pycsamt.agents.web``
    Optional Gradio interface for interactive agent usage.

``pycsamt.agents.__main__``
    Command-line entry point for lightweight agent demonstrations and utility
    commands.

``pycsamt.api.agents``
    The global provider, key, model, pricing, and budget configuration layer.
    See :doc:`llm_configuration`.


Related API
-----------

* :mod:`pycsamt.agents`
* :doc:`llm_configuration`
* :doc:`foundation_agents`
* :doc:`processing_agents`
* :doc:`inversion_agents`
* :doc:`ai_model_zoo_agents`
* :doc:`orchestration_output_agents`
* :doc:`../api/agents`
