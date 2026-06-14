Agent Catalogue
===============

The pyCSAMT agent catalogue is organized by workflow role.  This page is the
navigation map: use it to choose an agent family, then follow the agent links
to the detailed group pages.

The inventory below reflects the public classes exported by
:mod:`pycsamt.agents` and the implemented modules in ``pycsamt/agents/``.
Implementation modules such as ``web.py`` and ``__main__.py`` are documented
as interfaces, not agents, because they do not define agent classes.

Catalogue Groups
----------------

.. list-table::
   :header-rows: 1
   :widths: 28 46 26

   * - Group
     - Use it for
     - Detail page
   * - Foundation and survey intake
     - Configuration parsing, data loading, workflow infrastructure, and the
       standard result contract.
     - :doc:`foundation_agents`
   * - Processing and diagnostics
     - QC, static shift, tensor analysis, tipper analysis, rotation,
       frequency selection, and denoising.
     - :doc:`processing_agents`
   * - Forward and inversion workflows
     - Forward modelling, inversion-file preparation, backend execution,
       evaluation, and model comparison.
     - :doc:`inversion_agents`
   * - AI and model-zoo agents
     - Neural inversion, uncertainty, joint inversion, anomaly detection, and
       pre-trained model access.
     - :doc:`ai_model_zoo_agents`
   * - Orchestration, pipeline, and outputs
     - Natural-language routing, pyCSAMT pipeline execution, batch surveys,
       interpretation, maps, exports, reports, and reproducible scripts.
     - :doc:`orchestration_output_agents`

.. toctree::
   :maxdepth: 1
   :hidden:

   foundation_agents
   processing_agents
   inversion_agents
   ai_model_zoo_agents
   orchestration_output_agents

Foundation And Survey Intake
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 52 20

   * - Component
     - Main responsibility
     - Detail
   * - :ref:`BaseAgent <agent-base-agent>`
     - Shared execution base for LLM access, cost tracking, plotting helpers,
       JSON extraction, and common validation utilities.
     - :doc:`foundation_agents`
   * - :ref:`AgentResult <agent-result>`
     - Standard result object returned by every agent, with dict-like access
       to structured ``data``.
     - :doc:`foundation_agents`
   * - :ref:`ContextInputAgent <agent-context-input>`
     - Convert natural-language requests into structured workflow
       configuration with LLM or regex fallback.
     - :doc:`foundation_agents`
   * - :ref:`MTLoaderAgent <agent-mt-loader>`
     - Load EDI, AVG, J, existing ``Sites``, or EDI collections, then produce
       station-level quality summaries.
     - :doc:`foundation_agents`
   * - :ref:`AgentCoordinator <agent-coordinator-class>`
     - Chain agents into named workflows with step inputs, dry-run previews,
       checkpoints, and cost aggregation.
     - :doc:`foundation_agents`

Processing And Diagnostics
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 52 20

   * - Agent
     - Main responsibility
     - Detail
   * - :ref:`DataQCAgent <agent-data-qc>`
     - Assess frequency coverage, SNR-like indicators, dead bands, outliers,
       and station quality.
     - :doc:`processing_agents`
   * - :ref:`StaticShiftAgent <agent-static-shift>`
     - Detect and correct static shift using supported deterministic
       correction strategies.
     - :doc:`processing_agents`
   * - :ref:`PhaseAnalysisAgent <agent-phase-analysis>`
     - Generate phase-tensor, strike, dimensionality, Mohr, skew, and Argand
       diagnostics.
     - :doc:`processing_agents`
   * - :ref:`TensorRotationAgent <agent-tensor-rotation>`
     - Rotate impedance tensors into a target coordinate frame or strike
       reference.
     - :doc:`processing_agents`
   * - :ref:`TipperAnalysisAgent <agent-tipper-analysis>`
     - Analyze tipper amplitude, phase, and induction-arrow products.
     - :doc:`processing_agents`
   * - :ref:`FrequencyDecimationAgent <agent-frequency-decimation>`
     - Select stable, inversion-ready periods from dense or uneven frequency
       sampling.
     - :doc:`processing_agents`
   * - :ref:`DenoisingAgent <agent-denoising>`
     - Apply robust and AI-assisted denoising before analysis or inversion.
     - :doc:`processing_agents`

Forward And Inversion Workflows
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 52 20

   * - Agent
     - Main responsibility
     - Detail
   * - :ref:`ForwardModelAgent <agent-forward-model>`
     - Run 1-D, 2-D, or 3-D MT forward modelling from resistivity models.
     - :doc:`inversion_agents`
   * - :ref:`InversionPrepAgent <agent-inversion-prep>`
     - Prepare inversion input files through a general pre-inversion
       interface.
     - :doc:`inversion_agents`
   * - :ref:`Occam2DAgent <agent-occam2d>`
     - Write Occam2D data, mesh, startup, and related files.
     - :doc:`inversion_agents`
   * - :ref:`ModEmAgent <agent-modem>`
     - Write ModEM 3-D impedance data files.
     - :doc:`inversion_agents`
   * - :ref:`InversionBackendAgent <agent-inversion-backend>`
     - Drive pyCSAMT inversion backends from an agent workflow.
     - :doc:`inversion_agents`
   * - :ref:`InversionEvaluationAgent <agent-inversion-evaluation>`
     - Load inversion outputs and compute misfit, RMS, residual, and
       diagnostic products.
     - :doc:`inversion_agents`
   * - :ref:`InversionComparisonAgent <agent-inversion-comparison>`
     - Compare inversion sections, before/after results, or backend outputs.
     - :doc:`inversion_agents`

AI And Model-Zoo Agents
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 28 52 20

   * - Agent
     - Main responsibility
     - Detail
   * - :ref:`AIInversionAgent <agent-ai-inversion>`
     - Run end-to-end 1-D AI inversion from survey data.
     - :doc:`ai_model_zoo_agents`
   * - :ref:`Inv2DAgent <agent-inv2d>`
     - Run 2-D profile inversion with U-Net style models.
     - :doc:`ai_model_zoo_agents`
   * - :ref:`Inv3DAgent <agent-inv3d>`
     - Run 3-D spatial inversion using graph-based neural models.
     - :doc:`ai_model_zoo_agents`
   * - :ref:`EnsembleAgent <agent-ensemble>`
     - Estimate uncertainty using multiple inversion models or ensemble
       predictions.
     - :doc:`ai_model_zoo_agents`
   * - :ref:`JointInversionAgent <agent-joint-inversion>`
     - Run multi-modal inversion across MT, TEM, CSAMT, or other paired data.
     - :doc:`ai_model_zoo_agents`
   * - :ref:`AnomalyDetectionAgent <agent-anomaly-detection>`
     - Flag anomalous station-frequency samples or survey regions.
     - :doc:`ai_model_zoo_agents`
   * - :ref:`ModelZooAgent <agent-model-zoo>`
     - List, download, and use pre-trained EM inversion checkpoints.
     - :doc:`ai_model_zoo_agents`

Orchestration, Pipeline, And Outputs
------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 52 20

   * - Agent
     - Main responsibility
     - Detail
   * - :ref:`WorkflowOrchestratorAgent <agent-workflow-orchestrator>`
     - Classify natural-language requests and assemble matching workflows.
     - :doc:`orchestration_output_agents`
   * - :ref:`PipelineAgent <agent-pipeline>`
     - Bridge the pyCSAMT pipeline system into agent workflows.
     - :doc:`orchestration_output_agents`
   * - :ref:`BatchSurveyAgent <agent-batch-survey>`
     - Run workflows across multiple survey profiles or directories.
     - :doc:`orchestration_output_agents`
   * - :ref:`InterpretationAgent <agent-interpretation>`
     - Translate resistivity and inversion products into geological
       interpretation.
     - :doc:`orchestration_output_agents`
   * - :ref:`ResistivityMapAgent <agent-resistivity-map>`
     - Generate horizontal depth-slice maps from inversion outputs.
     - :doc:`orchestration_output_agents`
   * - :ref:`SensitivityAgent <agent-sensitivity>`
     - Estimate depth of investigation, vertical resolution, and sensitivity.
     - :doc:`orchestration_output_agents`
   * - :ref:`EDIExportAgent <agent-edi-export>`
     - Export processed survey objects back to EDI files.
     - :doc:`orchestration_output_agents`
   * - :ref:`ReportAgent <agent-report>`
     - Assemble workflow outputs into Markdown, HTML, or PDF reports.
     - :doc:`orchestration_output_agents`
   * - :ref:`CodeGenerationAgent <agent-code-generation>`
     - Generate a standalone Python script from a workflow configuration.
     - :doc:`orchestration_output_agents`

Choosing The Right Entry Point
------------------------------

Use ``ContextInputAgent`` when you need structured config from text.  Use
``MTLoaderAgent`` when the next step needs a ``Sites`` object.  Use
``AgentCoordinator`` when you know the workflow chain.  Use
``WorkflowOrchestratorAgent`` when the request itself should decide the chain.
Use direct processing agents when you are working interactively in a notebook
and want to inspect each intermediate result.

Related API
-----------

See :mod:`pycsamt.agents`, :doc:`../api/agents`, and the grouped pages linked
above.
