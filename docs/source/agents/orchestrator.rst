Workflow Orchestrator
=====================

:class:`pycsamt.agents.WorkflowOrchestratorAgent` is the high-level dispatcher
for agent workflows.  It accepts a natural-language request or an explicit
workflow configuration, chooses a workflow type, builds an
:class:`pycsamt.agents.AgentCoordinator`, and then either previews or executes
the selected chain.

This page is the detailed guide for the orchestrator itself.  The broader
catalogue page :doc:`orchestration_output_agents` lists the orchestrator
alongside related output agents such as ``PipelineAgent``, ``ReportAgent``,
``BatchSurveyAgent``, and ``CodeGenerationAgent``.

When To Use It
--------------

Use ``WorkflowOrchestratorAgent`` when the user can describe the objective
more naturally than they can name every processing step.

Good fits include:

.. list-table::
   :header-rows: 1
   :widths: 36 64

   * - Situation
     - Why the orchestrator helps
   * - Interactive notebooks
     - A short text request can become a reproducible agent chain.
   * - CLI or web interfaces
     - User text can be classified into an implemented workflow.
   * - Dry-run planning
     - The workflow can be previewed before reading, writing, or processing
       files.
   * - Repeated survey operations
     - The same high-level request can produce the same workflow structure.
   * - LLM-assisted workflows
     - When an API key is configured, the LLM can classify ambiguous requests
       before the rule-based fallback is used.

Use :class:`pycsamt.agents.AgentCoordinator` directly when you already know
the exact sequence of agents and want full control over step wiring.

How Routing Works
-----------------

The orchestrator follows a deterministic sequence:

1. Read ``request``, ``config``, ``dry_run``, ``output_dir``, and
   ``data_path`` from the input dictionary.
2. If ``config["workflow"]`` is provided, use it directly.
3. Otherwise, if an LLM API key is available, ask the LLM to return a JSON
   object containing ``workflow_type`` and ``reasoning``.
4. If no LLM result is available, classify the request with keyword rules.
5. If the workflow is still unknown, use ``default_workflow``.  The default is
   ``"qc"``.
6. Build an ``AgentCoordinator`` named ``orchestrated_<workflow_type>``.
7. Instantiate the agents needed by the selected workflow.
8. Execute the coordinator, or preview it when ``dry_run=True``.

The rule-based fallback is important.  It means the orchestrator is useful
even in no-LLM mode, especially for previews and deterministic tests.

Input Contract
--------------

The input to ``execute`` is a dictionary.

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - Key
     - Required
     - Meaning
   * - ``request``
     - Usually
     - Natural-language workflow request.  Used for LLM or keyword
       classification and path extraction.
   * - ``config``
     - No
     - Pre-built configuration.  If it contains ``workflow``, that workflow
       overrides natural-language classification.
   * - ``data_path``
     - Usually
     - Input EDI path or survey directory.  Overrides any path extracted from
       ``request``.
   * - ``output_dir``
     - No
     - Root directory for workflow outputs.  Defaults to
       ``"pycsamt_workflow_output"``.
   * - ``dry_run``
     - No
     - Preview the workflow without executing agent work.  Defaults to
       ``False``.

Output Contract
---------------

The orchestrator returns an :class:`pycsamt.agents.AgentResult`.  Its
``data`` dictionary contains:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Key
     - Meaning
   * - ``workflow_type``
     - Selected workflow name, such as ``"qc"`` or ``"phase_analysis"``.
   * - ``reasoning``
     - Short classification explanation from the LLM or keyword fallback.
   * - ``coordinator``
     - The built ``AgentCoordinator`` instance.
   * - ``result``
     - The result returned by the coordinator.
   * - ``steps``
     - List of step metadata dictionaries with step name, agent class, and
       description.

Basic Dry-Run Preview
---------------------

Start with dry-run mode whenever the request may write files, call LLMs, or
run expensive processing.

.. code-block:: python
   :linenos:

   from pycsamt.agents import WorkflowOrchestratorAgent

   agent = WorkflowOrchestratorAgent()

   result = agent.execute({
       "request": "QC the WILLY EDI files and prepare a short report",
       "data_path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_qc",
       "dry_run": True,
   })

   print(result.status)
   print(result["workflow_type"])
   print(result["reasoning"])

   for step in result["steps"]:
       print(step["name"], step["agent"], step["description"])

In dry-run mode, the returned coordinator result contains the plan and step
metadata.  The orchestrator result still reports the selected workflow and the
steps that would run.

Explicit Workflow Configuration
-------------------------------

If the workflow is already known, pass ``config={"workflow": ...}``.  This
skips classification and is the best choice for reproducible scripts.

.. code-block:: python
   :linenos:

   from pycsamt.agents import WorkflowOrchestratorAgent

   result = WorkflowOrchestratorAgent().execute({
       "config": {"workflow": "modem"},
       "data_path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_modem",
       "dry_run": True,
   })

   print(result["workflow_type"])
   print(result["steps"])

Use explicit configuration when the request text is user-facing prose but the
program already knows which workflow should run.

Supported Workflow Types
------------------------

The orchestrator has an internal workflow registry.  Each workflow maps to a
fixed list of agent steps.

.. list-table::
   :header-rows: 1
   :widths: 24 46 30

   * - Workflow
     - Typical chain
     - Keywords or use case
   * - ``qc``
     - ``MTLoaderAgent`` -> ``DataQCAgent`` -> ``StaticShiftAgent`` ->
       ``ReportAgent``
     - quality, cleaning, flagging, SNR, static shift
   * - ``phase_analysis``
     - ``MTLoaderAgent`` -> ``DataQCAgent`` -> ``StaticShiftAgent`` ->
       ``PhaseAnalysisAgent`` -> ``ReportAgent``
     - phase tensor, strike, dimensionality, Mohr, skew
   * - ``pre_inversion``
     - Load, QC, static shift, phase analysis, ``Occam2DAgent``,
       ``CodeGenerationAgent``
     - Occam2D, mesh, startup, inversion prep
   * - ``modem``
     - Load, QC, static shift, ``ModEmAgent``, report
     - ModEM, 3-D inversion, 3-D model
   * - ``ai_inversion``
     - Load, QC, ``DenoisingAgent``, ``AIInversionAgent``,
       ``InterpretationAgent``, report
     - neural, CNN, deep learning, AI inversion, Inv1D
   * - ``inv2d``
     - Load, QC, denoise, ``Inv2DAgent``, interpretation, report
     - U-Net, 2-D AI, profile inversion
   * - ``inv3d``
     - Load, QC, static shift, ``Inv3DAgent``, interpretation, report
     - GCN, graph convolutional, 3-D AI, spatial inversion
   * - ``ensemble_inversion``
     - Load, QC, denoise, ``EnsembleAgent``, interpretation, report
     - ensemble, uncertainty, confidence interval
   * - ``joint_inversion``
     - Load, QC, static shift, ``JointInversionAgent``, interpretation,
       report
     - joint inversion, multi-modal, MT+TEM, combined modality
   * - ``tipper``
     - Load, ``TipperAnalysisAgent``, report
     - tipper, induction arrows, Wiese, Parkinson
   * - ``sensitivity``
     - Load, QC, ``SensitivityAgent``, report
     - depth of investigation, DOI, resolution, sensitivity
   * - ``rotation``
     - Load, QC, phase analysis, ``TensorRotationAgent``
     - tensor rotation, strike rotation, coordinate frame
   * - ``full``
     - Load, QC, static shift, phase analysis, denoise, AI inversion,
       Occam2D preparation, report
     - full pipeline, complete processing, all steps

The keyword classifier also recognizes ``freq_decimation``, ``batch``, and
``comparison`` requests.  These are useful routing intents, but they should be
validated in dry-run mode before production use because they are more
specialized and depend on the availability of their downstream agents.

Running A Real Workflow
-----------------------

Remove ``dry_run`` when the workflow should execute.

.. code-block:: python
   :linenos:

   from pycsamt.agents import WorkflowOrchestratorAgent

   result = WorkflowOrchestratorAgent().execute({
       "request": "Run phase tensor analysis on the WILLY survey",
       "data_path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_phase",
   })

   print(result.summary)
   print(result["workflow_type"])
   print(result.cost_estimate_usd)

   coordinator_result = result["result"]
   print(coordinator_result.status)
   print(coordinator_result.summary)

For long or file-writing runs, keep the output directory stable.  The
orchestrator places coordinator checkpoints under
``<output_dir>/.checkpoints``.

LLM-Assisted Routing
--------------------

When the orchestrator is instantiated with an API key, it asks the configured
provider to classify the request.  The LLM must return JSON with
``workflow_type`` and ``reasoning``.  If this fails, the orchestrator falls
back to keyword classification.

.. code-block:: python
   :linenos:

   from pycsamt.agents import WorkflowOrchestratorAgent

   agent = WorkflowOrchestratorAgent(
       api_key="sk-ant-...",
       llm_provider="claude",
       model="claude-sonnet-4-6",
   )

   result = agent.execute({
       "request": (
           "The survey has noisy bands and I want an uncertainty-aware "
           "AI inversion report."
       ),
       "data_path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_uncertainty",
       "dry_run": True,
   })

   print(result["workflow_type"])
   print(result["reasoning"])

For global provider configuration, environment variables, budget caps, and
cost tracking, see :doc:`llm_configuration`.

No-LLM Routing
--------------

Without an API key, the orchestrator uses keyword matching.  This is the most
predictable mode for tests, examples, and CI.

.. code-block:: python
   :linenos:

   from pycsamt.agents import WorkflowOrchestratorAgent

   agent = WorkflowOrchestratorAgent()

   examples = [
       "run QC on the data",
       "compute phase tensor and strike analysis",
       "set up Occam2D mesh and startup file",
       "prepare ModEM 3D inversion",
       "ensemble uncertainty quantification",
   ]

   for request in examples:
       result = agent.execute({
           "request": request,
           "data_path": "/data/test",
           "dry_run": True,
       })
       print(request, "->", result["workflow_type"])

Handling Ambiguous Requests
---------------------------

By default, ambiguous text falls back to ``"qc"``.  You can change the
fallback workflow with ``default_workflow``.

.. code-block:: python
   :linenos:

   from pycsamt.agents import WorkflowOrchestratorAgent

   agent = WorkflowOrchestratorAgent(default_workflow="phase_analysis")

   result = agent.execute({
       "request": "prepare the usual diagnostic figures",
       "data_path": "/data/WILLY_EDIs",
       "dry_run": True,
   })

   print(result["workflow_type"])
   print(result.warnings)

Use explicit ``config={"workflow": ...}`` when ambiguity is not acceptable.

Inspecting The Planned Steps
----------------------------

The ``steps`` output is the easiest way to show a user what will happen.

.. code-block:: python
   :linenos:

   result = WorkflowOrchestratorAgent().execute({
       "request": "full pipeline for WILLY data",
       "data_path": "/data/WILLY_EDIs",
       "dry_run": True,
   })

   rows = []
   for index, step in enumerate(result["steps"], start=1):
       rows.append({
           "index": index,
           "name": step["name"],
           "agent": step["agent"],
           "description": step["description"],
       })

   for row in rows:
       print(
           f"{row['index']:02d}. {row['name']} "
           f"({row['agent']}): {row['description']}"
       )

This pattern is useful in a CLI, web UI, notebook, or report preamble.

Recommended Workflow
--------------------

For user-facing tools, use this operating pattern:

1. Run the orchestrator with ``dry_run=True``.
2. Display ``workflow_type``, ``reasoning``, and each step.
3. Let the user confirm the plan.
4. Re-run with the same ``request``, ``data_path``, and ``output_dir`` without
   ``dry_run``.
5. Inspect ``warnings`` and the nested coordinator ``result``.
6. Persist the output directory and generated reports.

Relationship To Other Pages
---------------------------

Use :doc:`orchestration_output_agents` for the catalogue of all orchestration
and output agents.  Use :doc:`coordinator` when you want to build the step
chain manually.  Use :doc:`agent_catalogue` when choosing the right agent
family.

Related API
-----------

See :class:`pycsamt.agents.WorkflowOrchestratorAgent`,
:class:`pycsamt.agents.AgentCoordinator`, :class:`pycsamt.agents.AgentResult`,
and :mod:`pycsamt.agents.orchestrator`.
