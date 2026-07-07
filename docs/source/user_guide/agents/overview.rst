Agent Overview
==============

``pycsamt.agents`` is the AI-assisted workflow layer of pyCSAMT.  It turns
common MT, AMT, and CSAMT work into composable agents: one agent loads survey
files, another performs quality control, another prepares inversion files,
another writes a report, and so on.

The agent layer is not only an LLM wrapper.  Many agents have deterministic
processing paths and can run without an API key.  LLM support is added where it
is useful: parsing natural-language requests, selecting workflow steps,
generating interpretation text, and producing report narratives.

Core Ideas
----------

The package is built around four concepts:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Concept
     - Role
   * - :class:`pycsamt.agents.BaseAgent`
     - Abstract base class for every agent.  It provides LLM access, cost
       tracking, JSON extraction, plotting helpers, and common validation
       helpers.
   * - :class:`pycsamt.agents.AgentResult`
     - Standard return object with ``status``, ``summary``, ``data``,
       ``warnings``, ``llm_interpretation``, elapsed time, and cost estimate.
   * - :class:`pycsamt.agents.AgentCoordinator`
     - Ordered workflow executor.  It chains agents, maps outputs between
       steps, supports dry-run previews, and checkpoints workflow state.
   * - :class:`pycsamt.agents.WorkflowOrchestratorAgent`
     - High-level entry point for natural-language workflow requests.  It
       classifies the request and builds the matching agent chain.

Imports are lazy.  Importing :mod:`pycsamt.agents` exposes the public names
without importing optional provider libraries such as ``anthropic``,
``openai``, ``google-generativeai``, ``torch``, or ``gradio`` until they are
actually needed.

When To Use Agents
------------------

Use agents when you want a workflow-oriented interface rather than individual
low-level function calls.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Task
     - Recommended entry point
   * - Parse a plain-English request into a structured config
     - ``ContextInputAgent``
   * - Load EDI, AVG, or J files and inspect station completeness
     - ``MTLoaderAgent``
   * - Run a reproducible multi-step processing chain
     - ``AgentCoordinator``
   * - Let pyCSAMT choose the workflow from a request
     - ``WorkflowOrchestratorAgent``
   * - Prepare inversion inputs
     - ``InversionPrepAgent``, ``Occam2DAgent``, or ``ModEmAgent``
   * - Run AI inversion or model-zoo workflows
     - ``AIInversionAgent``, ``Inv2DAgent``, ``Inv3DAgent``,
       ``ModelZooAgent``
   * - Generate final products
     - ``ReportAgent``, ``CodeGenerationAgent``, ``EDIExportAgent``

Installation
------------

The base package can be installed without LLM, GPU, or web-interface
dependencies:

.. code-block:: bash
   :linenos:

   pip install pycsamt

Install optional provider clients only when you plan to use them:

.. code-block:: bash
   :linenos:

   pip install anthropic
   pip install openai
   pip install google-generativeai

Install optional AI or web dependencies only for those features:

.. code-block:: bash
   :linenos:

   pip install torch
   pip install gradio

The Agent Result Contract
-------------------------

Every agent returns :class:`pycsamt.agents.AgentResult`.  The object exposes
structured fields and also supports dict-like access to ``data``.

.. code-block:: python
   :linenos:

   from pycsamt.agents import MTLoaderAgent

   result = MTLoaderAgent().execute({
       "path": "/data/WILLY_EDIs",
   })

   if result:
       print(result.status)
       print(result.summary)
       print(result["sites"])
       print(result.get("station_names", []))
       print(result.cost_estimate_usd)
   else:
       print(result.error)
       print(result.error_fix_hint)

Important fields:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Field
     - Meaning
   * - ``status``
     - ``"success"``, ``"failed"``, or ``"needs_review"``.
   * - ``summary``
     - Short human-readable description of the run.
   * - ``data``
     - Agent-specific structured outputs, such as ``sites``, ``qc_table``,
       ``workflow_type``, paths, figures, metrics, or generated files.
   * - ``warnings``
     - Non-fatal issues encountered during execution.
   * - ``llm_interpretation``
     - Optional LLM-generated interpretation text.  ``None`` in no-LLM mode.
   * - ``elapsed_seconds``
     - Wall-clock runtime measured by the agent.
   * - ``cost_estimate_usd``
     - Estimated LLM cost for that agent run.
   * - ``error`` and ``error_fix_hint``
     - Failure details and suggested remediation.

No-LLM Workflow Preview
-----------------------

No-LLM mode is useful for CI, reproducible examples, and workflow previews.
The context agent can fall back to rule-based parsing, and the coordinator can
preview a plan without executing expensive steps.

.. code-block:: python
   :linenos:

   from pycsamt.agents import ContextInputAgent, MTLoaderAgent, AgentCoordinator

   context = ContextInputAgent()
   loader = MTLoaderAgent()

   coord = AgentCoordinator("preview_qc")
   coord.add_step(
       "parse",
       context,
       description="Parse natural-language request into a config",
   )
   coord.add_step(
       "load",
       loader,
       input_fn=lambda results: {
           "path": (results["parse"].get("config") or {}).get("data_path", "")
       },
       description="Load survey files and create a station QC table",
   )

   result = coord.execute(
       {"request": "Load /data/WILLY_EDIs, QC them, period 0.001 to 10 s"},
       dry_run=True,
   )

   print(result.status)
   print(result["plan"])
   print(result.cost_estimate_usd)

Dry-run results are still :class:`~pycsamt.agents.AgentResult` objects.  The
``data`` payload contains a formatted ``plan`` and a list of ``steps``.

Direct Agent Execution
----------------------

Agents can also be used individually.  This is the clearest pattern when you
already know which operation you want.

.. code-block:: python
   :linenos:

   from pycsamt.agents import MTLoaderAgent, DataQCAgent

   load = MTLoaderAgent().execute({
       "path": "/data/WILLY_EDIs",
   })

   if load.status == "success":
       qc = DataQCAgent().execute({
           "sites": load["sites"],
           "output_dir": "/out/willy_qc",
       })
       print(qc.summary)
       print(qc.warnings)

This style is explicit and easy to debug.  It is also useful in notebooks
where intermediate objects such as ``Sites`` should be inspected between
steps.

Coordinated Workflows
---------------------

Use :class:`pycsamt.agents.AgentCoordinator` when several agents should run as
a named workflow.  Each step can receive either the original configuration or
a transformed view of earlier step results.

.. code-block:: python
   :linenos:

   from pycsamt.agents import (
       AgentCoordinator,
       MTLoaderAgent,
       DataQCAgent,
       StaticShiftAgent,
       ReportAgent,
   )

   coord = AgentCoordinator("static_shift_report")
   coord.add_step("load", MTLoaderAgent(), description="Load survey files")
   coord.add_step(
       "qc",
       DataQCAgent(),
       input_fn=lambda r: {
           "sites": r["load"]["sites"],
           "output_dir": "/out/willy/qc",
       },
       description="Run quality control",
   )
   coord.add_step(
       "static_shift",
       StaticShiftAgent(),
       input_fn=lambda r: {
           "sites": r["load"]["sites"],
           "qc_table": r["qc"].get("qc_table"),
           "output_dir": "/out/willy/static_shift",
       },
       description="Detect and correct static shift",
   )
   coord.add_step(
       "report",
       ReportAgent(),
       input_fn=lambda r: {
           "workflow_results": r,
           "output_dir": "/out/willy/report",
       },
       description="Assemble the workflow report",
   )

   result = coord.execute({
       "path": "/data/WILLY_EDIs",
   })

   print(result.summary)
   print(result.cost_estimate_usd)

Natural-Language Orchestration
------------------------------

Use :class:`pycsamt.agents.WorkflowOrchestratorAgent` when a user request
should determine the workflow.  The orchestrator can classify requests such as
quality control, phase-tensor analysis, Occam2D preparation, ModEM
preparation, AI inversion, 2-D inversion, 3-D inversion, ensemble inversion,
joint inversion, and full workflows.

.. code-block:: python
   :linenos:

   from pycsamt.agents import WorkflowOrchestratorAgent

   agent = WorkflowOrchestratorAgent()

   result = agent.execute({
       "request": "Run phase tensor analysis and generate a report",
       "data_path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_phase",
       "dry_run": True,
   })

   print(result["workflow_type"])
   print(result["steps"])

Switch ``dry_run`` to ``False`` or omit it when the workflow should actually
run.

LLM-Assisted Interpretation
---------------------------

Agents inherit provider settings from :data:`pycsamt.agents.AGENT_CONFIG`
unless an explicit per-agent key is supplied.

.. code-block:: python
   :linenos:

   from pycsamt.agents import configure_agents
   from pycsamt.agents import PhaseAnalysisAgent

   configure_agents(
       provider="claude",
       api_key="sk-ant-...",
       model="claude-sonnet-4-6",
   )

   result = PhaseAnalysisAgent().execute({
       "path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_phase",
   })

   print(result.summary)
   print(result.llm_interpretation)
   print(result.cost_estimate_usd)

For provider selection, environment variables, budgets, and custom pricing,
see :doc:`llm_configuration`.

AI And Model-Zoo Entry Points
-----------------------------

AI inversion agents are available as first-class workflow steps and direct
interfaces.  They are intended for workflows where trained models or model-zoo
checkpoints are available.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AIInversionAgent, ModelZooAgent

   zoo = ModelZooAgent()
   models = zoo.execute({"action": "list"})
   print(models["models"])

   inverter = AIInversionAgent.from_pretrained("mt1d-resnet-5layer-v1")
   result = inverter.execute({
       "path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_ai",
   })

   print(result["rms_global"])

Use ``Inv2DAgent``, ``Inv3DAgent``, ``EnsembleAgent``, and
``JointInversionAgent`` for specialized deep-learning workflows.

CLI And Web Interface
---------------------

The agent package also exposes a module CLI:

.. code-block:: bash
   :linenos:

   python -m pycsamt.agents preview "Load /data/WILLY EDIs, QC, PT analysis"
   python -m pycsamt.agents list
   python -m pycsamt.agents pricing
   python -m pycsamt.agents zoo
   python -m pycsamt.agents web --port=7860

The web interface requires ``gradio`` and can also be started from Python:

.. code-block:: python
   :linenos:

   from pycsamt.agents.web import launch

   launch()

Outputs And Reproducibility
---------------------------

Agent workflows are designed to be inspectable.  A robust workflow should:

1. Keep the original input path and output directory in the workflow config.
2. Preserve intermediate ``AgentResult`` objects while debugging.
3. Use ``dry_run=True`` before expensive or file-writing runs.
4. Inspect ``warnings`` and ``error_fix_hint`` before trusting final products.
5. Track ``cost_estimate_usd`` per agent and session-wide spend through
   ``AGENT_CONFIG``.
6. Use ``CodeGenerationAgent`` when a workflow should be converted into a
   reproducible standalone Python script.

Where To Go Next
----------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Page
     - Use it for
   * - :doc:`llm_configuration`
     - Provider setup, key resolution, budgets, and pricing.
   * - :doc:`agent_catalogue`
     - Choosing the right agent for each task.
   * - :doc:`coordinator`
     - Building explicit multi-agent workflows.
   * - :doc:`orchestrator`
     - Natural-language workflow classification and execution.
   * - :doc:`../api/agents`
     - Generated API reference for ``pycsamt.agents`` modules.

Related API
-----------

See :mod:`pycsamt.agents`, :class:`pycsamt.agents.BaseAgent`,
:class:`pycsamt.agents.AgentResult`,
:class:`pycsamt.agents.AgentCoordinator`,
:class:`pycsamt.agents.WorkflowOrchestratorAgent`, and
:data:`pycsamt.agents.AGENT_CONFIG`.
