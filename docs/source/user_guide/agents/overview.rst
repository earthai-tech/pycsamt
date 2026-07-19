Agent Overview
==============

``pycsamt.agents`` is the workflow layer of pyCSAMT. It packages common
MT, AMT, and CSAMT tasks into small executable agents: one agent loads a
survey, another checks data quality, another corrects static shift, another
prepares inversion files, another runs AI inversion, and another assembles a
report. The value is not that every operation becomes automatic; the value is
that each operation has the same execution shape, so a workflow can be
previewed, run, inspected, checkpointed, and repeated.

The agent layer is broader than LLM usage. Many agents are deterministic and
run without an API key. LLM support is added where text understanding or
narrative drafting is useful: parsing natural-language requests, choosing a
workflow route, generating draft interpretation text, or writing report prose.
The scientific data products remain structured outputs inside
:term:`AgentResult`, and those outputs should be reviewed directly.


Core ideas
----------

The package is built around four concepts:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Concept
     - Role
   * - :class:`pycsamt.agents.BaseAgent`
     - Abstract base class for every agent. It provides LLM access, cost
       tracking, JSON extraction, plotting helpers, and common validation
       helpers.
   * - :term:`AgentResult`
     - Standard return object with ``status``, ``summary``, structured
       ``data``, ``warnings``, optional ``llm_interpretation``, elapsed time,
       and cost estimate.
   * - :term:`Agent coordinator`
     - Ordered workflow executor. It chains agents, maps outputs between
       steps, supports :term:`dry run` previews, and checkpoints workflow
       state.
   * - :class:`pycsamt.agents.WorkflowOrchestratorAgent`
     - High-level entry point for natural-language workflow requests. It
       classifies the request and builds the matching agent chain.

Imports are lazy. Importing :mod:`pycsamt.agents` exposes public names without
importing optional provider libraries such as ``anthropic``, ``openai``,
``google-generativeai``, ``torch``, or ``gradio`` until a feature needs them.
This keeps basic processing and no-LLM workflows lightweight.


How the pieces fit
------------------

At the lowest level, an agent is just a class with one ``execute`` method:

.. math::

   R = A(I),

where :math:`I` is the input dictionary, :math:`A` is the agent, and
:math:`R` is the returned :term:`AgentResult`. A coordinated workflow repeats
that pattern:

.. code-block:: text

   config
     -> MTLoaderAgent
     -> DataQCAgent
     -> StaticShiftAgent
     -> PhaseAnalysisAgent
     -> ReportAgent

The step boundary matters. Each step should make clear which object it
consumes: raw loaded sites, corrected sites, selected periods, an inversion
model, a QC table, or a dictionary of previous results. That explicit data
flow is what makes an agent chain easier to audit than a long notebook cell.


When to use agents
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
     - ``InversionPrepAgent``, ``Occam2DAgent``, ``ModEmAgent``, or
       ``Mare2DEMAgent``
   * - Run AI inversion or model-zoo workflows
     - ``AIInversionAgent``, ``Inv2DAgent``, ``Inv3DAgent``,
       ``EnsembleAgent``, ``JointInversionAgent``, ``ModelZooAgent``
   * - Generate final products
     - ``ReportAgent``, ``CodeGenerationAgent``, ``EDIExportAgent``

Use lower-level pyCSAMT APIs when you need full control of an algorithm's
internal arrays, solver configuration, training loop, or plotting object. Use
agents when the built-in input and output contracts match the workflow you
want to run.


Installation
------------

The base package can be installed without LLM, GPU, or web-interface
dependencies:

.. code-block:: bash

   pip install pycsamt

Install optional provider clients only when you plan to use them:

.. code-block:: bash

   pip install anthropic
   pip install openai
   pip install google-generativeai

Install optional AI or web dependencies only for those features:

.. code-block:: bash

   pip install torch
   pip install gradio

Optional dependencies are intentionally lazy. A user can run loading, QC,
classical preparation, and many previews without installing a neural backend
or an LLM client.


The AgentResult contract
------------------------

Every agent returns :term:`AgentResult`. The object exposes structured fields
and also supports dict-like access to ``data``.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import AgentResult
   >>> result = AgentResult(
   ...     status="needs_review",
   ...     summary="QC finished with warnings.",
   ...     data={"n_stations": 28},
   ...     warnings=["2 stations need review"],
   ... )
   >>> bool(result)
   True
   >>> result.status
   'needs_review'
   >>> result.get("n_stations")
   28
   >>> result.warnings[0]
   '2 stations need review'

Only ``"failed"`` is false in a boolean test. That means
``if result: ...`` allows both ``"success"`` and ``"needs_review"`` to pass.
Use exact status checks when a review state should stop a production workflow.

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
     - Optional LLM-generated interpretation text. ``None`` in no-LLM mode.
   * - ``elapsed_seconds``
     - Wall-clock runtime measured by the agent.
   * - ``cost_estimate_usd``
     - Estimated LLM cost for that agent run.
   * - ``error`` and ``error_fix_hint``
     - Failure details and suggested remediation.


No-LLM request parsing
----------------------

``ContextInputAgent`` converts a plain-English request into a structured
workflow configuration. With LLM credentials it can ask the configured model;
without credentials it uses deterministic request parsing for common workflow
phrases, paths, period ranges, components, output directories, and inversion
codes.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import AGENT_CONFIG, ContextInputAgent
   >>> with AGENT_CONFIG.offline():
   ...     result = ContextInputAgent().execute({
   ...         "request": (
   ...             "Run phase tensor analysis on /data/WILLY_EDIs, "
   ...             "period 0.001 to 10 s, save to /out/willy_phase"
   ...         )
   ...     })
   >>> print(result.status)
   success
   >>> print(result["config"]["workflow"])
   phase_analysis
   >>> print(result["config"]["data_path"])
   /data/WILLY_EDIs
   >>> print(result["config"]["period_range"])
   [0.001, 10.0]
   >>> print(result["config"]["output_dir"])
   /out/willy_phase
   >>> print(result.cost_estimate_usd)
   0.0

This result is not yet a processed survey; it is a plan-ready configuration.
The warnings list should still be checked, especially when the referenced path
does not exist on the current machine or the request lacks an output
directory.


Dry-run workflow preview
------------------------

A :term:`dry run` lets the coordinator describe the chain before executing
any step. It is useful before long runs, file-writing steps, AI training, or
operations that depend on external executables.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import AgentCoordinator, AgentResult, BaseAgent
   >>>
   >>> class MiniAgent(BaseAgent):
   ...     def __init__(self, name, payload):
   ...         super().__init__(name)
   ...         self.payload = payload
   ...
   ...     def execute(self, input_data):
   ...         return AgentResult(
   ...             status="success",
   ...             summary=f"{self.name} complete",
   ...             data={**self.payload, "input": input_data},
   ...         )
   >>>
   >>> coord = AgentCoordinator(
   ...     "overview_demo",
   ...     checkpoint_dir="outputs/docs/overview_demo",
   ...     verbose=False,
   ... )
   >>> coord.add_step(
   ...     "load",
   ...     MiniAgent("load", {"sites": "SITES:28"}),
   ...     description="Load survey",
   ... )
   >>> coord.add_step(
   ...     "qc",
   ...     MiniAgent("qc", {"qc_score": 0.91}),
   ...     input_fn=lambda r: {"sites": r["load"]["sites"]},
   ...     description="Run QC",
   ... )
   >>> preview = coord.execute({"path": "/data/WILLY_EDIs"}, dry_run=True)
   >>> print(preview.status)
   success
   >>> print(preview.summary)
   Workflow preview: 2 steps.
   >>> [(s["name"], s["agent"], s["required"]) for s in preview["steps"]]
   [('load', 'MiniAgent', True), ('qc', 'MiniAgent', True)]

The preview result contains a formatted ``plan`` string and a structured
``steps`` list. No registered agent's ``execute`` method is called.


Direct agent execution
----------------------

Agents can also be used individually. This is the clearest pattern when you
already know which operation you want and want to inspect each intermediate
object before moving on.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import MTLoaderAgent, DataQCAgent
   >>> load = MTLoaderAgent().execute({
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ... })
   >>> if load.status == "success":
   ...     qc = DataQCAgent().execute({
   ...         "sites": load["sites"],
   ...         "output_dir": "outputs/willy_qc",
   ...     })
   ...     print(qc.status)
   ...     print(qc.summary)

This style is especially useful in notebooks. The ``Sites`` object returned by
the loader can be inspected, plotted, filtered, or passed to a different
agent before any larger workflow is committed.


Coordinated workflows
---------------------

Use :class:`pycsamt.agents.AgentCoordinator` when several agents should run as
a named workflow. Each step can receive either the original configuration or a
mapped view of earlier step results.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import (
   ...     AgentCoordinator,
   ...     MTLoaderAgent,
   ...     DataQCAgent,
   ...     StaticShiftAgent,
   ...     ReportAgent,
   ... )
   >>> config = {
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "output_dir": "outputs/willy_static_shift",
   ... }
   >>> coord = AgentCoordinator(
   ...     "static_shift_report",
   ...     checkpoint_dir=f"{config['output_dir']}/checkpoints",
   ... )
   >>> coord.add_step("load", MTLoaderAgent(), description="Load survey files")
   >>> coord.add_step(
   ...     "qc",
   ...     DataQCAgent(),
   ...     input_fn=lambda r: {
   ...         "sites": r["load"]["sites"],
   ...         "output_dir": f"{config['output_dir']}/qc",
   ...     },
   ...     description="Run quality control",
   ... )
   >>> coord.add_step(
   ...     "static_shift",
   ...     StaticShiftAgent(),
   ...     input_fn=lambda r: {
   ...         "sites": r["load"]["sites"],
   ...         "qc_table": r["qc"].get("qc_table"),
   ...         "output_dir": f"{config['output_dir']}/static_shift",
   ...     },
   ...     description="Detect and correct static shift",
   ... )
   >>> coord.add_step(
   ...     "report",
   ...     ReportAgent(),
   ...     input_fn=lambda r: {
   ...         "workflow_results": r,
   ...         "output_dir": f"{config['output_dir']}/report",
   ...     },
   ...     description="Assemble the workflow report",
   ...     required=False,
   ... )
   >>> preview = coord.execute(config, dry_run=True)
   >>> print([step["name"] for step in preview["steps"]])
   ['load', 'qc', 'static_shift', 'report']

The key idea is the ``input_fn`` boundary. The QC step consumes loaded sites.
The static-shift step consumes loaded sites plus any QC table the QC agent
returned. The report step consumes the dictionary of previous step results.
See :doc:`coordinator` for the full workflow contract.


Natural-language orchestration
------------------------------

Use :class:`pycsamt.agents.WorkflowOrchestratorAgent` when a user request
should determine the workflow. The orchestrator can classify requests such as
quality control, phase-tensor analysis, Occam2D preparation, ModEM
preparation, AI inversion, 2-D inversion, 3-D inversion, ensemble inversion,
joint inversion, and full workflows.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import AGENT_CONFIG, WorkflowOrchestratorAgent
   >>> with AGENT_CONFIG.offline():
   ...     result = WorkflowOrchestratorAgent().execute({
   ...         "request": "Run phase tensor and strike analysis",
   ...         "data_path": "data/AMT/WILLY_DATA/L18PLT",
   ...         "output_dir": "outputs/willy_phase",
   ...         "dry_run": True,
   ...     })
   >>> print(result["workflow_type"])
   phase_analysis
   >>> print([step["name"] for step in result["steps"]])
   ['load', 'qc', 'static_shift', 'phase_analysis', 'report']

Remove ``dry_run`` when the chain should execute. In executed mode, the
orchestrator returns both the built coordinator and the coordinator result so
the workflow can still be inspected step by step.


LLM-assisted interpretation
---------------------------

Agents inherit provider settings from :data:`pycsamt.agents.AGENT_CONFIG`
unless an explicit per-agent key is supplied.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import configure_agents, PhaseAnalysisAgent
   >>> configure_agents(
   ...     provider="claude",
   ...     api_key="sk-ant-...",
   ...     model="claude-sonnet-4-6",
   ... )
   >>> result = PhaseAnalysisAgent().execute({
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "output_dir": "outputs/willy_phase",
   ... })
   >>> print(result.summary)
   >>> print(result.llm_interpretation)
   >>> print(result.cost_estimate_usd)

``llm_interpretation`` is narrative assistance, not a substitute for the
structured result fields, diagnostic figures, or a geophysicist's review. For
provider selection, environment variables, budgets, and custom pricing, see
:doc:`llm_configuration`.


AI and model-zoo entry points
-----------------------------

AI inversion agents are available as first-class workflow steps and direct
interfaces. They are intended for workflows where a neural inverse model,
synthetic training configuration, or model-zoo checkpoint is part of the
analysis.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import ModelZooAgent
   >>> zoo = ModelZooAgent()
   >>> models = zoo.execute({"action": "list"})
   >>> print(models.status)
   success
   >>> print(models.summary)
   5 pre-trained models in zoo.
   >>> [item["name"] for item in models["details"][:2]]
   ['mt1d-resnet-5layer-v1', 'mt1d-cnn-5layer-v1']

Use ``Inv2DAgent``, ``Inv3DAgent``, ``EnsembleAgent``, and
``JointInversionAgent`` for specialized deep-learning workflows. See
:doc:`ai_model_zoo_agents` for the agent catalogue view and
:doc:`../ai_inversion/agents` for the detailed AI inversion workflow guide.


CLI and web interface
---------------------

The agent package also exposes a module CLI:

.. code-block:: bash

   python -m pycsamt.agents preview "Load /data/WILLY EDIs, QC, PT analysis"
   python -m pycsamt.agents list
   python -m pycsamt.agents pricing
   python -m pycsamt.agents zoo
   python -m pycsamt.agents web --port=7860

The web interface requires ``gradio`` and can also be started from Python:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents.web import launch
   >>> launch()

The CLI is useful for quick previews, catalogue checks, pricing checks, and
model-zoo inspection. The Python API remains the preferred interface for
auditable project workflows because it keeps structured results in memory and
lets the caller archive them explicitly.


Outputs and reproducibility
---------------------------

Agent workflows are designed to be inspectable. A robust workflow should:

1. Keep the original input path and output directory in the workflow config.
2. Preserve intermediate :term:`AgentResult` objects while debugging.
3. Use ``dry_run=True`` before expensive or file-writing runs.
4. Inspect ``warnings`` and ``error_fix_hint`` before trusting final products.
5. Track ``cost_estimate_usd`` per agent and session-wide spend through
   ``AGENT_CONFIG``.
6. Record checkpoint directories for coordinated workflows.
7. Use ``CodeGenerationAgent`` when a workflow should be converted into a
   reproducible standalone Python script.

For scientific results, reproducibility also means recording which object a
step consumed: raw survey, corrected survey, selected-period survey, AI
prediction, conventional inversion output, or full previous-step dictionary.
That record is what lets another user distinguish "the report was generated"
from "the report was generated from the reviewed corrected data."


Where to go next
----------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Page
     - Use it for
   * - :doc:`agent_catalogue`
     - Choosing the right agent for each task.
   * - :doc:`coordinator`
     - Building explicit multi-agent workflows.
   * - :doc:`orchestrator`
     - Natural-language workflow classification and execution.
   * - :doc:`llm_configuration`
     - Provider setup, key resolution, budgets, and pricing.
   * - :doc:`foundation_agents`
     - BaseAgent, AgentResult, ContextInputAgent, MTLoaderAgent, and
       AgentCoordinator basics.
   * - :doc:`processing_agents`
     - QC, static shift, phase analysis, tensor rotation, tipper analysis,
       frequency decimation, and denoising.
   * - :doc:`inversion_agents`
     - Forward modelling and classical inversion preparation.
   * - :doc:`ai_model_zoo_agents`
     - AI inversion, model-zoo, ensemble, joint inversion, and anomaly agents.


Related API
-----------

See :mod:`pycsamt.agents`, :class:`pycsamt.agents.BaseAgent`,
:class:`pycsamt.agents.AgentResult`,
:class:`pycsamt.agents.AgentCoordinator`,
:class:`pycsamt.agents.WorkflowOrchestratorAgent`, and
:data:`pycsamt.agents.AGENT_CONFIG`.
