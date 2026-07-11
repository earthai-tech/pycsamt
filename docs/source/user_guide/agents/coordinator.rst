Agent Coordinator
=================

:class:`pycsamt.agents.AgentCoordinator` is the explicit workflow runner for
the agent layer.  It executes a named sequence of agent steps, passes outputs
from earlier steps into later steps, writes checkpoints, supports dry-run
previews, and returns one workflow-level :class:`pycsamt.agents.AgentResult`.

Use the coordinator when the chain is known and should be reproducible.  Use
:class:`~pycsamt.agents.WorkflowOrchestratorAgent` when a natural-language
request should decide which chain to build.


What the coordinator solves
---------------------------

Many pyCSAMT workflows are naturally staged:

.. code-block:: text

   load data
   -> quality control
   -> correction
   -> tensor diagnostics
   -> inversion preparation
   -> interpretation and report

Each stage can be represented by one agent.  The coordinator provides the
glue:

* a stable workflow name;
* an ordered list of named steps;
* ``input_fn`` callbacks that map previous results into the next agent input;
* optional steps that do not abort the workflow when they fail;
* per-step checkpoints and a workflow summary file;
* dry-run previews for inspection before data are processed;
* workflow-level elapsed time, warnings, and LLM cost aggregation.


Core objects
------------

``WorkflowStep``
    Internal descriptor for one step.  It stores the step name, agent
    instance, optional input mapping function, description, and whether the
    step is required.

``AgentCoordinator``
    Public runner that registers steps with ``add_step(...)`` and runs them
    with ``execute(...)`` or previews them with ``preview(...)``.

``AgentResult``
    Standard result object returned by each step and by the coordinator
    itself.  The coordinator result stores per-step results in
    ``result.data``.


Minimal workflow
----------------

The first step usually receives the top-level workflow configuration directly.
Later steps usually receive selected outputs from earlier results.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AgentCoordinator, DataQCAgent, MTLoaderAgent

   coord = AgentCoordinator("mt_qc")

   coord.add_step(
       "load",
       MTLoaderAgent(),
       description="Load EDI files into a Sites object.",
   )

   coord.add_step(
       "qc",
       DataQCAgent(),
       input_fn=lambda results: {
           "sites": results["load"]["sites"],
       },
       description="Compute survey and station QC.",
   )

   workflow = coord.execute({
       "path": "data/edis",
   })

   print(workflow.status)
   print(workflow.summary)
   print(workflow["qc"].summary)

The coordinator's ``data`` dictionary is keyed by step name.  Each value is
the :class:`~pycsamt.agents.AgentResult` returned by that step.


Step registration
-----------------

Register steps with :meth:`AgentCoordinator.add_step
<pycsamt.agents.AgentCoordinator.add_step>`.

.. code-block:: python
   :linenos:

   coord.add_step(
       name="static_shift",
       agent=StaticShiftAgent(method="ama"),
       input_fn=lambda results: {
           "sites": results["load"]["sites"],
           "output_dir": "outputs/static_shift",
       },
       description="Detect and correct static-shift effects.",
       required=True,
   )

The parameters are:

``name``
    Unique identifier for the step.  It is used in ``results[name]``,
    checkpoint filenames, logs, and preview output.

``agent``
    A :class:`~pycsamt.agents.BaseAgent` instance.  Construct the agent before
    registration so it resolves the desired LLM provider and model.

``input_fn``
    Optional callable receiving the accumulated previous step results.  It
    must return the input dictionary for the current agent.  When omitted, the
    coordinator passes the original workflow config directly.

``description``
    Human-readable action text shown in dry-run previews and progress output.

``required``
    ``True`` by default.  If a required step fails, the workflow aborts.  If an
    optional step fails, the coordinator records the failure and continues.

Step names must be unique.  Registering the same name twice raises
``ValueError``.


Input mapping with ``input_fn``
-------------------------------

``input_fn`` is the most important part of a coordinator workflow.  It is the
bridge between one agent's output contract and the next agent's input
contract.

The callable receives only the accumulated step results:

.. code-block:: python
   :linenos:

   def static_shift_input(results):
       return {
           "sites": results["load"]["sites"],
           "method": "ama",
           "output_dir": "outputs/static_shift",
       }

If a later step also needs values from the original top-level config, capture
them with a closure:

.. code-block:: python
   :linenos:

   workflow_config = {
       "path": "data/edis",
       "output_dir": "outputs/willy",
       "period_range": (0.001, 10.0),
   }

   def qc_input(results):
       return {
           "sites": results["load"]["sites"],
           "period_range": workflow_config["period_range"],
           "output_dir": f"{workflow_config['output_dir']}/qc",
       }

   coord.add_step("qc", DataQCAgent(), input_fn=qc_input)
   result = coord.execute(workflow_config)

If an ``input_fn`` raises an exception in a required step, the workflow returns
a failed ``AgentResult`` immediately.  For optional steps, the coordinator
records a warning and continues.


Dry-run preview
---------------

Use ``dry_run=True`` or call ``preview(...)`` to inspect the workflow before
running agents.

.. code-block:: python
   :linenos:

   preview = coord.execute(
       {"path": "data/edis"},
       dry_run=True,
   )

   print(preview["plan"])
   print(preview["steps"])

Equivalent explicit form:

.. code-block:: python
   :linenos:

   preview = coord.preview({"path": "data/edis"})

The preview includes:

* workflow name;
* number of registered steps;
* formatted input config;
* step order;
* agent class for each step;
* LLM provider/model or ``"no-LLM"``;
* whether the step is required;
* the step description.

No agent ``execute(...)`` method is called during preview.


Checkpoints and resume
----------------------

The coordinator writes checkpoints after each step.  By default, checkpoints
go to:

.. code-block:: text

   pycsamt_agent_checkpoints/<workflow_name>/

For each step, two files are written when possible:

``<step>.pkl``
    Pickled :class:`~pycsamt.agents.AgentResult`, used by ``resume=True``.

``<step>.json``
    Human-readable sidecar containing status, summary, elapsed time, cost,
    warnings, and error information.

At the end of the workflow, the coordinator also writes:

``workflow_state.json``
    Summary of the workflow and status/cost metadata for all recorded steps.

Set a custom checkpoint directory when constructing the coordinator:

.. code-block:: python
   :linenos:

   coord = AgentCoordinator(
       "willy_qc",
       checkpoint_dir="outputs/willy/checkpoints",
   )

Resume from existing checkpoints:

.. code-block:: python
   :linenos:

   result = coord.execute(
       {"path": "data/edis"},
       resume=True,
   )

Clear checkpoints for a workflow:

.. code-block:: python
   :linenos:

   coord.reset_checkpoints()

Use checkpoints for long-running or expensive workflows.  For quick notebook
experiments, you can leave the default directory or reset it between runs.


Required and optional steps
---------------------------

Required steps abort the workflow when they fail:

.. code-block:: python
   :linenos:

   coord.add_step(
       "qc",
       DataQCAgent(),
       input_fn=lambda r: {"sites": r["load"]["sites"]},
       required=True,
   )

Optional steps let the workflow continue:

.. code-block:: python
   :linenos:

   coord.add_step(
       "report",
       ReportAgent(formats=["md", "html"]),
       input_fn=lambda r: {
           "results": r,
           "output_dir": "outputs/report",
       },
       description="Build a report if report dependencies are available.",
       required=False,
   )

Use optional steps for report generation, extra plots, secondary exports, or
experimental analyses.  Keep core loading, QC, correction, and inversion-input
steps required when later steps depend on them.


Workflow result structure
-------------------------

``execute(...)`` returns one :class:`~pycsamt.agents.AgentResult`.

.. code-block:: python
   :linenos:

   result = coord.execute({"path": "data/edis"})

   print(result.status)
   print(result.summary)
   print(result.elapsed_seconds)
   print(result.cost_estimate_usd)
   print(result.warnings)

   load_result = result["load"]
   qc_result = result["qc"]

The coordinator result has:

``status``
    ``"success"`` when every recorded step succeeded,
    ``"needs_review"`` when at least one optional or non-aborting step did not
    succeed, and ``"failed"`` when a required step aborted the workflow.

``summary``
    Human-readable workflow completion or failure summary.

``data``
    Dictionary of ``step_name -> AgentResult``.

``warnings``
    Combined warnings from step results plus coordinator warnings.

``cost_estimate_usd``
    Sum of per-step LLM cost estimates for steps run in this execution.


Complete QC and correction example
----------------------------------

The example below builds a practical survey workflow with loading, QC,
static-shift correction, phase diagnostics, EDI export, and report generation.

.. code-block:: python
   :linenos:

   from pycsamt.agents import (
       AgentCoordinator,
       DataQCAgent,
       EDIExportAgent,
       MTLoaderAgent,
       PhaseAnalysisAgent,
       ReportAgent,
       StaticShiftAgent,
   )

   config = {
       "path": "data/edis",
       "output_dir": "outputs/willy",
   }

   coord = AgentCoordinator(
       "willy_qc_correction",
       checkpoint_dir=f"{config['output_dir']}/checkpoints",
   )

   coord.add_step(
       "load",
       MTLoaderAgent(),
       description="Load survey files.",
   )

   coord.add_step(
       "qc",
       DataQCAgent(),
       input_fn=lambda r: {
           "sites": r["load"]["sites"],
           "output_dir": f"{config['output_dir']}/qc",
       },
       description="Run data quality control.",
   )

   coord.add_step(
       "static_shift",
       StaticShiftAgent(method="ama"),
       input_fn=lambda r: {
           "sites": r["load"]["sites"],
           "output_dir": f"{config['output_dir']}/static_shift",
       },
       description="Correct static-shift effects.",
   )

   coord.add_step(
       "phase",
       PhaseAnalysisAgent(),
       input_fn=lambda r: {
           "sites": r["static_shift"]["corrected_sites"],
           "output_dir": f"{config['output_dir']}/phase",
       },
       description="Compute phase tensor and strike diagnostics.",
   )

   coord.add_step(
       "export_edi",
       EDIExportAgent(),
       input_fn=lambda r: {
           "sites": r["static_shift"]["corrected_sites"],
           "output_dir": f"{config['output_dir']}/edis",
       },
       description="Export corrected EDI files.",
   )

   coord.add_step(
       "report",
       ReportAgent(formats=["md", "html"]),
       input_fn=lambda r: {
           "results": r,
           "output_dir": f"{config['output_dir']}/report",
       },
       description="Assemble workflow report.",
       required=False,
   )

   preview = coord.preview(config)
   print(preview["plan"])

   result = coord.execute(config)
   print(result.summary)

The exact output keys depend on each agent.  When composing a new chain, check
the agent-specific group pages:

* :doc:`foundation_agents`
* :doc:`processing_agents`
* :doc:`inversion_agents`
* :doc:`ai_model_zoo_agents`
* :doc:`orchestration_output_agents`


Inversion preparation example
-----------------------------

Use the same pattern to prepare inversion files after correction and frequency
selection.

.. code-block:: python
   :linenos:

   from pycsamt.agents import (
       AgentCoordinator,
       DataQCAgent,
       FrequencyDecimationAgent,
       MTLoaderAgent,
       Occam2DAgent,
       StaticShiftAgent,
   )

   coord = AgentCoordinator("occam2d_preparation")

   coord.add_step("load", MTLoaderAgent())

   coord.add_step(
       "qc",
       DataQCAgent(),
       input_fn=lambda r: {"sites": r["load"]["sites"]},
   )

   coord.add_step(
       "static_shift",
       StaticShiftAgent(method="ama"),
       input_fn=lambda r: {"sites": r["load"]["sites"]},
   )

   coord.add_step(
       "decimate",
       FrequencyDecimationAgent(),
       input_fn=lambda r: {
           "sites": r["static_shift"]["corrected_sites"],
           "period_range": [0.001, 10.0],
           "n_per_decade": 6,
       },
   )

   coord.add_step(
       "occam2d",
       Occam2DAgent(),
       input_fn=lambda r: {
           "sites": r["static_shift"]["corrected_sites"],
           "period_range": [0.001, 10.0],
           "output_dir": "outputs/occam2d",
       },
   )

   result = coord.execute({"path": "data/edis"})


Custom agents
-------------

Custom agents inherit from :class:`pycsamt.agents.BaseAgent` and return
:class:`pycsamt.agents.AgentResult`.

.. code-block:: python
   :linenos:

   import time

   from pycsamt.agents import AgentResult, BaseAgent

   class SurveyNoteAgent(BaseAgent):
       SYSTEM_PROMPT = "You are a careful MT/CSAMT survey reviewer."

       def __init__(self, *, api_key=None, model=None, llm_provider="claude"):
           super().__init__(
               "SurveyNoteAgent",
               api_key=api_key,
               model=model,
               llm_provider=llm_provider,
               section_preset="pseudosection",
           )

       def execute(self, input_data):
           self._last_cost = 0.0
           t0 = time.time()

           qc_summary = input_data.get("qc_summary", "")
           note = self.query_llm(
               f"Write a concise survey QC note: {qc_summary}",
               max_tokens=250,
           )

           return AgentResult(
               status="success",
               summary="Survey note generated.",
               data={"note": note or "LLM note unavailable."},
               llm_interpretation=note,
               elapsed_seconds=time.time() - t0,
               cost_estimate_usd=self._last_cost,
           )

Register the custom step like any built-in agent:

.. code-block:: python
   :linenos:

   coord.add_step(
       "survey_note",
       SurveyNoteAgent(),
       input_fn=lambda r: {
           "qc_summary": r["qc"].summary,
       },
       description="Generate a concise QC note.",
       required=False,
   )


Best practices
--------------

* Give every step a short, stable, lowercase name such as ``"load"``,
  ``"qc"``, ``"static_shift"``, or ``"occam2d"``.
* Keep the first step simple.  Usually it loads data from the top-level
  ``config``.
* Treat ``input_fn`` as an explicit contract between steps.  Avoid hiding
  complex processing inside the callback.
* Preview workflows before long runs or file-writing steps.
* Use a custom ``checkpoint_dir`` for project workflows so outputs are not
  mixed across surveys.
* Mark reports, extra figures, and secondary exports as optional when they
  should not block numerical outputs.
* Configure :data:`pycsamt.agents.AGENT_CONFIG` before constructing agents if
  the workflow should use LLM assistance.
* Prefer deterministic/no-LLM mode for unit tests and documentation examples.


Common mistakes
---------------

.. list-table::
   :header-rows: 1
   :widths: 32 34 34

   * - Symptom
     - Cause
     - Fix
   * - ``KeyError`` in ``input_fn``
     - The previous agent did not return the expected data key, or the step
       name is wrong.
     - Inspect ``result["previous_step"].data`` and update the mapping.
   * - Later steps use raw data instead of corrected data
     - The ``input_fn`` still points to ``results["load"]["sites"]``.
     - Point it to the correction step output, for example
       ``results["static_shift"]["corrected_sites"]``.
   * - Workflow aborts at a reporting step
     - The report step is marked required.
     - Use ``required=False`` for non-critical deliverables.
   * - Resume does not skip a step
     - The checkpoint directory or workflow name changed.
     - Reuse the same ``workflow_name`` and ``checkpoint_dir``.
   * - Agents do not use the expected LLM provider
     - Agents were constructed before ``AGENT_CONFIG`` was configured.
     - Configure LLM settings first, then instantiate agents.


Related API
-----------

* :class:`pycsamt.agents.AgentCoordinator`
* :class:`pycsamt.agents.WorkflowStep`
* :class:`pycsamt.agents.BaseAgent`
* :class:`pycsamt.agents.AgentResult`
* :doc:`llm_configuration`
* :doc:`agent_catalogue`
:html_theme.sidebar_secondary.remove:
