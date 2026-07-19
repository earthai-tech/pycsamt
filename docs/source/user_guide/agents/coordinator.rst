Agent Coordinator
=================

:class:`pycsamt.agents.AgentCoordinator` is the explicit workflow runner for
the agent layer. Use it when the chain is already known, when every transition
between steps should be visible, and when a run must leave enough state for a
reviewer to understand what happened without re-reading a notebook. It
registers named :term:`workflow step`\ s, executes them in order, passes
selected outputs forward through :term:`input mapping`, writes
:term:`workflow checkpoint`\ s, supports dry-run previews, aggregates warnings
and LLM cost, and returns one workflow-level :term:`AgentResult`.

The coordinator is deliberately less magical than
:class:`~pycsamt.agents.WorkflowOrchestratorAgent`. The orchestrator reads a
natural-language request and decides which chain to build. The coordinator is
for the moment after that decision: the project has a known processing order,
and the order itself is part of the reproducible record.


What the coordinator solves
---------------------------

Many pyCSAMT workflows are staged rather than monolithic:

.. code-block:: text

   load data
   -> quality control
   -> correction
   -> tensor diagnostics
   -> inversion preparation
   -> interpretation and report

Each stage has its own scientific contract. Loading returns station objects,
:term:`quality control` returns warnings and tables, correction returns a new
or modified survey object, inversion preparation writes solver files, and
reporting consumes the reviewed results. The coordinator does not blur those
contracts into one large function. It keeps a dictionary of named
:term:`AgentResult`\ s so the run can be inspected step by step:

.. math::

   R_k = A_k(I_k), \qquad
   I_k =
   \begin{cases}
      C, & \text{if no input mapping is supplied},\\
      g_k(R_1,\ldots,R_{k-1}), & \text{otherwise},
   \end{cases}

where :math:`C` is the top-level workflow configuration, :math:`A_k` is the
agent at step :math:`k`, :math:`g_k` is the step's ``input_fn``, and
:math:`R_k` is the returned :term:`AgentResult`. This small formulation is the
heart of reproducibility: a downstream step should receive exactly the fields
you name, not an implicit grab bag of global state.

The coordinator provides:

* a stable workflow name;
* an ordered list of named steps;
* explicit ``input_fn`` callbacks for step-to-step data flow;
* required and optional failure behavior;
* per-step checkpoints and a workflow summary file;
* dry-run previews for inspection before data are processed;
* workflow-level elapsed time, warnings, and LLM cost aggregation.


Core objects
------------

``WorkflowStep``
    Internal descriptor for one step. It stores the step name, agent instance,
    optional input mapping function, description, and whether the step is
    required.

``AgentCoordinator``
    Public runner that registers steps with ``add_step(...)`` and runs them
    with ``execute(...)`` or previews them with ``preview(...)``.

``AgentResult``
    Standard result object returned by each step and by the coordinator
    itself. The coordinator result stores per-step results in ``result.data``.

The workflow-level status is derived from the step statuses:

.. math::

   \mathrm{status}_{workflow} =
   \begin{cases}
      \mathrm{failed}, & \text{a required step fails},\\
      \mathrm{needs\_review}, & \text{at least one recorded step is not success},\\
      \mathrm{success}, & \text{all recorded steps succeed}.
   \end{cases}

This means ``success`` is still an execution statement, not a scientific
approval. A workflow can complete cleanly and still require geological,
statistical, or deliverable review.


Minimal workflow
----------------

The first step usually receives the top-level workflow configuration directly.
Later steps usually receive selected outputs from earlier results. The example
below uses tiny custom agents to expose the coordinator mechanics without
mixing in survey loading, plotting, or optional AI dependencies.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import AgentCoordinator, AgentResult, BaseAgent
   >>>
   >>> class ConstantAgent(BaseAgent):
   ...     def __init__(self, name, **payload):
   ...         super().__init__(name)
   ...         self.payload = payload
   ...
   ...     def execute(self, input_data):
   ...         return AgentResult(
   ...             status="success",
   ...             summary=f"{self.name} accepted {sorted(input_data)}.",
   ...             data=dict(self.payload, received=dict(input_data)),
   ...         )
   >>>
   >>> coord = AgentCoordinator(
   ...     "doc_demo",
   ...     checkpoint_dir="outputs/docs/doc_demo_checkpoints",
   ...     verbose=False,
   ... )
   >>> coord.add_step(
   ...     "load",
   ...     ConstantAgent("load", sites="SITES:3", n_stations=3),
   ...     description="Load EDI files into a Sites object.",
   ... )
   >>> coord.add_step(
   ...     "qc",
   ...     ConstantAgent("qc", qc_score=0.94),
   ...     input_fn=lambda r: {
   ...         "sites": r["load"]["sites"],
   ...         "min_score": 0.8,
   ...     },
   ...     description="Compute survey and station QC.",
   ... )
   >>> preview = coord.preview({"path": "data/AMT/WILLY_DATA/L18PLT"})
   >>> print(preview.status)
   success
   >>> print(preview.summary)
   Workflow preview: 2 steps.
   >>> print([step["name"] for step in preview["steps"]])
   ['load', 'qc']
   >>> result = coord.execute({"path": "data/AMT/WILLY_DATA/L18PLT"})
   >>> print(result.status)
   success
   >>> print(result["qc"].summary)
   qc accepted ['min_score', 'sites'].
   >>> print(result["qc"]["received"])
   {'sites': 'SITES:3', 'min_score': 0.8}

The important detail is not the fake ``sites`` value; it is the data flow.
The ``load`` step received the original config. The ``qc`` step received only
the two fields returned by its ``input_fn``. That is the pattern to preserve in
real chains: keep the first step simple, then make each transition explicit.


Step registration
-----------------

Register steps with :meth:`AgentCoordinator.add_step
<pycsamt.agents.AgentCoordinator.add_step>`.

.. code-block:: pycon
   :linenos:

   >>> coord.add_step(
   ...     name="static_shift",
   ...     agent=StaticShiftAgent(method="ama"),
   ...     input_fn=lambda results: {
   ...         "sites": results["load"]["sites"],
   ...         "output_dir": "outputs/static_shift",
   ...     },
   ...     description="Detect and correct static-shift effects.",
   ...     required=True,
   ... )

The parameters are:

``name``
    Unique identifier for the step. It is used in ``results[name]``,
    checkpoint filenames, logs, preview output, and report labels. Prefer
    short stable names such as ``"load"``, ``"qc"``, ``"static_shift"``, or
    ``"occam2d"``. Renaming a step changes the resume key.

``agent``
    A :class:`~pycsamt.agents.BaseAgent` instance. Construct the agent before
    registration so LLM provider, model, plotting preset, and constructor
    parameters are fixed before the run begins.

``input_fn``
    Optional callable receiving the accumulated previous step results. It must
    return the input dictionary for the current agent. When omitted, the
    coordinator passes the original workflow config directly.

``description``
    Human-readable action text shown in dry-run previews and progress output.

``required``
    ``True`` by default. If a required step fails, the workflow aborts. If an
    optional step fails, the coordinator records the failure and continues.

Step names must be unique. Registering the same name twice raises
``ValueError``.


Input mapping with ``input_fn``
-------------------------------

``input_fn`` is the most important part of a coordinator workflow. It is the
bridge between one agent's output contract and the next agent's input
contract. It should be short enough that a reviewer can read it as a wiring
diagram, not a hidden processing step.

.. code-block:: pycon
   :linenos:

   >>> def static_shift_input(results):
   ...     return {
   ...         "sites": results["load"]["sites"],
   ...         "method": "ama",
   ...         "output_dir": "outputs/static_shift",
   ...     }

If a later step also needs values from the original top-level config, keep
that config in the surrounding scope and read the required values inside the
mapping function:

.. code-block:: pycon
   :linenos:

   >>> workflow_config = {
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "output_dir": "outputs/willy",
   ...     "period_range": (0.001, 10.0),
   ... }
   >>>
   >>> def qc_input(results):
   ...     return {
   ...         "sites": results["load"]["sites"],
   ...         "period_range": workflow_config["period_range"],
   ...         "output_dir": f"{workflow_config['output_dir']}/qc",
   ...     }
   >>>
   >>> coord.add_step("qc", DataQCAgent(), input_fn=qc_input)

For reproducibility, keep these rules in mind:

* map from documented output keys, not from private attributes;
* pass corrected data to later steps only when the correction step is the
  intended source;
* keep constants such as ``period_range`` and output folders in the workflow
  config, not scattered across callbacks;
* inspect ``result["step_name"].data`` before wiring a new downstream step.

If an ``input_fn`` raises an exception in a required step, the workflow returns
a failed :term:`AgentResult` immediately. For optional steps, the coordinator
records a warning and continues.


Dry-run preview
---------------

Use ``dry_run=True`` or call ``preview(...)`` to inspect the workflow before
running agents.

.. code-block:: pycon
   :linenos:

   >>> preview = coord.execute(
   ...     {"path": "data/AMT/WILLY_DATA/L18PLT"},
   ...     dry_run=True,
   ... )
   >>> print(preview["steps"][0]["name"])
   load
   >>> print(preview["steps"][0]["llm"])
   no-LLM

Equivalent explicit form:

.. code-block:: pycon
   :linenos:

   >>> preview = coord.preview({"path": "data/AMT/WILLY_DATA/L18PLT"})

The preview includes workflow name, number of steps, formatted input config,
step order, agent class, LLM provider/model or ``"no-LLM"``, whether the step
is required, and the step description. No agent ``execute(...)`` method is
called during preview, so it is safe for expensive, file-writing, or
network-dependent chains.


Checkpoints and resume
----------------------

The coordinator writes checkpoints after each executed step. By default,
checkpoints go to:

.. code-block:: text

   pycsamt_agent_checkpoints/<workflow_name>/

For each step, two files are written when possible:

``<step>.pkl``
    Pickled :class:`~pycsamt.agents.AgentResult`, used by ``resume=True``.
    The checkpoint copy drops Matplotlib figures and containers of figures
    before pickling, because figure objects often contain unpicklable display
    closures. Keep figure files in ``output_dir`` for reporting.

``<step>.json``
    Human-readable sidecar containing status, summary, elapsed time, cost,
    warnings, and error information.

At the end of the workflow, the coordinator also writes:

``workflow_state.json``
    Summary of the workflow and status/cost metadata for all recorded steps.

Set a custom checkpoint directory when constructing the coordinator:

.. code-block:: pycon
   :linenos:

   >>> coord = AgentCoordinator(
   ...     "willy_qc",
   ...     checkpoint_dir="outputs/willy/checkpoints",
   ... )

Resume from existing checkpoints:

.. code-block:: pycon
   :linenos:

   >>> result = coord.execute(
   ...     {"path": "data/AMT/WILLY_DATA/L18PLT"},
   ...     resume=True,
   ... )

Resume is name-based. A step is skipped only when
``<checkpoint_dir>/<step_name>.pkl`` exists and can be loaded. If the workflow
name, checkpoint directory, or step name changes, the coordinator treats the
step as new work.

Clear checkpoints for a workflow:

.. code-block:: pycon
   :linenos:

   >>> coord.reset_checkpoints()

Use checkpoints for long-running or expensive workflows. For quick notebook
experiments, leave the default directory or reset it between runs. For formal
deliverables, checkpoints are only execution state; archive exported arrays,
figures, configuration files, and reports separately.


Required and optional steps
---------------------------

Required steps abort the workflow when they fail. The workflow result keeps
the partial results produced before the abort, which is important for
debugging, review, and restart decisions.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import AgentCoordinator, AgentResult, BaseAgent
   >>>
   >>> class FakeAgent(BaseAgent):
   ...     def __init__(self, name, result):
   ...         super().__init__(name)
   ...         self.result = result
   ...
   ...     def execute(self, input_data):
   ...         return self.result
   >>>
   >>> ok = AgentResult(
   ...     status="success",
   ...     summary="load ok",
   ...     data={"sites": "SITES:3"},
   ... )
   >>> bad = AgentResult.failed(
   ...     "qc failed",
   ...     hint="inspect station errors",
   ... )
   >>> coord = AgentCoordinator(
   ...     "doc_required_failure",
   ...     checkpoint_dir="outputs/docs/required_failure",
   ...     verbose=False,
   ... )
   >>> coord.add_step("load", FakeAgent("load", ok))
   >>> coord.add_step(
   ...     "qc",
   ...     FakeAgent("qc", bad),
   ...     input_fn=lambda r: {"sites": r["load"]["sites"]},
   ... )
   >>> coord.add_step(
   ...     "report",
   ...     FakeAgent("report", AgentResult("success", "report ok")),
   ... )
   >>> result = coord.execute({"path": "data/edis"})
   >>> print(result.status)
   failed
   >>> print(result.error)
   qc failed
   >>> print(result.error_fix_hint)
   inspect station errors
   >>> print(sorted(result.data))
   ['load', 'qc']

Optional steps let the workflow continue, but the workflow is marked
``needs_review`` because not all recorded steps succeeded.

.. code-block:: pycon
   :linenos:

   >>> coord = AgentCoordinator(
   ...     "doc_optional_failure",
   ...     checkpoint_dir="outputs/docs/optional_failure",
   ...     verbose=False,
   ... )
   >>> coord.add_step(
   ...     "optional_plot",
   ...     FakeAgent("optional", AgentResult.failed("optional plot failed")),
   ...     required=False,
   ... )
   >>> coord.add_step(
   ...     "report",
   ...     FakeAgent("report", AgentResult("success", "report ok")),
   ... )
   >>> result = coord.execute({"path": "data/edis"})
   >>> print(result.status)
   needs_review
   >>> print(result["optional_plot"].status)
   failed
   >>> print(result["report"].status)
   success

Use optional steps for reports, secondary figures, extra exports, and
experimental analyses. Keep loading, QC, correction, and inversion-input steps
required when later steps depend on them.


Workflow result structure
-------------------------

``execute(...)`` returns one :class:`~pycsamt.agents.AgentResult`.

.. code-block:: pycon
   :linenos:

   >>> result = coord.execute({"path": "data/AMT/WILLY_DATA/L18PLT"})
   >>> print(result.status)
   success
   >>> print(result.summary)
   Workflow 'doc_demo' complete: 2/2 steps succeeded in 0.5s ($0.000000).
   >>> print(result.elapsed_seconds >= 0)
   True
   >>> print(result.cost_estimate_usd)
   0.0
   >>> load_result = result["load"]
   >>> qc_result = result["qc"]

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

The workflow-level elapsed time is wall-clock time around the coordinator
execution. Per-step elapsed times remain inside each nested step result.


Complete QC and correction example
----------------------------------

The example below builds a practical survey workflow with loading, QC,
static-shift correction, phase diagnostics, EDI export, and report generation.
It is intentionally explicit: every downstream step states whether it consumes
raw loaded sites or corrected sites.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import (
   ...     AgentCoordinator,
   ...     DataQCAgent,
   ...     EDIExportAgent,
   ...     MTLoaderAgent,
   ...     PhaseAnalysisAgent,
   ...     ReportAgent,
   ...     StaticShiftAgent,
   ... )
   >>>
   >>> config = {
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "output_dir": "outputs/willy",
   ... }
   >>>
   >>> coord = AgentCoordinator(
   ...     "willy_qc_correction",
   ...     checkpoint_dir=f"{config['output_dir']}/checkpoints",
   ... )
   >>> coord.add_step("load", MTLoaderAgent(), description="Load survey files.")
   >>> coord.add_step(
   ...     "qc",
   ...     DataQCAgent(),
   ...     input_fn=lambda r: {
   ...         "sites": r["load"]["sites"],
   ...         "output_dir": f"{config['output_dir']}/qc",
   ...     },
   ...     description="Run data quality control.",
   ... )
   >>> coord.add_step(
   ...     "static_shift",
   ...     StaticShiftAgent(method="ama"),
   ...     input_fn=lambda r: {
   ...         "sites": r["load"]["sites"],
   ...         "output_dir": f"{config['output_dir']}/static_shift",
   ...     },
   ...     description="Correct static-shift effects.",
   ... )
   >>> coord.add_step(
   ...     "phase",
   ...     PhaseAnalysisAgent(),
   ...     input_fn=lambda r: {
   ...         "sites": r["static_shift"]["corrected_sites"],
   ...         "output_dir": f"{config['output_dir']}/phase",
   ...     },
   ...     description="Compute phase tensor and strike diagnostics.",
   ... )
   >>> coord.add_step(
   ...     "export_edi",
   ...     EDIExportAgent(),
   ...     input_fn=lambda r: {
   ...         "sites": r["static_shift"]["corrected_sites"],
   ...         "output_dir": f"{config['output_dir']}/edis",
   ...     },
   ...     description="Export corrected EDI files.",
   ... )
   >>> coord.add_step(
   ...     "report",
   ...     ReportAgent(formats=["md", "html"]),
   ...     input_fn=lambda r: {
   ...         "results": r,
   ...         "output_dir": f"{config['output_dir']}/report",
   ...     },
   ...     description="Assemble workflow report.",
   ...     required=False,
   ... )
   >>> preview = coord.preview(config)
   >>> print(preview.summary)
   Workflow preview: 6 steps.

On a real survey, inspect the generated QC tables, static-shift warnings,
phase-diagnostic figures, exported EDI paths, and report result separately.
The coordinator summary tells you whether the chain completed; the nested step
results tell you whether each scientific product is suitable for the next
stage.

The exact output keys depend on each agent. When composing a new chain, check
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

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import (
   ...     AgentCoordinator,
   ...     DataQCAgent,
   ...     FrequencyDecimationAgent,
   ...     MTLoaderAgent,
   ...     Occam2DAgent,
   ...     StaticShiftAgent,
   ... )
   >>>
   >>> config = {
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "period_range": [0.001, 10.0],
   ...     "output_dir": "outputs/occam2d",
   ... }
   >>> coord = AgentCoordinator("occam2d_preparation")
   >>> coord.add_step("load", MTLoaderAgent())
   >>> coord.add_step(
   ...     "qc",
   ...     DataQCAgent(),
   ...     input_fn=lambda r: {"sites": r["load"]["sites"]},
   ... )
   >>> coord.add_step(
   ...     "static_shift",
   ...     StaticShiftAgent(method="ama"),
   ...     input_fn=lambda r: {"sites": r["load"]["sites"]},
   ... )
   >>> coord.add_step(
   ...     "decimate",
   ...     FrequencyDecimationAgent(),
   ...     input_fn=lambda r: {
   ...         "sites": r["static_shift"]["corrected_sites"],
   ...         "period_range": config["period_range"],
   ...         "n_per_decade": 6,
   ...     },
   ... )
   >>> coord.add_step(
   ...     "occam2d",
   ...     Occam2DAgent(),
   ...     input_fn=lambda r: {
   ...         "sites": r["static_shift"]["corrected_sites"],
   ...         "period_range": config["period_range"],
   ...         "output_dir": config["output_dir"],
   ...     },
   ... )
   >>> preview = coord.execute(config, dry_run=True)
   >>> print([step["name"] for step in preview["steps"]])
   ['load', 'qc', 'static_shift', 'decimate', 'occam2d']

Notice that ``occam2d`` uses the corrected sites and the same period range
stored in ``config``. If frequency decimation returns a new site object or a
selected-period table in your workflow, wire that exact output key into the
Occam step instead of relying on the older config value.


Custom agents
-------------

Custom agents inherit from :class:`pycsamt.agents.BaseAgent` and return
:class:`pycsamt.agents.AgentResult`. They are useful for project-specific
checks, manifests, report annotations, or validation gates that are too local
to belong in pyCSAMT core.

.. code-block:: pycon
   :linenos:

   >>> import time
   >>> from pycsamt.agents import AgentResult, BaseAgent
   >>>
   >>> class SurveyNoteAgent(BaseAgent):
   ...     SYSTEM_PROMPT = "You are a careful MT/CSAMT survey reviewer."
   ...
   ...     def __init__(self, *, api_key=None, model=None, llm_provider="claude"):
   ...         super().__init__(
   ...             "SurveyNoteAgent",
   ...             api_key=api_key,
   ...             model=model,
   ...             llm_provider=llm_provider,
   ...             section_preset="pseudosection",
   ...         )
   ...
   ...     def execute(self, input_data):
   ...         self._last_cost = 0.0
   ...         t0 = time.time()
   ...         qc_summary = input_data.get("qc_summary", "")
   ...         note = self.query_llm(
   ...             f"Write a concise survey QC note: {qc_summary}",
   ...             max_tokens=250,
   ...         )
   ...         return AgentResult(
   ...             status="success",
   ...             summary="Survey note generated.",
   ...             data={"note": note or "LLM note unavailable."},
   ...             llm_interpretation=note,
   ...             elapsed_seconds=time.time() - t0,
   ...             cost_estimate_usd=self._last_cost,
   ...         )

Register the custom step like any built-in agent:

.. code-block:: pycon
   :linenos:

   >>> coord.add_step(
   ...     "survey_note",
   ...     SurveyNoteAgent(),
   ...     input_fn=lambda r: {
   ...         "qc_summary": r["qc"].summary,
   ...     },
   ...     description="Generate a concise QC note.",
   ...     required=False,
   ... )

For production workflows, keep deterministic validation gates separate from
optional LLM narrative steps. The workflow can then fail on a numerical or
data-quality condition while still treating generated prose as an optional
review aid.


Best practices
--------------

* Give every step a short, stable, lowercase name such as ``"load"``,
  ``"qc"``, ``"static_shift"``, or ``"occam2d"``.
* Keep the first step simple. Usually it loads data from the top-level
  ``config``.
* Treat ``input_fn`` as an explicit contract between steps. Avoid hiding
  complex processing inside the callback.
* Preview workflows before long runs or file-writing steps.
* Use a custom ``checkpoint_dir`` for project workflows so outputs are not
  mixed across surveys.
* Mark reports, extra figures, and secondary exports as optional when they
  should not block numerical outputs.
* Configure :data:`pycsamt.agents.AGENT_CONFIG` before constructing agents if
  the workflow should use LLM assistance.
* Prefer deterministic/no-LLM mode for tests, validation gates, and workflows
  that must be repeatable without credentials.
* Preserve output summaries, warnings, generated figure paths, and checkpoint
  locations in the final project record.


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
   * - Figures are missing after resume
     - Figures are stripped from pickled checkpoints because they may not be
       picklable.
     - Use saved figure paths from the original output directory or regenerate
       display figures.
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
