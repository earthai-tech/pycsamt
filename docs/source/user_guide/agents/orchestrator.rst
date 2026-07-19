Workflow Orchestrator
=====================

:class:`pycsamt.agents.WorkflowOrchestratorAgent` is the highest-level entry
point for the agent workflow system.  It turns either a natural-language
request or an explicit workflow configuration into a concrete
:term:`Agent coordinator` run: the workflow is selected, the required agents are
instantiated, step inputs are wired, a reproducibility plan is built, and the
chain is either previewed or executed.

The orchestrator is the right abstraction when a user says what they want in
survey language: "run phase tensor analysis", "prepare ModEM files", "do the
full AI workflow", or "make a QC report".  The coordinator is still responsible
for executing the ordered chain, but the orchestrator decides which chain should
exist in the first place.  In that sense, the orchestrator sits above
:class:`~pycsamt.agents.ContextInputAgent` and
:class:`~pycsamt.agents.AgentCoordinator`: context extraction helps interpret
the request, the orchestrator chooses the workflow, and the coordinator runs the
selected :term:`workflow steps <Workflow step>`.

.. seealso::

   :doc:`/applications/agent_master/index`
       The Agent Master application uses this orchestrator as the user-facing
       workflow dispatcher.

   :doc:`coordinator`
       The lower-level guide for building an ordered agent chain manually.

When To Use It
--------------

Use ``WorkflowOrchestratorAgent`` when the processing objective is clear but
the exact agent sequence should be inferred from the request.  It is especially
useful in notebooks, command-line tools, graphical assistants, batch dashboards,
and any interface where the user should not have to remember internal agent
class names.

For example, the request "compute phase tensor and strike analysis" resolves to
the ``phase_analysis`` workflow, which loads EDI data, checks quality, applies
static-shift correction, runs phase tensor analysis, and produces a report.  A
manual coordinator can do the same work, but only after the caller has already
chosen that five-step chain.  The orchestrator makes that choice reproducible
and inspectable before any heavy computation begins.

The central mapping can be viewed as a small routing function.  Let
:math:`r` be the request text, :math:`c` the optional configuration dictionary,
:math:`\mathcal{W}` the set of registered workflow identifiers, and
:math:`S_w = (s_1, \ldots, s_n)` the ordered step sequence for workflow
:math:`w \in \mathcal{W}`.  The orchestrator computes

.. math::

   w =
   \begin{cases}
   c_{\mathrm{workflow}}, & \text{if an explicit workflow is supplied},\\
   f_{\mathrm{LLM}}(r), & \text{if an API key is available and JSON parsing succeeds},\\
   f_{\mathrm{kw}}(r), & \text{otherwise}.
   \end{cases}

It then builds a coordinator :math:`C_w` from the registered sequence
:math:`S_w` and evaluates

.. math::

   R = C_w(x, d),

where :math:`x` is the resolved input data path and :math:`d` is the output
directory.  A :term:`dry run` evaluates the same route and builds the same
coordinator, but returns the preview instead of executing the agents.

Input Contract
--------------

The input to :meth:`~pycsamt.agents.WorkflowOrchestratorAgent.execute` is a
dictionary.  The orchestrator reads only a small set of top-level keys and
forwards the rest through the workflow configuration where relevant.

.. list-table::
   :header-rows: 1
   :widths: 24 18 58

   * - Key
     - Required
     - Meaning
   * - ``request``
     - Usually
     - Natural-language request used for workflow classification, reasoning,
       path extraction, and provenance.
   * - ``config``
     - No
     - Structured workflow configuration.  ``config["workflow"]`` has priority
       over request classification.
   * - ``data_path``
     - Usually
     - Survey directory or EDI path passed to the root load step.  If omitted,
       the orchestrator tries to extract a path from ``request``.
   * - ``output_dir``
     - No
     - Root directory for reports, figures, solver files, provenance files, and
       checkpoints.  Defaults to ``"pycsamt_workflow_output"``.
   * - ``dry_run``
     - No
     - When ``True``, preview the chain without running the processing agents.
   * - ``config["step_params"]``
     - No
     - Per-step parameter overrides.  Values under ``"load"`` are merged into
       the root execution config; values for later steps are injected through
       the step's :term:`input mapping`.
   * - ``config["checkpoint"]``
     - No
     - Checkpoint path forwarded to AI, PINN, hybrid, Inv2D, and Inv3D
       inversion steps when those steps are present.

Explicit configuration is the strongest reproducibility mode because it removes
route ambiguity:

.. code-block:: pycon

   >>> from contextlib import redirect_stdout
   >>> from io import StringIO
   >>> from pycsamt.agents import WorkflowOrchestratorAgent
   >>> from pycsamt.api.agents import AGENT_CONFIG
   >>> with AGENT_CONFIG.offline():
   ...     buffer = StringIO()
   ...     with redirect_stdout(buffer):
   ...         result = WorkflowOrchestratorAgent().execute({
   ...             "config": {"workflow": "modem"},
   ...             "data_path": "/data/WILLY_EDIs",
   ...             "output_dir": "/out/willy_modem",
   ...             "dry_run": True,
   ...         })
   >>> print(result.status)
   success
   >>> print(result["workflow_type"])
   modem
   >>> print([step["name"] for step in result["steps"]])
   ['load', 'qc', 'static_shift', 'modem', 'report']

Here, the request text is not needed because ``config["workflow"]`` already
selects the registered ``modem`` chain.  The dry run still creates the
coordinator and step metadata, so the caller can display the exact plan before
running a file-writing job.

Routing And Reproducibility
---------------------------

Routing follows a fixed priority order.  First, ``config["workflow"]`` is used
when present.  Second, if the orchestrator was created with an API key, the LLM
is asked to return a JSON object containing ``workflow_type`` and ``reasoning``.
Third, when no valid LLM route is available, pycsamt uses the shared keyword
classifier from :mod:`pycsamt.agents._workflows`.  This shared registry is also
used by :class:`~pycsamt.agents.ContextInputAgent`, which keeps text parsing and
workflow execution aligned.

No-LLM routing is deterministic.  The classifier scans ordered keyword groups,
so specific plotting, inversion, or solver phrases are matched before broad
phrases.  Compound requests that mention several processing families can route
to ``full`` instead of the first keyword match.  If no keyword is recognized,
the classifier returns ``qc`` with a reasoning message explaining the fallback.

.. code-block:: pycon

   >>> from contextlib import redirect_stdout
   >>> from io import StringIO
   >>> from pycsamt.agents import WorkflowOrchestratorAgent
   >>> from pycsamt.api.agents import AGENT_CONFIG
   >>> examples = [
   ...     "run QC on the data",
   ...     "compute phase tensor and strike analysis",
   ...     "set up Occam2D mesh and startup file",
   ...     "prepare ModEM 3D inversion",
   ...     "ensemble uncertainty quantification",
   ... ]
   >>> with AGENT_CONFIG.offline():
   ...     agent = WorkflowOrchestratorAgent()
   ...     for request in examples:
   ...         buffer = StringIO()
   ...         with redirect_stdout(buffer):
   ...             result = agent.execute({
   ...                 "request": request,
   ...                 "data_path": "/data/test",
   ...                 "dry_run": True,
   ...             })
   ...         print(f"{request} -> {result['workflow_type']}")
   run QC on the data -> qc
   compute phase tensor and strike analysis -> phase_analysis
   set up Occam2D mesh and startup file -> pre_inversion
   prepare ModEM 3D inversion -> modem
   ensemble uncertainty quantification -> ensemble_inversion

The ``default_workflow`` constructor argument is only a last-resort guard when
an unregistered route reaches the validation stage.  It should not be treated as
a general override for vague natural-language text, because no-keyword text
already resolves to ``qc`` in the deterministic classifier.  When the workflow
must not be inferred, pass ``config={"workflow": ...}``.

Output Contract
---------------

The orchestrator returns an :term:`AgentResult`.  Its status mirrors the nested
coordinator result, and its data dictionary exposes both the route and the
objects used to execute it.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Key
     - Meaning
   * - ``workflow_type``
     - Selected workflow identifier, such as ``"qc"``, ``"phase_analysis"``,
       ``"modem"``, or ``"full_ai_workflow"``.
   * - ``reasoning``
     - Short explanation from LLM routing or keyword routing.
   * - ``workflow_plan``
     - Structured plan used for validation and provenance.
   * - ``coordinator``
     - Built :class:`~pycsamt.agents.AgentCoordinator` instance.
   * - ``result``
     - Nested coordinator :term:`AgentResult`.
   * - ``steps``
     - List of dictionaries containing ``name``, ``agent``, and
       ``description`` for every registered step that will run.

For a phase-analysis dry run, the outer result tells you which route was
selected while the nested coordinator result tells you what the coordinator
would execute:

.. code-block:: pycon

   >>> from contextlib import redirect_stdout
   >>> from io import StringIO
   >>> from pycsamt.agents import WorkflowOrchestratorAgent
   >>> from pycsamt.api.agents import AGENT_CONFIG
   >>> with AGENT_CONFIG.offline():
   ...     buffer = StringIO()
   ...     with redirect_stdout(buffer):
   ...         result = WorkflowOrchestratorAgent().execute({
   ...             "request": "compute phase tensor and strike analysis",
   ...             "data_path": "/data/WILLY_EDIs",
   ...             "output_dir": "/out/willy_phase",
   ...             "dry_run": True,
   ...         })
   >>> print(result.status)
   success
   >>> print(result["workflow_type"])
   phase_analysis
   >>> print([step["name"] for step in result["steps"]])
   ['load', 'qc', 'static_shift', 'phase_analysis', 'report']
   >>> print(result["reasoning"])
   Matched keywords for workflow 'phase_analysis'.
   >>> print(result["result"].summary)
   Workflow preview: 5 steps.

The outer summary combines the route and the coordinator summary.  The outer
warnings combine route-level warnings, workflow-plan risk flags, and
coordinator warnings, so production callers should inspect ``result.warnings``
even when ``result.status == "success"``.

Supported Workflow Families
---------------------------

The registry maps each workflow identifier to an ordered chain.  The table
below groups the most important routes by the work they represent.

.. list-table::
   :header-rows: 1
   :widths: 24 42 34

   * - Workflow
     - Main chain
     - Typical request
   * - ``qc``
     - load -> QC -> static shift -> report
     - Clean the survey, flag poor bands, summarize quality.
   * - ``static_shift``
     - load -> QC -> static shift -> report
     - Detect and correct galvanic/static-shift effects.
   * - ``phase_analysis``
     - load -> QC -> static shift -> phase analysis -> report
     - Compute phase tensor, strike, skew, and dimensionality products.
   * - ``pre_inversion``
     - load -> QC -> static shift -> phase analysis -> Occam2D -> code
     - Prepare reproducible 2-D inversion inputs and a runnable script.
   * - ``modem``
     - load -> QC -> static shift -> ModEM files -> report
     - Prepare 3-D ModEM data, model, covariance, and control files.
   * - ``mare2dem``
     - load -> QC -> static shift -> MARE2DEM files -> report
     - Prepare 2.5-D MARE2DEM input files.
   * - ``ai_inversion`` / ``inv1d``
     - load -> QC -> denoise -> AI 1-D inversion -> interpretation -> report
     - Estimate 1-D resistivity models with EMInverter1D/CNN workflows.
   * - ``inv2d``
     - load -> QC -> denoise -> U-Net 2-D inversion -> interpretation -> report
     - Build a profile-scale 2-D AI section.
   * - ``inv3d``
     - load -> QC -> static shift -> GCN 3-D inversion -> interpretation -> report
     - Infer spatial 3-D structure from station-neighborhood information.
   * - ``ensemble_inversion``
     - load -> QC -> denoise -> ensemble inversion -> interpretation -> report
     - Quantify uncertainty with multiple 1-D model predictions.
   * - ``joint_inversion``
     - load -> QC -> static shift -> joint inversion -> interpretation -> report
     - Combine MT with complementary modalities or survey constraints.
   * - ``pinn_inversion``
     - load -> QC -> PINN inversion -> interpretation -> report
     - Use physics-informed optimization for resistivity estimation.
   * - ``hybrid_inversion``
     - load -> QC -> hybrid inversion -> interpretation -> report
     - Warm-start with AI, then refine with physics-guided inversion.
   * - ``tipper``
     - load -> tipper analysis -> report
     - Analyze induction arrows and tipper amplitudes.
   * - ``sensitivity``
     - load -> QC -> sensitivity/DOI -> report
     - Estimate depth of investigation and response sensitivity.
   * - ``rotation``
     - load -> QC -> phase analysis -> tensor rotation
     - Rotate tensors into a strike-consistent coordinate frame.
   * - ``freq_decimation``
     - load -> QC -> frequency decimation -> report
     - Select stable periods or reduce frequency density.
   * - ``full``
     - load -> QC -> static shift -> phase analysis -> denoise -> AI inversion
       -> Occam2D -> report
     - Run the classical end-to-end chain from QC to report.
   * - ``full_ai_workflow``
     - load -> QC -> static shift -> phase analysis -> denoise -> 1-D AI ->
       2-D AI -> code -> report
     - Run the extended AI-assisted chain.

The registry also contains focused routes such as ``forward``,
``inversion_eval``, ``interpretation``, ``report``, ``code_gen``, ``denoise``,
``batch``, and ``comparison``.  Plotting and correction tools are routed through
the shared workflow keyword table, and the Agent Master can use those routes to
dispatch more specialized single-purpose actions.

Running A Real Workflow
-----------------------

Remove ``dry_run`` only after the route and step list are acceptable.  A real
run can read EDI files, write reports and figures, create solver input files,
call configured LLM providers, and update :term:`workflow checkpoints
<Workflow checkpoint>` under ``<output_dir>/.checkpoints``.

.. code-block:: pycon

   >>> from pycsamt.agents import WorkflowOrchestratorAgent
   >>> result = WorkflowOrchestratorAgent().execute({
   ...     "request": "run phase tensor analysis on the WILLY survey",
   ...     "data_path": "/data/WILLY_EDIs",
   ...     "output_dir": "/out/willy_phase",
   ... })
   >>> print(result.summary)
   >>> print(result["workflow_type"])
   >>> print(result["result"].summary)

A non-preview run also writes provenance files into the output directory:
``workflow_plan.json`` records the validated plan, ``agent_trace.json`` records
the route and executed agents, ``environment.json`` records Python and package
versions, and ``output_manifest.json`` records generated files with hashes.  The
manifest is important for reproducibility because it ties the selected workflow
to the files produced by that run.

LLM-Assisted Routing
--------------------

When an API key is supplied, the orchestrator asks the configured provider to
classify the request before falling back to keyword routing.  The LLM response
is expected to be JSON with two fields: ``workflow_type`` and ``reasoning``.
The selected workflow still has to exist in the registry; an invented workflow
name is rejected rather than executed.

.. code-block:: pycon

   >>> from pycsamt.agents import WorkflowOrchestratorAgent
   >>> agent = WorkflowOrchestratorAgent(
   ...     api_key="sk-ant-...",
   ...     llm_provider="claude",
   ...     model="claude-sonnet-4-6",
   ... )
   >>> result = agent.execute({
   ...     "request": (
   ...         "The survey has noisy bands and I want an uncertainty-aware "
   ...         "AI inversion report."
   ...     ),
   ...     "data_path": "/data/WILLY_EDIs",
   ...     "output_dir": "/out/willy_uncertainty",
   ...     "dry_run": True,
   ... })
   >>> print(result["workflow_type"])
   >>> print(result["reasoning"])

This mode is useful for natural phrasing, but the reproducibility boundary is
still the same: preview first, inspect ``workflow_type`` and ``steps``, then
execute the same resolved request or an explicit ``config``.

Inspecting Planned Steps
------------------------

The ``steps`` list is the simplest representation to show in a CLI, notebook,
or application confirmation dialog.  Each row is the orchestrator's resolved
view of the workflow registry.

.. code-block:: pycon

   >>> from contextlib import redirect_stdout
   >>> from io import StringIO
   >>> from pycsamt.agents import WorkflowOrchestratorAgent
   >>> from pycsamt.api.agents import AGENT_CONFIG
   >>> with AGENT_CONFIG.offline():
   ...     buffer = StringIO()
   ...     with redirect_stdout(buffer):
   ...         result = WorkflowOrchestratorAgent().execute({
   ...             "request": "compute phase tensor and strike analysis",
   ...             "data_path": "/data/WILLY_EDIs",
   ...             "dry_run": True,
   ...         })
   >>> for index, step in enumerate(result["steps"], start=1):
   ...     print(f"{index:02d}. {step['name']} ({step['agent']})")
   01. load (MTLoaderAgent)
   02. qc (DataQCAgent)
   03. static_shift (StaticShiftAgent)
   04. phase_analysis (PhaseAnalysisAgent)
   05. report (ReportAgent)

For production interfaces, treat this preview as the confirmation boundary.
The user can see the route, the ordered agents, and the requested output
directory before the run touches survey data or creates files.

Validation, Warnings, And Failure Modes
---------------------------------------

Before a non-preview run executes, the orchestrator builds a workflow plan and
validates it.  Hard validation errors block execution and return a failed
:term:`AgentResult`; dry runs are allowed to return the plan even when the input
path does not exist, because their purpose is to expose the intended route.

There are three common failure classes:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Failure
     - Meaning
   * - Unknown workflow
     - The selected workflow identifier is not registered.  Use a known
       workflow name or pass an explicit configuration.
   * - Missing root agent
     - The first required agent in the chain could not be imported or
       instantiated.  The run stops immediately because downstream steps have
       no valid input.
   * - Plan validation error
     - A non-preview run failed the workflow-plan checks, for example because
       required execution information is missing.

Warnings should be read as reproducibility signals.  They can indicate missing
paths during preview, risk flags from the workflow plan, skipped optional
agents, or downstream coordinator warnings.  A successful status means the
workflow result is usable, not that every scientific assumption has been
validated.

Recommended Operating Pattern
-----------------------------

For user-facing tools, use a two-pass pattern.  First, run a dry preview and
show ``workflow_type``, ``reasoning``, ``steps``, ``warnings``, and
``output_dir``.  Second, after the user accepts the plan, run the same request
with the same path and output directory without ``dry_run``.  For scripted
production jobs, store the selected workflow in ``config["workflow"]`` once it
has been approved; that turns an interpreted request into an explicit,
repeatable workflow.

Related API
-----------

See :class:`pycsamt.agents.WorkflowOrchestratorAgent`,
:class:`pycsamt.agents.AgentCoordinator`,
:class:`pycsamt.agents.AgentResult`, :mod:`pycsamt.agents.orchestrator`, and
:mod:`pycsamt.agents._workflows`.
