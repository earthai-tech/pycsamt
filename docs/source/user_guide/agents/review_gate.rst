.. _agents-review-gate:

Review Gate
===========

A natural-language request reaches the agent loop as free text, and a
language model turns it into a workflow choice and a set of parameters.
Two things can go wrong on that path. The plan can be malformed (an
unknown workflow, a missing data path, a period band in the wrong
order), and the plan can be well formed but answer a different question
than the one the user asked. The first failure is caught by plan
validation. The second is harder: a misread request can name a valid,
registered workflow and pass every check on its parameters. pyCSAMT
therefore puts a :term:`review gate` between planning and execution in
:class:`~pycsamt.agents.WorkflowOrchestratorAgent`, made of three parts:

* **plan validation**, which blocks a plan with hard errors;
* a deterministic **route cross-check**, which compares the language
  model's route with a keyword route and turns a disagreement into a
  risk flag;
* **plan-hash approval**, which lets a plan that needs review run only
  when the caller returns the hash of that exact plan.

Writing :math:`P` for the validated plan, :math:`I` for the identity
(names and bytes) of its input files, :math:`E(P)` for the list of
validation errors and :math:`a` for the approval hash supplied by the
caller, a non-preview run executes only when

.. math::

   E(P) = \varnothing
   \quad\text{and}\quad
   \bigl(\neg\,\mathrm{review}(P) \;\lor\; a = H(P, I)\bigr),

where :math:`\mathrm{review}(P)` is the review flag of the plan and
:math:`H` is a SHA-256 digest of the plan and its inputs. A
:term:`dry run` is exempt from both conditions because it executes
nothing.

.. seealso::

   :doc:`orchestrator`
       The workflow dispatcher that hosts the gate.

   :doc:`plugins`
       The plugin interface through which another domain can supply its
       own validator and cross-check.


Plan Validation
---------------

Every request is converted into a
:class:`~pycsamt.agents._workflow_plan.WorkflowPlan` before any agent
runs. ``WorkflowPlan.validation_errors`` checks that the request text is
not empty, that the workflow identifier is one of the registered
workflows, that ``data_path`` exists, and that the core scientific
settings are coherent: ``period_range`` must hold two finite positive
bounds in increasing order, ``component`` and ``inversion_code`` must be
supported values, depths must be finite and positive, and counts such
as ``n_layers`` or ``epochs`` must be positive integers. The same checks
are applied to nested per-step settings under ``step_params``.

Any error blocks a non-preview run. The orchestrator returns a failed
:term:`AgentResult` whose ``error`` lists the problems and whose
``data["workflow_plan"]`` holds the rejected plan, so the interface can
show the user what was parsed. Nothing is written to the output
directory.

Validation is a bounded contract. It does not try to check every
solver-specific option, and it never repairs scientific intent: a
reversed period band is reported, not silently swapped.


Route Cross-Check
-----------------

When a language model chooses the workflow, pyCSAMT also classifies the
same text with the deterministic keyword table used by the offline
router. The two answers are compared at the level of the
:term:`analysis family`, because several workflow identifiers are
sub-tools of one analysis: ``strike``, ``pt_psection`` and
``phase_analysis`` all belong to the ``phase_analysis`` family, and the
static-shift correction methods ``corr_ss_ama`` or ``corr_ss_emap``
belong to ``static_shift``. An identifier outside every family is its
own family. The families are listed in
``pycsamt.agents._workflows.WORKFLOW_FAMILIES``.

The comparison runs at two levels:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Level
     - Behaviour on disagreement
   * - Router (``router.route_disagreement``)
     - :class:`~pycsamt.agents.router.IntentRouter` compares the language
       model's intent and workflow with the offline reading. A
       disagreement turns the request into a clarification question
       that names both readings, so the user chooses before anything is
       planned.
   * - Orchestrator (``orchestrator._route_cross_check``)
     - When the orchestrator itself asked the language model for the
       workflow, a family mismatch is added to the plan's
       ``risk_flags``. Any risk flag sets ``requires_human_review``, so
       the plan stops at plan-hash approval.

Two cases never constrain the language model. If the keyword table
matches nothing, there is no deterministic reading to compare with, and
the route stands. The keyword identifier ``code_gen`` describes an
output mode (write a script) and not an analysis, so it never counts as
a disagreement either. A workflow passed explicitly in ``config`` is not
cross-checked: the caller has already made the choice.

The built-in validator exposes both checks through the plugin
interface, which is the supported way to call them:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents.plugins import get_plugin
   >>> plugin = get_plugin()
   >>> plugin.registry.family("pt_psection"), plugin.registry.family("corr_ss_ama")
   ('phase_analysis', 'static_shift')
   >>> validator = plugin.validator
   >>> text = "estimate geoelectric strike for L26PLT"
   >>> validator.cross_check(text, "workflow", "strike") is None
   True
   >>> validator.route_flags(text, "phase_analysis")
   []
   >>> text = "denoise the survey data for L26PLT"
   >>> validator.cross_check(text, "workflow", "qc")
   "workflow 'denoise'"
   >>> validator.route_flags(text, "qc")
   ["LLM chose workflow 'qc' but the keyword route reads 'denoise'; confirm the workflow before it runs."]
   >>> validator.route_flags("tell me something about my survey", "inv2d")
   []

The first request agrees: the language model chose ``strike`` and the
keywords read ``phase_analysis``, which is the same family. The second
disagrees, because ``qc`` and ``denoise`` are different analyses; the
flag text is what the user sees next to the plan. The third request
matches no keyword, so the cross-check returns nothing and the model's
``inv2d`` choice is left to plan validation and the remaining review
conditions.


Plan-Hash Approval
------------------

A plan needs review when any of the following holds:

* the workflow executes an external inversion (``conventional_inversion``
  or ``field_inversion``);
* the configuration sets ``requires_human_review=True``;
* the plan carries at least one risk flag, including a route
  cross-check flag, a missing data path or an unrecognised workflow.

For such a plan the orchestrator computes
``WorkflowPlan.approval_digest()``: a SHA-256 hash of the serialised plan
(request text, workflow, data path, output directory, every parsed
parameter including ``step_params``, risk flags, provider and expected
outputs) together with the names and content hashes of the input files.
Without a matching ``approved_plan_hash`` in the request, the run stops
with status ``needs_review`` and returns the hash under
``data["approval_plan_hash"]``. A run record event ``review_gate`` with
``decision="blocked"`` or ``decision="approved"`` is written either way.

The caller approves by sending the same request again with that hash.
Because the hash covers the parameters and the input bytes, approval
cannot outlive a change: editing a period band, switching the language
model provider, or replacing an EDI file yields a new hash, and the old
approval no longer matches. A truthy flag such as ``approved=True`` is
not accepted.

The example below uses the public Baohuashan line bundled with the
repository. ``requires_human_review=True`` forces the gate on a QC run so
that the whole cycle is short; the coordinator's progress messages are
captured to keep the output readable.

.. code-block:: pycon
   :linenos:

   >>> import contextlib, io
   >>> from pycsamt.agents import WorkflowOrchestratorAgent

   >>> orch = WorkflowOrchestratorAgent(api_key=None)
   >>> request = {
   ...     "request": "run QC on line L22PLT",
   ...     "data_path": "data/AMT/WILLY_DATA/L22PLT",
   ...     "output_dir": "runs/review_gate_demo",
   ...     "config": {"workflow": "qc", "requires_human_review": True},
   ... }
   >>> first = orch.execute(request)
   >>> first.status
   'needs_review'
   >>> plan_hash = first.data["approval_plan_hash"]
   >>> len(plan_hash)
   64

   >>> changed = dict(
   ...     request,
   ...     config=dict(request["config"], period_range=[1e-3, 10.0]),
   ...     approved_plan_hash=plan_hash,
   ... )
   >>> orch.execute(changed).status
   'needs_review'

   >>> with contextlib.redirect_stdout(io.StringIO()):
   ...     approved = orch.execute(dict(request, approved_plan_hash=plan_hash))
   >>> approved.status
   'success'

The first call returns the plan without running any agent. The second
call reuses the hash after adding a period band; the plan is no longer
the one that was approved, so the gate stops it again and returns a new
hash. The third call sends the unchanged request with the original hash
and runs the QC chain. The hash value itself is not shown because it
depends on the absolute path and bytes of the input files.

An interface built on the orchestrator should follow the same cycle:
show the returned plan, its risk flags and its expected outputs, and
send back the hash of the plan that was displayed, never a hash
recomputed after the user edited something.


In The Agent Master Chat
------------------------

The Agent Master web chat (``pycsamt.app.agent_master``) follows that
cycle. Since 2.6.6:

* the chat picks the workflow from the parameter form, the keyword
  registry or the language model, and the orchestrator cross-check then
  compares that choice with the keyword route; a disagreement adds a risk
  flag, so the plan goes to the gate;
* a held plan is shown in the chat with its workflow, data, steps and
  risk flags and an **Approve and run** button. The button resends the
  plan that was displayed, with its ``approved_plan_hash``. The approval
  is kept in the server's memory only; it can be used once, and a restart
  of the app invalidates it;
* catalogue corrections (static shift, noise removal, rotation,
  coordinate and Stratagem methods) run as orchestrator workflows, one per
  correction id, with a ``load`` step and a ``correction`` step. They are
  validated, cross-checked, gated when flagged and recorded. The corrected
  data stay in memory for the post-processing dialog; the input files are
  not changed.

Plots and the other tools (strike, dimensionality, station validation,
coordinates, elevation, format conversion, batch plot export, the
frequency editor and the layered-model builder) still run directly and
are listed in the workflow history. None of them modifies the input
files: the converter and the batch export write new files to an export
folder, and the frequency editor, like a correction, returns an edited
copy that reaches the session only through the post-processing dialog.


Reading The Result
------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - ``status``
     - Meaning
   * - ``failed``
     - Plan validation found errors (listed in ``error``), or a
       workflow-specific precondition is missing; for
       ``field_inversion`` the orchestrator names every missing
       scientific setting and runs no step.
   * - ``needs_review``
     - The plan is valid but needs approval. Show
       ``data["workflow_plan"]`` and its ``risk_flags``, then resend the
       request with ``approved_plan_hash``.
   * - ``success``
     - The approved plan ran. ``success`` means the programmed chain
       completed; the scientific checks of each step still apply.

The gate decides whether a plan may run. It does not judge whether the
science is right: a QC report, an inversion misfit and the Occam2D
mesh-depth warning (see :ref:`models_occam2d`) remain the evidence for
that.
