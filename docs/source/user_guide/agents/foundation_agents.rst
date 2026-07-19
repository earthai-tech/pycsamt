Foundation And Survey Intake Agents
===================================

These components define the agent execution contract and the first steps of
most pyCSAMT workflows: parse intent, load data, and coordinate downstream
processing.  Every class documented here either *is* the contract --
``BaseAgent`` and ``AgentResult`` -- or is the first :term:`agent` a
workflow actually calls, so the order below doubles as the order most
workflows execute in: understand the contract, then intake.

.. _agent-base-agent:

BaseAgent
---------

``BaseAgent`` is the shared base class for all agent classes.  It resolves the
effective LLM configuration, exposes ``query_llm()``, tracks per-call cost,
provides JSON extraction helpers, and injects pyCSAMT plotting/style helpers
so figures produced by agents remain consistent with the rest of the package.
Every subclass in this catalogue inherits the same construction contract --
provider, key, and model resolved once at ``__init__`` -- so switching an
agent's LLM backend never means touching its ``execute`` logic; see
:doc:`llm_configuration` for how that resolution chooses between an explicit
argument, a stored key, an environment variable, and the global
``AGENT_CONFIG`` defaults.

``query_llm()`` is the one method most subclasses actually call.  It returns
``None`` immediately when no API key is configured -- the graceful
degradation every example on this page and in :doc:`agent_catalogue` relies
on -- and otherwise retries a failed call only when the error looks like a
rate limit, waiting longer each time:

.. math::

   d_i = 2^{\,i-1}\ \text{s}, \qquad i = 1, 2, 3,

so three attempts are spaced 1 s, 2 s, and 4 s apart before ``query_llm``
gives up and returns ``None``.  A non-rate-limit failure (a bad key, a
malformed request) is not retried at all, since waiting cannot fix it.  Every
successful call also prices itself against the resolved provider/model rate
table and folds the result into ``self._last_cost``:

.. math::

   \text{cost}_{\text{USD}} =
   \frac{n_{\text{in}} \cdot r_{\text{in}} + n_{\text{out}} \cdot r_{\text{out}}}
        {10^{6}},

with :math:`r_{\text{in}}`, :math:`r_{\text{out}}` the provider's USD-per-
million-token input/output rates and :math:`n_{\text{in}}`,
:math:`n_{\text{out}}` the token counts the provider itself reports.  For
``claude-sonnet-4-6`` (:math:`r_{\text{in}}=3.0`, :math:`r_{\text{out}}=15.0`
at the time of writing) a 500-input/200-output-token call costs
:math:`(500 \cdot 3.0 + 200 \cdot 15.0)/10^{6} = \$0.0045` -- small per call,
which is why :doc:`llm_configuration` also covers session-wide budgets
rather than relying on per-call intuition alone.

Beyond the LLM call itself, ``BaseAgent`` gives every subclass a few small,
reusable helpers worth knowing before writing a custom agent:
``extract_json()`` pulls the first balanced JSON object or array out of a
free-text LLM response (useful when a prompt asks for JSON but the model
wraps it in prose); ``require_keys()`` checks an input dictionary against a
list of required keys and logs the ones missing; and ``_save_figure()``
writes a Matplotlib figure through the shared ``PLOT_CONFIG`` so every
agent's figures share one save path and style.  Agents that ground an
answer or a plan in real package evidence -- rather than just running a
deterministic computation -- can also call the private ``_retrieve_context``
hook, which degrades to ``None`` exactly as silently as ``query_llm`` does
when :doc:`assistant_rag` is unavailable or returns nothing.

Use it when building custom agents:

.. code-block:: pycon

   >>> from pycsamt.agents import BaseAgent, AgentResult

   >>> class MyAgent(BaseAgent):
   ...     def __init__(self):
   ...         super().__init__("MyAgent")
   ...
   ...     def execute(self, input_data):
   ...         self._last_cost = 0.0
   ...         return AgentResult(
   ...             status="success",
   ...             summary="Custom step complete.",
   ...             data={"value": 1},
   ...             cost_estimate_usd=self._last_cost,
   ...         )

   >>> agent = MyAgent()
   >>> print(repr(agent))
   MyAgent(name='MyAgent', llm='claude/claude-sonnet-4-6')
   >>> print(agent.execute({}))
   AgentResult(status='success', elapsed=0.0s)

``llm='claude/claude-sonnet-4-6'`` in that ``repr`` is the *resolved
default*, not proof of a working connection -- with no API key configured,
``MyAgent`` above still ran end to end, because nothing in this minimal
``execute`` actually calls ``query_llm()``.  Two more helpers used the same
way, without needing a full agent:

.. code-block:: pycon

   >>> agent.extract_json('Sure, here is the plan: {"steps": 3, "ok": true} -- done.')
   {'steps': 3, 'ok': True}
   >>> agent.require_keys({"path": "x"}, "path", "sites")
   ['sites']

.. _agent-result:

AgentResult
-----------

``AgentResult`` is the standard return object for every :term:`agent`.  It
contains status, summary, structured data, warnings, optional LLM
interpretation, run time, cost estimate, and failure details -- the fields
``MyAgent`` above filled in by hand, and that every built-in agent fills in
the same way.

The result supports dict-like access to ``data``, and its truthiness tracks
``status``: only ``"failed"`` is falsy, so ``if result:`` is the idiomatic
success check rather than comparing ``status`` to a string.  Both branches,
run for real against :ref:`MTLoaderAgent <agent-mt-loader>` (introduced
next):

.. code-block:: pycon

   >>> from pycsamt.agents import MTLoaderAgent

   >>> result = MTLoaderAgent().execute({"path": "data/AMT/WILLY_DATA/L18PLT"})
   >>> if result:
   ...     print(result.summary)
   ...     print(result["n_stations"])
   ...     print(result.get("summary_stats"))
   ... else:
   ...     print(result.error)
   ...     print(result.error_fix_hint)
   Loaded 28 station(s) from L18PLT. Mean QC score 83/100. Period range: 9.62e-05–9.92e-01 s.
   28
   {'n_stations': 28, 'n_with_z': 28, 'n_with_coords': 28, 'n_with_tipper': 0, 'mean_qc_score': 82.96428571428571, 'min_qc_score': 81.0, 'median_n_freq': 53.0, 'global_t_min_s': 9.615384615384615e-05, 'global_t_max_s': 0.9920634920634921}

A path with nothing loadable in it takes the other branch:

.. code-block:: pycon

   >>> missing = MTLoaderAgent().execute({"path": "data/AMT/WILLY_DATA/DOES_NOT_EXIST"})
   >>> print(missing.status)
   failed
   >>> print(missing.error)
   No stations with valid impedance data found.
   >>> print(missing.error_fix_hint)
   Ensure the files contain ZXXR/ZXXI ... blocks or RES/PHS data.

Notice the error does not say "path not found" -- ``MTLoaderAgent`` does not
distinguish an absent directory from an existing one with no usable EDIs; both
collapse to the same downstream "no impedance data" failure, so a typo'd
path and a genuinely empty survey folder read identically from the
``AgentResult`` alone.  Checking ``n_stations`` in the input path with plain
``os``/``pathlib`` tools first, before handing it to the agent, is the more
diagnosable habit when that distinction matters.

.. _agent-context-input:

ContextInputAgent
-----------------

``ContextInputAgent`` converts natural-language workflow requests into a
structured configuration dictionary.  It can call an LLM when configured, but
also includes a regex fallback for common request patterns such as loading a
path, selecting a workflow, and extracting period ranges.  The regex path
patterns only match an absolute-looking path -- one starting with ``/`` or
``~`` -- so a request that only names a bare relative path (``data/edis``,
with no leading slash) will not populate ``data_path`` without an LLM key;
pass the path through a dedicated ``data_path`` config key in that case
instead of relying on the sentence alone.

Typical input keys:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``request``
     - Natural-language request.  Preferred key.
   * - ``text`` or ``prompt``
     - Alternative text keys accepted by the agent.

Typical output keys include ``config``, extracted path fields, workflow name,
period range, and validation warnings.

.. code-block:: pycon

   >>> from pycsamt.agents import ContextInputAgent

   >>> agent = ContextInputAgent()
   >>> result = agent.execute({
   ...     "request": "Load EDIs from /data/WILLY, QC them, period 0.001 to 10 s",
   ... })

   >>> print(result["config"])
   {'workflow': 'qc', 'data_path': '/data/WILLY', 'period_range': [0.001, 10.0], 'component': 'xy', 'station': None, 'inversion_code': None, 'output_dir': 'pycsamt_agent_output', 'verbose': True}
   >>> print(result.llm_interpretation)
   None

The leading ``/`` on ``/data/WILLY`` above is what let the regex fallback
extract ``data_path`` at all -- ``llm_interpretation`` stays ``None``
because no LLM key is configured, and the request was still parsed into a
complete, workflow-ready config regardless.  ``"period 0.001 to 10 s"`` in
plain English became ``period_range: [0.001, 10.0]`` and ``"QC"`` became
``workflow: 'qc'``: both are matched by the same kind of fixed
pattern-to-field table as the :term:`query expansion` used in
:doc:`assistant_rag`, not by anything resembling language understanding, so
requests that stray far from these phrasings are exactly where an LLM key
starts to earn its cost.

Use this agent at the start of natural-language workflows or before building
an ``AgentCoordinator`` from user text.

.. _agent-mt-loader:

MTLoaderAgent
-------------

``MTLoaderAgent`` loads MT, AMT, or CSAMT survey data into a pyCSAMT ``Sites``
object.  It accepts a single file, a directory, a list of paths, an existing
``Sites`` object, or an EDI collection supported by the lower-level loading
utilities.  ``sites`` is what every processing, inversion, and reporting
agent downstream actually consumes -- everything else in ``data`` is a
summary of that same object, computed once so later steps and reports do
not each recompute it.

Typical output keys:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``sites``
     - Loaded ``Sites`` object for downstream agents.
   * - ``station_names``
     - Ordered station names.
   * - ``n_stations``
     - Number of loaded stations.
   * - ``quality_table``
     - Per-station loading and quality summary.
   * - ``summary_stats``
     - Survey-level loading and quality statistics.

The success run from :ref:`AgentResult <agent-result>` above already showed
``n_stations`` and ``summary_stats``; the remaining two output keys round
out the picture:

.. code-block:: pycon

   >>> result = MTLoaderAgent().execute({"path": "data/AMT/WILLY_DATA/L18PLT"})
   >>> result["station_names"][:5]
   ['18-001A', '18-002U', '18-003A', '18-004A', '18-005U']
   >>> result["quality_table"].columns.tolist()
   ['station', 'has_z', 'has_tipper', 'has_coords', 'n_freq', 't_min_s', 't_max_s', 'snr_proxy', 'qc_score']
   >>> len(result["sites"])
   28

``quality_table`` is the per-station row behind the survey-level averages in
``summary_stats`` -- ``qc_score`` here is a loading-stage completeness and
SNR check, a coarser, faster signal than the dedicated confidence ratio
:ref:`DataQCAgent <agent-data-qc>` computes once a ``Sites`` object exists;
treat a low score here as a reason to inspect a station before it enters
QC, not as a substitute for QC itself.

Use this agent whenever a workflow needs a validated ``Sites`` object before
QC, static-shift correction, phase analysis, or inversion preparation.

.. _agent-coordinator-class:

AgentCoordinator
----------------

``AgentCoordinator`` is not a ``BaseAgent`` subclass, but it is the central
workflow runner for explicit multi-agent pipelines.  It registers ordered
steps, maps previous results into the next step with ``input_fn``, supports
:term:`dry run` previews, writes per-step checkpoints, and aggregates cost.
:doc:`coordinator` is the full guide -- input mapping, resume, required
versus optional steps, and a complete QC-and-correction chain; this section
only closes the loop between the two agents just introduced.

.. code-block:: pycon

   >>> from pycsamt.agents import AgentCoordinator, ContextInputAgent, MTLoaderAgent

   >>> coord = AgentCoordinator("load_preview")
   >>> coord.add_step("parse", ContextInputAgent())
   >>> coord.add_step(
   ...     "load",
   ...     MTLoaderAgent(),
   ...     input_fn=lambda r: {
   ...         "path": (r["parse"].get("config") or {}).get("data_path", "")
   ...     },
   ... )

   >>> result = coord.execute(
   ...     {"request": "Load /data/WILLY_EDIs"},
   ...     dry_run=True,
   ... )
   Workflow: load_preview
   Steps   : 2
   Config  : {
     "request": "Load /data/WILLY_EDIs"
   }

   ────────────────────────────────────────────────────────────
      1. [parse]
          Agent  : ContextInputAgent
          LLM    : claude/claude-sonnet-4-6
          Action : Run ContextInputAgent
      2. [load]
          Agent  : MTLoaderAgent
          LLM    : claude/claude-sonnet-4-6
          Action : Run MTLoaderAgent
   ────────────────────────────────────────────────────────────

   >>> print(result["plan"])
   Workflow: load_preview
   Steps   : 2
   Config  : {
     "request": "Load /data/WILLY_EDIs"
   }

   ────────────────────────────────────────────────────────────
      1. [parse]
          Agent  : ContextInputAgent
          LLM    : claude/claude-sonnet-4-6
          Action : Run ContextInputAgent
      2. [load]
          Agent  : MTLoaderAgent
          LLM    : claude/claude-sonnet-4-6
          Action : Run MTLoaderAgent
   ────────────────────────────────────────────────────────────

That banner appears twice on purpose here, not as a copy-paste slip: the
coordinator prints it once, automatically, the moment ``execute(dry_run=True)``
builds the preview, and ``result["plan"]`` stores that exact same text for
``print()``, logging, or a report to reuse later without re-running
anything.  The two ``load`` steps' generic ``"Run MTLoaderAgent"``
description is the default -- pass ``description="..."`` to ``add_step``
for a preview that reads like the :doc:`coordinator` examples do.

Use the coordinator when the chain is known and should be reproducible.
