Foundation And Survey Intake Agents
===================================

These components define the agent execution contract and the first steps of
most pyCSAMT workflows: parse intent, load data, and coordinate downstream
processing.

.. _agent-base-agent:

BaseAgent
---------

``BaseAgent`` is the shared base class for all agent classes.  It resolves the
effective LLM configuration, exposes ``query_llm()``, tracks per-call cost,
provides JSON extraction helpers, and injects pyCSAMT plotting/style helpers
so figures produced by agents remain consistent with the rest of the package.

Use it when building custom agents:

.. code-block:: python
   :linenos:

   from pycsamt.agents import BaseAgent, AgentResult

   class MyAgent(BaseAgent):
       def __init__(self):
           super().__init__("MyAgent")

       def execute(self, input_data):
           self._last_cost = 0.0
           return AgentResult(
               status="success",
               summary="Custom step complete.",
               data={"value": 1},
               cost_estimate_usd=self._last_cost,
           )

.. _agent-result:

AgentResult
-----------

``AgentResult`` is the standard return object for every agent.  It contains
status, summary, structured data, warnings, optional LLM interpretation, run
time, cost estimate, and failure details.

The result supports dict-like access to ``data``:

.. code-block:: python
   :linenos:

   result = agent.execute({"path": "/data/EDIs"})

   if result:
       print(result.summary)
       print(result["sites"])
       print(result.get("qc_table"))
   else:
       print(result.error)
       print(result.error_fix_hint)

.. _agent-context-input:

ContextInputAgent
-----------------

``ContextInputAgent`` converts natural-language workflow requests into a
structured configuration dictionary.  It can call an LLM when configured, but
also includes a regex fallback for common request patterns such as loading a
path, selecting a workflow, and extracting period ranges.

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

.. code-block:: python
   :linenos:

   from pycsamt.agents import ContextInputAgent

   agent = ContextInputAgent()
   result = agent.execute({
       "request": "Load EDIs from /data/WILLY, QC them, period 0.001 to 10 s",
   })

   print(result["config"])
   print(result.llm_interpretation)

Use this agent at the start of natural-language workflows or before building
an ``AgentCoordinator`` from user text.

.. _agent-mt-loader:

MTLoaderAgent
-------------

``MTLoaderAgent`` loads MT, AMT, or CSAMT survey data into a pyCSAMT ``Sites``
object.  It accepts a single file, a directory, a list of paths, an existing
``Sites`` object, or an EDI collection supported by the lower-level loading
utilities.

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

.. code-block:: python
   :linenos:

   from pycsamt.agents import MTLoaderAgent

   result = MTLoaderAgent().execute({
       "path": "/data/WILLY_EDIs",
   })

   print(result["n_stations"])
   print(result["quality_table"])

Use this agent whenever a workflow needs a validated ``Sites`` object before
QC, static-shift correction, phase analysis, or inversion preparation.

.. _agent-coordinator-class:

AgentCoordinator
----------------

``AgentCoordinator`` is not a ``BaseAgent`` subclass, but it is the central
workflow runner for explicit multi-agent pipelines.  It registers ordered
steps, maps previous results into the next step with ``input_fn``, supports
dry-run previews, and aggregates cost.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AgentCoordinator, ContextInputAgent, MTLoaderAgent

   coord = AgentCoordinator("load_preview")
   coord.add_step("parse", ContextInputAgent())
   coord.add_step(
       "load",
       MTLoaderAgent(),
       input_fn=lambda r: {
           "path": (r["parse"].get("config") or {}).get("data_path", "")
       },
   )

   result = coord.execute(
       {"request": "Load /data/WILLY_EDIs"},
       dry_run=True,
   )

   print(result["plan"])

Use the coordinator when the chain is known and should be reproducible.
:html_theme.sidebar_secondary.remove:
