Agent Coordinator
=================

:class:`pycsamt.agents.AgentCoordinator` executes an ordered list of agents as
a named workflow.  It provides step registration, dependency mapping,
checkpointing, dry-run previews, and workflow-level cost aggregation.

Basic pattern
-------------

Each step has a unique name, an agent instance, and an optional ``input_fn``.
The ``input_fn`` receives the accumulated previous results and returns the
input dictionary for the next agent.

.. code-block:: python

   from pycsamt.agents import AgentCoordinator, MTLoaderAgent, DataQCAgent

   coord = AgentCoordinator("mt_qc")
   coord.add_step("load", MTLoaderAgent(), description="Load EDIs")
   coord.add_step(
       "qc",
       DataQCAgent(),
       input_fn=lambda r: {"sites": r["load"]["sites"]},
       description="Compute station QC",
   )

   result = coord.execute({"path": "/data/EDIs"})
   print(result.summary)

Dry-run preview
---------------

Use dry-run execution when a workflow should be inspected before files are
read, written, or processed.

.. code-block:: python

   result = coord.execute({"path": "/data/EDIs"}, dry_run=True)
   print(result["plan"])

Workflow recipes
----------------

Common agent chains include:

.. code-block:: text

   ContextInputAgent -> MTLoaderAgent -> DataQCAgent
   MTLoaderAgent -> DataQCAgent -> StaticShiftAgent -> ReportAgent
   MTLoaderAgent -> DataQCAgent -> StaticShiftAgent -> PhaseAnalysisAgent
   MTLoaderAgent -> DataQCAgent -> StaticShiftAgent -> Occam2DAgent
   MTLoaderAgent -> DataQCAgent -> DenoisingAgent -> AIInversionAgent

Custom agents
-------------

Custom agents inherit from :class:`pycsamt.agents.BaseAgent` and return
:class:`pycsamt.agents.AgentResult`.

.. code-block:: python

   from pycsamt.agents import BaseAgent, AgentResult

   class MyCustomAgent(BaseAgent):
       SYSTEM_PROMPT = "You are a custom MT processing expert."

       def __init__(self, *, api_key=None, model=None, llm_provider="claude"):
           super().__init__(
               "MyCustomAgent",
               api_key=api_key,
               model=model,
               llm_provider=llm_provider,
               section_preset="pseudosection",
           )

       def execute(self, input_data):
           sites = input_data.get("sites")
           interpretation = self.query_llm(
               f"Interpret this survey object: {sites}",
               max_tokens=200,
           )
           return AgentResult(
               status="success",
               summary="Custom agent complete.",
               data={"result": 42},
               warnings=[],
               llm_interpretation=interpretation,
               elapsed_seconds=0.1,
               cost_estimate_usd=self._last_cost,
           )

Register the custom step like any built-in agent:

.. code-block:: python

   from pycsamt.agents import AgentCoordinator, MTLoaderAgent

   coord = AgentCoordinator("custom_workflow")
   coord.add_step("load", MTLoaderAgent(), description="Load EDIs")
   coord.add_step(
       "custom",
       MyCustomAgent(),
       input_fn=lambda r: {"sites": r["load"]["sites"]},
       description="Run custom survey analysis",
   )

   result = coord.execute({"path": "/data/EDIs"})

Related API
-----------

See :class:`pycsamt.agents.AgentCoordinator`,
:class:`pycsamt.agents.WorkflowStep`, :class:`pycsamt.agents.BaseAgent`, and
:class:`pycsamt.agents.AgentResult`.
