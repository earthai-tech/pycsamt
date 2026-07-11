.. _agents_gallery:

Agents
------

:mod:`pycsamt.agents` provides small, composable agents for interpreting a
request, planning a workflow, running scientific tools, and coordinating
multi-step processing.  The examples in this section run entirely **offline**:
each one executes inside
:meth:`AGENT_CONFIG.offline() <pycsamt.api.agents.AgentConfig.offline>`, so the
agents use deterministic rule-based logic and local numerical code — no API
key, no network request, and zero cost — even when provider keys are present in
the environment.

The sequence progresses from the smallest building blocks to the high-level
:class:`~pycsamt.agents.AgentMaster` front door, and every example ends with a
figure:

* route a plain-language request to an intent (confidence chart);
* extract a structured workflow configuration (parsed-config table);
* inspect the offline model catalogue (architecture/depth chart);
* run a local 1-D forward model (earth model + synthetic sounding);
* validate a workflow plan before execution (pass/fail report card);
* preview a coordinated multi-agent chain (flow diagram);
* plan a complete workflow through ``AgentMaster`` (pipeline diagram).
