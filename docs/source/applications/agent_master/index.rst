.. _applications-agent-master:

Agent Master
============

Agent Master is the **conversational** application of pyCSAMT. You describe a
CSAMT / AMT / MT task in plain language and it plans the workflow, asks for any
parameters it needs, runs the right pyCSAMT agents, and returns figures,
reports, reproducible code, and a full provenance trace. Its tagline says it
plainly: *describe your workflow — I handle the rest.*

It complements the other surfaces rather than replacing them. The
:doc:`desktop GUI <../desktop/index>`, :doc:`web app <../web/index>`, and
:doc:`MapView <../mapview/index>` are point-and-click; Agent Master is the
chat-driven front end over the same agents and the same scientific core, for
when it is faster to *say* what you want than to click through to it.

.. figure:: ../../_static/applications/agent_master/agent-master-walkthrough.gif
   :alt: Animated walkthrough of an Agent Master session, from welcome to results
   :class: pycsamt-screenshot

   A full Agent Master session at a glance: launch, configure the LLM, load
   EDI data, ask *what can you do?*, run a workflow (with line selection and
   generated figures), generate a validated script, and get a report with a
   sensitivity analysis. Each step below is documented in detail on the
   following pages.

A good first pass, in order, is **Overview → Installation → LLM
Configuration → Welcome & Chat → Workflows & Agents → Tools, Memory &
Outputs → Agent Master & Other Apps → Troubleshooting**.

.. toctree::
   :numbered: 4
   :maxdepth: 1
   :class: pycsamt-guide-toc

   overview
   installation
   welcome_and_chat
   workflows_and_agents
   llm_configuration
   tools_memory_outputs
   handoff_from_apps
   troubleshooting
