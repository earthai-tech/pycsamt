.. _applications-agent-master-handoff:

Agent Master And The Other Apps
===============================

The conversational, agent-driven experience is not unique to this app — it
resurfaces inside the desktop GUI and the web app as well, because all three
draw on the same :mod:`pycsamt.agents` package: the same :term:`Agent`
classes, and, for full-workflow requests, the same
:class:`~pycsamt.agents.WorkflowOrchestratorAgent` that plans a chain from a
plain-language request.  This page works out how the standalone Agent Master
relates to those other agent surfaces, so you can pick the right one for a
given moment of work and move between them without losing anything.

One Set Of Agents, Several Front Ends
-------------------------------------

Nothing scientific changes between surfaces — routing, step execution, and the
:term:`AgentResult` returned by every step are identical wherever the request
originates.  What changes is framing: how much of the screen the conversation
gets, and what sits next to it.

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Surface
     - What it is
   * - **Agent Master** (this app)
     - The dedicated, full-window conversational app — launched with
       ``pycsamt-agent``.  Best when the chat *is* the workflow.
   * - **Web app → AI Agents**
     - A page inside the :doc:`web app <../web/index>` with an **Agent Runner**
       and a **Chat** tab, alongside all the point-and-click pages.
   * - **Desktop GUI**
     - The :doc:`desktop app <../desktop/index>` opens an **Agent Master**
       panel window (``Ctrl+Shift+A``) alongside its Profile, Map, QC, and
       Inversion panels, so agent review sits next to the manual diagnostics
       it is commenting on.

Because the underlying :term:`agent coordinator` and orchestrator are what
actually run — see :doc:`/user_guide/agents/orchestrator` for how a
plain-language request becomes an ordered agent chain — a survey processed
through one surface behaves identically in another, and the products of a
run (figures, reports, generated code) are interchangeable regardless of
which window produced them.

Choosing A Surface
------------------

The right choice usually follows from what is already open, not from any
difference in capability:

* Use **Agent Master** when you want a distraction-free conversation that plans
  and runs a whole workflow — and when you want generated, reproducible code or
  a report from a description.
* Use the **web app's AI Agents** page when you are already working in the
  browser studio and want to run an agent or ask a question without leaving it
  (see :doc:`../web/processing_pages`).
* Use the **desktop GUI**'s Agent Master panel when your work is mostly manual
  review — QC plots, correction previews, inversion diagnostics — and you want
  agent commentary next to the evidence rather than in a separate app. As the
  desktop guide puts it, treat the agent as an explanation aid: let it point
  you back to the profile or QC view, then make the processing decision from
  the data itself, not from the agent's confidence alone.

Moving Between Them
-------------------

Agent Master and the web app each start their own local server, and the
desktop GUI runs as its own process again, so there is no shared in-memory
session for a running agent to hand off — closing one does not hand a live
conversation to another.  What crosses freely is the same thing every surface
already reads and writes: files.

* Load the **same EDI folder** in whichever surface you move to. The
  :term:`EDI` format and the loading behaviour are identical everywhere, so a
  survey that quality-checks cleanly in the web app opens the same way in
  Agent Master or the desktop GUI.
* Carry **generated scripts** and **exported figures** across freely — a
  script produced by Agent Master's code generation is a standalone pyCSAMT
  program with no dependency on the app that wrote it (see
  :doc:`tools_memory_outputs`), and an exported figure is just an image or
  HTML file.
* Keep the **survey folder** alongside any **saved session** so the working
  context can be rebuilt on the next surface even though the conversation
  itself does not travel with you.

Because each surface binds its own port, you are not limited to one at a
time either: Agent Master listens on ``8765``, the web app on ``8050``, and
MapView on ``8770`` by default, so you can leave several running and move
between browser tabs as the task shifts from chat-driven planning to
point-and-click review.

Failure Recovery
-----------------

The same underlying stack means most handoff problems trace back to one of
three places — the server, the model provider, or a mid-workflow parameter —
regardless of which surface surfaced the failure.

* If Agent Master is not reachable, confirm the ``pycsamt-agent`` process is
  still running and that its terminal shows ``http://127.0.0.1:8765`` (see
  :doc:`installation`); a closed terminal is the most common cause, since it
  is what hosts the server.
* If a request needs a model and none is configured, add a provider and key in
  :doc:`Settings <llm_configuration>` and resend the request — you can still
  browse the interface without a key, but the language-driven planning and
  chat that make Agent Master useful will not run until one is set.
* If a run fails partway through, read the **steps trace** for the failing
  step, adjust the parameters in the in-chat form, and run again (see
  :doc:`workflows_and_agents`).  Because every step returns its own
  :term:`AgentResult`, the trace tells you exactly which step stopped the
  chain rather than leaving you to re-run the whole workflow blind.

Next Steps
----------

* :doc:`../web/processing_pages` -- the web app's AI Agents page.
* :doc:`troubleshooting` -- when a surface will not start or connect.
