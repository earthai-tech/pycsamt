.. _applications-agent-master-handoff:

Agent Master And The Other Apps
===============================

The conversational, agent-driven experience appears in more than one place in
pyCSAMT.  This page explains how the standalone Agent Master relates to the
agent surfaces inside the desktop GUI and the web app, so you can pick the right
one and move between them.

One Set Of Agents, Several Front Ends
-------------------------------------

Agent Master, the web app's **AI Agents** page, and the desktop GUI's agent
panel are all front ends over the **same** pyCSAMT agents and orchestrator.
They differ only in framing:

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
     - The :doc:`desktop app <../desktop/index>` exposes the agents from within
       its native workflow windows.

Because they share the agents and the scientific core, a survey processed in
one surface behaves identically in another, and the products (figures, reports,
generated code) are interchangeable.

Choosing A Surface
------------------

* Use **Agent Master** when you want a distraction-free conversation that plans
  and runs a whole workflow — and when you want generated, reproducible code or
  a report from a description.
* Use the **web app's AI Agents** page when you are already working in the
  browser studio and want to run an agent or ask a question without leaving it
  (see :doc:`../web/processing_pages`).
* Use the **desktop GUI** when your work is mostly manual review and you want
  agent help in context.

Moving Between Them
-------------------

The surfaces do not share a live in-memory session, but they *do* share data
and outputs, so the handoff is through files:

* load the **same EDI folder** in whichever surface you move to — loading is
  identical everywhere;
* carry **generated scripts** and **exported figures** across freely;
* keep the survey folder alongside any saved session so the context can be
  rebuilt.

Each surface runs as its own local server on its own port — Agent Master on
``8765``, the web app on ``8050``, MapView on ``8770`` — so you can even keep
more than one open at once and move between browser tabs.

Failure Recovery
----------------

* If Agent Master is not reachable, confirm the ``pycsamt-agent`` process is
  running and that its terminal shows ``http://127.0.0.1:8765`` (see
  :doc:`installation`).
* If a request needs a model and none is configured, add a provider and key in
  :doc:`Settings <llm_configuration>` and resend.
* If a run fails mid-workflow, read the step trace, adjust the parameters in the
  in-chat form, and run again (see :doc:`workflows_and_agents`).

Next Steps
----------

* :doc:`../web/processing_pages` -- the web app's AI Agents page.
* :doc:`troubleshooting` -- when a surface will not start or connect.
