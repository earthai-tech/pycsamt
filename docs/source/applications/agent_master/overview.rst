.. _applications_agent_master_overview:

Overview
========

What Agent Master Does
--------------------------

* **Understands plain-language requests** and routes them to the correct
  workflow — QC, static-shift, denoising, phase/strike analysis, tipper,
  rotation, frequency decimation, sensitivity, inversion prep, full pipelines,
  reports, and code generation.
* **Asks for parameters** when a workflow needs them, through small in-chat
  forms, and **asks which lines** to process when several are loaded.
* **Runs the real pyCSAMT agents** — the same code as the GUI and web app — and
  streams progress step by step.
* **Returns products**: figures (collected in a Figures panel), reports,
  validated standalone Python scripts, and a step-by-step trace with timing and
  cost for provenance.
* **Uses your choice of LLM** — Claude, OpenAI, Gemini, or DeepSeek — with the
  API key kept on your machine.

When To Use It
------------------

Reach for Agent Master when:

* you want to run a multi-step workflow by describing the goal, not by
  operating each page;
* you want a **reproducible script** for a workflow you just ran;
* you want a quick **report** or a guided **inversion preparation**;
* you are learning what pyCSAMT can do and want to ask it directly.

For dense manual review — station-by-station QC, correction previews, map and
3-D inspection — the :doc:`desktop <../desktop/index>`, :doc:`web
<../web/index>`, and :doc:`MapView <../mapview/index>` surfaces are still the
right tools. Results are interchangeable across all of them.

Run Agent Master
--------------------

.. code-block:: bash

   pycsamt-agent

The launcher starts a local server on ``http://127.0.0.1:8765`` and opens the
app in your browser. The module form works everywhere the package is
importable:

.. code-block:: bash

   python -m pycsamt.app.agent_master

See :doc:`installation` for host, port, and browser options, and
:doc:`llm_configuration` for connecting a model.

The Interface At A Glance
-----------------------------

.. figure:: ../../_static/applications/agent_master/agent-master-interface.png
   :alt: The Agent Master interface: top bar, history sidebar, chat area, input
   :class: pycsamt-screenshot

   The main window: the command bar (**Load EDI**, **Save**, theme, help,
   **Settings**), the **History** sidebar (Chat / Session tabs, pinned prompts,
   and a Figures panel), the chat area with suggestion chips, and the
   natural-language input at the bottom.

.. seealso::

   :doc:`/user_guide/agents/orchestrator`
       How the orchestrating agent behind this app plans and runs
       multi-step workflows — the guide to the brain, while this section
       documents the app surface.
