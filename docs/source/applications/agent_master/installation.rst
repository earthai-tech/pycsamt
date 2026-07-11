.. _applications-agent-master-installation:

Installation And Launch
=======================

Agent Master ships with pyCSAMT and uses the same optional application extra as
the desktop, web, and MapView apps — one install covers every application
surface.

Install The App Extra
---------------------

.. code-block:: bash

   pip install "pycsamt[app]"

For development from a source checkout:

.. code-block:: bash

   pip install -e ".[app,dev]"

The extra brings the Dash/Plotly stack the app renders with; the scientific
core and the agents are part of the base install.  To actually *drive* the
LLM-backed features you also need a model provider and API key — see
:doc:`llm_configuration`.

Launch
------

.. code-block:: bash

   pycsamt-agent

The launcher starts a local server on ``http://127.0.0.1:8765`` and opens the
app in your default browser.  The module form launches the same app:

.. code-block:: bash

   python -m pycsamt.app.agent_master

Leave the terminal open while you work — it hosts the server.  Press ``Ctrl+C``
to stop.

Options
-------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Option
     - Effect
   * - ``--host HOST``
     - Bind address (default ``127.0.0.1``, local-only).  Use ``0.0.0.0`` to
       reach the app from another machine on a trusted network.
   * - ``--port PORT``
     - HTTP port (default ``8765``).
   * - ``--no-browser``
     - Start the server without opening a browser tab.
   * - ``--debug``
     - Enable Dash debug / hot-reload — useful only when developing the app
       itself.

Examples:

.. code-block:: bash

   pycsamt-agent --port 8770          # serve on a different port
   pycsamt-agent --no-browser         # headless / remote host
   pycsamt-agent --host 0.0.0.0       # allow LAN access (trusted networks only)

The Welcome Screen
------------------

On first load, Agent Master shows a welcome screen with a single call to
action.

.. figure:: ../../_static/applications/agent_master/welcome-page.png
   :alt: The Agent Master welcome screen
   :class: pycsamt-screenshot

   The welcome screen summarises the four things Agent Master does — **Load
   EDI**, **Chat & Plan**, **AI Inversion**, and **Reports** — and starts the
   session with **Start Agent Master**.

Verify
------

After launch, confirm that:

* the terminal reports the server is running on ``http://127.0.0.1:8765``;
* the browser opens to the welcome screen (or copy the URL manually on a
  headless host);
* **Settings** (the gear icon) opens and shows the LLM provider options.

If the browser does not open, the page is blank, or the port is busy, see
:doc:`troubleshooting`.

Next Steps
----------

* :doc:`llm_configuration` -- connect Claude, OpenAI, Gemini, or DeepSeek.
* :doc:`welcome_and_chat` -- the interface, loading data, and chatting.
