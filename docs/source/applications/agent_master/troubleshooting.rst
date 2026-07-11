.. _applications-agent-master-troubleshooting:

Troubleshooting
===============

Common problems when launching and using Agent Master, grouped by where they
occur.

Launch
------

**Dash / app dependencies missing.**
    Install the application extra and relaunch::

        pip install "pycsamt[app]"

**Port 8765 already in use.**
    Another process (or a previous Agent Master) holds the port.  Launch on a
    different one::

        pycsamt-agent --port 8770

**The browser does not open.**
    On headless or remote hosts there may be no browser to open.  Launch with
    ``--no-browser`` and browse to the printed URL (``http://127.0.0.1:8765``
    by default) yourself.  From another machine on a trusted network, use
    ``--host 0.0.0.0`` and the server's address.

**The terminal exits / nothing serves.**
    The terminal running ``pycsamt-agent`` hosts the server — keep it open.
    Closing it stops the app.

Model Provider
--------------

**A request needs a model / "no API key".**
    The language-driven parts need a provider and key.  Open **Settings**,
    choose a **provider**, paste the **API key**, set the **model** and the
    **active provider**, then **Save API keys** (see
    :doc:`llm_configuration`).  Resend the request.

**Provider connection errors.**
    Check the key, the selected model id, and network access to the provider.
    A proxy or firewall may block the provider's endpoint.  Switching the
    **active provider** to one you have a working key for is a quick workaround.

Data And Workflows
------------------

**Nothing runs / the assistant asks for data.**
    Load an EDI survey first with **Load EDI** — every workflow operates on a
    loaded survey (see :doc:`welcome_and_chat`).

**A workflow keeps asking which lines to process.**
    That is expected when several lines are loaded: choose lines in the
    **Select line(s) to process** panel and **Run selected**, or **Run all
    lines** (see :doc:`workflows_and_agents`).

**A workflow fails partway.**
    Read the **steps trace** in the response for the failing step, adjust the
    parameters in the in-chat form (method, thresholds, period limits), and run
    again.

**Generated code will not import.**
    Agent Master validates generated scripts before showing them, so imports
    should resolve.  If a script still fails, confirm the pyCSAMT version in the
    environment you run it in matches the one that generated it, and that the
    data path in the script header points at your survey.

Session And Browser
-------------------

**A new tab starts empty.**
    Sessions are per browser context.  Use **Save** and the **Session** tab to
    persist and restore a conversation rather than relying on a second tab.

**Stale state after an update.**
    Start a **New Chat**, or reload the page; re-load the survey if the loaded
    badge looks wrong.

When You Need More Detail
-------------------------

Relaunch with dev tools for tracebacks in the terminal and browser::

   pycsamt-agent --debug

Use ``--debug`` only for local diagnosis; run without it for normal use.

See Also
--------

* :doc:`installation` -- dependencies and launch options.
* :doc:`llm_configuration` -- provider and key setup.
* :doc:`workflows_and_agents` -- line selection and parameter forms.
