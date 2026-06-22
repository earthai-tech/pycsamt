.. _web_application:

Web Application
===============

The pyCSAMT web application is a Dash/Plotly interface for interactive
survey work in a browser. It mirrors the desktop application's main
workflow zones while using web-native layout, callbacks, cache handling,
and Plotly figures.

.. note::

   This page is a stub. A full deployment and user guide will be added
   later after the desktop and pipeline documentation are complete.

Run The Web App
---------------

Install the application extra:

.. code-block:: bash
   :linenos:

   pip install "pycsamt[app]"

Start the web app launcher:

.. code-block:: bash
   :linenos:

   pycsamt-web

The launcher starts the local Dash server and opens the app in your default
browser. By default it binds to ``127.0.0.1`` and prefers port ``8050``:

.. code-block:: text
   :linenos:

   http://127.0.0.1:8050

If port ``8050`` is already in use, ``pycsamt-web`` automatically chooses a
free port and opens the browser at the correct address.

Advanced launch options:

.. code-block:: bash
   :linenos:

   pycsamt-web --no-browser
   pycsamt-web --port 8060
   pycsamt-web --port 0
   pycsamt-web --host 0.0.0.0
   pycsamt-web --debug

Use ``--host 0.0.0.0`` only when you intentionally want the app to be reachable
from another machine on the local network.

Current Architecture
--------------------

The web application lives under :mod:`pycsamt.app.web`.

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Package area
     - Role
   * - ``app.py``
     - Dash application factory, Flask asset routes, callback registration,
       and the ``pycsamt-web`` launcher.
   * - ``layout.py``
     - Top-level browser layout and navigation composition.
   * - ``pages/``
     - Page modules for home, QC, correction, forward modelling,
       inversion, interpretation, pipeline, TDEM, agents, and advanced
       tools.
   * - ``callbacks/``
     - Dash callbacks for data loading, maps, profiles, frequency controls,
       QC, corrections, forward modelling, inversion, interpretation,
       pipelines, TDEM, agents, navigation, themes, and session state.
   * - ``cache.py``
     - Cache and optional diskcache integration for long-running or
       background-capable operations.
   * - ``assets/``
     - Web-specific CSS and JavaScript assets.

Feature Areas To Document
-------------------------

The detailed web manual should later cover:

* local development startup and debug mode;
* browser layout and navigation;
* loading survey data and maintaining session state;
* station maps and profile panels;
* QC, correction, forward, inversion, interpretation, pipeline, and TDEM
  pages;
* agent page behavior and long-running callbacks;
* cache configuration and diskcache behavior;
* running behind a local network address or deployment server;
* differences between the web and desktop application surfaces.

Relationship To The Desktop App
-------------------------------

The web app is not a separate scientific workflow implementation. Its
callbacks should delegate to the same pyCSAMT processing, pipeline, agent,
and controller logic used by the rest of the project whenever possible.
This keeps numerical behavior consistent across Python scripts, CLI runs,
desktop sessions, and browser sessions.

Related API
-----------

The application API stub is available at :doc:`../api/app`.
