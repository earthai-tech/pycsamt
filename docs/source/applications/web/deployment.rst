.. _applications-web-deployment:

Deployment Notes
================

The web app is designed for **local and internal-network** use: one person on
their own machine, or a small team sharing a running instance on a trusted
network.  This page covers the supported patterns and the security
expectations that come with them.

.. important::

   The pyCSAMT web app is a research and review tool, not a hardened public web
   service.  It has no built-in authentication or multi-tenant isolation.  Run
   it on a machine or network you trust, and do not expose it directly to the
   public internet.

Local-Only (Default)
--------------------

The default binds to the loopback address, so the app is reachable only from
the same machine:

.. code-block:: bash

   pycsamt-web

This is the right mode for single-user work.  Nothing is exposed to the
network, and the browser opens automatically to ``http://127.0.0.1:8050`` (or
the next free port).

Internal-Network Access
-----------------------

To let colleagues on the same trusted network open the app, bind to all
interfaces:

.. code-block:: bash

   pycsamt-web --host 0.0.0.0 --port 8060

Others then reach it at ``http://<your-machine-ip>:8060``.  Use this on a
trusted LAN only.  Anyone who can reach the address can load data and run
workflows — there is no login.

Choosing A Port
---------------

* ``--port N`` sets the preferred port (default ``8050``).
* If the preferred port is busy, the launcher selects a free port automatically
  and prints it — unless you pass ``--strict-port``, which fails instead.
* ``--port 0`` always asks the operating system for a free port.

.. code-block:: bash

   pycsamt-web --port 8050 --strict-port   # fixed port, fail if taken
   pycsamt-web --port 0                     # always a free port

``--strict-port`` is the right choice behind a reverse proxy or when a firewall
rule expects a specific, stable port.

Behind A Reverse Proxy
----------------------

For a more permanent internal deployment, run the app on a fixed host and port
with a strict port, and place a reverse proxy (for example nginx) in front of
it to handle TLS and access control:

.. code-block:: bash

   pycsamt-web --host 127.0.0.1 --port 8060 --strict-port --no-browser

Bind to ``127.0.0.1`` and let the proxy be the only thing that listens on the
external interface.  Terminate TLS and enforce authentication at the proxy —
the app itself provides neither.

.. note::

   Do **not** use ``--debug`` for shared or proxied deployments.  Debug mode
   enables developer tools and exposes callback tracebacks in the browser,
   which leak internal detail.  Keep it for local development only.

Cache And Temporary Files
-------------------------

Uploaded data and intermediate results are held in the app's cache while the
survey is loaded.  Two practical consequences:

* uploaded files are written to a server-side cache/temporary area so the
  callbacks can read them — treat the host machine as holding a copy of any
  data you load;
* large multi-line surveys use correspondingly more memory and cache space, so
  size the host machine to the surveys you expect to load.

Because uploads land on the server, only run a shared instance on a machine
that is allowed to hold the survey data being reviewed.

Where State Lives
-----------------

* **Browser session** — auto-saved workflow state and the AI **API key** live
  in the user's browser (``localStorage``), per browser and per machine.  They
  are *not* stored on the server.
* **Downloaded session JSON** — the portable copy of a session, produced by the
  **Session** drawer, that a user can carry to another machine.
* **Server cache** — the loaded survey and intermediate products for the
  running process.

This split matters for shared deployments: each user's browser holds their own
key and view, while the loaded survey is shared by the running server process.

Security Expectations
---------------------

* Run on trusted machines and networks only; there is no authentication.
* Do not expose the raw app to the public internet; front it with a proxy that
  provides TLS and authentication if remote access is required.
* Uploaded survey data is cached server-side — host it accordingly.
* AI API keys are per-user in the browser; they are not written to the server
  and are not part of shared exports.
* Keep ``--debug`` off outside local development.

Framework
---------

For reference when planning a deployment, the web layer is built on **Dash 4**,
**Dash Bootstrap Components 2**, and **Plotly 6**, over NumPy, SciPy,
Matplotlib, and pyproj (the same stack reported in the **About** dialog).  It
runs anywhere that stack installs — Windows, macOS, and Linux.

Next Steps
----------

* :doc:`installation` -- the full launch option list.
* :doc:`troubleshooting` -- port, upload, and callback problems.
