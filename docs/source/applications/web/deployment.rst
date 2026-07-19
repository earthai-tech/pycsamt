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

``pycsamt-web`` is a :term:`console script` around ``pycsamt.app.web.app.main``
— the same entry point ``python -m pycsamt.app.web`` calls — so both forms
launch identically.  This is the right mode for single-user work: the
:term:`local server` binds only to ``127.0.0.1``, nothing is exposed to the
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
it to handle :term:`TLS` and access control:

.. code-block:: bash

   pycsamt-web --host 127.0.0.1 --port 8060 --strict-port --no-browser

Bind to ``127.0.0.1`` and let the proxy be the only thing that listens on the
external interface.  Terminate TLS and enforce authentication at the proxy —
the app itself provides neither.

.. note::

   Do **not** use ``--debug`` for shared or proxied deployments.  It enables
   :term:`Dash debug mode`, which exposes callback diagnostics and reload
   behaviour, including tracebacks in the browser — helpful while developing,
   but a detail leak on anything reachable by someone else.  Keep it for
   local development only.

Cache And Temporary Files
-------------------------

Once the browser has parsed and uploaded a survey folder, the resulting
station objects are cached server-side in a ``diskcache.Cache`` directory at
``~/.pycsamt/web_cache`` (or ``$PYCSAMT_DATA/web_cache`` if you set the
:term:`environment variable` ``PYCSAMT_DATA`` to relocate pyCSAMT's data
root), keyed by a per-tab session id and expiring automatically after 24
hours.  Two practical consequences follow directly from that:

* the host machine holds a real, on-disk copy of every survey loaded through
  it — for as long as the cache entry lives — so only run a shared instance
  on a machine that is allowed to hold the data being reviewed, and consider
  where ``PYCSAMT_DATA`` points if the default under the server's home
  directory is not an acceptable location;
* large multi-line surveys, and more concurrent browser tabs, use
  correspondingly more disk and memory, so size the host to the surveys and
  the number of simultaneous users you expect — cached entries are not
  deleted when a tab closes, only when they expire, so a busy shared instance
  accumulates cache until the 24-hour sweep catches up.

If the optional ``diskcache`` dependency is unavailable, the app falls back
to an in-memory cache instead.  That keeps a single-process development
server working, but it does not survive a restart and cannot back the
background-callback manager, so treat it as a degraded mode rather than a
deployment option.

Where State Lives
-----------------

* **Browser, per tab (``sessionStorage``)** — the session id that keys the
  server-side cache above.  It is a UUID generated fresh for each browser
  tab, so closing the tab (not just reloading it) severs the link to that
  tab's cached survey; a new tab gets a new id and an empty cache entry even
  on the same machine.
* **Browser, per machine (``localStorage``)** — auto-saved workflow state,
  the theme choice, and the AI :term:`API key` persist across restarts of
  the same browser on the same machine.  None of this is written to the
  server.
* **Downloaded session JSON** — the portable copy of a session, produced by
  the **Session** drawer, that a user can carry to another machine.  Because
  it is workflow state rather than the cache entry itself, restoring it does
  not require the original tab or its session id.
* **Server disk cache** — the loaded survey and intermediate products
  described above, shared by every tab that presents the same session id and
  isolated between tabs that do not.

This split matters for shared deployments: each user's browser holds only
their own key, theme, and view, while the actual survey data being reviewed
sits on the server, addressed by an id that lives only as long as its tab
does.

Security Expectations
---------------------

* Run on trusted machines and networks only; there is no authentication.
* Do not expose the raw app to the public internet; front it with a proxy that
  provides :term:`TLS` and authentication if remote access is required —
  the app terminates neither itself.
* Uploaded survey data is cached to disk under ``~/.pycsamt/web_cache`` (or
  ``PYCSAMT_DATA``) for up to 24 hours — host and secure that path
  accordingly, the same as any other location holding a copy of the data.
* The AI :term:`API key` is a :term:`secret` kept per-user in the browser's
  ``localStorage``; it is never written to the server and is not part of
  shared exports or the downloaded session JSON.
* Keep ``--debug`` off outside local development.

Framework
---------

For reference when planning a deployment, the web layer is built on
:term:`Dash` 4, Dash Bootstrap Components 2, and :term:`Plotly` 6, over
NumPy, SciPy, Matplotlib, and pyproj (the same stack reported in the
**About** dialog).  It runs anywhere that stack installs — Windows, macOS,
and Linux.

Next Steps
----------

* :doc:`installation` -- the full launch option list.
* :doc:`troubleshooting` -- port, upload, and callback problems.
