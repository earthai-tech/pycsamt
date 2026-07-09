Troubleshooting
===============

Port already in use
-------------------

MapView defaults to port ``8770``. If another process (or a previous
MapView instance) holds it, launch on a different port:

.. code-block:: bash

   pycsamt-mapview --port 8780

No browser opens
----------------

Some environments (WSL, remote shells, servers) cannot open a browser.
Start with ``--no-browser`` and browse to the printed URL manually —
``http://127.0.0.1:8770`` by default. From another machine, launch with
``--host 0.0.0.0`` and use the server's address.

The 3-D view is blank
---------------------

The 3-D scene needs WebGL. If the station map works but the 3-D canvas
stays empty, enable hardware acceleration in the browser (or try
another browser) — remote-desktop sessions often disable WebGL.

Large surveys feel slow
-----------------------

* Toggle survey lines off in the line picker while you focus on one
  profile.
* Turn *Smooth sections* off, or lower the *Section resolution* — the
  interpolation grid dominates the figure cost.
* Hide station labels in dense maps.

Data does not load
------------------

* The EDI tab expects station files with coordinates in their headers;
  files without any location can be placed only after loading stations
  that do have coordinates.
* Inversion results without coordinates need *match coordinates from
  already-loaded EDI stations* enabled — load the survey's EDIs first.
* Sessions are per browser tab: a second tab starts empty by design.

Where nothing else helps
------------------------

Run with ``--debug`` for the Dash traceback in the terminal, and check
the browser console. The app never modifies your data folders, so a
restart is always safe.
