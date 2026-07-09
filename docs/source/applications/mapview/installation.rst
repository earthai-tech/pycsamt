Installation And Launch
=======================

MapView ships with the same optional application extra as the desktop
and web apps — one install covers all pyCSAMT application surfaces:

.. code-block:: bash

   pip install "pycsamt[app]"

For development from a source checkout:

.. code-block:: bash

   pip install -e ".[app,dev]"

The extra brings the Dash/Plotly stack the app renders with; the
scientific core (readers, ``pycsamt.map``) is part of the base install.

Launching
---------

.. code-block:: bash

   pycsamt-mapview

which is equivalent to:

.. code-block:: bash

   python -m pycsamt.app.mapview

By default the server binds ``127.0.0.1:8770`` and opens your browser.
All options:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Effect
   * - ``--port 8780``
     - Serve on another port (default ``8770``).
   * - ``--host 0.0.0.0``
     - Bind all interfaces, e.g. to reach the app from another machine
       on the network.
   * - ``--data <folder>``
     - Preload an EDI folder (one or more survey lines) at startup.
   * - ``--no-browser``
     - Start the server without opening a browser tab.
   * - ``--debug``
     - Dash debug mode with hot reload — useful when developing the app
       itself.

Verify
------

After launch, the welcome panel offers to load data; if the page is
blank, browse to ``http://127.0.0.1:8770`` manually (see
:doc:`troubleshooting`).
