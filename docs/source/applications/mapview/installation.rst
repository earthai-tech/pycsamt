.. _applications-mapview-installation:

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

Preloading Data
---------------

``--data`` is the quickest way to open straight onto a survey — point it at a
folder of EDI files (or a parent folder with one subfolder per line):

.. code-block:: bash

   pycsamt-mapview --data path/to/edi_folder

Behind the scenes this builds a :class:`~pycsamt.map.MapView` from the folder
and seeds the app with it, so the map is populated the moment the browser
opens.  You can still load or add more lines from the UI afterwards.

Verify
------

After launch, the app opens on the :ref:`welcome screen
<applications-mapview>` with a **Load Survey Lines - Start** button.  If the
browser does not open (headless or remote hosts), browse to
``http://127.0.0.1:8770`` manually.  If the page is blank or the port is busy,
see :doc:`troubleshooting`.
