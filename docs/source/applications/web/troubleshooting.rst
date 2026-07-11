.. _applications-web-troubleshooting:

Troubleshooting
===============

This page collects the common problems when launching and using the web app,
grouped by where they occur: launch, loading, pages and callbacks, plotting,
agents, and the browser.  For deployment-specific concerns (ports behind a
proxy, cache location, security), see :doc:`deployment`.

Launch Problems
---------------

**Dash / DBC not installed.**
    The launcher stops with a message telling you to install ``pycsamt[app]``.
    Install the extra and relaunch::

        pip install "pycsamt[app]"

**Port already in use.**
    By default the launcher selects the next free port automatically and prints
    it (``requested port 8050 was busy; using 8051``) — use the printed URL.
    If you passed ``--strict-port``, it fails instead: choose another port with
    ``--port`` or drop ``--strict-port``.  To always get a free port, use
    ``--port 0``.

**The browser did not open.**
    The server still started — copy the ``http://…`` URL from the terminal.
    If you launched with ``--no-browser``, that is expected.  On a headless or
    remote host there may be no browser to open; use ``--no-browser`` and reach
    the URL yourself.

**Nothing happens / the terminal exits.**
    Confirm the terminal that runs ``pycsamt-web`` stays open — it hosts the
    server.  Closing it stops the app.

Loading Problems
----------------

**"No data selected" or nothing loads.**
    Confirm you pointed at a folder that actually contains ``.edi``, ``.avg``,
    or ``.j`` files.  Browse-folder selection is **recursive**, so a survey
    root with one subfolder per line works well; an empty or wrong folder loads
    nothing.  Review **Preflight details** in the loader before pressing
    **Load Survey**.

**Lines are wrong or missing.**
    Line grouping comes from folder names by default.  If stations landed in
    the wrong line, open the **Lines** drawer and switch line assignment to
    *Auto-detect from ID* (then **Detect Lines**) or *Edit / Rename* (then
    **Apply Renames**).  See :doc:`loading_and_sessions`.

**A page ignores some lines.**
    That is by design: only **Active** lines flow through QC, Correction, and
    Profile.  Check the **Lines** drawer and activate the lines you need
    (**All** activates every line).

**Upload is slow or the browser struggles with a big survey.**
    Very large multi-line surveys take time to read, parse, and hold in cache.
    Load the lines you need, or load a subset and use **+Lines** to add more.
    On a shared host, size the machine to the surveys you expect (see
    :doc:`deployment`).

Page And Callback Problems
--------------------------

**A page is blank or a plot did not appear.**
    Most pages do not draw until you press the primary action — **Plot**,
    **Generate**, **Run …**, or **Preview**.  Press it after setting the
    controls.  If it still does not draw, check that a survey is loaded and that
    the page's **Lines/Stations** pickers are not empty.

**A processing step reports an error in the log.**
    Read the **Log** tab.  Failures are surfaced explicitly — for example
    ``ERROR: SVD did not converge`` from a static-shift solve.  Adjust the
    method's parameters (window, kernel, station selection) and re-run the step.

**A control seems stuck.**
    Reload the page in the browser.  Because sessions auto-save to the browser,
    a reload restores your place without losing the loaded survey; if it does
    not, restore your downloaded session JSON from the **Session** drawer.

Plotting Problems
-----------------

**Stations plot in the wrong place on the map.**
    Set the correct **Coordinate System** in Map View.  Standard EDI files are
    usually geographic latitude/longitude, but project files may use UTM or
    another EPSG code.

**Contours look blocky or misleading.**
    Contours interpolate over station locations.  Reduce the number of levels,
    switch back to the station overview, and confirm that station spacing
    supports interpolation.

**Basemap tiles do not appear.**
    Basemap tiles need optional geospatial dependencies and network access.
    The station and contour maps work without tiles.

**Profile tabs are empty.**
    Select a station with valid impedance data and enough frequency samples.
    Some tabs — the tipper view in particular — require the matching data to
    exist in the loaded files.

Agent Problems
--------------

**Agents or chat do not respond / "no API key".**
    LLM-backed agents need a provider and API key.  Open **Settings**, choose a
    **provider**, paste the **API key**, pick a **model**, and use **Test
    Connection**.  The key is stored in the browser only.

**"Test Connection" fails.**
    Check the key, the selected model, and network access to the provider.  A
    proxy or firewall may block the provider's endpoint.

**Non-LLM (processing) agents still work without a key.**
    Processing agents such as QC Quicklook, Dimensionality, and Static Shift do
    not need an LLM; only the LLM-backed agents and chat do.

Browser Problems
----------------

**Theme or settings look stale after an update.**
    The theme, API key, and last view live in the browser's ``localStorage``.
    If they look inconsistent after upgrading the app, a hard refresh (or
    clearing site data for the app's origin) resets browser state.  You can
    re-import a session with the **Session** drawer afterwards.

**Downloads do not appear.**
    Figure and data downloads use the browser's normal download flow.  Check
    the browser's download bar and its download folder; if a download prompt is
    blocked, allow downloads for the app's address.

When To Use Debug Mode
----------------------

If a problem is reproducible and you need detail, relaunch with dev tools:

.. code-block:: bash

   pycsamt-web --debug

Debug mode exposes callback tracebacks in the browser.  Use it locally to
diagnose a problem, then run without ``--debug`` for normal and shared use (see
:doc:`deployment`).

Still Stuck?
------------

* Re-check :doc:`installation` for dependencies and launch options.
* Re-check :doc:`loading_and_sessions` for line assignment and sessions.
* Confirm the survey is healthy on the Home dashboard (see :doc:`navigation`).
