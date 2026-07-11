.. _applications-web-exports:

Exports And Reproducibility
===========================

The web app produces two kinds of output: **figures** you save for reports, and
**data products** (corrected data, inversion figures, and session files) you
reuse elsewhere.  Because every page delegates to the pyCSAMT package, anything
you build interactively can also be regenerated in code — the browser is a
convenient surface, not a separate result.

Saving Figures
--------------

**Interactive Plotly figures** — the Map View, 3-D scene, and many section
plots are Plotly figures.  Each has a **modebar** in the top-right corner:

* the **camera** button downloads the current view as a **PNG**;
* pan, zoom, and autoscale let you frame the figure before you snapshot it;
* double-click resets the view.

**Matplotlib page figures** — QC, correction, advanced, forward, interpretation,
and results figures are rendered as images.  Where a page provides an explicit
**Export** / **Export PNG** button, use it rather than a screen capture: it
writes the figure at full resolution through the app's download mechanism.

.. figure:: ../../_static/applications/web/results-section.png
   :alt: The Results View page with an Export PNG button beneath the plot controls
   :class: pycsamt-screenshot

   Pages that render matplotlib figures — like Results View — provide an
   **Export PNG** button so you save the figure itself, not a screenshot.

Page And Tool Exports
---------------------

Several pages and tools export more than a single figure:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Surface
     - What it exports
   * - **Correction page**
     - **Export Corrected Data** writes the corrected response as a downloadable
       table (CSV) reflecting the current correction chain.
   * - **3D Map page**
     - The **Export** panel saves the 3-D scene as a **PNG** or a standalone
       interactive **HTML** file.
   * - **Results View page**
     - **Export PNG** saves the current results tab (section, convergence,
       response, depth map, …).
   * - **Interpretation page**
     - The **Export** category saves interpretation products.
   * - **Pipeline page**
     - The **Export** step writes processed EDIs and products to an export
       folder you choose.
   * - **Batch Export Plots** (Tools)
     - Exports many figures at once for the loaded survey.

Batch Export
------------

For anything more than a few figures, use **Tools → Batch Export Plots**.  It
generates a set of plots for the loaded survey in one pass, which is faster and
more consistent than exporting page by page — useful when preparing a full
figure set for a report.

Session Files
-------------

The **Session** drawer (see :doc:`loading_and_sessions`) downloads a session
JSON that captures your workflow state and settings.  It is the reproducibility
anchor for the browser:

* **Download JSON** archives the current session or shares it with a colleague;
* **Restore Session** re-opens it in any browser running the app;
* the **What is included?** link lists exactly what the file records.

The session JSON stores workflow state, not the raw EDI files.  Keep the
original survey folder alongside the session so the survey can be reloaded when
the session is restored on another machine.

The **Settings** drawer can likewise download its configuration, so computation
defaults and provider settings can be carried between machines.

.. note::

   The AI provider **API key** lives in the browser's ``localStorage`` and is
   deliberately not part of shared exports.  Each user enters their own key in
   **Settings** (see :doc:`navigation`).

Reproducing Web Output In Python
--------------------------------

The web app is a view over the package, so its output is reproducible in code:

* interactive maps and 3-D scenes correspond to the ``pycsamt.map`` façade, so
  a view you build in Map View or 3D Map can be regenerated in a script;
* corrections applied in the Correction chain correspond to the same catalogue
  functions available from Python, so a corrected survey can be rebuilt
  headlessly;
* the **Pipeline** page mirrors the scripted ``load → QC → edit → correct →
  strike → export`` workflow, so a pipeline you tune interactively can be run
  from the Python API on new data.

This parity is what makes the exports trustworthy: a figure in a report can be
tied back to a reproducible computation rather than a one-off screen capture.

Practical Export Habits
-----------------------

* Prefer a page's **Export** button over a screenshot — it saves the figure at
  full resolution.
* Give exported files descriptive names (line, quantity, and settings), for
  example ``L34_phase_tensor_section.png`` or ``modem_run_section.png``.
* Keep the **survey folder**, the **session JSON**, and the exported figures
  together, so any view can be regenerated later.
* For a full figure set, use **Batch Export Plots** rather than exporting each
  page by hand.

Next Steps
----------

* :doc:`deployment` -- where exports and uploads live when the app is shared.
* :doc:`troubleshooting` -- when an export or download does not appear.
