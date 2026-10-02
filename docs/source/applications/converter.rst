.. _applications-converter:

Format Converter
=================

Every classical solver writes its own file layout, EDI and EMTF-XML disagree
about how a transfer function should be serialized, and a borehole, a
geology legend, a fault trace, or a target list each need their own small
JSON container before they can sit next to a resistivity model in
:doc:`Map View </applications/mapview/index>` or the desktop's 3-D panel.
**pyCSAMT Format Studio** — the standalone converter app, ``pycsamt.app.converter``
— is the one small tool whose entire job is turning all of that into a file
another pyCSAMT surface can open directly.

It is deliberately not a scaled-down desktop GUI. The full
:doc:`desktop suite </applications/desktop/index>` loads a survey, processes
it, models it, and only then offers a **Format Converter** tool (from its
**Tools** menu) that exports an *already-loaded* survey to EDI/CSV/JSON —
EMTF-XML, PCBH, PCGL, PCGS, and PCPT are not reachable from the desktop app
at all. This standalone app skips the survey entirely: point it at a file or
folder, choose an output, click the one orange button, and it is done. That
focus is also what makes it freezable — see
:ref:`applications-converter-compile` — into a double-click binary you can
hand to a field engineer or office colleague who does not want a Python
environment at all.

.. figure:: /_static/applications/converter/converter-inversion-to-pcsf.png
   :alt: pyCSAMT Format Studio open on the Inversion to PCSF/PCSM page
   :class: pycsamt-screenshot
   :align: center
   :width: 92%

   The Inversion → PCSF/PCSM page after converting a real ModEM 3-D run. The
   sidebar lists all nine tools in the same order as the **Tools** menu; the
   source field auto-detected a ModEM working directory ("Detected: solver
   (modem)"); the **Advanced options** panel exposes solver-specific
   overrides two-up so they never force horizontal scrolling; and the shared
   log dock at the bottom — one dock every page reports to, not nine
   separate ones — shows the result of the *previous* action taken on the
   currently-loaded file (here, a follow-up **Info** read-back: a
   :math:`41\times 50\times 288` grid, 125 stations, resistivity spanning
   four orders of magnitude).

What it converts, end to end:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Page
     - Converts
   * - Inversion → PCSF/PCSM
     - Occam2D / ModEM 3-D / MARE2DEM working directory, or an AI/DL
       ``.npz``/``.npy`` array bundle → :term:`PCSF` / :term:`PCSM`.
   * - PCSF ⇄ PCSM / Validate / Info
     - Lossless transcoding between the two, plus structural validation and
       a full summary of an existing file.
   * - EDI ⇄ EMTF-XML
     - Historical SEG :term:`EDI` transfer functions ⇄ IRIS/EMTF-XML, one
       file or a whole directory, in either direction.
   * - Build PCBH
     - CSV / relational CSV directory / XLSX / LAS 2.0 log → a PCBH borehole
       document.
   * - Build PCGL
     - CSV of named resistivity units → a PCGL geology legend.
   * - Build PCGS
     - CSV(s) of planar/linear structural measurements and fault traces → a
       PCGS structural-evidence document.
   * - Build PCPT
     - CSV/XLSX of named points → a PCPT points-of-interest document.
   * - Batch queue
     - A mixed folder of EDI / EMTF-XML / inversion sources, converted
       sequentially with a per-row status.
   * - Settings
     - Persisted defaults (EPSG, UTM zone, ModEM conventions, EDI/XML
       policy, theme) applied to every new job.

Every page is a thin widget around one function in
:mod:`pycsamt.app.converter.jobs` — plain, Qt-free Python that reuses exactly
the engine the ``pycsamt format`` command line drives
(:mod:`pycsamt.format.convert_engine`, :mod:`pycsamt.emtf`,
:mod:`pycsamt.format.borehole`, :mod:`pycsamt.format.geology`,
:mod:`pycsamt.format.structure`, :mod:`pycsamt.format.pointset`). This page
documents both faces side by side: click through the GUI, or drive the same
conversion from a script or CI job with the CLI command shown in each
section. :doc:`/cli/format` already has the exhaustive option reference for
``convert`` / ``detect`` / ``info`` / ``validate``; this page links there
rather than repeating it, and gives full treatment to the six sub-commands
(``edi-to-xml``, ``xml-to-edi``, ``build-pcbh``, ``build-pcgl``,
``build-pcgs``, ``build-pcpt``) that are not documented anywhere else yet.

Install And Launch
-------------------

The converter app ships with the same optional application extra as the
full desktop suite — one install covers every pyCSAMT GUI surface:

.. code-block:: bash

   pip install "pycsamt[app]"

For development from a source checkout:

.. code-block:: bash

   pip install -e ".[app,dev]"

Launch it with the dedicated entry point:

.. code-block:: bash

   pycsamt-converter

which is equivalent to:

.. code-block:: bash

   python -m pycsamt.app.converter

Either form opens a window titled **pyCSAMT Format Studio** at a modest
920×620 by default (minimum 760×520) — small enough to sit next to a file
browser, not a window that swallows the whole screen. Every page scrolls
internally if its content genuinely needs more room than the window offers,
so resizing smaller never clips a field.

Unlike the full desktop app, the converter app never needs survey data
loaded first — there is no station list, no map, no EDI folder to open at
startup. You pick a source path on whichever page you need and go. Nothing
here writes into your data folders except the output path you explicitly
choose, and every builder/convert action refuses to overwrite an existing
file unless its **Overwrite existing output** checkbox is ticked (CLI:
``--overwrite``).

.. _applications-converter-compile:

Compiling A Standalone Build
------------------------------

The converter app never *needs* a binary — ``pycsamt-converter`` /
``python -m pycsamt.app.converter`` work directly once
``pip install "pycsamt[app]"`` is done. Packaging exists for one specific
case: handing a double-click-able build to someone who does not want a
Python environment at all. It is the *only* pyCSAMT application shipped this
way — the full desktop suite, the web app, Map View, and the agent app are
not packaged as frozen binaries, precisely because this app's narrow scope
(convert files, nothing else) is what keeps the freeze small enough to be
worth building.

Build it with `PyInstaller <https://pyinstaller.org/>`__ from a Python
environment that has ``pycsamt`` installed together with its ``app`` extra,
plus PyInstaller itself:

.. code-block:: bash

   pip install -e ".[app]"
   pip install pyinstaller

Then, from that environment and the repository root, on Windows:

.. code-block:: powershell

   powershell -ExecutionPolicy Bypass -File packaging\pyinstaller\build_converter.ps1

or the equivalent on Linux/macOS:

.. code-block:: bash

   bash packaging/pyinstaller/build_converter.sh

Both scripts wrap ``pyinstaller --noconfirm --clean
packaging/pyinstaller/pycsamt_converter.spec`` — run that directly for more
control over PyInstaller's own flags.

The build is a **onedir** build, not a onefile: the app lands in
``dist/pycsamt-converter/`` as ``pycsamt-converter.exe`` (or
``pycsamt-converter`` on Linux/macOS) plus an ``_internal/`` folder of
dependencies. Onedir starts faster and is easier to debug than onefile's
self-extracting-to-temp approach — the whole folder is what you zip up and
hand to someone, not just the executable.

Expect roughly 600 MB — that is the scientific Python stack (NumPy, SciPy,
pandas, Matplotlib, h5py, pyproj, openpyxl) plus Qt/PySide6, not the
converter app's own code. The spec explicitly excludes PyTorch, TensorFlow,
and the other pyCSAMT apps (``pycsamt.ai``,
``pycsamt.app.desktop``/``web``/``mapview``/``agent_master``): none of those
are reachable from :mod:`pycsamt.app.converter.jobs`, so leaving them out
costs no functionality while keeping the frozen build a few hundred MB
smaller. If a future page needs one of them, drop it from the ``excludes``
list in the spec rather than fighting a ``ModuleNotFoundError`` at runtime.

Before handing out a build, smoke-test it: launch
``dist/pycsamt-converter/pycsamt-converter.exe`` on a clean machine, then
convert one real file on each page — an inversion working directory or
``.npz`` on the Inversion page, a PCSF↔PCSM transcode, a real ``.edi`` on the
EDI↔XML page, a CSV on each of the four builder pages, and a mixed batch
queue. A page that silently does nothing usually means an uncaught
traceback the windowed (``console=False``) build swallowed; rebuild with
``console=True`` in the spec's ``EXE(...)`` call temporarily to see it.
``pytest pycsamt/app/converter/tests/test_jobs.py`` already covers every
conversion function's correctness in-process — the frozen-build smoke test
above only checks that PyInstaller's dependency bundling did not drop
something the tests do not exercise (import resolution, Qt plugin loading,
data files).

.. note::

   The build is not code-signed, so Windows SmartScreen / macOS Gatekeeper
   will warn on first launch — mention this to whoever receives it. First
   launch is also slower than later ones while the OS caches the bundled
   DLLs; this is normal, not a hang. Rebuild after any dependency bump —
   there is no CI job producing this binary automatically yet, so an old
   build silently stops covering a new format/emtf submodule until someone
   rebuilds it.

Finding Your Way Around
-------------------------

Every page shares the same frame, visible in the screenshot above:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Chrome
     - What it does
   * - Menu bar
     - **File** (Settings, Quit) · **Tools** — one checkable entry per page,
       ``Ctrl+1`` .. ``Ctrl+9``, kept in sync with the sidebar selection ·
       **View** (toggle the log dock, Light/Dark theme) · **Help**
       (documentation, GitHub, about).
   * - Sidebar
     - The same nine tools, in the same order as the **Tools** menu; click
       or use the shortcut, either updates the other.
   * - Log dock
     - One shared, dockable panel at the bottom (capped around 140 px) that
       every page reports progress and results to — not a separate log
       embedded in each page. Toggle it from **View → Log**.
   * - Status bar
     - A quiet "Ready ●" indicator; most feedback lives in the log dock
       instead.

Theme (**View → Theme → Light/Dark**) persists across launches through the
same :class:`~pycsamt.app.converter.settings.ConverterSettings` that backs
the Settings page. Primary actions — Convert, Build, Transcode, Run all,
Save settings — are always the single orange button on a page; anything
secondary (Validate, Info, Browse…) stays neutral, so the one action that
actually writes a file is never ambiguous.

.. _applications-converter-source-detection:

Source Detection
-------------------

Two pages — Inversion and the batch queue — rely on the same detector,
:func:`pycsamt.format.detect_source`. It never
loads heavy array data; it looks at file names, extensions, and a handful of
solver signatures (an ``*.iter`` file for Occam2D, ``Modular_NLCG.log`` /
``ModEM.inv`` for ModEM, a ``*.poly`` PSLG for MARE2DEM, resistivity/key
names for an AI ``.npz``) and reports what it found. :doc:`/cli/format`
documents every signature and the full category table
(``pcsf`` / ``pcsm`` / ``solver`` / ``ai_arrays`` / ``unknown``); on the
Inversion page, the same result appears live under the source field as soon
as you drop a path in, for example:

.. code-block:: text

   Detected: solver (modem) — Detected modem working directory.

When a directory carries more than one solver's signature, detection
reports medium confidence and picks one — use **Force solver** in Advanced
options (CLI: ``--solver``) to override it rather than renaming files.

Inversion → PCSF/PCSM
------------------------

This page takes whatever an inversion produced and turns it into one
self-describing :term:`PCSF` (or its ASCII sibling, :term:`PCSM`) file that
:doc:`Map View </applications/mapview/index>`, the desktop's 3-D panel, the
web map view, or a plain analysis script can open without knowing which
solver produced it. See :doc:`/user_guide/models/pcsf_format` for the format
itself — canonical resistivity is always linear :math:`\Omega\,\mathrm{m}`,
with the source backend's native encoding preserved alongside it as
provenance, never assumed.

**Source** accepts a solver working directory (or a single signature file
inside one, such as ``run/ITER17.iter`` or ``demo.poly`` — the detector
resolves it to that file's folder), a ``.npz``/``.npy`` AI array bundle, or
an existing ``.pcsf``/``.pcsm`` file, and shows the live detection result
described above. **Output folder** and **Output format** (``pcsf`` /
``pcsm`` / ``pcsm.gz``) control where the result lands.

**Advanced options** lays its short scalar fields out two per row —
nine fields become five grid rows instead of nine — with two full-width
DropZones and two checkboxes below:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Field
     - Meaning
   * - Force solver
     - Override auto-detection (``occam2d`` / ``modem`` / ``mare2dem``).
       Needed only when a directory carries more than one solver's
       signature.
   * - Occam2D iteration
     - Which ``.iter`` to convert; the special value **final (default)**
       (spin value ``0``) picks the last one automatically.
   * - EPSG / UTM zone
     - Projection used together with **Topography source** and, for AI
       ``grid2d``/``grid3d`` sources, to place lon/lat station coordinates.
   * - AI array encoding
     - ``auto`` / ``linear`` / ``log10`` / ``ln`` — overrides the
       key-name-based encoding guess for an AI ``.npz`` (a ``log10_rho`` or
       ``ln_rho`` key is detected automatically; this forces it when the
       array key is ambiguous, e.g. plain ``resistivity``).
   * - ModEM station-Z
     - How to read the ``.dat`` file's station Z column. ``auto`` flips a
       positive-down depth to an elevation only when the file's own
       comments say so; ``depth_down`` forces that flip; ``elevation``
       trusts the column as-is. Get this wrong and every station plots at
       the wrong depth relative to topography.
   * - ModEM air threshold (Ω·m)
     - Cells above this resistivity are masked as above-topography air fill
       rather than real model resistivity (default :math:`10^{8}`); pass
       ``0`` to disable masking entirely. This exists because ModEM's air
       layers are typically filled with an enormous but finite resistivity
       (not NaN), which would otherwise dominate any colour scale a viewer
       tries to build from the array.
   * - Created by / Description
     - Free-text values written into the PCSF ``created_by`` /
       ``description`` attributes.

Below the grid, **Topography source** accepts a ``.bln``/``.csv``/``.stn``
file or an EDI directory and is attached per station (Occam2D and ModEM
only) via :mod:`pycsamt.format.topo_source`'s name/positional matching;
**MARE2DEM .poly** lets you rebuild the mesh from a specific PSLG file
instead of the one auto-discovered next to the run (MARE2DEM conversion
needs the ``triangle`` package — see :ref:`applications-converter-troubleshooting`
if it is missing). **Write PCSM block in log10(rho)** only affects the
human-readable text block for ``--to pcsm``; the canonical stored
resistivity stays linear either way.

Under the hood this is one call to :mod:`pycsamt.format.convert_engine`'s
adapters, exactly what ``pycsamt format convert`` runs:

.. code-block:: console

   pycsamt format convert data/modem/run01/ run01.pcsm.gz \
       --topo data/AMT/WILLY_DATA --epsg 32648

See :doc:`/cli/format` for the complete ``convert`` option table (every
field above has a matching flag) and more worked examples across all three
solvers and AI arrays.

PCSF ⇄ PCSM / Validate / Info
--------------------------------

One source field, three actions. **Transcode** (the orange button)
round-trips a ``.pcsf``/``.pcsm``/``.pcsm.gz`` file to the other encoding —
**Transcode to** and **Output file** only matter for this action.
**Validate** and **Info** read the same source field and ignore the output
field entirely; they exist to answer "is this file sound?" and "what is
actually in this file?" without writing anything.

.. figure:: /_static/applications/converter/converter-pcsf-pcsm-validate.png
   :alt: PCSF/PCSM page after running Validate on a real ModEM-derived file
   :class: pycsamt-screenshot
   :align: center
   :width: 92%

   **Validate** run against the same ``27-freque-watex-data-04.pcsf`` from
   the Inversion example above. Four checks run in order — header, load,
   schema, round-trip — and every one has to pass for the file to report
   **VALID**; the round-trip step specifically re-writes the model to a
   temporary file of the same encoding, reads it back, and compares the
   resistivity array element-for-element, which catches a writer silently
   dropping precision that a header-only check would miss.

``Validate``'s four steps and ``Info``'s full field-by-field summary
(geometry, every stored resistivity-like array with its NaN count, station
table, topography kind, provenance) are documented in full, with more
examples, under :doc:`/cli/format`:

.. code-block:: console

   pycsamt format info model.pcsf
   pycsamt format validate model.pcsm --no-roundtrip

EDI ⇄ EMTF-XML
-----------------

Losslessly converts between historical SEG :term:`EDI` transfer functions
and IRIS/EMTF-XML, one file or an entire directory at a time. See
:doc:`/user_guide/emtf/edi_interop` for exactly what maps to what; the short
version is that both formats store the same impedance/tipper tensors and
uncertainties, but EDI is a fixed-width ASCII block format from the 1980s
and EMTF-XML is a modern, richly-annotated schema, so some EMTF-only
metadata (arbitrary transfer-function types, full covariance matrices,
several provenance blocks) has no EDI destination — the app never silently
drops it without telling you (see **EDI/XML data-loss policy** below).

**Direction** switches the whole form between **EDI → EMTF-XML** and
**EMTF-XML → EDI**; the options below it change to match:

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - Option
     - Direction
     - Meaning
   * - Prefer EDI SPECTRA blocks
     - EDI → XML
     - When an EDI carries raw ``SPECTRA`` cross-power blocks *and*
       impedance/tipper, prefer re-deriving the transfer function from
       SPECTRA (closer to the original field measurement) over trusting the
       impedance already written in the file.
   * - EDI/XML data-loss policy
     - XML → EDI
     - ``warn`` (default) logs what could not be represented and continues;
       ``raise`` aborts the conversion instead; ``ignore`` converts
       silently. This is exactly the gap described above — an EMTF document
       with no lossless EDI equivalent for some of its content.
   * - Strict XML parsing
     - XML → EDI
     - Reject a malformed ``<Period>`` (or similar) outright versus drop
       just that element and continue with a shorter document.

.. figure:: /_static/applications/converter/converter-edi-xml-convert.png
   :alt: EDI to EMTF-XML page after converting a folder of 25 EDI files
   :class: pycsamt-screenshot
   :align: center
   :width: 92%

   A real folder conversion in progress: 25 EDI files converted one at a
   time, each logged as it completes, ending in a one-line summary. The same
   per-item progress reporting backs the batch queue below.

Both directions accept a single file or a whole folder for **Source**, and
neither is documented anywhere else in the CLI reference yet, so here is
the full picture. EDI → XML:

.. code-block:: console

   pycsamt format edi-to-xml SOURCE [--prefer-spectra/--no-prefer-spectra] \
       [-o OUTPUT_DIR] [--overwrite] [-v]

Run against the bundled Broken Hill sample survey, it converts every
``.edi`` in the folder and reports each one as it finishes:

.. code-block:: console

   pycsamt format edi-to-xml data/MT/broken-hill/edis -o xml_out

which prints:

.. code-block:: text

                             format edi-to-xml
   ┌─────────────────────────────────┬─────────────────────────────────┐
   │ data/MT/broken-hill/edis/BH_10_… │ xml_out/BH_10_imp_rev.xml       │
   │ data/MT/broken-hill/edis/BH_11_… │ xml_out/BH_11_imp.xml           │
   │                              …   │                              …  │
   └─────────────────────────────────┴─────────────────────────────────┘

   ✓ 21 file(s) written to xml_out

XML → EDI takes the same shape, with the loss-policy and strictness flags
instead:

.. code-block:: console

   pycsamt format xml-to-edi SOURCE [--on-loss {warn,raise,ignore}] \
       [--strict/--permissive] [-o OUTPUT_DIR] [--overwrite] [-v]

Feeding the XML directory produced above straight back through this
direction round-trips it to EDI, aborting instead of warning on any
irrecoverable EMTF content this time:

.. code-block:: console

   pycsamt format xml-to-edi xml_out -o edi_out --on-loss raise

.. dropdown:: Data and model citation
   :animate: fade-in
   :color: secondary

   The Broken Hill EDI files used in the examples above are the real,
   bundled sample survey under ``data/MT/broken-hill/edis`` — cite both the
   article and the data release if you reuse them:

   AlQahtani, Y., Ozaydin, S., Chatzaras, V., Rey, P. F., & Passos, T.
   (2026). Why does the Broken Hill deposit sit in resistive crust?
   Magnetotelluric evidence for metamorphic decoupling of a world-class
   mineral system. *Journal of Geophysical Research: Solid Earth*,
   **131**, e2026JB035666. `<https://doi.org/10.1029/2026JB035666>`__

   AlQahtani, Y., Ozaydin, S., Chatzaras, V., Rey, P. F., & Passos, T.
   (2026). Open data for the article "Why does the Broken Hill deposit
   sit in resistive crust?". *Zenodo*.
   `<https://doi.org/10.5281/zenodo.21272106>`__

.. _applications-converter-pcbh:

Build PCBH (Borehole)
------------------------

Builds a :doc:`PCBH </user_guide/geology/pcbh_format>` document — one or more
boreholes with collar, interval, and lithology data — from a combined CSV, a
directory of relational tables (with an ``import.yaml`` manifest), an XLSX
workbook, or a single LAS 2.0 log. **Source kind** defaults to ``auto``
(inferred from the path: a directory is ``csv-dir``, ``.xlsx``/``.xlsm`` is
``xlsx``, ``.las`` is ``las``, anything else is ``csv``); force it
explicitly when a source needs a specific reader.

.. figure:: /_static/applications/converter/converter-pcbh-build.png
   :alt: PCBH page after building a two-borehole document from CSV
   :class: pycsamt-screenshot
   :align: center
   :width: 92%

   A combined collar+interval CSV — two boreholes, two lithology intervals
   each — built successfully: 4 rows read, 4 accepted, 0 rejected, but
   **11 issues** logged anyway. Issues are non-fatal warnings (missing
   optional fields, unusual but not invalid values); the document still
   writes. Check the log for what they were before treating the output as
   final — a large issue count on a small file is worth reading.

The **LAS collar** panel only appears when **Source kind** is (or resolves
to) ``las`` — a single wireline log has no collar coordinates of its own, so
**Collar id**, **X**, **Y**, **Z**, and **Horizontal CRS** are required and
the build refuses to run without all five.

.. important::

   Borehole ``kind`` is a controlled vocabulary, not free text: ``water``,
   ``mining_exploration``, ``mining_production``, ``geotechnical``,
   ``environmental``, ``petroleum``, ``geothermal``, ``scientific``,
   ``monitoring``, or ``unknown``. A CSV value outside that list — a plain
   ``"exploration"`` instead of ``"mining_exploration"``, for example — is
   rejected per-row with ``unsupported borehole kind '...'`` rather than
   silently accepted. ``status`` is similarly controlled: ``planned``,
   ``drilling``, ``completed``, ``suspended``, ``abandoned``,
   ``decommissioned``, ``unknown``.

The equivalent CLI command:

.. code-block:: console

   pycsamt format build-pcbh SOURCE [-o OUT.pcbh.json] [--from {auto,csv,csv-dir,xlsx,las}] \
       [--collar-id ID --x X --y Y --z Z --crs CRS] [--document-id ID] \
       [--created-by TEXT] [--overwrite]

Against the same combined CSV as the screenshot above (``--from`` inferred
as ``csv`` automatically):

.. code-block:: console

   pycsamt format build-pcbh boreholes.csv

which prints the same read/accepted/rejected and issue counts the GUI's log
dock shows:

.. code-block:: text

                     format build-pcbh
   ┌──────────────────────────────┬──────────────────────┐
   │ wrote                        │ boreholes.pcbh.json  │
   │ document_id                  │ csv:boreholes         │
   │ boreholes                    │ 2                      │
   │ rows read/accepted/rejected  │ 4/4/0                  │
   │ issues                       │ 11                     │
   └──────────────────────────────┴──────────────────────┘

   ✓ boreholes.pcbh.json

A LAS log needs the collar flags on the command line, since a LAS file has
no location of its own:

.. code-block:: console

   pycsamt format build-pcbh ZK01.las --collar-id ZK01 \
       --x 512340 --y 2894210 --z 118.4 --crs EPSG:32650

Build PCGL (Geology Legend)
------------------------------

A PCGL document is a small, named lookup table: geological unit names mapped
to a resistivity range, so a rendered depth slice or fence panel can be
coloured and labelled by rock type instead of by raw ohm-metres alone. It is
the same schema :meth:`~pycsamt.geology.lithology.RockDatabase.from_csv`
accepts, wrapped in PCGL's document envelope.

.. figure:: /_static/applications/converter/converter-pcgl-build.png
   :alt: PCGL page after building a four-unit geology legend from CSV
   :class: pycsamt-screenshot
   :align: center
   :width: 92%

   Four resistivity units built into one legend from a plain CSV — the page
   itself states the accepted columns right above the source field, so
   there is no need to look them up elsewhere before starting.

**Source** columns: ``name``, ``rho_min``, ``rho_max`` required;
``color``, ``description``, ``code``, ``source``, ``pattern_id``,
``pattern_source`` optional. **Title** and **Document id** are free text
(document id defaults to a stable value derived from the source filename).

The CLI equivalent:

.. code-block:: console

   pycsamt format build-pcgl SOURCE.csv [-o OUT.pcgl.json] [--title TEXT] \
       [--document-id ID] [--created-by TEXT] [--overwrite]

Against the same ``units.csv`` as the screenshot above:

.. code-block:: console

   pycsamt format build-pcgl units.csv --title "Site A geology legend"

which prints:

.. code-block:: text

              format build-pcgl
   ┌─────────────┬────────────────────────┐
   │ wrote       │ units.pcgl.json         │
   │ document_id │ pcgl:units              │
   │ title       │ Site A geology legend   │
   │ entries     │ 4                       │
   └─────────────┴────────────────────────┘

   ✓ units.pcgl.json

.. _applications-converter-pcgs:

Build PCGS (Structure)
-------------------------

A PCGS document holds structural field evidence — planar measurements
(bedding, foliation), linear measurements (fold axes, lineations), and fault
traces — each in its own optional CSV, combined into one document. At least
one of the three is required; any combination of the other two is fine.

.. figure:: /_static/applications/converter/converter-pcgs-build.png
   :alt: PCGS page after building a structural-evidence document from three CSVs
   :class: pycsamt-screenshot
   :align: center
   :width: 92%

   All three optional CSVs supplied at once — 3 planar, 2 linear, and 2
   fault measurements combined into one PCGS document. Leaving any of the
   three DropZones empty is fine; the hint text above them spells out each
   CSV's required columns so you never have to leave the page to check.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - CSV
     - Required columns
   * - Planar
     - ``x``, ``kind``, ``strike_deg``, ``dip_deg``, ``dip_direction_deg``
   * - Linear
     - ``x``, ``kind``, ``trend_deg``, ``plunge_deg``
   * - Faults
     - ``x``, ``dip_deg``, ``downthrown_side``

.. warning::

   ``downthrown_side`` accepts exactly ``left`` or ``right`` — compass
   directions such as ``east``/``west`` are rejected outright
   (``downthrown_side must be one of ('left', 'right'), got '...'``).
   "Left"/"right" here means relative to the profile direction implied by
   increasing ``x``, not a geographic compass bearing; translate a mapped
   compass direction to left/right before building the CSV.

The CLI equivalent takes the same three optional paths:

.. code-block:: console

   pycsamt format build-pcgs [--planar PATH] [--linear PATH] [--faults PATH] \
       -o OUT.pcgs.json [--title TEXT] [--document-id ID] [--created-by TEXT] [--overwrite]

Against the same three CSVs as the screenshot above:

.. code-block:: console

   pycsamt format build-pcgs --planar planar.csv --linear linear.csv --faults faults.csv \
       -o structure.pcgs.json --title "Site A structural evidence"

which prints:

.. code-block:: text

                     format build-pcgs
   ┌─────────────┬─────────────────────────────────┐
   │ wrote       │ structure.pcgs.json               │
   │ document_id │ pcgs:2026-09-15T20:13:16.233925Z │
   │ planar      │ 3                                  │
   │ linear      │ 2                                  │
   │ faults      │ 2                                  │
   └─────────────┴─────────────────────────────────┘

   ✓ structure.pcgs.json

The ``document_id`` shown above is the generated default — a UTC
timestamp — used whenever ``--document-id`` is omitted; pass an explicit id
for anything you intend to reference elsewhere later.

Build PCPT (Points)
----------------------

A PCPT document is a flat list of named points — drill targets, geophysical
anomalies, sample locations — each with coordinates and an optional note, for
overlay on a map or 3-D scene alongside a resistivity model.

.. figure:: /_static/applications/converter/converter-pcpt-build.png
   :alt: PCPT page after building a three-point target list from CSV
   :class: pycsamt-screenshot
   :align: center
   :width: 92%

   Three targets built from a CSV with an explicit **CRS**. When the source
   has no CRS column and none is given, the builder infers ``EPSG:4326`` if
   any point carries longitude/latitude, or a local unprojected grid
   otherwise — passing **CRS** explicitly, as here, avoids relying on that
   inference for anything that matters.

**Source** columns: ``name``, ``x``, ``y``\ [, ``z``] at minimum (CSV or
XLSX); **Sheet** and a header-row override apply to XLSX only. **CRS** is
free text (``EPSG:4326``, ``EPSG:32650``, …).

The CLI equivalent:

.. code-block:: console

   pycsamt format build-pcpt SOURCE [-o OUT.pcpt.json] [--sheet NAME_OR_INDEX] \
       [--header-row N] [--crs CRS] [--document-id ID] [--overwrite]

Against the same ``targets.csv`` as the screenshot above:

.. code-block:: console

   pycsamt format build-pcpt targets.csv --crs EPSG:32650

which prints:

.. code-block:: text

            format build-pcpt
   ┌─────────────┬────────────────────┐
   │ wrote       │ targets.pcpt.json  │
   │ document_id │ pcpt:targets       │
   │ points      │ 3                  │
   │ crs         │ EPSG:32650         │
   └─────────────┴────────────────────┘

   ✓ targets.pcpt.json

Batch Queue
--------------

Converts a mixed pile of files — inversion results/AI arrays, ``.pcsf``,
``.pcsm``, ``.edi``, EMTF-XML — one after another with a per-row status,
instead of repeating the Inversion or EDI ⇄ XML page once per file.

.. figure:: /_static/applications/converter/converter-batch-queue-run.png
   :alt: Batch queue after converting four Broken Hill EDI files to EMTF-XML
   :class: pycsamt-screenshot
   :align: center
   :width: 92%

   Four real EDI files added via **Add folder…**, auto-classified as
   ``edi`` in the **Kind** column, run, and marked **done** one at a time —
   the log dock mirrors the same per-item progress shown in the table, and
   the closing summary line ("4 done, 0 failed, 0 skipped") is the fastest
   way to confirm a large batch actually finished cleanly.

**Add files…**/**Add folder…** classify each path automatically via
:func:`~pycsamt.app.converter.jobs.classify_batch_item` — ``.edi`` →
``edi``, ``.xml`` → ``xml``, anything :func:`~pycsamt.format.detect_source`
calls ``solver``/``ai_arrays``/``pcsf``/``pcsm`` → ``inversion``, everything
else → ``unknown`` (queued but always skipped). Each row converts
independently: one failure is recorded in **Status** and does not stop the
rest of the queue, which matters for a folder where most files are good but
a handful are corrupt or from an unrelated survey.

**Inversion target format** (``pcsf``/``pcsm``) applies only to
``inversion``-kind rows; an ``edi`` row always converts to EMTF-XML and an
``xml`` row always converts to EDI — there is no direction choice per row,
unlike the dedicated EDI ⇄ XML page. **PCBH/PCGL/PCGS/PCPT builders are not
queued here** — each needs its own structured per-file options (collar
coordinates, planar/linear/fault CSV groupings) that do not fit one generic
row; use their own dedicated pages, one file at a time, instead.

There is no single CLI command equivalent to the batch queue — script the
same behaviour with a shell loop over the per-format commands documented
above, for example:

.. code-block:: bash

   for f in survey_edi/*.edi; do
       pycsamt format edi-to-xml "$f" -o survey_xml/
   done

Settings
-----------

Persisted defaults applied to every *new* job on every page — not a global
override you have to fight per-run. Any value set explicitly on a page
(EPSG, UTM zone, a checkbox) still wins over the default for that one job;
Settings only fills in what you leave blank.

.. figure:: /_static/applications/converter/converter-settings-saved.png
   :alt: Settings page after saving persisted converter defaults
   :class: pycsamt-screenshot
   :align: center
   :width: 92%

   Defaults set once — an EPSG/UTM zone pair, the ModEM conventions, and the
   EDI/XML policy — apply to the Inversion and EDI ⇄ XML pages from the next
   job onward without retyping them.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Setting
     - Applies to
   * - Default EPSG / Default UTM zone
     - Inversion page's ``--topo``/lon-lat projection, when left blank there.
   * - ModEM station-Z convention / air threshold
     - Inversion page's ModEM-specific fields, as their own defaults.
   * - EDI/XML data-loss policy
     - EDI ⇄ EMTF-XML page's XML → EDI direction.
   * - Prefer EDI SPECTRA blocks / Strict EMTF-XML parsing
     - EDI ⇄ EMTF-XML page's two direction-specific checkboxes.
   * - Write PCSM blocks in log10(rho) by default
     - Inversion and PCSF ⇄ PCSM pages' log10-view checkboxes.
   * - Open output folder after each conversion
     - Every builder/convert action, once a job finishes.
   * - Overwrite existing output without asking
     - Every page's overwrite protection — the app-wide equivalent of
       ticking **Overwrite existing output** on every single page.

Settings persist through Qt's ``QSettings`` under the
``earthai-tech`` / ``pycsamt-converter`` org/app pair — an INI file under
the OS user-config directory on Linux, or the native registry on Windows,
resolved automatically by :mod:`pycsamt.app.converter.settings`. There is no
CLI equivalent: a script instead passes the same options as flags on every
invocation, or wraps them in a shell alias/wrapper for repeated use.

.. _applications-converter-troubleshooting:

Troubleshooting
-------------------

``… exists — pass --overwrite to replace it`` / a page refuses to run
   The target file already exists and **Overwrite existing output** is
   unticked (CLI: pass ``--overwrite``). This is deliberate — nothing in
   the app replaces a file you did not explicitly allow it to.

``Nothing to convert — … looks like: …``
   :func:`~pycsamt.format.detect_source` returned ``unknown`` for the
   Inversion page's source. Check the live "Detected: …" hint under the
   source field, and force **Force solver** if it is a non-standard solver
   directory — see :ref:`applications-converter-source-detection`.

``MARE2DEM → PCSF needs the 'triangle' package``
   MARE2DEM conversion rebuilds the run's triangulation from its ``.poly``
   PSLG and needs the optional ``triangle`` package:
   ``pip install triangle``.

``unsupported borehole kind '...'`` / ``downthrown_side must be one of (...)``
   A CSV value falls outside PCBH's or PCGS's controlled vocabulary — see
   the tables under :ref:`Build PCBH <applications-converter-pcbh>` and
   :ref:`Build PCGS <applications-converter-pcgs>` above for the exact
   accepted values.

``--from las needs --collar-id, --x, --y, --z and --crs``
   A LAS 2.0 log has no collar location of its own; all five must be
   supplied (GUI: the **LAS collar** panel; CLI: the matching flags).

"Pick a source file or folder first." / "Pick an output … first." in the log
   A required field is still empty — every page validates its own required
   fields locally before starting a background job, rather than letting a
   background worker fail on missing input.

A page seems to hang
   Every conversion runs on a background ``QThread``
   (:class:`~pycsamt.app.converter.workers.ConversionWorker`); the log
   dock's progress bar switches to indeterminate mode for jobs with no
   per-item progress (most builders) and shows real ``[n/total]`` steps for
   ones that do (EDI ⇄ XML, batch queue). A large EDI/XML folder or a big
   ModEM grid genuinely takes time — watch the log dock rather than the
   window title for confirmation it is still working.

Frozen-build-specific issues (blank window, missing DLL, slow first launch) are
covered separately under :ref:`applications-converter-compile`.

.. seealso::

   :doc:`/cli/format`
       The complete ``pycsamt format`` command reference — full option
       tables for ``convert``/``detect``/``info``/``validate``, more worked
       examples across every solver, and the Python-API equivalents this
       whole app (and CLI) sits on top of.

   :doc:`/user_guide/models/pcsf_format`
       What PCSF/PCSM actually are, the canonical-resistivity invariant,
       and how every solver's native encoding maps into it.

   :doc:`/user_guide/geology/pcbh_format`
       The full PCBH schema — collars, trajectories, intervals, structures,
       samples — behind the Build PCBH page.

   :doc:`/user_guide/emtf/edi_interop`
       Exactly what does and does not survive an EDI ⇄ EMTF-XML round trip.

   :doc:`/applications/desktop/index`
       The full survey-processing suite this app deliberately leaves out —
       QC, corrections, forward modelling, inversion, interpretation.
