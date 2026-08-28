Format Commands
===============

``pycsamt format`` is the command group for the **PCSF / PCSM**
inversion-result format. PCSF (``.pcsf``, an HDF5 container) is
pyCSAMT's backend-neutral representation of a finished resistivity
model; PCSM (``.pcsm``, or gzip-compressed ``.pcsm.gz``) is its
lossless, hand-editable ASCII sibling.

Every classical solver (Occam2D, ModEM 3-D, MARE2DEM) and any AI/DL
inversion converts *to* PCSF. The group has one job: take whatever an
inversion produced — a solver working directory, a ``.npz`` of
predicted arrays, or an existing ``.pcsf`` / ``.pcsm`` — and turn it
into a single self-describing file that every downstream viewer
(``pycsamt map``, the desktop 3-D panel, the web map view, a plain
analysis script) can open without knowing which backend produced it.

Command Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Command
     - Purpose
     - Main output
   * - ``pycsamt format convert``
     - Convert any inversion result to ``.pcsf`` or ``.pcsm``.
     - One PCSF/PCSM file plus a written-file summary.
   * - ``pycsamt format detect``
     - Report what a file or folder is, and how ``convert`` would
       treat it.
     - A classification table (or JSON).
   * - ``pycsamt format info``
     - Full summary of an existing ``.pcsf`` / ``.pcsm`` file.
     - Geometry, resistivity statistics, stations, topography,
       provenance, history.
   * - ``pycsamt format validate``
     - Structural + round-trip check of a ``.pcsf`` / ``.pcsm`` file.
     - Per-step pass/fail; non-zero exit on any failure.

Source Detection
----------------

``convert`` and ``detect`` share one detector
(:func:`pycsamt.format.detect_source`). It never loads heavy array
data — it looks at file names, extensions, and a handful of solver
signatures, then reports a ``category``:

.. list-table::
   :header-rows: 1
   :widths: 18 30 52

   * - ``category``
     - Recognised from
     - Converted through
   * - ``pcsf``
     - a ``.pcsf`` file
     - direct transcode to ``.pcsm``
   * - ``pcsm``
     - a ``.pcsm`` / ``.pcsm.gz`` file
     - direct transcode to ``.pcsf``
   * - ``solver``
     - an Occam2D / ModEM / MARE2DEM working directory (or a
       signature file inside one)
     - :mod:`pycsamt.format.adapters` (``occam2d_to_pcsf``,
       ``modem3d_to_pcsf``, ``mare2dem_to_pcsf``)
   * - ``ai_arrays``
     - a ``.npz`` / ``.npy`` array bundle
     - :mod:`pycsamt.format.adapters.generic` (``grid2d_to_pcsf`` /
       ``grid3d_to_pcsf`` / ``mesh_to_pcsf``)
   * - ``unknown``
     - anything else
     - nothing — ``convert`` refuses

Solver directory signatures:

* **Occam2D** — any ``*.iter`` file, or names such as ``Occam2DMesh``,
  ``Occam2DModel``, ``OccamStartup``, ``OccamDataFile.dat``.
* **ModEM** — ``Modular_NLCG.log``, ``ModEM.inv``, or a ``*.rho`` model
  file next to a ``*.dat`` data file.
* **MARE2DEM** — a ``*.poly`` PSLG next to a ``*.resistivity`` file, or
  ``mare2dem.settings``.

When a directory carries more than one signature, the detector reports
``medium`` confidence and picks one; pass ``--solver`` to force the
choice. Pointing ``convert`` at a single file inside a run (for
example ``run/ITER17.iter`` or ``demo.poly``) resolves to that file's
folder automatically.

AI array bundles
~~~~~~~~~~~~~~~~~

For a ``.npz``, the detector matches keys case-insensitively:

* **resistivity** — ``resistivity``, ``rho``, ``model``, ``prediction``,
  ``pred``, ``output`` … (linear ohm-m). A ``log10_rho`` / ``ln_rho``
  key is also recognised and sets the encoding automatically.
* **grid2d** — coordinate arrays ``x`` and ``z`` (plus optional
  ``x_nodes`` / ``z_nodes``), or a 2-D resistivity array.
* **grid3d** — ``x``, ``y``, ``z``, or a 3-D resistivity array.
* **mesh** — ``nodes`` (or ``vertices`` / ``points``) together with
  ``connectivity`` (or ``triangles`` / ``elements``), plus optional
  ``region_ids``.
* optional ``uncertainty`` and ``sensitivity`` arrays are carried
  through when present.

A bare ``.npy`` is accepted as linear resistivity on a unit-spaced
grid (``low`` confidence) — prefer an ``.npz`` with real coordinates.

Detect
------

Usage:

.. code-block:: console

   pycsamt format detect PATH [--solver {occam2d,modem,mare2dem}] [-f {text,json}]

``detect`` performs no conversion and writes nothing. It exits non-zero
when the source is not convertible, which makes it a cheap pre-flight
check in scripts.

.. code-block:: console

   pycsamt format detect data/occam2D
   pycsamt format detect data/mare2dem/demo_mt_inversion
   pycsamt format detect unet_prediction.npz -f json

Convert
-------

Usage:

.. code-block:: console

   pycsamt format convert SOURCE [TARGET] [OPTIONS]

If ``TARGET`` is omitted, the output name is derived from the source
and written into ``--output-dir`` (default: the current directory).
The output format comes from ``TARGET``'s extension, or from ``--to``
(default: ``pcsf``).

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Option
     - Default
     - Meaning
   * - ``--to {pcsf,pcsm,pcsm.gz}``
     - ``pcsf``
     - Output format when ``TARGET`` has no telling extension.
   * - ``--solver {occam2d,modem,mare2dem}``
     - auto
     - Force the source backend instead of fingerprinting it.
   * - ``--iteration INT``
     - final
     - Occam2D iteration index to convert.
   * - ``--topo PATH``
     - none
     - Topography source (``.bln`` / ``.csv`` / ``.stn`` / EDI
       directory) attached per station. Occam2D and ModEM only.
   * - ``--epsg INT`` / ``--utm-zone ZONE``
     - none
     - Projection for ``--topo`` and station lon/lat.
   * - ``--encoding {linear,log10,ln}``
     - auto
     - Encoding of the AI resistivity array. Overrides key-name
       detection.
   * - ``--poly PATH``
     - auto
     - MARE2DEM ``.poly`` PSLG to rebuild the mesh from.
   * - ``--origin X,Y[,Z]``
     - none
     - Real-world grid origin for AI ``grid2d`` / ``grid3d`` sources.
   * - ``--azimuth DEG``
     - none
     - Profile azimuth (``grid2d``) or volume rotation (``grid3d``).
   * - ``--created-by TEXT`` / ``--description TEXT``
     - CLI defaults
     - Values written into the PCSF attributes.
   * - ``--log10-view``
     - false
     - For ``--to pcsm``: render the human-readable block in
       ``log10(rho)``. The canonical stored resistivity stays linear.
   * - ``--dry-run``
     - false
     - Detect and print the plan without writing anything.
   * - ``--overwrite``
     - false
     - Replace an existing output file.
   * - ``-f {text,json}``
     - ``text``
     - Console output format for the result summary.

Canonical resistivity in a PCSF file is **always linear ohm-m**,
regardless of the source. A backend's native encoding
(``log10`` for Occam2D, ``ln`` for ModEM, a ``--encoding`` for an AI
array) is preserved separately as provenance.

Examples
~~~~~~~~

Occam2D working directory to PCSF:

.. code-block:: console

   pycsamt format convert data/occam2D occam.pcsf

ModEM 3-D directory to a compressed PCSM, with topography from the EDI
survey:

.. code-block:: console

   pycsamt format convert data/modem/run01/ run01.pcsm.gz \
       --topo data/AMT/WILLY_DATA --epsg 32648

MARE2DEM run (the mesh is rebuilt from the run's ``.poly`` PSLG — this
needs the ``triangle`` package):

.. code-block:: console

   pycsamt format convert data/mare2dem/demo_mt_inversion mare.pcsf

An AI/DL prediction saved as ``resistivity`` / ``x`` / ``z`` arrays:

.. code-block:: console

   pycsamt format convert unet_prediction.npz unet.pcsf --encoding log10

Transcode between the two encodings (lossless both ways):

.. code-block:: console

   pycsamt format convert occam.pcsf occam.pcsm --log10-view
   pycsamt format convert occam.pcsm roundtrip.pcsf

Preview without writing:

.. code-block:: console

   pycsamt format convert data/mare2dem/demo_mt_inversion --dry-run

Info
----

Usage:

.. code-block:: console

   pycsamt format info FILE [-f {text,json}]

``info`` loads the file and prints geometry (axis extents, mesh sizes,
origin/rotation), resistivity statistics for every stored array
(``resistivity``, ``resistivity_native``, ``resistivity_by_region``,
``resistivity_by_node``, ``uncertainty``, ``sensitivity``), the station
table, topography kind, inversion history keys, and — for AI results —
the embedded ``model_provenance`` block.

.. code-block:: console

   pycsamt format info model.pcsf
   pycsamt format info model.pcsm.gz -f json

Validate
--------

Usage:

.. code-block:: console

   pycsamt format validate FILE [--no-roundtrip] [-f {text,json}]

``validate`` runs four checks and exits non-zero if any fail:

1. **header** — the geometry kind can be peeked;
2. **load** — the file parses into a :class:`~pycsamt.format.PCSFModel`;
3. **schema** — ``PCSFModel.validate()`` passes;
4. **roundtrip** — re-writing the model to a temp file of the same
   encoding and reading it back reproduces the resistivity array
   (skip with ``--no-roundtrip``).

.. code-block:: console

   pycsamt format validate model.pcsf
   pycsamt format validate model.pcsm --no-roundtrip -f json

Common Failures
---------------

``Nothing to convert — … looks like: …``
    The detector returned ``unknown``. Run ``pycsamt format detect`` on
    the same path to see why, and pass ``--solver`` if it is a
    non-standard solver directory.

``No .poly PSLG found``
    A MARE2DEM directory has no polygon mesh file. Pass ``--poly PATH``.

``MARE2DEM → PCSF needs the 'triangle' package``
    Install it with ``pip install triangle``. It rebuilds the run's
    triangulation from the ``.poly`` PSLG.

``… exists — pass --overwrite to replace it``
    The target file is already there.

``… is already PCSF; give a different TARGET or --to the other format``
    ``convert`` will not transcode a file onto itself.

Python Equivalents
------------------

The CLI is a thin layer over :mod:`pycsamt.format`:

.. code-block:: python

   from pycsamt.format import detect_source, write_pcsf
   from pycsamt.format.adapters import occam2d_to_pcsf
   from pycsamt.models.occam2d.results import InversionResult

   sk = detect_source("data/occam2D")           # -> SourceKind(category="solver", ...)
   result = InversionResult(workdir=sk.path)
   model = occam2d_to_pcsf(result, created_by="me")
   write_pcsf(model, "occam.pcsf")

For an AI result, skip the solver classes entirely:

.. code-block:: python

   import numpy as np
   from pycsamt.format import write_pcsm
   from pycsamt.format.adapters.generic import grid2d_to_pcsf

   d = np.load("unet_prediction.npz")
   model = grid2d_to_pcsf(d["resistivity"], d["x"], d["z"], encoding="log10")
   write_pcsm(model, "unet.pcsm")
