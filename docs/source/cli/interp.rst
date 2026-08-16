Interpretation Commands
=======================

``pycsamt interp`` turns electrical resistivity models into geological
interpretation products. The command group has three jobs:

* inspect or query the built-in rock-physics database;
* classify an inversion result into pseudo-stratigraphic station logs;
* export those interpreted products to common geoscience formats.

The current CLI implementation supports Occam2D work directories for
``classify`` and ``export``. The lower-level interpretation API can be
extended to other model containers, but this CLI page documents only the
commands that are wired today.

Command Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Command
     - Purpose
   * - ``pycsamt interp rocks``
     - List the built-in rock database or classify one resistivity value.
   * - ``pycsamt interp classify WORKDIR``
     - Load an Occam2D inversion result and print pseudo-stratigraphic logs.
   * - ``pycsamt interp export WORKDIR``
     - Export interpreted logs or grids to ``xyz``, ``las``, ``csv``, ``vtk``,
       ``surfer_grid``, or ``surfer_xyz``.

Global Conventions
------------------

``--format {text,json,csv}``
    Output format for commands that print tables. This option is used by
    ``rocks`` and ``classify``.

``-v`` / ``--verbose``
    Increase logging. For ``classify`` and ``export``, verbose mode reports
    solver detection, loaded iteration information, and station counts.

``--no-color``
    Disable colored logs.

Resistivity values passed to ``interp`` are linear Ohm.m values unless a
field explicitly says otherwise. Internal model arrays may be log-scaled;
the CLI handles that through ``ResistivityModel.from_occam2d`` and the
calibrator.

``interp rocks``
----------------

``pycsamt interp rocks`` is the quickest way to inspect the lithology
catalog used by the interpretation layer. It uses
``pycsamt.geology.lithology.RockDatabase.default()``.

List the catalog:

.. code-block:: console

   pycsamt interp rocks
   pycsamt interp rocks --format json
   pycsamt interp rocks --format csv > rock_catalog.csv

Classify one resistivity value:

.. code-block:: console

   pycsamt interp rocks --rho 250
   pycsamt interp rocks --rho 0.05 --method overlap
   pycsamt interp rocks --rho 5000 --format json

Catalog Output
~~~~~~~~~~~~~~

When ``--rho`` is omitted, the command prints all database entries sorted
by ``rho_min``. Output columns are:

``code``
    Short rock/unit code.

``name``
    Human-readable rock or hydrogeological unit name.

``rho_min`` / ``rho_max``
    Linear resistivity range in Ohm.m.

``rho_mid``
    Geometric middle of the range, rounded for display.

``description``
    Short interpretation note.

Single-Value Output
~~~~~~~~~~~~~~~~~~~

When ``--rho`` is supplied, it must be positive. The result contains:

``query_rho_Ohm_m``
    The value you asked to classify.

``code`` / ``name``
    Matched database entry.

``rho_min_Ohm_m`` / ``rho_max_Ohm_m`` / ``rho_mid_Ohm_m``
    Resistivity range for the matched entry.

``description``
    Interpretation note for the matched entry.

Classification Methods
~~~~~~~~~~~~~~~~~~~~~~

``nearest``
    Default. Selects the database entry with the closest geometric mean
    resistivity.

``overlap``
    Selects the first database entry whose resistivity range brackets the
    query value. This is useful when you want range membership instead of
    nearest-center matching.

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Option
     - Default
     - Meaning
   * - ``--rho``, ``-r``
     - none
     - Positive linear resistivity value in Ohm.m. Omit to list the catalog.
   * - ``--method``
     - ``nearest``
     - One of ``nearest`` or ``overlap``.
   * - ``--format``
     - ``text``
     - Print text, JSON, or CSV.
   * - ``-v``, ``--verbose``
     - ``0``
     - Increase logging verbosity.
   * - ``--no-color``
     - false
     - Disable colored logs.

``interp classify``
-------------------

``pycsamt interp classify WORKDIR`` loads an inversion work directory,
converts the inversion result to a ``ResistivityModel``, classifies the
model against the built-in rock database, and prints a flat layer table.

Basic use:

.. code-block:: console

   pycsamt interp classify occam_run/
   pycsamt interp classify occam_run/ --format csv
   pycsamt interp classify occam_run/ --format json

Select an iteration or save a CSV copy:

.. code-block:: console

   pycsamt interp classify occam_run/ --iteration 10
   pycsamt interp classify occam_run/ --merge-tol 0.1 --output logs.csv
   pycsamt interp classify occam_run/ --max-depth 1000 --format json

Workflow
~~~~~~~~

The command does the following:

1. Resolve the solver. With ``--solver auto``, a workdir is considered
   Occam2D when it contains at least one Occam signature file:
   ``Occam2DMesh``, ``Occam2DModel``, ``OccamData.dat``, or
   ``OccamStartup``.
2. Load ``pycsamt.models.occam2d.results.InversionResult`` from the
   workdir. If ``--iteration`` is omitted, the result loader uses its
   default final/last iteration behavior.
3. Convert the native result with ``ResistivityModel.from_occam2d``.
4. Optionally crop by ``--max-depth`` before classification.
5. Fit ``ModelCalibrator`` with the requested ``--ptol``.
6. Build one stratigraphic log per station using
   ``stratigraphic_logs(merge_tolerance=...)``.
7. Emit the flattened layer table, and optionally save it to
   ``--output`` as CSV.

Layer Table
~~~~~~~~~~~

``classify`` emits one row per interpreted layer. Columns are:

``station``
    Station name from the resistivity model.

``depth_top_m`` / ``depth_bot_m``
    Top and bottom depth in metres, rounded to two decimals.

``thickness_m``
    ``depth_bot_m - depth_top_m``.

``rock``
    Matched rock/unit name.

``rho_min`` / ``rho_max``
    Resistivity range for the matched rock entry.

``description``
    Rock database description.

Interpretation Controls
~~~~~~~~~~~~~~~~~~~~~~~

``--ptol``
    Borehole-calibration tolerance fraction passed to ``ModelCalibrator``.
    The current CLI path fits without external boreholes, but the parameter
    is preserved for the calibrator interface. Valid range is ``0`` to ``1``.

``--merge-tol``
    Log10-rho difference threshold for merging adjacent cells into one
    stratigraphic layer. Smaller values preserve more thin layers; larger
    values produce smoother logs. Minimum is ``0``.

``--max-depth``
    Clips the model before classification by keeping only cells with
    ``z_centers <= max_depth``.

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Option
     - Default
     - Meaning
   * - ``WORKDIR``
     - required
     - Existing inversion work directory.
   * - ``--solver``
     - ``auto``
     - ``auto`` or ``occam2d``. Only Occam2D is currently supported.
   * - ``--iteration``, ``-i``
     - last/default
     - Positive iteration number to load.
   * - ``--ptol``
     - ``0.10``
     - Calibrator tolerance fraction, from ``0`` to ``1``.
   * - ``--merge-tol``
     - ``0.2``
     - Log10-rho threshold for adjacent-layer merging.
   * - ``--max-depth``
     - none
     - Maximum depth in metres to classify.
   * - ``--output``, ``-o``
     - none
     - Save the layer table to this CSV file.
   * - ``--format``
     - ``text``
     - Print text, JSON, or CSV to stdout.
   * - ``-v``, ``--verbose``
     - ``0``
     - Print detection/classification metadata to stderr.
   * - ``--no-color``
     - false
     - Disable colored logs.

Important Boundaries
~~~~~~~~~~~~~~~~~~~~

``classify`` currently rejects unknown solver directories when
``--solver auto`` cannot detect Occam2D. It also rejects any effective
solver other than ``occam2d``. ModEM and MARE2DEM interpretation can be
added later, but they are not part of this CLI command today.

``interp export``
-----------------

``pycsamt interp export WORKDIR`` runs the same load/classify pipeline as
``classify`` and writes one or more files to ``--output-dir``.

The export format is required:

.. code-block:: console

   pycsamt interp export occam_run/ --format xyz --output-dir exports/
   pycsamt interp export occam_run/ --format las --output-dir las_logs/
   pycsamt interp export occam_run/ --format csv --output-dir exports/
   pycsamt interp export occam_run/ --format vtk --output-dir vtk_out/
   pycsamt interp export occam_run/ --format surfer_grid --output-dir surfer/
   pycsamt interp export occam_run/ --format surfer_xyz --output-dir surfer/

Overwrite existing output:

.. code-block:: console

   pycsamt interp export occam_run/ --format csv --output-dir exports/ --overwrite

Formats
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 14 28 58

   * - Format
     - File name(s)
     - Purpose
   * - ``xyz``
     - ``profile.xyz``
     - Oasis Montaj / GIS-compatible interpreted profile.
   * - ``las``
     - ``<station>.las``
     - One CWLS LAS 2.0-style well log per station.
   * - ``csv``
     - ``layers.csv``
     - Flat interpreted layer table for all stations.
   * - ``vtk``
     - ``model.vtk``
     - ASCII rectilinear grid for ParaView/QGIS-style inspection.
   * - ``surfer_grid``
     - ``model.grd``
     - Golden Software Surfer DSAA regular grid — resampled onto a
       uniform grid and opens directly in Surfer, no further gridding
       needed.
   * - ``surfer_xyz``
     - ``model.dat``
     - Golden Software Surfer scattered ``X Y Z`` points at the model's
       real cell centres — no resampling; Surfer grids it with
       whatever method the user chooses.

Export Pipeline
~~~~~~~~~~~~~~~

For all formats, the command:

1. detects or validates Occam2D;
2. loads ``InversionResult(workdir, iteration=iteration)``;
3. converts to ``ResistivityModel.from_occam2d``;
4. optionally crops by ``--max-depth``;
5. fits ``ModelCalibrator(ptol=ptol)``;
6. builds stratigraphic logs with ``merge_tolerance=merge_tol``;
7. dispatches to ``pycsamt.interp.export`` writer functions:
   ``to_oasis_montaj_xyz``, ``to_las``, ``to_csv``, ``to_vtk``,
   ``to_surfer_grid``, or ``to_surfer_xyz``.

The ``surfer_grid``/``surfer_xyz`` writers accept *any* 2-D inversion
result — not just Occam2D — via
:meth:`~pycsamt.interp.ResistivityModel.from_any`; the CLI command
above is Occam2D-only today only because ``interp export`` as a whole
is (``--solver`` only accepts ``auto``/``occam2d``, above), not because
the writers themselves are limited. Called directly from Python, they
also accept ModEM, the backend-neutral
:class:`~pycsamt.inversion.results.InversionResult`, an AI agent
result, or a raw array, and an optional real ``elevation`` to drape
real terrain instead of a flat datum.

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Option
     - Default
     - Meaning
   * - ``WORKDIR``
     - required
     - Existing inversion work directory.
   * - ``--solver``
     - ``auto``
     - ``auto`` or ``occam2d``. Only Occam2D is currently supported.
   * - ``--iteration``, ``-i``
     - last/default
     - Positive iteration number to load.
   * - ``--format``, ``-f``
     - required
     - One of ``xyz``, ``las``, ``csv``, ``vtk``, ``surfer_grid``, or
       ``surfer_xyz``.
   * - ``--output-dir``, ``-o``
     - current directory
     - Directory where output files are written.
   * - ``--ptol``
     - ``0.10``
     - Calibrator tolerance fraction, from ``0`` to ``1``.
   * - ``--merge-tol``
     - ``0.2``
     - Log10-rho threshold for adjacent-layer merging.
   * - ``--max-depth``
     - none
     - Maximum depth in metres to export.
   * - ``--overwrite``
     - false
     - Replace existing output files.
   * - ``-v``, ``--verbose``
     - ``0``
     - Print classification/export metadata to stderr.
   * - ``--no-color``
     - false
     - Disable colored logs.

Overwrite Rules
~~~~~~~~~~~~~~~

The exporter refuses to replace an existing output file unless
``--overwrite`` is supplied. This applies to:

* ``profile.xyz`` for ``xyz``;
* each ``<station>.las`` file for ``las``;
* ``layers.csv`` for ``csv``;
* ``model.vtk`` for ``vtk``.

Typical Workflows
-----------------

Start by learning the resistivity-to-rock mapping:

.. code-block:: console

   pycsamt interp rocks --rho 250
   pycsamt interp rocks --rho 250 --format json

Inspect an Occam2D result before exporting:

.. code-block:: console

   pycsamt interp classify occam_run/ --merge-tol 0.15
   pycsamt interp classify occam_run/ --merge-tol 0.15 --output layers.csv

Export station logs and a 3-D-viewable grid:

.. code-block:: console

   pycsamt interp export occam_run/ --format las --output-dir las_logs/
   pycsamt interp export occam_run/ --format vtk --output-dir vtk/

Clip to the shallow section for deliverables:

.. code-block:: console

   pycsamt interp classify occam_run/ --max-depth 1000 --format csv > shallow_layers.csv
   pycsamt interp export occam_run/ --max-depth 1000 --format csv --output-dir shallow/

Troubleshooting
---------------

``Resistivity must be a positive value``
    ``interp rocks --rho`` received ``0`` or a negative value. Use a
    positive linear Ohm.m value.

``Could not detect an inversion solver``
    The work directory does not contain an Occam2D signature file. Pass
    ``--solver occam2d`` only if the directory is truly an Occam2D result
    with files that the result loader can parse.

``Only occam2d is currently available``
    ``classify`` and ``export`` do not yet support the detected/requested
    solver.

``already exists. Pass --overwrite to replace it.``
    The export target exists. Choose a new ``--output-dir`` or add
    ``--overwrite``.

``Error loading inversion result`` / ``Error loading / classifying model``
    The solver was accepted, but the result loader, model conversion, depth
    crop, or calibrator failed. Check that the workdir has a readable
    Occam2D result and that ``--iteration`` exists.

API Equivalents
---------------

Rock lookup:

.. code-block:: python

   from pycsamt.geology.lithology import RockDatabase

   db = RockDatabase.default()
   entry = db.classify(250.0, method="nearest")
   print(entry.name, entry.rho_min, entry.rho_max)

Classify an Occam2D model:

.. code-block:: python

   from pycsamt.models.occam2d.results import InversionResult
   from pycsamt.interp import ResistivityModel, ModelCalibrator

   result = InversionResult("occam_run")
   model = ResistivityModel.from_occam2d(result)
   calibrator = ModelCalibrator(ptol=0.10).fit(model)
   logs = calibrator.stratigraphic_logs(merge_tolerance=0.2)

Export logs:

.. code-block:: python

   from pycsamt.interp import export

   export.to_csv(logs, "layers.csv")
   export.to_oasis_montaj_xyz(logs, "profile.xyz")
   for log in logs:
       export.to_las(log, f"{log.station_name}.las")
   export.to_vtk(model, "model.vtk")

