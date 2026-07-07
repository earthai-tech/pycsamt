.. _models_mare2dem:

MARE2DEM
========

:mod:`pycsamt.models.mare2dem` provides the pyCSAMT integration layer for
MARE2DEM, a 2.5-D finite-element electromagnetic modelling and inversion
code. MARE2DEM supports magnetotelluric (MT) and controlled-source
electromagnetic (CSEM) workflows, adaptive triangular meshes, topography, and
MPI execution.

pyCSAMT does not vendor the compiled MARE2DEM executable. Instead, it provides
the tools needed to manage a MARE2DEM project around that executable:

* configuration templates for source, binary, MPI, and file-name settings;
* source download/build/location helpers;
* native ``.emdata``, ``.resistivity``, ``.poly``, ``.settings``, log, group
  RMS, and data-group readers/writers;
* survey builders for MT and CSEM synthetic or prepared data files;
* ZMM-to-MARE2DEM MT data conversion;
* geometry utilities for topography, UTM conversion, profile projection,
  polygon simplification, triangle-region assignment, and area-of-interest
  estimation;
* input builders, runners, result loaders, plotting helpers, model-difference
  utilities, merge tools, and synthetic-noise helpers.

This page is a practical guide to those pieces. It focuses on how a user should
prepare, run, inspect, and archive a MARE2DEM project from pyCSAMT.

When To Use MARE2DEM
--------------------

Use the MARE2DEM integration when a project needs a native 2.5-D finite-element
workflow rather than a pyCSAMT built-in inversion. Common reasons include:

* a survey geometry that is naturally profile-based but not adequately treated
  by a simple 1-D or 2-D approximation;
* MT, CSEM, or combined data that must be represented in MARE2DEM's native
  ``.emdata`` format;
* topography, seafloor, receiver, or transmitter geometry that should be
  explicitly represented;
* adaptive triangular discretization through Triangle/PSLG ``.poly`` geometry;
* existing MARE2DEM project files that need Python-side validation,
  conversion, plotting, or result loading;
* an HPC workflow where native input files are prepared locally, run on a
  cluster, and loaded back into pyCSAMT afterward.

MARE2DEM is not the fastest path for every inversion. If you only need a
high-level backend-neutral 2-D run, start with :doc:`choosing_backend`.
MARE2DEM is most valuable when the native file set and engine-specific control
are part of the scientific workflow.

Package Map
-----------

The public MARE2DEM package surface is intentionally broad. It includes low
level file readers, high level lifecycle helpers, and project utilities.

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Area
     - Main objects
     - Purpose
   * - Configuration
     - ``Mare2DEMConfig``
     - Stores source location, compiler overrides, binary/MPI settings,
       inversion controls, initial model value, and default file names.
   * - Source and binary
     - ``SourceManager``
     - Locates, downloads, builds, and reports status for the external
       MARE2DEM source tree and executable.
   * - Build inputs
     - ``InputBuilder``
     - Writes a MARE2DEM run directory from an existing ``.emdata`` file, an
       ``EMDataFile`` object, or MT/CSEM survey configuration objects.
   * - Execute
     - ``Mare2DEMRunner``
     - Builds the command line, handles MPI options, runs the executable from
       the working directory, and can return an ``InversionResult``.
   * - Load results
     - ``InversionResult``, ``Mare2DEMLog``, ``GroupRMSLog``
     - Scans output directories, parses iteration logs, exposes RMS history,
       convergence state, final model, observed data, and response data.
   * - Native I/O
     - ``read_emdata``, ``write_emdata``, ``read_resistivity``,
       ``write_resistivity``, ``read_poly``, ``write_poly``,
       ``write_settings``
     - Reads and writes MARE2DEM-native project files.
   * - Data management
     - ``MTSurveyConfig``, ``CSEMSurveyConfig``, ``make_data_file``,
       ``read_zmm``, ``make_mt_data_from_zmm``, ``merge_data_files``
     - Builds, converts, and merges MT/CSEM data products.
   * - Geometry
     - ``parse_topo``, ``lonlat_to_utm``, ``project_onto_line``,
       ``get_line_orientation``, ``simplify_poly``, ``get_triangle_regions``
     - Handles topography, survey projection, coordinate conversion, polygon
       cleanup, and triangle-region utilities.
   * - QC and interpretation
     - ``NoiseConfig``, ``add_synthetic_noise``, ``diff_resistivity``,
       ``PlotConvergence``, ``PlotSurveyLayout``, ``plot_poly``
     - Supports synthetic tests, model comparison, convergence review, survey
       layout QC, and mesh/geometry inspection.

Configuration
-------------

``Mare2DEMConfig`` is the source-of-truth object for MARE2DEM runs. It is a
plain dataclass, so it can be created in Python, written as a template, edited
outside Python, and loaded again.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import Mare2DEMConfig

   cfg = Mare2DEMConfig(
       source_dir="/opt/mare2dem/source",
       binary="MARE2DEM",
       use_mpi=True,
       n_procs=16,
       mpi_command="mpirun",
       max_iterations=120,
       target_rms=1.0,
       initial_rho=10.0,
       data_file="line12.emdata",
       resistivity_file="line12.resistivity",
       settings_file="line12.settings",
   )

   cfg.to_template("runs/line12/mare2dem.yml")
   loaded = Mare2DEMConfig.from_file("runs/line12/mare2dem.yml")

The configuration groups five concerns.

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Concern
     - Fields
     - Meaning
   * - Source management
     - ``source_dir``, ``fc_compiler``, ``cc_compiler``
     - Where the Fortran source lives and which compiler commands should be
       used when building.
   * - Binary and MPI
     - ``binary``, ``use_mpi``, ``n_procs``, ``mpi_command``
     - How pyCSAMT should locate and launch the executable.
   * - Inversion control
     - ``max_iterations``, ``target_rms``
     - Iteration limit and normalized RMS target written into the native
       model/settings files.
   * - Initial model
     - ``initial_rho``
     - Starting homogeneous half-space resistivity, in ohm metres.
   * - File names
     - ``data_file``, ``resistivity_file``, ``settings_file``
     - Native file names used by builders and runners.

The ``resistivity_stem`` property returns the stem of
``resistivity_file``. MARE2DEM receives this stem on the command line and then
derives related filenames from it. For example, ``line12.resistivity`` is
passed as ``line12``.

Source And Binary Management
----------------------------

``SourceManager`` handles the external source tree. It resolves a source
directory in this order:

1. the explicit ``source_dir`` argument passed to ``SourceManager``;
2. ``Mare2DEMConfig.source_dir``;
3. the ``PYCSAMT_MARE2DEM_SOURCE`` environment variable;
4. the package ``_source/`` directory when it is writable, which is common in
   editable development installs;
5. a platform user-data directory.

Use ``status`` first, before trying to download or build anything.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import Mare2DEMConfig, SourceManager

   cfg = Mare2DEMConfig(source_dir="/opt/mare2dem/source")
   source = SourceManager(config=cfg, verbose=1)

   status = source.status()
   print(status["source_dir"])
   print(status["downloaded"])
   print(status["built"])
   print(status["binary_path"])

When source code is not present, ``download`` can use Git or a source archive.
When source code is present but the binary is missing, ``build`` compiles the
external code.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import SourceManager

   source = SourceManager(source_dir="/opt/mare2dem/source", verbose=1)

   # Network access and compiler availability are environment-dependent.
   # source.download(method="auto")
   # source.build(clean_first=False)

   source.print_status()

Compilation is system-dependent. In many environments MARE2DEM requires MPI
Fortran/C tooling and Intel MKL/ScaLAPACK/BLACS support. On Windows, build from
WSL or another Unix-like environment rather than a native Windows shell.

Binary Resolution
-----------------

``Mare2DEMRunner`` resolves the executable through:

1. ``Mare2DEMConfig.binary`` on ``PATH``;
2. ``<source_dir>/<binary>``;
3. the platform user-data binary location.

Use ``runner.command`` as a dry-run check before launching an inversion.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import Mare2DEMConfig, Mare2DEMRunner

   cfg = Mare2DEMConfig(
       binary="MARE2DEM",
       use_mpi=True,
       n_procs=8,
       mpi_command="mpirun",
       resistivity_file="line12.resistivity",
   )

   runner = Mare2DEMRunner("runs/line12/native", config=cfg)
   print(runner.command("line12"))

For cluster workflows, put the command string and the loaded module list in the
run provenance file. The same native directory can then be executed by a job
scheduler and loaded later with ``InversionResult``.

Native Files
------------

MARE2DEM projects revolve around native text files. pyCSAMT treats these as
first-class scientific records, not temporary build artifacts.

.. list-table::
   :header-rows: 1
   :widths: 22 28 50

   * - File
     - Reader/writer
     - Role
   * - ``.emdata``
     - ``read_emdata``, ``write_emdata``
     - Observed or synthetic MT/CSEM/DC data, receiver/transmitter metadata,
       UTM origin, frequencies, and data rows.
   * - ``*_MARE2DEM.emdata`` or ``.resp``
     - ``read_emdata``
     - Predicted response data produced by the engine.
   * - ``.resistivity``
     - ``read_resistivity``, ``write_resistivity``
     - Resistivity parameters, free/fixed flags, bounds, prejudice values, and
       references to data, settings, and polygon files.
   * - ``.poly``
     - ``read_poly``, ``write_poly``
     - Triangle PSLG geometry: nodes, segments, holes, and regions.
   * - ``.settings``
     - ``write_settings``
     - Parallel decomposition and inversion settings.
   * - ``.emdata_group``
     - ``read_data_group``, ``write_data_group``
     - Group definitions used for grouped RMS diagnostics.
   * - Group RMS logs
     - ``read_group_rms_log``
     - Per-group RMS evolution, useful for diagnosing which data families are
       controlling the inversion.
   * - Iteration logs
     - ``Mare2DEMLog``
     - Iteration number, RMS misfit, roughness, Lagrange multiplier, and
       convergence state.

The validation helpers classify common MARE2DEM files by suffix and naming
convention.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import detect_file_type, is_response_file

   print(detect_file_type("line12.emdata"))
   print(detect_file_type("line12.resistivity"))
   print(is_response_file("line12_MARE2DEM.emdata"))

Build A Run Directory
---------------------

``InputBuilder`` writes a minimal MARE2DEM input set. It can start from an
existing ``.emdata`` file, an in-memory ``EMDataFile``, or MT/CSEM survey
configuration objects.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import InputBuilder, Mare2DEMConfig

   cfg = Mare2DEMConfig(
       initial_rho=10.0,
       max_iterations=120,
       target_rms=1.0,
       data_file="line12.emdata",
       resistivity_file="line12.resistivity",
       settings_file="line12.settings",
   )

   builder = InputBuilder(config=cfg, verbose=1)
   files = builder.build(
       "prepared/line12.emdata",
       workdir="runs/line12/native",
   )

   print(files["data"])
   print(files["model"])
   print(files["settings"])

The builder writes:

* the data file, copied or generated into the run directory;
* a homogeneous starting ``.resistivity`` file based on ``initial_rho``;
* a ``.settings`` file.

For production inversions, the generated starting model is often only the first
step. Review the ``.resistivity`` and ``.poly`` inputs before launching the
external code.

Create MT Data
--------------

``MTSurveyConfig`` and ``make_data_file`` are useful for synthetic tests or
for constructing simple MT-native files from survey arrays.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.models.mare2dem import MTSurveyConfig, make_data_file

   mt = MTSurveyConfig(
       frequencies=np.logspace(-3, 3, 25),
       rx_y=np.linspace(-6000.0, 6000.0, 31),
       rx_type="land",
       lTE=True,
       lTM=True,
       lTipper=False,
   )

   em = make_data_file(
       "runs/line12/native/line12.emdata",
       topo=0.0,
       mt=mt,
   )

   print(em.n_data)

For real MT processing, pyCSAMT also includes ZMM readers and converters.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.models.mare2dem import make_mt_data_from_zmm, read_zmm

   station = read_zmm("zmm/station001.zmm")
   print(station.name)

   make_mt_data_from_zmm(
       sorted(Path("zmm").glob("*.zmm")),
       "runs/line12/native/line12.emdata",
       error_floor=0.05,
   )

Check the generated data rows carefully. Conversion utilities help with file
mechanics, but the interpreter still owns station selection, component choice,
frequency band selection, and error-floor policy.

Create CSEM Data
----------------

``CSEMSurveyConfig`` builds a controlled-source survey with transmitter and
receiver layout parameters.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.models.mare2dem import CSEMSurveyConfig, make_data_file

   csem = CSEMSurveyConfig(
       frequencies=np.array([0.25, 0.5, 1.0, 2.0]),
       rx_y=np.linspace(-4000.0, 4000.0, 17),
       tx_y=np.array([-2500.0, 0.0, 2500.0]),
       rx_type="marine",
       tx_type="marine",
       lEx=True,
       lEy=True,
       lBx=True,
       lBy=True,
   )

   make_data_file(
       "runs/csem_line/native/csem_line.emdata",
       topo=-1000.0,
       csem=csem,
   )

For CSEM projects, geometry QC is essential. Plot receiver/transmitter
locations, confirm signs and offsets, and record whether coordinates are local
profile coordinates, UTM coordinates, or a transformed system.

Merge And Noise Utilities
-------------------------

MARE2DEM workflows often require combining data families or creating synthetic
observations from a forward response. pyCSAMT provides utilities for both.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import merge_data_files

   merged = merge_data_files(
       ["runs/mt/native/mt.emdata", "runs/csem/native/csem.emdata"],
       "runs/joint/native/joint.emdata",
   )

   print(merged.n_data)

Synthetic noise is useful for controlled tests and algorithm validation.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import (
       NoiseConfig,
       add_synthetic_noise,
       read_emdata,
       write_emdata,
   )

   response = read_emdata("runs/forward/native/forward_MARE2DEM.emdata")
   noise = NoiseConfig(
       mt_rel_noise=0.05,
       mt_abs_noise_tipper=0.01,
   )

   synthetic = add_synthetic_noise(response, noise)
   write_emdata(synthetic, "runs/synthetic/native/synthetic.emdata")

Do not treat noisy synthetic data as field data. Keep synthetic sources,
random seeds, and noise parameters in the provenance notes.

Geometry And Topography
-----------------------

The MARE2DEM integration includes geometry helpers because geometry mistakes
are one of the easiest ways to produce convincing but wrong models.

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Helper
     - Use
   * - ``parse_topo``, ``topo_depth``, ``topo_slope``
     - Read and interpolate topography/seafloor profiles and compute local
       slopes.
   * - ``lonlat_to_utm``, ``utm_to_lonlat``
     - Convert coordinates with ``pyproj`` when available and a pure-Python
       WGS-84 fallback otherwise.
   * - ``get_line_orientation``, ``project_onto_line``
     - Estimate profile orientation and project stations onto a survey line.
   * - ``dp_simplify``, ``simplify_poly``
     - Simplify polylines and remove collinear polygon nodes.
   * - ``get_intersections``, ``do_rects_overlap``
     - Find segment intersections with bounding-box pre-filtering.
   * - ``triangle_centroids``, ``triangle_areas``, ``get_centroids``
     - Compute area-weighted centroids and triangle geometry summaries.
   * - ``estimate_area_of_interest``
     - Estimate a practical modelling area from survey geometry.
   * - ``get_triangle_regions``
     - Assign finite-element regions for Triangle-based meshes.

Example profile projection:

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.models.mare2dem import (
       get_line_orientation,
       lonlat_to_utm,
       project_onto_line,
   )

   lon = np.array([11.501, 11.507, 11.514])
   lat = np.array([3.842, 3.845, 3.849])

   east, north, zone, hemisphere = lonlat_to_utm(lon, lat)
   azimuth = get_line_orientation(north, east)
   cross_profile, along_profile = project_onto_line(
       north,
       east,
       north[0],
       east[0],
       azimuth,
   )

   print(zone, hemisphere)
   print(along_profile)

For production documentation, include a small figure or table showing original
coordinates, projected distances, and topography values. That record is often
as important as the inversion result.

Grid And Mesh Utilities
-----------------------

The package includes helpers for constructing MARE2DEM geometry and models
from gridded resistivity information.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import grid_to_mare2dem

   files = grid_to_mare2dem(
       grid="models/resistivity_grid.npy",
       out_prefix="runs/grid_model/native/line12",
   )

   print(files["poly"])
   print(files["resistivity"])

The exact grid format and metadata requirements depend on the upstream grid
object. Always inspect the written ``.poly`` and ``.resistivity`` files before
using them in an inversion.

Run MARE2DEM
------------

Once native files are prepared and the executable is available, use
``Mare2DEMRunner`` to launch the run.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import (
       InversionResult,
       Mare2DEMConfig,
       Mare2DEMRunner,
   )

   cfg = Mare2DEMConfig.from_file("runs/line12/mare2dem.yml")
   runner = Mare2DEMRunner("runs/line12/native", config=cfg, verbose=1)

   # Check the command before running it.
   print(runner.command(cfg.resistivity_stem))

   # Run locally. On a cluster, submit the same command through the scheduler.
   # result = runner.run(cfg.resistivity_stem, timeout=None)

   # Load an existing completed run.
   result = InversionResult("runs/line12/native", config=cfg)
   print(result.summary())

``run`` returns an ``InversionResult`` by default. Set ``load_result=False`` if
you only want process execution and will load results later.

.. code-block:: python
   :linenos:

   result = runner.run(
       cfg.resistivity_stem,
       use_mpi=True,
       n_procs=16,
       extra_args=None,
       timeout=None,
       load_result=True,
   )

For long inversions, prefer scheduler-managed execution. Build the native
directory with pyCSAMT, submit the command on the cluster, then use
``InversionResult`` after the files are complete.

Inspect Results
---------------

``InversionResult`` scans a run directory for recognized output files and
exposes the main products.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import InversionResult

   result = InversionResult("runs/line12/native")

   print(result.converged)
   print(result.final_rms)
   print(result.n_iterations)
   print(result.model)
   print(result.data)
   print(result.response)

The log parser exposes the RMS history, roughness history, and regularization
multiplier history.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import Mare2DEMLog

   log = Mare2DEMLog("runs/line12/native/OccamLog.2012.0")

   print(log.final_rms)
   print(log.converged)
   print(log.rms_history())
   print(log.roughness_history())

Group RMS logs should be reviewed when the inversion combines data types,
components, stations, or source groups.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import read_group_rms_log

   group_log = read_group_rms_log("runs/line12/native/line12.group_rms.log")
   print(group_log.headers)
   print(group_log.rms_log[-1])

A low total RMS can hide a poor fit to one group. Inspect total RMS, group RMS,
response residuals, and model roughness together.

Plotting And QC
---------------

The plotting helpers are designed for quality control and interpretation.

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Plot helper
     - Use
   * - ``PlotConvergence``
     - RMS and convergence history from logs.
   * - ``PlotSurveyLayout``
     - Receiver/transmitter layout in map/profile coordinates.
   * - ``PlotRxParams``
     - Receiver geometry and parameter checks.
   * - ``PlotTxParams``
     - Transmitter geometry and parameter checks.
   * - ``plot_poly``
     - PSLG/``.poly`` geometry inspection.
   * - ``PlotModel``
     - Model-section visualization support.
   * - ``PlotResponse``
     - Observed/predicted response inspection.

Example convergence review:

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import Mare2DEMLog, PlotConvergence

   log = Mare2DEMLog("runs/line12/native/OccamLog.2012.0")
   plotter = PlotConvergence(log)
   plotter.plot()

Before presenting a final model, review at least:

* station and transmitter geometry;
* topography or bathymetry representation;
* mesh or polygon geometry;
* RMS history and convergence status;
* group RMS history;
* observed/predicted response fits;
* model smoothness, bounds, and fixed/free parameter behavior.

Model Comparison
----------------

``diff_resistivity`` compares two ``.resistivity`` files, typically to inspect
how a model changed between iterations, parameter choices, or preprocessing
decisions.

.. code-block:: python
   :linenos:

   from pycsamt.models.mare2dem import diff_resistivity

   diff = diff_resistivity(
       "runs/line12_v01/native/line12.resistivity",
       "runs/line12_v02/native/line12.resistivity",
       "runs/comparison/line12_v02_minus_v01.resistivity",
   )

   print(diff)

Use model differences as diagnostics, not as standalone geological evidence.
Pair them with data-fit changes and model regularization changes.

Recommended Run Layout
----------------------

Keep source/build material separate from scientific run directories.

.. code-block:: text
   :linenos:

   external/
     mare2dem/
       source/
       build_notes.md

   runs/
     survey_alpha/
       line12_mare2dem_v01/
         mare2dem.yml
         provenance.yml
         native/
           line12.emdata
           line12.poly
           line12.resistivity
           line12.settings
           OccamLog.2012.0
           line12_MARE2DEM.emdata
           line12.group_rms.log
         qc/
           survey_layout.png
           convergence.png
           response_fit.png
         exports/
           line12_snapshot.zip

Use versioned run names. MARE2DEM uses simple native filenames, so running
several experiments in one directory can easily mix stale outputs with new
inputs.

Pre-Run Checklist
-----------------

Before launch:

* load ``Mare2DEMConfig`` from the edited template;
* confirm source directory and executable resolution;
* check ``runner.command`` and record it;
* confirm MPI process count and scheduler settings;
* verify that ``.emdata``, ``.resistivity``, ``.settings``, and ``.poly`` files
  exist when required;
* inspect station, transmitter, and profile coordinates;
* inspect topography or seafloor representation;
* confirm frequency band, components, and error floors;
* move old output files out of the run directory;
* record pyCSAMT version, MARE2DEM source/binary path, and compiler/MPI
  context.

Post-Run Checklist
------------------

After completion:

* read the main log before plotting;
* check ``result.converged``, ``result.final_rms``, and iteration count;
* review RMS and roughness histories;
* review group RMS where groups are used;
* confirm response-file timestamps match the intended run;
* compare observed and predicted responses;
* inspect the final model against bounds and fixed/free parameter flags;
* archive native input and output files with the configuration and provenance.

Common Mistakes
---------------

Using a stale response file
    If the binary fails, an old ``*_MARE2DEM.emdata`` file can remain in the
    directory. Check timestamps and logs before loading results.

Passing the wrong stem
    MARE2DEM receives the resistivity stem, not usually the full project path.
    ``Mare2DEMRunner`` normalizes the stem, but users should still keep
    filenames consistent.

Mixing coordinate systems
    Longitude/latitude, UTM, and local profile distance are different records.
    State which system is used in every generated native file.

Treating total RMS as enough
    Total RMS can hide poor fits to one component, station range, source, or
    data family. Review group RMS and response residuals.

Ignoring build provenance
    MARE2DEM behavior can depend on source revision, compiler, MKL/ScaLAPACK,
    MPI runtime, and cluster environment. Record them.

Next Steps
----------

* :doc:`configuration_and_io` for source-of-truth configuration and archive
  practice.
* :doc:`choosing_backend` for deciding when MARE2DEM is the right integration.
* :ref:`inversion_concepts` for regularized inversion and misfit concepts.
* :doc:`../api/models` for generated API reference pages.
