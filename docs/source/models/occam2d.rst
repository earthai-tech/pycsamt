.. _models_occam2d:

Occam2D
=======

:mod:`pycsamt.models.occam2d` provides the pyCSAMT interface to an
Occam2DMT-style smooth 2-D inversion workflow. It prepares native input files,
can run the external Fortran executable, loads iteration and response outputs,
and provides plotting helpers for models, responses, pseudosections, misfit
curves, station-level residuals, and 1-D station extracts.

Occam inversion is deliberately conservative. It seeks the smoothest model that
fits the data to an acceptable normalized RMS misfit. For a residual vector
:math:`r`, the usual RMS diagnostic is

.. math::

   \mathrm{RMS} =
   \sqrt{\frac{1}{N}\sum_{i=1}^{N} r_i^2}.

Occam iteration files store model parameters as log10 resistivity,
:math:`m = \log_{10}(\rho)`. pyCSAMT preserves that convention in
``InversionResult.rho_2d`` and converts to physical resistivity only when a
plot or downstream calculation needs it.

When To Use Occam2D
-------------------

Occam2D is a good first production engine when:

* stations follow a profile and can be represented by 2-D geometry;
* the dominant resistivity structure is expected to be approximately 2-D;
* a smooth model is scientifically appropriate;
* TE/TM apparent resistivity and phase are the primary inversion data;
* the deliverable must include native Occam data, mesh, model, startup,
  iteration, response, and log files;
* the interpreter wants a reproducible external-code workflow rather than only
  an in-memory inversion result.

Occam2D is not a blocky geology engine. A sharp geological contact may appear
as a smooth gradient because the regularization intentionally spreads structure
unless the data require a sharper transition.

Package Map
-----------

The Occam2D package is organized around the native inversion lifecycle.

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Area
     - Main objects
     - Purpose
   * - Configuration
     - ``OccamConfig``
     - Stores data selection, error floors, mesh options, startup controls,
       file names, and binary discovery name.
   * - Input construction
     - ``InputBuilder``
     - Builds data, mesh, model, and startup files from a survey source.
   * - Native data
     - ``OccamData``
     - Reads/writes Occam data files and builds data rows from EDI/Sites-like
       sources.
   * - Mesh and model
     - ``OccamMesh``, ``OccamModel``
     - Build/read finite-element mesh geometry and the model-parameter mapping.
   * - Startup and iterations
     - ``OccamStartup``, ``OccamIter``
     - Represent the iteration-zero startup file and non-zero iteration output
       files.
   * - Execution
     - ``OccamRunner``
     - Finds or compiles the executable, patches startup controls when
       requested, launches the solver, and captures stdout/stderr logs.
   * - Results
     - ``InversionResult``
     - Scans a completed run, loads mesh/model/data/log/iteration/response
       files, reconstructs ``rho_2d``, and exports ``iter2dat``.
   * - Diagnostics
     - ``OccamResponse``, ``OccamLog``
     - Parse response residuals and convergence history.
   * - Plotting
     - ``PlotModel``, ``PlotResponse``, ``PlotPseudo``, ``PlotMisfit``,
       ``PlotSounding1D``, ``PlotSiteMisfit``, ``PlotResponseGrid``,
       ``PlotStation1DFit``
     - Visual QC and interpretation views.
   * - Validation
     - ``detect_file_type`` and ``is_*`` helpers
     - Recognize data, mesh, model, startup, iteration, response, and log files.

Configuration
-------------

``OccamConfig`` is the source-of-truth object for a native Occam2D run. It can
be created in Python, written to a template, edited, and loaded again.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import OccamConfig

   cfg = OccamConfig(
       modes=["TE", "TM"],
       error_floor_rho=0.05,
       error_floor_phase=0.5,
       freq_min=0.1,
       freq_max=1000.0,
       n_layers=36,
       n_airlayers=5,
       cell_size_horizontal=75.0,
       cell_size_vertical_top=10.0,
       depth_scale=1.18,
       n_padding_x=8,
       max_iterations=80,
       target_misfit=1.0,
       initial_rho=100.0,
       data_file="OccamDataFile.dat",
       mesh_file="Occam2DMesh",
       model_file="Occam2DModel",
       startup_file="Startup",
       binary_name="Occam2D",
   )

   cfg.to_template("runs/profile_a/occam2d.yml")
   loaded = OccamConfig.from_file("runs/profile_a/occam2d.yml")

The configuration groups four concerns.

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Concern
     - Fields
     - Meaning
   * - Data selection
     - ``modes``, ``error_floor_rho``, ``error_floor_phase``,
       ``freq_min``, ``freq_max``
     - Which data rows are written and how minimum uncertainties are enforced.
   * - Mesh geometry
     - ``n_layers``, ``n_airlayers``, ``cell_size_horizontal``,
       ``cell_size_vertical_top``, ``depth_scale``, ``n_padding_x``
     - Finite-element discretization near stations, at depth, and near lateral
       boundaries.
   * - Startup controls
     - ``max_iterations``, ``target_misfit``, ``roughness_type``,
       ``diagonal_penalties``, ``stepsize_cut_count``, ``debug_level``,
       ``initial_rho``, ``lagrange_start``
     - Values written to the startup control file.
   * - Files and binary
     - ``data_file``, ``mesh_file``, ``model_file``, ``startup_file``,
       ``binary_name``
     - Native file names inside the run directory and executable name used by
       the runner.

Use strict loading for project work. Unknown keys usually mean spelling
mistakes or stale configuration files.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import OccamConfig

   cfg = OccamConfig.from_file("runs/profile_a/occam2d.yml")

   # Migration only: ignore retired or unknown keys while cleaning old files.
   migrated = OccamConfig.from_file(
       "runs/profile_a/old_occam2d.yml",
       strict=False,
   )

Native Files
------------

Occam2D projects should be archived as native file sets. The final image alone
is not enough to reproduce the inversion.

.. list-table::
   :header-rows: 1
   :widths: 22 26 52

   * - File
     - Object
     - Role
   * - ``OccamDataFile.dat``
     - ``OccamData``
     - Observed data rows, station names, offsets, frequencies, Occam type
       codes, datum values, and uncertainties.
   * - ``Occam2DMesh``
     - ``OccamMesh``
     - Finite-element mesh geometry, including air layers, earth layers,
       horizontal cells, and padding cells.
   * - ``Occam2DModel``
     - ``OccamModel``
     - Mapping from mesh cells to inversion parameters.
   * - ``Startup``
     - ``OccamStartup``
     - Iteration-zero control file passed to the executable: data/model/mesh
       names, target misfit, iteration controls, starting parameters.
   * - ``*.iter``
     - ``OccamIter``
     - Non-zero iteration output files containing log10-resistivity parameter
       values and iteration diagnostics.
   * - ``*.resp``
     - ``OccamResponse``
     - Modeled responses and residual information for an iteration.
   * - ``*.logfile`` or ``LogFile*``
     - ``OccamLog``
     - Convergence history and run-level diagnostic messages.
   * - ``occam_stdout.log`` and ``occam_stderr.log``
     - ``OccamRunner``
     - Captured process streams from pyCSAMT-launched runs.

The validation helpers can classify files when scanning a run directory.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d.validation import detect_file_type

   for path in [
       "runs/profile_a/native/OccamDataFile.dat",
       "runs/profile_a/native/Occam2DMesh",
       "runs/profile_a/native/Startup",
   ]:
       print(path, detect_file_type(path))

Build Input Files
-----------------

``InputBuilder`` constructs the four files required before an external Occam2D
run can start.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import InputBuilder, OccamConfig
   from pycsamt.site import Sites

   sites = Sites.from_dir("data/edi/profile_a")

   cfg = OccamConfig(
       modes=["TE", "TM"],
       freq_min=0.1,
       freq_max=1000.0,
       error_floor_rho=0.05,
       error_floor_phase=0.5,
       n_layers=32,
       target_misfit=1.0,
   )

   builder = InputBuilder(
       sites,
       workdir="runs/profile_a/native",
       config=cfg,
       verbose=1,
   )
   builder.build(title="Profile A Occam2D inversion")

   print(builder.summary())

The build chain is fixed:

1. ``OccamData.from_edi`` converts the survey source into Occam data rows.
2. ``OccamMesh.from_data`` builds a mesh from station offsets and mesh options.
3. ``OccamModel.from_mesh`` maps mesh cells to inversion parameters.
4. ``OccamStartup.from_model`` writes the initial parameter vector and run
   controls.

One-shot overrides passed to ``build`` update the stored configuration before
files are written.

.. code-block:: python
   :linenos:

   builder.build(
       modes=["TM"],
       n_layers=40,
       cell_size=50.0,
       error_floor_rho=0.07,
       freq_min=0.2,
       freq_max=500.0,
       title="Profile A TM-only sensitivity run",
   )

Because these overrides persist on ``builder.config``, write the resulting
configuration to the run directory if the build is retained.

Data Rows And Type Codes
------------------------

Occam2D data files written by pyCSAMT use TE/TM apparent resistivity and phase
rows. In the configuration:

* ``"TE"`` selects the :math:`Z_{xy}` component;
* ``"TM"`` selects the :math:`Z_{yx}` component;
* apparent resistivity is stored as :math:`\log_{10}(\rho_a)`;
* phase is stored in degrees;
* ``error_floor_rho`` is relative, for example ``0.05`` for five percent;
* ``error_floor_phase`` is absolute, in degrees.

Inspect the data object before running.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import OccamData

   data = OccamData.read("runs/profile_a/native/OccamDataFile.dat")

   print(data.n_sites)
   print(data.n_frequencies)
   print(data.n_data)
   print(data.site_names)
   print(data.offsets)

Station order and offsets matter because the mesh and pseudosections are built
around that profile geometry.

Mesh And Model Review
---------------------

The mesh and model determine what kind of smoothness the inversion can express.
Before launching a long external run, inspect:

* horizontal cell width near stations;
* number of padding cells on each side;
* top cell thickness and depth growth factor;
* number of air layers and earth layers;
* total number of model parameters;
* whether the mesh is far wider and deeper than the interpreted target.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import OccamMesh, OccamModel

   mesh = OccamMesh.read("runs/profile_a/native/Occam2DMesh")
   model = OccamModel.read("runs/profile_a/native/Occam2DModel")

   print(mesh.n_xcells, mesh.n_zcells)
   print(model.n_params)

Fine cells can improve near-surface representation, but they also increase
runtime and may exaggerate the apparent resolution of poorly constrained
structure. Padding moves boundaries away from the profile but also increases
mesh size.

Run Occam2D
-----------

``OccamRunner`` executes a prepared native directory. It does not build input
files; use ``InputBuilder`` first when starting from EDI data.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import OccamRunner

   runner = OccamRunner(
       workdir="runs/profile_a/native",
       binary_path="/usr/local/bin/Occam2D",
       startup_file="Startup",
       verbose=1,
   )

   runner.discover_binary(auto_compile=False)

   # Run only when the executable and native files are ready.
   # exit_code = runner.run(max_iter=80, target_misfit=1.0)

Binary discovery follows this order:

1. explicit ``binary_path``;
2. ``Occam2D`` or ``Occam2D.exe`` in the run directory;
3. executable on ``PATH``;
4. bundled ``_source`` directory, if ``auto_compile=True``.

Automatic compilation uses the bundled Fortran source and a compiler such as
``gfortran`` through ``make``. For reproducible production work, prefer an
explicit binary path and record compiler provenance separately.

``run`` can patch the startup file in place when ``max_iter`` or
``target_misfit`` is supplied. Archive the startup file that was actually run,
not only the template that created it.

Asynchronous execution is available for scripts that need to poll an external
process.

.. code-block:: python
   :linenos:

   runner = OccamRunner("runs/profile_a/native")

   # process = runner.run_async(auto_compile=False)
   # while runner.is_running:
   #     ...
   # exit_code = runner.wait()

For HPC usage, build and validate the native directory locally, then submit the
equivalent ``Occam2D Startup`` command through the scheduler. Load the completed
directory afterward with ``InversionResult``.

Backend-Neutral Occam2D Runs
----------------------------

The backend-neutral inversion API can drive the same native workflow through
``backend="occam2d"``.

.. code-block:: python
   :linenos:

   from pycsamt.inversion import InversionConfig, run_inversion

   cfg = InversionConfig(
       method="mt",
       dimension="2d",
       backend="occam2d",
       data="data/edi/profile_a",
       workdir="runs/profile_a",
       run_external=False,
       backend_options={
           "occam_config": {
               "modes": ["TE", "TM"],
               "n_layers": 32,
               "target_misfit": 1.0,
           },
       },
   )

   result = run_inversion(cfg)

With ``run_external=False``, pyCSAMT prepares or validates the native directory
without requiring the external binary to run. This is useful for documentation,
cluster workflows, and dry-run checks.

Load Results
------------

``InversionResult`` scans a completed run directory. If ``iteration`` is not
specified, it loads the highest numbered ``.iter`` file. It tries to match the
corresponding ``.resp`` file and reconstructs a log10-resistivity grid from the
mesh, model, and iteration vector.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import InversionResult

   result = InversionResult("runs/profile_a/native")

   print(result.summary())
   print(result.final_rms)
   print(result.n_iterations)
   print(result.rho_2d.shape if result.rho_2d is not None else None)

   selected = InversionResult("runs/profile_a/native", iteration=12)
   selected.iter2dat("runs/profile_a/exports/profile_a_iter12.dat")

The loader is tolerant of missing optional files. Missing logs or response
files leave the corresponding attributes as ``None``. A missing run directory
raises ``NotADirectoryError``.

Response And Misfit Diagnostics
-------------------------------

``OccamResponse`` reads modeled responses and residuals. Use it to understand
which sites, frequencies, or components are controlling the misfit.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import OccamResponse

   response = OccamResponse.read("runs/profile_a/native/RESP12.resp")

   print(response.rms)
   print(response.misfit_per_site())
   print(response.misfit_per_frequency())

Weighted residuals depend on the error column in the Occam data file. If error
floors are too small, the inversion may chase noise. If they are too large, the
model may stop before fitting reliable signal.

Log And Convergence
-------------------

``OccamLog`` parses convergence history from Occam log files. Use the log
together with the selected iteration file; do not judge a run only by the final
model image.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import OccamLog

   log = OccamLog.read("runs/profile_a/native/LogFile")

   print(log.converged)
   print(log.rms[-1] if log.rms.size else None)
   print(log.n_iter)

Review:

* starting RMS and final RMS;
* whether the run reached target misfit;
* whether roughness changes stabilize;
* whether the best-looking model corresponds to a sensible iteration;
* whether response residuals improve where the data are trustworthy.

Plotting And QC
---------------

The plotting helpers are designed to replace common Occam2DMT MATLAB
post-processing views.

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Plot helper
     - Use
   * - ``PlotModel`` or ``result.plot_model()``
     - Plot the reconstructed 2-D resistivity section.
   * - ``PlotResponse`` or ``result.plot_response()``
     - Compare observed and modeled responses.
   * - ``PlotPseudo`` or ``result.plot_pseudo()``
     - Plot observed-data pseudosections.
   * - ``PlotMisfit`` or ``result.plot_misfit()``
     - Plot RMS/convergence metrics by iteration.
   * - ``PlotSounding1D``
     - Extract station-centered 1-D profiles from the 2-D model.
   * - ``PlotSiteMisfit``
     - Plot per-site residual diagnostics.
   * - ``PlotResponseGrid``
     - Inspect response behavior across sites/frequencies.
   * - ``PlotStation1DFit`` and ``plot_station_1d_fit``
     - Review station-level 1-D fit style diagnostics.

Example:

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import InversionResult

   result = InversionResult("runs/profile_a/native")

   result.plot_misfit()
   result.plot_pseudo(mode="TE", data_type="rho")
   result.plot_response(site=0)
   result.plot_model()

Before interpretation, compare the section against residual plots. A coherent
anomaly with poor response fit is not reliable geological evidence.

Recommended Run Layout
----------------------

Keep one run directory per scientific experiment.

.. code-block:: text
   :linenos:

   runs/
     profile_a_occam2d_v01/
       occam2d.yml
       provenance.yml
       native/
         OccamDataFile.dat
         Occam2DMesh
         Occam2DModel
         Startup
         ITER12.iter
         RESP12.resp
         LogFile
         occam_stdout.log
         occam_stderr.log
       qc/
         pseudosection.png
         convergence.png
         response_fit_site_014.png
         section_iter12.png
       exports/
         profile_a_iter12.dat
         run_snapshot.zip

Avoid launching new experiments in an old output directory. Occam-style native
files often have simple names, and stale ``.iter`` or ``.resp`` files can be
mistaken for fresh output.

Pre-Run Checklist
-----------------

Before launching:

* load ``OccamConfig`` from the edited template;
* confirm station order and profile offsets;
* inspect selected modes and frequency band;
* confirm apparent-resistivity and phase error floors;
* inspect mesh dimensions, padding, and depth growth;
* inspect model parameter count;
* confirm startup file references the intended data, mesh, and model files;
* confirm the executable path and compiler provenance;
* move old iteration, response, and log files out of the native directory;
* record the exact command and runtime environment.

Post-Run Checklist
------------------

After completion:

* read stdout/stderr logs if pyCSAMT launched the run;
* load the run with ``InversionResult``;
* confirm the selected iteration number and final RMS;
* inspect convergence history;
* inspect per-site and per-frequency misfit;
* compare observed and modeled responses;
* plot the pseudosection and final model together;
* export ``iter2dat`` only after confirming the selected iteration;
* archive the native input and output files with the configuration.

Common Mistakes
---------------

Interpreting smooth gradients too literally
    Occam regularization intentionally smooths structure. A gradient may
    represent a sharper geological contact that is not resolved by the data.

Ignoring station geometry
    The data, mesh, pseudosection, and model section all depend on station
    ordering and offsets. Check them before running.

Using unrealistic error floors
    Too-small floors can force the inversion to fit noise; too-large floors can
    underfit useful signal.

Mixing old and new output files
    If a run fails, old ``.iter`` or ``.resp`` files can remain. Check
    timestamps and log files before loading.

Forgetting startup patching
    Passing ``max_iter`` or ``target_misfit`` to ``OccamRunner.run`` modifies
    ``Startup`` in place. Archive the modified file.

Next Steps
----------

* :doc:`configuration_and_io` for source-of-truth configuration and native file
  archive practice.
* :doc:`choosing_backend` for deciding when Occam2D is the right model
  integration.
* :doc:`../tutorials/prepare_occam2d_inversion` for a practical Occam2D
  workflow.
* :ref:`inversion_concepts` for Occam-style objective functions and
  regularization.
* :doc:`../api/models` for generated API details.
