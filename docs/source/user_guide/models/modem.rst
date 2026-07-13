.. _models_modem:

ModEM
=====

:mod:`pycsamt.models.modem` is the ModEM integration layer in pyCSAMT v2. It
does not reimplement the ModEM solver. Instead, it gives Python objects for
preparing ModEM input files, checking native file types, launching an external
ModEM executable, loading completed runs, and plotting the results.

ModEM projects are file-oriented. A reproducible run is normally defined by an
observed data file, an initial model file, an inversion-control file, and, for
3-D inversions, a covariance file. The solver then produces iteration models,
predicted responses, and logs. pyCSAMT keeps those files explicit so a project
can move between Python, command-line ModEM, and archived native run folders
without losing provenance.

When To Use ModEM
-----------------

Use ModEM when the survey, interpretation target, or legacy project requires a
native ModEM workflow.

Good ModEM candidates include:

* MT or AMT surveys with stations distributed over an area, not only one line;
* geological settings where off-profile structure is expected to affect the
  response;
* projects that need 3-D inversion, covariance masks, or region-specific
  smoothing;
* existing ModEM projects that pyCSAMT should validate, plot, or archive;
* forward-response checks against a known ModEM model;
* backend-neutral pyCSAMT inversion workflows where the numerical solver is an
  external ModEM executable.

For simple 2-D CSAMT profile inversions, :doc:`occam2d` may be more direct.
For large 3-D MT or AMT modelling, ModEM is usually the more appropriate
native backend.

Dimensionality
--------------

ModEM supports both 2-D and 3-D modes in pyCSAMT. The mode controls the model
class, file set, executable name, and whether a covariance file is expected.

.. list-table::
   :header-rows: 1
   :widths: 18 26 26 30

   * - Mode
     - Main model object
     - Native model file
     - Typical use
   * - ``"2d"``
     - :class:`pycsamt.models.modem.ModEmModel2D`
     - ``.rho``
     - Profile inversion, TE/TM or selected impedance components.
   * - ``"3d"``
     - :class:`pycsamt.models.modem.ModEmModel3D`
     - ``.ws`` for starting models; ``.rho`` or ``.ws`` products
     - Area surveys, full impedance tensors, and 3-D covariance control.

.. list-table::
   :header-rows: 1
   :widths: 22 28 25 25

   * - Mode
     - Default executable
     - Covariance file
     - Default builder output
   * - ``"2d"``
     - ``Mod2DMT``
     - Not created by :class:`~pycsamt.models.modem.InputBuilder`
     - ``data``, ``model``, ``control``
   * - ``"3d"``
     - ``Mod3DMT``
     - Created from the 3-D model grid
     - ``data``, ``model``, ``covariance``, ``control``

Package Map
-----------

The public ModEM API is grouped around native file roles.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Object or function group
     - Role
   * - :class:`pycsamt.models.modem.ModEmConfig`
     - Stores data settings, dimensionality, grid controls, covariance
       controls, inversion-control values, file names, executable names, and
       MPI settings.
   * - :class:`pycsamt.models.modem.InputBuilder`
     - Builds a consistent ModEM input directory from EDI-like station objects
       or a populated :class:`~pycsamt.models.modem.ModEmData`.
   * - :class:`pycsamt.models.modem.ModEmData`
     - Reads, writes, and builds ModEM observed or predicted data files.
   * - :class:`pycsamt.models.modem.ModEmModel2D`
     - Reads, writes, and creates 2-D ModEM resistivity models.
   * - :class:`pycsamt.models.modem.ModEmModel3D`
     - Reads, writes, and creates 3-D ModEM resistivity models.
   * - :class:`pycsamt.models.modem.ModEmCovariance`
     - Represents 3-D smoothing coefficients, region masks, and smoothing
       exceptions.
   * - :class:`pycsamt.models.modem.ModEmControl`
     - Reads and writes the ``.inv`` nonlinear inversion-control file.
   * - :class:`pycsamt.models.modem.ModEmRunner`
     - Assembles and runs external ModEM inversion or forward commands.
   * - :class:`pycsamt.models.modem.InversionResult`
     - Scans a run directory and loads models, responses, covariance, control,
       and logs.
   * - :class:`pycsamt.models.modem.ModEmLog`
     - Parses iteration number, RMS, objective value, model norm, lambda, and
       alpha from a ModEM log.
   * - ``PlotMisfit``, ``PlotModel2D``, ``PlotModel3D``,
       ``PlotResponse``, ``PlotPseudo``
     - Matplotlib diagnostics for convergence, model inspection, station
       responses, and pseudo-sections.
   * - ``detect_file_type`` and ``is_*`` validators
     - Identify ModEM data, model, covariance, control, and log files before
       reading or routing them.
   * - ``read_z3d_old``, ``write_z3d_old``, ``convert_z3d``,
       ``read_mackie3d``, ``write_meshtools3d``, and related utilities
     - Convert older impedance, Mackie, MeshTools, and interpolation-oriented
       file formats.

Configuration
-------------

Most workflows start with :class:`~pycsamt.models.modem.ModEmConfig`. The same
configuration object is passed to the builder, runner, control file, model
factory, covariance factory, and result loader.

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import ModEmConfig

   cfg = ModEmConfig(
       mode="3d",
       component_type="Full_Impedance",
       error_floor_z=0.05,
       initial_rho=100.0,
       nx=24,
       ny=24,
       nz=36,
       cell_size_h=500.0,
       cell_size_v_top=10.0,
       depth_scale=1.18,
       n_padding_xy=8,
       smooth_x=0.2,
       smooth_y=0.2,
       smooth_z=0.1,
       n_smooth_iter=2,
       max_iterations=80,
       target_rms=1.05,
       binary_3d="Mod3DMT",
       use_mpi=True,
       n_procs=16,
   )

   cfg.write_template("modem_config.ini")

The template can be edited and loaded again:

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import ModEmConfig

   cfg = ModEmConfig.from_file("modem_config.ini", strict=True)
   print(cfg.mode, cfg.binary_name)

Important settings are grouped as follows.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Setting group
     - Main fields
   * - Dimensionality
     - ``mode`` controls ``"2d"`` or ``"3d"`` behavior. ``cfg.is_3d`` and
       ``cfg.binary_name`` are convenience properties.
   * - Data block
     - ``component_type``, ``sign_convention``, ``units``,
       ``error_floor_z``, ``error_floor_z_floor``, ``freq_min``,
       ``freq_max``.
   * - 2-D grid
     - ``nx_2d``, ``nz_2d``, ``n_airlayers_2d``, ``cell_size_h_2d``,
       ``cell_size_v_top_2d``, ``depth_scale_2d``, ``n_padding_x_2d``.
   * - 3-D grid
     - ``nx``, ``ny``, ``nz``, ``n_airlayers``, ``cell_size_h``,
       ``cell_size_v_top``, ``depth_scale``, ``n_padding_xy``.
   * - Regularization
     - ``smooth_x``, ``smooth_y``, ``smooth_z``, ``n_smooth_iter``.
   * - Nonlinear inversion
     - ``max_iterations``, ``target_rms``, ``initial_lambda``,
       ``lambda_divisor``, ``initial_alpha``, ``rms_diff_tol``,
       ``lambda_exit``.
   * - File names
     - ``data_file``, ``model_file``, ``covariance_file``, ``control_file``,
       ``log_file``, ``output_stem``.
   * - Execution
     - ``binary_2d``, ``binary_3d``, ``use_mpi``, ``n_procs``,
       ``mpi_command``.

Native Files
------------

A ModEM run folder should be understandable without Python. pyCSAMT therefore
writes and reads the same native file roles that the executable uses.

.. list-table::
   :header-rows: 1
   :widths: 22 22 56

   * - File role
     - Python object
     - Notes
   * - Observed data
     - ``ModEmData``
     - ASCII data blocks grouped by component type. Built from EDI-like
       station objects or read from existing ModEM files.
   * - 2-D model
     - ``ModEmModel2D``
     - ``.rho`` resistivity model with horizontal and vertical cell widths.
   * - 3-D model
     - ``ModEmModel3D``
     - ``.ws`` starting model or 3-D ``.rho``/``.ws`` products loaded from a
       run directory.
   * - Covariance
     - ``ModEmCovariance``
     - 3-D earth-only smoothing file. It stores smoothing arrays, exceptions,
       and integer masks.
   * - Control
     - ``ModEmControl``
     - ``.inv`` file with output stem, lambda controls, target RMS, and
       iteration limits.
   * - Log
     - ``ModEmLog``
     - Parsed convergence history from ModEM text logs.
   * - Predicted data
     - ``ModEmData``
     - Loaded by :class:`~pycsamt.models.modem.InversionResult` as
       ``data_pred`` when a predicted response file is present.

Validate files before routing them into a workflow:

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.models.modem import detect_file_type

   for path in Path("runs/modem_3d").iterdir():
       kind = detect_file_type(path)
       if kind.value != "unknown":
           print(path.name, kind.value)

Data Files And Components
-------------------------

:class:`~pycsamt.models.modem.ModEmData` stores observed or predicted response
rows. pyCSAMT can read existing ModEM data files or build new ones from
station objects.

When building from EDI-like objects, each station should provide:

* a station name;
* station coordinates or latitude/longitude/elevation metadata;
* positive frequencies;
* complex impedance values shaped like ``(n_frequency, 2, 2)``;
* impedance errors or uncertainty estimates.

The configuration selects which response components are written. Supported
``component_type`` values include:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - ``component_type``
     - Components written
   * - ``"Full_Impedance"``
     - ``ZXX``, ``ZXY``, ``ZYX``, ``ZYY``
   * - ``"Off_Diagonal_Impedance"``
     - ``ZXY``, ``ZYX``
   * - ``"TE_Impedance"``
     - ``TE``
   * - ``"TM_Impedance"``
     - ``TM``
   * - ``"Determinant_Impedance"``
     - ``ZDet``
   * - ``"Full_Vertical_Components"``
     - ``HZX``, ``HZY``
   * - ``"Phase_Tensor"``
     - ``PTxx``, ``PTxy``, ``PTyx``, ``PTyy``

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import ModEmConfig, ModEmData

   cfg = ModEmConfig(
       mode="3d",
       component_type="Off_Diagonal_Impedance",
       error_floor_z=0.05,
       freq_min=1e-3,
       freq_max=1e3,
   )

   data = ModEmData.from_edi(stations, config=cfg)
   data.write("runs/modem_3d/ModEMData.dat")

   print(data.n_sites)
   print(data.n_periods)
   print(data.component_types)

.. important::

   ModEM data files depend on consistent local coordinates. The pyCSAMT writer
   stores local ``x``, ``y``, and ``z`` coordinates used by the model builder.
   Confirm projection, station order, and units before inversion. Do not mix
   metres, kilometres, and geographic degrees inside the same run folder.

Build A 3-D Input Set
---------------------

:class:`~pycsamt.models.modem.InputBuilder` creates the standard 3-D input
set: observed data, starting model, covariance, and control file.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.models.modem import InputBuilder, ModEmConfig

   workdir = Path("runs/modem_3d/native")
   cfg = ModEmConfig(
       mode="3d",
       component_type="Full_Impedance",
       initial_rho=100.0,
       nx=28,
       ny=28,
       nz=38,
       n_airlayers=5,
       cell_size_h=500.0,
       cell_size_v_top=10.0,
       depth_scale=1.18,
       n_padding_xy=8,
       smooth_x=0.2,
       smooth_y=0.2,
       smooth_z=0.1,
       n_smooth_iter=2,
   )

   builder = InputBuilder(config=cfg)
   files = builder.build(
       stations,
       workdir=workdir,
       data_filename=cfg.data_file,
       model_filename="m0.ws",
       cov_filename=cfg.covariance_file,
       ctrl_filename=cfg.control_file,
   )

   for role, path in files.items():
       print(role, path)

The returned mapping contains ``data``, ``model``, ``covariance``, and
``control``. The builder also keeps the populated objects on
``builder.data``, ``builder.model``, ``builder.covariance``, and
``builder.control`` for inspection.

Build A 2-D Input Set
---------------------

The 2-D builder path writes data, a 2-D half-space model, and a control file.
It does not create a 3-D covariance file.

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import InputBuilder, ModEmConfig

   cfg = ModEmConfig(
       mode="2d",
       component_type="TE_Impedance",
       initial_rho=100.0,
       nx_2d=120,
       nz_2d=60,
       n_airlayers_2d=5,
       cell_size_h_2d=100.0,
       cell_size_v_top_2d=10.0,
       depth_scale_2d=1.18,
       n_padding_x_2d=8,
       max_iterations=60,
   )

   files = InputBuilder(config=cfg).build(
       profile_stations,
       workdir="runs/modem_2d/native",
       data_filename="ModEMData.dat",
       model_filename="m0.rho",
       ctrl_filename="ModEM.inv",
   )

   assert "covariance" not in files

Build From Existing Data
------------------------

Use ``build_from_data`` when the observed data object has already been read,
filtered, or edited. In this mode the builder creates the model, covariance
when needed, and control file. The caller is responsible for writing the data
file itself.

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import InputBuilder, ModEmConfig, ModEmData

   cfg = ModEmConfig(mode="3d")
   data = ModEmData.read("survey/ModEMData.dat")
   data.write("runs/modem_3d/native/ModEMData.dat")

   files = InputBuilder(config=cfg).build_from_data(
       data,
       workdir="runs/modem_3d/native",
       model_filename="m0.ws",
       cov_filename="ModEM.cov",
       ctrl_filename="ModEM.inv",
   )

Models
------

The model objects store cell widths and resistivity values. Internally,
resistivity is stored in natural-log units because ModEM solves for log
resistivity. The ``rho_linear`` property returns resistivity in linear
``ohm m`` units for plotting and interpretation.

Create and inspect a 3-D starting model:

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import ModEmData, ModEmModel3D, ModEmConfig

   cfg = ModEmConfig(mode="3d", initial_rho=100.0)
   data = ModEmData.read("runs/modem_3d/native/ModEMData.dat")

   model = ModEmModel3D.halfspace(data, config=cfg)
   model.write("runs/modem_3d/native/m0.ws")

   print(model.shape)
   print(model.n_air)
   print(model.rho_linear.min(), model.rho_linear.max())

Create and inspect a 2-D starting model:

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import ModEmData, ModEmModel2D, ModEmConfig

   cfg = ModEmConfig(mode="2d", initial_rho=100.0)
   data = ModEmData.read("runs/modem_2d/native/ModEMData.dat")

   model = ModEmModel2D.halfspace(data, config=cfg)
   model.write("runs/modem_2d/native/m0.rho")

   print(model.nx, model.nz)
   print(model.x_nodes[-1], model.z_nodes[-1])

The half-space factories are intentionally conservative. They are useful for a
first run, but a production inversion usually deserves a deliberate mesh
review: station spacing, padding, first-layer thickness, air layers, expected
skin depths, and target depth all matter.

Covariance
----------

The covariance file is central to 3-D ModEM interpretation. It controls the
model regularization term through smoothing coefficients and integer masks.
pyCSAMT creates a uniform active earth region by default, then lets advanced
users edit masks and exceptions before writing the file.

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import ModEmCovariance

   cov = ModEmCovariance.from_model(model, config=cfg)
   cov.exceptions.append((1, 2, 0.0))  # turn off smoothing across two regions
   cov.write("runs/modem_3d/native/ModEM.cov")

Mask values follow the ModEM convention:

* ``0`` is reserved for air;
* ``9`` is reserved for ocean;
* ``1`` through ``8`` are user-defined earth regions.

``ModEmCovariance.from_model`` excludes air layers from the covariance grid.
The number of covariance layers is therefore ``model.nz - model.n_air``.

Control Files
-------------

:class:`~pycsamt.models.modem.ModEmControl` stores the nonlinear inversion
settings written to the ``.inv`` file.

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import ModEmConfig, ModEmControl

   cfg = ModEmConfig(
       output_stem="ModEM_out",
       initial_lambda=10.0,
       lambda_divisor=100.0,
       initial_alpha=10.0,
       rms_diff_tol=5e-4,
       target_rms=1.05,
       lambda_exit=1e-4,
       max_iterations=100,
   )

   control = ModEmControl.from_config(cfg)
   control.write("runs/modem_3d/native/ModEM.inv")

The control file does not replace data-quality assessment. A low target RMS is
meaningful only when uncertainty floors, component selection, and bad-period
masking are realistic.

Run ModEM
---------

:class:`~pycsamt.models.modem.ModEmRunner` assembles the external command and
can execute it with :func:`subprocess.run`. The executable is resolved from
``PATH`` or from the run directory and local ``_source/2D`` or ``_source/3D``
subdirectories.

Always inspect the command first:

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import ModEmConfig, ModEmRunner

   cfg = ModEmConfig(
       mode="3d",
       binary_3d="Mod3DMT",
       use_mpi=True,
       n_procs=16,
       mpi_command="mpirun",
   )

   runner = ModEmRunner("runs/modem_3d/native", config=cfg)
   command = runner.command(
       "m0.ws",
       "ModEMData.dat",
       "ModEM.inv",
       covariance="ModEM.cov",
   )
   print(command)

Run the inversion only after file paths, executable names, and MPI settings
are correct:

.. code-block:: python
   :linenos:

   result = runner.run(
       "m0.ws",
       "ModEMData.dat",
       "ModEM.inv",
       covariance="ModEM.cov",
       timeout=24 * 3600,
       load_result=True,
   )

For a forward response check, call ``run_forward``:

.. code-block:: python
   :linenos:

   runner.run_forward(
       "m0.ws",
       "ModEMData.dat",
       timeout=3600,
       load_result=False,
   )

.. warning::

   The runner is a subprocess wrapper around an external executable. It cannot
   make a physically poor mesh, incorrect component selection, or inconsistent
   coordinate system valid. Treat the generated command as a reproducibility
   aid, not as a scientific approval stamp.

Load Results
------------

:class:`~pycsamt.models.modem.InversionResult` scans a completed run directory
and loads what it finds: logs, controls, covariance, data files, and iteration
models.

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import InversionResult

   result = InversionResult("runs/modem_3d/native", config=cfg)

   print(result.mode)
   print(result.n_iter)
   print(result.final_rms)
   print(result.best_rms)
   print(result.iteration_numbers)

   final_model = result.model_final
   observed = result.data_obs
   predicted = result.data_pred

The result loader recognizes common ModEM output naming patterns, including
numbered ``Modular_NLCG`` products. The lowest numbered response can be used
as an observed-data fallback and the highest numbered response can be used as
the predicted response when explicit filenames are not available.

Log Diagnostics
---------------

Use :class:`~pycsamt.models.modem.ModEmLog` when the convergence history is
the primary diagnostic.

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import ModEmLog

   log = ModEmLog.read("runs/modem_3d/native/Modular_NLCG.log")

   print(log.n_iter)
   print(log.final_rms)
   print(log.best_iter)
   print(log.rms)
   print(log.lagrange)

Review more than the final RMS. Sudden RMS stalls, unstable lambda changes,
or a best iteration far earlier than the final iteration can indicate
overfitting, inconsistent errors, or a regularization setting that should be
revisited.

Plotting
--------

The ModEM plotters operate on
:class:`~pycsamt.models.modem.InversionResult` objects and return Matplotlib
figures.

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import (
       InversionResult,
       PlotMisfit,
       PlotModel3D,
       PlotPseudo,
       PlotResponse,
   )

   result = InversionResult("runs/modem_3d/native")

   fig = PlotMisfit(result=result).plot()
   fig.savefig("runs/modem_3d/figures/rms.png", dpi=200)

   fig = PlotModel3D(
       result=result,
       depths=[500, 1000, 2000, 4000],
       rho_min=1.0,
       rho_max=1000.0,
   ).plot()
   fig.savefig("runs/modem_3d/figures/model_slices.png", dpi=200)

   fig = PlotResponse(
       result=result,
       stations=["S001", "S002"],
       max_stations=2,
   ).plot()
   fig.savefig("runs/modem_3d/figures/responses.png", dpi=200)

   fig = PlotPseudo(result=result, component="ZXY").plot()
   fig.savefig("runs/modem_3d/figures/pseudo_zxy.png", dpi=200)

For 2-D results, use :class:`pycsamt.models.modem.PlotModel2D`:

.. code-block:: python
   :linenos:

   from pycsamt.models.modem import InversionResult, PlotModel2D

   result = InversionResult("runs/modem_2d/native")
   fig = PlotModel2D(result=result, depth_max=5000.0).plot()
   fig.savefig("runs/modem_2d/figures/model_section.png", dpi=200)

Conversion And Utility Tools
----------------------------

The ModEM package also exposes utility functions for existing projects and
format conversion.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Utility group
     - Examples
   * - Impedance files
     - ``ZBlock``, ``ImpedanceFile``, ``read_z3d_old``,
       ``write_z3d_old``, ``read_z2d_old``, ``write_z2d_old``,
       ``write_z3d_list``, ``write_z2d_list``, ``convert_z3d``,
       ``convert_z2d``.
   * - Mackie formats
     - ``read_mackie2d``, ``write_mackie2d``, ``read_mackie3d``,
       ``write_mackie3d``.
   * - MeshTools export
     - ``write_meshtools3d``.
   * - Interpolation
     - ``interp_model3d``, ``interp_z3d``.
   * - Units and transforms
     - ``skin_depth``, ``imp_units_factor``, ``loge_to_log10``,
       ``log10_to_loge``, ``loge_to_linear``, ``linear_to_loge``.

These helpers are most useful when a project arrives with older ModEM,
Mackie, or impedance-list files and pyCSAMT is being used as a bridge into a
clean v2 run directory.

Backend-Neutral Workflows
-------------------------

The native objects above are the most explicit way to work with ModEM. pyCSAMT
also exposes ModEM through the backend-neutral inversion interface.

.. code-block:: python
   :linenos:

   from pycsamt.inversion.config import InversionConfig
   from pycsamt.inversion.backends.modem import ModEMBackend

   inv_cfg = InversionConfig(
       method="mt",
       dimension="3d",
       backend="modem",
       workdir="runs/modem_3d/native",
       run_external=False,
       backend_options={
           "config": {
               "component_type": "Full_Impedance",
               "initial_rho": 100.0,
               "binary_3d": "Mod3DMT",
               "use_mpi": True,
               "n_procs": 16,
           },
           "files": {
               "data": "ModEMData.dat",
               "model": "m0.ws",
               "control": "ModEM.inv",
               "covariance": "ModEM.cov",
           },
           "runner": {
               "timeout": 24 * 3600,
               "load_result": True,
           },
       },
   )

   result = ModEMBackend(inv_cfg).run()
   print(result.status)
   print(result.metadata["command"])

With ``run_external=False``, the backend prepares or validates the run folder
and reports the command that would be executed. Set ``run_external=True`` only
when the external ModEM executable is installed and the run directory has been
reviewed.

Recommended Run Layout
----------------------

A stable ModEM project layout separates native inputs, figures, and notes:

.. code-block:: text
   :linenos:

   runs/
     modem_3d/
       README.md
       config/
         modem_config.ini
       native/
         ModEMData.dat
         m0.ws
         ModEM.cov
         ModEM.inv
         Modular_NLCG.log
         Modular_NLCG_001.rho
         Modular_NLCG_001.dat
       figures/
         rms.png
         model_slices.png
         responses.png
       exports/
         final_model.vtk

Keep hand-edited files under version control when possible. Large model
outputs and response grids can be archived separately if repository size is a
concern.

Pre-Run Checklist
-----------------

Before starting a ModEM inversion, verify:

* station coordinates use one local coordinate system and metre units;
* periods and frequencies are in the intended range;
* impedance sign convention and units match the executable expectation;
* component selection is appropriate for the survey dimensionality;
* error floors are realistic and bad periods have been removed or masked;
* horizontal cell size reflects station spacing and target resolution;
* vertical first-layer thickness and depth growth resolve shallow structure;
* padding extends far enough from the station footprint;
* air layers and ocean or inactive masks are correct;
* covariance smoothing values and exceptions are scientifically justified;
* the dry-run command points to the intended executable and files;
* MPI process count matches the machine and executable build.

Post-Run Checklist
------------------

After a run finishes, review:

* RMS history and best iteration, not only final RMS;
* response fits by station, period, and component;
* residual concentration around specific stations or period bands;
* model updates near boundaries, air layers, and padding cells;
* sensitivity of interpretation to error floors and covariance settings;
* consistency between final model features and known geology;
* reproducibility of the run directory, configuration, command, and code
  version.

Common Mistakes
---------------

``The command is ready but the executable cannot be found.``
   Check ``binary_2d`` or ``binary_3d`` in :class:`ModEmConfig`. The runner
   searches ``PATH``, the working directory, and local ``_source`` folders.

``The 3-D run is missing a covariance file.``
   Build the input set with ``mode="3d"`` or create
   :class:`ModEmCovariance` from the 3-D model and pass it to the runner.

``The model loads but plotted station positions look wrong.``
   Recheck coordinate origin, projection, and units before trusting any
   response or model diagnostics.

``The RMS is low but the model has unrealistic detail.``
   Inspect uncertainty floors, removed periods, smoothing values, covariance
   masks, and response residuals. A visually detailed model is not necessarily
   better constrained.

``The result loader did not find a predicted response.``
   Confirm the output stem in the control file and inspect the run directory
   for numbered ``Modular_NLCG`` response files.

Next Steps
----------

* :doc:`choosing_backend` explains when ModEM is preferable to other model
  backends.
* :doc:`configuration_and_io` gives the shared model-backend configuration
  and file-layout policy.
* :doc:`occam2d` documents the 2-D Occam-style alternative.
* :ref:`inversion_concepts` introduces misfit, regularization, and inversion
  diagnostics.
* :doc:`../../api/models` links to generated API pages for the ModEM objects.
