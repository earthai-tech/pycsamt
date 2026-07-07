.. _models_configuration_and_io:

Configuration And File I/O
==========================

Model integrations are only useful when a run can be repeated, inspected, and
shared. In pyCSAMT, reproducibility is built around two records:

* an editable pyCSAMT configuration file that describes the scientific intent
  of the run;
* the native files required by the selected modelling or inversion engine.

The configuration file answers questions such as "which backend was selected",
"where are the data", "what is the working directory", and "which numerical
controls should be used". The native files answer questions such as "what exact
mesh did the engine see", "what startup file was passed to the binary", "what
response file came back", and "which log records the final misfit".

This page documents that contract. It is intentionally practical: if a run is
sent to a colleague, moved to a cluster, or archived for a report, these are the
files and conventions that should move with it.

The Three-Layer Record
----------------------

A well-organized model run should keep three layers separate.

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Layer
     - Purpose
     - Typical contents
   * - pyCSAMT configuration
     - Stores the editable parameters used to build, run, or load the model.
     - ``InversionConfig``, ``OccamConfig``, ``ModEmConfig``,
       ``Mare2DEMConfig`` templates in ``.py``, ``.json``, ``.yml``, or
       ``.yaml`` format.
   * - Native engine files
     - Stores the files read and written by the external modelling or inversion
       code.
     - Occam2D data, mesh, model, startup, iteration, response, and log files;
       ModEM data, model, covariance, control, response, and log files;
       MARE2DEM ``.emdata``, ``.poly``, ``.resistivity``, ``.settings``,
       ``.EMResp``, and logs.
   * - Derived outputs
     - Stores products created after the engine has completed.
     - Parsed result objects, CSV summaries, figures, GeoJSON/VTK/NPZ exports,
       archive snapshots, and quality-control reports.

The most common mistake is to preserve only the final figure or final model
array. That is not enough. For scientific work, keep the configuration and the
native directory together.

Configuration Template Formats
------------------------------

The shared helper :mod:`pycsamt.models.config_io` writes and reads editable
configuration templates. The same mechanism is used by the high-level
inversion API and by model-specific integrations.

.. code-block:: python
   :linenos:

   from pycsamt.inversion import InversionConfig

   cfg = InversionConfig(
       method="mt",
       dimension="2d",
       backend="occam2d",
       data="data/edis",
       workdir="runs/occam2d_profile_a",
       run_external=False,
   )

   path = cfg.write_template("runs/configs/occam2d_profile_a.yml")
   loaded = InversionConfig.from_file(path)

The model-specific configuration classes follow the same pattern.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import OccamConfig
   from pycsamt.models.modem import ModEmConfig
   from pycsamt.models.mare2dem import Mare2DEMConfig

   OccamConfig.write_template("runs/configs/occam2d_native.py")
   ModEmConfig.write_template("runs/configs/modem_3d.json")
   Mare2DEMConfig.write_template("runs/configs/mare2dem.yml")

   occam = OccamConfig.from_file("runs/configs/occam2d_native.py")
   modem = ModEmConfig.from_file("runs/configs/modem_3d.json")
   mare2dem = Mare2DEMConfig.from_file("runs/configs/mare2dem.yml")

The writer chooses the format from the file suffix. If no suffix is supplied,
``.py`` is used by default.

.. list-table::
   :header-rows: 1
   :widths: 16 34 30 20

   * - Format
     - Best use
     - How it is read
     - Notes
   * - ``.py``
     - Human-edited templates with comments and Python literal values.
     - The file must define a ``CONFIG`` dictionary. pyCSAMT parses the value
       with :func:`ast.literal_eval`.
     - The file is not executed as a script.
   * - ``.json``
     - Machine-generated configs, validation snapshots, and interchange with
       other tools.
     - JSON is loaded as an object. If a top-level ``"config"`` object exists,
       that object is used.
     - Generated JSON includes a ``"_schema"`` documentation block that is
       ignored when loading.
   * - ``.yml`` or ``.yaml``
     - User-facing project files with comments and compact syntax.
     - YAML is read with ``yaml.safe_load``.
     - Requires PyYAML when reading YAML files.

The schema entries behind each template are instances of
:class:`pycsamt.models.config_io.ConfigParameter`. They provide parameter
names, groups, and descriptions so long configuration files remain navigable.
For JSON, the descriptions are written into ``"_schema"`` metadata. For Python
and YAML, they are written as comments.

Strict Loading
--------------

Configuration loading is strict by default. Unknown keys raise an error, which
is useful because most unknown keys are misspellings, obsolete parameters, or
copy-paste mistakes from another backend.

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import OccamConfig

   # Recommended for normal project work.
   cfg = OccamConfig.from_file("runs/configs/occam2d_native.yml")

   # Useful when migrating old files and intentionally ignoring retired keys.
   migrated = OccamConfig.from_file(
       "runs/configs/old_occam2d_native.yml",
       strict=False,
   )

Use ``strict=False`` only as a migration tool. After the file has been cleaned,
write a fresh template and return to strict loading.

Backend-Neutral Versus Native Configuration
-------------------------------------------

pyCSAMT has two related but different configuration levels.

``InversionConfig`` is the backend-neutral entry point. It describes the
workflow: method, dimensionality, backend, input data, working directory,
external execution policy, output paths, and common inversion controls.

Native configuration classes describe one engine family in more detail:
``OccamConfig`` for Occam2D files, ``ModEmConfig`` for ModEM files, and
``Mare2DEMConfig`` for MARE2DEM source/build/run settings.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.inversion import InversionConfig
   from pycsamt.models.modem import ModEmConfig

   native_dir = Path("runs/modem_regional_v01/native")
   workflow = InversionConfig(
       method="mt",
       dimension="3d",
       backend="modem",
       data="data/edis",
       workdir="runs/modem_regional_v01",
       run_external=False,
   )

   native = ModEmConfig(
       data_file="ModEM_Data.dat",
       model_file="ModEM_Model.rho",
       covariance_file="covariance.cov",
       control_file="control.inv",
   )

   workflow.write_template("runs/modem_regional_v01/inversion.yml")
   native.to_template(native_dir / "modem.yml")

Keeping both files is useful. The first explains why this backend was selected;
the second explains exactly how the engine-facing files were named and
organized.

Native Files By Engine
----------------------

Native files are not hidden implementation details. They are the audit trail of
the modelling engine. pyCSAMT classes should therefore preserve paths, parsed
values, warnings, and engine metadata where possible.

.. list-table::
   :header-rows: 1
   :widths: 18 32 28 22

   * - Engine
     - Input responsibility
     - Output responsibility
     - pyCSAMT modules
   * - Occam2D
     - Data file, mesh file, model file, and startup file.
     - Iteration files, response files, log files, misfit summaries.
     - ``pycsamt.models.occam2d.data``, ``mesh``, ``model``, ``startup``,
       ``response``, ``log``, ``runner``.
   * - ModEM
     - Data file, 2-D or 3-D model file, covariance file, and control file.
     - Response files, updated model iterations, run status, log records.
     - ``pycsamt.models.modem.data``, ``model2d``, ``model3d``,
       ``covariance``, ``control``, ``log``, ``runner``.
   * - MARE2DEM
     - ``.emdata``, ``.poly``, ``.resistivity``, ``.settings``, source and
       executable configuration.
     - ``.EMResp``, inversion logs, group-RMS files, updated resistivity
       products.
     - ``pycsamt.models.mare2dem.iotools``, ``builder``, ``runner``,
       ``log``, ``source``, ``validation``.

Occam2D files are compact and profile-oriented. ModEM files are usually larger
and more sensitive to coordinate conventions, dimensionality, and covariance
settings. MARE2DEM projects often include a broader environment record because
source checkout, build settings, and executable provenance matter.

Recommended Directory Layouts
-----------------------------

Use a separate working directory for each scientific run. Avoid reusing the
same directory for exploratory attempts unless the previous outputs have been
archived or removed.

For a high-level backend-neutral run:

.. code-block:: text
   :linenos:

   runs/
     survey_alpha/
       occam2d_profile_a_v01/
         inversion.yml
         inputs/
           edi/
           station_table.csv
         native/
           occam2d.yml
           data.dat
           mesh
           model
           startup
         outputs/
           response.dat
           iter_08
           run.log
         figures/
           apparent_resistivity.png
           section.png
         exports/
           result.npz
           result.csv
           run_snapshot.zip

For a direct native-engine workflow:

.. code-block:: text
   :linenos:

   runs/
     survey_alpha/
       modem_regional_v03/
         modem.yml
         native/
           ModEM_Data.dat
           ModEM_Model.rho
           covariance.cov
           control.inv
           responses/
           models/
           logs/
         qc/
           data_coverage.csv
           rms_history.csv
         exports/
           model.vtk
           stations.geojson

For MARE2DEM, keep source/build material separate from project runs.

.. code-block:: text
   :linenos:

   external/
     mare2dem/
       source/
       build/
       bin/

   runs/
     survey_alpha/
       mare2dem_line_12_v02/
         mare2dem.yml
         native/
           line12.emdata
           line12.poly
           line12.resistivity
           line12.settings
         outputs/
           line12.EMResp
           inversion.log
           group_rms.txt
         exports/
           line12_archive.zip

The directory names should encode the experiment, not only the engine. Names
such as ``line_12_static_shift_corrected_v04`` are more useful than
``run_new`` when comparing alternatives six months later.

Builder, Runner, Loader Contract
--------------------------------

Model integrations should keep three responsibilities distinct.

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Component
     - What it does
     - What it should avoid
   * - Builder
     - Creates or updates native input files from pyCSAMT objects and config
       values.
     - Silently launching a long external inversion.
   * - Runner
     - Launches the external executable with explicit paths, working directory,
       process settings, and logging.
     - Guessing missing scientific parameters.
   * - Loader
     - Reads completed native outputs into pyCSAMT result objects.
     - Modifying native outputs in place.

This contract makes cluster workflows easier. A user can build native files on
a laptop, run the external engine on a cluster, then load the completed outputs
back into pyCSAMT without changing the scientific configuration.

.. code-block:: python
   :linenos:

   from pycsamt.inversion import InversionConfig, run_inversion

   cfg = InversionConfig.from_file("runs/occam2d_profile_a_v01/inversion.yml")

   # Build and validate the run directory, but do not launch the external code.
   cfg.run_external = False
   result = run_inversion(cfg)

   # The native directory can now be transferred to another machine if needed.

For external binaries, prefer explicit execution over implicit execution. The
configuration should say where the binary is, which working directory is used,
and whether pyCSAMT should launch the process or only prepare/load files.

Validation Before A Run
-----------------------

Validation should happen before an external binary is launched. It is cheaper
to catch a bad coordinate convention or missing startup file before a long
inversion starts.

.. list-table::
   :header-rows: 1
   :widths: 28 42 30

   * - Check
     - Why it matters
     - Example failure
   * - Required files exist
     - External engines often fail late or write cryptic logs.
     - ``startup`` points to a mesh file that was never written.
   * - Dimensionality matches data
     - A 2-D profile, 3-D grid, and 1-D sounding have different assumptions.
     - ModEM 3-D config is used with a profile-only station set.
   * - Units are explicit
     - Resistivity, conductivity, distance, depth, frequency, and period must
       remain consistent.
     - Mesh spacing is interpreted as metres when the source table was in
       kilometres.
   * - Coordinate frame is known
     - Native engines may expect local profile coordinates rather than
       longitude/latitude.
     - Station ordering is correct but profile offsets are reversed.
   * - Old outputs are isolated
     - Stale response files can be mistaken for new results.
     - A loader reads yesterday's ``.EMResp`` after today's run failed.
   * - Executable provenance is recorded
     - Results can vary with engine version, compilation flags, and MPI setup.
     - A report cannot identify which MARE2DEM binary produced the final model.

When a validation helper is available, use it before running the engine.
Validation modules exist for the engine-specific packages, including Occam2D,
ModEM, and MARE2DEM.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.models.occam2d import OccamConfig

   native_dir = Path("runs/occam2d_profile_a_v01/native")
   cfg = OccamConfig.from_file("runs/occam2d_profile_a_v01/native/occam2d.yml")

   required = [
       cfg.data_file,
       cfg.mesh_file,
       cfg.model_file,
       cfg.startup_file,
   ]
   missing = [
       name for name in required
       if not (native_dir / name).exists()
   ]
   if missing:
       raise FileNotFoundError(f"Missing Occam2D native files: {missing}")

The exact field names depend on the native configuration class, but the pattern
is the same: load the config, resolve paths from the working directory, and
fail before launching the engine.

Provenance To Keep With Every Run
---------------------------------

A model directory should make the run understandable without relying on memory
or notebook state. At minimum, keep:

* the edited pyCSAMT configuration file;
* the engine-specific native configuration file, when one was used;
* the native input files supplied to the external code;
* the native output files produced by the external code;
* logs from pyCSAMT and from the external engine;
* the pyCSAMT version or source revision used to create the run;
* the external executable path and version information when available;
* the command line, MPI settings, or scheduler script used to run the engine;
* a short note describing the scientific purpose of the run.

A small provenance file is often enough.

.. code-block:: yaml
   :linenos:

   project: survey_alpha
   run_id: occam2d_profile_a_v01
   created_by: pycsamt
   pycsamt_version: 2.x
   configuration:
     workflow: inversion.yml
     native: native/occam2d.yml
   engine:
     name: Occam2D
     executable: /opt/occam2d/bin/occam2d
     launched_by_pycsamt: false
   data:
     source: inputs/edi
     station_table: inputs/station_table.csv
   notes: >
     Initial 2-D profile inversion after static-shift review.

Prefer plain text formats for provenance. They survive version control,
archives, and long-term project storage better than notebook-only metadata.

Archiving Results
-----------------

The inversion export helpers can include native files in a portable archive
when the result object carries native-file metadata.

.. code-block:: python
   :linenos:

   from pycsamt.inversion.export import to_archive

   archive = to_archive(
       result,
       "runs/occam2d_profile_a_v01/exports/run_snapshot.zip",
       include_native=True,
   )

The archive should not replace the run directory while work is active. Treat it
as a snapshot for transfer, publication support, or long-term storage. During
active interpretation, keep the full directory tree available so native files
can be inspected directly.

Common Mistakes
---------------

Do not mix configuration levels
    ``InversionConfig`` chooses the workflow and backend. ``OccamConfig``,
    ``ModEmConfig``, and ``Mare2DEMConfig`` describe native engine details.
    Keeping those concerns separate makes documentation, testing, and
    migration easier.

Do not edit generated native outputs by hand
    If a response, iteration, or log file is edited manually, the run is no
    longer a clean record of the external engine. Write a derived file instead.

Do not run new experiments in an old output directory
    External codes often reuse simple file names. A stale response file can
    make a failed run look successful.

Do not ignore unknown configuration keys
    Unknown keys should usually fail. Use ``strict=False`` only for controlled
    migration of old templates.

Do not archive only figures
    Figures are interpretation products. They are not enough to reproduce the
    inversion.

Pre-Run Checklist
-----------------

Before launching or submitting a run:

* generate a fresh configuration template;
* edit the template instead of changing values only inside a notebook;
* load the edited template with strict validation;
* confirm data paths and working directory paths;
* inspect station coordinates, profile direction, and units;
* confirm the selected dimensionality and backend;
* confirm the external executable and runtime policy;
* write or refresh native input files;
* move old outputs away from the run directory;
* record provenance for the intended experiment.

Post-Run Checklist
------------------

After the external code finishes:

* inspect the engine log before plotting results;
* confirm that output timestamps match the intended run;
* parse responses and models with the engine-specific loader;
* preserve native output paths inside the result metadata where possible;
* export compact summaries for downstream analysis;
* create an archive snapshot if the run will be shared;
* write a short interpretation note while the run context is still fresh.

Next Steps
----------

* :doc:`choosing_backend` explains how to decide between backend integrations.
* :doc:`occam2d` documents Occam2D profile-oriented files and workflows.
* :doc:`modem` documents ModEM 2-D/3-D files and execution conventions.
* :doc:`mare2dem` documents MARE2DEM source, files, geometry, and logs.
