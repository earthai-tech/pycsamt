.. _models_choosing_backend:

Choosing A Model Backend
========================

Choosing :term:`Occam2D`, :term:`ModEM`, :term:`MARE2DEM`, SimPEG, pyGIMLi,
or the built-in pyCSAMT :term:`inversion backend` is a scientific decision
first and a software decision second. A :term:`backend` is not only a solver.
It encodes assumptions about :term:`dimensionality`, source physics,
:term:`mesh` geometry, :term:`regularization`, file formats, and how much
native-engine control the user must manage.

This page helps you decide which path to use before you start preparing
files. It complements :ref:`inversion_concepts`, which explains the
:term:`objective function`, :term:`regularization`, :term:`RMS misfit`, and
:term:`non-uniqueness`.

Two Different Choices
---------------------

pyCSAMT has two related but different model choices:

* **inversion backend** - selected by
  :class:`pycsamt.inversion.config.InversionConfig` through ``backend``;
* **model integration** - a direct interface under :mod:`pycsamt.models`
  that prepares native files, runs external programs, and loads results.

The backend-neutral inversion API currently accepts:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.inversion.backends import available_backends
   >>> available_backends()
   ['builtin', 'modem', 'occam2d', 'pygimli', 'simpeg']

The :mod:`pycsamt.models` documentation covers:

* :mod:`pycsamt.models.occam2d`;
* :mod:`pycsamt.models.modem`;
* :mod:`pycsamt.models.mare2dem`;
* shared configuration and native file I/O.

At the time of writing, :term:`MARE2DEM` is a rich :term:`model integration`,
not a value for ``InversionConfig(backend=...)``. Use the MARE2DEM package
directly for source management, input building, execution, and output loading.

Quick Recommendation
--------------------

.. list-table::
   :header-rows: 1
   :widths: 18 26 28 28

   * - Situation
     - Recommended path
     - Why
     - Main caution
   * - First quality check or smoke test
     - ``backend="builtin"``
     - Runs without external binaries and validates the broad data flow.
     - Use for checks, not final production interpretation.
   * - Single station or layered setting
     - ``builtin``, ``pygimli``, or ``simpeg``
     - 1-D models are fast and useful for sounding-scale interpretation.
     - Lateral structure can be mapped into false layers.
   * - Profile-scale MT/AMT/CSAMT
     - Occam2D
     - Smooth 2-D sections are robust for profile interpretation.
     - Requires approximately 2-D geology and correct mode/strike handling.
   * - Area survey or strong 3-D structure
     - ModEM
     - Supports 3-D model files, covariance controls, and ModEM-native runs.
     - Mesh, covariance, coordinates, and error floors require care.
   * - MT/CSEM with 2.5-D finite elements
     - MARE2DEM
     - Handles MARE2DEM-native ``.emdata``, ``.poly``, and adaptive FEM
       workflows.
     - Requires external source/binary management and MPI/compiler setup.
   * - Python-native research workflow
     - SimPEG or pyGIMLi
     - Useful when optional scientific inversion libraries are installed.
     - API compatibility and dependency versions matter.

Backend Matrix
--------------

.. list-table::
   :header-rows: 1
   :widths: 17 14 17 20 18 14

   * - Path
     - Dimensionality
     - Methods
     - Main outputs
     - Execution style
     - pyCSAMT route
   * - Built-in
     - 1-D, stitched 2-D
     - MT, AMT, CSAMT, TDEM
     - Common :class:`pycsamt.inversion.results.InversionResult`
     - Pure Python local execution
     - ``InversionConfig``
   * - Occam2D
     - 2-D
     - MT, AMT, CSAMT-style profiles
     - Occam data, mesh, model, response, log, common result adapter
     - External executable or load-only mode
     - ``InversionConfig`` or ``models.occam2d``
   * - ModEM
     - 2-D, 3-D
     - Primarily MT
     - Data, model, covariance, control, response, log
     - External executable or file parsing
     - ``InversionConfig`` or ``models.modem``
   * - MARE2DEM
     - 2.5-D
     - MT, CSEM, DC-style native files
     - ``.emdata``, ``.poly``, ``.resistivity``, ``.settings``, logs
     - External source/binary and MPI
     - ``models.mare2dem``
   * - SimPEG
     - Backend dependent
     - MT/AMT/CSAMT/TDEM paths in pyCSAMT
     - Common inversion result
     - Optional Python dependency
     - ``InversionConfig``
   * - pyGIMLi
     - 1-D, stitched 2-D
     - MT/AMT/CSAMT/TDEM paths in pyCSAMT
     - Common inversion result
     - Optional Python dependency
     - ``InversionConfig``

Decision Axis 1: Dimensionality
-------------------------------

Dimensionality is the most important choice because it defines what the
forward problem is allowed to vary. In a :term:`1D` model, electrical
resistivity is only a function of depth,
:math:`\rho(\mathbf{x})=\rho(z)`. In a :term:`2-D` profile model it varies
across the section and with depth,
:math:`\rho(\mathbf{x})=\rho(x,z)`, while being assumed constant along
:term:`geoelectric strike`. In a :term:`3D` model it can vary in all spatial
directions, :math:`\rho(\mathbf{x})=\rho(x,y,z)`. The more freedom the model
has, the more data support, computation, and prior control it needs.

Use a **1-D** path when:

* the target is approximately layered;
* stations are sparse;
* the goal is sounding-scale reconnaissance;
* you need fast sensitivity checks;
* lateral structure is not the main interpretation target.

Use a **2-D** path when:

* stations follow a profile;
* geology is approximately constant along strike;
* TE/TM mode interpretation is meaningful;
* off-profile effects are limited or can be discussed;
* a section is the main deliverable.

Use a **3-D** path when:

* stations cover an area;
* geology varies in all directions;
* tipper, phase tensor, induction arrows, or residuals suggest 3-D behavior;
* 2-D inversions leave systematic station/component misfits;
* the project can support the larger computation and modelling effort.

Use a **2.5-D** path such as :term:`MARE2DEM` when:

* the model is a 2-D section but sources/receivers or fields need richer
  finite-element treatment;
* MT/CSEM project files are already MARE2DEM-native;
* topography and mesh geometry must be represented explicitly.

Decision Axis 2: Data Type
--------------------------

Once dimensionality is defensible, match the data to the physics that produced
it. Passive :term:`MT`/:term:`AMT` workflows usually interpret the impedance
relation

.. math::

   \mathbf{E}(\omega) = \mathbf{Z}(\omega)\,\mathbf{H}(\omega),

where :math:`\mathbf{Z}` is the complex :term:`impedance tensor` at angular
frequency :math:`\omega=2\pi f`. :term:`Apparent resistivity` and
:term:`phase` then reduce each usable impedance component to
:math:`\rho_a=(\mu_0\omega)^{-1}|Z|^2` and
:math:`\phi=\tan^{-1}(\operatorname{Im}Z/\operatorname{Re}Z)`. In
:term:`CSAMT` and :term:`CSEM`, the transmitter is explicit, so the same
far-field impedance interpretation is only appropriate after source effects
have been checked. In :term:`TDEM`, the measured signal is a transient decay
rather than a steady frequency-domain impedance, so a time-gate inversion and
an EDI-style conversion answer different questions.

.. list-table::
   :header-rows: 1
   :widths: 22 32 46

   * - Data type
     - Suitable paths
     - Notes
   * - MT
     - Built-in, Occam2D, ModEM, SimPEG, pyGIMLi, MARE2DEM
     - Choose by dimensionality and deliverable requirements.
   * - AMT
     - Built-in, Occam2D-style profile workflows, SimPEG, pyGIMLi
     - Often profile-oriented; static shift and cultural noise are common.
   * - CSAMT
     - Built-in checks, Occam2D-style profile workflows, MARE2DEM for
       controlled-source native projects
     - Confirm far-field assumptions and source effects before choosing an
       MT-style inversion.
   * - TDEM
     - Built-in, pyGIMLi, SimPEG, TDEM conversion tools
     - Direct decay inversion and EDI conversion are different choices.
   * - CSEM
     - MARE2DEM
     - Use MARE2DEM-native data and geometry when controlled-source physics
       is central.

Decision Axis 3: Operational Control
------------------------------------

Choose the backend-neutral :mod:`pycsamt.inversion` route when:

* you want a common configuration object;
* you need one result interface across several backend experiments;
* you are preparing or loading external projects rather than manually tuning
  every native file;
* pipeline or agent workflows will orchestrate the run.

Choose the direct :mod:`pycsamt.models` route when:

* native files are the deliverable;
* you need full control over mesh, model, covariance, startup, or settings;
* you need to run, reload, compare, or plot native outputs independently;
* the external executable is managed outside the Python workflow;
* the project requires detailed provenance at the file level.

Both routes eventually compare observed data :math:`\mathbf{d}_{obs}` with a
prediction from a :term:`forward operator`, :math:`F(\mathbf{m})`, for model
parameters :math:`\mathbf{m}`. A common smooth inversion minimizes

.. math::

   \Phi(\mathbf{m}) =
   \left\|\mathbf{W}_d\left(\mathbf{d}_{obs} - F(\mathbf{m})\right)\right\|_2^2
   + \lambda^2\left\|\mathbf{W}_m\left(\mathbf{m}-\mathbf{m}_{ref}\right)\right\|_2^2,

where :math:`\mathbf{W}_d` weights data by their errors,
:math:`\mathbf{W}_m` controls model roughness or departure from a reference
model :math:`\mathbf{m}_{ref}`, and :math:`\lambda` balances fit against
smoothness. Backend choice decides the form of :math:`F`, the available mesh
and model controls, and how transparently those choices are written to native
files.

Typical File Responsibility
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 18 34 48

   * - Engine
     - User must understand
     - pyCSAMT helps with
   * - Occam2D
     - Data file, mesh, model, startup, response, log, target RMS.
     - Builders, validation, runner, result loader, response/model plots.
   * - ModEM
     - Data components, model grid, covariance, control file, response/logs.
     - Data/model/covariance/control classes, runner, iotools, plots.
   * - MARE2DEM
     - Source/binary, ``.emdata``, ``.poly``, ``.resistivity``,
       ``.settings``, MPI, topography and geometry.
     - Source manager, file parsers/writers, geometry utilities, runner,
       result loader, plots.

External Executables
--------------------

External engines should be treated as scientific dependencies. Before a
production run, document:

* exact executable name and path;
* source version or build provenance;
* compiler and MPI details where relevant;
* environment variables and library paths;
* run command and working directory;
* whether pyCSAMT executed the binary or only prepared/loaded files.

In :class:`pycsamt.inversion.config.InversionConfig`, ``run_external=False``
is a conservative default. It lets pyCSAMT prepare or validate a working
directory without launching an external program unexpectedly.

Workflow Map
------------

.. code-block:: text
   :linenos:

   Field data
      |
      v
   QC, coordinates, units, static shift, source checks
      |
      v
   Choose modelling assumption
      |
      +--> 1-D: builtin / pyGIMLi / SimPEG
      |
      +--> smooth 2-D profile: Occam2D
      |
      +--> 3-D MT: ModEM
      |
      +--> 2.5-D MT/CSEM FEM: MARE2DEM
      |
      v
   Build native files or backend-neutral config
      |
      v
   Validate before execution
      |
      v
   Run or load external result
      |
      v
   Plot responses, residuals, model, convergence
      |
      v
   Interpret with uncertainty and archive provenance

Connection To InversionConfig
-----------------------------

The backend-neutral route uses :class:`pycsamt.inversion.config.InversionConfig`.
The example below keeps ``run_external=False`` so it can be used to validate
configuration and file responsibility before any external executable is
launched.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.inversion import InversionConfig, InversionWorkflow
   >>> cfg = InversionConfig(
   ...     method="mt",
   ...     dimension="2d",
   ...     backend="occam2d",
   ...     data="data/3edis",
   ...     workdir="runs/occam2d_profile",
   ...     run_external=False,
   ...     backend_options={
   ...         "files": {
   ...             "data": "OccamData.dat",
   ...             "mesh": "OccamMesh",
   ...             "model": "OccamModel",
   ...         },
   ...     },
   ... )
   >>> cfg.backend
   'occam2d'
   >>> cfg.run_external
   False
   >>> cfg.workdir
   'runs/occam2d_profile'
   >>> # result = InversionWorkflow(cfg).run()

This style is useful for pipeline-friendly experiments. For native-project
work, instantiate the model package directly.

Direct Model Integration Example
--------------------------------

Direct model use gives you the native builder/runner/result objects:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.models.occam2d import (
   ...     InputBuilder,
   ...     InversionResult,
   ...     OccamConfig,
   ...     OccamRunner,
   ... )
   >>> from pycsamt.emtools import ensure_sites
   >>> sites = ensure_sites("data/3edis", recursive=True, verbose=0)
   >>> type(sites).__name__, len(sites)
   ('Sites', 3)
   >>> cfg = OccamConfig()
   >>> builder = InputBuilder(sites, workdir="runs/occam2d_profile", config=cfg)
   >>> builder.__class__.__name__
   'InputBuilder'
   >>> runner = OccamRunner("runs/occam2d_profile", config=cfg)
   >>> runner.__class__.__name__
   'OccamRunner'
   >>> # builder.build()
   >>> # runner.run(target_misfit=1.0)
   >>> # result = InversionResult("runs/occam2d_profile", config=cfg)

Use the same pattern mentally for ModEM and MARE2DEM: configuration,
builder, optional runner, result loader, plots, archive.

Risk Checklist
--------------

Before committing to a backend, check:

* **Dimensionality** - Do data diagnostics support 1-D, 2-D, 2.5-D, or 3-D?
* **Coordinates** - Are station coordinates, elevations, and profile axes
  correct?
* **Units** - Are impedance, apparent resistivity, phase, time gates, or CSEM
  values in expected units?
* **Errors** - Are :term:`error floor`\ s realistic, and are bad samples
  masked?
* **Source physics** - Are CSAMT/CSEM/TDEM source assumptions represented?
* **Topography** - Does the backend need explicit topography or air cells?
* **Executable** - Can the binary run reproducibly on the target machine?
* **Native files** - Can another user inspect and rerun the generated project?
* **Deliverable** - Does the final report require smooth sections, 3-D
  volumes, response fits, native logs, or hydrogeological handoff?

Common Mistakes
---------------

Avoid these backend-selection mistakes:

* using 3-D inversion because it looks more advanced, without enough spatial
  coverage;
* using a 1-D backend for a profile with obvious lateral contacts;
* applying an MT-style backend to CSAMT data without checking source effects;
* treating MARE2DEM as an ``InversionConfig`` backend when using the direct
  model integration is required;
* hiding native file choices inside an undocumented script;
* running an external binary before validating data, mesh, and coordinates;
* judging the backend only by final :term:`RMS misfit` instead of residual
  structure and geological plausibility.

Recommended Path By Deliverable
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 30 36

   * - Deliverable
     - Recommended start
     - Why
   * - Fast reconnaissance model
     - Built-in or pyGIMLi 1-D
     - Quick, simple, and useful for checking data scale.
   * - Publication-quality 2-D profile section
     - Occam2D plus residual plots
     - Smooth model, target RMS, native convergence evidence.
   * - 3-D MT conductivity volume
     - ModEM
     - Mature 3-D file ecosystem and covariance controls.
   * - MT/CSEM finite-element project
     - MARE2DEM
     - Native support for 2.5-D FEM geometry and controlled-source files.
   * - Pipeline report with several backend trials
     - ``InversionConfig`` plus pipeline outputs
     - Common result vocabulary and reproducible run directories.

Next Steps
----------

* :doc:`occam2d` for smooth 2-D profile inversion.
* :doc:`modem` for ModEM 2-D/3-D workflows.
* :doc:`mare2dem` for MARE2DEM MT/CSEM projects.
* :doc:`configuration_and_io` for native file and run-directory practice.
* :ref:`inversion_concepts` for objective functions and model assumptions.
