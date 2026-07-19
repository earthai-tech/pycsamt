.. _user_guide_models_overview:

Overview
========

The :mod:`pycsamt.models` layer is for projects where the modelling engine
itself is part of the scientific record. It sits close to native files,
external executables, solver logs, and review figures, while still letting the
rest of pyCSAMT load, compare, and archive results in a consistent way.

Treat this page as the entrance to a reproducible modelling chain. The
question is not only "which program can invert these data?", but "which set of
assumptions can another geophysicist inspect, rerun, and defend?" That means
the :term:`mesh`, :term:`starting model`, :term:`regularization`,
:term:`error floor`, coordinate system, source treatment, executable identity,
and final :term:`residual` structure all belong in the project record.

At a glance
-----------

.. list-table::
   :header-rows: 1
   :widths: 18 22 25 35

   * - Integration
     - Typical geometry
     - Typical methods
     - Good starting point
   * - :term:`Occam2D`
     - Smooth 2-D profile
     - :term:`MT`, :term:`AMT`, :term:`CSAMT`
     - A line-oriented survey where the :term:`2-D` assumption is defensible
       and a smooth resistivity section is the intended product.
   * - :term:`ModEM`
     - 2-D or 3-D grid
     - Primarily :term:`MT` and :term:`AMT`
     - A station grid, important off-profile structure,
       :term:`covariance`-controlled smoothing, or an existing native ModEM
       project.
   * - :term:`MARE2DEM`
     - 2.5-D finite-element section
     - :term:`MT` and :term:`CSEM`
     - :term:`topography`, CSEM source/receiver geometry, adaptive triangular
       finite-element meshes, or an existing MARE2DEM workflow.

Selection is a scientific decision, not only a file-format choice. Review
survey geometry, :term:`dimensionality`, :term:`geoelectric strike`,
topography, source configuration, target scale, available compute resources,
and the assumptions of each solver before committing to a
:term:`inversion backend` or direct :term:`model integration`.

The practical boundary is easiest to see through the model function. For a
1-D sounding, resistivity is written as :math:`\rho(z)` and each station can be
reviewed almost independently. For a 2-D profile, the model becomes
:math:`\rho(x,z)` and the interpretation depends on a defensible strike
direction. For a 3-D volume, :math:`\rho(x,y,z)` can explain off-profile
structure, but only if station coverage, errors, and computation are strong
enough to support the extra freedom. A 2.5-D finite-element workflow keeps a
2-D earth section while allowing richer source, receiver, and field behaviour
than a simple plane-wave profile assumption.

Shared workflow
---------------

Despite different native formats, a responsible engine workflow follows the
same broad sequence. The order matters because every later modelling choice
inherits the units, coordinates, masks, and assumptions accepted at the start:

#. review the processed observations and dimensionality evidence;
#. choose the backend and model dimension;
#. create an isolated, versioned run directory;
#. configure data components, errors, mesh, starting model, and regularization;
#. build native files without overwriting source observations;
#. validate shapes, units, coordinates, file references, and executable setup;
#. inspect or dry-run the assembled command;
#. launch the external executable only when explicitly intended;
#. monitor logs, iteration history, residuals, and solver termination;
#. load final and intermediate products back into pyCSAMT;
#. compare observed and predicted responses before interpreting the model;
#. archive native files, configuration, command, executable identity, and
   review products together.

The integration page for each engine expands this sequence gradually with
step-by-step examples grounded in its real :mod:`pycsamt.models` modules. The
overview here is intentionally engine-neutral: it defines what must remain
visible no matter whether the final files are Occam data/startup files, ModEM
data/model/covariance/control files, or MARE2DEM ``.emdata``/``.poly``/
``.settings`` projects.

What the inversion is minimizing
--------------------------------

All three native integrations ultimately compare observed data
:math:`\mathbf{d}_{obs}` with predicted data
:math:`\mathbf{d}_{pred}=F(\mathbf{m})`, where :math:`F` is the
:term:`forward operator` and :math:`\mathbf{m}` is the vector of model
parameters. In an electromagnetic inversion, :math:`\mathbf{d}` may contain
complex impedance components, :term:`apparent resistivity`, :term:`phase`,
:term:`tipper` values, or controlled-source fields. A standard weighted data
misfit is

.. math::

   \Phi_d =
   \left\|\mathbf{W}_d\left(\mathbf{d}_{obs}-F(\mathbf{m})\right)\right\|_2^2,

with :math:`\mathbf{W}_d` normally built from inverse data standard
deviations. This is why error floors are scientific parameters, not clerical
settings: if a datum is assigned an uncertainty :math:`\sigma_i`, its
contribution scales like :math:`(d_i-F_i(\mathbf{m}))^2/\sigma_i^2`.
Underestimating :math:`\sigma_i` can force the inversion to chase noise or
source effects; overestimating it can hide real structure.

Because EM inverse problems are usually :term:`non-uniqueness`\ -limited, the
misfit is paired with a model penalty:

.. math::

   \Phi(\mathbf{m}) = \Phi_d + \lambda^2 \Phi_m,

where :math:`\Phi_m` expresses smoothness, damping, or departure from a
reference model and :math:`\lambda` controls the trade-off. Occam-style
workflows emphasize the smoothest model that reaches an acceptable target
misfit. ModEM exposes much of this behaviour through covariance and control
files. MARE2DEM wraps the same appraisal problem in a finite-element geometry
where source and topography choices become especially visible.

Model integrations versus the common inversion API
--------------------------------------------------

pyCSAMT separates backend-neutral orchestration from native engine tooling:

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Layer
     - Purpose
     - Use it when
   * - :mod:`pycsamt.inversion`
     - Common configuration, workflow, result handling, plotting, and export.
     - You want one interface across built-in and external backends.
   * - :mod:`pycsamt.models`
     - Native builders, file objects, validators, runners, result loaders, and
       engine-specific plots.
     - You need direct control of Occam2D, ModEM, or MARE2DEM project files.
   * - :mod:`pycsamt.pipeline`
     - Repeatable multi-step processing and project orchestration.
     - The same preparation or review chain must be rerun consistently.
   * - :mod:`pycsamt.agents`
     - Assisted workflow routing, configuration guidance, and reporting.
     - You want agent-supported orchestration over the underlying science APIs.

The layers complement one another. The app and agent interfaces do not replace
the native science objects when detailed engine control is required. A good
rule is to use :mod:`pycsamt.inversion` when comparability across backends is
the main goal, and :mod:`pycsamt.models` when native files are the evidence
that must survive review.

Engine assumptions in plain terms
---------------------------------

:term:`Occam2D` is a profile tool. It is strongest when stations lie along a
line and geological variation along strike is small enough that the 2-D
assumption is meaningful. Its smoothness is a feature, not a cosmetic choice:
the model is deliberately restrained so features that remain after reasonable
regularization have more interpretive weight.

:term:`ModEM` is the better starting point when a 3-D conductivity volume is
scientifically required, or when an existing ModEM project already controls
data components, grid geometry, covariance, and inversion settings. The price
of that freedom is bookkeeping. Coordinate consistency, inactive cells,
topography handling, and component-specific error floors must be reviewed
before the final volume is interpreted.

:term:`MARE2DEM` is useful when a 2-D section is still the modelling target but
the source and receiver physics need more explicit treatment. In CSEM and
controlled-source MT work, the transmitter is not just metadata: source
position, orientation, frequency, current, and receiver geometry are part of
the forward problem. MARE2DEM therefore belongs in the workflow when that
native geometry needs to be preserved rather than approximated away.

Run-directory principles
------------------------

Use one directory per backend, dataset, configuration, and revision. A useful
pattern is:

.. code-block:: text

   inversion_runs/
   `-- <line>_<backend>_<date>_<revision>/
       |-- source/          # immutable prepared observations
       |-- configuration/   # templates and resolved settings
       |-- native/          # engine input and output files
       |-- logs/            # command, stdout, stderr, solver logs
       |-- review/          # residuals, convergence, response plots
       `-- manifest.yml     # provenance and file roles

Do not mix files from competing meshes, error floors, or starting models in a
single native directory. Engine outputs often have conventional names and can
silently become inconsistent when copied between runs.

The run directory should answer four questions without relying on memory:

* What exact observations entered the inversion?
* Which model space, mesh, priors, errors, and source assumptions were used?
* Which command and executable produced the result?
* Which residuals, convergence traces, and comparison plots justify accepting
  or rejecting the result?

If any one of those answers lives only in a notebook cell, desktop GUI state,
or field notebook margin, the run is not yet reproducible.

Minimum review evidence
-----------------------

Before passing a model into :doc:`../interpretation/index`, retain:

* source data and preprocessing identifiers;
* station coordinates, profile direction, components, and units;
* backend and executable identity;
* resolved configuration and command;
* native data, mesh/model, control/startup, covariance, and response files as
  applicable;
* starting and final models;
* convergence and RMS history;
* residuals by station, period/frequency, and component;
* observed-versus-predicted responses;
* sensitivity or resolution evidence appropriate to the backend;
* failed or alternative runs relevant to model appraisal;
* reviewer notes, limitations, and accepted interpretation depth.

Never select a final model from RMS alone. A low global RMS can hide structured
misfit, unrealistic parameters, error-floor problems, or geometry that violates
the assumed dimension. The stronger check is residual anatomy: ask whether
misfit clusters by station, component, period/frequency band, line segment, or
data type. A model with slightly higher RMS but random, explainable residuals
can be more defensible than a lower-RMS model that repeatedly misses the same
component or depth range.

The accepted interpretation depth should also be stated separately from the
maximum plotted depth. A model cell may be drawn because it exists in the mesh,
not because the data have useful sensitivity there. As a practical review
habit, keep the plotted model, response fits, residual section, convergence
curve, and sensitivity or resolution evidence together in the same review
folder.

Related reading
---------------

* :doc:`choosing_backend` for deciding between the common inversion API and
  the direct model integrations.
* :doc:`configuration_and_io` for native file and run-directory practice.
* :doc:`../../theory/inversion_concepts` for objective functions, errors,
  regularization, RMS, resolution, non-uniqueness, and model appraisal.
* :doc:`../../api/models` for generated class and function signatures.
