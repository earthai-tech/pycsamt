.. _models:

Classical model integrations
============================

The :mod:`pycsamt.models` package connects pyCSAMT to three established
electromagnetic modelling and inversion ecosystems: Occam2D, ModEM, and
MARE2DEM. These integrations keep native files and external executables
visible while providing Python tools to build inputs, validate projects, run
solvers when explicitly requested, load results, and create review figures.

This section is the classical engine layer of the
:doc:`../inversion/index` guide. It is separate from
:doc:`../forward/index`, which focuses on predicting synthetic responses, and
from :doc:`../ai_inversion/index`, which covers learned inversion methods.

.. admonition:: External solvers
   :class: important

   pyCSAMT provides integration code and project tooling. Availability of an
   external solver, its executable license, compilation requirements, and
   computational resources remain the user's responsibility. External
   execution is explicit rather than automatic.

Choose an integration
---------------------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Occam2D
      :link: occam2d
      :link-type: doc
      :img-top: ../../_static/icons/model.svg
      :class-card: pycsamt-card sd-text-center

      Prepare and review smooth 2-D MT, AMT, and CSAMT profile inversions with
      explicit data, mesh, model, startup, response, and iteration files.

   .. grid-item-card:: ModEM
      :link: modem
      :link-type: doc
      :img-top: ../../_static/icons/model-builder.svg
      :class-card: pycsamt-card sd-text-center

      Build and inspect ModEM 2-D or 3-D projects, including data, model,
      covariance, control, response, iteration, and native run products.

   .. grid-item-card:: MARE2DEM
      :link: mare2dem
      :link-type: doc
      :img-top: ../../_static/icons/external-validation.svg
      :class-card: pycsamt-card sd-text-center

      Work with 2.5-D finite-element MT and CSEM projects, triangular meshes,
      topography, native data/settings files, and MPI-oriented execution.

At a glance
-----------

.. list-table::
   :header-rows: 1
   :widths: 18 22 25 35

   * - Integration
     - Typical geometry
     - Typical methods
     - Good starting point
   * - Occam2D
     - Smooth 2-D profile
     - MT, AMT, CSAMT
     - A line-oriented survey where the 2-D assumption is defensible and a
       smooth resistivity section is the intended product.
   * - ModEM
     - 2-D or 3-D grid
     - Primarily MT and AMT
     - A station grid, important off-profile structure, covariance-controlled
       smoothing, or an existing native ModEM project.
   * - MARE2DEM
     - 2.5-D finite-element section
     - MT and CSEM
     - Topography, CSEM source/receiver geometry, adaptive triangular meshes,
       or an existing MARE2DEM workflow.

Selection is a scientific decision, not only a file-format choice. Review
survey geometry, dimensionality, strike, topography, source configuration,
target scale, available compute resources, and the assumptions of each solver
before committing to a backend.

Shared workflow
---------------

Despite different native formats, a responsible engine workflow follows the
same broad sequence:

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

The integration page for each engine will expand this sequence gradually with
step-by-step examples grounded in its real :mod:`pycsamt.models` modules.

Supporting guides
-----------------

.. grid:: 1 1 2 2
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Choose a backend
      :link: choosing_backend
      :link-type: doc
      :img-top: ../../_static/icons/user-guide-inversion.svg
      :class-card: pycsamt-card sd-text-center

      Compare geometry, physics, native files, runtime requirements, and
      expected outputs before creating a project.

   .. grid-item-card:: Configuration and I/O
      :link: configuration_and_io
      :link-type: doc
      :img-top: ../../_static/icons/controls.svg
      :class-card: pycsamt-card sd-text-center

      Organize reproducible run folders, shared configuration, validation,
      executable paths, file naming, and provenance.

   .. grid-item-card:: Inversion concepts
      :link: ../../theory/inversion_concepts
      :link-type: doc
      :img-top: ../../_static/icons/inversion-concepts.svg
      :class-card: pycsamt-card sd-text-center

      Review objective functions, errors, regularization, RMS, resolution,
      non-uniqueness, and model appraisal.

   .. grid-item-card:: API reference
      :link: ../../api/models
      :link-type: doc
      :img-top: ../../_static/icons/api-reference-icon-braces.svg
      :class-card: pycsamt-card sd-text-center

      Consult generated class and function signatures after choosing the
      appropriate narrative workflow.

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
the native science objects when detailed engine control is required.

Run-directory principles
------------------------

Use one directory per backend, dataset, configuration, and revision. A useful
pattern is:

.. code-block:: text

   inversion_runs/
   └── <line>_<backend>_<date>_<revision>/
       ├── source/          # immutable prepared observations
       ├── configuration/   # templates and resolved settings
       ├── native/          # engine input and output files
       ├── logs/            # command, stdout, stderr, solver logs
       ├── review/          # residuals, convergence, response plots
       └── manifest.yml     # provenance and file roles

Do not mix files from competing meshes, error floors, or starting models in a
single native directory. Engine outputs often have conventional names and can
silently become inconsistent when copied between runs.

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
the assumed dimension.

Contents
--------

.. toctree::
   :maxdepth: 1

   choosing_backend
   configuration_and_io
   occam2d
   modem
   mare2dem
