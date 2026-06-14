.. _forward:

Forward Modelling
=================

:mod:`pycsamt.forward` contains pyCSAMT's in-process forward modelling and
synthetic data generation tools. These tools are used to test physical
assumptions, generate training data for AI inversion, compare model responses,
and prepare intuition before running expensive inversion workflows.

Forward modelling answers the question:

   If this resistivity model were true, what response would the survey record?

That question is essential for inversion. Inversion searches for a model that
matches observed data, while forward modelling computes the data predicted by
one known model. A reliable inversion workflow therefore depends on forward
models for synthetic tests, sensitivity checks, data-feature design, and
quality control.

Forward And Models
------------------

The :mod:`pycsamt.forward` and :mod:`pycsamt.models` packages are related but
serve different purposes.

.. list-table::
   :header-rows: 1
   :widths: 26 37 37

   * - Package
     - Main purpose
     - Typical use
   * - :mod:`pycsamt.forward`
     - In-process 1-D, 2-D, quasi-3-D forward responses, synthetic models,
       noise, datasets, and plotting.
     - Generate synthetic data, test model sensitivity, train AI inversion
       networks, and inspect response behaviour.
   * - :mod:`pycsamt.models`
     - Integration with external modelling and inversion engines such as
       Occam2D, ModEM, and MARE2DEM.
     - Prepare native files, run external binaries, load outputs, and plot
       engine-specific results.

Use ``forward`` when you need a controlled synthetic response. Use ``models``
when you need to prepare or load a native external-engine project.

Package Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Object or module
     - Role
   * - ``ForwardConfig``
     - Configuration for 1-D MT, CSAMT, and TEM dataset generation.
   * - ``ForwardConfig2D``
     - Configuration for 2-D MT finite-difference runs.
   * - ``ForwardConfig3D``
     - Configuration for quasi-3-D MT forward runs.
   * - ``LayeredModel``
     - 1-D layered earth model used by MT1D, CSAMT1D, TEM1D, and synthetic
       dataset generation.
   * - ``Grid2D`` and ``Grid3D``
     - Non-uniform resistivity grids with station layouts.
   * - ``MT1DForward``, ``CSAMT1DForward``, ``TEM1DForward``
     - 1-D electromagnetic forward solvers.
   * - ``MT2DForward``
     - 2-D MT finite-difference forward solver returning TE and TM responses.
   * - ``MT3DForward``
     - Quasi-3-D MT forward solver returning impedance tensor components.
   * - ``ForwardDataset`` and ``SurveyDataset3D``
     - Containers for synthetic datasets saved as compressed ``.npz`` files.
   * - ``GaussianNoise``, ``MultiplicativeNoise``, ``FieldRealisticNoise``
     - Noise models used to make synthetic responses more field-like.
   * - ``plot_*`` functions
     - Diagnostic plots for 1-D responses, 2-D models, pseudo-sections, 3-D
       maps, and tensor components.

Recommended Reading Order
-------------------------

For a user who wants to understand the forward layer:

#. Read :doc:`concepts`.
#. Read :doc:`configuration`.
#. Read :doc:`solvers_and_grids`.
#. Read :doc:`synthetic_datasets`.
#. Read :doc:`plotting`.
#. Read :doc:`forward_to_inversion` to connect synthetic tests with inversion.

For AI-assisted inversion or training datasets:

#. Start with :doc:`synthetic_datasets`.
#. Review :doc:`configuration` for reproducible source-of-truth files.
#. Use :doc:`plotting` to inspect generated examples.
#. Continue to :doc:`../agents/ai_model_zoo_agents` and
   :doc:`../user_guide/ai_inversion`.

Contents
--------

.. toctree::
   :maxdepth: 1

   concepts
   configuration
   solvers_and_grids
   synthetic_datasets
   plotting
   forward_to_inversion

