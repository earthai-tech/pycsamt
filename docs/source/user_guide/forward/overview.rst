.. _forward_overview:

Overview
========

Forward And Model Integrations
---------------------------------

The :mod:`pycsamt.forward` and :mod:`pycsamt.models` packages are related, but
they should not be treated as the same layer.

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Package
     - Main responsibility
     - Typical use
   * - :mod:`pycsamt.forward`
     - Compute controlled synthetic responses inside Python using pyCSAMT
       model containers, solvers, dataset generators, noise models, and
       plotting helpers.
     - Test survey design, generate training data, validate inversion
       assumptions, compare response behaviour, and build reproducible
       synthetic examples.
   * - :mod:`pycsamt.models`
     - Integrate external modelling and inversion engines such as Occam2D,
       ModEM, and MARE2DEM.
     - Prepare native engine files, run external binaries, load engine
       outputs, and manage backend-specific inversion projects.

Use :mod:`pycsamt.forward` when you need a controlled synthetic response from
a known model. Use :mod:`pycsamt.models` when you need to prepare, run, or load
a native external-engine project. The two layers meet when a forward synthetic
experiment is converted into an inversion benchmark.

Package Map
-------------

The public forward API is built around a small set of object families.

.. list-table::
   :header-rows: 1
   :widths: 28 34 38

   * - Object family
     - Main objects
     - Role
   * - Configuration
     - ``ForwardConfig``, ``ForwardConfig2D``, ``ForwardConfig3D``
     - Store reproducible frequency/time axes, model ranges, station layouts,
       solver options, noise settings, and output choices.
   * - Model containers
     - ``LayeredModel``, ``Grid2D``, ``Grid3D``
     - Represent 1-D layered earths, 2-D profile grids, and quasi-3-D
       resistivity volumes.
   * - Solvers
     - ``MT1DForward``, ``CSAMT1DForward``, ``TEM1DForward``, ``MT2DForward``,
       ``MT3DForward``
     - Compute predicted electromagnetic responses from model containers.
   * - Response containers
     - ``ForwardResponse``, ``ForwardResponse2D``, ``ForwardResponse3D``
     - Hold predicted apparent resistivity, phase, transient values,
       impedance components, station geometry, and array conversion helpers.
   * - Noise models
     - ``GaussianNoise``, ``MultiplicativeNoise``, ``FieldRealisticNoise``
     - Convert ideal synthetic responses into more realistic observations.
   * - Dataset containers
     - ``ForwardDataset``, ``SurveyDataset3D``
     - Store synthetic feature and target arrays, metadata, splits, and
       compressed ``.npz`` archives.
   * - Plotting helpers
     - ``plot_response_1d``, ``plot_model_2d``, ``plot_pseudosection_2d``,
       ``plot_model_3d``, ``plot_response_map_3d``, and related functions
     - Inspect models, responses, tensor components, noisy examples, and
       generated datasets before they are trusted downstream.

Core Workflow
---------------

A typical forward modelling workflow has six stages.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.forward import LayeredModel, MT1DForward, plot_response_and_model_1d

   freqs = np.logspace(-2, 3, 40)
   model = LayeredModel(
       resistivity=[100.0, 15.0, 800.0],
       thickness=[250.0, 700.0],
   )

   solver = MT1DForward(freqs)
   response = solver.run(model)

   fig = plot_response_and_model_1d(response, model)
   fig.savefig("runs/forward/mt1d_response.png", dpi=200)

For production work, the same idea should be driven by a configuration file
and archived with outputs:

#. define the physical question and target dimensionality;
#. create a ``ForwardConfig*`` template;
#. build a model or grid from the configuration;
#. run the selected solver;
#. apply a documented noise model when simulating observations;
#. plot and archive the model, response, configuration, and metadata.

Choosing A Path
------------------

Different users enter the forward section with different goals.

.. list-table::
   :header-rows: 1
   :widths: 30 34 36

   * - Goal
     - Start here
     - Then read
   * - Understand the physics and vocabulary
     - :doc:`concepts`
     - :doc:`solvers_and_grids`, :doc:`plotting`
   * - Build a reproducible synthetic run
     - :doc:`configuration`
     - :doc:`solvers_and_grids`, :doc:`forward_to_inversion`
   * - Generate training data for AI inversion
     - :doc:`synthetic_datasets`
     - :doc:`configuration`, :doc:`plotting`,
       :doc:`../agents/ai_model_zoo_agents`
   * - Compare 1-D, 2-D, and quasi-3-D behaviour
     - :doc:`solvers_and_grids`
     - :doc:`concepts`, :doc:`plotting`
   * - Prepare an inversion benchmark
     - :doc:`forward_to_inversion`
     - :doc:`../models/index`, :doc:`../pipeline/index`
   * - Diagnose a generated dataset or response
     - :doc:`plotting`
     - :doc:`synthetic_datasets`, :doc:`forward_to_inversion`

Relationship To Theory
-------------------------

Forward modelling is practical, but it is not detached from theory. When a
plot or synthetic response looks surprising, the relevant background pages are:

* :doc:`../../theory/csamt_amt_mt_overview` for the distinction between CSAMT,
  AMT, MT, and TDEM survey assumptions;
* :doc:`../../theory/impedance_tensor` for impedance, apparent resistivity, phase,
  and tensor notation;
* :doc:`../../theory/static_shift` for shallow distortion effects that can make
  synthetic and field curves disagree;
* :doc:`../../theory/inversion_concepts` for how forward responses are used inside
  objective functions and regularized inversion.
