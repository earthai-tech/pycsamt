.. _forward_concepts:

Forward Modelling Concepts
==========================

Forward modelling is the controlled side of electromagnetic interpretation.
The model is known, the survey geometry is known, and the solver computes the
response. Inversion reverses this relationship: the response is known and the
model is estimated.

In pyCSAMT, forward modelling is used for five related tasks:

* computing synthetic MT, CSAMT, and TEM responses from known models;
* testing whether a target should be visible at a planned frequency range,
  time range, and station spacing;
* generating training data for AI inversion;
* validating inversion settings with synthetic recovery tests;
* producing diagnostic plots that explain why an inversion succeeds or fails.

The forward package is intentionally separate from
:mod:`pycsamt.models`. :mod:`pycsamt.forward` contains in-process solvers,
grids, synthetic models, noise models, and datasets. :mod:`pycsamt.models`
wraps external engines such as Occam2D, ModEM, and MARE2DEM.

The Forward Problem
-------------------

For an electromagnetic forward problem, pyCSAMT usually works with:

* a resistivity model :math:`\rho(\mathbf{x})`;
* a source assumption, such as a plane-wave MT source, a controlled-source
  approximation, or a TEM transmitter loop;
* receiver locations and sampled frequencies or time gates;
* a solver that implements a physical approximation;
* response quantities such as impedance, apparent resistivity, phase, or
  transient decay.

The basic relationship is:

.. math::

   d_\mathrm{pred} = F(m),

where :math:`m` is the model, :math:`F` is the forward operator, and
:math:`d_\mathrm{pred}` is the predicted data vector.

The same expression appears inside inversion:

.. math::

   m_\mathrm{inv} =
   \arg\min_m
   \left[
      \| W_d(F(m) - d_\mathrm{obs}) \|_2^2
      + \lambda R(m)
   \right].

The forward operator :math:`F` is therefore not a side detail. It defines how
models are compared to data. If the forward assumption is wrong, the inversion
can fit the data while still producing a misleading model.

Physical Inputs
---------------

A useful forward experiment should make every physical input explicit.

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Input
     - pyCSAMT object or setting
     - Why it matters
   * - Resistivity model
     - ``LayeredModel``, ``Grid2D``, ``Grid3D``
     - Defines where current can flow and where fields diffuse.
   * - Dimensionality
     - 1-D, 2-D, quasi-3-D
     - Controls which geological structures the solver can represent.
   * - Frequency grid
     - ``freqs`` or ``ForwardConfig.freq_grid()``
     - Sets the period range and depth sensitivity for MT/CSAMT.
   * - Time gates
     - ``times`` or ``ForwardConfig.time_grid()``
     - Sets the transient sampling for TEM.
   * - Source assumption
     - MT plane wave, CSAMT source offset, TEM loop
     - Determines whether the response is natural-source, controlled-source,
       or transient.
   * - Station layout
     - ``x_stations`` or ``stations_xy``
     - Controls lateral sampling and pseudo-section geometry.
   * - Noise model
     - ``GaussianNoise``, ``MultiplicativeNoise``,
       ``FieldRealisticNoise``
     - Makes synthetic data closer to field conditions.

The minimum reproducible forward run should record the model, solver, sampled
axis, station geometry, noise model, random seed, and output file.

Dimensionality
--------------

The forward package exposes three practical dimensionality levels. Each is a
deliberate approximation, not merely a larger array shape.

.. list-table::
   :header-rows: 1
   :widths: 16 30 30 24

   * - Level
     - Objects
     - Best use
     - Main limitation
   * - 1-D
     - ``LayeredModel``, ``MT1DForward``, ``CSAMT1DForward``,
       ``TEM1DForward``
     - Fast soundings, method comparison, synthetic AI datasets, and
       sanity checks.
     - No lateral structure.
   * - 2-D
     - ``Grid2D``, ``MT2DForward``
     - Survey-line responses, TE/TM comparison, anomaly tests, and
       pseudo-sections.
     - Assumes structure is constant along strike.
   * - Quasi-3-D
     - ``Grid3D``, ``MT3DForward``
     - Area-style station layouts, approximate tensor responses, and
       survey-scale synthetic experiments.
     - Approximate 3-D physics; not a replacement for a full 3-D production
       solver.

Use the simplest dimensionality that can represent the question. If the
question is "Can a conductive layer be detected?", 1-D may be enough. If the
question is "Can stations across a profile detect a buried block?", use 2-D.
If the question involves area coverage, tensor components, or spatially
correlated training data, quasi-3-D is more appropriate.

1-D Forward Models
------------------

The 1-D model is a stack of horizontal layers. Resistivity varies with depth
only:

.. math::

   \rho(\mathbf{x}) = \rho(z).

:class:`pycsamt.forward.LayeredModel` stores layer resistivities and
thicknesses. The final layer is the halfspace and has no thickness value.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.forward import LayeredModel, MT1DForward

   model = LayeredModel(
       resistivity=[100.0, 10.0, 500.0],
       thickness=[300.0, 800.0],
       name="conductive_middle_layer",
   )

   freqs = np.logspace(-3, 4, 40)
   response = MT1DForward(freqs).run(model)

   print(response.rho_a.shape)
   print(response.phase.shape)

1-D forward modelling is fast because each frequency or time gate only needs
to propagate fields through a layer stack. That makes it ideal for:

* synthetic training datasets with thousands or millions of samples;
* checking apparent resistivity and phase behaviour;
* designing frequency grids;
* comparing MT, CSAMT, and TEM sensitivity;
* teaching and debugging inversion concepts.

It is not appropriate when station-to-station variation is caused by lateral
geology rather than noise or near-surface distortion.

MT, CSAMT, And TEM Assumptions
------------------------------

The 1-D solvers share the same model container but represent different field
methods.

.. list-table::
   :header-rows: 1
   :widths: 24 26 26 24

   * - Solver
     - Source assumption
     - Sampled axis
     - Main output
   * - ``MT1DForward``
     - Natural plane wave
     - Frequency :math:`f`
     - Impedance, :math:`\rho_a`, phase.
   * - ``CSAMT1DForward``
     - Controlled-source AMT with far-field MT-like behaviour and optional
       near-field correction.
     - Frequency :math:`f`
     - Impedance, :math:`\rho_a`, phase.
   * - ``TEM1DForward``
     - Central-loop step-off transient.
     - Time :math:`t`
     - :math:`\partial B_z / \partial t`.

For MT and CSAMT, angular frequency is:

.. math::

   \omega = 2 \pi f.

For a 1-D MT response, apparent resistivity is computed from surface
impedance:

.. math::

   \rho_a = \frac{|Z|^2}{\mu_0 \omega}.

The impedance phase is:

.. math::

   \phi = \tan^{-1}\left(\frac{\operatorname{Im}(Z)}
                              {\operatorname{Re}(Z)}\right).

TEM is different. It observes a transient decay after transmitter current
changes. A useful scale estimate is:

.. math::

   z(t) \propto \sqrt{\frac{\rho t}{\mu_0}},

where later times and more resistive earth tend to sample greater depths.
This is only a sensitivity scale, not a direct depth conversion.

2-D Forward Modelling
---------------------

In a 2-D model, resistivity varies with horizontal distance and depth but is
constant along strike:

.. math::

   \rho(\mathbf{x}) = \rho(x, z).

:class:`pycsamt.forward.Grid2D` stores cell widths, depths, resistivity, and
station positions. :class:`pycsamt.forward.MT2DForward` solves both TE and TM
polarization modes on the grid.

.. code-block:: python
   :linenos:

   from pycsamt.forward import Grid2D, MT2DForward

   grid = Grid2D.with_anomaly(
       bg_rho=500.0,
       anomaly_rho=5.0,
       anomaly_bounds=(2000.0, 6000.0, 300.0, 1500.0),
       nx=50,
       nz=35,
       x_max=10_000.0,
       z_max=6000.0,
       n_pad=8,
       n_stations=16,
   )

   response = MT2DForward(
       freqs=[1.0, 10.0, 100.0],
       grid=grid,
       verbose=False,
   ).run()

   features = response.to_feature_array(mode="both")

The TE mode uses electric field parallel to strike. In pyCSAMT's 2-D
response object, this corresponds to ``ZXY``. The TM mode uses magnetic field
parallel to strike and corresponds to ``ZYX`` by convention.

.. list-table::
   :header-rows: 1
   :widths: 18 32 50

   * - Mode
     - Response component
     - Interpretation note
   * - TE
     - ``zxy``, ``rho_a_te``, ``phase_te``
     - Often sensitive to conductive structures connected along strike.
   * - TM
     - ``zyx``, ``rho_a_tm``, ``phase_tm``
     - Often responds strongly to charge accumulation across resistivity
       contrasts.

The TE/TM distinction is important for inversion. A 2-D inversion that fits
one mode but not the other may be telling you that the dimensionality
assumption, strike direction, static shift correction, or error model needs
review.

Quasi-3-D Forward Modelling
---------------------------

The quasi-3-D solver is designed for survey-scale synthetic experiments and
AI training where station layouts are two-dimensional but a full production
3-D solver would be too expensive.

:class:`pycsamt.forward.MT3DForward` extracts orthogonal 2-D slices from a
:class:`pycsamt.forward.Grid3D`, solves 2-D responses, and assembles an
approximate tensor response. It can return ``xy``, ``yx``, ``xx``, and ``yy``
apparent resistivity and phase arrays.

.. code-block:: python
   :linenos:

   from pycsamt.forward import Grid3D, MT3DForward

   grid = Grid3D.block_anomaly(
       bg_rho=500.0,
       anomaly_rho=20.0,
       bounds=(2000.0, 6000.0, 2000.0, 6000.0, 300.0, 1500.0),
       nx=20,
       ny=20,
       nz=15,
       x_max=8000.0,
       y_max=8000.0,
       z_max=4000.0,
       nx_stations=5,
       ny_stations=5,
   )

   response = MT3DForward(
       freqs=[1.0, 10.0, 100.0],
       grid=grid,
       method="quasi3d",
       verbose=False,
   ).run()

   x = response.to_feature_array(components="xy_yx")

The quasi-3-D result is useful for:

* generating spatially varying synthetic responses;
* testing station-layout effects;
* building AI datasets with multi-station features;
* checking tensor-style plotting and response maps;
* developing workflows before moving to ModEM or MARE2DEM.

It should not be presented as a full 3-D Maxwell solver result for final
geological interpretation. For production 3-D inversion and reporting, use an
appropriate external backend such as :doc:`../models/modem` or
:doc:`../models/mare2dem`.

Frequency, Period, And Depth Sensitivity
----------------------------------------

Frequency-domain EM responses are commonly plotted against period
:math:`T = 1/f`. Lower frequencies, or longer periods, generally sample
greater depths. A useful plane-wave skin-depth scale is:

.. math::

   \delta \approx 503 \sqrt{\frac{\rho}{f}},

where :math:`\delta` is in metres when :math:`\rho` is in :math:`\Omega m`
and :math:`f` is in Hz.

This is only a scale estimate. Actual sensitivity depends on:

* dimensionality;
* source geometry;
* resistivity contrast;
* station spacing;
* topography and near-surface structure;
* error floors and noise;
* which response component is used.

Still, the skin-depth estimate is useful when choosing frequency ranges. If
the lowest frequency is too high, a deep target may not affect the response
enough to recover. If the highest frequency is too low, shallow structure may
be poorly resolved.

Response Quantities
-------------------

Forward responses are returned in physical arrays and can also be flattened
into feature arrays for AI and inversion experiments.

.. list-table::
   :header-rows: 1
   :widths: 24 35 41

   * - Quantity
     - Where it appears
     - Use
   * - Complex impedance ``z``
     - ``ForwardResponse.z``, ``zxy``, ``zyx``, tensor components
     - Native frequency-domain EM response.
   * - Apparent resistivity ``rho_a``
     - 1-D, 2-D, quasi-3-D response objects
     - Easy to inspect across period and station.
   * - Phase
     - 1-D, 2-D, quasi-3-D response objects
     - Adds independent information about conductivity contrasts and
       inductive behaviour.
   * - ``dBz_dt``
     - ``TEM1DForward``
     - Time-domain transient decay.
   * - Feature array
     - ``to_array`` or ``to_feature_array``
     - Machine-learning or inversion-ready data vector.

For MT/CSAMT training data, a common feature vector is:

.. code-block:: text
   :linenos:

   [log10(rho_a_0), ..., log10(rho_a_n),
    phase_0, ..., phase_n]

For 2-D and quasi-3-D responses, the feature array is usually organized by
station, so the result has shape ``(n_stations, n_features)``.

Synthetic Priors
----------------

A synthetic dataset is only as useful as its model prior. If the prior
contains unrealistic resistivity ranges, layer counts, depths, or lateral
correlation, the trained model may perform well on synthetic validation data
and poorly on field data.

``LayeredModel`` supports several constructors:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Constructor
     - Purpose
   * - ``LayeredModel.random``
     - Draws random layer resistivities and thicknesses from broad bounds.
   * - ``LayeredModel.blocky``
     - Builds sharper layered contrasts.
   * - ``LayeredModel.smooth``
     - Builds smoother vertical variation.
   * - ``LayeredModel.from_geology``
     - Uses named geological priors such as sedimentary, crystalline,
       geothermal, marine, and permafrost settings.

Use geological priors when the goal is AI training for a known target class.
Use broader random priors when the goal is method exploration or robustness
testing.

Noise And Realism
-----------------

Noise-free synthetic data are useful for checking solver behaviour, but they
are usually too clean for training robust AI models or evaluating field
workflows. pyCSAMT includes Gaussian, multiplicative, and field-realistic
noise models so synthetic responses can include controlled uncertainty.

.. list-table::
   :header-rows: 1
   :widths: 26 36 38

   * - Noise model
     - Meaning
     - Good use
   * - ``GaussianNoise``
     - Adds independent perturbations to response values.
     - Baseline uncertainty tests.
   * - ``MultiplicativeNoise``
     - Applies log-space relative perturbations.
     - Responses spanning orders of magnitude.
   * - ``FieldRealisticNoise``
     - Adds frequency-dependent noise, dead-band effects, and power-line
       style contamination.
     - MT/CSAMT training data and field-like stress tests.

Noise should be documented as part of the experiment:

* noise type;
* noise level;
* random seed;
* whether phase was included;
* any response clipping or log transform;
* the final saved dataset path.

Synthetic Recovery
------------------

The most useful forward experiment is often a synthetic recovery test. The
workflow is:

#. create a known model;
#. compute a clean response;
#. add controlled noise;
#. invert the synthetic response;
#. compare recovered and true models;
#. repeat with different frequencies, noise levels, and regularization.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.forward import FieldRealisticNoise, LayeredModel, MT1DForward

   true_model = LayeredModel(
       resistivity=[200.0, 20.0, 800.0],
       thickness=[400.0, 1200.0],
   )

   freqs = np.logspace(-3, 4, 40)
   clean = MT1DForward(freqs).run(true_model)
   noisy = FieldRealisticNoise(base_level=0.05).apply(clean, seed=0)

   features = noisy.to_array(include_phase=True)

This test does not prove that a field inversion is correct. It proves that,
under controlled assumptions, the inversion workflow can recover a known
model from a response similar to the one being studied.

Common Failure Modes
--------------------

Forward modelling is often where mistakes are easiest to catch.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Symptom
     - Likely cause
   * - Apparent resistivity is flat for every model.
     - Frequency range may not sample the target depth, or the model contrast
       is too weak.
   * - Phase looks unrealistic or unstable.
     - Sign convention, component choice, or numerical resolution may need
       review.
   * - 2-D stations show no lateral variation.
     - The anomaly may be too deep, too small, outside the station aperture,
       or masked by grid/padding choices.
   * - Quasi-3-D tensor components look over-interpreted.
     - Remember that the method is approximate and assembled from 2-D
       profiles.
   * - AI training accuracy is high but field performance is poor.
     - Synthetic priors, noise model, or feature layout may not match field
       conditions.
   * - Inversion recovers the wrong model in a synthetic test.
     - The problem may be non-unique, under-sampled, over-regularized, or
       using incompatible data errors.

Good Practice Checklist
-----------------------

Before using a forward response in inversion, training, or reporting, verify:

* the dimensionality matches the question;
* resistivity values and thicknesses are physically plausible;
* sampled frequencies or times cover the expected target sensitivity;
* station spacing resolves the lateral scale of interest;
* padding and grid resolution are adequate;
* response units and component labels are clear;
* noise settings are recorded and reproducible;
* feature arrays include the intended quantities;
* plots have been reviewed for obvious physical or numerical issues.

Next Steps
----------

* :doc:`configuration` explains source-of-truth configuration files.
* :doc:`solvers_and_grids` gives concrete solver and grid examples.
* :doc:`synthetic_datasets` explains dataset generation and noise.
* :doc:`plotting` shows how to inspect responses and models.
* :doc:`forward_to_inversion` connects forward tests to inversion workflows.
* :doc:`../theory/inversion_concepts` provides the inversion background.
