.. _ai_inversion_concepts:

AI inversion concepts
=====================

AI inversion estimates an earth model from electromagnetic observations with
a learned or differentiable computational model. In pyCSAMT, this includes
three related approaches:

* **supervised surrogate inversion**, which learns from synthetic
  response–model pairs;
* **physics-informed inversion**, which optimizes model parameters through a
  differentiable forward-physics loss without labelled earth models;
* **hybrid inversion**, which uses an AI estimate as a starting point or prior
  for subsequent physics-based refinement.

These approaches can make repeated inversion fast, provide useful screening
models, and support uncertainty experiments. They do not remove the
non-uniqueness of electromagnetic inversion. Their results remain conditional
on data quality, parameterization, physics, training coverage, regularization,
and validation evidence.

.. admonition:: Core principle
   :class: important

   A neural network does not learn "the subsurface." It learns a relationship
   defined by its inputs, targets, forward simulator, parameter ranges, loss,
   architecture, and training examples. A field prediction is defensible only
   to the extent that this learned relationship represents the survey.

Where AI inversion fits
-----------------------

The principal EM workflow is:

.. figure:: ../../images/user_guide/ai_inversion/workflow.png
   :alt: Relationship among field data, forward modeling, AI inversion, validation, and interpretation
   :align: center
   :width: 92%

   Forward modeling supplies controlled response–model relationships. AI
   inversion predicts candidate models. Response reconstruction, uncertainty,
   classical baselines, and independent evidence determine whether those
   candidates are suitable for interpretation.

The stages have distinct roles:

Forward modeling
   Predicts observations :math:`\mathbf d` from an earth model
   :math:`\mathbf m` through :math:`\mathbf d=F(\mathbf m)`.

Classical inversion
   Searches for a model by repeatedly evaluating physics and minimizing an
   objective such as data misfit plus regularization.

Supervised AI inversion
   Learns an approximation :math:`G_\theta` to the inverse relationship from
   examples and predicts :math:`\hat{\mathbf m}=G_\theta(\mathbf d)`.

Physics-informed inversion
   Uses differentiable physics during optimization, often minimizing
   :math:`\|F(\mathbf m_\theta)-\mathbf d_{obs}\|` plus constraints.

Interpretation
   Connects the reviewed resistivity model to geology or hydrogeology. It is
   not part of network prediction and requires additional evidence.

Why EM inversion remains difficult
----------------------------------

Electromagnetic inversion is nonlinear and ill posed. Different resistivity
models can produce responses that are indistinguishable within measurement
error. Sensitivity changes with period, depth, component, geometry, and
conductivity structure. Finite bandwidth limits depth resolution, while noise,
distortion, source effects, and dimensionality violations introduce additional
ambiguity.

Classical inversion manages this with regularization, error models, starting
models, constraints, and model appraisal. Supervised AI manages it implicitly
through the training distribution, target representation, loss, architecture,
and data augmentation. Those choices are forms of prior information even when
they are not called regularization.

Consequently:

* a single response need not have a single correct target model;
* a low supervised test error can coexist with field failure;
* a smooth prediction may reflect training targets rather than geology;
* network confidence can be high outside the supported domain;
* fast inference changes computational cost, not physical information content.

Three AI inversion families
---------------------------

Supervised surrogate inversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The supervised workflow samples earth models, calculates synthetic responses,
and fits a network to recover the sampled parameters. pyCSAMT provides:

* :class:`pycsamt.ai.inversion.EMInverter1D`;
* :class:`pycsamt.ai.inversion.EMInverter2D`;
* :class:`pycsamt.ai.inversion.GCNInverter3D`;
* :class:`pycsamt.ai.inversion.JointInverter`;
* :class:`pycsamt.ai.inversion.EnsembleInverter`.

Advantages include rapid inference after training, consistent execution across
many stations, and direct access to synthetic ground truth during development.
The principal risk is the simulation-to-field domain gap.

Physics-informed inversion
~~~~~~~~~~~~~~~~~~~~~~~~~~

pyCSAMT exposes:

* :class:`pycsamt.ai.inversion.PINNInverter1D`;
* :class:`pycsamt.ai.inversion.PINNInverter2D`;
* :class:`pycsamt.ai.inversion.PINNInverter3D`.

These approaches optimize through a physics loss and do not require labelled
training models in the same way as supervised inversion. They still depend on
the selected forward operator, model parameterization, regularization weights,
boundary assumptions, optimizer, and stopping rule. "Physics-informed" does
not mean exact physics or unique recovery.

Hybrid inversion
~~~~~~~~~~~~~~~~

pyCSAMT provides:

* :class:`pycsamt.ai.inversion.HybridInverter1D`;
* :class:`pycsamt.ai.inversion.HybridInverter2D`;
* :class:`pycsamt.ai.inversion.HybridInverter3D`.

A hybrid workflow combines learned speed with physics-based correction. For
example, an AI result can initialize an iterative refinement. This can reduce
runtime or starting-model sensitivity, but refinement cannot rescue an
incompatible parameterization or missing physics automatically.

Choosing 1-D, 2-D, or 3-D
-------------------------

Dimension describes the model and information-sharing assumption, not merely
the shape of an output array.

.. list-table::
   :header-rows: 1
   :widths: 15 27 27 31

   * - Level
     - Representation
     - Appropriate context
     - Principal risk
   * - 1-D
     - One layered model per sounding.
     - Locally layered earth, sparse surveys, rapid station screening, or a
       baseline before higher-dimensional analysis.
     - Lateral structure is mapped into misleading layer parameters.
   * - 2-D
     - A depth–station section predicted from an ordered profile panel.
     - A line survey whose strike and dimensionality evidence support a 2-D
       approximation.
     - Learned lateral continuity is confused with 2-D EM physics.
   * - Graph 3-D
     - Layered station parameters share information through a spatial graph.
     - Multi-line or areal surveys with reviewed coordinates and meaningful
       neighborhood relationships.
     - Graph smoothing is presented as a full 3-D numerical EM inversion.

:class:`~pycsamt.ai.inversion.GCNInverter3D` is graph-context inversion. It
predicts layered parameters at spatial stations and shares information through
adjacency. The label 3-D should not be interpreted as proof that a full 3-D EM
forward operator was used in every supervised training example.

Select dimension using strike, tensor dimensionality, tipper behavior, survey
layout, target geometry, station spacing, and classical response evidence—not
because a higher number sounds more advanced.

Model parameterization
----------------------

The network can only recover parameters represented in its output.

Layered 1-D target
~~~~~~~~~~~~~~~~~~

A common target contains ``n_layers`` resistivities and ``n_layers - 1``
thicknesses:

.. code-block:: text

   y = [log10(rho_1), ..., log10(rho_L),
        log10(h_1), ..., log10(h_(L-1))]

Resistivity is commonly expressed in ohm metres and thickness in metres before
log10 transformation. The bottom layer is a half-space and has no thickness.
Increasing layer count increases flexibility but also ambiguity and training
difficulty.

Fixed-grid 2-D target
~~~~~~~~~~~~~~~~~~~~~

A profile inverter may predict ``log10(rho)`` on a fixed
``(n_depth, n_stations)`` grid. Depth discretization is then part of the model
prior. A visually sharp boundary cannot be more resolved than the information
in the responses and training examples.

Graph target
~~~~~~~~~~~~

Graph inversion commonly predicts the layered target at every station. A
later volume or depth slice interpolates these predictions spatially. Separate
network output, graph smoothing, and post-prediction gridding in reports.

Parameterization decisions include:

* linear versus logarithmic variables;
* bounds and probability distributions;
* number and thickness of layers or depth cells;
* fixed versus predicted interfaces;
* anisotropy or isotropy;
* components and modes represented;
* topography, source, and distortion parameters;
* whether properties vary independently or through geological rules.

Data representation
-------------------

The input representation determines what the network can use.

Possible features include:

* log10 apparent resistivity;
* phase;
* real and imaginary impedance components;
* TE/TM or tensor-component combinations;
* TDEM decay values;
* frequency or period coordinates;
* masks, errors, and quality indicators;
* station coordinates or adjacency.

The same observed data can yield different learned problems depending on
feature choice. Apparent resistivity and phase are derived from impedance;
treating all of them as independent channels can duplicate information.

pyCSAMT provides public bridge utilities rather than requiring invented field
feature helpers:

* :class:`pycsamt.ai.inversion.SiteObs1D` and
  :class:`pycsamt.ai.inversion.SiteObs2D`;
* :func:`pycsamt.ai.inversion.sites_to_obs_1d`;
* :func:`pycsamt.ai.inversion.sites_to_obs_2d`;
* :func:`pycsamt.ai.inversion.sites_to_features_1d`;
* :func:`pycsamt.ai.inversion.sites_to_panel_2d`;
* :func:`pycsamt.ai.inversion.sites_to_coords_3d`;
* :func:`pycsamt.ai.inversion.obs_to_features_1d`.

These utilities help enforce the contract between :class:`pycsamt.site.Sites`
and AI arrays. Their outputs still need review for station order, component,
frequency coverage, masking, and finite values.

The feature contract
--------------------

Training and inference must use the same:

* feature names and order;
* units and transformations;
* frequency or period grid;
* component convention;
* interpolation method;
* missing-value and mask behavior;
* normalization statistics;
* station ordering and padding policy;
* coordinate units and adjacency rules.

A checkpoint without its feature contract is incomplete. A matrix with the
right shape but wrong channel order can produce plausible yet meaningless
predictions.

Training distribution as prior
------------------------------

The synthetic model distribution defines which earth structures the network
expects. It therefore acts as an implicit prior.

.. figure:: ../../images/user_guide/ai_inversion/training_distribution.png
   :alt: Comparison of synthetic training coverage and field observations
   :align: center
   :width: 90%

   Field responses outside the synthetic envelope require extrapolation.
   Network behavior there is not validated by ordinary held-out synthetic
   performance.

Design choices include:

Geology
   Resistivity ranges, layer counts, thicknesses, correlations, interfaces,
   faults, lateral structures, anisotropy, and rare targets.

Survey
   Frequency range, sampling density, station count, spacing, line layout,
   components, source geometry, and coordinate uncertainty.

Noise and nuisance effects
   Heteroscedastic errors, missing bands, outliers, static shift, source
   effects, cultural noise, distortion, and processing variability.

Physics
   Forward solver, dimensionality, boundary conditions, mesh accuracy, and
   approximations.

Class balance
   Frequency of common backgrounds versus rare conductive or resistive targets.

Sampling parameters independently from broad uniform ranges is convenient but
may create geologically impossible examples. Conversely, an overly narrow
geological generator can produce excellent test metrics and poor field
coverage.

Simulation-to-field domain gap
------------------------------

The domain gap is the difference between synthetic training examples and field
observations. Sources include:

* incomplete forward physics;
* wrong dimensionality;
* underestimated noise and missing data;
* static shift and galvanic distortion;
* source and near-field effects;
* topography and anisotropy;
* coordinate and instrument errors;
* processing differences;
* geology outside sampled priors.

Good synthetic validation does not measure this gap. Diagnose it by comparing
feature distributions, response envelopes, residual patterns, classical
inversions, and independent field evidence. Fine-tuning on field labels can
help only when trustworthy labels exist and leakage is controlled.

Supervised learning objective
-----------------------------

A supervised inverter minimizes a model-space loss over examples:

.. math::

   \mathcal L_{model}(\theta)
   = \frac{1}{N}\sum_i
   \ell\!\left(G_\theta(\mathbf d_i),\mathbf m_i\right).

This rewards resemblance to the selected synthetic target. It does not
guarantee that the predicted model reconstructs the field response. A stronger
workflow also evaluates a physics-space diagnostic:

.. math::

   \mathcal L_{data}
   = \left\|W_d\left(F(\hat{\mathbf m})-mathbf d_{obs}\right)\right\|^2.

The forward operator, error weights, components, and residual space must be
stated. An unweighted RMS of log apparent resistivity is not equivalent to a
full complex-impedance likelihood.

Physics-informed objective
--------------------------

A PINN-style objective can combine data fit and regularization:

.. math::

   \mathcal L
   = \mathcal L_{data}
   + \lambda_v\mathcal R_v
   + \lambda_l\mathcal R_l
   + \lambda_g\mathcal R_g,

where vertical, lateral, or graph regularizers apply according to dimension.
The weights determine the balance between fit and structure. They must be
treated as inversion parameters and tested, not hidden as neural-network
details.

Architectures encode assumptions
--------------------------------

Fully connected network
   Treats the feature vector globally. It is simple but does not explicitly
   encode local frequency structure.

1-D convolution
   Learns local patterns along an ordered feature or frequency axis. Channel
   arrangement and padding affect meaning.

Residual network
   Uses skip connections to train deeper transformations and often provides a
   strong station-level baseline.

U-Net
   Combines local convolution and multiscale skip connections for profile
   panels. Its receptive field encourages lateral continuity but does not
   implement a classical 2-D EM solver.

Graph convolutional network
   Propagates information along an adjacency graph. Radius, coordinate system,
   normalization, and graph connectivity become inversion assumptions.

Architecture comparison is valid only when datasets, splits, preprocessing,
target scale, training budget, and selection procedure are controlled.

Training, validation, and test separation
-----------------------------------------

Use distinct roles:

Training set
   Updates network parameters.

Validation set
   Selects epochs, hyperparameters, and checkpoints.

Calibration set
   Fits conformal or posterior calibration without changing the base network.

Synthetic test set
   Estimates final performance on unseen synthetic cases.

Field validation set
   Tests transfer using observations or independent constraints not used in
   training, tuning, or calibration.

Random row splitting can leak nearly identical models or profiles across sets.
Prefer group splitting by geological scenario, profile family, random generator
seed, or survey realization. For spatial field data, nearby stations are not
independent test samples when the model shares spatial information.

.. figure:: ../../images/user_guide/ai_inversion/training_convergence.png
   :alt: Training and validation loss showing checkpoint selection
   :align: center
   :width: 82%

   Training loss measures optimization on seen examples. Validation divergence
   can indicate overfitting, but parallel curves do not prove field transfer.

Evaluation hierarchy
--------------------

Evaluate at several levels:

Model-space metrics
   Error in log resistivity, thickness, boundary depth, or complete sections on
   synthetic cases with known truth.

Response-space metrics
   Difference between observed responses and responses recomputed from the
   prediction.

Structural metrics
   Boundary location, anomaly continuity, target detection, volume overlap, or
   graph consistency where scientifically relevant.

Calibration metrics
   Interval coverage, reliability, and sharpness on held-out calibration/test
   data.

Out-of-distribution diagnostics
   Distance or coverage relative to synthetic inputs and latent
   representations.

Field baselines
   Bostick-style transforms, classical 1-D/2-D/3-D inversion, boreholes,
   mapped geology, and other geophysics.

Decision metrics
   Whether target ranking, boundary depth, or classification remains correct
   under accepted uncertainty and scenarios.

No single metric is sufficient. A low model-space error can hide response
mismatch, and a low response residual can correspond to the wrong model because
the inverse problem is non-unique.

.. figure:: ../../images/user_guide/ai_inversion/predicted_section.png
   :alt: Predicted AI resistivity section displayed with residual and review diagnostics
   :align: center
   :width: 92%

   Review predictions together with coverage, response reconstruction,
   residuals, uncertainty, and classical or geological evidence.

Uncertainty concepts
--------------------

Useful distinctions are:

Aleatoric uncertainty
   Variation associated with observation noise or irreducible ambiguity.

Epistemic uncertainty
   Uncertainty in learned parameters due to finite or incomplete training data.

Distributional uncertainty
   Risk that the field input comes from a different distribution than training.

Inverse non-uniqueness
   Multiple earth models fit the observations. A deterministic network can
   collapse these possibilities into one conditional estimate.

Structural uncertainty
   Wrong dimension, parameterization, forward physics, architecture, or
   geological assumptions.

pyCSAMT exposes :class:`pycsamt.ai.inversion.EnsembleInverter`,
:class:`pycsamt.ai.inversion.ConformalPredictor`,
:class:`pycsamt.ai.inversion.PosteriorCalibrator`, and MC-dropout prediction for
the graph inverter. Each covers only part of this taxonomy.

.. figure:: ../../images/user_guide/ai_inversion/uncertainty_profile.png
   :alt: Predictive interval widening with depth
   :align: center
   :width: 84%

   Predictive spread should be interpreted conditionally on the ensemble,
   calibration set, or dropout model that produced it.

Conformal coverage applies under exchangeability assumptions between
calibration and future examples. Synthetic calibration does not guarantee the
same coverage on shifted field data. Narrow intervals can be confidently
wrong when structural uncertainty is omitted.

Configuration objects
---------------------

:class:`pycsamt.ai.inversion.InversionConfig` describes an AI inverter, while
:class:`pycsamt.ai.inversion.RunConfig` coordinates dataset and training
configuration. Configuration-first work is preferable to undocumented notebook
state:

.. code-block:: python

   from pycsamt.ai.inversion import RunConfig

   RunConfig.write_template("ai_inversion.yml")
   config = RunConfig.from_file("ai_inversion.yml")
   config.validate()
   print(config.summary())

Validation checks internal configuration consistency. It cannot establish that
the selected ranges, noise, physics, or architecture represent a particular
field survey.

The Sites bridge
----------------

Use the canonical loader and public bridge utilities:

.. code-block:: python

   from pycsamt.emtools._core import ensure_sites
   from pycsamt.ai.inversion import (
       sites_to_features_1d,
       sites_to_obs_1d,
   )

   sites = ensure_sites(
       "data/AMT/WILLY_DATA/L18",
       recursive=True,
       verbose=0,
   )

   observations = sites_to_obs_1d(sites)
   X, frequencies_hz, station_names = sites_to_features_1d(
       sites,
       comp="xy",
       n_freqs=32,
   )

Before using these arrays, inspect their documented return type and shape in
the installed version, station names, frequency grid, components, and missing
data. Do not invent a feature helper whose transformation differs from the
checkpoint contract.

Agents and automation
---------------------

The agents described in :doc:`agents` orchestrate standard workflows around
the lower-level objects. They are useful for screening and repeatable task
execution, but their defaults encode specific synthetic-data and diagnostic
choices. A successful :class:`pycsamt.agents.AgentResult` means the programmed
workflow completed; it is not a scientific acceptance decision.

Optional LLM text is separate from numerical inversion. It may summarize
structured outputs, but it does not alter the recovered model and must be
reviewed before reporting.

Pretrained checkpoints
----------------------

A checkpoint is usable only with its full model contract:

* architecture and layer/depth parameterization;
* solver and method;
* feature order, units, frequency grid, and normalization;
* training priors and noise distribution;
* dataset split and performance metrics;
* software/backend version;
* checkpoint checksum and model-card limitations.

Registry metadata can be inspected through
:func:`pycsamt.ai._zoo.list_pretrained` or
:class:`pycsamt.agents.ModelZooAgent`. Checkpoint availability is separate
from registry availability. Never claim pretrained inference if execution
silently fell back to new training.

Reproducibility and provenance
------------------------------

A reproducible AI inversion record includes:

.. code-block:: text

   ai_inversion_run/
   ├── manifest.yml
   ├── environment/
   │   └── dependencies.txt
   ├── configuration/
   │   ├── dataset.yml
   │   ├── model.yml
   │   └── training.yml
   ├── datasets/
   │   ├── metadata.json
   │   └── split_indices.npz
   ├── checkpoints/
   │   ├── best_model.*
   │   └── checksums.sha256
   ├── predictions/
   │   ├── field_predictions.npz
   │   └── station_order.csv
   ├── diagnostics/
   │   ├── training_history.csv
   │   ├── response_residuals.csv
   │   └── coverage_report.csv
   ├── figures/
   └── review/
       └── model_card.md

Record random seeds, backend, hardware, precision, preprocessing code, dataset
generator version, station ordering, graph construction, model selection, and
all failed or excluded cases. A checkpoint alone is not reproducible.

Scientific acceptance framework
-------------------------------

Before field interpretation, confirm:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Question
     - Required evidence
   * - Is the dimension defensible?
     - Strike, tensor, tipper, survey geometry, and classical diagnostics.
   * - Is the parameterization adequate?
     - Target depth, expected structures, layer/grid resolution, and bounds.
   * - Is the feature contract exact?
     - Channel order, units, grid, masks, normalization, and station order.
   * - Does training cover the field?
     - Response envelopes, parameter support, nuisance effects, and explicit
       out-of-distribution tests.
   * - Was evaluation independent?
     - Leakage-resistant training, validation, calibration, synthetic test,
       and field-validation roles.
   * - Does the model reconstruct data?
     - Forward responses, error-aware residuals, station/frequency/component
       patterns, and failure cases.
   * - Is uncertainty calibrated?
     - Coverage and sharpness on appropriate held-out data plus distributional
       and structural limitations.
   * - Do baselines agree?
     - Classical inversion and independent geological or borehole evidence,
       with mismatches reported.
   * - Is the result reproducible?
     - Dataset, configuration, code, environment, checkpoint, checksum,
       predictions, and review record.

Common conceptual mistakes
--------------------------

Avoid these misunderstandings:

* assuming synthetic ground truth is the unique inverse solution;
* treating a neural network as free of regularization or prior assumptions;
* equating training interpolation with field generalization;
* selecting 2-D or 3-D solely for visual sophistication;
* confusing U-Net continuity or graph smoothing with numerical EM physics;
* checking array shape but not feature semantics;
* using random splits that leak related synthetic profiles;
* selecting a checkpoint on the test set;
* reporting model-space error without response reconstruction;
* interpreting dropout or ensemble spread as total uncertainty;
* assuming physics-informed optimization eliminates non-uniqueness;
* treating fast prediction as increased depth resolution;
* using an undocumented pretrained checkpoint;
* accepting an automated or LLM narrative as geological validation.

Next steps
----------

Continue in this order:

* :doc:`data_preparation` to define synthetic and field data contracts;
* :doc:`model_selection` to choose dimension and architecture;
* :doc:`training` to fit and preserve a model correctly;
* :doc:`validation` to establish acceptance evidence;
* :doc:`inference` to apply an approved checkpoint;
* :doc:`uncertainty` to assess predictive calibration and domain shift;
* :doc:`pinn_2d` for physics-informed profile inversion;
* :doc:`agents` for standard orchestration;
* :doc:`reporting` for model cards and release packages.

Documentation figures
---------------------

The figures on this page are generated by
``docs/scripts/generate_ai_inversion_figures.py``. They illustrate validation
logic and do not represent a field result or a production training run.
