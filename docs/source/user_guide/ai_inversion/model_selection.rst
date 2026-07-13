.. _ai_inversion_model_selection:

AI model selection
==================

Model selection chooses the scientific representation, learning strategy, and
architecture used for AI inversion. It is not simply a search for the network
with the lowest validation loss. The selected model must match survey physics,
geometry, data volume, target parameterization, deployment constraints, and the
evidence required for field acceptance.

pyCSAMT supports supervised 1-D, profile 2-D, graph-context 3-D, joint,
ensemble, physics-informed, and hybrid workflows. Each answers a different
question and encodes different priors.

.. admonition:: Prefer the simplest defensible model
   :class: important

   Added architectural complexity is justified only when it improves a
   predeclared scientific or operational criterion on independent data. A more
   elaborate network can hide stronger assumptions, increase leakage risk, and
   make deployment less reproducible without adding recoverable information.

Selection workflow
------------------

#. define the decision, output, and acceptance metrics;
#. determine defensible physical dimension from the survey;
#. freeze the feature and target contracts;
#. establish non-AI and simple AI baselines;
#. identify candidate model families rather than arbitrary hyperparameters;
#. define leakage-resistant training, validation, calibration, and test roles;
#. train candidates under comparable budgets;
#. compare model-space, response-space, calibration, and field metrics;
#. test robustness, domain shift, compute, and persistence;
#. select once using the declared rule and confirm on an untouched test set;
#. document rejected candidates and the reason for the final choice.

1. Define the selection objective
---------------------------------

Different objectives can select different models. Examples include:

* minimize median log-resistivity error on unseen synthetic geology;
* recover conductor-top depth within an operational tolerance;
* minimize error-weighted field response residuals;
* maintain calibrated 90% intervals with useful sharpness;
* rank drilling targets consistently across accepted scenarios;
* run within memory and latency limits on the deployment machine;
* remain stable under missing frequency bands or station dropout.

Declare a primary metric, secondary constraints, and unacceptable failure
modes before comparing candidates. Do not select using whichever metric looks
best after training.

2. Decide dimension before architecture
---------------------------------------

Dimension is a geophysical assumption.

Choose 1-D when
~~~~~~~~~~~~~~~

* local geology is approximately layered;
* station-by-station screening is the goal;
* survey spacing is too sparse for stable lateral learning;
* no reviewed profile or areal geometry is available;
* a transparent baseline is needed.

Use :class:`pycsamt.ai.inversion.EMInverter1D`, a 1-D PINN, or a 1-D hybrid.

Choose profile 2-D when
~~~~~~~~~~~~~~~~~~~~~~~

* stations form a meaningful ordered line;
* tensor, strike, and dimensionality evidence support a 2-D approximation;
* lateral continuity is central to the decision;
* training examples represent profile geometry and missingness;
* a fixed station and depth representation is acceptable.

Use :class:`pycsamt.ai.inversion.EMInverter2D`,
:class:`pycsamt.ai.inversion.PINNInverter2D`, or
:class:`pycsamt.ai.inversion.HybridInverter2D`.

Choose graph-context 3-D when
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* stations form several lines or an areal layout;
* projected coordinates and station order are reliable;
* neighborhood relationships carry useful information;
* training surveys represent graph layout and correlation scales;
* node-wise layered outputs are suitable.

Use :class:`pycsamt.ai.inversion.GCNInverter3D`,
:class:`pycsamt.ai.inversion.PINNInverter3D`, or
:class:`pycsamt.ai.inversion.HybridInverter3D`.

The GCN shares information among station nodes. It does not automatically mean
that supervised examples were generated with a full numerical 3-D EM solver.
Do not select it merely to label the product "3-D."

3. Compare model families
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 19 22 28 31

   * - Family
     - Data requirement
     - Strength
     - Principal limitation
   * - Supervised 1-D
     - Synthetic labelled soundings.
     - Fast, simple, portable checkpoint support, strong baseline.
     - Ignores lateral structure and depends on synthetic coverage.
   * - Supervised U-Net 2-D
     - Labelled profile panels and fixed section targets.
     - Learns multiscale lateral patterns over a complete line.
     - Continuity is architectural; current class lacks explicit public
       save/load parity with the 1-D class.
   * - Supervised GCN
     - Labelled multi-station surveys and fixed graph convention.
     - Uses irregular spatial relationships and supports dropout spread.
     - Sensitive to graph construction; pseudo-3-D training may use 1-D
       forward responses.
   * - Joint supervised
     - Aligned, correlated modalities with common model targets.
     - Can combine complementary sensitivity.
     - Alignment, modality imbalance, and cross-property assumptions increase
       complexity and leakage risk.
   * - Ensemble
     - Several fitted compatible estimators plus calibration data.
     - Measures member disagreement and supports calibrated intervals.
     - Multiplies training/storage cost and misses shared structural bias.
   * - PINN
     - Observations and differentiable physics; no labelled model targets.
     - Direct response-based optimization and explicit regularization.
     - Slower per survey; sensitive to optimizer, physics, and loss weights.
   * - Hybrid
     - AI initial estimate plus physics-refinement workflow.
     - Combines fast initialization with response-space refinement.
     - Two-stage provenance and failure behavior are more complex.

4. Choose the 1-D architecture
------------------------------

:class:`pycsamt.ai.inversion.EMInverter1D` supports ``"fcn"``, ``"cnn1d"``,
and ``"resnet"``:

.. code-block:: python

   from pycsamt.ai.inversion import EMInverter1D

   candidate = EMInverter1D(
       arch="resnet",
       n_layers=5,
       solver="mt1d",
       include_phase=True,
       log_thickness=True,
   )

FCN
~~~

A fully connected network treats the complete feature vector globally.

Prefer it when:

* the feature vector is short and fixed;
* a simple capacity baseline is needed;
* channel adjacency should not be assumed meaningful;
* training data are limited relative to larger models.

Its weakness is that it does not explicitly exploit ordered frequency-local
patterns.

CNN1D
~~~~~

A one-dimensional convolution shares filters along the arranged input axis.

Prefer it when:

* features have a meaningful local ordering;
* the same transformation should detect patterns across frequency;
* robustness to local shifts or patterns is useful.

Check the actual block layout. Concatenated resistivity and phase blocks create
a boundary that is not a physical frequency adjacency. Architecture benefit
depends on how the network implementation consumes that vector.

ResNet
~~~~~~

Residual connections support deeper mappings and stable gradient flow.

Prefer it as a strong default candidate when:

* enough synthetic examples cover the parameter space;
* FCN underfits important nonlinear structure;
* compute and checkpoint complexity remain acceptable.

Do not assume ResNet must win because it is deeper. Compare it against FCN and
CNN1D under the same data, split, optimization budget, and selection metric.

5. Choose layer count and target complexity
-------------------------------------------

For a layered model, output size is commonly ``2 * n_layers - 1``. Increasing
``n_layers`` adds resistivity and interface parameters, expanding ambiguity.

Select layer count using:

* available frequency/time bandwidth;
* expected number of resolvable electrical units;
* target depth and station spacing;
* sensitivity and classical inversion evidence;
* synthetic recovery by depth;
* stability across noise and starting distributions;
* the decision's required detail.

Compare candidate layer counts by response reconstruction and boundary
stability, not only target-vector MSE. A five-layer prediction can fit a
three-unit earth by splitting units artificially; a three-layer model can smear
a thin critical conductor.

Use a fixed layer count for straightforward supervised outputs. Variable-layer
targets from ``generate_dataset`` are NaN-padded and require a proven mask-aware
strategy; they do not automatically make the network variable-output.

6. Select solver and feature content
------------------------------------

For 1-D configuration, supported solver names are ``"mt1d"``, ``"csamt1d"``,
and ``"tem1d"``. The solver must agree with synthetic generation and field
representation.

Phase
~~~~~

Set ``include_phase=True`` when phase is available, quality-controlled, and
represented consistently in synthetic data. Phase adds information but also
exposes phase convention, wrapping, and noise differences.

Do not include a phase channel filled with constants merely to satisfy a model
trained with phase. Either use a validated missing-channel strategy or select
a model trained without phase.

Components
~~~~~~~~~~

Select ``xy``, ``yx``, TE/TM, or multiple tensor channels based on survey
orientation and dimensionality. More channels can improve identifiability only
when their physics, units, errors, and ordering are modeled consistently.

Frequency grid
~~~~~~~~~~~~~~

Grid extent controls sensitivity and deployment compatibility. A larger
``n_freqs`` increases input resolution and compute but cannot create
information absent from the field observations. Prefer a common range with
adequate coverage across required stations.

7. Select the 2-D U-Net configuration
-------------------------------------

:class:`pycsamt.ai.inversion.EMInverter2D` maps
``(n_components, n_freqs, n_stations)`` panels to
``(n_depth, n_stations)`` sections.

.. code-block:: python

   from pycsamt.ai.inversion import EMInverter2D

   candidate_2d = EMInverter2D(
       n_components=4,
       n_freqs=32,
       n_stations=48,
       n_depth=64,
       arch="unet",
       dropout=0.20,
   )

Selection variables include:

``n_components``
   Two channels can represent TE log resistivity and phase; four add TM log
   resistivity and phase in the public field bridge.

``n_freqs``
   Controls vertical input sampling and U-Net input height.

``n_stations``
   Fixes profile width. Field profiles with different station counts need a
   validated resampling, crop, pad, or separate model.

``n_depth``
   Fixes output depth sampling. More cells improve display resolution, not
   necessarily geophysical resolution.

``dropout``
   Regularizes training. Excessive dropout can underfit; insufficient dropout
   can overfit synthetic profile patterns.

``channels`` or U-Net depth
   Control capacity and receptive field. The implementation constrains pooling
   depth according to the smaller input dimension.

Select using complete-profile holdouts and boundary/response metrics. Never
split stations from one synthetic profile across training and validation.

8. Select graph architecture and adjacency
------------------------------------------

:class:`pycsamt.ai.inversion.GCNInverter3D` requires per-station features,
layer count, hidden widths, dropout, and graph context:

.. code-block:: python

   from pycsamt.ai.inversion import GCNInverter3D

   candidate_3d = GCNInverter3D(
       n_features=64,
       n_layers=5,
       hidden=(256, 128, 64),
       dropout=0.10,
   )

Adjacency is as important as network width. If neither adjacency nor
coordinates is passed, the implementation warns and uses identity adjacency,
removing inter-station coupling.

Compare radius or graph policies using:

* node degree and connected components;
* survey spacing distribution;
* known structural boundaries;
* recovery across correlation-length scenarios;
* station dropout robustness;
* performance on irregular layouts;
* uncertainty near isolated nodes and survey edges.

A radius that is too small reduces the GCN toward station-wise prediction. A
radius that is too large shares information across unrelated structures.
Select radius as a spatial prior in coordinate units, not as a generic neural
hyperparameter.

9. Decide whether joint inversion is justified
----------------------------------------------

:class:`pycsamt.ai.inversion.JointInverter` fuses two or more feature matrices
into one layered target:

.. code-block:: python

   from pycsamt.ai.inversion import JointInverter

   joint = JointInverter(
       n_features_list=(64, 25),
       n_layers=5,
       growth_rate=32,
       n_dense_layers=6,
       hidden_dim=256,
       dropout=0.20,
       log_thickness=True,
   )

Joint inversion is justified when modalities:

* observe compatible locations and support volumes;
* share a defensible common target or cross-property relationship;
* are aligned without label leakage;
* provide complementary rather than redundant information;
* have realistic joint synthetic examples;
* remain available under the deployment policy.

Compare against each modality alone. A fused network that improves aggregate
loss may rely almost entirely on one modality or learn spurious correlations.
Test missing-modality behavior explicitly.

10. Decide between supervised, PINN, and hybrid
-----------------------------------------------

Choose supervised when
~~~~~~~~~~~~~~~~~~~~~~

* a representative synthetic generator exists;
* repeated low-latency inference is important;
* target labels and parameterization are trustworthy;
* domain coverage can be tested.

Choose PINN when
~~~~~~~~~~~~~~~~

* labelled model generation is a bottleneck;
* differentiable physics matches the observations;
* per-survey optimization cost is acceptable;
* data-fit and regularization need direct control;
* convergence can be reviewed for every run.

Choose hybrid when
~~~~~~~~~~~~~~~~~~

* an AI prior provides a useful initialization;
* field response refinement is required;
* supervised domain shift is a concern;
* the extra two-stage complexity is justified by independent improvement.

Compare final response fit, structural stability, runtime, and failure rate.
Do not compare supervised inference latency alone with a PINN's complete
optimization and conclude that the latter is inferior without considering the
different task.

11. Decide whether an ensemble is required
------------------------------------------

:class:`pycsamt.ai.inversion.EnsembleInverter` trains or combines compatible
base estimators.

Use an ensemble when:

* uncertainty or stability is part of the decision;
* member diversity can be created through seeds, data, or architectures;
* calibration and test sets are large enough;
* multiplied training and storage cost are acceptable.

Member disagreement measures only differences among members. If every member
uses the same incorrect simulator or narrow prior, the ensemble can agree and
still be wrong. Compare raw spread, calibrated coverage, interval sharpness,
and field domain diagnostics.

12. Establish baselines
-----------------------

Every selection study should include:

Simple statistical baseline
   Mean target, nearest synthetic neighbor, or simple regression where
   appropriate. This reveals whether the network learns more than dataset bias.

Simple neural baseline
   Usually FCN or compact 1-D network.

Dimensional baseline
   Independent 1-D predictions before U-Net or graph context.

Physics baseline
   Forward reconstruction and a classical inversion appropriate to the survey,
   such as built-in 1-D, Occam2D, ModEM, or MARE2DEM.

Operational baseline
   Current project method, including runtime and reviewer effort.

A candidate should improve a declared criterion without unacceptable losses in
calibration, robustness, interpretability, or deployability.

13. Define a fair comparison protocol
-------------------------------------

Hold constant across candidates:

* dataset version and split indices;
* feature and target transformations;
* training and validation roles;
* number of repeated seeds;
* early-stopping policy;
* optimization budget where comparable;
* metric implementation;
* test and field evaluation sets;
* hardware or normalized resource accounting.

Architecture-specific tuning is legitimate, but the tuning budget should be
comparable and documented. Report the distribution across repeated runs rather
than only the best seed.

Use nested selection when possible:

* inner validation chooses architecture and hyperparameters;
* calibration fits uncertainty after model choice;
* untouched synthetic test estimates final in-domain performance;
* withheld field evidence assesses transfer.

Do not choose architecture, layer count, and reporting threshold on the same
test observations.

14. Compare multiple metric families
------------------------------------

Model-space
   Log-resistivity error, thickness error, boundary-depth error, section error,
   and performance by depth or geology.

Response-space
   Error-weighted impedance, apparent-resistivity, phase, or decay residuals
   after forward reconstruction.

Structural
   Conductor detection, boundary overlap, continuity, false targets, and edge
   behavior.

Uncertainty
   Coverage, interval width, reliability, ensemble spread, and failure under
   domain shift.

Robustness
   Noise, missing bands, station dropout, static shift, coordinate perturbation,
   and out-of-prior geology.

Operational
   Training time, inference time, peak memory, artifact size, persistence,
   backend availability, and reproducibility.

Field
   Classical response fit, borehole or geological agreement, and stability of
   the actual decision.

Aggregate averages can hide failure in deep layers, rare conductive targets,
or sparse stations. Report distributions and scenario-stratified results.

15. Consider persistence and deployment
---------------------------------------

Current public persistence support differs among classes:

* :class:`~pycsamt.ai.inversion.EMInverter1D` exposes ``save()`` and ``load()``;
* :class:`~pycsamt.ai.inversion.EnsembleInverter` saves and loads member
  checkpoints and ensemble metadata, but not attached uncertainty calibrators;
* the 2-D and graph inverter classes do not currently expose the same explicit
  public persistence pair.

If deployment requires a portable approved checkpoint, this can favor 1-D or a
project-tested persistence implementation. Do not select a model that cannot be
reliably restored in the target environment.

Also evaluate:

* CPU/GPU availability;
* PyTorch or TensorFlow backend compatibility;
* memory for profile/graph batch sizes;
* offline or network restrictions;
* checkpoint security and trusted loading;
* latency and throughput;
* monitoring and rollback requirements.

16. Use InversionConfig for 1-D candidates
------------------------------------------

:class:`pycsamt.ai.inversion.InversionConfig` captures supported 1-D settings:

.. code-block:: python

   from pycsamt.ai.inversion import InversionConfig

   config = InversionConfig(
       arch="resnet",
       n_layers=5,
       solver="mt1d",
       include_phase=True,
       log_thickness=True,
       augment_noise=0.02,
       epochs=150,
       batch_size=256,
       lr=1e-3,
       patience=20,
       val_frac=0.15,
       grad_clip=1.0,
       seed=42,
       checkpoint_dir="checkpoints",
       checkpoint_name="mt1d_resnet_5l",
   )
   config.validate()

   inverter = config.to_inverter()
   inverter.fit(dataset, **config.to_fit_kwargs())

``validate()`` checks supported architecture, solver, device, and numerical
ranges. It does not evaluate scientific suitability.

``weight_decay`` and ``min_delta`` are stored in the configuration for
documentation but are not currently included by ``to_fit_kwargs()`` because
``EMInverter1D.fit`` does not expose them. Do not assume every configuration
field affects training.

17. Record the selection decision
---------------------------------

A model-selection record should contain:

.. code-block:: text

   model_selection/
   ├── selection_plan.yml
   ├── dataset_reference.yml
   ├── split_indices.npz
   ├── candidates.csv
   ├── metrics_by_seed.csv
   ├── metrics_by_scenario.csv
   ├── compute_profile.csv
   ├── rejected_candidates.md
   ├── selected_model.yml
   └── review_signoff.md

For every candidate record:

* family, dimension, architecture, and parameter count;
* feature and target contract IDs;
* training configuration and seeds;
* validation/calibration/test metrics;
* response reconstruction results;
* domain-shift and robustness tests;
* runtime, memory, artifact size, and persistence method;
* known limitations and failure cases;
* selection status and reason.

Complete selection example
--------------------------

The following skeleton runs supported 1-D architectures repeatedly on one
dataset. Because the current ``fit(..., seed=...)`` uses the seed for its
internal train/validation split as well as training randomness, these runs
measure combined split-and-training sensitivity rather than a strictly fixed
split comparison. A controlled study should supply pre-separated arrays or a
project harness with frozen indices. Final selection must use project metrics,
not only the validation score stored by the inverter:

.. code-block:: python

   import csv
   from pathlib import Path

   from pycsamt.ai.inversion import InversionConfig

   output = Path("model_selection")
   output.mkdir(parents=True, exist_ok=True)

   rows = []
   for arch in ("fcn", "cnn1d", "resnet"):
       for seed in (11, 22, 33):
           cfg = InversionConfig(
               arch=arch,
               n_layers=5,
               solver="mt1d",
               include_phase=True,
               epochs=100,
               batch_size=256,
               lr=1e-3,
               patience=15,
               val_frac=0.15,
               seed=seed,
               checkpoint_dir=None,
               verbose=False,
           )
           cfg.validate()
           model = cfg.to_inverter()
           model.fit(dataset, **cfg.to_fit_kwargs())

           rows.append({
               "architecture": arch,
               "seed": seed,
               "best_validation_loss": model._meta.get(
                   "best_val_loss", float("nan")
               ),
           })

   with (output / "candidate_validation.csv").open(
       "w", newline="", encoding="utf-8"
   ) as stream:
       writer = csv.DictWriter(stream, fieldnames=rows[0])
       writer.writeheader()
       writer.writerows(rows)

This example reads a private metadata attribute for compact illustration.
Production selection should use stable public histories or independently
computed metrics and should never make a decision solely from private internals.

Selection checklist
-------------------

.. list-table::
   :header-rows: 1
   :widths: 31 69

   * - Check
     - Evidence
   * - Objective declared first
     - Primary metric, constraints, failure conditions, and decision tolerance.
   * - Dimension is geophysical
     - Survey layout, strike, tensor/tipper diagnostics, and target geometry.
   * - Parameterization is resolvable
     - Layer/depth tests, bandwidth, sensitivity, and boundary stability.
   * - Candidates are meaningful
     - Simple baseline, supervised/PINN/hybrid rationale, and architecture
       assumptions.
   * - Comparison is fair
     - Frozen data and splits, comparable tuning budget, repeated seeds, and
       untouched test set.
   * - Metrics are multidimensional
     - Model, response, structure, uncertainty, robustness, field, and compute.
   * - Domain shift is tested
     - Missingness, noise, nuisance effects, geology, geometry, and field
       coverage.
   * - Deployment is feasible
     - Persistence, backend, hardware, security, latency, and rollback.
   * - Decision is auditable
     - All candidates, rejected alternatives, declared selection rule, reviewer,
       and signoff.

Common mistakes
---------------

Avoid these errors:

* choosing dimension from architecture availability rather than geophysics;
* selecting the most complex model by default;
* comparing candidates trained on different random splits;
* tuning on the test or field-validation set;
* reporting only the best seed;
* increasing output depth cells and claiming improved physical resolution;
* using variable-layer NaN targets without a mask-aware strategy;
* confusing U-Net continuity with 2-D forward physics;
* using identity adjacency and describing the result as spatial GCN inversion;
* selecting graph radius without inspecting connectivity;
* fusing modalities without testing each modality alone;
* interpreting ensemble agreement as absence of shared bias;
* comparing only validation loss and ignoring response reconstruction;
* selecting a model that cannot be restored in production;
* assuming every ``InversionConfig`` field is used by ``to_fit_kwargs()``.

Next steps
----------

Continue with:

* :doc:`training` to fit the selected candidate reproducibly;
* :doc:`validation` to run the final acceptance protocol;
* :doc:`inference` to deploy the approved model contract;
* :doc:`uncertainty` to select and calibrate predictive intervals;
* :doc:`pinn_2d` for detailed physics-informed profile choices;
* :doc:`reporting` for the model-selection record and model card.
