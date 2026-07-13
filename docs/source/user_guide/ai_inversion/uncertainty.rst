.. _ai_inversion_uncertainty:

AI inversion uncertainty
========================

An AI inversion result is not complete when it contains only one resistivity
model.  Electromagnetic responses are non-unique, field data are noisy, the
training distribution is finite, and several processing and modelling choices
remain uncertain.  Uncertainty analysis describes how strongly the result is
supported, where alternatives remain plausible, and when the learned model is
being used outside its evidence base.

This page explains the uncertainty tools currently available in pyCSAMT and
the limits of their interpretation.  It assumes that training and model
selection have already followed :doc:`training` and :doc:`model_selection`.

What uncertainty should answer
------------------------------

A useful analysis should answer four separate questions:

#. How much do independently trained predictors disagree?
#. Do reported intervals contain held-out truth at approximately their stated
   rate?
#. How sensitive is the result to realistic perturbations of data and workflow
   choices?
#. Is the field observation sufficiently similar to the data on which the
   model and its calibration were established?

No single standard-deviation array answers all four.  Report the method and
the source represented by every uncertainty quantity.

Sources of uncertainty
----------------------

Aleatoric uncertainty
   Variability associated with observations: measurement noise, cultural
   interference, missing frequencies, imperfect station positioning, and
   unresolved small-scale structure.  It cannot generally be removed by
   training a larger network.

Epistemic uncertainty
   Uncertainty in the learned mapping caused by finite training coverage,
   architecture, optimization, and model parameters.  Deep-ensemble spread
   is primarily a practical proxy for this source.

Inverse-problem non-uniqueness
   Distinct conductivity structures may reproduce nearly the same EM
   response.  A narrow neural prediction does not prove that the physical
   inverse problem is unique.

Model-form uncertainty
   The synthetic earth family, dimensionality, forward solver, anisotropy
   assumptions, layer parameterization, and regularization may all be
   incomplete.

Workflow uncertainty
   Quality-control thresholds, static-shift correction, component selection,
   interpolation, frequency filtering, and coordinate handling can change the
   inferred model.

Domain-shift uncertainty
   Field inputs may differ from the synthetic training and calibration
   distributions.  Calibrated intervals do not automatically remain valid
   after this shift.

Use a source register in the final report.  For each source, state whether it
was quantified, tested qualitatively, assumed negligible, or left unresolved.

Deep ensembles
--------------

:class:`pycsamt.ai.inversion.EnsembleInverter` trains several copies of a
compatible base inverter.  Its mean is the point prediction and its
inter-member standard deviation measures disagreement among the members.

.. code-block:: python

   from pycsamt.ai.inversion import EMInverter1D, EnsembleInverter

   base = EMInverter1D(
       n_features=X_train.shape[1],
       n_layers=n_layers,
       solver="pytorch",
       device="auto",
   )
   ensemble = EnsembleInverter(
       base_estimator=base,
       n_estimators=5,
       seed=42,
   )
   ensemble.fit(
       X_train,
       y_train,
       epochs=250,
       batch_size=256,
       patience=25,
       val_frac=0.15,
       verbose=True,
   )

   mean, raw_std = ensemble.predict_with_uncertainty(X_test)
   quantiles = ensemble.predict_quantiles(
       X_test, q=(0.05, 0.25, 0.50, 0.75, 0.95)
   )

The returned arrays have shape ``(n_samples, n_parameters)``.  Before
calibration, ``raw_std`` is the sample standard deviation of member
predictions.  Quantiles are empirical quantiles across members; with only five
members, extreme quantiles are necessarily coarse.

.. caution:: Ensemble spread is not a confidence interval, a Bayesian
   posterior, or a complete measure of non-uniqueness.  Members share the same
   architecture, synthetic generator, target parameterization, and usually
   the same systematic errors.  They can agree closely and still all be wrong.

In the current implementation, member seeds affect both the internal data
split and stochastic optimization.  The spread therefore combines split and
optimization variability.  If those sources must be separated, use an
experimental training protocol with a fixed split and controlled
initializations.

Choosing the ensemble size
~~~~~~~~~~~~~~~~~~~~~~~~~~

More members provide a better estimate of model disagreement but increase
training and inference cost roughly in proportion to member count.  Evaluate
stability rather than selecting a count by convention:

#. train an initial ensemble;
#. recompute mean, standard deviation, and interval coverage using increasing
   numbers of members;
#. inspect whether rankings of uncertain stations and layers stabilize;
#. stop only when the remaining variation is acceptable for the decision.

Calibration data are a separate resource
----------------------------------------

Use four conceptual partitions when calibrated uncertainty is required:

``training``
   Fits the ensemble members.

``validation``
   Controls early stopping and hyperparameter selection.

``calibration``
   Learns how raw prediction errors relate to ensemble spread.

``test``
   Evaluates point accuracy and coverage after all choices are frozen.

The calibration set must not overlap training or validation, and the test set
must not be reused to recalibrate a disappointing result.  Split related noise
realizations, profiles, and surveys as groups so one parent earth model cannot
appear on both sides of a boundary.

Calibration also requires representative examples.  A large calibration set
from the wrong geology or acquisition design does not validate intervals for
the field survey.

Conformal prediction intervals
------------------------------

Calling :meth:`EnsembleInverter.calibrate
<pycsamt.ai.inversion.EnsembleInverter.calibrate>` attaches a split-conformal
predictor and a posterior calibrator:

.. code-block:: python

   ensemble.calibrate(X_cal, y_cal, alpha=0.10)

   center, lower, upper = ensemble.predict_intervals(
       X_test, alpha=0.10
   )

Here ``alpha=0.10`` requests nominal 90% coverage.  The conformal predictor
computes normalized residual scores on the held-out calibration set and uses
a finite-sample corrected quantile to scale ensemble standard deviations.
pyCSAMT takes the maximum normalized residual across output parameters for
each calibration sample.  Its diagnostic coverage therefore checks whether
all target parameters for a sample fall inside their bands simultaneously.

The statistical coverage statement depends on exchangeability: calibration
and future cases must behave like draws from the same distribution, and the
calibration cases must not have influenced training or model selection.  It is
a marginal repeated-sample statement, not a guarantee that every individual
station, layer, or geological subgroup has 90% coverage.  Distribution shift,
serial dependence, spatial grouping, or adaptive reuse of the calibration set
can invalidate it.

.. note:: An interval can achieve nominal coverage by becoming very wide.
   Always report interval width together with coverage and point accuracy.

Direct conformal use
~~~~~~~~~~~~~~~~~~~~

The lower-level :class:`pycsamt.ai.inversion.calibration.ConformalPredictor`
can wrap any fitted predictor exposing
``predict_with_uncertainty(X) -> (mean, std)``:

.. code-block:: python

   from pycsamt.ai.inversion.calibration import ConformalPredictor

   conformal = ConformalPredictor(ensemble, alpha=0.10)
   conformal.calibrate(X_cal, y_cal)
   center, lower, upper = conformal.predict_intervals(X_new)

Use either this explicit wrapper or ``ensemble.calibrate()`` in a workflow,
not two independent calibrations on the same data without a stated reason.

Coverage diagnostics
--------------------

Evaluate calibration on the untouched test set:

.. code-block:: python

   diagnostics = ensemble.coverage_diagnostics(
       X_test,
       y_test,
       alphas=(0.50, 0.20, 0.10, 0.05),
   )

   for alpha, actual in diagnostics.items():
       nominal = 1.0 - alpha
       print(f"nominal={nominal:.2f}, actual={actual:.3f}")

Plot nominal coverage against empirical coverage.  Values below the diagonal
indicate under-coverage; values far above it can indicate unnecessarily broad
intervals.  Add uncertainty bars or bootstrap intervals around empirical
coverage when the test set is small.

The convenience method ``ensemble.coverage(X, y, n_sigma=1.96)`` measures the
fraction of individual finite target entries within ``mean +/- 1.96 * std``.
This is different from conformal simultaneous sample coverage and should be
labelled explicitly.  The familiar 95% interpretation of ``1.96`` assumes an
appropriate Gaussian error model; check it empirically rather than assuming
it.

Go beyond one aggregate number.  Calculate coverage and interval width by:

* resistivity versus thickness parameters;
* layer or depth range;
* noise level and missing-data fraction;
* geological family and resistivity contrast;
* frequency coverage and station spacing;
* in-distribution versus stressed or shifted cases.

Small subgroups give noisy estimates, but they can reveal dangerous failures
hidden by acceptable global coverage.

Calibrated posterior samples
----------------------------

After ``ensemble.calibrate()``, posterior-like draws can be generated with
the fitted monotone recalibration map:

.. code-block:: python

   import numpy as np

   rng = np.random.default_rng(2026)
   draws = ensemble.predict_posterior(
       X_new,
       n_samples=1000,
       rng=rng,
   )
   # draws: (n_draws, n_cases, n_parameters)
   median = np.median(draws, axis=0)
   lo, hi = np.quantile(draws, (0.05, 0.95), axis=0)

The :class:`pycsamt.ai.inversion.calibration.PosteriorCalibrator` learns a
global monotone correction from standardized calibration residuals and a
global scale correction for raw standard deviations.  These samples are
useful for propagation and visualization, but they remain conditional on the
ensemble, calibration distribution, and chosen parameterization.  Do not call
them a complete geological posterior without qualifying those assumptions.

Respect physical parameter domains.  If sampling or symmetric intervals
produce non-positive resistivity or thickness, do not silently clip values
and continue as though the distribution were unchanged.  Prefer a
positive-domain parameterization such as logarithmic resistivity or thickness
where supported, or document the transformation and truncation explicitly.

From parameter uncertainty to model uncertainty
------------------------------------------------

For a layered 1-D prediction with ``L`` layers, the first ``L`` outputs are
resistivities and the remaining ``L - 1`` are interface thicknesses.  Convert
each draw to cumulative interface depths before making a depth plot.  Taking
quantiles of thickness first and then accumulating them is generally not the
same as accumulating every joint draw and then taking depth quantiles.

Preserve complete draws so correlations between parameters are retained.
Independent error bars lose trade-offs such as a conductive layer becoming
thicker while its resistivity increases.  These correlations are often central
to EM equivalence.

Forward-response uncertainty
----------------------------

An uncertainty band in parameter space should be tested in observation space:

#. select calibrated draws or representative ensemble members;
#. run each model through the same forward solver and frequency grid;
#. compare predicted apparent resistivity and phase with held-out observations;
#. summarize response quantiles and normalized residuals by frequency and
   component.

A wide range of models that all reproduce the response demonstrates
non-uniqueness.  Conversely, narrow parameter bands whose forward responses
miss the data reveal miscalibration, domain shift, or a forward-model mismatch.
Keep measurement error separate from predictive model spread when plotting the
comparison.

Sensitivity and perturbation tests
----------------------------------

Ensemble uncertainty should be supplemented with controlled perturbations.
Repeat inference after changing one scientifically plausible factor at a time:

* add noise consistent with estimated impedance uncertainty;
* mask selected frequencies or components;
* vary static-shift factors within defensible bounds;
* perturb station coordinates or topography within survey uncertainty;
* compare accepted interpolation and quality-control settings;
* change graph radius or adjacency construction for GCN inversion;
* compare plausible layer counts or target bounds;
* compare architectures or synthetic geological priors.

Store each scenario name, input version, seed, and resulting model.  The spread
across workflow scenarios represents a different source from inter-member
spread and should normally be reported separately.  Combining all sources by
root-sum-of-squares is justified only when independence and scale assumptions
are defensible.

Out-of-distribution checks
--------------------------

Calibration cannot rescue a model applied far outside its training support.
Before interpreting intervals, compare field inputs with training and
calibration inputs using quantities meaningful to the inversion:

* frequency range and missing-frequency pattern;
* apparent-resistivity and phase ranges by component;
* impedance uncertainty and signal-to-noise distribution;
* station count, spacing, and profile length;
* feature-space distance after the frozen training transform;
* QC flags, dimensionality indicators, and phase-tensor behavior.

No universal distance threshold is provided by pyCSAMT.  Establish alert
thresholds on held-out, deliberately shifted synthetic experiments.  When an
input is flagged, show the point prediction only as an exploratory result,
mark its calibrated interval as unsupported, and prefer additional modelling
or classical inversion.

Uncertainty for 2-D, graph, joint, and hybrid models
---------------------------------------------------

The same principles apply beyond 1-D, but the uncertainty object changes:

2-D inversion
   Evaluate interval width and coverage by depth and lateral position.  Use
   profile-level splits and inspect boundary artifacts introduced by resizing
   or padding.

Graph inversion
   Vary graph construction and coordinate uncertainty.  Report isolated nodes
   and distinguish uncertainty from weak connectivity from uncertainty in the
   learned weights.

Joint inversion
   Perturb or remove each modality.  Calibration is not transferable when a
   modality is absent, reordered, or drawn from a different noise regime.

PINN and hybrid inversion
   Explore data-error realizations, initial models, physics-loss weights,
   regularization strengths, and optimizer seeds.  Optimization variability is
   not a substitute for a posterior.  See :doc:`pinn_2d`.

Current high-level calibrated ensemble utilities are centered on compatible
ensemble predictors such as 1-D inversion.  Do not imply equivalent calibrated
interval support for every architecture unless the actual predictor exposes
the required uncertainty interface and has been tested.

Persistence and reproducibility
-------------------------------

Save ensemble members after training:

.. code-block:: python

   ensemble.save("checkpoints/mt1d_ensemble")

   from pycsamt.ai.inversion import EnsembleInverter
   restored = EnsembleInverter.load("checkpoints/mt1d_ensemble")

.. warning:: The current ensemble serialization preserves its members and
   ensemble metadata, but it does not serialize the attached conformal or
   posterior calibrators.  After loading, calibrated interval and posterior
   methods require recalibration from the original, versioned calibration set.
   Preserve calibration-set identity and calibration settings as first-class
   experiment artifacts.

For every uncertainty product, record:

* member count and member seeds;
* training, validation, calibration, and test partition identifiers;
* target transform and units;
* requested ``alpha`` or quantiles;
* calibration sample count and subgroup composition;
* software environment and model checkpoint identity;
* random generator seed used for posterior draws;
* domain-shift and perturbation tests performed.

Reporting uncertainty responsibly
----------------------------------

Report a central estimate, an interval, and its meaning together.  For
example: "median resistivity with a split-conformal 90% simultaneous parameter
band, calibrated on held-out synthetic models from the stated generator."
Avoid labels such as "90% confidence" unless the construction and repeated-
sampling interpretation genuinely support that term.

Every figure should state whether color or shading represents inter-member
standard deviation, empirical member quantiles, conformal limits, calibrated
posterior quantiles, or scenario sensitivity.  Include:

* coverage versus nominal level;
* interval width or sharpness;
* point-error metrics in physical units;
* forward-response misfit;
* subgroup and depth-dependent diagnostics;
* out-of-distribution flags;
* unresolved sources and known assumptions.

Common mistakes
---------------

* interpreting ensemble standard deviation as the total geological
  uncertainty;
* calibrating on training data or repeatedly tuning against the test set;
* claiming conformal guarantees under untested field-domain shift;
* reporting coverage without interval width;
* pooling all layers and geological regimes into one reassuring average;
* treating parameter-wise error bars as independent;
* clipping non-physical samples without documenting the changed distribution;
* assuming calibration objects survive ensemble serialization;
* using many posterior draws to hide a small or unrepresentative calibration
  set--Monte Carlo sample count does not create new information.

Decision checklist
------------------

Before using an AI inversion result for interpretation, verify that:

* uncertainty sources and intended decisions are explicitly named;
* calibration and test data are independent and representative;
* empirical coverage and interval width are both acceptable;
* diagnostics are resolved by parameter, depth, and relevant subgroup;
* parameter draws reproduce observations through forward modelling;
* perturbation and domain-shift tests have been completed;
* physical bounds and parameter correlations are preserved;
* the model, partitions, calibration recipe, and random seeds are reproducible;
* limitations are carried into :doc:`reporting` and the interpretation rather
  than reduced to a single confidence number.

Uncertainty is evidence about the reliability of an inference, not decoration
around a preferred model.  When the evidence is weak or the input is outside
the validated domain, the correct result is an explicit warning and a request
for additional data or modelling.
