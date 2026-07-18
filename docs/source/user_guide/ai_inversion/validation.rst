.. _ai_inversion_validation:

AI inversion validation
=======================

Validation determines whether an :term:`AI inversion` model is fit for a
stated use, not merely whether it can predict a synthetic test array.  A
defensible result must be accurate in parameter space, reproduce the
electromagnetic response, remain stable under realistic perturbations, behave
honestly when uncertain, and add value relative to a simpler
:term:`baseline model` or to an established classical inversion.

This page provides a validation protocol for models trained as described in
:doc:`training`.  Uncertainty calibration is treated in :doc:`uncertainty`,
while the final evidence package is assembled in :doc:`reporting`.

Validation is a claim with a scope
----------------------------------

Write the intended-use statement before computing metrics.  It is the domain
condition attached to every number that follows.  If a fitted model is
represented by :math:`g_\theta` and the measured feature vector by
:math:`\mathbf{x}`, validation is not the unconditional claim that
:math:`g_\theta(\mathbf{x})` is correct everywhere.  It is the narrower claim
that, for samples drawn from a declared domain :math:`\mathcal{D}`, the
prediction :math:`\hat{\mathbf{m}} = g_\theta(\mathbf{x})` is accurate enough
for a named decision and fails visibly outside that domain.  The statement
should name:

* survey method and components, such as :term:`AMT`, :term:`MT`,
  :term:`apparent resistivity`, and :term:`phase`;
* frequency range, station spacing, and expected data-quality range;
* geological families and resistivity contrasts represented;
* output :term:`dimensionality`, layer count, depth range, and parameter units;
* whether the result is for screening, initialization, interpretation, or a
  decision requiring quantitative accuracy;
* conditions under which the model must abstain or be reviewed.

A model validated for five-layer synthetic :term:`1D` earths is not thereby
validated for a strongly :term:`3D` field setting.  Every acceptance statement
should be read as "validated for this declared domain under these tests."

The evaluation partitions
-------------------------

Keep the roles of data partitions distinct:

``training``
   Updates model weights.

``validation``
   Controls early stopping, architecture choice, and hyperparameters.

``calibration``
   Fits conformal or posterior calibration when required.

``test``
   Measures final performance after all choices are frozen.

``challenge``
   Optional deliberately shifted cases used to identify failure boundaries.

The :term:`calibration set` and :term:`challenge set` have different jobs.  A
calibration set adjusts uncertainty statements after the model architecture is
chosen; a challenge set maps the failure boundary and should not quietly become
a source of tuning decisions.  The test set is not a second validation set.  If
its results cause another round of tuning, it has become development data and a
new untouched test set is required.

Split by independent geological realization or survey, not merely by row.
Noise variants from one :term:`forward model`, overlapping profile windows, or
nodes from one synthetic graph must remain in the same partition.  Record the
split manifest and verify that parent identifiers do not cross boundaries.  Any
path by which held-out information influences model selection, preprocessing,
or threshold setting is :term:`validation leakage`; it usually makes the final
score optimistic in a way that is difficult to repair after the fact.

.. important:: The current high-level supervised trainers estimate
   normalization statistics before making their internal training/validation
   split.  This limitation affects the internal validation curve.  It does not
   justify exposing the external test set: keep test and calibration arrays
   entirely outside ``fit`` and document the preprocessing behavior in the
   evidence report.

Freeze the artifact before testing
----------------------------------

Before opening the test set, preserve:

* model :term:`checkpoint` and checksum;
* configuration and random seeds;
* dataset and split identifiers;
* preprocessing and target transforms;
* pyCSAMT, backend, Python, and dependency versions;
* selected epoch and selection rationale;
* acceptance thresholds written without reference to test results.

Reload the checkpoint in a fresh process and reproduce a small reference
prediction.  This catches missing normalizers, incompatible shapes, and
serialization assumptions before scientific evaluation begins.

Parameter-space evaluation
--------------------------

For 1-D inversion, request predictions in the same representation as the
stored targets.  The default output of
:meth:`pycsamt.ai.inversion.EMInverter1D.predict` contains logarithmic
resistivity and, when ``log_thickness=True``, logarithmic thickness.

:term:`Model-space metric` values are easiest to audit when each target
quantity keeps its own units.  For sample :math:`i` and layer :math:`\ell`,
let :math:`y^{\rho}_{i\ell}=\log_{10}\rho_{i\ell}` and
:math:`\hat{y}^{\rho}_{i\ell}` be the predicted log-resistivity.  The
log-space error is
:math:`e^{\rho}_{i\ell}=\hat{y}^{\rho}_{i\ell}-y^{\rho}_{i\ell}`, while the
corresponding multiplicative physical error is
:math:`10^{|e^{\rho}_{i\ell}|}`.  If :math:`\mathcal{M}` is the finite-value
mask for the entries included in a metric, the masked RMSE is

.. math::

   \operatorname{RMSE}
   =
   \sqrt{
     \frac{1}{|\mathcal{M}|}
     \sum_{(i,j)\in\mathcal{M}}
     \left(\hat{y}_{ij}-y_{ij}\right)^2
   }.

A trained workflow would call ``inverter.predict(X_test, as_log_rho=True)``.
The miniature example below isolates the metric semantics with fixed arrays so
the reported output is reproducible:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.ai.training.metrics import layer_rmse, summarise
   >>>
   >>> n_layers = 3
   >>> y_test = np.array([
   ...     [2.0, 2.5, 3.0, 60.0, 120.0],
   ...     [2.3, 2.8, 3.2, 80.0, 160.0],
   ...     [1.9, 2.4, 2.9, np.nan, 140.0],
   ... ])
   >>> y_pred = np.array([
   ...     [2.1, 2.4, 3.1, 65.0, 110.0],
   ...     [2.2, 2.9, 3.1, 90.0, 150.0],
   ...     [2.0, 2.5, 3.0, 75.0, 135.0],
   ... ])
   >>>
   >>> overall = summarise(y_test, y_pred, n_layers=n_layers)
   >>> per_parameter = layer_rmse(y_test, y_pred)
   >>> {name: round(value, 4) for name, value in overall.items()}
   {'rmse': 5.4778, 'mae': 3.3933, 'r2': 0.9944, 'relative_rmse': 0.0692, 'depth_rmse': 0.0577}
   >>> print("resistivity RMSE by layer:", np.round(per_parameter[:n_layers], 4))
   resistivity RMSE by layer: [0.1 0.1 0.1]
   >>> print("thickness RMSE by interface:", np.round(per_parameter[n_layers:], 4))
   thickness RMSE by interface: [ 7.0711 10.    ]
   >>> print("finite target values:", int(np.isfinite(y_test).sum()))
   finite target values: 15

The available metric helpers ignore non-finite entries.  Report how many
values each metric used; a score based on a small surviving fraction can look
misleadingly good.

Interpret metrics in their actual space:

RMSE and MAE
   Express typical error magnitude.  In log10 resistivity, an absolute error
   of 0.3 corresponds to approximately a factor of two, not 0.3 ohm-m.

R-squared
   Compares residual variation with target variation.  It can be negative and
   can look strong when broad target ranges dominate scientifically important
   local errors.

Relative RMSE
   Divides errors by target magnitude.  It becomes unstable near zero and is
   awkward for log-valued targets, so do not use it without explaining the
   representation.

Depth-weighted RMSE
   Summarizes resistivity errors with layer-index weighting.  The current
   helper gives deeper layers less influence; it does not replace explicit
   depth-resolved results.

The first scalar RMSE in the example is intentionally poor as a headline
number: it mixes log-resistivity entries and metre-valued thickness entries.
Never combine resistivity and thickness into one headline score without also
showing their separate physical-unit errors.  Convert predictions back to
ohm-m and metres for decision-facing tables:

.. code-block:: pycon

   >>> rho_true = 10.0 ** y_test[:, :n_layers]
   >>> rho_pred = 10.0 ** y_pred[:, :n_layers]
   >>> h_true = y_test[:, n_layers:]
   >>> h_pred = y_pred[:, n_layers:]
   >>>
   >>> rho_factor_error = np.maximum(rho_pred / rho_true, rho_true / rho_pred)
   >>> thickness_abs_error = np.abs(h_pred - h_true)
   >>> print("median rho factor error:", np.round(np.nanmedian(rho_factor_error, axis=0), 3))
   median rho factor error: [1.259 1.259 1.259]
   >>> print("median thickness abs error m:", np.round(np.nanmedian(thickness_abs_error, axis=0), 3))
   median thickness abs error m: [ 5. 10.]
   >>> print("worst rho factor error:", round(float(np.nanmax(rho_factor_error)), 3))
   worst rho factor error: 1.259

Report medians, upper quantiles, and worst credible cases in addition to means.
Aggregate metrics should be broken down by layer, cumulative interface depth,
resistivity contrast, conductor thickness, and data quality.

Model geometry and boundary recovery
------------------------------------

Layered targets need geometry-aware diagnostics.  Compute interface depths by
cumulatively summing thicknesses for every sample,
:math:`z_{ik}=\sum_{\ell=1}^{k}h_{i\ell}`, where :math:`h_{i\ell}` is layer
thickness and :math:`z_{ik}` is the depth to interface :math:`k`.  Then
evaluate:

* interface-depth absolute and relative error;
* whether important conductive or resistive units are detected;
* top and base depth error for target units;
* resistivity-thickness trade-offs;
* false layers, merged layers, and missed thin layers;
* performance as interfaces approach the investigation limit.

For 2-D and 3-D outputs, report cell-wise error together with structural
measures: lateral boundary position, depth extent, connected-body recovery,
and smoothing or edge artifacts.  A low image RMSE can coexist with a badly
placed target boundary.

Response-space validation
-------------------------

The strongest physical check is to forward-model the predicted earth and
compare its response with the input observation.  If :math:`F` is the
:term:`forward operator`, :math:`\hat{\mathbf{m}}_i` is the predicted earth
model, :math:`\mathbf{d}_i` is the observed response vector, and
:math:`\sigma_{ij}` is the observational standard error for component
:math:`j`, the normalized residual is

.. math::

   r_{ij}
   =
   \frac{F_j(\hat{\mathbf{m}}_i)-d_{ij}}{\sigma_{ij}}.

The normalized :term:`response-space metric` is then

.. math::

   \operatorname{NRMS}
   =
   \sqrt{
     \frac{1}{|\mathcal{M}|}
     \sum_{(i,j)\in\mathcal{M}} r_{ij}^{2}
   }.

For a trained 1-D MT model, ``inverter.predict_models(X_test)`` can be passed
to :class:`pycsamt.forward.MT1DForward`.  The compact example below uses
already reconstructed responses so the residual definition is explicit:

.. code-block:: pycon

   >>> import numpy as np
   >>>
   >>> observed = np.array([
   ...     [2.0, 2.2, 2.4, 45.0, 47.0, 49.0],
   ...     [1.8, 2.1, 2.5, 43.0, 46.0, 50.0],
   ... ])
   >>> reconstructed = np.array([
   ...     [2.0, 2.1, 2.4, 44.0, 47.5, 49.5],
   ...     [1.9, 2.0, 2.6, 42.5, 46.0, 50.0],
   ... ])
   >>> sigma = np.array([
   ...     [0.1, 0.2, 0.2, 2.0, 2.0, 3.0],
   ...     [0.1, 0.2, 0.2, 1.0, 2.0, 3.0],
   ... ])
   >>>
   >>> residual = (reconstructed - observed) / sigma
   >>> nrms = np.sqrt(np.mean(residual ** 2))
   >>> inside_one_sigma = np.mean(np.abs(residual) <= 1.0)
   >>> print("normalized RMS:", round(float(nrms), 3))
   normalized RMS: 0.961
   >>> print("fraction within 1 sigma:", round(float(inside_one_sigma), 3))
   fraction within 1 sigma: 0.667
   >>> print("per-feature mean residual:", np.round(residual.mean(axis=0), 3))
   per-feature mean residual: [ 0.    -0.5    0.    -0.5    0.25   0.167]

In the full workflow, ``predict_models`` converts the network output to
:class:`pycsamt.forward.synthetic.LayeredModel` objects and enforces positive
minimum values.  Preserve the raw network output as well; otherwise automatic
flooring can hide invalid predictions.

Compare forward responses with the unnormalized test observations using:

* log-apparent-resistivity residual by frequency;
* phase residual in degrees, with the adopted phase convention stated;
* component-specific residuals where applicable;
* normalized residual using observational standard errors;
* fraction of residuals inside declared error bounds;
* systematic frequency trends and station-correlated misfit.

Do not reduce the comparison immediately to one scalar.  A similar total
misfit can conceal a narrow frequency band that controls the target depth.
When training and validation use one forward solver, repeat a subset with an
independent trusted implementation if possible; otherwise the neural model and
its validation can share the same simulator bias.

Baselines and ablations
-----------------------

:term:`AI inversion` validation requires comparisons that answer whether the
complexity adds value.  Use the same held-out cases, preprocessing, units, and
metrics for all methods.  Relevant :term:`baseline model` choices include:

* a constant predictor based only on training-target medians;
* a simple nearest-neighbour or low-capacity regression baseline;
* the selected classical Occam, ModEM, or MARE2DEM workflow where dimensionality
  and data support permit;
* an AI-initialized classical or hybrid refinement;
* alternative AI architectures selected before test evaluation.

Compare accuracy, forward misfit, uncertainty calibration, runtime, memory,
failure rate, and analyst effort.  A model that is faster but less accurate may
still be valuable as an initializer; state that narrower role rather than
claiming replacement of classical inversion.

:term:`Ablation study` results determine which inputs and mechanisms matter.
Repeat evaluation after removing :term:`phase`, individual components,
auxiliary modalities, graph edges, augmentation, or a physics-loss term.  An
unchanged score may reveal that a claimed information source is being ignored.

Robustness and stress testing
-----------------------------

Construct challenge sets before field deployment.  Vary one factor at a time
and then test realistic combinations:

* noise amplitude, outliers, and coherent cultural interference;
* missing frequencies, shortened bandwidth, and irregular sampling;
* static shift and residual distortion;
* station spacing, profile length, and coordinate error;
* resistivity and thickness near or beyond training bounds;
* thin conductors, sharp contrasts, anisotropy, and 2-D/3-D structures;
* component loss or modality dropout;
* graph radius, connectivity, and isolated stations;
* alternative forward-model or discretization settings.

Plot performance against stress severity.  Define the :term:`operating
envelope` at the point where error, coverage, or failure rate crosses a
predeclared limit.  Stress tests are not expected to all pass; their purpose is
to reveal when the model should abstain.

Uncertainty validation
----------------------

For an :term:`ensemble inversion`, point accuracy and uncertainty quality are
separate axes.  On an untouched test set, examine:

* empirical versus nominal conformal coverage;
* interval width or sharpness;
* coverage by layer, depth, geology, and noise regime;
* whether larger predicted spread corresponds to larger actual error;
* forward-response coverage after propagating parameter draws;
* behavior on deliberately shifted challenge cases.

For interval estimates :math:`[L_{ij}(\alpha), U_{ij}(\alpha)]` designed to
cover target entry :math:`y_{ij}` at nominal probability :math:`1-\alpha`, the
entry-wise empirical coverage is

.. math::

   \widehat{C}(\alpha)
   =
   \frac{1}{|\mathcal{M}|}
   \sum_{(i,j)\in\mathcal{M}}
   \mathbf{1}\{L_{ij}(\alpha)\le y_{ij}\le U_{ij}(\alpha)\}.

The mean interval width,
:math:`|\mathcal{M}|^{-1}\sum_{(i,j)\in\mathcal{M}}
\left(U_{ij}(\alpha)-L_{ij}(\alpha)\right)`, must be reported beside
coverage, because overly wide intervals can cover well while still being
scientifically unhelpful.

.. code-block:: pycon

   >>> import numpy as np
   >>>
   >>> y_test = np.array([
   ...     [2.0, 2.5, 3.0],
   ...     [2.1, 2.6, 2.9],
   ... ])
   >>> lower = np.array([
   ...     [1.85, 2.35, 2.85],
   ...     [1.95, 2.45, 2.75],
   ... ])
   >>> upper = np.array([
   ...     [2.15, 2.65, 3.15],
   ...     [2.25, 2.75, 3.05],
   ... ])
   >>> diagnostics = {
   ...     alpha: float(np.mean((lower <= y_test) & (y_test <= upper)))
   ...     for alpha in (0.20, 0.10, 0.05)
   ... }
   >>> mean_width = np.mean(upper - lower, axis=0)
   >>> diagnostics
   {0.2: 1.0, 0.1: 1.0, 0.05: 1.0}
   >>> print("mean 90% interval width:", np.round(mean_width, 3))
   mean 90% interval width: [0.3 0.3 0.3]

Coverage measured on individual entries with ``ensemble.coverage`` is not the
same diagnostic as the conformal simultaneous sample coverage.  Label the
definition used.  See :doc:`uncertainty` for the :term:`exchangeability`
assumptions and :term:`calibration set` requirements.

Field-data validation
---------------------

Synthetic truth enables parameter scoring; field data rarely provide complete
truth.  Field validation therefore combines several weaker but independent
lines of evidence:

Data compatibility
   Confirm component order, units, phase convention, frequency support,
   missing-value handling, station geometry, and preprocessing match the
   validated :term:`feature contract`.

Quality control
   Review impedance validity, uncertainty, coherence or quality indicators,
   outliers, static shift, phase tensor, dimensionality, and strike.  Load EDI
   data through :func:`pycsamt.emtools._core.ensure_sites`; do not bypass the
   canonical data model with ad-hoc arrays unless their provenance is retained.

Distribution support
   Compare field feature ranges and distances with training and calibration
   distributions.  Flag :term:`domain gap` and
   :term:`out-of-distribution diagnostic` results rather than interpreting a
   narrow ensemble spread as safety.

Forward consistency
   Forward-model predicted structures and compare with the measured response
   and observational errors.

Method agreement
   Compare with an appropriate classical inversion and with alternative
   parameterizations.  Investigate disagreements instead of averaging them
   away.

External evidence
   Compare interfaces and conductors with boreholes, logs, geology, hydrochemical
   information, seismic constraints, or known structures that were not used to
   tune the model.

Spatial coherence
   Adjacent stations should vary in geologically plausible ways, but do not use
   smoothness alone as proof: a biased model can be smoothly wrong.

Blind wells or withheld geotechnical observations provide especially valuable
checks.  Keep them hidden until the workflow and interpretation rules are
frozen.

Validation by inversion family
------------------------------

1-D models
   Validate resistivity, thickness, cumulative interface depth, response
   misfit, and lateral consistency across independently inverted stations.
   Diagnose where the 1-D assumption fails.

2-D U-Net models
   Split at profile or parent-model level.  Validate depth/lateral geometry,
   station-boundary effects, resizing behavior, and profile forward responses.

GCN models
   Split at survey level.  Verify node order and adjacency, then stratify error
   by node degree, edge distance, boundary location, and connected component.
   Compare against the identity-graph ablation.

Joint models
   Validate row alignment and each modality independently.  Test missing,
   corrupted, and shuffled modalities and compare against the strongest
   single-modality baseline.

PINN and hybrid models
   Validate total and component losses, initial-model sensitivity, physics
   residual, observed-data misfit, regularization sensitivity, and convergence
   from multiple seeds.  See :doc:`hybrid` and :doc:`pinn_2d`.

Failure analysis
----------------

Do not report only representative successes.  Create a failure table with one
row per failed or high-error case and include:

* case and parent-group identifier;
* input quality and out-of-distribution indicators;
* true and predicted parameter summaries where truth exists;
* response residuals by frequency and component;
* uncertainty interval and whether it covered the truth;
* architecture, seed, and preprocessing version;
* likely failure category and supporting evidence;
* disposition: accepted limitation, abstention rule, data correction, or model
  redevelopment.

Inspect the worst cases by a metric chosen before opening them.  Avoid deleting
difficult cases unless a reproducible data-quality rule, independent of model
error, requires exclusion.

Acceptance criteria
-------------------

Define thresholds from the intended scientific decision, observational error,
baseline performance, and acceptable risk.  Avoid universal thresholds such
as "R-squared above 0.9" without context.  A promotion gate can require all of
the following:

* parameter error below declared limits overall and in critical subgroups;
* forward-response residual consistent with the data-error model;
* improvement over the agreed baseline by a declared margin;
* calibrated coverage close to nominal without excessive interval width;
* no severe degradation inside the declared noise and acquisition envelope;
* an out-of-distribution or QC rule that catches unsupported inputs;
* acceptable failure, latency, and resource rates;
* reproducible checkpoint reload and prediction;
* documented review of worst cases and unresolved risks.

Use three outcomes rather than forcing pass or fail:

``accepted``
   All mandatory gates pass for the stated intended use.

``conditional``
   The model is allowed only inside a narrower operating envelope or with
   mandatory classical/analyst review.

``rejected``
   Evidence is insufficient or a critical gate fails.  Return to data design,
   model selection, or training rather than weakening the threshold after
   seeing results.

Statistical reporting
---------------------

Point estimates alone hide sampling variability.  Report sample counts and
confidence intervals for aggregate metrics, preferably using bootstrap units
that respect the grouping structure.  For example, resample whole geological
realizations or surveys rather than correlated rows.  Use paired resampling
when comparing two methods on the same cases.

Report results across training seeds, but distinguish variation across seeds
from uncertainty of the finite test sample.  Do not select the best seed on the
test set.  If many architectures, groups, or metrics are inspected, acknowledge
the resulting multiplicity and prioritize predeclared primary endpoints.

Minimum validation record
-------------------------

The validation artifact should contain:

* intended-use and exclusion statements;
* immutable model, configuration, environment, and data identifiers;
* partition manifest with grouping and leakage checks;
* parameter metrics in model and physical units;
* depth-, layer-, geology-, and quality-stratified results;
* forward-response residual figures and tables;
* baseline and ablation comparisons;
* robustness, domain-shift, and uncertainty diagnostics;
* field and independent-evidence comparisons;
* complete failure analysis;
* predeclared thresholds and final accepted, conditional, or rejected status.

Validation checklist
--------------------

Before signing off, confirm that:

* no test or calibration case influenced training or model selection;
* parent earths, profiles, and surveys do not cross partitions;
* metrics use stated units, transforms, masks, and denominators;
* physical geometry and forward responses were evaluated;
* comparisons use identical held-out cases and preprocessing;
* uncertainty coverage and width were both measured;
* challenge tests define a credible operating envelope;
* field inputs were checked for QC, dimensionality, and domain support;
* failures and negative evidence remain visible;
* the conclusion is limited to the tested intended use.

Validation is complete when another analyst can reproduce the evidence and
reach the same promotion decision without relying on undocumented judgment.
Carry that evidence into :doc:`reporting` and, for field interpretation, the
broader :doc:`../interpretation/index` workflow.
