.. _ai_inversion_inference:

AI inversion inference
======================

Inference applies a fitted AI inversion model to observations that were not
used to update its parameters. Although the network call may take only
milliseconds, a reliable inference workflow includes checkpoint verification,
exact preprocessing replay, input-domain assessment, output decoding, forward
response reconstruction, uncertainty, scientific review, and controlled
export.

This guide assumes that the dataset and model have already passed the
procedures in :doc:`data_preparation`, :doc:`training`, and :doc:`validation`.
Inference is not the stage at which an unvalidated checkpoint becomes valid.

.. admonition:: Compatibility before prediction
   :class: important

   Never run a checkpoint on an array merely because its shape matches. The
   feature meaning, order, units, grid, normalization, target parameterization,
   backend, and model metadata must match the approved training contract.

Inference workflow
------------------

#. identify the approved model and deployment policy;
#. verify checkpoint integrity and compatibility;
#. load and QC field observations;
#. reproduce the feature contract exactly;
#. apply saved missing-value and normalization transformations;
#. gate unsupported or out-of-distribution inputs;
#. run prediction without changing model state;
#. decode outputs into documented physical units;
#. reconstruct responses and inspect residuals;
#. calculate calibrated uncertainty where available;
#. apply acceptance, review, or rejection rules;
#. export predictions with station order and provenance.

1. Define the inference unit
----------------------------

The unit passed to the model depends on architecture:

1-D
   One row per station. Stations can be batched, but each prediction is
   independent unless an ensemble or later spatial operation is applied.

2-D
   One complete ordered profile panel per batch item. Station count, order,
   channel count, frequency grid, and depth target must match training.

Graph 3-D
   One survey graph or a batch of graphs sharing compatible node and feature
   conventions. Predictions depend on feature rows and adjacency together.

PINN or hybrid
   Observations are generally attached before optimization or refinement.
   Their ``predict()`` methods may not accept the same deployment array as a
   supervised surrogate; follow the class-specific fitted-state contract.

Do not split a profile or graph into convenient batches if that changes the
spatial context learned by the model.

2. Assemble the approved model package
--------------------------------------

An inference package should contain:

* checkpoint or fitted-model artifact;
* checksum and model identifier;
* architecture and backend information;
* feature schema and exact frequency/time grid;
* saved normalization or preprocessing state;
* output schema, transformations, and units;
* training-distribution summary;
* validation and calibration metrics;
* intended methods, geometry, and operating domain;
* known failure modes and rejection thresholds;
* software and dependency versions;
* approval status and reviewer.

A standalone weights file is insufficient. If the feature contract cannot be
reconstructed unambiguously, the checkpoint should not be deployed.

3. Verify integrity and compatibility
-------------------------------------

Before loading, verify the artifact checksum using the project's approved
integrity tooling. Then compare model metadata with the inference request:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Compatibility field
     - Required match
   * - Method and solver
     - MT, CSAMT, TEM, or other documented physics and response convention.
   * - Parameterization
     - Layer count, depth grid, thickness representation, and output order.
   * - Features
     - Components, apparent resistivity/phase or impedance channels, block
       order, transformations, and masks.
   * - Sampling grid
     - Exact frequencies or times, order, unit, and interpolation rule.
   * - Normalization
     - Saved training statistics and transformation version.
   * - Geometry
     - Station count/order for 2-D; coordinates and graph policy for 3-D.
   * - Backend
     - Compatible PyTorch/TensorFlow model construction and weight format.
   * - Operating domain
     - Supported ranges, noise, missingness, geology, and survey layout.

Reject incompatibility rather than modifying field arrays until the model runs.
Any deliberate adapter becomes a new preprocessing version that requires
validation.

4. Load a 1-D checkpoint
------------------------

:class:`pycsamt.ai.inversion.EMInverter1D` explicitly saves weights,
normalizers, hyperparameters, history, and selected metadata:

.. code-block:: python

   from pathlib import Path
   from pycsamt.ai.inversion import EMInverter1D

   checkpoint = Path("checkpoints/mt1d_resnet_5layer.pkl")
   inverter = EMInverter1D.load(checkpoint)

Loading restores the backend recorded in the checkpoint and rebuilds the
network. The compatible deep-learning backend must be installed. Treat a
checkpoint as executable model content and load only trusted artifacts.

Registry checkpoints can be requested with:

.. code-block:: python

   inverter = EMInverter1D.from_pretrained(
       "mt1d-resnet-5layer-v1",
       cache_dir="model_cache",
   )

Registry presence does not guarantee that weights are currently downloadable.
Preserve registry metadata and the downloaded file checksum. Do not replace a
failed download with new training while continuing to label the result as the
registry model.

5. Prepare 1-D field features
-----------------------------

Use the public bridge and the checkpoint's exact grid:

.. code-block:: python

   import numpy as np

   from pycsamt.ai.inversion import sites_to_features_1d
   from pycsamt.emtools._core import ensure_sites

   sites = ensure_sites(
       "data/AMT/WILLY_DATA/L18",
       recursive=True,
       verbose=0,
   )

   X_field, frequencies_hz, station_names = sites_to_features_1d(
       sites,
       comp="xy",
       n_freqs=32,
       freq_min=1e-3,
       freq_max=1e3,
   )

The bridge returns the block layout ``[log10(rho_a), phase_deg]``. Frequency
endpoints and ``n_freqs`` must come from model metadata, not memory or a nearby
example.

Check the matrix before prediction:

.. code-block:: python

   print("Feature shape:", X_field.shape)
   print("Stations:", station_names)
   print("Non-finite fraction:", np.mean(~np.isfinite(X_field), axis=1))

``sites_to_features_1d`` leaves values outside a station's observed range as
``nan``. Apply only the imputation/mask policy fitted and validated with the
checkpoint. The inverter's stored normalizer does not define a missing-value
policy by itself.

6. Replay preprocessing exactly
-------------------------------

For :class:`~pycsamt.ai.inversion.EMInverter1D`, ``predict()`` applies the
normalizer restored from the fitted object. Therefore supply features in the
same pre-normalized representation used by training. Do not normalize them a
second time.

For custom models, distinguish:

Raw transformation
   Conversion from physical observations into log resistivity, phase, masks,
   or other feature channels.

Grid transformation
   Sorting and interpolation onto the approved frequency/time/station axes.

Missing-value transformation
   Masking or imputation using the trained policy.

Statistical transformation
   Scaling with training-only means, standard deviations, or other fitted
   parameters.

Model input adaptation
   Batch and channel dimensions expected by the selected backend.

Version this pipeline as one unit. A change to any stage creates a new
deployment configuration and may invalidate prior validation.

7. Gate out-of-domain inputs
----------------------------

Run domain checks before predictions are visible to the interpreter. At
minimum, compare each field row against training feature percentiles:

.. code-block:: python

   train_low = np.load("model_package/training_p01.npy")
   train_high = np.load("model_package/training_p99.npy")

   outside = (X_field < train_low) | (X_field > train_high)
   outside_fraction = np.nanmean(outside, axis=1)

   review_mask = outside_fraction > 0.10
   for name, fraction, review in zip(
       station_names, outside_fraction, review_mask
   ):
       print(name, fraction, "REVIEW" if review else "within marginal gate")

This is only a marginal gate. Strongly correlated feature patterns can remain
out of domain even when every value lies inside its individual range. Where
available, also apply multivariate distance, latent-space, density, ensemble
disagreement, missingness, and geometry checks.

Define actions in advance:

``accept_for_prediction``
   Input falls inside validated operating conditions.

``predict_with_review``
   Limited departure exists, prediction is retained for expert evaluation, and
   the departure is visible in outputs.

``reject``
   Input is incompatible or materially outside the validated domain.

Do not use the visual plausibility of the prediction to override an input gate
without documenting a new review decision.

8. Run 1-D prediction
---------------------

Raw parameter vectors
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   y_pred = inverter.predict(
       X_field,
       as_log_rho=True,
   )
   print(y_pred.shape)

The output contains resistivity parameters first and thickness parameters
after them. With ``as_log_rho=True``, resistivity remains log10 ohm metres.
Thickness behavior depends on the fitted inverter's ``log_thickness`` setting.
Do not assume the entire output vector has one unit.

With ``as_log_rho=False``, the method converts resistivity to linear ohm metres
and converts thickness to linear metres when ``log_thickness`` is enabled.

LayeredModel output
~~~~~~~~~~~~~~~~~~~

Prefer :meth:`pycsamt.ai.inversion.EMInverter1D.predict_models` when downstream
code expects physical layered models:

.. code-block:: python

   models = inverter.predict_models(X_field)

   for name, model in zip(station_names, models):
       if model is None:
           print(name, "prediction could not be decoded")
           continue
       print(name, model.resistivity, model.thickness)

This method back-transforms resistivity and thickness, clips resistivity to a
positive minimum and thickness to at least one metre, and can return ``None``
when model construction fails. Clipping should be reported; a boundary value
may signal unsupported output rather than a physical estimate.

Single synthetic response
~~~~~~~~~~~~~~~~~~~~~~~~~

``predict_response(response)`` is a convenience for a compatible
``ForwardResponse``. Field :class:`~pycsamt.site.Sites` should use the public
bridge so station identity and interpolation remain explicit.

9. Run 2-D profile prediction
-----------------------------

Create the field panel with the contract used during training:

.. code-block:: python

   from pycsamt.ai.inversion import sites_to_panel_2d

   X_profile, frequencies_hz, station_names = sites_to_panel_2d(
       ordered_sites,
       n_freqs=32,
       n_components=4,
       comp_te="xy",
       comp_tm="yx",
       freq_min=1e-3,
       freq_max=1e3,
   )

   section = inverter_2d.predict(
       X_profile,
       as_log_rho=True,
   )

The input shape is
``(n_profiles, n_components, n_freqs, n_stations)``. The output shape is
``(n_profiles, n_depth, n_stations)``.

The 2-D inverter applies its learned input and target normalization internally.
The field station count must match the inverter's configured station axis.
Station ordering is preserved by the bridge and must already follow reviewed
profile chainage.

Use ``as_log_rho=False`` for linear ohm metres. Retain the fixed depth grid
from model metadata; the predicted array alone does not contain depth
coordinates.

.. warning::

   The current :class:`~pycsamt.ai.inversion.EMInverter2D` class exposes fit
   and predict behavior but does not provide the same explicit public
   ``save()/load()`` pair as ``EMInverter1D``. Preserve and restore 2-D fitted
   models only through a project-tested mechanism, and document backend and
   normalization state. Do not imply portable checkpoint support that the
   class does not expose.

10. Run graph 3-D prediction
----------------------------

Graph prediction requires features plus adjacency or coordinates:

.. code-block:: python

   from pycsamt.ai.inversion import (
       sites_to_coords_3d,
       sites_to_features_1d,
   )

   X_nodes, frequencies_hz, station_names = sites_to_features_1d(
       ordered_sites,
       comp="xy",
       n_freqs=32,
       freq_min=1e-3,
       freq_max=1e3,
   )
   coords_m = sites_to_coords_3d(ordered_sites)

   graph_prediction = inverter_3d.predict(
       X_nodes,
       coords=coords_m,
       radius=3000.0,
       as_log_rho=True,
   )

The field feature order, coordinate order, and station-name order must be
identical. Use authoritative projected coordinates where possible; the helper
can fall back to artificial uniform spacing when coordinates are unavailable.

If the fitted inverter stored an approved adjacency matrix, omit new geometry
only when field nodes and order are exactly the same. Otherwise pass a reviewed
adjacency explicitly:

.. code-block:: python

   graph_prediction = inverter_3d.predict(
       X_nodes,
       adjacency=approved_adjacency,
       as_log_rho=True,
   )

Inspect graph degree and disconnected nodes. Changing ``radius`` changes the
model context and is not a harmless inference option.

As with the 2-D class, the current graph inverter does not expose the explicit
public ``save()/load()`` pair provided by ``EMInverter1D``. Deployment requires
a tested project persistence mechanism that preserves weights, normalizers,
adjacency policy, backend, and configuration.

11. Predict graph uncertainty
-----------------------------

Use MC dropout where supported:

.. code-block:: python

   mean, standard_deviation = inverter_3d.predict_with_uncertainty(
       X_nodes,
       coords=coords_m,
       radius=3000.0,
       n_mc=50,
   )

Verify the exact signature in the installed API when passing adjacency versus
coordinates. The standard deviation is in output-parameter space and should be
reported with the same transformation as the mean. It captures stochastic
dropout variation, not total inversion or domain uncertainty.

12. Ensemble inference
----------------------

:class:`pycsamt.ai.inversion.EnsembleInverter` can return a mean, spread,
quantiles, intervals, or calibrated posterior draws:

.. code-block:: python

   from pycsamt.ai.inversion import EnsembleInverter

   ensemble = EnsembleInverter.load("checkpoints/mt1d_ensemble")

   mean = ensemble.predict(X_field)
   mean, std = ensemble.predict_with_uncertainty(X_field)
   quantiles = ensemble.predict_quantiles(
       X_field,
       q=(0.05, 0.50, 0.95),
   )

An ensemble has explicit ``save()`` and ``load()`` support. Its members must
share a compatible input and output contract. The current serialization stores
member checkpoints, estimator count, and seeds; it does **not** serialize the
attached conformal or posterior calibrators.

After loading, restore calibration through the project's reviewed calibration
data and call ``calibrate()`` again before requesting intervals. When the
in-memory ensemble has been calibrated on a separate calibration set:

.. code-block:: python

   center, lower, upper = ensemble.predict_intervals(
       X_field,
       alpha=0.10,
   )

   posterior_draws = ensemble.predict_posterior(
       X_field,
       n_samples=500,
       rng=np.random.default_rng(42),
   )

Conformal marginal coverage assumes calibration and deployment examples are
exchangeable. A synthetic calibration set does not guarantee 90% coverage for
shifted field data. Preserve the calibration dataset ID and empirical test
coverage with the inference output. Do not assume a loaded ensemble remains
calibrated merely because the ensemble directory was saved after calibration.

13. PINN and hybrid inference
-----------------------------

PINN and hybrid classes do not follow the same stateless ``predict(X)`` pattern
as all supervised surrogates. For example,
:class:`pycsamt.ai.inversion.PINNInverter1D.predict` returns fitted layered
models after observation-specific optimization, and
:class:`pycsamt.ai.inversion.HybridInverter1D.predict` returns the refined
models associated with its configured observations.

Treat these as inversion runs rather than pure checkpoint deployment:

* preserve observations and optimizer configuration;
* verify convergence and loss components;
* record regularization and initialization;
* inspect fitted responses and per-station residuals;
* do not reuse a result for different observations without rerunning the
  documented workflow.

See :doc:`pinn_2d` for the complete profile workflow.

14. Decode outputs safely
-------------------------

Every exported array needs a schema. Record:

* array shape and axis names;
* station order;
* depth/layer coordinates;
* parameter block boundaries;
* log10 versus linear values;
* resistivity and thickness units;
* masked, rejected, or clipped predictions;
* model/checkpoint identifier;
* preprocessing and inference version.

Validate decoded properties:

.. code-block:: python

   for name, model in zip(station_names, models):
       if model is None:
           continue
       if not np.all(np.isfinite(model.resistivity)):
           raise ValueError(f"{name}: non-finite resistivity")
       if not np.all(model.resistivity > 0):
           raise ValueError(f"{name}: non-positive resistivity")
       if not np.all(model.thickness > 0):
           raise ValueError(f"{name}: non-positive thickness")

Positive and finite values are necessary, not sufficient. Also compare with
training bounds and flag predictions near limits.

15. Reconstruct forward responses
---------------------------------

A predicted model should be passed through the relevant forward operator and
compared with field observations. For layered MT:

.. code-block:: python

   from pycsamt.forward import MT1DForward

   reconstructed = []
   for model in models:
       if model is None:
           reconstructed.append(None)
           continue
       reconstructed.append(
           MT1DForward(freqs=frequencies_hz).run(model)
       )

Calculate residuals in a clearly defined space. A robust report states:

* apparent resistivity, phase, complex impedance, or decay values used;
* linear or log residual transformation;
* observational errors and weights;
* interpolation grid;
* components included;
* missing-data mask;
* global and station/frequency summaries.

Do not accept a model solely because it resembles expected geology. Conversely,
response agreement alone cannot establish uniqueness.

16. Apply acceptance rules
--------------------------

Combine gates rather than relying on one RMS:

Input compatibility
   Exact feature contract, finite values, frequency support, geometry, and
   accepted missingness.

Domain support
   Training-envelope, multivariate, latent, or ensemble diagnostics.

Output validity
   Finite physical parameters, supported bounds, no unexplained clipping, and
   correct shape.

Response fit
   Error-aware reconstructed response with no structured residual pattern.

Uncertainty
   Calibrated interval or model spread suitable for the decision.

Scientific consistency
   Dimensionality, classical baseline, boreholes, geology, and neighboring
   observations.

Assign ``accepted``, ``needs_review``, or ``rejected`` per station/profile and
preserve the reason. Do not remove rejected stations from exports without a
rejection table.

17. Batch and resource behavior
-------------------------------

For 1-D models, batch size affects memory and throughput but should not change
station semantics. Confirm numerical consistency between a small batch and the
production batch.

For 2-D and graph models, one batch item contains a complete profile or survey.
Padding multiple geometries to one size requires a trained mask policy.

Run inference in evaluation mode. The provided supervised prediction methods
handle backend evaluation internally. MC dropout intentionally reintroduces
stochastic behavior for uncertainty; use a documented sample count and seed
where the API permits.

Record device, backend, precision, batch size, elapsed time, and memory-related
fallbacks. Numerical differences across devices should be evaluated against an
approved tolerance.

18. Export an inference record
------------------------------

A reproducible inference directory can contain:

.. code-block:: text

   inference/L18_mt1d_resnet_v001/
   ├── manifest.yml
   ├── checkpoint_reference.yml
   ├── input/
   │   ├── station_inventory.csv
   │   ├── field_features.npz
   │   └── domain_diagnostics.csv
   ├── predictions/
   │   ├── parameter_vectors.npz
   │   ├── layered_models.csv
   │   └── status_by_station.csv
   ├── responses/
   │   └── reconstructed_residuals.csv
   ├── uncertainty/
   │   └── prediction_intervals.npz
   ├── figures/
   └── review/
       └── inference_review.md

The manifest should record checkpoint checksum, feature contract, input survey
and QC IDs, preprocessing version, domain thresholds, software/backend, device,
random settings, output schema, accepted/rejected stations, reviewer, status,
and date.

Complete 1-D inference example
-----------------------------

.. code-block:: python

   from pathlib import Path
   import json
   import numpy as np

   from pycsamt.ai.inversion import (
       EMInverter1D,
       sites_to_features_1d,
   )
   from pycsamt.emtools._core import ensure_sites
   from pycsamt.forward import MT1DForward

   root = Path("inference/L18_mt1d_resnet_v001")
   root.mkdir(parents=True, exist_ok=True)

   inverter = EMInverter1D.load(
       "checkpoints/mt1d_resnet_5layer.pkl"
   )

   sites = ensure_sites(
       "data/AMT/WILLY_DATA/L18",
       recursive=True,
       verbose=0,
   )
   X, freqs, names = sites_to_features_1d(
       sites,
       comp="xy",
       n_freqs=32,
       freq_min=1e-3,
       freq_max=1e3,
   )

   # Replace this gate with the approved model-package policy.
   finite_fraction = np.mean(np.isfinite(X), axis=1)
   accepted = finite_fraction == 1.0

   models = [None] * len(names)
   if accepted.any():
       accepted_models = inverter.predict_models(X[accepted])
       for index, model in zip(np.flatnonzero(accepted), accepted_models):
           models[index] = model

   responses = []
   forward = MT1DForward(freqs=freqs)
   for model in models:
       responses.append(None if model is None else forward.run(model))

   rho = np.full((len(names), inverter.n_layers), np.nan)
   thickness = np.full((len(names), inverter.n_layers - 1), np.nan)
   for index, model in enumerate(models):
       if model is not None:
           rho[index] = model.resistivity
           thickness[index] = model.thickness

   np.savez_compressed(
       root / "predictions.npz",
       station_names=np.asarray(names),
       frequencies_hz=freqs,
       resistivity_ohm_m=rho,
       thickness_m=thickness,
       accepted=accepted,
   )

   manifest = {
       "checkpoint": "mt1d_resnet_5layer.pkl",
       "component": "xy",
       "n_freqs": 32,
       "feature_layout": "log10_rho_then_phase_deg",
       "accepted_stations": int(accepted.sum()),
       "rejected_stations": int((~accepted).sum()),
   }
   (root / "manifest.json").write_text(
       json.dumps(manifest, indent=2),
       encoding="utf-8",
   )

This example deliberately rejects rows with any missing feature. A production
model may use a different approved policy, but it must not improvise one during
deployment.

Review checklist
----------------

.. list-table::
   :header-rows: 1
   :widths: 31 69

   * - Check
     - Required evidence
   * - Model is approved
     - Identifier, checksum, model card, validation, calibration, and status.
   * - Contract matches
     - Method, features, grid, normalization, output schema, and geometry.
   * - Inputs are traceable
     - Survey, QC, processing, station order, coordinates, and exclusions.
   * - Missingness follows policy
     - Masks/imputation identical to validation and visible in results.
   * - Domain gate runs first
     - Thresholds, diagnostics, accepted/review/rejected status, and reasons.
   * - Prediction is immutable
     - Evaluation mode, no fitting, no normalization update, and recorded
       backend/device.
   * - Outputs decode correctly
     - Log/linear conversions, units, layers, depth, station axes, and clipping.
   * - Responses are reconstructed
     - Forward operator, components, errors, residual definition, and patterns.
   * - Uncertainty is conditional
     - Method, calibration set, coverage, domain-shift limitation, and sample
       count.
   * - Release is auditable
     - Manifest, arrays, status table, checkpoint reference, reviewer, and date.

Common mistakes
---------------

Avoid these errors:

* loading an untrusted checkpoint;
* checking only feature count rather than feature meaning;
* rebuilding normalization from field data;
* silently filling NaNs with convenient values;
* predicting before out-of-domain screening;
* treating all output columns as log10 values;
* losing station order when exporting predictions;
* changing profile station count or graph radius at deployment;
* calling graph-context output a full numerical 3-D inversion;
* treating MC dropout or ensemble spread as total uncertainty;
* claiming conformal field coverage from synthetic calibration alone;
* accepting geological-looking models without response reconstruction;
* exporting predictions without rejected-input records;
* implying portable 2-D/3-D checkpoint support that the public classes do not
  currently expose.

Next steps
----------

Continue with:

* :doc:`validation` for response, baseline, field, and acceptance testing;
* :doc:`uncertainty` for calibration, ensembles, and distribution shift;
* :doc:`reporting` for model cards and controlled prediction packages;
* :doc:`agents` for standard EDI-to-prediction orchestration;
* :doc:`pinn_2d` for observation-specific physics-informed inference.
