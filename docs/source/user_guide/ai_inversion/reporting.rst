.. _ai_inversion_reporting:

AI inversion reporting
======================

AI inversion reporting documents how a prediction was produced, where it is
valid, how it was tested, and what decisions it may support. A resistivity
section and a training-loss curve are not a complete report. The deliverable
must connect the field survey to a versioned dataset, feature contract, model
selection record, checkpoint, inference run, response reconstruction,
uncertainty analysis, and scientific review.

This page focuses on AI-specific reporting. For the later geological and
hydrogeophysical interpretation package, also use
:doc:`../interpretation/reporting`.

.. admonition:: Separate prediction from interpretation
   :class: important

   A network output is a model prediction. Geological labels, aquifer targets,
   and drilling recommendations are interpretations built from that prediction
   plus independent evidence. Report them in separate sections and preserve the
   transition from numerical result to professional judgment.

Reporting objectives
--------------------

An AI inversion package should allow a reviewer to answer:

* Which field observations and processing revision were used?
* Which feature transformation was applied?
* Which synthetic or field data trained the model?
* How were train, validation, calibration, and test examples separated?
* Why was this model family and architecture selected?
* Which exact checkpoint generated the prediction?
* Was each field input inside the validated operating domain?
* Does the predicted model reconstruct the observations?
* Which uncertainty sources were quantified and which were omitted?
* Which stations, profiles, or regions were accepted, flagged, or rejected?
* What classical and independent evidence supports the result?
* Who reviewed and approved the release?

Recommended reporting workflow
------------------------------

#. assign stable IDs to dataset, model, checkpoint, and inference run;
#. freeze the structured numerical outputs before writing narrative;
#. generate dataset and model cards;
#. report model-selection evidence and rejected candidates;
#. report training, calibration, synthetic test, and field validation results;
#. report input-domain status and per-station inference decisions;
#. report reconstructed response residuals and uncertainty;
#. state limitations, prohibited uses, and follow-up requirements;
#. assemble machine-readable artifacts, figures, and narrative;
#. perform scientific, ML, technical, and editorial review;
#. publish an immutable approved revision with a manifest and checksums.

1. Define report status and audience
------------------------------------

Use explicit status labels:

``experimental``
   Research output; architecture or dataset is still changing.

``working``
   Internal project result not ready for independent review.

``validation_candidate``
   Frozen model submitted to the declared validation protocol.

``review``
   Complete candidate package awaiting scientific or technical approval.

``approved``
   Versioned model and inference package accepted for its stated intended use.

``superseded``
   Retained for audit but no longer authoritative.

An approved model is approved only for the documented method, feature
contract, survey conditions, output, and decision. Approval does not transfer
automatically to another line, region, component, or acquisition system.

Tailor the summary to the audience, while keeping one authoritative technical
appendix. A manager may need risk and decision language; a geophysicist needs
residuals and dimensionality; an ML reviewer needs split and calibration
evidence; a downstream interpreter needs units, axes, masks, and uncertainty.

2. Use linked identifiers
-------------------------

Use immutable identifiers for:

``survey_id``
   Original acquisition and station inventory.

``processing_run_id``
   QC, masks, corrections, components, and frequency selection.

``dataset_id``
   Synthetic/field examples, generator revision, feature/target contracts, and
   split indices.

``model_selection_id``
   Candidate comparison and declared selection rule.

``training_run_id``
   Architecture, optimizer, seed, environment, history, and checkpoint output.

``checkpoint_id``
   Exact artifact and checksum.

``calibration_id``
   Conformal/posterior calibration dataset and state, where applicable.

``inference_run_id``
   Field features, domain gates, predictions, residuals, and statuses.

``report_package_id``
   This complete reviewed release.

A filename such as ``best_model.pkl`` is not a stable identifier. Include a
human-readable alias only in addition to a unique run/checkpoint ID.

3. Recommended package structure
--------------------------------

.. code-block:: text

   L18_ai_inversion_20260712_r01/
   ├── README.md
   ├── manifest.yml
   ├── CHANGELOG.md
   ├── cards/
   │   ├── dataset_card.md
   │   └── model_card.md
   ├── provenance/
   │   ├── environment.txt
   │   ├── source_revisions.yml
   │   └── checksums.sha256
   ├── configuration/
   │   ├── data_generation.yml
   │   ├── feature_contract.yml
   │   ├── model.yml
   │   ├── training.yml
   │   └── inference.yml
   ├── checkpoints/
   │   └── checkpoint_reference.yml
   ├── metrics/
   │   ├── model_selection.csv
   │   ├── training_history.csv
   │   ├── synthetic_test.csv
   │   ├── calibration.csv
   │   ├── field_residuals.csv
   │   └── robustness.csv
   ├── predictions/
   │   ├── station_status.csv
   │   ├── predicted_models.npz
   │   └── output_schema.yml
   ├── uncertainty/
   │   ├── intervals.npz
   │   └── uncertainty_summary.csv
   ├── figures/
   ├── report/
   │   ├── technical_report.md
   │   └── technical_report.html
   └── review/
       ├── comments.csv
       └── approval.yml

Store authoritative structured outputs separately from rendered reports. A PDF
or HTML file is convenient for reading but unsuitable as the only source of
prediction values.

4. Write a dataset card
-----------------------

The dataset card should describe both synthetic and field inputs.

Purpose and scope
   Intended inversion task, method, dimension, parameterization, target depth,
   and prohibited uses.

Synthetic model distribution
   Resistivity/thickness ranges, dependencies, geology scenarios, lateral
   correlation, rare targets, rejected samples, and sample count.

Forward physics
   Solver, dimension, frequencies/times, source geometry, numerical settings,
   approximations, and code version.

Observation effects
   Noise types and levels, missingness, outliers, static shift, source effects,
   coordinate perturbations, and augmentation probabilities.

Feature contract
   Names, order, units, transformations, grid, component convention, masks,
   imputation, and normalization.

Target contract
   Parameter order, shapes, units, log transformations, layer/depth axes, and
   padding/masks.

Splits
   Grouping policy, exact counts and indices, leakage controls, calibration
   role, and untouched test set.

Field data
   Survey/processing IDs, stations, components, frequency coverage, exclusions,
   corrections, and data-governance constraints.

Known gaps
   Physics, geology, noise, or geometry not represented.

The dataset card must travel with the checkpoint. A model cannot be reviewed
without knowing the distribution it learned.

5. Write a model card
---------------------

A model card should contain:

Model identity
   Name, version, checkpoint ID/checksum, code revision, backend, and status.

Intended use
   Supported method, dimension, survey conditions, output, users, and decision.

Out-of-scope use
   Unsupported components, grids, geology, geometry, depth, noise, or safety-
   critical decisions.

Architecture
   FCN/CNN/ResNet/U-Net/GCN/PINN/hybrid/joint/ensemble, parameterization,
   important hyperparameters, and architectural assumptions.

Training
   Dataset ID, optimizer, epochs, early stopping, regularization, seeds,
   hardware, precision, runtime, and selected checkpoint rule.

Selection
   Candidate set, primary metric, constraints, repeated-run statistics, and why
   rejected alternatives were not chosen.

Evaluation
   Synthetic model-space, response-space, structural, robustness, calibration,
   field, and operational metrics.

Uncertainty
   Method, calibration data, coverage, interval width, distribution-shift
   limitations, and omitted sources.

Ethical and operational considerations
   Data sensitivity, automation limits, human review, monitoring, and rollback.

Limitations
   Known failure cases, depth/resolution cautions, physics approximations, and
   domain boundaries.

6. Report model selection
-------------------------

Do not report only the winning architecture. Include a compact candidate table:

.. list-table::
   :header-rows: 1
   :widths: 16 16 18 18 16 16

   * - Candidate
     - Dimension
     - Synthetic test
     - Response fit
     - Calibration
     - Decision
   * - FCN-5L
     - 1-D
     - Reported metric
     - Reported metric
     - N/A or metric
     - Rejected/accepted reason
   * - ResNet-5L
     - 1-D
     - Reported metric
     - Reported metric
     - N/A or metric
     - Selected reason
   * - PINN-2D
     - pseudo-2-D
     - N/A or scenario
     - Reported metric
     - Scenario spread
     - Baseline/rejected reason

State whether repeated seeds changed internal train/validation splits. Report
means, spreads, failures, and tuning budget rather than the best run alone. See
:doc:`model_selection` for the complete protocol.

7. Report training behavior
---------------------------

For supervised models, retain per-epoch training and validation histories where
available. The inverter classes store training history internally, but there is
not one uniform public ``history_`` property across every class. Agent results
may expose ``train_history`` or ``loss_df``; PINN classes expose methods such as
``convergence_curve()``. Use the documented interface for the selected class
and store the underlying numeric table, not only the plotted curve.

Report:

* requested and completed epochs;
* batch size, optimizer, learning rate, and scheduling;
* early-stopping rule and selected epoch;
* training/validation losses and their definitions;
* augmentation and regularization;
* all seeds and number of repeated runs;
* backend, device, precision, and runtime;
* interrupted, failed, or restarted runs;
* checkpoint selection rule.

For PINN/hybrid runs, report observation-specific optimization, total loss,
physics and regularization terms where available, learning-rate/epoch
sensitivity, initialization, and convergence limitations. A PINN convergence
curve is not a supervised validation curve.

8. Report test metrics in context
---------------------------------

Model-space metrics
~~~~~~~~~~~~~~~~~~~

State whether values represent log10 resistivity, linear resistivity,
thickness, boundary depth, or a section. Include units, aggregation, and depth
or scenario stratification.

Response-space metrics
~~~~~~~~~~~~~~~~~~~~~~

State:

* forward operator;
* apparent resistivity, phase, impedance, or decay channels;
* log/linear residual space;
* observational errors and weights;
* frequency/time mask;
* station/component aggregation.

Do not label an unweighted log-apparent-resistivity RMS as a full impedance RMS.

Structural and decision metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Report conductor-top error, anomaly detection, section overlap, target ranking,
or other decision metrics only when their definitions and tolerances were
declared before model selection.

Performance slices
~~~~~~~~~~~~~~~~~~

Include results by:

* resistivity and thickness range;
* layer/depth;
* geology scenario;
* noise and missingness level;
* station density or graph degree;
* in-domain versus near-boundary inputs;
* known rare or difficult targets.

One aggregate test score can conceal operationally unacceptable failure.

9. Report field-domain status
-----------------------------

Every inference station or profile should have a status table:

.. code-block:: text

   station,feature_finite_pct,ood_score,domain_status,prediction_status,reason
   S00,100.0,0.12,inside,accepted,
   S01,93.8,0.41,near_boundary,needs_review,missing_high_frequency_band
   S02,71.9,0.88,outside,rejected,unsupported_frequency_coverage

Document the domain method and threshold. A marginal percentile envelope,
latent distance, density estimate, or ensemble disagreement score has different
meaning. Do not call all of them probability.

Keep rejected rows. Removing them makes the prediction set look more complete
and breaks traceability to the survey inventory.

10. Report predictions with an output schema
---------------------------------------------

For layered 1-D or graph output, include:

* station identifier and order;
* layer index;
* resistivity in linear ohm metres and/or explicitly named log10 form;
* thickness in metres;
* cumulative interface depth;
* accepted/review/rejected status;
* clipping or bound flags;
* checkpoint and inference IDs.

For 2-D output, include:

* profile and station axes;
* physical chainage;
* layer or depth axis;
* whether thickness varies by station;
* resistivity transformation and unit;
* missing/rejected columns;
* whether the section is supervised U-Net, pseudo-2-D PINN, or hybrid.

For graph output, additionally include coordinates, CRS, adjacency/radius,
degree, connected component, and whether a later volume was interpolated from
node predictions.

NPZ is useful for trusted arrays, while CSV is useful for audit. Accompany both
with a schema. Do not require a reviewer to infer array axes from code.

11. Report response reconstruction
----------------------------------

Preserve reconstructed responses and residuals per station, frequency/time,
and component. A useful table contains:

.. code-block:: text

   station,frequency_hz,component,rho_obs_ohm_m,rho_pred_ohm_m,
   phase_obs_deg,phase_pred_deg,valid,error_rho,error_phase

Summaries should include distributions, not only one global mean. Plot residual
heatmaps and curves with consistent scales. Flag structured residuals associated
with dimensionality, source effects, bad bands, station gaps, or out-of-domain
inputs.

For PINN 2-D, explicitly state that response physics are per-station MT1D with
lateral model smoothing. For supervised 2-D and graph workflows, state the
forward physics used to generate training examples separately from the forward
operator used for field reconstruction.

12. Report uncertainty and calibration
--------------------------------------

State the method:

Ensemble spread
   Member construction, number of members, seeds/architectures, raw standard
   deviation, and shared-bias limitation.

MC dropout
   Dropout rate, stochastic pass count, parameter space, and missing structural
   uncertainty.

Conformal intervals
   Calibration set ID/size, ``alpha``, empirical untouched-test coverage,
   interval width, and exchangeability limitation.

Posterior calibration
   Calibration inputs, base mean/spread, draw count, random generator, and
   interpretation of samples.

Scenario spread
   PINN/hybrid/classical configurations, data perturbations, accepted scenario
   criteria, and aggregation grid.

Report calibration and sharpness together. A very wide interval can achieve
coverage without supporting a useful decision. A narrow interval can be wrong
under field domain shift.

.. warning::

   ``EnsembleInverter.save()`` preserves members, estimator count, and seeds,
   but does not serialize attached conformal or posterior calibrators. A report
   must reference the calibration dataset and restoration procedure separately;
   loading the ensemble alone does not restore calibrated inference.

13. Report classical and independent validation
-----------------------------------------------

Include comparisons with an appropriate baseline:

* built-in or conventional 1-D inversion;
* Occam2D for defensible line geometry;
* ModEM or MARE2DEM where appropriate;
* boreholes, laboratory measurements, mapped geology, or other geophysics;
* withheld field responses or sites.

State which evidence was used in training, tuning, calibration, interpretation,
and final validation. Evidence cannot be called independent after it influenced
model or threshold selection.

Report disagreement rather than averaging it away. Explain whether it suggests
non-uniqueness, wrong dimension, domain shift, data error, or geological
ambiguity.

14. Report figures responsibly
-------------------------------

Minimum useful figures include:

* synthetic parameter and response distributions;
* field-versus-training coverage;
* training/validation or PINN convergence;
* synthetic test predictions versus truth;
* observed and reconstructed field responses;
* residuals by station/frequency/component;
* prediction section or layered profiles;
* uncertainty intervals or maps;
* classical/model comparison;
* station status and rejected inputs.

Every figure should identify model/checkpoint, dataset, run, axes, units,
transformations, missing data, status, and color scale. Use fixed scales across
model comparisons. A smooth color image must not hide station spacing or
rejected columns.

For varying-thickness layered sections, do not plot layer index as depth. Build
interface depths and resample to a physical grid while retaining original
layered values.

15. Use ReportAgent carefully
-----------------------------

:class:`pycsamt.agents.ReportAgent` can assemble selected agent results into a
survey-level Markdown report and optional HTML:

.. code-block:: python

   from pycsamt.agents import ReportAgent

   reporter = ReportAgent(
       report_title="L18 AI inversion review",
       formats=["md", "html"],
   )

   report_result = reporter.execute({
       "results": {
           "load": load_result,
           "qc": qc_result,
           "ai_inversion": ai_result,
       },
       "output_dir": "report/agent_draft",
   })

   print(report_result["report_path_md"])
   print(report_result.warnings)

The current implementation has predefined narrative sections for loading, QC,
static shift, phase analysis, forward modeling, and recommendations. It copies
figure paths from ``AgentResult`` objects, but it does not provide a dedicated
AI model-card section or serialize full prediction/uncertainty arrays.

Markdown is always attempted. HTML requires the optional ``markdown`` package.
Although the module documentation mentions PDF as an intended format, the
current ``execute()`` implementation does not generate or return a PDF path.
Do not promise PDF output from this agent without an additional tested step.

Use ``ReportAgent`` as a draft survey-report assembler, not as the complete AI
governance package described on this page.

16. Handle optional LLM narrative
---------------------------------

AI agents can return ``llm_interpretation`` when configured with a provider.
Report generated text only after human verification.

The narrative must not introduce:

* station counts or metrics absent from structured outputs;
* geological labels without independent evidence;
* claims of uniqueness or certainty;
* unsupported depth precision;
* hidden recommendations not approved by the reviewer;
* confidential data sent to an unauthorized external service.

Store generated text separately from approved conclusions, along with provider,
model, prompt context policy, timestamp, estimated cost, reviewer, and edits.
The structured numerical artifacts remain authoritative.

17. Write the narrative report
------------------------------

A detailed AI inversion report should contain:

Executive summary
   Intended decision, selected model, accepted coverage, principal prediction,
   uncertainty, limitations, and next action.

Survey and processing
   Acquisition, station inventory, QC, components, masks, corrections,
   dimensionality, and source run IDs.

AI inversion concept
   Supervised/PINN/hybrid approach, physical dimension, parameterization, and
   relationship to classical inversion.

Data preparation
   Synthetic generator, forward physics, noise, features, targets, splits,
   normalization, and field-domain comparison.

Model selection and training
   Candidates, declared metrics, selected architecture, hyperparameters,
   repeated runs, checkpoint, and compute environment.

Validation
   Synthetic test, response reconstruction, robustness, calibration, field
   evidence, failure cases, and acceptance thresholds.

Inference
   Field feature replay, domain gates, station statuses, predictions, clipping,
   residuals, and runtime.

Uncertainty
   Quantified sources, coverage and intervals, distribution shift, scenario
   spread, and omitted sources.

Interpretation
   Separate evidence-based geological meaning, alternatives, and confidence.

Limitations and intended use
   Unsupported conditions, prohibited uses, monitoring triggers, and review
   requirements.

Conclusions and recommendations
   Decision-focused conclusions conditional on the documented evidence.

18. Build a machine-readable manifest
-------------------------------------

.. code-block:: yaml

   schema_version: 1
   report_package_id: L18_ai_inversion_20260712_r01
   status: review
   survey_id: willy_L18_acquisition_v01
   processing_run_id: L18_processing_r03
   dataset_id: mt1d_5layer_v001
   model_selection_id: L18_model_selection_r01
   training_run_id: mt1d_resnet_seed42_r01
   checkpoint:
     id: mt1d_resnet_5l_ckpt_r01
     path: external_or_controlled_checkpoint_reference
     sha256: replace_with_verified_checksum
   feature_contract:
     component: xy
     n_frequencies: 32
     frequency_unit: Hz
     layout: log10_rho_block_then_phase_deg_block
   inference_run_id: L18_inference_r01
   prediction:
     resistivity_unit: ohm_m
     thickness_unit: m
     accepted_stations: 24
     needs_review_stations: 3
     rejected_stations: 2
   uncertainty:
     method: ensemble_conformal
     interval: P05_P95
     calibration_id: mt1d_calibration_v001
   software:
     package: pycsamt
     version: 2.0.0
   review:
     scientific_status: pending
     ml_status: pending

Use relative paths and checksums where files are inside the package. Do not
place API keys, credentials, personal information, or machine-specific absolute
paths in the manifest.

19. Validate the report package
-------------------------------

Structural validation
   Required cards, configuration, metrics, predictions, manifest, and review
   records exist and are non-empty.

Schema validation
   Array axes, CSV fields, units, transformations, null values, and IDs match
   the declared schema.

Numerical validation
   Prediction arrays match figures; log/linear conversion round-trips;
   thicknesses and resistivities are positive; intervals are ordered; station
   counts reconcile with the source inventory.

Provenance validation
   Checkpoint checksum, dataset/model IDs, code version, and configuration
   references resolve correctly.

Scientific validation
   Conclusions match structured results, response residuals, uncertainty,
   classical comparisons, and accepted/rejected statuses.

Reproduction validation
   A clean approved environment can load the artifact, replay preprocessing,
   reproduce a small reference prediction, and match within tolerance.

Security validation
   Checkpoints are trusted, external-service use was authorized, and the package
   contains no credentials or unintended sensitive information.

20. Review and approval roles
-----------------------------

Geophysical review
   Physics, dimensionality, response fit, classical comparison, depth, and
   geological plausibility.

ML review
   Dataset, leakage, architecture selection, training, calibration, robustness,
   domain shift, and checkpoint reproducibility.

Software review
   API contracts, units, shapes, persistence, environment, tests, and package
   integrity.

Interpretation review
   Independent evidence, alternatives, confidence, and decision language.

Governance approval
   Intended use, release status, privacy, revision, and downstream notification.

Record comments and dispositions. A material change to data, preprocessing,
checkpoint, calibration, domain gate, or conclusion requires a new revision and
revalidation of dependent artifacts.

21. Monitor deployed models
---------------------------

Reporting does not end at initial approval. Define monitoring for:

* feature-distribution and missingness drift;
* new acquisition or processing versions;
* station geometry outside the validated domain;
* response residual degradation;
* uncertainty or ensemble-disagreement increase;
* repeated review/rejection patterns;
* new boreholes or classical inversions that contradict predictions;
* backend or dependency changes affecting numerical output.

Specify retraining, recalibration, rollback, and reapproval triggers. Preserve
the prior approved package when issuing a replacement.

Complete reporting skeleton
---------------------------

The following code writes core structured inference products. It assumes the
arrays and statuses have already been validated:

.. code-block:: python

   from pathlib import Path
   import csv
   import json
   import numpy as np

   root = Path("L18_ai_inversion_20260712_r01")
   predictions_dir = root / "predictions"
   metrics_dir = root / "metrics"
   predictions_dir.mkdir(parents=True, exist_ok=True)
   metrics_dir.mkdir(parents=True, exist_ok=True)

   np.savez_compressed(
       predictions_dir / "predicted_models.npz",
       station_names=np.asarray(station_names),
       resistivity_ohm_m=resistivity_ohm_m,
       thickness_m=thickness_m,
       accepted=accepted_mask,
       needs_review=review_mask,
   )

   status_rows = []
   for index, name in enumerate(station_names):
       if accepted_mask[index]:
           status = "accepted"
       elif review_mask[index]:
           status = "needs_review"
       else:
           status = "rejected"
       status_rows.append({
           "station": name,
           "prediction_status": status,
           "reason": status_reasons[index],
       })

   with (predictions_dir / "station_status.csv").open(
       "w", newline="", encoding="utf-8"
   ) as stream:
       writer = csv.DictWriter(
           stream,
           fieldnames=["station", "prediction_status", "reason"],
       )
       writer.writeheader()
       writer.writerows(status_rows)

   manifest = {
       "report_package_id": "L18_ai_inversion_20260712_r01",
       "status": "review",
       "dataset_id": "mt1d_5layer_v001",
       "checkpoint_id": "mt1d_resnet_5l_ckpt_r01",
       "inference_run_id": "L18_inference_r01",
       "n_stations": len(station_names),
       "accepted": int(np.sum(accepted_mask)),
       "needs_review": int(np.sum(review_mask)),
       "rejected": int(np.sum(~accepted_mask & ~review_mask)),
       "resistivity_unit": "ohm_m",
       "thickness_unit": "m",
   }
   (root / "manifest.json").write_text(
       json.dumps(manifest, indent=2),
       encoding="utf-8",
   )

This skeleton does not replace dataset/model cards, checksums, response
residuals, uncertainty, figures, narrative, or approval records.

Release checklist
-----------------

.. list-table::
   :header-rows: 1
   :widths: 31 69

   * - Check
     - Acceptance evidence
   * - Scope and status are explicit
     - Intended use, prohibited use, audience, revision, and approval state.
   * - IDs connect the lineage
     - Survey, processing, dataset, selection, training, checkpoint,
       calibration, inference, and report IDs.
   * - Dataset card is complete
     - Priors, physics, noise, features, targets, splits, field data, and gaps.
   * - Model card is complete
     - Architecture, training, selection, metrics, limitations, and deployment.
   * - Checkpoint is identifiable
     - Trusted path/reference, checksum, backend, compatibility, and status.
   * - Metrics are contextual
     - Definitions, units, aggregation, scenarios, repeated runs, and failures.
   * - Field domain is visible
     - Per-station/profile score, threshold, accepted/review/rejected status.
   * - Responses are reconstructed
     - Forward operator, component, weighting, residual tables, and patterns.
   * - Uncertainty is calibrated conditionally
     - Method, calibration ID, coverage, sharpness, shift, and omitted sources.
   * - Interpretation is separate
     - Evidence, alternatives, confidence, and responsible reviewer.
   * - Package is reproducible
     - Configurations, environment, arrays, schemas, manifest, and reference
       prediction.
   * - Release is controlled
     - Checksums, comments, approvals, changelog, immutable archive, and
       monitoring triggers.

Common reporting mistakes
-------------------------

Avoid these errors:

* reporting only a prediction figure and global RMS;
* omitting the dataset or feature contract from the model record;
* naming one run ``best`` without reporting the selection process;
* presenting synthetic test accuracy as field validation;
* hiding rejected stations or model failures;
* confusing log10 and linear resistivity or thickness;
* losing station order or graph geometry in exports;
* calling U-Net/GCN/PINN output full 2-D/3-D physics without qualification;
* reporting ensemble or dropout spread as total uncertainty;
* claiming conformal field coverage from shifted synthetic calibration;
* assuming a saved ensemble retains uncertainty calibrators;
* treating LLM text as an authoritative scientific result;
* promising PDF output from ``ReportAgent`` when its current execution path
  produces Markdown and optional HTML only;
* changing preprocessing or thresholds without revising the model package;
* approving a checkpoint without a trusted checksum and reproducibility test.

Next steps
----------

Use this page with:

* :doc:`data_preparation` for the dataset card source;
* :doc:`model_selection` for candidate and decision evidence;
* :doc:`training` for training-run provenance;
* :doc:`validation` for acceptance metrics and field evidence;
* :doc:`inference` for station-level status and structured predictions;
* :doc:`uncertainty` for calibration and distribution shift;
* :doc:`agents` for ``AgentResult`` and report orchestration;
* :doc:`../interpretation/reporting` for downstream interpretation delivery.
