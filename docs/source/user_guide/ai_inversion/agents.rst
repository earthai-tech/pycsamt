.. _ai_inversion_agents:

AI inversion agents
===================

The agents in :mod:`pycsamt.agents` provide task-oriented orchestration around
the AI inversion classes in :mod:`pycsamt.ai.inversion`. They can load EDI
data, construct model inputs, train or load an inverter, predict resistivity,
compute selected diagnostics, create figures, save outputs, and return a
standard :class:`pycsamt.agents.AgentResult`.

Agents do not replace the science API. Use the lower-level inverter classes
when you need complete control of datasets, network construction, training
loops, loss functions, checkpoints, or research experiments. Use agents when
the built-in workflow matches the task and a consistent result contract is
valuable.

.. admonition:: Human review remains required
   :class: important

   An agent can coordinate calculations but cannot establish that the
   training distribution represents the field site, that a predicted model is
   unique, or that geological interpretation is justified. Review input data,
   training provenance, forward response fit, warnings, uncertainty, and
   independent evidence before using a result operationally.

Agent map
---------

.. list-table::
   :header-rows: 1
   :widths: 23 25 52

   * - Agent
     - Main engine
     - Intended use
   * - :class:`pycsamt.agents.AIInversionAgent`
     - :class:`pycsamt.ai.inversion.EMInverter1D`
     - Independent layered 1-D prediction at each station, trained from
       synthetic MT responses or loaded from a checkpoint.
   * - :class:`pycsamt.agents.Inv2DAgent`
     - :class:`pycsamt.ai.inversion.EMInverter2D`
     - U-Net profile inversion using the complete station–frequency panel to
       predict a laterally coherent section.
   * - :class:`pycsamt.agents.Inv3DAgent`
     - :class:`pycsamt.ai.inversion.GCNInverter3D`
     - Graph-based multi-station prediction using coordinates and spatial
       adjacency, with optional MC-dropout spread.
   * - :class:`pycsamt.agents.EnsembleAgent`
     - :class:`pycsamt.ai.inversion.EnsembleInverter`
     - Multiple 1-D estimators, prediction intervals, and optional conformal
       calibration.
   * - :class:`pycsamt.agents.PINNInversionAgent`
     - PINN 1-D, 2-D, or 3-D inverter
     - Physics-informed optimization without labelled training targets.
   * - :class:`pycsamt.agents.ModelZooAgent`
     - pyCSAMT model registry and ``EMInverter1D``
     - List model metadata, obtain a released checkpoint, or run the zoo
       prediction shortcut.

Additional agents such as ``HybridInversionAgent``, ``JointInversionAgent``,
``AnomalyDetectionAgent``, ``SensitivityAgent``,
``InversionComparisonAgent``, and ``InversionEvaluationAgent`` support more
specialized review and combination workflows. Start with one primary agent and
add these only when the scientific question requires them.

Common execution contract
-------------------------

All primary agents use:

.. code-block:: python

   result = agent.execute({
       "path": "data/AMT/WILLY_DATA/L18",
       "output_dir": "outputs/ai_inversion/L18",
   })

``path`` or ``sites``
   Required observed input. ``path`` is passed through the canonical
   :func:`pycsamt.emtools._core.ensure_sites` loader. ``sites`` accepts an
   already loaded :class:`pycsamt.site.Sites` object.

``output_dir``
   Optional output directory for supported checkpoints and figures. If it is
   omitted, in-memory results can still be returned.

Some agents accept additional execution overrides, such as ``period_range``,
``coords``, ``adjacency``, ``epochs``, or ``n_train_samples``. Constructor
parameters remain the clearest way to define a reproducible workflow; runtime
overrides should be recorded with the result.

AgentResult
-----------

Every execution returns :class:`pycsamt.agents.AgentResult`, not a bare array:

.. code-block:: python

   if result.status == "failed":
       raise RuntimeError(
           f"{result.error}\nSuggested fix: {result.error_fix_hint}"
       )

   print(result.status)
   print(result.summary)
   print(result.warnings)
   print(result.elapsed_seconds)

   rms = result.get("rms_global")

The standard fields are:

``status``
   ``"success"``, ``"failed"``, or ``"needs_review"``. A truth test is false
   only for ``"failed"``; therefore check the exact status when
   ``"needs_review"`` must stop a production workflow.

``summary``
   Short human-readable execution summary.

``data``
   Agent-specific arrays, inverter objects, figures, and paths. Dictionary-like
   access is available through ``result["key"]`` and ``result.get(...)``.

``warnings``
   Non-fatal issues. A successful status does not make these optional reading.

``error`` and ``error_fix_hint``
   Failure detail and suggested remediation.

``llm_interpretation``
   Optional generated narrative when an API key is configured. This text is
   commentary, not a validated scientific conclusion.

``elapsed_seconds`` and ``cost_estimate_usd``
   Execution duration and estimated LLM cost. Neural-network compute cost is
   not represented by the LLM cost field.

Prepare data before running an agent
------------------------------------

Load and inspect the same data independently before delegating inversion:

.. code-block:: python

   from pycsamt.emtools._core import ensure_sites

   sites = ensure_sites(
       "data/AMT/WILLY_DATA/L18",
       recursive=True,
       verbose=0,
   )
   print("Usable stations:", len(sites))

Review QC, impedance components, frequencies, coordinates, dimensionality,
static shift, and processing provenance first. The 1-D agent includes an
important fast-fail guard: it verifies that at least one station can produce a
finite impedance feature vector before spending time generating synthetic data
and training. This guard catches empty or corrupted inputs, but it is not a
complete QC assessment.

Deep-learning backend
---------------------

The inversion agents lazily require PyTorch or TensorFlow. If neither backend
is available, they return a failed ``AgentResult`` with an installation hint
instead of failing package import.

Verify the project environment before a long run. Also record backend name,
version, hardware, random seeds where exposed, and numerical precision. Results
can differ across backends and devices even when high-level parameters match.

AIInversionAgent: 1-D workflow
------------------------------

:class:`pycsamt.agents.AIInversionAgent` performs this built-in sequence:

#. load observed stations;
#. interpolate impedance-derived features onto a common frequency grid;
#. generate synthetic layered MT training examples unless a checkpoint loads;
#. fit :class:`~pycsamt.ai.inversion.EMInverter1D`;
#. predict layer resistivities and thicknesses for each usable station;
#. forward-model each prediction and calculate a log-resistivity RMS where
   possible;
#. create convergence and section figures;
#. optionally save the inverter and figures.

Train from synthetic data
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pycsamt.agents import AIInversionAgent

   agent = AIInversionAgent(
       arch="resnet",
       n_layers=5,
       n_train_samples=10_000,
       epochs=100,
   )

   result = agent.execute({
       "sites": sites,
       "output_dir": "outputs/ai_inversion/L18_1d",
   })

Supported architecture names are ``"resnet"``, ``"cnn1d"``, and ``"fcn"``.
The default workflow uses 40 log-spaced frequencies from :math:`10^{-4}` to
:math:`10^3` Hz, 2,000 synthetic examples, five layers, and 30 epochs. Those
defaults are convenient workflow settings, not evidence of production
adequacy.

The generated training set currently uses the MT 1-D solver, 3% noise, seed
42, and one generation job. If a project requires different priors, correlated
noise, missing-data patterns, other solvers, or explicit train/validation/test
control, use :class:`pycsamt.ai.inversion.EMInverter1D` directly as documented
in :doc:`data_preparation` and :doc:`training`.

Run from a local checkpoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   agent = AIInversionAgent(
       pretrained="checkpoints/L18_resnet.pkl",
   )
   result = agent.execute({"sites": sites})

If checkpoint loading fails, the agent adds a warning and falls back to fresh
training. Inspect warnings so an expensive fallback is not mistaken for use of
the approved checkpoint.

Inspect 1-D outputs
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   predictions = result["predictions"]
   rms_per_station = result["rms_per_station"]
   rms_global = result["rms_global"]
   history = result["train_history"]
   first_model = result["best_model"]

``predictions`` maps station names to log10 resistivity arrays. The
``best_model`` label means the first successfully predicted station in the
current implementation; it is not selected as the statistically best station.
It contains linear resistivity and thickness arrays plus log resistivity.

The forward RMS is computed in log10 apparent-resistivity space after
interpolation. It is useful as a response-fit diagnostic, but it is not a
complete impedance likelihood, phase fit, uncertainty-weighted RMS, or proof of
geological correctness.

Inv2DAgent: profile workflow
----------------------------

:class:`pycsamt.agents.Inv2DAgent` assembles a station–frequency input panel,
generates synthetic profile examples, trains a U-Net-style inverter, and
predicts the complete profile at once:

.. code-block:: python

   from pycsamt.agents import Inv2DAgent

   agent = Inv2DAgent(
       n_depth=40,
       n_freqs=32,
       n_components=4,
       arch="unet",
       n_train_profiles=500,
       n_stations_per_profile=20,
       epochs=80,
   )

   result = agent.execute({
       "sites": sites,
       "period_range": [1e-3, 1e3],
       "output_dir": "outputs/ai_inversion/L18_2d",
   })

Important outputs are ``pred_section`` with shape
``(n_depth, n_stations)``, ``depths_km``, ``station_names``, ``rms_global``,
the fitted ``inverter``, and figure dictionaries.

The learned lateral continuity comes from the training construction and network
architecture; it is not equivalent to a conventional 2-D EM forward operator
or proof that the field structure is two-dimensional. Validate the section
against dimensionality evidence and classical 2-D inversion where feasible.

Inv3DAgent: spatial graph workflow
----------------------------------

:class:`pycsamt.agents.Inv3DAgent` represents stations as graph nodes. Edges
connect stations within ``radius`` unless a normalized adjacency matrix is
supplied:

.. code-block:: python

   import numpy as np
   from pycsamt.agents import Inv3DAgent

   coords_m = np.array([
       [0.0, 0.0],
       [500.0, 0.0],
       [0.0, 500.0],
       [500.0, 500.0],
   ])

   agent = Inv3DAgent(
       n_layers=5,
       n_freqs=32,
       n_train_profiles=300,
       epochs=80,
       radius=1000.0,
       hidden=(256, 128, 64),
       dropout=0.1,
       n_mc=50,
   )

   result = agent.execute({
       "sites": sites,
       "coords": coords_m,
       "output_dir": "outputs/ai_inversion/survey_3d",
   })

Pass projected coordinates in metres. When coordinates are omitted, the agent
attempts to derive planar positions from EDI latitude/longitude. For controlled
production work, supply reviewed coordinates explicitly and ensure their order
matches the loaded stations.

Outputs include ``pred_rho`` in log10 resistivity, ``pred_thick`` in log10
metres, optional ``pred_uncertainty``, ``depths_km``, coordinates, adjacency,
edge count, global RMS, figures, and the inverter.

``pred_uncertainty`` is MC-dropout standard deviation across stochastic
passes. It represents one model-based uncertainty component; it does not
include training-prior misspecification, coordinate error, inversion
non-uniqueness, or field domain shift.

Check graph connectivity explicitly:

.. code-block:: python

   adjacency = result["adjacency"]
   degree = (adjacency > 0).sum(axis=1) - 1
   print("Graph degree by station:", degree)

Disconnected or weakly connected stations cannot receive the intended spatial
context. A very large radius can oversmooth across unrelated structures.

EnsembleAgent: uncertainty-aware 1-D workflow
---------------------------------------------

:class:`pycsamt.agents.EnsembleAgent` trains independent 1-D estimators with
different seeds and returns prediction intervals:

.. code-block:: python

   from pycsamt.agents import EnsembleAgent

   agent = EnsembleAgent(
       n_estimators=5,
       arch="resnet",
       n_layers=5,
       n_train_samples=10_000,
       epochs=100,
       calibrate=True,
   )
   result = agent.execute({
       "sites": sites,
       "output_dir": "outputs/ai_inversion/L18_ensemble",
   })

Inspect ``pred_mean``, ``pred_std``, ``pred_lo``, ``pred_hi``, ``coverage``,
and ``rms_global``. Bounds represent the agent's ensemble and optional
calibration procedure, conditional on its synthetic distribution. Empirical
coverage on held-out synthetic examples is not automatically field coverage.

PINNInversionAgent
------------------

:class:`pycsamt.agents.PINNInversionAgent` supports dimensions 1, 2, and 3 and
optimizes a physics-informed loss without labelled model targets:

.. code-block:: python

   from pycsamt.agents import PINNInversionAgent

   agent = PINNInversionAgent(
       dim=2,
       n_layers=10,
       depth_max=2000.0,
       smoothness_weight=0.01,
       lateral_weight=0.005,
       epochs=300,
       lr=1e-2,
       solver="mt1d",
   )
   result = agent.execute({
       "sites": sites,
       "output_dir": "outputs/ai_inversion/L18_pinn2d",
   })

Outputs include ``section``, optional layered ``models`` for 1-D,
``rms_per_station``, ``rms_global``, loss and residual dataframes when
available, figures, and the fitted inverter.

Physics-informed does not mean assumption-free. Solver dimension,
regularization weights, depth parameterization, optimizer, learning rate, and
stopping behavior all condition the result. See :doc:`pinn_2d` for the full
scientific workflow.

ModelZooAgent
-------------

List registry metadata without downloading weights:

.. code-block:: python

   from pycsamt.agents import ModelZooAgent

   zoo = ModelZooAgent()
   listed = zoo.execute({"action": "list"})
   for item in listed["details"]:
       print(item["name"], item["arch"], item["n_layers"])

Download and prediction actions require ``model_name``:

.. code-block:: python

   downloaded = zoo.execute({
       "action": "download",
       "model_name": "mt1d-resnet-5layer-v1",
   })

   predicted = zoo.execute({
       "action": "predict",
       "model_name": "mt1d-resnet-5layer-v1",
       "sites": sites,
       "output_dir": "outputs/ai_inversion/zoo_prediction",
   })

Check exact status and ``checkpoint_path``. When weights are not released or
available, download can return ``"needs_review"``. The prediction action may
fall back to on-the-fly training; confirm warnings and provenance before
describing a result as pretrained.

Checkpoint and output policy
----------------------------

An approved AI run should preserve:

* input survey and QC identifiers;
* exact agent class and constructor parameters;
* runtime payload and overrides;
* backend, dependency versions, hardware, and seeds;
* synthetic dataset configuration or checkpoint identity and checksum;
* training history and stopping behavior;
* prediction arrays in their documented units;
* station ordering, frequencies, coordinates, and adjacency where applicable;
* RMS definition and per-station diagnostics;
* ensemble or dropout uncertainty settings;
* all warnings and failures;
* figure paths and a serialized machine-readable result summary;
* reviewer, validation evidence, status, and date.

Do not rely only on the live ``inverter`` object stored in ``AgentResult.data``.
Save the supported checkpoint plus a plain configuration and manifest that can
be inspected without unpickling arbitrary objects.

Optional LLM interpretation
---------------------------

Passing ``api_key``, ``model``, and ``llm_provider`` enables an optional text
interpretation. Without an API key, science execution still runs and
``llm_interpretation`` remains ``None``.

Treat generated text as a draft. Verify every station count, RMS, depth,
geological claim, and recommendation against structured outputs. Do not send
sensitive project data to an external provider unless authorized by project
policy.

Failure handling
----------------

Use a consistent guard:

.. code-block:: python

   result = agent.execute(payload)

   if result.status == "failed":
       print("Failure:", result.error)
       print("Remediation:", result.error_fix_hint)
   elif result.status == "needs_review":
       print("Review required:", result.summary)
       for warning in result.warnings:
           print(" -", warning)
   else:
       for warning in result.warnings:
           print("Warning:", warning)
       # Continue to scientific validation.

Common failure causes include no installed deep-learning backend, no input
path or sites, unusable impedance, insufficient valid stations, incompatible
frequency coverage, disconnected graphs, checkpoint incompatibility, resource
exhaustion, and figure-save errors.

Agent success means the programmed workflow returned. It does not mean the
model passed scientific acceptance criteria.

Choosing the right interface
----------------------------

Use an agent when:

* the built-in synthetic training assumptions match the project;
* standard EDI-to-result orchestration is desired;
* consistent figures and ``AgentResult`` outputs are useful;
* a baseline or screening workflow is being established.

Use :mod:`pycsamt.ai.inversion` directly when:

* training priors or noise models must be customized;
* train, validation, calibration, and test splits require precise control;
* network or loss functions are being researched;
* checkpoint and optimizer behavior must be managed explicitly;
* custom metrics, callbacks, or distributed training are required;
* agent fallbacks are inappropriate for controlled production.

Use the workflow orchestrator only when multiple reviewed steps must be chained.
Understand each individual result contract before hiding it inside a larger
pipeline.

Review checklist
----------------

.. list-table::
   :header-rows: 1
   :widths: 31 69

   * - Check
     - Evidence
   * - Correct agent selected
     - Survey geometry, dimensionality, target, and architecture rationale.
   * - Input data reviewed
     - QC, components, frequency coverage, coordinates, processing, and usable
       station count.
   * - Training provenance retained
     - Priors, noise, solver, sample count, split, seed, epochs, and history.
   * - Checkpoint identity proven
     - Model name/path, metadata, checksum, compatibility, and no silent
       fallback.
   * - Result contract inspected
     - Exact status, warnings, units, shapes, station order, figures, and paths.
   * - Response fit reviewed
     - Per-station diagnostics, RMS definition, phase/component limitations,
       and structured residuals.
   * - Spatial assumptions checked
     - 2-D dimensionality or 3-D coordinates, radius, adjacency, and
       connectivity.
   * - Uncertainty interpreted conditionally
     - Ensemble/dropout method, calibration evidence, domain shift, and omitted
       sources.
   * - Independent validation completed
     - Classical baseline, synthetic holdout, field response, boreholes, and
       geological consistency.
   * - Release is auditable
     - Configuration, checkpoint, outputs, manifest, reviewer, status, and
       limitations.

Common mistakes
---------------

Avoid these errors:

* treating agent orchestration as a substitute for AI model validation;
* using default synthetic priors without checking field representativeness;
* ignoring a successful result's warnings;
* assuming ``best_model`` means lowest-error station;
* reporting the 1-D log-resistivity RMS as a full impedance RMS;
* calling U-Net lateral continuity proof of 2-D earth structure;
* using latitude/longitude-derived graph positions without reviewing units and
  station order;
* treating MC dropout or ensemble spread as total uncertainty;
* claiming pretrained inference after a fallback training run;
* trusting generated LLM interpretation without checking structured outputs;
* preserving figures but not the configuration, checkpoint, or data contract;
* continuing automatically after ``needs_review`` in a controlled workflow.

Next steps
----------

Continue with:

* :doc:`data_preparation` for representative training and field datasets;
* :doc:`training` for controlled lower-level model fitting;
* :doc:`inference` for checkpoint compatibility and field prediction;
* :doc:`validation` for acceptance tests and classical baselines;
* :doc:`uncertainty` for predictive calibration and domain shift;
* :doc:`reporting` for model cards and release packages;
* :doc:`../agents/index` for the complete pyCSAMT agent architecture.
