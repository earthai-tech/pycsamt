AI And Model-Zoo Agents
=======================

The AI and model-zoo agents are the neural inversion layer of
:mod:`pycsamt.agents`. They turn a reviewed survey object into structured
:term:`AgentResult` outputs by loading stations, building impedance-derived
:term:`feature vector`\ s, training or loading a neural model, predicting a
resistivity model, computing diagnostic misfit, saving figures when requested,
and preserving warnings. This page gives the agent catalogue view; the longer
walk-through, with full run discussion and visual outputs, is in
:doc:`../ai_inversion/agents`.

These agents are powerful because they make baseline AI inversion repeatable,
but the repeatability comes from recording the full calculation. A useful run
identifies the input path or ``Sites`` object, station ordering, frequency
grid, feature definition, synthetic training settings, checkpoint identity,
backend, random seeds where exposed, output directory, warnings, and figure
paths. The prediction itself is usually stored in logarithmic resistivity:
:math:`m=\log_{10}\rho`, where :math:`\rho` is in ohm metres. That convention
keeps the learning target numerically stable and avoids treating a change from
10 to 100 ohm m as the same absolute error as a change from 1010 to 1100 ohm m.

Before choosing an AI agent, decide what kind of spatial assumption is
scientifically defensible. ``AIInversionAgent`` predicts one layered column per
station; ``Inv2DAgent`` predicts a profile section from a station-frequency
panel; ``Inv3DAgent`` lets neighbouring stations exchange information through
a graph; ``EnsembleAgent`` repeats the 1-D workflow to estimate predictive
spread; ``JointInversionAgent`` combines multiple modalities; and
``ModelZooAgent`` handles named checkpoints. None of them removes
:term:`non-uniqueness`, replaces :term:`quality control`, or proves that the
training distribution represents the field site.

.. _agent-ai-inversion:

AIInversionAgent
----------------

``AIInversionAgent`` runs end-to-end 1-D AI inversion. It is the right first
agent when the survey can be screened station by station and a
:term:`layered model` is an acceptable starting approximation. The agent builds
features from the observed :term:`impedance tensor`, normally
:math:`\log_{10}\rho_a` and phase for the ``xy`` component on a common
frequency grid, trains :class:`pycsamt.ai.inversion.EMInverter1D` on synthetic
1-D responses unless a checkpoint is supplied, predicts layer resistivities
and thicknesses, then forward-models the prediction to obtain a response-fit
RMS.

Mathematically, the learned mapping can be read as
:math:`\hat{\mathbf{m}}_s=f_\theta(\mathbf{x}_s)`, where station :math:`s`
has feature vector :math:`\mathbf{x}_s` and predicted model
:math:`\hat{\mathbf{m}}_s=(\log_{10}\rho_1,\ldots,\log_{10}\rho_L)`. The
diagnostic RMS is computed in log apparent-resistivity space after the
predicted layered model has been passed through the 1-D
:term:`forward operator`:

.. math::

   \mathrm{RMS}_s =
   \sqrt{\frac{1}{N_s}\sum_{i=1}^{N_s}
   \left(\log_{10}\rho_{a,i}^{\mathrm{obs}} -
         \log_{10}\rho_{a,i}^{\mathrm{pred}}\right)^2}.

This is a practical response-fit diagnostic. It is not a fully weighted
:term:`RMS misfit` unless the calculation is explicitly normalized by data
uncertainties.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import AIInversionAgent
   >>> agent = AIInversionAgent(
   ...     arch="resnet",
   ...     n_layers=5,
   ...     n_train_samples=2000,
   ...     epochs=30,
   ... )
   >>> result = agent.execute({
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "output_dir": "outputs/agents/ai1d",
   ... })
   >>> print(result.status)
   success
   >>> print(round(result.get("rms_global"), 3))
   3.9
   >>> print(sorted(result.get("figure_paths", {})))
   ['ai_section', 'convergence']

The ``ai_section`` figure is the main station-by-layer view of the predicted
model, while ``convergence`` shows whether training stopped because the model
continued to improve or because validation loss had flattened. Treat both as
part of the result: a low RMS is more convincing when the convergence curve is
stable and the section does not show isolated station artefacts.

.. _agent-inv2d:

Inv2DAgent
----------

``Inv2DAgent`` performs 2-D profile inversion using a U-Net style model. It
uses the whole station-frequency panel rather than treating stations as
independent columns, so it is useful when lateral continuity along a survey
line is part of the modelling question. The input tensor has shape similar to
:math:`(C,F,S)`, where :math:`C` is the feature channel count, :math:`F` the
frequency count, and :math:`S` the number of stations used by the profile
window. The predicted section is
:math:`\hat M\in\mathbb{R}^{D\times S}`, with :math:`D` depth cells.

The U-Net smoothness is learned from synthetic profile examples and from the
architecture. It should not be confused with the roughness penalty in a
classical 2-D inversion objective. The agent's data-space check maps predicted
depth cells back toward observed apparent resistivity using the Bostick-style
sampling depth

.. math::

   d_B(T) = 503\sqrt{\rho_a(T)\,T},

where :math:`T=1/f` is period. This is a quick diagnostic bridge between the
section and the observed curve, not a replacement for a full 2-D EM forward
solve.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import Inv2DAgent
   >>> result = Inv2DAgent(
   ...     n_depth=40,
   ...     n_freqs=32,
   ...     n_train_profiles=200,
   ...     n_stations_per_profile=10,
   ...     epochs=30,
   ... ).execute({
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "output_dir": "outputs/agents/inv2d",
   ... })
   >>> print(result.status)
   success
   >>> print(result.get("pred_section").shape)
   (40, 10)
   >>> print(result.get("station_names")[:3])
   ['18-001A', '18-002U', '18-003A']

The predicted section should be read as a learned profile image. Laterally
smooth features may reflect coherent structure, but they may also reflect the
U-Net training prior and station-window size. Compare the section with QC,
phase-tensor dimensionality, station spacing, and any classical 2-D inversion
before treating a continuous band as a geological boundary.

.. _agent-inv3d:

Inv3DAgent
----------

``Inv3DAgent`` performs spatial AI inversion with a graph-convolutional model.
Use it when station geometry matters and when neighbouring sites should inform
each other. The graph is controlled by either a supplied normalized adjacency
matrix or by a radius search on station coordinates. With raw adjacency
:math:`A`, the graph convolution uses self-loops and symmetric normalization,

.. math::

   \tilde A = A + I, \qquad
   \hat A = \tilde D^{-1/2}\tilde A\tilde D^{-1/2},

where :math:`\tilde D` is the degree matrix of :math:`\tilde A`. A larger
radius increases graph density, which can help isolated stations but can also
smooth across unrelated structures. Always inspect the station coordinates,
edge count, and degree distribution before interpreting the section.

If ``n_mc`` is positive, the agent also estimates :term:`Monte Carlo dropout`
spread. For layer :math:`\ell` and :math:`M` stochastic inference passes,

.. math::

   \sigma_\ell =
   \sqrt{\frac{1}{M}\sum_{m=1}^{M}
   \left(\hat y_\ell^{(m)}-\bar y_\ell\right)^2},
   \qquad
   \bar y_\ell=\frac{1}{M}\sum_{m=1}^{M}\hat y_\ell^{(m)}.

This is an :term:`epistemic uncertainty` diagnostic for the trained network,
not a complete uncertainty budget for acquisition error, prior mismatch, or
geological ambiguity.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import Inv3DAgent
   >>> result = Inv3DAgent(
   ...     n_layers=5,
   ...     n_freqs=32,
   ...     n_train_profiles=150,
   ...     epochs=30,
   ...     radius=250.0,
   ...     n_mc=20,
   ... ).execute({
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "output_dir": "outputs/agents/inv3d",
   ... })
   >>> print(result.status)
   success
   >>> print(result.get("pred_rho").shape)
   (28, 5)
   >>> degree = (result.get("adjacency") > 0).sum(axis=1) - 1
   >>> print(degree[:6])
   [1 2 2 3 2 3]

Review the graph products together. The resistivity section shows the
along-profile prediction, the depth slices show how the same prediction varies
spatially at selected layers, and the uncertainty map shows where dropout
passes disagree. A resistive or conductive feature is more credible when it is
consistent across these views and does not coincide with weak graph
connectivity or high uncertainty.

.. _agent-ensemble:

EnsembleAgent
-------------

``EnsembleAgent`` runs repeated 1-D neural inversions and summarizes prediction
spread. It is useful when a single network prediction is too brittle for the
decision being made, but it should be read carefully: an ensemble estimates
variation among the fitted members, not all uncertainty in the inverse
problem.

For a predicted parameter :math:`j`, the ensemble mean and standard deviation
are

.. math::

   \bar m_j = \frac{1}{K}\sum_{k=1}^{K}\hat m_{k,j},
   \qquad
   s_j = \sqrt{\frac{1}{K-1}\sum_{k=1}^{K}
   \left(\hat m_{k,j}-\bar m_j\right)^2}.

When conformal calibration is enabled, the interval is widened by a quantile
of calibration residuals. The guarantee only applies when the calibration
examples and future inputs are exchangeable; it does not automatically become
field coverage.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import EnsembleAgent
   >>> result = EnsembleAgent(
   ...     n_estimators=5,
   ...     n_train_samples=2000,
   ...     epochs=30,
   ...     calibrate=True,
   ... ).execute({
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "output_dir": "outputs/agents/ensemble",
   ... })
   >>> print(result.status)
   success
   >>> print(round(result.get("coverage"), 3))
   0.0

Low or zero empirical coverage is a scientific warning. It usually means the
ensemble members agree with each other more than they agree with the
calibration targets, so the reported intervals are too narrow for the tested
distribution.

.. _agent-joint-inversion:

JointInversionAgent
-------------------

``JointInversionAgent`` runs multi-modal inversion, for example MT with TEM or
CSAMT with a supporting data type. Use it only when the modalities are meant
to constrain the same subsurface target and their coordinate systems, depth
support, and preprocessing histories can be reconciled. In a joint objective,
the data terms normally take the form

.. math::

   \Phi(\mathbf{m}) =
   \sum_{q=1}^{Q}
   \left\|\mathbf{W}_q
   \left(\mathbf{d}_q-\mathbf{F}_q(\mathbf{m})\right)\right\|_2^2
   + \lambda R(\mathbf{m}),

where modality :math:`q` has observed data :math:`\mathbf{d}_q`, forward map
:math:`\mathbf{F}_q`, weighting matrix :math:`\mathbf{W}_q`, and shared model
:math:`\mathbf{m}`. The weights are part of the interpretation; a weakly
weighted modality may barely influence the result, while an over-weighted one
can dominate the inversion.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import JointInversionAgent
   >>> result = JointInversionAgent(modalities=["mt", "tem"]).execute({
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "secondary_path": "data/TEM/WILLY_LINE",
   ...     "output_dir": "outputs/agents/joint",
   ... })
   >>> print(result.status)
   success

Record the source and preprocessing chain for every modality. A joint result
without modality provenance is difficult to reproduce and easy to overstate.

.. _agent-anomaly-detection:

AnomalyDetectionAgent
---------------------

``AnomalyDetectionAgent`` flags anomalous data patterns before expensive
training or inversion. Use it for station-frequency screening, survey triage,
and manual review queues. A typical anomaly score compares a feature vector
with a learned or robust reference distribution. In standardized form, a
simple score can be written

.. math::

   a_s = \left\|\Sigma^{-1/2}(\mathbf{x}_s-\boldsymbol{\mu})\right\|_2,

where :math:`\boldsymbol{\mu}` and :math:`\Sigma` describe the reference
feature population. The practical interpretation is simple: high-scoring
stations deserve review before they are allowed to shape a neural inversion.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import AnomalyDetectionAgent
   >>> anomalies = AnomalyDetectionAgent().execute({
   ...     "path": "data/AMT/WILLY_DATA/L18PLT",
   ...     "output_dir": "outputs/agents/anomalies",
   ... })
   >>> print(anomalies.status)
   success
   >>> print(sorted(anomalies.data))
   ['anomaly_table', 'figures', 'figure_paths', 'summary']

Use the anomaly table and figures as a review queue. High-score stations or
frequency bands should be checked against raw impedance curves, acquisition
notes, dead-band behaviour, and neighbouring stations before they are used in
training or inversion.

.. _agent-model-zoo:

ModelZooAgent
-------------

``ModelZooAgent`` lists available model metadata, downloads checkpoints, and
runs predictions where supported. Listing the registry is local and
reproducible; downloading depends on whether the public weights have been
released and whether the cache already contains the checkpoint.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents import ModelZooAgent
   >>> zoo = ModelZooAgent()
   >>> models = zoo.execute({"action": "list"})
   >>> print(models.status)
   success
   >>> print(models.summary)
   5 pre-trained models in zoo.
   >>> [item["name"] for item in models["details"]]
   ['mt1d-resnet-5layer-v1', 'mt1d-cnn-5layer-v1', 'mt1d-resnet-7layer-v1', 'csamt1d-resnet-5layer-v1', 'tem1d-fcn-5layer-v1']

When a checkpoint is present, a zoo prediction is just a named version of the
1-D workflow with metadata attached. When it is absent, the current
implementation falls back to on-the-fly training and records the missing
checkpoint in warnings. Check ``status``, ``warnings``, and
``checkpoint_path`` before calling a result pretrained.

.. code-block:: pycon
   :linenos:

   >>> checkpoint = zoo.execute({
   ...     "action": "download",
   ...     "model_name": "mt1d-resnet-5layer-v1",
   ... })
   >>> print(checkpoint.status)
   needs_review
   >>> print(checkpoint.get("checkpoint_path"))
   None

AI Workflow Pattern
-------------------

.. code-block:: text

   MTLoaderAgent -> DataQCAgent -> DenoisingAgent
   -> AIInversionAgent or Inv2DAgent or Inv3DAgent
   -> InterpretationAgent -> ReportAgent

Add ``AnomalyDetectionAgent`` before inversion when the survey contains
questionable stations or frequency bands. Add ``EnsembleAgent`` when
uncertainty is part of the objective. Add ``JointInversionAgent`` when a
secondary modality is available and traceable. Add ``ModelZooAgent`` when a
named checkpoint is part of the reproducible workflow.

The human-readable flow is as important as the code order: first establish
that the data are usable, then choose the spatial assumption, then state the
neural target and feature contract, then run the agent, then preserve outputs,
figures, warnings, and provenance. That sequence makes the AI layer auditable
instead of merely fast.
