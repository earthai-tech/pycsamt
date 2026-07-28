.. _ai_inversion_losses:

Staged training objective
=========================

:doc:`concepts` introduces a supervised :term:`objective function` in two
pieces, a :term:`model-space metric` :eq:`eq-ai-model-loss` and a
:term:`response-space metric` :eq:`eq-ai-data-loss`, and a PINN-style
combination with generic vertical, lateral, and graph
:term:`regularization` terms :eq:`eq-ai-pinn-objective`.
:mod:`pycsamt.ai.losses` is where that idea becomes five concrete,
independently testable functions rather than one monolithic training
step:

.. math::
   :label: eq-ai-losses-staged

   L = w_m L_{\mathrm{model}} + \lambda_x L_{\mathrm{grad}_x}
       + \lambda_z L_{\mathrm{grad}_z} + \lambda_{tv} L_{TV}
       + \lambda_d L_{\mathrm{response}},

with :math:`L_{\mathrm{model}}` from :mod:`~pycsamt.ai.losses.model`,
the gradient and total-variation terms from
:mod:`~pycsamt.ai.losses.spatial`, and :math:`L_{\mathrm{response}}`
from :mod:`~pycsamt.ai.losses.response`. Two further terms sit outside
:eq:`eq-ai-losses-staged` rather than being folded into it:
:mod:`~pycsamt.ai.losses.boundary` anchors cells the training data has
no sensitivity to, and :mod:`~pycsamt.ai.losses.uncertainty` scores a
network's declared confidence, which is not something a single scalar
data-fit term can express. Every function returns a small, immutable
result dataclass — value together with ``kind``, ``reduction``,
``n_valid``, and ``weight_sum`` — rather than a bare float, so a
near-zero loss can always be told apart from a term that quietly had
nothing left to evaluate after masking. All five submodules operate on
plain NumPy arrays; nothing here imports PyTorch or TensorFlow.

Weighting a residual by what it means
--------------------------------------

:func:`~pycsamt.ai.losses.model.model_l1_loss`,
:func:`~pycsamt.ai.losses.model.model_l2_loss`, and
:func:`~pycsamt.ai.losses.model.model_huber_loss` share one masked,
weighted reduction over a predicted and a true grid — L1 and L2 in the
usual sense, and Huber blending the two,

.. math::
   :label: eq-ai-losses-huber

   \ell_\delta(r) =
   \begin{cases}
     \tfrac12 r^2, & |r| \le \delta \\
     \delta\left(|r| - \tfrac12\delta\right), & |r| > \delta,
   \end{cases}

quadratic near zero and linear beyond the transition point
:math:`\delta`, so a handful of badly recovered cells cannot dominate
the gradient the way a pure L2 term would:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.ai.losses import model_l2_loss, model_huber_loss

   >>> model_l2_loss(np.array([1.0, 3.0]), np.array([1.0, 1.0])).value
   2.0
   >>> small = model_huber_loss(np.array([0.5]), np.array([0.0]), delta=1.0).value
   >>> large = model_huber_loss(np.array([5.0]), np.array([0.0]), delta=1.0).value
   >>> round(small, 3), round(large, 3)
   (0.125, 4.5)

A cell can be excluded three ways — a non-finite value in either
array, an explicit ``valid`` mask, or a zero ``weights`` entry — and
all three compose, since a masked-out cell should never silently
survive because a weight happened to be positive.
:func:`~pycsamt.ai.losses.model.depth_weights` builds one common
weight pattern directly: shallow cells outweigh deep ones,
proportional to :math:`1/(1+\text{depth index})` and normalized to sum
to one, matching the intuition that near-surface structure is both
easier to recover and more heavily sampled by the data than depth is:

.. code-block:: pycon

   >>> from pycsamt.ai.losses import ModelLoss, depth_weights
   >>> depth_weights(3)
   array([0.54545455, 0.27272727, 0.18181818])

   >>> loss = ModelLoss.with_depth_weights(n_depth=3, dimension=2, kind="l2")
   >>> loss.weights.shape
   (3, 1)
   >>> true = np.array([[1.0, 1.0, 1.0], [2.0, 2.2, 2.0], [3.0, 3.0, 3.5]])
   >>> pred = np.array([[1.1, 0.9, 1.0], [2.3, 2.0, 2.1], [2.5, 3.4, 3.6]])
   >>> loss(pred, true)
   ModelLossResult(value=0.04181818181818182, kind='l2', reduction='mean', n_valid=9, weight_sum=3.0)

The depth-weighted value is well below the unweighted
``model_l2_loss(pred, true).value`` of ``0.0644``, because most of
this section's error sits in the poorly recovered bottom row, which
:class:`~pycsamt.ai.losses.model.ModelLoss`'s depth weighting
deliberately discounts. :class:`~pycsamt.ai.losses.model.ModelLoss`
bundles a kind, reduction, and optional fixed weights into one
reusable callable, so a training loop configures it once outside the
loop rather than re-passing the same keyword arguments on every batch.

Penalizing structure, not just misfit
----------------------------------------

A network can minimize :math:`L_{\mathrm{model}}` and still produce a
section that is pixel-noisy between neighbouring cells that the data
cannot individually distinguish.
:func:`~pycsamt.ai.losses.spatial.gradient_smoothness_loss` penalizes
the first difference along one grid axis,

.. code-block:: pycon

   >>> from pycsamt.ai.losses import gradient_smoothness_loss, total_variation_loss
   >>> grid = np.array([[0.0, 1.0, 3.0], [0.0, 0.0, 0.0]])
   >>> gradient_smoothness_loss(grid, axis=1, kind="l1")
   SpatialLossResult(value=0.75, kind='l1', label='grad_axis1', reduction='mean', n_valid=4, weight_sum=4.0)

and :func:`~pycsamt.ai.losses.spatial.total_variation_loss` sums that
same penalty over every axis at once — the standard *anisotropic*
total-variation definition (a per-axis sum of directional gradients),
not an isotropic pointwise gradient norm. On the same predicted
section used above, the three terms of :eq:`eq-ai-losses-staged`'s
spatial part read differently depending on which direction the
structure actually varies in:

.. code-block:: pycon

   >>> gradient_smoothness_loss(pred, axis=0, kind="l2").value  # depth
   1.3516666666666666
   >>> gradient_smoothness_loss(pred, axis=1, kind="l2").value  # lateral
   0.16666666666666663
   >>> total_variation_loss(pred, kind="l1").value
   0.6916666666666668

The depth gradient is nearly ten times the lateral one here, because
``pred`` genuinely varies with depth by design — a reminder that
:math:`\lambda_z` should not, in general, be tuned to the same value as
:math:`\lambda_x` for a survey where vertical layering is expected and
lateral continuity is the assumption worth enforcing.
:class:`~pycsamt.ai.losses.spatial.SpatialLoss` combines all three into
:eq:`eq-ai-losses-staged`'s :math:`\lambda_x L_{\mathrm{grad}_x} +
\lambda_z L_{\mathrm{grad}_z} + \lambda_{tv} L_{TV}` for a canonical
2-D ``(z, x)`` grid, and returns the combined float directly rather
than a result dataclass, since it is a sum of differently-labelled
terms with no single ``kind`` left to report:

.. code-block:: pycon

   >>> from pycsamt.ai.losses import SpatialLoss
   >>> spatial = SpatialLoss(lambda_x=1.0, lambda_z=0.5, lambda_tv=0.1, kind="l2")
   >>> spatial(pred)
   0.9116666666666666

which is exactly ``1.0 * 0.1667 + 0.5 * 1.3517 + 0.1 * 0.6917`` from
the three terms computed above — :class:`SpatialLoss` is bookkeeping
over the same functions, not a separate implementation.

Anchoring what the training data cannot see
------------------------------------------------

Not every cell in a predicted section is constrained by the response
data at all. Air cells above topography, and the outer padding of a
mesh built wide enough to keep boundary effects away from the survey
footprint, have no sensitivity in the forward problem — a network is
free to put anything there unless something tells it otherwise.
:func:`~pycsamt.ai.losses.boundary.boundary_condition_loss` is a
masked data-fit loss against an explicit required ``target`` value,
restricted to a caller-supplied ``boundary_mask``: there is no implicit
default target, since an agent should never silently hide a physical
assumption like "air is very resistive" inside a keyword default.
A common source for that mask is
:meth:`~pycsamt.ai.geology.topography.TopographicSurface.air_mask`:

.. code-block:: pycon

   >>> from pycsamt.ai.losses import BoundaryLoss
   >>> air_mask = np.array([[True, True, True], [False, False, False], [False, False, False]])
   >>> boundary = BoundaryLoss(kind="l1")
   >>> boundary(pred, boundary_mask=air_mask, target=1.0)
   ModelLossResult(value=0.0666666666666667, kind='l1', reduction='mean', n_valid=3, weight_sum=3.0)

``target=1.0`` here stands in for a fixed log-resistivity air value;
the small residual reflects that the top row of ``pred`` was already
close to it. Internally
:func:`~pycsamt.ai.losses.boundary.boundary_condition_loss` reuses the
exact same L1/L2/Huber machinery as
:mod:`~pycsamt.ai.losses.model` against a broadcast target array, so
everything said about robust penalties and masking above applies here
unchanged — a boundary constraint is a data-fit loss against a target
that happens to be known in advance rather than observed.

Closing the loop back through the forward solver
------------------------------------------------------

:math:`L_{\mathrm{model}}` and the spatial terms only ever look at
resistivity cells; nothing yet checks whether the recovered section,
run back through a solver, actually reproduces the data.
:func:`~pycsamt.ai.losses.response.response_residual_loss` closes that
loop directly on complex impedance,

.. code-block:: pycon

   >>> from pycsamt.ai.losses import response_residual_loss
   >>> z_pred = np.array([1 + 1j, 2 + 2j])
   >>> z_obs = np.array([1 + 1j, 0 + 0j])
   >>> response_residual_loss(z_pred, z_obs, kind="l2").value
   4.0
   >>> response_residual_loss(z_pred, z_obs, errors=np.array([1.0, 2.0]), kind="l2").value
   1.0

dividing each residual by its declared standard error before the
penalty when one is supplied, matching the usual normalized-RMS EM
data-misfit convention — an observation sitting at its
:term:`error floor` contributes less to the loss than one whose small
declared uncertainty claims it should already be matched almost
exactly.
:func:`~pycsamt.ai.losses.response.response_loss_from_contracts` skips
the manual array bookkeeping and takes a
:class:`~pycsamt.forward.maxwell.contracts.ForwardResult` (see
:doc:`forward_physics`) and a
:class:`~pycsamt.ai.data.contracts.SurveyData` (see
:doc:`data_contracts`) directly, refusing to run if their stations,
components, or frequencies are not identically ordered — silently
interpolating or reordering either axis before a misfit calculation is
exactly the kind of survey-matching mistake this contract exists to
rule out:

.. code-block:: pycon

   >>> from pycsamt.ai.data import SurveyData
   >>> from pycsamt.forward.maxwell import ForwardResult, SolverDiagnostics
   >>> from pycsamt.ai.losses import response_loss_from_contracts

   >>> z = np.array([[[1 + 1j]]])
   >>> observed = SurveyData(z, [10.0], ["S1"], ["zxy"], [[0, 0]])
   >>> diagnostics = SolverDiagnostics([[True]], [[1]], [[0.0]], 0.01)
   >>> forward = ForwardResult(
   ...     "a" * 64, [10.0], ["S1"], ["zxy"], z, None, "demo", "1", diagnostics,
   ... )
   >>> response_loss_from_contracts(forward, observed).value
   0.0

This is the same pairing the two-family split in :doc:`domain_gap`
exists for: a network trained only against clean synthetic responses
has no reason to expect the noise, static shift, and distortion an
:math:`L_{\mathrm{response}}` term will encounter once ``observed``
comes from a real survey rather than another
:func:`~pycsamt.ai.training.dataset2d.generate_2d_maxwell_dataset`
sample.

Making declared uncertainty answerable to evidence
----------------------------------------------------------

A model that predicts a resistivity value without a credible error bar
is only half a result. :func:`~pycsamt.ai.losses.uncertainty.gaussian_nll_loss`
trains a heteroscedastic :term:`aleatoric uncertainty` head alongside
the mean prediction, with each cell contributing

.. math::
   :label: eq-ai-losses-gaussian-nll

   \ell_i = \tfrac12\left(\frac{(\hat y_i - y_i)^2}{\sigma_i^2}
   + \log\sigma_i^2 + \log 2\pi\right),
   \qquad \sigma_i^2 = \exp(\log\sigma_i^2),

parameterized through a predicted log-variance rather than the
variance itself, so the network output stays unconstrained in sign
while :math:`\sigma_i^2` is guaranteed positive after the exponential:

.. code-block:: pycon

   >>> from pycsamt.ai.losses import gaussian_nll_loss, calibration_loss
   >>> pred_u = np.array([1.0, 2.2, 2.9, 4.5])
   >>> true_u = np.array([1.1, 2.0, 3.0, 4.0])
   >>> log_var = np.log(np.array([0.05, 0.1, 0.2, 0.3]))
   >>> gaussian_nll_loss(pred_u, true_u, log_var)
   UncertaintyLossResult(value=0.09038918945783034, kind='gaussian_nll', reduction='mean', n_valid=4, weight_sum=4.0)

A low, well-behaved log-variance still needs checking against reality:
declaring :math:`\sigma` honestly is not the same as being right about
it. :func:`~pycsamt.ai.losses.uncertainty.calibration_loss` scores that
separately, penalizing the gap between empirical coverage measured on a
held-out :term:`calibration set` and the nominal confidence level each
coverage value was supposed to hit:

.. code-block:: pycon

   >>> coverage = np.array([0.42, 0.68, 0.83, 0.94])
   >>> nominal = np.array([0.5, 0.7, 0.8, 0.9])
   >>> calibration_loss(coverage, nominal)
   UncertaintyLossResult(value=0.0023249999999999968, kind='calibration', reduction='mean', n_valid=4, weight_sum=4.0)

The two stay separate functions rather than one combined call
deliberately: the Gaussian NLL is a per-cell term evaluated and
backpropagated every training step, while calibration summarizes
predictive intervals across many held-out realizations and is checked
periodically, not differentiated through.

Wiring the spatial terms into a trainable network
--------------------------------------------------------

Every function above operates on NumPy arrays, which is what keeps
:mod:`pycsamt.ai.losses` importable without an optional deep-learning
backend and lets every example on this page run without one — but a
NumPy array cannot participate in backpropagation.
:class:`~pycsamt.ai.inversion.EMInverter2D`'s ``fit`` method exposes
``lambda_x``, ``lambda_z``, and ``lambda_tv`` keywords that add the
lateral-gradient, depth-gradient, and total-variation terms to its
plain MSE data fit during PyTorch training, mirroring the math of
:eq:`eq-ai-losses-staged`'s spatial part exactly but computed directly
on normalized ``torch`` tensors rather than by calling
:mod:`~pycsamt.ai.losses.spatial` itself, since a training step needs
gradients to flow through the penalty and not just past it. All three
weights default to zero, which reproduces plain MSE training
unchanged, and the staged terms are currently only implemented for the
PyTorch backend: passing a nonzero weight while TensorFlow is active
raises ``NotImplementedError`` rather than silently ignoring the
request. :mod:`~pycsamt.ai.losses.response`,
:mod:`~pycsamt.ai.losses.boundary`, and
:mod:`~pycsamt.ai.losses.uncertainty` are not yet wired into that same
training step; they remain available today for offline scoring —
exactly the checks :doc:`scientific_validation` runs after training,
independently of whichever loss shaped it.
