.. _ai_inversion_scientific_validation:

Recovery, residual, and OOD diagnostics
=======================================

:doc:`losses` shapes what a network is trained to minimize;
:mod:`pycsamt.ai.validation` is where the result gets checked
afterwards, by someone other than the optimizer. It splits that
checking into four independent areas rather than one aggregate score:
a :term:`model-space metric` for known-truth synthetic recovery
(:mod:`~pycsamt.ai.validation.recovery`), a
:term:`response-space metric` broken down by where it comes from
(:mod:`~pycsamt.ai.validation.residuals`), whether a declared
predictive uncertainty deserves trust
(:mod:`~pycsamt.ai.validation.calibration`), and an
:term:`out-of-distribution diagnostic` for inputs the training
support never covered (:mod:`~pycsamt.ai.validation.ood`). Plotting
stays in :mod:`pycsamt.ai.plot`; every number a validation figure
draws comes from here. The older, still-active :doc:`validation`
page documents ``Inv2DAgent``/``Inv3DAgent``'s own separate,
built-in validation protocol; this package is the newer, modular
library described in :doc:`roadmap`, and — as the closing section
below shows — the two are already partly connected rather than
entirely separate.

Scoring recovery when the truth is known
-------------------------------------------

A field survey has no ground truth to compare against, but a
synthetic realization from :func:`~pycsamt.ai.training.dataset2d.generate_2d_maxwell_dataset`
does, which is exactly what makes
:func:`~pycsamt.ai.validation.recovery.recovery_report` possible:
masked RMSE and MAE reusing :mod:`~pycsamt.ai.losses.model`'s own
reductions, an :math:`R^2` against the true section's own variance,
a windowed structural similarity index, and a per-depth breakdown.
Standing in for a trained network's prediction with one of this
package's own realizations, lightly perturbed and made deliberately
worse at the bottom row:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.ai.geology import GeologyGrid
   >>> from pycsamt.ai.training.dataset2d import (
   ...     Maxwell2DDatasetConfig, generate_2d_maxwell_dataset,
   ... )
   >>> from pycsamt.ai.validation import recovery_report

   >>> grid = GeologyGrid.regular_2d(nx=10, nz=8, dx_m=300, dz_m=150)
   >>> config = Maxwell2DDatasetConfig(
   ...     dataset_id="sv-demo", grid=grid,
   ...     correlation_length_x_m=(600.0, 1200.0),
   ...     correlation_length_z_m=(200.0, 400.0),
   ...     frequencies_hz=np.logspace(0, 4, 12),
   ...     station_x_m=list(np.linspace(500.0, 2500.0, 10)),
   ...     n_realizations=6, seed=7,
   ...     log_resistivity_mean=2.0, log_resistivity_std=0.4,
   ...     validation_fraction=0.0, test_fraction=0.0,
   ... )
   >>> dataset = generate_2d_maxwell_dataset(config)
   >>> true_log_rho = np.log10(dataset.samples[0].resistivity_ohm_m)

   >>> rng = np.random.default_rng(11)
   >>> pred_log_rho = true_log_rho + rng.normal(scale=0.08, size=true_log_rho.shape)
   >>> pred_log_rho[-1, :] += 0.35  # a network that degrades at depth

   >>> report = recovery_report(pred_log_rho, true_log_rho, compute_ssim=True, ssim_window=3)
   >>> round(report.rmse, 4), round(report.mae, 4), round(report.r2, 4)
   (0.1533, 0.0971, 0.8531)
   >>> round(report.ssim, 4)
   0.8553
   >>> np.round(report.depth_rmse, 4)
   array([0.0728, 0.0538, 0.0872, 0.0522, 0.0717, 0.0793, 0.0411, 0.3954])

The global RMSE alone (0.153) hides exactly what
``report.depth_rmse`` shows plainly: seven of eight layers sit near
0.05-0.09, and the bottom layer alone is at 0.395 — the injected
depth-dependent bias would be invisible in one aggregate number.
:func:`~pycsamt.ai.validation.recovery.structural_similarity`, used
inside ``report.ssim`` above, follows the windowed
luminance/contrast/structure formulation of Wang et al. (2004),

.. math::
   :label: eq-ai-scival-ssim

   \mathrm{SSIM}(x, y) =
   \frac{(2\mu_x\mu_y + c_1)(2\sigma_{xy} + c_2)}
        {(\mu_x^2 + \mu_y^2 + c_1)(\sigma_x^2 + \sigma_y^2 + c_2)},

with :math:`\mu`, :math:`\sigma^2`, and :math:`\sigma_{xy}` the
windowed mean, variance, and covariance and :math:`c_1`, :math:`c_2`
small stabilizing constants proportional to the data range — it has
no defined masking rule, so it is skipped (``ssim=None`` in
:class:`~pycsamt.ai.validation.recovery.RecoveryReport`) whenever the
grid is partially masked or smaller than the requested window, rather
than silently computed over an incomplete comparison.
:func:`~pycsamt.ai.validation.recovery.depth_profile_rmse` and
:func:`~pycsamt.ai.validation.recovery.depth_profile_mae` are also
available standalone, for a report that only needs the per-layer
breakdown without the rest.

Breaking a response residual down by where it comes from
------------------------------------------------------------

Recovery only applies where a true model exists.
:mod:`~pycsamt.ai.validation.residuals` instead breaks
:mod:`~pycsamt.ai.losses.response`'s complex-impedance residual down
by station, frequency, and component — the diagnostic that still
works on a field survey with no synthetic truth at all, since it
compares a forward-simulated response against an *observed* one
rather than a predicted model against a true one. A residual that is
large everywhere means something different from one localized to a
handful of stations, and only the second shape is visible once it is
broken down:

.. code-block:: pycon

   >>> from pycsamt.ai.validation import response_residual_report

   >>> rng = np.random.default_rng(5)
   >>> n_sta, n_freq, n_comp = 4, 3, 2
   >>> z_true = (
   ...     rng.normal(scale=40, size=(n_sta, n_freq, n_comp))
   ...     + 1j * rng.normal(scale=40, size=(n_sta, n_freq, n_comp))
   ...     + (60 + 40j)
   ... )
   >>> z_pred = z_true.copy()
   >>> z_pred[2:, :, 1] += (15 + 10j)  # bias only zyx at stations S3/S4

   >>> report = response_residual_report(
   ...     z_pred, z_true, kind="l2",
   ...     station_names=["S1", "S2", "S3", "S4"],
   ...     frequencies_hz=[1.0, 10.0, 100.0],
   ...     components=["zxy", "zyx"],
   ... )
   >>> report.overall.value
   81.25000000000001
   >>> np.round(report.by_station, 3)
   array([  0. ,   0. , 162.5, 162.5])
   >>> np.round(report.by_frequency, 3)
   array([81.25, 81.25, 81.25])
   >>> np.round(report.by_component, 3)
   array([  0. , 162.5])

The overall value of 81.25 gives no hint that the true residual is
zero at two of four stations; ``by_station`` and ``by_component``
localize it exactly to ``S3``/``S4`` and ``zyx``, while
``by_frequency`` stays flat because the injected bias does not depend
on frequency — three different projections of the same underlying
array, each answering a different question about the same misfit.
:func:`~pycsamt.ai.validation.residuals.response_residual_report_from_contracts`
skips the manual array handling the same way
:func:`~pycsamt.ai.losses.response.response_loss_from_contracts` does
(see :doc:`losses`), taking a
:class:`~pycsamt.forward.maxwell.contracts.ForwardResult` and a
:class:`~pycsamt.ai.data.contracts.SurveyData` directly and refusing
to run on station, component, or frequency sets that are not
identically ordered.

Checking whether declared uncertainty deserves trust
----------------------------------------------------------

A network's declared confidence is only useful if it is honest.
:func:`~pycsamt.ai.validation.calibration.reliability_curve` checks
that directly: at each of several nominal confidence levels, what
fraction of a held-out :term:`calibration set`'s true values actually
fall inside the predicted ``mean +/- z(level) * std`` interval — the
same Gaussian parameterization
:func:`~pycsamt.ai.losses.uncertainty.gaussian_nll_loss` trains
against, so a network trained with that loss can be checked with the
metric that matches it. Two models with identical predicted means but
different declared standard deviations make the point concrete: one
whose declared spread matches the data's real scatter, and one that
is simply overconfident:

.. code-block:: pycon

   >>> from pycsamt.ai.validation import reliability_curve

   >>> rng = np.random.default_rng(21)
   >>> true = rng.normal(loc=0.0, scale=1.0, size=2000)
   >>> mean = np.zeros(2000)

   >>> curve_ok = reliability_curve(true, mean, np.ones(2000), levels=[0.5, 0.8, 0.9, 0.95])
   >>> np.round(curve_ok.coverage, 3)
   array([0.514, 0.821, 0.908, 0.961])
   >>> round(curve_ok.calibration.value, 5), curve_ok.sharpness
   (0.0002, 1.0)

   >>> curve_bad = reliability_curve(true, mean, np.full(2000, 0.5), levels=[0.5, 0.8, 0.9, 0.95])
   >>> np.round(curve_bad.coverage, 3)
   array([0.27 , 0.496, 0.606, 0.681])
   >>> round(curve_bad.calibration.value, 5), curve_bad.sharpness
   (0.07603, 0.5)

Sharpness alone would call the second model better: a declared
standard deviation of 0.5 is a tighter, more confident interval than
1.0. Its coverage tells the opposite story — a nominal 90% interval
only actually contains the true value 61% of the time, because the
declared spread understates the data's real scatter by exactly the
factor it was shrunk by.
:func:`~pycsamt.ai.validation.calibration.predictive_sharpness` is
therefore documented, and should be read, as meaningful only
alongside good calibration, never as a quality signal by itself.

Flagging predictions the training distribution never covered
--------------------------------------------------------------------

None of the checks above say anything about an input a network was
never trained near in the first place.
:func:`~pycsamt.ai.validation.ood.ood_score` measures that directly
against a reference set of training feature vectors, either by
Mahalanobis distance from the reference mean and covariance,

.. math::
   :label: eq-ai-scival-mahalanobis

   d(\mathbf x) = \sqrt{(\mathbf x - \boldsymbol\mu)^{\!\top}
   \boldsymbol\Sigma^{-1} (\mathbf x - \boldsymbol\mu)},

which assumes an elliptical reference distribution and needs more
reference samples than features to invert
:math:`\boldsymbol\Sigma`, or by plain Euclidean distance to the
:math:`k`-th nearest reference point, which makes no distributional
assumption at all.
:func:`~pycsamt.ai.validation.ood.flag_out_of_distribution` turns
either score into a flag by comparing against a quantile of the
reference set's own leave-one-out self-scores — "how unusual is a
typical reference point" — rather than an arbitrary fixed number.
Using each realization's own depth- and lateral-roughness (the
standard deviation of its first difference along each axis) as a
two-feature summary of geological character, a realization from the
same :class:`~pycsamt.ai.training.dataset2d.Maxwell2DDatasetConfig`
used above sits well inside the reference set, while one generated
with correlation lengths ten times shorter — a visibly rougher,
structurally different field — does not:

.. code-block:: pycon

   >>> from pycsamt.ai.validation import flag_out_of_distribution

   >>> reference = np.array([
   ...     [np.diff(np.log10(s.resistivity_ohm_m), axis=0).std(),
   ...      np.diff(np.log10(s.resistivity_ohm_m), axis=1).std()]
   ...     for s in dataset.samples
   ... ])
   >>> short_config = Maxwell2DDatasetConfig(
   ...     dataset_id="sv-ood-shortcorr", grid=grid,
   ...     correlation_length_x_m=(60.0, 90.0),
   ...     correlation_length_z_m=(20.0, 40.0),
   ...     frequencies_hz=np.logspace(0, 4, 12),
   ...     station_x_m=list(np.linspace(500.0, 2500.0, 10)),
   ...     n_realizations=1, seed=100,
   ...     log_resistivity_mean=2.0, log_resistivity_std=0.4,
   ...     validation_fraction=0.0, test_fraction=0.0,
   ... )
   >>> short_dataset = generate_2d_maxwell_dataset(short_config)
   >>> short_rho = np.log10(short_dataset.samples[0].resistivity_ohm_m)
   >>> x = np.array([
   ...     reference[0],
   ...     [np.diff(short_rho, axis=0).std(), np.diff(short_rho, axis=1).std()],
   ... ])

   >>> report = flag_out_of_distribution(x, reference, method="mahalanobis", quantile=0.9)
   >>> np.round(report.scores, 3)
   array([ 0.952, 33.178])
   >>> round(report.threshold, 3)
   1.524
   >>> report.flagged.tolist()
   [False, True]

A Mahalanobis distance of 33.2 against a threshold of 1.5 is not a
subtle case — the short-correlation field's roughness lies far
outside anything the six reference realizations sampled, which is
exactly the situation :doc:`domain_gap` warns a network never sees
during training and :func:`flag_out_of_distribution` exists to catch
before a confident-looking, unsupported map reaches a report.

What actually reaches an agent's result today
----------------------------------------------------

:func:`~pycsamt.ai.validation.recovery.recovery_report` is not only
documented here — it already runs inside ``Inv2DAgent``'s
``physics="mt2d"`` path. When a held-out synthetic test or validation
partition is available, the agent predicts on it, compares against
its known-truth resistivity with exactly this function, and folds the
averaged RMSE/MAE/:math:`R^2` into its result, precisely because a
field survey has no ground truth to run the same check against — a
synthetic partition is the only place this class of check can run at
all. :mod:`~pycsamt.ai.validation.residuals`,
:mod:`~pycsamt.ai.validation.calibration`, and
:mod:`~pycsamt.ai.validation.ood` are not wired into that same
automatic result yet; they are available today the same way
:doc:`losses`'s response, boundary, and uncertainty terms are —
correct, tested, and ready to be called explicitly on a
:class:`~pycsamt.forward.maxwell.contracts.ForwardResult`,
:class:`~pycsamt.ai.data.contracts.SurveyData`, or reference feature
set, rather than produced automatically by every inversion run.
