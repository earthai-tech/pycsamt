# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Calibrated uncertainty quantification for EM neural-network inverters.

Two complementary post-hoc calibration strategies are provided, both of
which convert the *empirical* uncertainty from
:class:`~pycsamt.ai.inversion.ensemble.EnsembleInverter` or
:class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D` into *statistically
rigorous* intervals:

1. :class:`ConformalPredictor` — split conformal prediction
   (Vovk et al. 2005; Angelopoulos & Bates 2023).
   Guarantees exact marginal coverage :math:`P(y \\in \\hat{C}(x)) \\ge 1-\alpha`
   without distributional assumptions.

2. :class:`PosteriorCalibrator` — Gaussianization normalising flow
   (Kuleshov et al. 2018).
   Fits a monotone recalibration map :math:`g: [0,1] \\to [0,1]` via the
   Pool Adjacent Violators (PAVA) algorithm, then inverts it through the
   standard-normal quantile function to draw calibrated posterior samples.
   This constitutes a one-dimensional normalising flow that transforms the
   empirical residual distribution to a standard Gaussian — the asymptotic
   posterior.

Typical usage
-------------
>>> from pycsamt.ai.inversion import EnsembleInverter
>>> from pycsamt.ai.inversion.calibration import ConformalPredictor, PosteriorCalibrator
>>> ens = EnsembleInverter(base_estimator=inv, n_estimators=5)
>>> ens.fit(ds_train)                                           # doctest: +SKIP
>>> cp = ConformalPredictor(ens).calibrate(X_cal, y_cal)       # doctest: +SKIP
>>> center, lo, hi = cp.predict_intervals(X_test, alpha=0.10)  # doctest: +SKIP
>>> pc = PosteriorCalibrator().fit(y_cal, y_pred_cal, sigma_cal) # doctest: +SKIP
>>> samples = pc.predict_posterior(y_pred_test, sigma_test)     # doctest: +SKIP

References
----------
Vovk V., Gammerman A. & Shafer G. (2005) *Algorithmic Learning in a
Random World*. Springer.

Angelopoulos A. N. & Bates S. (2023) Conformal prediction: A gentle
introduction. *Foundations and Trends in Machine Learning*, **16**, 494–591.

Kuleshov V., Fenner N. & Ermon S. (2018) Accurate uncertainties for deep
learning using calibrated regression. *Proceedings of ICML*, PMLR 80,
2796–2804.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

__all__ = [
    "ConformalPredictor",
    "PosteriorCalibrator",
]


# ---------------------------------------------------------------------------
# Internal: Pool Adjacent Violators (isotonic regression, non-decreasing)
# ---------------------------------------------------------------------------


def _pava(y: np.ndarray) -> np.ndarray:
    """
    Non-decreasing isotonic regression via the Pool Adjacent Violators
    Algorithm (PAVA).

    Parameters
    ----------
    y : ndarray, shape (n,)
        Values to isotonically regress.

    Returns
    -------
    y_iso : ndarray, shape (n,)
        Isotonically non-decreasing version of *y*.
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n == 0:
        return y.copy()
    # Blocks: list of (mean, count)
    blocks: list = []
    for v in y:
        blocks.append([v, 1])
        while len(blocks) >= 2 and blocks[-2][0] > blocks[-1][0]:
            m1, c1 = blocks.pop()
            m2, c2 = blocks.pop()
            merged_mean = (m2 * c2 + m1 * c1) / (c1 + c2)
            blocks.append([merged_mean, c1 + c2])
    out = np.empty(n, dtype=float)
    idx = 0
    for mean, count in blocks:
        out[idx : idx + count] = mean
        idx += count
    return out


# ---------------------------------------------------------------------------
# ConformalPredictor
# ---------------------------------------------------------------------------


class ConformalPredictor:
    r"""
    Split conformal prediction wrapper for calibrated regression intervals.

    Given a fitted base predictor :math:`f(\mathbf{x})` that also exposes
    per-output uncertainty estimates :math:`\boldsymbol{\sigma}(\mathbf{x})`
    (e.g. :class:`~pycsamt.ai.inversion.ensemble.EnsembleInverter`), this
    class fits *normalised nonconformity scores* on a held-out calibration
    set and uses them to produce prediction bands with provable marginal
    coverage.

    **Normalised nonconformity score** (per sample, over all output dims):

    .. math::

        s_i = \max_j \frac{|y_{ij} - f(\mathbf{x}_i)_j|}
                          {\sigma_{ij} + \varepsilon}

    **Calibration quantile** (finite-sample correction, Vovk et al. 2005):

    .. math::

        \hat{q} = \text{Quantile}_{1-\alpha + \frac{1}{n_{\rm cal}+1}}
                  (s_1, \ldots, s_{n_{\rm cal}})

    **Prediction interval** for a new point :math:`\mathbf{x}`:

    .. math::

        \hat{C}_j(\mathbf{x}) =
        \bigl[f(\mathbf{x})_j - \hat{q}\,\sigma_j(\mathbf{x}),\;
              f(\mathbf{x})_j + \hat{q}\,\sigma_j(\mathbf{x})\bigr]

    **Marginal coverage guarantee** (Vovk et al. 2005):

    .. math::

        P\!\bigl(y_j \in \hat{C}_j(\mathbf{x})\bigr) \;\ge\; 1 - \alpha
        \quad \forall j

    Parameters
    ----------
    base_predictor : any
        A fitted object exposing
        ``predict_with_uncertainty(X) -> (mean, std)`` where both outputs
        have shape ``(n_samples, n_params)``.
    alpha : float
        Default coverage level ``1-alpha``.  Can be overridden per-call.
    eps : float
        Floor added to :math:`\boldsymbol{\sigma}` to prevent division by zero.

    Examples
    --------
    >>> cp = ConformalPredictor(ens_inverter)         # doctest: +SKIP
    >>> cp.calibrate(X_cal, y_cal)                    # doctest: +SKIP
    >>> center, lo, hi = cp.predict_intervals(X_new)  # doctest: +SKIP
    """

    def __init__(
        self,
        base_predictor,
        alpha: float = 0.10,
        eps: float = 1e-10,
    ) -> None:
        self.base_predictor = base_predictor
        self.alpha = float(alpha)
        self.eps = float(eps)

        self._scores: np.ndarray | None = None
        self._alpha_fit: float | None = None
        self._is_calibrated: bool = False

    # ── calibration ──────────────────────────────────────────────────────────

    def calibrate(
        self,
        X_cal: np.ndarray,
        y_cal: np.ndarray,
        alpha: float | None = None,
    ) -> ConformalPredictor:
        """
        Fit nonconformity scores on a held-out calibration set.

        Parameters
        ----------
        X_cal : ndarray, shape (n_cal, n_features)
        y_cal : ndarray, shape (n_cal, n_params)
        alpha : float or None
            Coverage level for this calibration run; defaults to
            ``self.alpha``.

        Returns
        -------
        self
        """
        alpha = float(alpha) if alpha is not None else self.alpha
        X_cal = np.asarray(X_cal, dtype=float)
        y_cal = np.asarray(y_cal, dtype=float)

        mean, sigma = self.base_predictor.predict_with_uncertainty(X_cal)
        mean = np.asarray(mean, dtype=float)
        sigma = np.asarray(sigma, dtype=float)

        # normalised residuals: max over output dimensions
        scores = np.max(
            np.abs(y_cal - mean) / (sigma + self.eps), axis=-1
        )  # (n_cal,)

        self._scores = np.sort(scores)
        self._alpha_fit = alpha
        self._is_calibrated = True
        return self

    # ── prediction ───────────────────────────────────────────────────────────

    def _q_hat(self, alpha: float) -> float:
        """Finite-sample conformal quantile at level ``1-alpha``."""
        n = len(self._scores)
        # Vovk (2005) finite-sample correction: ceil((n+1)*(1-alpha))/n
        q_level = np.ceil((n + 1) * (1.0 - alpha)) / n
        q_level = min(q_level, 1.0)
        return float(np.quantile(self._scores, q_level))

    def predict_intervals(
        self,
        X: np.ndarray,
        alpha: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return calibrated prediction intervals.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
        alpha : float or None
            Coverage level; defaults to the value used in :meth:`calibrate`.

        Returns
        -------
        center : ndarray, shape (n_samples, n_params)
        lower  : ndarray, shape (n_samples, n_params)
        upper  : ndarray, shape (n_samples, n_params)
            The interval :math:`[\\text{center} - \\hat{q}\\,\\sigma,\\;
            \\text{center} + \\hat{q}\\,\\sigma]` satisfies the marginal
            coverage guarantee.
        """
        self._check_calibrated()
        alpha = float(alpha) if alpha is not None else self._alpha_fit
        q = self._q_hat(alpha)

        mean, sigma = self.base_predictor.predict_with_uncertainty(X)
        mean = np.asarray(mean, dtype=float)
        sigma = np.asarray(sigma, dtype=float)

        half = q * sigma
        return mean, mean - half, mean + half

    # ── diagnostics ──────────────────────────────────────────────────────────

    def coverage(
        self,
        X_test: np.ndarray,
        y_test: np.ndarray,
        alpha: float | None = None,
    ) -> float:
        """
        Empirical marginal coverage on a test set.

        A well-calibrated predictor returns a value close to ``1 - alpha``.

        Parameters
        ----------
        X_test : ndarray
        y_test : ndarray, shape (n_test, n_params)
        alpha : float or None

        Returns
        -------
        coverage : float ∈ [0, 1]
        """
        self._check_calibrated()
        alpha = float(alpha) if alpha is not None else self._alpha_fit
        _, lo, hi = self.predict_intervals(X_test, alpha=alpha)
        y_test = np.asarray(y_test, dtype=float)
        covered = np.all((y_test >= lo) & (y_test <= hi), axis=-1)
        return float(covered.mean())

    def coverage_diagnostics(
        self,
        X_test: np.ndarray,
        y_test: np.ndarray,
        alphas: Sequence[float] | None = None,
    ) -> dict[float, float]:
        """
        Coverage table across a range of significance levels.

        Useful for a reliability diagram: plot ``actual_coverage`` vs
        ``nominal_coverage = 1 - alpha``.

        Parameters
        ----------
        X_test : ndarray
        y_test : ndarray
        alphas : sequence of float or None
            Significance levels to evaluate.  Defaults to 25 values
            uniformly spaced in ``[0.02, 0.50]``.

        Returns
        -------
        diag : dict  {alpha: actual_coverage}
        """
        self._check_calibrated()
        if alphas is None:
            alphas = np.linspace(0.02, 0.50, 25)
        return {
            float(a): self.coverage(X_test, y_test, alpha=a) for a in alphas
        }

    # ── misc ─────────────────────────────────────────────────────────────────

    def _check_calibrated(self) -> None:
        if not self._is_calibrated:
            raise RuntimeError(
                "ConformalPredictor is not calibrated.  Call calibrate() first."
            )

    def __repr__(self) -> str:
        status = "calibrated" if self._is_calibrated else "uncalibrated"
        n = len(self._scores) if self._scores is not None else 0
        return f"ConformalPredictor(alpha={self.alpha}, {status}, n_cal={n})"


# ---------------------------------------------------------------------------
# PosteriorCalibrator
# ---------------------------------------------------------------------------


class PosteriorCalibrator:
    r"""
    Gaussianization normalising flow for calibrated posterior sampling.

    Implements the regression calibration scheme of Kuleshov et al. (2018):
    a monotone recalibration function :math:`g: [0,1] \to [0,1]` is
    learned from calibration-set residuals via isotonic regression (PAVA),
    then composed with the standard-normal quantile function to produce
    draws from a calibrated posterior.

    **Step 1 — fit the recalibration map** on
    :math:`(y_{\rm cal}, \hat{y}_{\rm cal}, \sigma_{\rm cal})`:

    .. math::

        p_i = \Phi\!\left(\frac{y_i - \hat{y}_i}{\sigma_i}\right),
        \qquad
        g = \text{PAVA}\!\left(p_{(i)},\, \frac{i}{n_{\rm cal}}\right)

    where :math:`p_{(i)}` are the sorted predicted CDF values and
    :math:`i/n_{\rm cal}` are the empirical CDF targets.

    **Step 2 — calibrated posterior samples**:

    .. math::

        u \sim \mathcal{U}(0,1),\quad
        p^* = g^{-1}(u),\quad
        \theta^* = \hat{y} + \sigma \cdot \Phi^{-1}(p^*)

    The map :math:`g^{-1}` (inverse PAVA) is computed via linear
    interpolation.  Because the composition
    :math:`\Phi^{-1} \circ g^{-1} \circ \Phi` is a monotone bijection
    :math:`\mathbb{R} \to \mathbb{R}`, this constitutes a one-dimensional
    *normalising flow* that transforms the empirical residual distribution
    into a standard Gaussian.

    Parameters
    ----------
    n_bins : int
        Number of equally spaced CDF bins used as knots for the
        interpolated inverse map.  Higher values give a smoother flow at
        the cost of slight over-smoothing in data-sparse tails.

    Examples
    --------
    >>> pc = PosteriorCalibrator()                                  # doctest: +SKIP
    >>> pc.fit(y_cal, y_pred_cal, sigma_cal)                       # doctest: +SKIP
    >>> samples = pc.predict_posterior(y_pred_test, sigma_test)    # doctest: +SKIP
    >>> cal_std = pc.calibrated_std(sigma_raw_test)                # doctest: +SKIP

    References
    ----------
    Kuleshov V., Fenner N. & Ermon S. (2018) Accurate uncertainties for
    deep learning using calibrated regression.  *ICML*, PMLR 80, 2796–2804.
    """

    def __init__(self, n_bins: int = 100) -> None:
        self.n_bins = int(n_bins)

        # fitted quantities
        self._p_knots: np.ndarray | None = None  # sorted predicted CDF
        self._q_targets: np.ndarray | None = None  # empirical CDF targets
        self._sigma_scale: float | None = None  # global variance scale
        self._is_fitted: bool = False

    # ── fitting ──────────────────────────────────────────────────────────────

    def fit(
        self,
        y_cal: np.ndarray,
        y_pred_cal: np.ndarray,
        sigma_cal: np.ndarray,
    ) -> PosteriorCalibrator:
        """
        Fit the Gaussianization recalibration map.

        Parameters
        ----------
        y_cal : ndarray, shape (n_cal, n_params)
            True parameter vectors on the calibration set.
        y_pred_cal : ndarray, shape (n_cal, n_params)
            Point predictions (e.g. ensemble mean).
        sigma_cal : ndarray, shape (n_cal, n_params)
            Predicted standard deviations (e.g. ensemble std).

        Returns
        -------
        self
        """
        from scipy.special import (  # fast normal CDF / quantile
            ndtr,
        )

        y_cal = np.asarray(y_cal, dtype=float)
        y_pred_cal = np.asarray(y_pred_cal, dtype=float)
        sigma_cal = np.asarray(sigma_cal, dtype=float)

        eps = 1e-10
        # Standardised residuals: one scalar per (sample, param)
        z_cal = (y_cal - y_pred_cal) / (sigma_cal + eps)  # (n_cal, n_params)

        # Predicted CDF values — treat each (sample, param) as one observation
        p_pred = ndtr(z_cal).ravel()  # U[0,1] if calibrated

        # Sort and assign empirical CDF targets
        order = np.argsort(p_pred)
        p_sorted = p_pred[order]
        n_obs = len(p_sorted)
        p_target = (np.arange(1, n_obs + 1)) / n_obs  # empirical CDF

        # Isotonic regression: fit g such that g(p_sorted) ≈ p_target (monotone)
        g_vals = _pava(p_target)

        # Build interpolation knots at n_bins uniform p values
        p_knots = np.linspace(0.0, 1.0, self.n_bins + 2)
        g_knots = np.interp(p_knots, p_sorted, g_vals, left=0.0, right=1.0)
        g_knots = np.clip(g_knots, 1e-12, 1.0 - 1e-12)

        self._p_knots = p_knots
        self._g_knots = g_knots

        # Inverse map: from empirical CDF value to predicted CDF value
        # g^{-1}(u): linear interpolation on (g_knots, p_knots)
        # Sort by g_knots for invertibility (should already be monotone)
        sort_g = np.argsort(g_knots)
        self._g_inv_x = g_knots[sort_g]
        self._g_inv_y = p_knots[sort_g]

        # Global sigma scaling factor: ratio std(z_cal * sigma_scale) → 1
        # Equivalent to learning a global variance correction
        z_cal_flat = z_cal.ravel()
        self._sigma_scale = (
            float(np.std(z_cal_flat)) if np.std(z_cal_flat) > 0 else 1.0
        )

        self._is_fitted = True
        return self

    # ── prediction ───────────────────────────────────────────────────────────

    def predict_posterior(
        self,
        y_pred: np.ndarray,
        sigma_pred: np.ndarray,
        n_samples: int = 500,
        rng: np.random.Generator | None = None,
    ) -> np.ndarray:
        """
        Draw samples from the calibrated posterior via inverse-CDF sampling.

        Parameters
        ----------
        y_pred : ndarray, shape (n_samples_pred, n_params)
            Point predictions (ensemble mean).
        sigma_pred : ndarray, shape (n_samples_pred, n_params)
            Raw (uncalibrated) uncertainty estimates.
        n_samples : int
            Number of posterior draws per input point.
        rng : numpy.random.Generator or None
            Random number generator for reproducibility.

        Returns
        -------
        draws : ndarray, shape (n_samples, n_samples_pred, n_params)
            Calibrated posterior samples.  ``draws.mean(axis=0) ≈ y_pred``
            and ``draws.std(axis=0) ≈ calibrated_sigma``.
        """
        from scipy.special import ndtri

        self._check_fitted()
        if rng is None:
            rng = np.random.default_rng()

        y_pred = np.asarray(y_pred, dtype=float)
        sigma_pred = np.asarray(sigma_pred, dtype=float)
        n_pts, n_p = y_pred.shape

        # Sample u ~ U(0,1) then invert the calibration map g^{-1}
        u = rng.uniform(0.0, 1.0, size=(n_samples, n_pts, n_p))
        u_clipped = np.clip(u, 1e-12, 1.0 - 1e-12)
        p_star = np.interp(
            u_clipped.ravel(), self._g_inv_x, self._g_inv_y
        ).reshape(u.shape)
        p_star = np.clip(p_star, 1e-12, 1.0 - 1e-12)

        # Transform to standardised residual via inverse-normal
        z_star = ndtri(p_star)  # (n_samp, n_pts, n_p)

        # Calibrated sigma (apply global scale)
        sigma_cal = sigma_pred / self._sigma_scale  # (n_pts, n_p)

        # Posterior draws
        draws = y_pred[None, :, :] + sigma_cal[None, :, :] * z_star
        return draws

    def calibrated_std(self, sigma_raw: np.ndarray) -> np.ndarray:
        """
        Apply the global variance recalibration to raw sigma estimates.

        Parameters
        ----------
        sigma_raw : ndarray

        Returns
        -------
        sigma_cal : ndarray
            Recalibrated standard deviation such that coverage is
            approximately nominal on the calibration distribution.
        """
        self._check_fitted()
        return np.asarray(sigma_raw, dtype=float) / self._sigma_scale

    def calibration_error(
        self,
        y_test: np.ndarray,
        y_pred_test: np.ndarray,
        sigma_test: np.ndarray,
    ) -> float:
        """
        Mean absolute calibration error (MACE, Kuleshov et al. 2018).

        For a perfectly calibrated predictor, the fraction of test samples
        falling below the :math:`p`-quantile of the predicted distribution
        should equal :math:`p`.  MACE measures the average deviation:

        .. math::

            \\text{MACE} = \\frac{1}{B} \\sum_{b=1}^{B}
                |\\hat{F}(q_b) - q_b|

        where :math:`B` is the number of CDF evaluation points.

        Parameters
        ----------
        y_test : ndarray, shape (n_test, n_params)
        y_pred_test : ndarray, shape (n_test, n_params)
        sigma_test : ndarray, shape (n_test, n_params)

        Returns
        -------
        mace : float ∈ [0, 1]
            0 = perfectly calibrated, 1 = maximally mis-calibrated.
        """
        from scipy.special import ndtr

        self._check_fitted()
        y_test = np.asarray(y_test, dtype=float)
        y_pred_test = np.asarray(y_pred_test, dtype=float)
        sigma_test = np.asarray(sigma_test, dtype=float)

        sigma_cal = self.calibrated_std(sigma_test)
        z = (y_test - y_pred_test) / (sigma_cal + 1e-10)
        p_obs = ndtr(z).ravel()

        p_levels = np.linspace(0.0, 1.0, self.n_bins + 2)
        acal = np.array([np.mean(p_obs <= p) for p in p_levels])
        return float(np.mean(np.abs(acal - p_levels)))

    # ── misc ─────────────────────────────────────────────────────────────────

    def _check_fitted(self) -> None:
        if not self._is_fitted:
            raise RuntimeError(
                "PosteriorCalibrator is not fitted.  Call fit() first."
            )

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return f"PosteriorCalibrator(n_bins={self.n_bins}, {status})"
