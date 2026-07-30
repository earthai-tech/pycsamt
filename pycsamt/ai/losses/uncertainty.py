# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Uncertainty-aware losses for calibrated inversion outputs.

These losses support the M9 uncertainty milestone in the
AI-inversion plan: a heteroscedastic aleatoric training term
(:func:`gaussian_nll_loss`) and a coverage-calibration diagnostic
(:func:`calibration_loss`).  The two are kept separate rather than
fused into one combined call: the NLL is a per-cell term evaluated
every training step, while calibration summarizes predictive
intervals across many held-out realizations and is evaluated
periodically, not differentiated through.

All functions operate on plain NumPy arrays so the module stays
importable without an optional deep-learning backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ._common import _reduce, _weight_array

__all__ = [
    "UncertaintyLossResult",
    "gaussian_nll_loss",
    "calibration_loss",
    "UncertaintyLoss",
]

_KINDS = ("gaussian_nll", "calibration")
_REDUCTIONS = ("mean", "sum")
_LOG_TWO_PI = float(np.log(2.0 * np.pi))


@dataclass(frozen=True)
class UncertaintyLossResult:
    """Immutable scalar result of an uncertainty-aware loss.

    Parameters
    ----------
    value : float
        Reduced loss value. ``nan`` when no cell was included and
        ``reduction="mean"``.
    kind : {"gaussian_nll", "calibration"}
        Loss family that produced ``value``.
    reduction : {"mean", "sum"}
        Reduction applied over valid, weighted cells.
    n_valid : int
        Number of cells included after masking.
    weight_sum : float
        Sum of weights over included cells. Equals ``n_valid`` when
        no explicit weights were supplied.

    Examples
    --------
    >>> import numpy as np
    >>> pred = np.array([0.0, 1.0])
    >>> true = np.array([0.0, 0.0])
    >>> log_var = np.array([0.0, 0.0])
    >>> result = gaussian_nll_loss(pred, true, log_var)
    >>> result.kind, result.n_valid
    ('gaussian_nll', 2)
    """

    value: float
    kind: str
    reduction: str
    n_valid: int
    weight_sum: float

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise ValueError(f"kind must be one of {_KINDS}.")
        if self.reduction not in _REDUCTIONS:
            raise ValueError(f"reduction must be one of {_REDUCTIONS}.")
        if self.n_valid < 0:
            raise ValueError("n_valid must be non-negative.")
        object.__setattr__(self, "value", float(self.value))
        object.__setattr__(self, "weight_sum", float(self.weight_sum))
        object.__setattr__(self, "n_valid", int(self.n_valid))


def gaussian_nll_loss(
    y_pred: Any,
    y_true: Any,
    log_variance: Any,
    *,
    valid: Any | None = None,
    weights: Any | None = None,
    reduction: str = "mean",
) -> UncertaintyLossResult:
    """Compute a heteroscedastic Gaussian negative log-likelihood.

    Each cell contributes
    ``0.5 * ((y_pred - y_true)**2 / variance + log_variance
    + log(2*pi))`` with ``variance = exp(log_variance)``.
    Parameterizing the log-variance rather than the variance itself
    keeps it unconstrained in sign while ``variance`` stays positive.

    Parameters
    ----------
    y_pred : array-like
        Predicted mean values.
    y_true : array-like
        True values, same shape as ``y_pred``.
    log_variance : array-like
        Predicted log-variance, same shape as ``y_pred``. Any finite
        real value is valid.
    valid : array-like of bool or None, optional
        Explicit cell mask, combined with finite-value masking of
        all three inputs.
    weights : array-like or None, optional
        Non-negative per-cell weights broadcastable to ``y_pred``.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over included cells.

    Returns
    -------
    UncertaintyLossResult
        Reduced negative log-likelihood.

    Examples
    --------
    >>> import numpy as np
    >>> pred = np.array([0.0, 1.0])
    >>> true = np.array([0.0, 0.0])
    >>> log_var = np.array([0.0, 0.0])
    >>> round(gaussian_nll_loss(pred, true, log_var).value, 6)
    1.168939
    """
    pred = np.asarray(y_pred, dtype=float)
    true = np.asarray(y_true, dtype=float)
    log_var = np.asarray(log_variance, dtype=float)
    if pred.shape != true.shape or pred.shape != log_var.shape:
        raise ValueError(
            "y_pred, y_true, and log_variance must share one shape."
        )
    if pred.size == 0:
        raise ValueError("y_pred, y_true, and log_variance must not be empty.")
    mask = np.isfinite(pred) & np.isfinite(true) & np.isfinite(log_var)
    if valid is not None:
        supplied = np.asarray(valid, dtype=bool)
        if supplied.shape != pred.shape:
            raise ValueError("valid must have the same shape as y_pred.")
        mask &= supplied
    safe_log_var = np.where(mask, log_var, 0.0)
    variance = np.exp(safe_log_var)
    difference = np.where(mask, pred - true, 0.0)
    per_cell = 0.5 * (
        np.square(difference) / variance + safe_log_var + _LOG_TWO_PI
    )
    per_cell = np.where(mask, per_cell, 0.0)
    weight_array = _weight_array(weights, pred.shape, "weights")
    active_weights = (
        mask.astype(float) if weight_array is None else weight_array * mask
    )
    value, weight_sum = _reduce(per_cell, active_weights, reduction)
    return UncertaintyLossResult(
        value=value,
        kind="gaussian_nll",
        reduction=reduction,
        n_valid=int(np.count_nonzero(mask)),
        weight_sum=weight_sum,
    )


def calibration_loss(
    coverage: Any,
    nominal_levels: Any,
    *,
    kind: str = "l2",
    valid: Any | None = None,
    weights: Any | None = None,
    reduction: str = "mean",
) -> UncertaintyLossResult:
    """Penalize deviation between empirical and nominal coverage.

    Parameters
    ----------
    coverage : array-like
        Empirical coverage observed at each nominal level, in
        ``[0, 1]``.
    nominal_levels : array-like
        Declared confidence levels in ``[0, 1]``, same shape as
        ``coverage``.
    kind : {"l1", "l2"}, default="l2"
        Elementwise penalty applied to each calibration residual.
    valid : array-like of bool or None, optional
        Explicit level mask, combined with finite-value masking of
        both inputs.
    weights : array-like or None, optional
        Non-negative per-level weights broadcastable to ``coverage``.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over included levels.

    Returns
    -------
    UncertaintyLossResult
        Reduced calibration penalty.

    Examples
    --------
    >>> import numpy as np
    >>> coverage = np.array([0.4, 0.9])
    >>> nominal = np.array([0.5, 0.8])
    >>> round(calibration_loss(coverage, nominal).value, 6)
    0.01
    """
    if kind not in ("l1", "l2"):
        raise ValueError("kind must be 'l1' or 'l2'.")
    cov = np.asarray(coverage, dtype=float)
    nominal = np.asarray(nominal_levels, dtype=float)
    if cov.shape != nominal.shape:
        raise ValueError("coverage and nominal_levels must share one shape.")
    if cov.size == 0:
        raise ValueError("coverage and nominal_levels must not be empty.")
    if np.any((cov[np.isfinite(cov)] < 0.0) | (cov[np.isfinite(cov)] > 1.0)):
        raise ValueError("coverage must be within [0, 1].")
    finite_nominal = nominal[np.isfinite(nominal)]
    if np.any((finite_nominal < 0.0) | (finite_nominal > 1.0)):
        raise ValueError("nominal_levels must be within [0, 1].")
    mask = np.isfinite(cov) & np.isfinite(nominal)
    if valid is not None:
        supplied = np.asarray(valid, dtype=bool)
        if supplied.shape != cov.shape:
            raise ValueError("valid must have the same shape as coverage.")
        mask &= supplied
    difference = np.where(mask, cov - nominal, 0.0)
    per_cell = np.abs(difference) if kind == "l1" else np.square(difference)
    weight_array = _weight_array(weights, cov.shape, "weights")
    active_weights = (
        mask.astype(float) if weight_array is None else weight_array * mask
    )
    value, weight_sum = _reduce(per_cell, active_weights, reduction)
    return UncertaintyLossResult(
        value=value,
        kind="calibration",
        reduction=reduction,
        n_valid=int(np.count_nonzero(mask)),
        weight_sum=weight_sum,
    )


@dataclass(frozen=True)
class UncertaintyLoss:
    """Configurable, callable heteroscedastic Gaussian NLL loss.

    Wraps :func:`gaussian_nll_loss` for reuse across training
    batches. Calibration is scored separately with
    :func:`calibration_loss` over binned coverage, since it
    summarizes predictive intervals across many held-out
    realizations rather than acting as a per-cell training term.

    Parameters
    ----------
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over included cells.

    Examples
    --------
    >>> import numpy as np
    >>> loss = UncertaintyLoss()
    >>> pred = np.array([0.0, 1.0])
    >>> true = np.array([0.0, 0.0])
    >>> log_var = np.array([0.0, 0.0])
    >>> round(loss(pred, true, log_var).value, 6)
    1.168939
    """

    reduction: str = "mean"

    def __post_init__(self) -> None:
        if self.reduction not in _REDUCTIONS:
            raise ValueError(f"reduction must be one of {_REDUCTIONS}.")

    def __call__(
        self,
        y_pred: Any,
        y_true: Any,
        log_variance: Any,
        *,
        valid: Any | None = None,
        weights: Any | None = None,
    ) -> UncertaintyLossResult:
        """Evaluate the configured Gaussian NLL on one batch."""
        return gaussian_nll_loss(
            y_pred,
            y_true,
            log_variance,
            valid=valid,
            weights=weights,
            reduction=self.reduction,
        )
