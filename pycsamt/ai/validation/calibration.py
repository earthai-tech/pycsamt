# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Uncertainty-calibration diagnostics for predictive intervals.

These diagnostics support the Uncertainty row of the validation
matrix in the AI-inversion plan: "calibration, coverage, sharpness,
OOD sensitivity".  Predictive distributions are assumed Gaussian,
i.e. a cell's predictive interval at nominal level ``p`` is
``mean +/- z(p) * std`` with ``z(p)`` the two-sided normal quantile;
this matches how ``gaussian_nll_loss`` in
:mod:`pycsamt.ai.losses.uncertainty` parameterizes aleatoric
uncertainty.

All functions operate on plain NumPy arrays so the module stays
importable without an optional deep-learning backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.stats import norm

from ..losses.uncertainty import UncertaintyLossResult, calibration_loss

__all__ = [
    "ReliabilityCurve",
    "reliability_curve",
    "empirical_coverage",
    "predictive_sharpness",
]

_DEFAULT_LEVELS = (0.5, 0.8, 0.9, 0.95, 0.99)


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    """Return a read-only copy of *value* as an ndarray."""
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _validate_inputs(
    y_true: Any, y_pred_mean: Any, y_pred_std: Any
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return y_true/y_pred_mean/y_pred_std as float arrays of one
    non-empty shape.
    """
    true = np.asarray(y_true, dtype=float)
    mean = np.asarray(y_pred_mean, dtype=float)
    std = np.asarray(y_pred_std, dtype=float)
    if true.shape != mean.shape or true.shape != std.shape:
        raise ValueError(
            "y_true, y_pred_mean, and y_pred_std must share one "
            "shape."
        )
    if true.size == 0:
        raise ValueError(
            "y_true, y_pred_mean, and y_pred_std must not be empty."
        )
    return true, mean, std


def _validate_levels(levels: Any | None) -> np.ndarray:
    """Return validated nominal confidence levels in ``(0, 1)``."""
    array = np.asarray(
        _DEFAULT_LEVELS if levels is None else levels, dtype=float
    )
    if array.ndim != 1 or array.size == 0:
        raise ValueError("levels must be a non-empty 1-D array.")
    if (
        not np.all(np.isfinite(array))
        or np.any(array <= 0.0)
        or np.any(array >= 1.0)
    ):
        raise ValueError("levels must be finite values within (0, 1).")
    return array


def _coverage_mask(
    true: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    valid: Any | None,
) -> np.ndarray:
    """Combine finite/positive-std masking with an optional mask."""
    mask = (
        np.isfinite(true)
        & np.isfinite(mean)
        & np.isfinite(std)
        & (std > 0.0)
    )
    if valid is not None:
        supplied = np.asarray(valid, dtype=bool)
        if supplied.shape != true.shape:
            raise ValueError("valid must have the same shape as y_true.")
        mask &= supplied
    return mask


@dataclass(frozen=True)
class ReliabilityCurve:
    """Immutable calibration report for Gaussian predictive
    intervals.

    Parameters
    ----------
    levels : ndarray
        Nominal confidence levels in ``(0, 1)``.
    coverage : ndarray
        Empirical coverage at each level, same shape as ``levels``.
    calibration : UncertaintyLossResult
        Reduced deviation between ``coverage`` and ``levels``, from
        :func:`~pycsamt.ai.losses.uncertainty.calibration_loss`.
    sharpness : float
        Mean predictive standard deviation over included cells.
        Lower is sharper (more confident); meaningful only alongside
        good calibration.
    n_valid : int
        Number of cells included after masking.
    shape : tuple of int
        Shape of the compared ``y_true`` array.

    Examples
    --------
    >>> import numpy as np
    >>> true = np.array([0.0, 0.0, 0.0, 0.0, 10.0])
    >>> mean = np.zeros(5)
    >>> std = np.ones(5)
    >>> curve = reliability_curve(true, mean, std)
    >>> curve.n_valid
    5
    """

    levels: np.ndarray
    coverage: np.ndarray
    calibration: UncertaintyLossResult
    sharpness: float
    n_valid: int
    shape: tuple[int, ...]

    def __post_init__(self) -> None:
        levels = np.asarray(self.levels, dtype=float)
        coverage = np.asarray(self.coverage, dtype=float)
        if levels.shape != coverage.shape:
            raise ValueError("levels and coverage must share one shape.")
        if not isinstance(self.calibration, UncertaintyLossResult):
            raise TypeError(
                "calibration must be an UncertaintyLossResult."
            )
        object.__setattr__(self, "levels", _readonly(levels))
        object.__setattr__(self, "coverage", _readonly(coverage))
        object.__setattr__(self, "sharpness", float(self.sharpness))
        object.__setattr__(self, "n_valid", int(self.n_valid))
        object.__setattr__(
            self, "shape", tuple(int(size) for size in self.shape)
        )


def empirical_coverage(
    y_true: Any,
    y_pred_mean: Any,
    y_pred_std: Any,
    *,
    levels: Any | None = None,
    valid: Any | None = None,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Compute empirical coverage of Gaussian predictive intervals.

    Parameters
    ----------
    y_true : array-like
        True values.
    y_pred_mean, y_pred_std : array-like
        Predicted mean and positive standard deviation, same shape
        as ``y_true``.
    levels : array-like or None, optional
        Nominal confidence levels in ``(0, 1)``. Defaults to
        ``(0.5, 0.8, 0.9, 0.95, 0.99)``.
    valid : array-like of bool or None, optional
        Explicit cell mask, combined with finite-value masking of
        all three inputs and with ``y_pred_std > 0``.

    Returns
    -------
    levels : ndarray
        The validated nominal levels.
    coverage : ndarray, same shape as ``levels``
        Fraction of included cells whose true value falls inside the
        ``mean +/- z(level) * std`` interval.
    n_valid : int
        Number of cells included after masking.

    Examples
    --------
    >>> import numpy as np
    >>> true = np.array([0.0, 0.0, 0.0, 0.0, 10.0])
    >>> mean = np.zeros(5)
    >>> std = np.ones(5)
    >>> levels, coverage, n_valid = empirical_coverage(
    ...     true, mean, std, levels=[0.5]
    ... )
    >>> coverage.tolist()
    [0.8]
    """
    true, mean, std = _validate_inputs(y_true, y_pred_mean, y_pred_std)
    mask = _coverage_mask(true, mean, std, valid)
    if not np.any(mask):
        raise ValueError("no valid cell to evaluate.")
    levels_arr = _validate_levels(levels)
    residual = np.abs(true[mask] - mean[mask]) / std[mask]
    coverage = np.array(
        [
            float(np.mean(residual <= norm.ppf(0.5 + level / 2.0)))
            for level in levels_arr
        ]
    )
    return levels_arr, coverage, int(np.count_nonzero(mask))


def predictive_sharpness(
    y_pred_std: Any, *, valid: Any | None = None
) -> float:
    """Return the masked mean predictive standard deviation.

    Sharpness summarizes how confident a model's predictive
    distribution is, independent of correctness; it is only a
    meaningful quality signal alongside good calibration.

    Parameters
    ----------
    y_pred_std : array-like
        Predicted positive standard deviation.
    valid : array-like of bool or None, optional
        Explicit cell mask, combined with finite-value masking and
        with ``y_pred_std > 0``.

    Returns
    -------
    float
        Mean predictive standard deviation over included cells.

    Examples
    --------
    >>> import numpy as np
    >>> predictive_sharpness(np.array([1.0, 2.0, 3.0]))
    2.0
    """
    std = np.asarray(y_pred_std, dtype=float)
    if std.size == 0:
        raise ValueError("y_pred_std must not be empty.")
    mask = np.isfinite(std) & (std > 0.0)
    if valid is not None:
        supplied = np.asarray(valid, dtype=bool)
        if supplied.shape != std.shape:
            raise ValueError(
                "valid must have the same shape as y_pred_std."
            )
        mask &= supplied
    if not np.any(mask):
        raise ValueError("no valid cell to evaluate.")
    return float(np.mean(std[mask]))


def reliability_curve(
    y_true: Any,
    y_pred_mean: Any,
    y_pred_std: Any,
    *,
    levels: Any | None = None,
    valid: Any | None = None,
    kind: str = "l2",
    reduction: str = "mean",
) -> ReliabilityCurve:
    """Build a full calibration report for Gaussian predictive
    intervals.

    Parameters
    ----------
    y_true : array-like
        True values.
    y_pred_mean, y_pred_std : array-like
        Predicted mean and positive standard deviation, same shape
        as ``y_true``.
    levels : array-like or None, optional
        Nominal confidence levels in ``(0, 1)``. Defaults to
        ``(0.5, 0.8, 0.9, 0.95, 0.99)``.
    valid : array-like of bool or None, optional
        Explicit cell mask, combined with finite-value masking of
        all three inputs and with ``y_pred_std > 0``.
    kind : {"l1", "l2"}, default="l2"
        Elementwise penalty forwarded to
        :func:`~pycsamt.ai.losses.uncertainty.calibration_loss`.
    reduction : {"mean", "sum"}, default="mean"
        Reduction forwarded to
        :func:`~pycsamt.ai.losses.uncertainty.calibration_loss`.

    Returns
    -------
    ReliabilityCurve
        Combined coverage, calibration penalty, and sharpness.

    Examples
    --------
    >>> import numpy as np
    >>> true = np.array([0.0, 0.0, 0.0, 0.0, 10.0])
    >>> mean = np.zeros(5)
    >>> std = np.ones(5)
    >>> curve = reliability_curve(true, mean, std, levels=[0.5])
    >>> curve.coverage.tolist(), curve.sharpness
    ([0.8], 1.0)
    """
    true, mean, std = _validate_inputs(y_true, y_pred_mean, y_pred_std)
    mask = _coverage_mask(true, mean, std, valid)
    if not np.any(mask):
        raise ValueError("no valid cell to evaluate.")
    levels_arr, coverage, n_valid = empirical_coverage(
        true, mean, std, levels=levels, valid=mask
    )
    sharpness = predictive_sharpness(std, valid=mask)
    calibration = calibration_loss(
        coverage, levels_arr, kind=kind, reduction=reduction
    )
    return ReliabilityCurve(
        levels=levels_arr,
        coverage=coverage,
        calibration=calibration,
        sharpness=sharpness,
        n_valid=n_valid,
        shape=true.shape,
    )
