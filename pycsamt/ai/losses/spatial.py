# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Spatial regularization losses on predicted resistivity grids.

These losses implement the ``L_grad_x``, ``L_grad_z``, and ``L_TV``
terms of the staged inversion objective in the AI-inversion plan::

    L = w_m * L_model + lambda_x * L_grad_x + lambda_z * L_grad_z
        + lambda_tv * L_TV + lambda_d * L_response

Inputs are predicted model arrays on a canonical geological grid,
e.g. ``(z, x)`` or ``(z, y, x)`` from
:class:`~pycsamt.ai.geology.GeologyGrid`.  Only forward-difference
pairs where both endpoints are finite, user-valid, and non-zero
weight are included.

All functions operate on plain NumPy arrays so the module stays
importable without an optional deep-learning backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ._common import _reduce, _weight_array

__all__ = [
    "SpatialLossResult",
    "gradient_smoothness_loss",
    "total_variation_loss",
    "SpatialLoss",
]

_KINDS = ("l1", "l2")
_REDUCTIONS = ("mean", "sum")


@dataclass(frozen=True)
class SpatialLossResult:
    """Immutable scalar result of a spatial regularization loss.

    Parameters
    ----------
    value : float
        Reduced penalty value. ``nan`` when no difference was
        included and ``reduction="mean"``.
    kind : {"l1", "l2"}
        Elementwise penalty applied to each spatial difference.
    label : str
        Loss identity, e.g. ``"grad_axis0"`` or ``"tv"``.
    reduction : {"mean", "sum"}
        Reduction applied over valid, weighted differences.
    n_valid : int
        Number of differences included after masking.
    weight_sum : float
        Sum of weights over included differences.

    Examples
    --------
    >>> import numpy as np
    >>> grid = np.array([[0.0, 1.0, 3.0], [0.0, 0.0, 0.0]])
    >>> result = gradient_smoothness_loss(grid, axis=1)
    >>> result.label
    'grad_axis1'
    """

    value: float
    kind: str
    label: str
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


def _normalize_axis(ndim: int, axis: int) -> int:
    """Return *axis* normalized into ``[0, ndim)`` or raise."""
    normalized = axis + ndim if axis < 0 else axis
    if not (0 <= normalized < ndim):
        raise ValueError(f"axis must be in [-{ndim}, {ndim}).")
    return normalized


def _axis_pair(array: np.ndarray, axis: int) -> tuple[np.ndarray, np.ndarray]:
    """Return the front/back views shifted by one cell along *axis*."""
    front = [slice(None)] * array.ndim
    back = [slice(None)] * array.ndim
    front[axis] = slice(1, None)
    back[axis] = slice(None, -1)
    return array[tuple(front)], array[tuple(back)]


def _spatial_difference_loss(
    array: np.ndarray,
    *,
    axis: int,
    kind: str,
    valid: Any | None,
    weights: Any | None,
    reduction: str,
    label: str,
) -> SpatialLossResult:
    """Compute one masked, weighted first-difference penalty.

    Assumes *axis* is already a valid, normalized array axis.
    """
    if kind not in _KINDS:
        raise ValueError(f"kind must be one of {_KINDS}.")
    if array.shape[axis] < 2:
        raise ValueError(
            f"axis {axis} must have at least two cells to difference."
        )
    mask = np.isfinite(array)
    if valid is not None:
        supplied = np.asarray(valid, dtype=bool)
        if supplied.shape != array.shape:
            raise ValueError("valid must have the same shape as y_pred.")
        mask &= supplied
    front_values, back_values = _axis_pair(array, axis)
    front_mask, back_mask = _axis_pair(mask, axis)
    pair_mask = front_mask & back_mask
    difference = np.where(pair_mask, front_values - back_values, 0.0)
    per_cell = np.abs(difference) if kind == "l1" else np.square(difference)
    weight_array = _weight_array(weights, array.shape, "weights")
    if weight_array is None:
        active_weights = pair_mask.astype(float)
    else:
        front_weights, back_weights = _axis_pair(weight_array, axis)
        active_weights = np.where(
            pair_mask, np.minimum(front_weights, back_weights), 0.0
        )
    value, weight_sum = _reduce(per_cell, active_weights, reduction)
    return SpatialLossResult(
        value=value,
        kind=kind,
        label=label,
        reduction=reduction,
        n_valid=int(np.count_nonzero(pair_mask)),
        weight_sum=weight_sum,
    )


def gradient_smoothness_loss(
    y_pred: Any,
    *,
    axis: int,
    kind: str = "l2",
    valid: Any | None = None,
    weights: Any | None = None,
    reduction: str = "mean",
) -> SpatialLossResult:
    """Penalize first-difference magnitude along one grid axis.

    Parameters
    ----------
    y_pred : array-like
        Predicted model values on a canonical geological grid.
    axis : int
        Grid axis along which to difference, e.g. ``0`` for depth
        or ``-1`` for the horizontal direction. Negative axes are
        supported.
    kind : {"l1", "l2"}, default="l2"
        Elementwise penalty applied to each difference.
    valid : array-like of bool or None, optional
        Cell mask applied before differencing. A difference is kept
        only if both of its endpoint cells are valid and finite.
    weights : array-like or None, optional
        Non-negative per-cell weights broadcastable to ``y_pred``.
        A difference is weighted by the minimum of its two endpoint
        weights.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over included differences.

    Returns
    -------
    SpatialLossResult
        Reduced smoothness penalty, e.g. ``L_grad_x`` or ``L_grad_z``.

    Examples
    --------
    >>> import numpy as np
    >>> grid = np.array([[0.0, 1.0, 3.0], [0.0, 0.0, 0.0]])
    >>> gradient_smoothness_loss(grid, axis=1, kind="l1").value
    0.75
    """
    array = np.asarray(y_pred, dtype=float)
    if array.ndim == 0 or array.size == 0:
        raise ValueError("y_pred must be a non-empty array.")
    normalized = _normalize_axis(array.ndim, axis)
    return _spatial_difference_loss(
        array,
        axis=normalized,
        kind=kind,
        valid=valid,
        weights=weights,
        reduction=reduction,
        label=f"grad_axis{normalized}",
    )


def total_variation_loss(
    y_pred: Any,
    *,
    kind: str = "l1",
    valid: Any | None = None,
    weights: Any | None = None,
    reduction: str = "mean",
) -> SpatialLossResult:
    """Penalize anisotropic total variation over every spatial axis.

    Computes :func:`gradient_smoothness_loss` along each axis of
    ``y_pred`` and combines them, matching the standard anisotropic
    total-variation definition (a per-axis sum of directional
    gradients, as opposed to an isotropic pointwise gradient norm).

    Parameters
    ----------
    y_pred : array-like
        Predicted model values on a canonical geological grid.
    kind : {"l1", "l2"}, default="l1"
        Elementwise penalty applied to each difference.
    valid : array-like of bool or None, optional
        Cell mask shared by every axis.
    weights : array-like or None, optional
        Non-negative per-cell weights shared by every axis.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over all included differences from every
        axis combined.

    Returns
    -------
    SpatialLossResult
        Reduced total-variation penalty, ``L_TV``, labeled ``"tv"``.

    Examples
    --------
    >>> import numpy as np
    >>> grid = np.array([[0.0, 1.0], [0.0, 3.0]])
    >>> total_variation_loss(grid).value
    1.5
    """
    array = np.asarray(y_pred, dtype=float)
    if array.ndim == 0 or array.size == 0:
        raise ValueError("y_pred must be a non-empty array.")
    total = 0.0
    weight_total = 0.0
    n_valid = 0
    for axis in range(array.ndim):
        axis_result = _spatial_difference_loss(
            array,
            axis=axis,
            kind=kind,
            valid=valid,
            weights=weights,
            reduction="sum",
            label=f"grad_axis{axis}",
        )
        total += axis_result.value
        weight_total += axis_result.weight_sum
        n_valid += axis_result.n_valid
    if reduction == "sum":
        value = total
    elif reduction == "mean":
        value = total / weight_total if weight_total > 0 else float("nan")
    else:
        raise ValueError(f"reduction must be one of {_REDUCTIONS}.")
    return SpatialLossResult(
        value=value,
        kind=kind,
        label="tv",
        reduction=reduction,
        n_valid=n_valid,
        weight_sum=weight_total,
    )


@dataclass(frozen=True)
class SpatialLoss:
    """Configurable combination of gradient and TV regularizers.

    Combines the ``lambda_x * L_grad_x + lambda_z * L_grad_z +
    lambda_tv * L_TV`` terms of the staged inversion objective for a
    canonical 2-D ``(z, x)`` grid, where depth is axis 0 and the
    horizontal direction is axis 1.

    Parameters
    ----------
    lambda_x : float, default=1.0
        Weight applied to the horizontal-gradient term.
    lambda_z : float, default=1.0
        Weight applied to the depth-gradient term.
    lambda_tv : float, default=0.0
        Weight applied to the total-variation term.
    kind : {"l1", "l2"}, default="l2"
        Elementwise penalty shared by the two gradient terms. The
        total-variation term always uses ``"l1"``, matching its
        standard definition.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied within each enabled term.

    Examples
    --------
    >>> import numpy as np
    >>> grid = np.array([[0.0, 1.0], [0.0, 3.0]])
    >>> loss = SpatialLoss(
    ...     lambda_x=1.0, lambda_z=0.0, lambda_tv=0.0, kind="l1"
    ... )
    >>> loss(grid)
    2.0
    """

    lambda_x: float = 1.0
    lambda_z: float = 1.0
    lambda_tv: float = 0.0
    kind: str = "l2"
    reduction: str = "mean"

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise ValueError(f"kind must be one of {_KINDS}.")
        if self.reduction not in _REDUCTIONS:
            raise ValueError(f"reduction must be one of {_REDUCTIONS}.")
        for name in ("lambda_x", "lambda_z", "lambda_tv"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and non-negative.")
            object.__setattr__(self, name, value)

    def __call__(
        self,
        y_pred: Any,
        *,
        valid: Any | None = None,
        weights: Any | None = None,
    ) -> float:
        """Evaluate the weighted sum of the configured spatial terms.

        Parameters
        ----------
        y_pred : array-like, shape (n_z, n_x)
            Predicted model values on a canonical 2-D grid.
        valid : array-like of bool or None, optional
            Cell mask shared by every enabled term.
        weights : array-like or None, optional
            Non-negative per-cell weights shared by every enabled
            term.

        Returns
        -------
        float
            ``lambda_x * L_grad_x + lambda_z * L_grad_z +
            lambda_tv * L_TV``.
        """
        array = np.asarray(y_pred, dtype=float)
        if array.ndim != 2:
            raise ValueError("y_pred must be a 2-D (z, x) grid.")
        total = 0.0
        if self.lambda_z > 0:
            total += (
                self.lambda_z
                * gradient_smoothness_loss(
                    array,
                    axis=0,
                    kind=self.kind,
                    valid=valid,
                    weights=weights,
                    reduction=self.reduction,
                ).value
            )
        if self.lambda_x > 0:
            total += (
                self.lambda_x
                * gradient_smoothness_loss(
                    array,
                    axis=1,
                    kind=self.kind,
                    valid=valid,
                    weights=weights,
                    reduction=self.reduction,
                ).value
            )
        if self.lambda_tv > 0:
            total += (
                self.lambda_tv
                * total_variation_loss(
                    array,
                    valid=valid,
                    weights=weights,
                    reduction=self.reduction,
                ).value
            )
        return total
