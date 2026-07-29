# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Masked, weighted data-fit losses on canonical resistivity grids.

These losses implement the ``L_model`` term of the staged inversion
objective in the AI-inversion plan::

    L = w_m * L_model + lambda_x * L_grad_x + lambda_z * L_grad_z
        + lambda_tv * L_TV + lambda_d * L_response

Inputs are predicted and true model arrays that share one canonical
geological-grid shape, e.g. ``(z, x)`` or ``(z, y, x)`` from
:class:`~pycsamt.ai.geology.GeologyGrid`.  Cells may be excluded with
an explicit boolean mask, a non-finite value in either array, or a
zero weight.

All functions operate on plain NumPy arrays so the module stays
importable without an optional deep-learning backend.  Differentiable
use inside a training loop is left to the caller's autograd tensors,
which support the same elementwise arithmetic used here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ._common import _paired_arrays, _reduce, _validity_mask, _weight_array

__all__ = [
    "ModelLossResult",
    "ModelLoss",
    "model_l1_loss",
    "model_l2_loss",
    "model_huber_loss",
    "depth_weights",
]

_KINDS = ("l1", "l2", "huber")
_REDUCTIONS = ("mean", "sum")


@dataclass(frozen=True)
class ModelLossResult:
    """Immutable scalar result of a model data-fit loss.

    Parameters
    ----------
    value : float
        Reduced loss value. ``nan`` when no cell was included and
        ``reduction="mean"``.
    kind : {"l1", "l2", "huber"}
        Elementwise loss family that produced ``value``.
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
    >>> result = model_l2_loss(
    ...     np.array([1.0, 2.0]), np.array([1.0, 0.0])
    ... )
    >>> result.value, result.n_valid
    (2.0, 2)
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


def _model_loss(
    y_pred: Any,
    y_true: Any,
    *,
    kind: str,
    delta: float = 1.0,
    valid: Any | None = None,
    weights: Any | None = None,
    reduction: str = "mean",
) -> ModelLossResult:
    """Compute one masked, weighted elementwise model loss."""
    if kind not in _KINDS:
        raise ValueError(f"kind must be one of {_KINDS}.")
    pred, true = _paired_arrays(y_pred, y_true)
    mask = _validity_mask(pred, true, valid)
    cell_weights = _weight_array(weights, pred.shape, "weights")
    difference = np.where(mask, pred - true, 0.0)
    if kind == "l1":
        per_cell = np.abs(difference)
    elif kind == "l2":
        per_cell = np.square(difference)
    else:
        if not np.isfinite(delta) or delta <= 0:
            raise ValueError("delta must be finite and positive.")
        absolute = np.abs(difference)
        quadratic = np.minimum(absolute, delta)
        linear = absolute - quadratic
        per_cell = 0.5 * np.square(quadratic) + delta * linear
    active_weights = (
        mask.astype(float) if cell_weights is None else cell_weights * mask
    )
    value, weight_sum = _reduce(per_cell, active_weights, reduction)
    return ModelLossResult(
        value=value,
        kind=kind,
        reduction=reduction,
        n_valid=int(np.count_nonzero(mask)),
        weight_sum=weight_sum,
    )


def model_l1_loss(
    y_pred: Any,
    y_true: Any,
    *,
    valid: Any | None = None,
    weights: Any | None = None,
    reduction: str = "mean",
) -> ModelLossResult:
    """Masked, weighted mean/summed absolute error.

    Parameters
    ----------
    y_pred, y_true : array-like
        Predicted and true model values sharing one shape, typically
        log-resistivity on a canonical geological grid.
    valid : array-like of bool or None, optional
        Explicit cell mask, combined with finite-value masking of
        both inputs.
    weights : array-like or None, optional
        Non-negative per-cell weights broadcastable to the shared
        shape, e.g. from :func:`depth_weights`.
    reduction : {"mean", "sum"}, default="mean"
        Whether to divide by the total weight or return the raw sum.

    Returns
    -------
    ModelLossResult
        Reduced ``L1`` loss with provenance.

    Examples
    --------
    >>> import numpy as np
    >>> model_l1_loss(
    ...     np.array([1.0, 3.0]), np.array([1.0, 1.0])
    ... ).value
    1.0
    """
    return _model_loss(
        y_pred,
        y_true,
        kind="l1",
        valid=valid,
        weights=weights,
        reduction=reduction,
    )


def model_l2_loss(
    y_pred: Any,
    y_true: Any,
    *,
    valid: Any | None = None,
    weights: Any | None = None,
    reduction: str = "mean",
) -> ModelLossResult:
    """Masked, weighted mean/summed squared error.

    Parameters
    ----------
    y_pred, y_true : array-like
        Predicted and true model values sharing one shape.
    valid : array-like of bool or None, optional
        Explicit cell mask, combined with finite-value masking of
        both inputs.
    weights : array-like or None, optional
        Non-negative per-cell weights broadcastable to the shared
        shape.
    reduction : {"mean", "sum"}, default="mean"
        Whether to divide by the total weight or return the raw sum.

    Returns
    -------
    ModelLossResult
        Reduced ``L2`` loss with provenance.

    Examples
    --------
    >>> import numpy as np
    >>> model_l2_loss(
    ...     np.array([1.0, 3.0]), np.array([1.0, 1.0])
    ... ).value
    2.0
    """
    return _model_loss(
        y_pred,
        y_true,
        kind="l2",
        valid=valid,
        weights=weights,
        reduction=reduction,
    )


def model_huber_loss(
    y_pred: Any,
    y_true: Any,
    *,
    delta: float = 1.0,
    valid: Any | None = None,
    weights: Any | None = None,
    reduction: str = "mean",
) -> ModelLossResult:
    """Masked, weighted Huber loss, robust to outlier cells.

    Parameters
    ----------
    y_pred, y_true : array-like
        Predicted and true model values sharing one shape.
    delta : float, default=1.0
        Positive transition point between the quadratic and linear
        regimes.
    valid : array-like of bool or None, optional
        Explicit cell mask, combined with finite-value masking of
        both inputs.
    weights : array-like or None, optional
        Non-negative per-cell weights broadcastable to the shared
        shape.
    reduction : {"mean", "sum"}, default="mean"
        Whether to divide by the total weight or return the raw sum.

    Returns
    -------
    ModelLossResult
        Reduced Huber loss with provenance.

    Examples
    --------
    >>> import numpy as np
    >>> small = model_huber_loss(
    ...     np.array([0.5]), np.array([0.0]), delta=1.0
    ... ).value
    >>> large = model_huber_loss(
    ...     np.array([5.0]), np.array([0.0]), delta=1.0
    ... ).value
    >>> round(small, 3), round(large, 3)
    (0.125, 4.5)
    """
    return _model_loss(
        y_pred,
        y_true,
        kind="huber",
        delta=delta,
        valid=valid,
        weights=weights,
        reduction=reduction,
    )


def depth_weights(n_depth: int) -> np.ndarray:
    """Return inverse-depth weights normalized to sum to one.

    Shallower cells (small index) receive more weight than deeper
    ones, matching the intuition that shallow structure is easier to
    recover and should not dominate a training loss.

    Parameters
    ----------
    n_depth : int
        Number of depth cells, at least one.

    Returns
    -------
    ndarray, shape (n_depth,)
        Weights proportional to ``1 / (1 + depth_index)``.

    Examples
    --------
    >>> import numpy as np
    >>> weights = depth_weights(2)
    >>> np.round(weights, 6)
    array([0.666667, 0.333333])
    """
    if (
        not isinstance(n_depth, int)
        or isinstance(n_depth, bool)
        or n_depth < 1
    ):
        raise ValueError("n_depth must be a positive integer.")
    raw = 1.0 / (1.0 + np.arange(n_depth, dtype=float))
    return raw / raw.sum()


@dataclass(frozen=True)
class ModelLoss:
    """Configurable, callable masked model data-fit loss.

    Bundles a loss family, reduction, and optional fixed per-cell
    weights so the same configuration can be reused across batches in
    a training loop.

    Parameters
    ----------
    kind : {"l1", "l2", "huber"}, default="l2"
        Elementwise loss family.
    delta : float, default=1.0
        Huber transition point. Ignored unless ``kind="huber"``.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over valid, weighted cells.
    weights : ndarray or None, optional
        Fixed per-cell weights reused on every call, broadcastable to
        the input shape. A ``weights`` argument passed directly to
        :meth:`__call__` overrides this default for that call only.

    Examples
    --------
    >>> import numpy as np
    >>> loss = ModelLoss(kind="l1")
    >>> loss(np.array([1.0, 3.0]), np.array([1.0, 1.0])).value
    1.0
    """

    kind: str = "l2"
    delta: float = 1.0
    reduction: str = "mean"
    weights: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise ValueError(f"kind must be one of {_KINDS}.")
        if self.reduction not in _REDUCTIONS:
            raise ValueError(f"reduction must be one of {_REDUCTIONS}.")
        if self.weights is not None:
            weights = np.asarray(self.weights, dtype=float)
            if not np.all(np.isfinite(weights)) or np.any(weights < 0):
                raise ValueError(
                    "weights must be finite and non-negative."
                )
            object.__setattr__(self, "weights", weights)

    def __call__(
        self,
        y_pred: Any,
        y_true: Any,
        *,
        valid: Any | None = None,
        weights: Any | None = None,
    ) -> ModelLossResult:
        """Evaluate the configured loss on one pair of model grids."""
        return _model_loss(
            y_pred,
            y_true,
            kind=self.kind,
            delta=self.delta,
            valid=valid,
            weights=self.weights if weights is None else weights,
            reduction=self.reduction,
        )

    @classmethod
    def with_depth_weights(
        cls,
        n_depth: int,
        *,
        dimension: int = 2,
        kind: str = "l2",
        delta: float = 1.0,
        reduction: str = "mean",
    ) -> ModelLoss:
        """Build a loss weighted by inverse depth on the leading axis.

        Parameters
        ----------
        n_depth : int
            Number of depth cells along the grid's leading axis.
        dimension : {2, 3}, default=2
            Grid rank, matching
            :attr:`~pycsamt.ai.geology.GeologyGrid.dimension`: 2 for
            ``(z, x)`` grids, 3 for ``(z, y, x)`` grids. Used only to
            reshape the depth weights for broadcasting.
        kind, delta, reduction
            Forwarded to the constructor.

        Returns
        -------
        ModelLoss
            Loss whose stored weights broadcast against grids shaped
            ``(n_depth, ...)``.

        Examples
        --------
        >>> loss = ModelLoss.with_depth_weights(3, dimension=2)
        >>> loss.weights.shape
        (3, 1)
        """
        if dimension not in (2, 3):
            raise ValueError("dimension must be 2 or 3.")
        shape = (n_depth,) + (1,) * (dimension - 1)
        return cls(
            kind=kind,
            delta=delta,
            reduction=reduction,
            weights=depth_weights(n_depth).reshape(shape),
        )
