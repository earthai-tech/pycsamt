# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Boundary-condition losses on predicted resistivity grids.

Boundary constraints anchor prediction cells that the training data
has no sensitivity to, e.g. air cells above topography or the outer
mesh padding, to an explicit required value instead of letting the
network hallucinate structure there.  A boundary constraint is
therefore a masked data-fit loss between the prediction and a
required ``target``, restricted to a caller-supplied
``boundary_mask``.  There is no implicit default target: callers must
state the physically motivated value explicitly, e.g. a fixed air
resistivity, matching the plan's requirement that agents contain no
hidden physics.

A common ``boundary_mask`` source is
:meth:`~pycsamt.ai.geology.topography.TopographicSurface.air_mask`.

All functions operate on plain NumPy arrays so the module stays
importable without an optional deep-learning backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .model import ModelLossResult, _model_loss

__all__ = [
    "boundary_condition_loss",
    "BoundaryLoss",
]

_KINDS = ("l1", "l2", "huber")
_REDUCTIONS = ("mean", "sum")


def _boundary_target(target: Any, shape: tuple[int, ...]) -> np.ndarray:
    """Broadcast a scalar or array boundary target to *shape*."""
    try:
        array = np.broadcast_to(
            np.asarray(target, dtype=float), shape
        ).astype(float)
    except ValueError as error:
        raise ValueError(
            f"target must be broadcastable to shape {shape}."
        ) from error
    if not np.all(np.isfinite(array)):
        raise ValueError("target must be finite.")
    return array


def boundary_condition_loss(
    y_pred: Any,
    *,
    boundary_mask: Any,
    target: Any,
    kind: str = "l2",
    delta: float = 1.0,
    valid: Any | None = None,
    weights: Any | None = None,
    reduction: str = "mean",
) -> ModelLossResult:
    """Penalize predicted values that violate a boundary constraint.

    Parameters
    ----------
    y_pred : array-like
        Predicted model values on a canonical geological grid.
    boundary_mask : array-like of bool, same shape as ``y_pred``
        Cells subject to the boundary constraint, e.g. air cells
        above topography or the outer mesh padding. At least one
        cell must be selected.
    target : float or array-like
        Required value on ``boundary_mask`` cells, e.g. a fixed air
        resistivity. A scalar is broadcast to the grid shape.
    kind : {"l1", "l2", "huber"}, default="l2"
        Elementwise penalty, as in
        :func:`~pycsamt.ai.losses.model.model_l2_loss`.
    delta : float, default=1.0
        Huber transition point. Ignored unless ``kind="huber"``.
    valid : array-like of bool or None, optional
        Additional cell mask combined with ``boundary_mask`` and
        with finite-value masking of ``y_pred``.
    weights : array-like or None, optional
        Non-negative per-cell weights broadcastable to the grid
        shape.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over included boundary cells.

    Returns
    -------
    ModelLossResult
        Reduced boundary-condition penalty.

    Examples
    --------
    >>> import numpy as np
    >>> grid = np.array([[1.0, 1.0], [3.0, 3.0]])
    >>> air = np.array([[True, True], [False, False]])
    >>> boundary_condition_loss(
    ...     grid, boundary_mask=air, target=0.0, kind="l1"
    ... ).value
    1.0
    """
    array = np.asarray(y_pred, dtype=float)
    mask = np.asarray(boundary_mask, dtype=bool)
    if mask.shape != array.shape:
        raise ValueError(
            "boundary_mask must have the same shape as y_pred."
        )
    if not np.any(mask):
        raise ValueError("boundary_mask must select at least one cell.")
    if valid is not None:
        supplied = np.asarray(valid, dtype=bool)
        if supplied.shape != array.shape:
            raise ValueError("valid must have the same shape as y_pred.")
        mask = mask & supplied
    target_array = _boundary_target(target, array.shape)
    return _model_loss(
        array,
        target_array,
        kind=kind,
        delta=delta,
        valid=mask,
        weights=weights,
        reduction=reduction,
    )


@dataclass(frozen=True)
class BoundaryLoss:
    """Configurable, callable boundary-condition penalty.

    Parameters
    ----------
    kind : {"l1", "l2", "huber"}, default="l2"
        Elementwise penalty applied to each boundary-cell residual.
    delta : float, default=1.0
        Huber transition point. Ignored unless ``kind="huber"``.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over included boundary cells.

    Examples
    --------
    >>> import numpy as np
    >>> loss = BoundaryLoss(kind="l1")
    >>> grid = np.array([[1.0, 1.0], [3.0, 3.0]])
    >>> air = np.array([[True, True], [False, False]])
    >>> loss(grid, boundary_mask=air, target=0.0).value
    1.0
    """

    kind: str = "l2"
    delta: float = 1.0
    reduction: str = "mean"

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise ValueError(f"kind must be one of {_KINDS}.")
        if self.reduction not in _REDUCTIONS:
            raise ValueError(f"reduction must be one of {_REDUCTIONS}.")

    def __call__(
        self,
        y_pred: Any,
        *,
        boundary_mask: Any,
        target: Any,
        valid: Any | None = None,
        weights: Any | None = None,
    ) -> ModelLossResult:
        """Evaluate the configured boundary-condition penalty."""
        return boundary_condition_loss(
            y_pred,
            boundary_mask=boundary_mask,
            target=target,
            kind=self.kind,
            delta=self.delta,
            valid=valid,
            weights=weights,
            reduction=self.reduction,
        )
