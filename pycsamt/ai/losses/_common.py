# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Private array helpers shared by :mod:`pycsamt.ai.losses`
submodules.

Every helper operates on plain NumPy arrays so the package stays
importable without an optional deep-learning backend.
"""

from __future__ import annotations

from typing import Any

import numpy as np

__all__: list[str] = []


def _paired_arrays(
    y_pred: Any,
    y_true: Any,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``y_pred``/``y_true`` as float arrays of one shape."""
    pred = np.asarray(y_pred, dtype=float)
    true = np.asarray(y_true, dtype=float)
    if pred.shape != true.shape:
        raise ValueError("y_pred and y_true must share one shape.")
    if pred.size == 0:
        raise ValueError("y_pred and y_true must not be empty.")
    return pred, true


def _validity_mask(
    pred: np.ndarray,
    true: np.ndarray,
    valid: Any | None,
) -> np.ndarray:
    """Combine finiteness of both arrays with an optional mask."""
    mask = np.isfinite(pred) & np.isfinite(true)
    if valid is not None:
        supplied = np.asarray(valid, dtype=bool)
        if supplied.shape != mask.shape:
            raise ValueError(
                "valid must have the same shape as y_pred/y_true."
            )
        mask &= supplied
    return mask


def _weight_array(
    weights: Any,
    shape: tuple[int, ...],
    name: str,
) -> np.ndarray | None:
    """Validate optional non-negative *weights* broadcast to *shape*."""
    if weights is None:
        return None
    try:
        array = np.broadcast_to(
            np.asarray(weights, dtype=float), shape
        ).astype(float)
    except ValueError as error:
        raise ValueError(
            f"{name} must be broadcastable to shape {shape}."
        ) from error
    if not np.all(np.isfinite(array)) or np.any(array < 0):
        raise ValueError(f"{name} must be finite and non-negative.")
    return array


def _reduce(
    per_cell: np.ndarray,
    active_weights: np.ndarray,
    reduction: str,
) -> tuple[float, float]:
    """Reduce weighted per-cell values to ``(value, weight_sum)``."""
    if reduction not in {"mean", "sum"}:
        raise ValueError("reduction must be 'mean' or 'sum'.")
    weight_sum = float(np.sum(active_weights))
    total = float(np.sum(per_cell * active_weights))
    if reduction == "sum":
        return total, weight_sum
    if weight_sum <= 0:
        return float("nan"), weight_sum
    return total / weight_sum, weight_sum
