# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Evaluation metrics for 1-D EM inversion networks.

All functions accept numpy arrays and ignore NaN entries (corresponding
to padding in variable-depth datasets).

Metrics
-------
``rmse``           Root mean square error on log₁₀(ρ).
``mae``            Mean absolute error.
``r2``             Coefficient of determination R².
``relative_rmse``  RMSE normalised by |y_true|.
``depth_rmse``     RMSE weighted by inverse depth (shallower layers
                   are easier to recover and receive less weight).
``layer_rmse``     Per-layer RMSE vector.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

__all__ = [
    "rmse",
    "mae",
    "r2",
    "relative_rmse",
    "depth_rmse",
    "layer_rmse",
    "masked_mse_loss",
    "summarise",
]


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _valid(y_true: np.ndarray, y_pred: np.ndarray):
    """Return (yt, yp) with NaN rows removed (any NaN in a row discards it)."""
    mask = np.isfinite(y_true).all(axis=-1) & np.isfinite(y_pred).all(axis=-1)
    return y_true[mask], y_pred[mask]


# ─────────────────────────────────────────────────────────────────────────────
# Numpy metrics
# ─────────────────────────────────────────────────────────────────────────────

def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """
    Root mean square error, ignoring NaN.

    Operates element-wise on ``y_true`` and ``y_pred``.
    Both arrays are assumed to be in log₁₀(Ω·m) or normalised space.
    """
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    if mask.sum() == 0:
        return float("nan")
    return float(np.sqrt(np.mean((y_true[mask] - y_pred[mask]) ** 2)))


def mae(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Mean absolute error, ignoring NaN."""
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    if mask.sum() == 0:
        return float("nan")
    return float(np.mean(np.abs(y_true[mask] - y_pred[mask])))


def r2(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Coefficient of determination R², ignoring NaN."""
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    if mask.sum() == 0:
        return float("nan")
    yt, yp = y_true[mask], y_pred[mask]
    ss_res = np.sum((yt - yp) ** 2)
    ss_tot = np.sum((yt - np.mean(yt)) ** 2)
    return float(1.0 - ss_res / (ss_tot + 1e-12))


def relative_rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """
    Normalised RMSE: ``sqrt(mean((y_true - y_pred)² / y_true²))``.

    Useful when comparing models with very different resistivity ranges.
    """
    mask = np.isfinite(y_true) & np.isfinite(y_pred) & (np.abs(y_true) > 1e-12)
    if mask.sum() == 0:
        return float("nan")
    yt, yp = y_true[mask], y_pred[mask]
    return float(np.sqrt(np.mean(((yt - yp) / yt) ** 2)))


def depth_rmse(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    n_layers: int,
    *,
    depth_weight: bool = True,
) -> float:
    """
    RMSE on the resistivity sub-vector only (first ``n_layers`` columns),
    optionally weighted by layer index so that deeper (harder) layers
    have less influence on the score.

    Parameters
    ----------
    y_true, y_pred : ndarray, shape (n_samples, n_params)
    n_layers : int
        Number of resistivity values in the parameter vector.
    depth_weight : bool
        If True, weight by ``1 / (1 + layer_index)``.
    """
    yt = y_true[:, :n_layers]
    yp = y_pred[:, :n_layers]
    mask = np.isfinite(yt) & np.isfinite(yp)
    if mask.sum() == 0:
        return float("nan")
    diff_sq = (yt - yp) ** 2
    if depth_weight:
        w = 1.0 / (1.0 + np.arange(n_layers, dtype=float))
        w /= w.sum()
        diff_sq = diff_sq * w[None, :]
    return float(np.sqrt(np.nanmean(diff_sq)))


def layer_rmse(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    """
    Per-column RMSE vector.

    Returns
    -------
    rmse_per_col : ndarray, shape (n_params,)
        RMSE for each parameter column independently.
    """
    n_cols = y_true.shape[1] if y_true.ndim > 1 else 1
    out = np.empty(n_cols)
    for j in range(n_cols):
        out[j] = rmse(y_true[:, j], y_pred[:, j])
    return out


def summarise(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    n_layers: Optional[int] = None,
) -> dict:
    """
    Compute all scalar metrics and return as a dict.

    Returns
    -------
    metrics : dict with keys
        ``'rmse'``, ``'mae'``, ``'r2'``, ``'relative_rmse'``,
        ``'depth_rmse'`` (if n_layers given).
    """
    d = {
        "rmse": rmse(y_true, y_pred),
        "mae": mae(y_true, y_pred),
        "r2": r2(y_true, y_pred),
        "relative_rmse": relative_rmse(y_true, y_pred),
    }
    if n_layers is not None:
        d["depth_rmse"] = depth_rmse(y_true, y_pred, n_layers)
    return d


# ─────────────────────────────────────────────────────────────────────────────
# PyTorch masked loss
# ─────────────────────────────────────────────────────────────────────────────

def masked_mse_loss(pred, target):
    """
    MSE loss that ignores ``NaN`` (padding) entries in *target*.

    Parameters
    ----------
    pred : torch.Tensor, shape (batch, n_out)
    target : torch.Tensor, shape (batch, n_out)

    Returns
    -------
    loss : torch.Tensor  (scalar)
    """
    mask = torch_isfinite(target)
    if not mask.any():
        return pred.sum() * 0.0   # differentiable zero
    diff = (pred[mask] - target[mask]) ** 2
    return diff.mean()


def torch_isfinite(t):
    """Return a bool mask without importing torch at module level."""
    try:
        import torch
        return torch.isfinite(t)
    except ImportError:
        raise ImportError("PyTorch required for masked_mse_loss")
