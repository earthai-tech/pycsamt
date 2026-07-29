# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Synthetic-recovery diagnostics for known-truth geological grids.

These diagnostics implement the M0 baseline metrics named in the
AI-inversion plan: log-resistivity MAE/RMSE, structural similarity,
and recovery broken down by depth.  Inputs are predicted and true
model arrays sharing one canonical geological-grid shape, e.g.
``(z, x)`` or ``(z, y, x)`` from
:class:`~pycsamt.ai.geology.GeologyGrid`.  Recovery diagnostics are
only meaningful when the true model is known, i.e. on synthetic
data; field-survey validation instead relies on
:mod:`~pycsamt.ai.losses.response` residuals.

All functions operate on plain NumPy arrays so the module stays
importable without an optional deep-learning backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.ndimage import uniform_filter

from ..losses.model import model_l1_loss, model_l2_loss

__all__ = [
    "RecoveryReport",
    "recovery_report",
    "structural_similarity",
    "depth_profile_rmse",
    "depth_profile_mae",
]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    """Return a read-only copy of *value* as an ndarray."""
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _validate_pair(y_pred: Any, y_true: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return y_pred/y_true as float arrays of one 2-D/3-D shape."""
    pred = np.asarray(y_pred, dtype=float)
    true = np.asarray(y_true, dtype=float)
    if pred.shape != true.shape:
        raise ValueError("y_pred and y_true must share one shape.")
    if pred.ndim not in (2, 3):
        raise ValueError("y_pred and y_true must be a 2-D or 3-D grid.")
    if pred.size == 0:
        raise ValueError("y_pred and y_true must not be empty.")
    return pred, true


def _validate_valid(
    valid: Any | None, shape: tuple[int, ...]
) -> np.ndarray | None:
    """Return an optional boolean mask validated against *shape*."""
    if valid is None:
        return None
    array = np.asarray(valid, dtype=bool)
    if array.shape != shape:
        raise ValueError("valid must have the same shape as y_pred.")
    return array


def _depth_profile(
    pred: np.ndarray,
    true: np.ndarray,
    *,
    axis: int,
    mask: np.ndarray,
    kind: str,
) -> np.ndarray:
    """Return per-layer masked RMSE (``kind="l2"``) or MAE along axis."""
    n_layers = pred.shape[axis]
    out = np.full(n_layers, np.nan)
    for layer in range(n_layers):
        index = [slice(None)] * pred.ndim
        index[axis] = layer
        index = tuple(index)
        layer_mask = mask[index]
        if not np.any(layer_mask):
            continue
        difference = pred[index][layer_mask] - true[index][layer_mask]
        if kind == "l1":
            out[layer] = float(np.mean(np.abs(difference)))
        else:
            out[layer] = float(np.sqrt(np.mean(np.square(difference))))
    return out


@dataclass(frozen=True)
class RecoveryReport:
    """Immutable synthetic-recovery diagnostics for one grid pair.

    Parameters
    ----------
    rmse, mae : float
        Global masked root-mean-square and mean-absolute error.
    r2 : float
        Coefficient of determination. ``nan`` when the true values
        are numerically constant, making R² undefined.
    ssim : float or None
        Structural similarity index (Wang et al., 2004), or ``None``
        when it was not requested or not computable, e.g. the grid
        is partially masked or smaller than the requested window.
    depth_rmse, depth_mae : ndarray
        Per-layer RMSE/MAE along the requested depth axis.
    n_valid : int
        Number of cells included after masking.
    shape : tuple of int
        Shape of the compared grids.

    Examples
    --------
    >>> import numpy as np
    >>> pred = np.array([[1.0, 2.0], [3.0, 4.0]])
    >>> true = np.array([[1.0, 2.0], [3.0, 6.0]])
    >>> report = recovery_report(pred, true, compute_ssim=False)
    >>> report.n_valid, report.shape
    (4, (2, 2))
    """

    rmse: float
    mae: float
    r2: float
    ssim: float | None
    depth_rmse: np.ndarray
    depth_mae: np.ndarray
    n_valid: int
    shape: tuple[int, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "rmse", float(self.rmse))
        object.__setattr__(self, "mae", float(self.mae))
        object.__setattr__(self, "r2", float(self.r2))
        object.__setattr__(
            self, "ssim", None if self.ssim is None else float(self.ssim)
        )
        object.__setattr__(self, "n_valid", int(self.n_valid))
        object.__setattr__(
            self, "shape", tuple(int(size) for size in self.shape)
        )
        object.__setattr__(
            self, "depth_rmse", _readonly(self.depth_rmse, float)
        )
        object.__setattr__(
            self, "depth_mae", _readonly(self.depth_mae, float)
        )


def structural_similarity(
    y_pred: Any,
    y_true: Any,
    *,
    window: int = 7,
    data_range: float | None = None,
) -> float:
    """Compute the mean structural similarity index (SSIM).

    Uses the windowed luminance/contrast/structure formulation of
    Wang et al. (2004) with a uniform (box) window, evaluated on the
    interior of the grid to avoid boundary-filter artifacts.

    Parameters
    ----------
    y_pred, y_true : array-like
        Fully finite 2-D or 3-D grids sharing one shape. SSIM has no
        defined masking rule, so both must be complete.
    window : int, default=7
        Positive odd window size, no larger than the smallest grid
        axis.
    data_range : float or None, optional
        Dynamic range of the compared values. Defaults to the range
        spanned by the combined ``y_pred``/``y_true`` values.

    Returns
    -------
    float
        Mean SSIM over the interior window positions, at most 1.0
        for identical grids.

    Examples
    --------
    >>> import numpy as np
    >>> grid = np.arange(64, dtype=float).reshape(8, 8)
    >>> structural_similarity(grid, grid, window=3)
    1.0
    """
    pred, true = _validate_pair(y_pred, y_true)
    if not np.all(np.isfinite(pred)) or not np.all(np.isfinite(true)):
        raise ValueError(
            "structural_similarity requires fully finite grids; mask "
            "or exclude cells before calling it."
        )
    if (
        not isinstance(window, int)
        or isinstance(window, bool)
        or window < 1
        or window % 2 == 0
    ):
        raise ValueError("window must be a positive odd integer.")
    if window > min(pred.shape):
        raise ValueError(
            "window must not exceed the smallest grid axis "
            f"{min(pred.shape)}."
        )
    span = data_range
    if span is None:
        span = float(max(pred.max(), true.max()) - min(pred.min(), true.min()))
    if not np.isfinite(span) or span <= 0:
        raise ValueError("data_range must be finite and positive.")
    mu_x = uniform_filter(pred, size=window)
    mu_y = uniform_filter(true, size=window)
    sigma_x2 = uniform_filter(pred * pred, size=window) - mu_x * mu_x
    sigma_y2 = uniform_filter(true * true, size=window) - mu_y * mu_y
    sigma_xy = uniform_filter(pred * true, size=window) - mu_x * mu_y
    c1 = (0.01 * span) ** 2
    c2 = (0.03 * span) ** 2
    numerator = (2 * mu_x * mu_y + c1) * (2 * sigma_xy + c2)
    denominator = (mu_x**2 + mu_y**2 + c1) * (sigma_x2 + sigma_y2 + c2)
    ssim_map = numerator / denominator
    crop = window // 2
    interior = tuple(
        slice(crop, size - crop) for size in ssim_map.shape
    )
    cropped = ssim_map[interior]
    if cropped.size == 0:
        raise ValueError("window is too large for this grid shape.")
    return float(np.mean(cropped))


def depth_profile_rmse(
    y_pred: Any,
    y_true: Any,
    *,
    axis: int = 0,
    valid: Any | None = None,
) -> np.ndarray:
    """Return per-layer masked RMSE along one grid axis.

    Parameters
    ----------
    y_pred, y_true : array-like
        Predicted and true model values sharing one 2-D or 3-D grid
        shape.
    axis : int, default=0
        Grid axis to break down by, typically depth (``z``).
    valid : array-like of bool or None, optional
        Explicit cell mask, combined with finite-value masking of
        both inputs.

    Returns
    -------
    ndarray, shape (y_pred.shape[axis],)
        RMSE for each layer; ``nan`` for a layer with no valid
        cells.

    Examples
    --------
    >>> import numpy as np
    >>> pred = np.array([[1.0, 2.0], [3.0, 4.0]])
    >>> true = np.array([[1.0, 2.0], [3.0, 6.0]])
    >>> depth_profile_rmse(pred, true).tolist()
    [0.0, 1.4142135623730951]
    """
    pred, true = _validate_pair(y_pred, y_true)
    mask = np.isfinite(pred) & np.isfinite(true)
    supplied = _validate_valid(valid, pred.shape)
    if supplied is not None:
        mask &= supplied
    if axis < 0:
        axis += pred.ndim
    if not (0 <= axis < pred.ndim):
        raise ValueError(f"axis must be in [-{pred.ndim}, {pred.ndim}).")
    return _depth_profile(pred, true, axis=axis, mask=mask, kind="l2")


def depth_profile_mae(
    y_pred: Any,
    y_true: Any,
    *,
    axis: int = 0,
    valid: Any | None = None,
) -> np.ndarray:
    """Return per-layer masked MAE along one grid axis.

    Parameters
    ----------
    y_pred, y_true : array-like
        Predicted and true model values sharing one 2-D or 3-D grid
        shape.
    axis : int, default=0
        Grid axis to break down by, typically depth (``z``).
    valid : array-like of bool or None, optional
        Explicit cell mask, combined with finite-value masking of
        both inputs.

    Returns
    -------
    ndarray, shape (y_pred.shape[axis],)
        MAE for each layer; ``nan`` for a layer with no valid cells.

    Examples
    --------
    >>> import numpy as np
    >>> pred = np.array([[1.0, 2.0], [3.0, 4.0]])
    >>> true = np.array([[1.0, 2.0], [3.0, 6.0]])
    >>> depth_profile_mae(pred, true).tolist()
    [0.0, 1.0]
    """
    pred, true = _validate_pair(y_pred, y_true)
    mask = np.isfinite(pred) & np.isfinite(true)
    supplied = _validate_valid(valid, pred.shape)
    if supplied is not None:
        mask &= supplied
    if axis < 0:
        axis += pred.ndim
    if not (0 <= axis < pred.ndim):
        raise ValueError(f"axis must be in [-{pred.ndim}, {pred.ndim}).")
    return _depth_profile(pred, true, axis=axis, mask=mask, kind="l1")


def recovery_report(
    y_pred: Any,
    y_true: Any,
    *,
    valid: Any | None = None,
    depth_axis: int = 0,
    compute_ssim: bool = True,
    ssim_window: int = 7,
) -> RecoveryReport:
    """Summarize synthetic-recovery quality for one grid pair.

    Parameters
    ----------
    y_pred, y_true : array-like
        Predicted and true model values sharing one 2-D or 3-D grid
        shape.
    valid : array-like of bool or None, optional
        Explicit cell mask, combined with finite-value masking of
        both inputs.
    depth_axis : int, default=0
        Grid axis passed to :func:`depth_profile_rmse` and
        :func:`depth_profile_mae`.
    compute_ssim : bool, default=True
        Attempt :func:`structural_similarity`. Skipped (``ssim`` is
        ``None``) when the grid is partially masked or smaller than
        ``ssim_window``, since SSIM has no defined masking rule.
    ssim_window : int, default=7
        Window forwarded to :func:`structural_similarity`.

    Returns
    -------
    RecoveryReport
        Combined recovery diagnostics.

    Raises
    ------
    ValueError
        If shapes mismatch, inputs are empty, or no cell is valid.

    Examples
    --------
    >>> import numpy as np
    >>> pred = np.array([[1.0, 2.0], [3.0, 4.0]])
    >>> true = np.array([[1.0, 2.0], [3.0, 6.0]])
    >>> report = recovery_report(pred, true, compute_ssim=False)
    >>> round(report.rmse, 6), round(report.mae, 6)
    (1.0, 0.5)
    >>> round(report.r2, 6)
    0.714286
    """
    pred, true = _validate_pair(y_pred, y_true)
    supplied = _validate_valid(valid, pred.shape)
    mask = np.isfinite(pred) & np.isfinite(true)
    if supplied is not None:
        mask &= supplied
    if not np.any(mask):
        raise ValueError("no valid cell to compare.")
    rmse = float(np.sqrt(model_l2_loss(pred, true, valid=mask).value))
    mae = float(model_l1_loss(pred, true, valid=mask).value)
    yt = true[mask]
    yp = pred[mask]
    ss_res = float(np.sum(np.square(yt - yp)))
    ss_tot = float(np.sum(np.square(yt - np.mean(yt))))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    ssim = None
    if compute_ssim and np.all(mask) and min(pred.shape) >= ssim_window:
        ssim = structural_similarity(pred, true, window=ssim_window)
    return RecoveryReport(
        rmse=rmse,
        mae=mae,
        r2=r2,
        ssim=ssim,
        depth_rmse=_depth_profile(
            pred, true, axis=depth_axis, mask=mask, kind="l2"
        ),
        depth_mae=_depth_profile(
            pred, true, axis=depth_axis, mask=mask, kind="l1"
        ),
        n_valid=int(np.count_nonzero(mask)),
        shape=pred.shape,
    )
