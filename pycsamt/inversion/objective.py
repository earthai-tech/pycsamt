# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Objective-function helpers for inversion backends."""

from __future__ import annotations

import numpy as np

__all__ = ["weighted_rms", "relative_errors"]


def relative_errors(values, floor: float = 0.05) -> np.ndarray:
    """Return absolute errors from a relative floor."""
    arr = np.asarray(values, dtype=float)
    return np.maximum(np.abs(arr) * float(floor), 1e-12)


def weighted_rms(observed, predicted, errors=None) -> float:
    """Return normalized weighted RMS misfit."""
    obs = np.asarray(observed, dtype=float)
    pred = np.asarray(predicted, dtype=float)
    if errors is None:
        err = np.ones_like(obs, dtype=float)
    else:
        err = np.asarray(errors, dtype=float)
    mask = np.isfinite(obs) & np.isfinite(pred) & np.isfinite(err) & (err > 0)
    if not np.any(mask):
        return float("nan")
    resid = (pred[mask] - obs[mask]) / err[mask]
    return float(np.sqrt(np.mean(resid ** 2)))
