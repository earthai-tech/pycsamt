# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Observation-reliability factors for two-dimensional inversion.

This module keeps measurement quality separate from compatibility with a
2-D physical model. For phase-tensor skew angle :math:`\beta` and a frozen
scale :math:`\beta_0`, dimensionality reliability is

.. math::

    c_{2D} = \exp[-(|\beta| / \beta_0)^2].

The final observation reliability is the bounded product of measurement and
dimensionality factors. Keeping both arrays allows explicit ablation of the
two effects.
"""

from __future__ import annotations

import numpy as np

from ...compat.sklearn import validate_params
from .schema import (
    DIMENSIONALITY_RELIABILITY_SCHEMA,
    RELIABILITY_COMBINATION_SCHEMA,
)

__all__ = [
    "combine_observation_reliability",
    "dimensionality_reliability",
]


@validate_params(DIMENSIONALITY_RELIABILITY_SCHEMA)
def dimensionality_reliability(
    beta_deg: np.ndarray,
    *,
    beta_scale_deg: float = 5.0,
    minimum: float = 0.0,
) -> np.ndarray:
    r"""Convert phase-tensor skew into compatibility with 2-D physics.

    Parameters
    ----------
    beta_deg : array-like of float
        Phase-tensor skew angles in degrees. Non-finite entries receive the
        configured minimum reliability.
    beta_scale_deg : float, default=5.0
        Positive skew scale :math:`\beta_0`. Reliability equals ``exp(-1)``
        at this absolute skew.
    minimum : float, default=0.0
        Lower bound in the closed interval ``[0, 1]``.

    Returns
    -------
    numpy.ndarray of float
        Reliability values with the same shape as ``beta_deg``.

    Raises
    ------
    ValueError
        If the input is empty or scalar controls are invalid.

    Examples
    --------
    >>> dimensionality_reliability([0.0, 5.0]).round(6).tolist()
    [1.0, 0.367879]
    """
    beta = np.asarray(beta_deg, dtype=float)
    scale = float(beta_scale_deg)
    floor = float(minimum)
    if beta.size == 0:
        raise ValueError("beta_deg must contain at least one value")
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("beta_scale_deg must be finite and positive")
    if not np.isfinite(floor) or not 0 <= floor <= 1:
        raise ValueError("minimum must lie in [0, 1]")
    result = np.full(beta.shape, floor, dtype=float)
    finite = np.isfinite(beta)
    result[finite] = np.exp(-np.square(np.abs(beta[finite]) / scale))
    return np.clip(result, floor, 1.0)


@validate_params(RELIABILITY_COMBINATION_SCHEMA)
def combine_observation_reliability(
    measurement: np.ndarray,
    dimensionality: np.ndarray,
    *,
    minimum: float = 0.0,
) -> np.ndarray:
    """Return bounded measurement-by-dimensionality reliability.

    Parameters
    ----------
    measurement, dimensionality : array-like of float
        Reliability factors in ``[0, 1]``. Inputs must be broadcastable to
        one common shape.
    minimum : float, default=0.0
        Lower bound applied after multiplication.

    Returns
    -------
    numpy.ndarray of float
        Broadcast product clipped to ``[minimum, 1]``.

    Raises
    ------
    ValueError
        If inputs are empty, non-finite, outside ``[0, 1]``, or cannot be
        broadcast together.

    Examples
    --------
    >>> combine_observation_reliability(
    ...     [0.8, 0.5], [0.5, 0.2]
    ... ).tolist()
    [0.4, 0.1]
    """
    noise = np.asarray(measurement, dtype=float)
    physics = np.asarray(dimensionality, dtype=float)
    floor = float(minimum)
    if noise.size == 0 or physics.size == 0:
        raise ValueError("reliability inputs must not be empty")
    if not np.isfinite(floor) or not 0 <= floor <= 1:
        raise ValueError("minimum must lie in [0, 1]")
    try:
        noise, physics = np.broadcast_arrays(noise, physics)
    except ValueError as exc:
        raise ValueError("reliability inputs are not broadcastable") from exc
    for name, values in (
        ("measurement", noise),
        ("dimensionality", physics),
    ):
        if not np.all(np.isfinite(values)) or np.any(
            (values < 0) | (values > 1)
        ):
            raise ValueError(f"{name} must be finite and lie in [0, 1]")
    return np.clip(noise * physics, floor, 1.0)
