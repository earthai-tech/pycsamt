# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Gradient-boosted DUHI prejudice weights (experimental).

Motivation
----------
A family-specific DUHI-paper diagnostic (``multiple_body``: several
separate, alternating conductive/resistive compact bodies) found that
the M4-vs-M2 point-accuracy gap tracks how much *roughness* the
converged Occam model needs to satisfy data and prejudice jointly
(Spearman-like correlation 0.55-0.91 across families, tightest for
``multiple_body``). Occam trades roughness against misfit and against
the DUHI prejudice term through a *single global* Lagrange multiplier
per iteration -- it cannot apply a different smoothness locally. A
geology made of several separate compact anomalies is exactly the case
that most often needs non-uniform local smoothness, which the global
trade-off cannot supply, so the prejudice term is too weak relative to
roughness right where the AI ensemble has already resolved a sharp,
real boundary between two adjacent bodies.

:func:`pycsamt.ai.inversion.edge_roughness.compute_roughness_exceptions`
addresses this from the roughness side (relaxing the penalty between
AI-confident adjacent cells) and was found safe but weakly effective in
aggregate, and untested on the specific families/realizations that fail
worst. This module addresses the same trade-off from the *prejudice*
side instead: it locally strengthens the existing point-value prejudice
weight (already spatially varying through AI predictive uncertainty,
see :class:`pycsamt.ai.inversion.duhi2d.DUHIInverter2D`) wherever the
AI mean itself predicts a sharp local gradient, so the prejudice term
has more leverage to hold that boundary's shape against the roughness
penalty's pull toward smoothness. No Fortran changes are required; the
result is written through the same
:class:`pycsamt.models.occam2d.OccamPrejudice` sparse file already used
by :class:`pycsamt.ai.inversion.duhi2d.DUHIInverter2D`.

This is a genuinely different lever than roughness exceptions: it never
changes the roughness penalty itself, only how strongly the model is
pulled toward the AI's own predicted value at cells that sit on a real
AI-resolved boundary.

Entry points
------------
``gradient_boost_multiplier(local_gradient, gradient_scale, boost_gain, boost_cap)``
    Per-parameter boost formula.
``boost_prejudice_weights(model, ai_mean_parameters, prejudice_weights, ...)``
    Apply the boost over an entire Occam model.
"""

from __future__ import annotations

import numpy as np

from ..inversion.mapping2d import build_occam_parameter_adjacency
from ...models.occam2d import OccamModel

__all__ = ["gradient_boost_multiplier", "boost_prejudice_weights"]


def gradient_boost_multiplier(
    local_gradient: float,
    *,
    gradient_scale: float,
    boost_gain: float,
    boost_cap: float,
) -> float:
    r"""Return the prejudice-weight multiplier for one parameter.

    The multiplier grows with the parameter's local AI-mean gradient
    (the largest absolute log10-resistivity difference to any spatial
    neighbour):

    .. math::

        \mu
        =
        1 + g \, \min\!\left(
            \frac{|\nabla_{\mathrm{AI}}|}{\sigma_{\nabla}},\; c
        \right),

    where :math:`|\nabla_{\mathrm{AI}}|` is ``local_gradient``,
    :math:`\sigma_{\nabla}` is ``gradient_scale``, :math:`g` is
    ``boost_gain``, and :math:`c` is ``boost_cap``. A flat
    neighbourhood (``local_gradient == 0``) returns ``1.0`` (no
    change); a strong AI-resolved boundary approaches
    ``1 + boost_gain * boost_cap``.

    Parameters
    ----------
    local_gradient : float
        Non-negative largest absolute AI-mean difference between this
        parameter and any spatial neighbour, in log10 resistivity.
    gradient_scale : float
        Positive calibration scale. Gradients much smaller than this
        are treated as smooth background; gradients much larger
        saturate toward ``boost_cap``.
    boost_gain : float
        Non-negative gain controlling how strongly a saturated
        gradient amplifies the prejudice weight.
    boost_cap : float
        Non-negative ceiling on the normalized gradient term, keeping
        the multiplier bounded even for extreme gradients.

    Returns
    -------
    float
        Multiplier greater than or equal to ``1.0``.

    Raises
    ------
    ValueError
        Raised when ``local_gradient`` is negative or non-finite, or
        when ``gradient_scale`` is not positive and finite, or when
        ``boost_gain``/``boost_cap`` are negative or non-finite.

    See Also
    --------
    boost_prejudice_weights
        Applies this formula over every Occam parameter.

    Examples
    --------
    A flat neighbourhood leaves the weight unchanged:

    >>> from pycsamt.ai.inversion.gradient_prejudice import (
    ...     gradient_boost_multiplier,
    ... )
    >>> gradient_boost_multiplier(
    ...     0.0, gradient_scale=0.3, boost_gain=2.0, boost_cap=5.0,
    ... )
    1.0

    A strong AI-resolved boundary saturates toward the cap:

    >>> round(gradient_boost_multiplier(
    ...     10.0, gradient_scale=0.3, boost_gain=2.0, boost_cap=5.0,
    ... ), 3)
    11.0
    """
    if not (np.isfinite(local_gradient) and local_gradient >= 0):
        raise ValueError("local_gradient must be finite and non-negative")
    if not (np.isfinite(gradient_scale) and gradient_scale > 0):
        raise ValueError("gradient_scale must be finite and positive")
    if not (np.isfinite(boost_gain) and boost_gain >= 0):
        raise ValueError("boost_gain must be finite and non-negative")
    if not (np.isfinite(boost_cap) and boost_cap >= 0):
        raise ValueError("boost_cap must be finite and non-negative")

    normalized = min(local_gradient / gradient_scale, boost_cap)
    return float(1.0 + boost_gain * normalized)


def boost_prejudice_weights(
    model: OccamModel,
    ai_mean_parameters: np.ndarray,
    prejudice_weights: np.ndarray,
    *,
    gradient_scale: float = 0.3,
    boost_gain: float = 2.0,
    boost_cap: float = 5.0,
) -> np.ndarray:
    r"""Locally boost DUHI prejudice weights at AI-resolved boundaries.

    For every Occam parameter, computes the largest absolute AI-mean
    difference to any spatial neighbour (from
    :func:`pycsamt.ai.inversion.mapping2d.build_occam_parameter_adjacency`)
    and multiplies the existing prejudice weight by
    :func:`gradient_boost_multiplier` evaluated at that local gradient.
    Parameters with no neighbours (a degenerate single-parameter model)
    are left unmodified.

    Parameters
    ----------
    model : OccamModel
        Populated Occam model definition, matching the parameter order
        of ``ai_mean_parameters`` and ``prejudice_weights``.
    ai_mean_parameters : array-like of float, shape (model.n_params,)
        Mapped AI mean in Occam parameter order, e.g.
        ``DUHIInverter2D.ai_mean_parameters`` after
        :meth:`DUHIInverter2D.prepare`.
    prejudice_weights : array-like of float, shape (model.n_params,)
        Existing non-negative prejudice weights to boost, e.g.
        ``DUHIInverter2D.prejudice_weights``.
    gradient_scale : float, default 0.3
        Forwarded to :func:`gradient_boost_multiplier`.
    boost_gain : float, default 2.0
        Forwarded to :func:`gradient_boost_multiplier`.
    boost_cap : float, default 5.0
        Forwarded to :func:`gradient_boost_multiplier`.

    Returns
    -------
    numpy.ndarray of float, shape (model.n_params,)
        Boosted prejudice weights, always greater than or equal to the
        input weights element-wise.

    Raises
    ------
    ValueError
        Raised when ``ai_mean_parameters`` or ``prejudice_weights``
        length does not match ``model.n_params``, when either array is
        non-finite, or when ``prejudice_weights`` contains negative
        values.

    See Also
    --------
    gradient_boost_multiplier
        Per-parameter boost formula.
    pycsamt.ai.inversion.edge_roughness.compute_roughness_exceptions
        Addresses the same roughness/prejudice trade-off from the
        roughness side instead.
    pycsamt.models.occam2d.OccamPrejudice.from_dense
        Consumes the returned weights to rebuild a sparse prejudice
        file.

    Examples
    --------
    >>> from pycsamt.ai.inversion.gradient_prejudice import (
    ...     boost_prejudice_weights,
    ... )
    >>> from pycsamt.models.occam2d import OccamModel
    >>> model = OccamModel.read("occam_run/Occam2DModel")  # doctest: +SKIP
    >>> boosted = boost_prejudice_weights(  # doctest: +SKIP
    ...     model, ai_mean_parameters, prejudice_weights,
    ... )
    """
    mean = np.asarray(ai_mean_parameters, dtype=float).reshape(-1)
    weights = np.asarray(prejudice_weights, dtype=float).reshape(-1)
    n = int(model.n_params)
    if mean.size != n:
        raise ValueError(
            "ai_mean_parameters length does not match model.n_params"
        )
    if weights.size != n:
        raise ValueError(
            "prejudice_weights length does not match model.n_params"
        )
    if not np.all(np.isfinite(mean)):
        raise ValueError("ai_mean_parameters must be finite")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("prejudice_weights must be finite and non-negative")

    pairs = build_occam_parameter_adjacency(model)
    local_gradient = np.zeros(n, dtype=float)
    for i, j in pairs:
        diff = abs(mean[i - 1] - mean[j - 1])
        local_gradient[i - 1] = max(local_gradient[i - 1], diff)
        local_gradient[j - 1] = max(local_gradient[j - 1], diff)

    boosted = weights.copy()
    for idx in range(n):
        multiplier = gradient_boost_multiplier(
            local_gradient[idx],
            gradient_scale=gradient_scale,
            boost_gain=boost_gain,
            boost_cap=boost_cap,
        )
        boosted[idx] = weights[idx] * multiplier
    return boosted
