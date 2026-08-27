# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Uncertainty-adaptive roughness exceptions for DUHI (experimental).

Motivation
----------
Two controlled diagnostics (documented in the DUHI-paper reproduction
workspace, not part of the manuscript's frozen results) established
that the gap between DUHI (M4) and the direct AI point estimate (M2)
is not explained by either of the two cheapest hypotheses:

* blanket-disabling Occam's roughness penalty (``roughness_type=0``)
  does not help on average and is catastrophic for AI-untrustworthy
  geology, so roughness is not simply "fighting" the AI prejudice
  term everywhere; and
* denser synthetic surveys make classical recovery *worse*, not
  better, ruling out survey sparsity as the dominant limitation.

Both point toward a representational mismatch instead: Occam's
uniform first-difference roughness penalty cannot express the sharp,
piecewise structure the AI ensemble often recovers correctly, and a
uniform on/off toggle cannot distinguish where that structure should
be trusted from where it should not.

This module computes a spatially *selective* relaxation instead,
using the native OCCAM 3.0 roughness-penalty "exceptions" mechanism
(:attr:`pycsamt.models.occam2d.OccamModel.exceptions`, already read by
the bundled compiled solver but never previously populated by
PyCSAMT). For each pair of spatially adjacent Occam parameters, the
roughness penalty between them is relaxed only when the AI ensemble is
confident about *both* cells; it is left at its default strength
whenever either cell is uncertain. No Fortran changes are required.

Entry points
------------
``roughness_exception_weight(std_i, std_j, sigma_ai_floor, roughness_weight_floor)``
    Per-pair relaxation formula.
``compute_roughness_exceptions(model, ai_std_parameters, ...)``
    Build the full exception list for one Occam model.
"""

from __future__ import annotations

import numpy as np

from ..inversion.mapping2d import build_occam_parameter_adjacency
from ...models.occam2d import OccamModel

__all__ = ["roughness_exception_weight", "compute_roughness_exceptions"]


def roughness_exception_weight(
    std_i: float,
    std_j: float,
    *,
    sigma_ai_floor: float,
    roughness_weight_floor: float,
) -> float:
    r"""Return the roughness-penalty multiplier for one adjacent pair.

    Each parameter's AI confidence is

    .. math::

        \gamma
        =
        \frac{1}{1 + (\sigma_{\mathrm{AI}} / \sigma_0)^2}
        \in (0, 1],

    which approaches one as the predictive standard deviation
    :math:`\sigma_{\mathrm{AI}}` vanishes and approaches zero as it
    grows large relative to the floor :math:`\sigma_0`. The pair
    confidence is the weaker of its two endpoints,
    :math:`\gamma_{\mathrm{pair}} = \min(\gamma_i, \gamma_j)`, and the
    returned multiplier is

    .. math::

        \mathrm{expen}
        =
        1 - (1 - w_{\min})\,\gamma_{\mathrm{pair}},

    which equals ``1.0`` (the unmodified default penalty) when either
    endpoint is uncertain and approaches ``roughness_weight_floor``
    only when both endpoints are simultaneously confident.

    Parameters
    ----------
    std_i, std_j : float
        Non-negative mapped AI predictive standard deviations (log10
        resistivity) at the two adjacent Occam parameters.
    sigma_ai_floor : float
        Positive uncertainty scale. Values much smaller than this
        floor are treated as confident; values much larger are treated
        as uncertain. Matches the floor already used by
        :class:`pycsamt.ai.inversion.duhi2d.DUHIInverter2D` for the
        model-space prejudice weight, so both mechanisms respond to
        the same notion of "confident."
    roughness_weight_floor : float
        Smallest allowed multiplier, in ``[0, 1]``. A blanket
        diagnostic (uniformly disabling roughness) was found to be
        harmful; this floor keeps the relaxation partial rather than
        eliminating the penalty outright even in the most confident
        regions.

    Returns
    -------
    float
        Multiplier in ``[roughness_weight_floor, 1]``.

    Raises
    ------
    ValueError
        Raised when ``std_i``/``std_j`` are negative or non-finite,
        ``sigma_ai_floor`` is not positive and finite, or
        ``roughness_weight_floor`` is outside ``[0, 1]``.

    See Also
    --------
    compute_roughness_exceptions
        Applies this formula over every adjacent parameter pair.

    Examples
    --------
    Confident on both sides -- relaxed toward the floor:

    >>> from pycsamt.ai.inversion.edge_roughness import (
    ...     roughness_exception_weight,
    ... )
    >>> round(roughness_exception_weight(
    ...     0.01, 0.01, sigma_ai_floor=0.1, roughness_weight_floor=0.3,
    ... ), 3)
    0.307

    One uncertain endpoint -- effectively unchanged:

    >>> round(roughness_exception_weight(
    ...     0.01, 5.0, sigma_ai_floor=0.1, roughness_weight_floor=0.3,
    ... ), 3)
    1.0
    """
    if not (np.isfinite(std_i) and std_i >= 0):
        raise ValueError("std_i must be finite and non-negative")
    if not (np.isfinite(std_j) and std_j >= 0):
        raise ValueError("std_j must be finite and non-negative")
    if not (np.isfinite(sigma_ai_floor) and sigma_ai_floor > 0):
        raise ValueError("sigma_ai_floor must be finite and positive")
    if not (np.isfinite(roughness_weight_floor) and 0.0 <= roughness_weight_floor <= 1.0):
        raise ValueError("roughness_weight_floor must lie in [0, 1]")

    confidence_i = 1.0 / (1.0 + (std_i / sigma_ai_floor) ** 2)
    confidence_j = 1.0 / (1.0 + (std_j / sigma_ai_floor) ** 2)
    pair_confidence = min(confidence_i, confidence_j)
    return float(1.0 - (1.0 - roughness_weight_floor) * pair_confidence)


def compute_roughness_exceptions(
    model: OccamModel,
    ai_std_parameters: np.ndarray,
    *,
    sigma_ai_floor: float = 0.1,
    roughness_weight_floor: float = 0.3,
    active_threshold: float = 0.999,
) -> list[tuple[int, int, float]]:
    r"""Build native Occam roughness exceptions from mapped AI std.

    Computes :func:`roughness_exception_weight` for every spatially
    adjacent Occam parameter pair (from
    :func:`pycsamt.ai.inversion.mapping2d.build_occam_parameter_adjacency`)
    and keeps only pairs whose multiplier differs meaningfully from
    the default of one, so the written exception list stays sparse.

    Parameters
    ----------
    model : OccamModel
        Populated Occam model definition, matching the parameter order
        of ``ai_std_parameters``.
    ai_std_parameters : array-like of float, shape (model.n_params,)
        Mapped AI predictive standard deviation in Occam parameter
        order, e.g. ``DUHIInverter2D.ai_std_parameters`` after
        :meth:`DUHIInverter2D.prepare`.
    sigma_ai_floor : float, default 0.1
        Forwarded to :func:`roughness_exception_weight`.
    roughness_weight_floor : float, default 0.3
        Forwarded to :func:`roughness_exception_weight`.
    active_threshold : float, default 0.999
        Pairs whose multiplier is at least this close to ``1.0`` are
        omitted, keeping the exception list sparse (matching how
        :meth:`OccamPrejudice.from_dense` omits zero-weight records).

    Returns
    -------
    list of (int, int, float)
        One-based ``(brick_i, brick_j, expen)`` triples suitable for
        :attr:`pycsamt.models.occam2d.OccamModel.exceptions`.

    Raises
    ------
    ValueError
        Raised when ``ai_std_parameters`` length does not match
        ``model.n_params``, or contains negative or non-finite values.

    See Also
    --------
    roughness_exception_weight
        Per-pair relaxation formula.
    pycsamt.ai.inversion.mapping2d.build_occam_parameter_adjacency
        Supplies the candidate pairs.
    pycsamt.models.occam2d.OccamModel.exceptions
        Consumes the returned triples.

    Examples
    --------
    >>> from pycsamt.ai.inversion.edge_roughness import (
    ...     compute_roughness_exceptions,
    ... )
    >>> from pycsamt.models.occam2d import OccamModel
    >>> model = OccamModel.read("occam_run/Occam2DModel")  # doctest: +SKIP
    >>> exceptions = compute_roughness_exceptions(  # doctest: +SKIP
    ...     model, ai_std_parameters,
    ... )
    >>> model.exceptions = exceptions  # doctest: +SKIP
    >>> model.write("occam_run/Occam2DModel")  # doctest: +SKIP
    """
    std = np.asarray(ai_std_parameters, dtype=float).reshape(-1)
    if std.size != int(model.n_params):
        raise ValueError(
            "ai_std_parameters length does not match model.n_params"
        )
    if not np.all(np.isfinite(std)) or np.any(std < 0):
        raise ValueError("ai_std_parameters must be finite and non-negative")
    if not (np.isfinite(active_threshold) and 0.0 <= active_threshold < 1.0):
        raise ValueError("active_threshold must lie in [0, 1)")

    pairs = build_occam_parameter_adjacency(model)
    exceptions: list[tuple[int, int, float]] = []
    for i, j in pairs:
        expen = roughness_exception_weight(
            std[i - 1],
            std[j - 1],
            sigma_ai_floor=sigma_ai_floor,
            roughness_weight_floor=roughness_weight_floor,
        )
        if expen <= active_threshold:
            exceptions.append((i, j, expen))
    return exceptions
