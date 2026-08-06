# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regularization controls for inversion workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject
from .doc import _inversion_param_docs

__all__ = [
    "Regularization",
    "regularization_from_config",
    "regularization_residual",
    "regularization_weight",
    "pygimli_lambda",
]


@dataclass
class Regularization(PyCSAMTObject, MetadataMixin):
    kind: str = "smooth"
    alpha_s: float = 1.0
    alpha_x: float = 1.0
    alpha_z: float = 1.0
    reference_weight: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.kind = str(self.kind).lower()
        self.alpha_s = float(self.alpha_s)
        self.alpha_x = float(self.alpha_x)
        self.alpha_z = float(self.alpha_z)
        self.reference_weight = float(self.reference_weight)
        self.validate()

    def validate(self) -> None:
        """Validate regularization kind and non-negative weights.

        Raises
        ------
        ValueError
            If ``kind`` is not one of ``"none"``, ``"smooth"``,
            ``"damped"``, or ``"blocky"``, or if any alpha/reference weight is
            negative.

        Examples
        --------
        >>> from pycsamt.inversion.regularization import Regularization
        >>> Regularization(kind="smooth", alpha_x=1.0).validate() is None
        True
        """
        if self.kind not in {"none", "smooth", "damped", "blocky"}:
            raise ValueError(
                "regularization kind must be none/smooth/damped/blocky."
            )
        for name in ("alpha_s", "alpha_x", "alpha_z", "reference_weight"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be non-negative.")


Regularization.__doc__ = rf"""
Backend-neutral regularization settings.

``Regularization`` stores the model-structure penalty vocabulary shared by
built-in, SimPEG, and pyGIMLi inversion paths. The residual helpers operate on
model-parameter arrays, usually log-domain resistivity values, and return
unweighted structure residuals that can be multiplied by a global scalar such
as :math:`\sqrt{{\lambda}}`.

For a smooth model :math:`m`, the penalty contains finite differences such as

.. math::

    \phi_m = \alpha_x \|\nabla_x m\|_2^2
             + \alpha_z \|\nabla_z m\|_2^2 .

For damped/reference regularization, a smallness term
:math:`\alpha_s\|m-m_{{ref}}\|_2^2` is added. For ``"blocky"`` models, the
finite differences are normalized by
:math:`\sqrt{{(\Delta m)^2 + \epsilon^2}}` to reduce sensitivity to sharp
edges.

Parameters
----------
{_inversion_param_docs.regularization.kind}
{_inversion_param_docs.regularization.alpha_s}
{_inversion_param_docs.regularization.alpha_x}
{_inversion_param_docs.regularization.alpha_z}
{_inversion_param_docs.regularization.reference_weight}
{_inversion_param_docs.regularization.metadata}

Notes
-----
This class stores relative penalty settings only. The scalar objective weight
is read separately by :func:`regularization_weight`, while pyGIMLi-specific
``lambda`` handling is read by :func:`pygimli_lambda`.

{_inversion_param_docs.regularization.residual_examples}

{_inversion_param_docs.regularization.references}
"""


def regularization_from_config(cfg: Any) -> Regularization:
    """Build backend-neutral controls from a config object.

    Parameters
    ----------
    cfg : object
        Inversion config-like object exposing ``regularization`` and optional
        ``backend_options``. Recognized backend options include ``alpha_s``,
        ``alpha_model``, ``alpha_x``, ``alpha_lateral``, ``alpha_z``,
        ``alpha_vertical``, ``reference_weight``, and
        ``regularization_metadata``.

    Returns
    -------
    Regularization
        Parsed regularization settings.

    Examples
    --------
    >>> from pycsamt.inversion.config import InversionConfig
    >>> from pycsamt.inversion.regularization import regularization_from_config
    >>> cfg = InversionConfig(
    ...     regularization="blocky",
    ...     backend_options={"alpha_x": 2.0, "alpha_z": 0.5},
    ... )
    >>> reg = regularization_from_config(cfg)
    >>> (reg.kind, reg.alpha_x, reg.alpha_z)
    ('blocky', 2.0, 0.5)
    """
    opts = dict(getattr(cfg, "backend_options", {}) or {})
    metadata = dict(opts.get("regularization_metadata", {}))
    return Regularization(
        kind=getattr(cfg, "regularization", "smooth"),
        alpha_s=float(opts.get("alpha_s", opts.get("alpha_model", 1.0))),
        alpha_x=float(opts.get("alpha_x", opts.get("alpha_lateral", 1.0))),
        alpha_z=float(opts.get("alpha_z", opts.get("alpha_vertical", 1.0))),
        reference_weight=float(opts.get("reference_weight", 0.0)),
        metadata=metadata,
    )


def regularization_weight(cfg: Any, *, default: float = 0.0) -> float:
    """Return the scalar penalty weight for least-squares backends.

    Parameters
    ----------
    cfg : object
        Inversion config-like object with optional ``backend_options``.
    default : float, default 0.0
        Value returned when no shared weight is configured.

    Returns
    -------
    float
        Value from ``backend_options["regularization_weight"]`` or
        ``backend_options["reg_weight"]``.

    Examples
    --------
    >>> from pycsamt.inversion.config import InversionConfig
    >>> from pycsamt.inversion.regularization import regularization_weight
    >>> cfg = InversionConfig(backend_options={"regularization_weight": 3.0})
    >>> regularization_weight(cfg)
    3.0
    """
    opts = dict(getattr(cfg, "backend_options", {}) or {})
    return float(
        opts.get("regularization_weight", opts.get("reg_weight", default))
    )


def pygimli_lambda(cfg: Any, *, default: float = 20.0) -> float:
    """Return the pyGIMLi lambda value using shared option names.

    Parameters
    ----------
    cfg : object
        Inversion config-like object with optional ``backend_options``.
    default : float, default 20.0
        Fallback value when no lambda-like option is configured.

    Returns
    -------
    float
        Value from ``lam``, ``lambda``, or ``regularization_weight`` in that
        priority order.

    Examples
    --------
    >>> from pycsamt.inversion.config import InversionConfig
    >>> from pycsamt.inversion.regularization import pygimli_lambda
    >>> cfg = InversionConfig(backend_options={"lam": 7.0})
    >>> pygimli_lambda(cfg)
    7.0
    """
    opts = dict(getattr(cfg, "backend_options", {}) or {})
    return float(
        opts.get(
            "lam",
            opts.get("lambda", opts.get("regularization_weight", default)),
        )
    )


def regularization_residual(
    values: Any,
    *,
    reference: Any | None = None,
    regularization: Regularization | None = None,
    blocky_eps: float = 1e-2,
    axes: tuple[str, ...] | None = None,
) -> np.ndarray:
    """Return unweighted residual terms for smooth/damped/blocky penalties.

    ``values`` are normally log-domain model parameters. The returned vector
    already contains the square-root alpha factors; callers can multiply by
    a global scalar such as ``sqrt(regularization_weight)``.

    Parameters
    ----------
    values : array-like
        Model parameter array. One-dimensional arrays are treated as depth
        profiles. Two-dimensional arrays are treated as ``(z, x)`` sections.
    reference : array-like, optional
        Reference model with the same shape as ``values`` or broadcastable to
        that shape.
    regularization : Regularization, optional
        Regularization settings. If omitted, ``Regularization()`` is used.
    blocky_eps : float, default 1e-2
        Stabilization value in the blocky normalized-gradient penalty.
    axes : tuple of {"x", "z"}, optional
        Axes on which to compute roughness. Defaults to ``("z",)`` for 1-D
        arrays and ``("x", "z")`` for 2-D or higher arrays.

    Returns
    -------
    ndarray
        Concatenated residual vector. Empty when ``kind="none"`` or all active
        penalty terms have zero weight.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.inversion.regularization import Regularization
    >>> from pycsamt.inversion.regularization import regularization_residual
    >>> reg = Regularization(kind="damped", alpha_s=1.0, alpha_z=1.0)
    >>> residual = regularization_residual(
    ...     np.array([1.0, 2.0, 4.0]), reference=np.ones(3), regularization=reg
    ... )
    >>> residual.size
    5

    References
    ----------
    .. [regularization-residual-1] Tikhonov, A. N. and Arsenin, V. Y. (1977). *Solutions of
       Ill-Posed Problems*. Winston.
    .. [regularization-residual-2] Farquharson, C. G. and Oldenburg, D. W. (1998). Non-linear
       inversion using general measures of data misfit and model structure.
       *Geophysical Journal International*, 134(1), 213-227.
    """
    reg = regularization or Regularization()
    arr = np.asarray(values, dtype=float)
    if reg.kind == "none":
        return np.array([], dtype=float)
    if reference is None:
        ref = np.zeros_like(arr)
    else:
        ref = np.asarray(reference, dtype=float)
        if ref.shape != arr.shape:
            ref = np.broadcast_to(ref, arr.shape)

    parts: list[np.ndarray] = []
    alpha_s = reg.alpha_s
    if reg.kind == "damped" and alpha_s == 0.0:
        alpha_s = 1.0
    if alpha_s > 0.0:
        ref_weight = max(
            float(reg.reference_weight), 1.0 if reference is not None else 0.0
        )
        scale = np.sqrt(
            alpha_s * max(ref_weight, 1.0 if reg.kind == "damped" else 0.0)
        )
        if scale > 0.0:
            parts.append(scale * (arr - ref).reshape(-1))

    if reg.kind in {"smooth", "damped", "blocky"}:
        selected_axes = axes or _default_axes(arr.ndim)
        axis_map = {"x": -1, "z": 0}
        alpha_map = {"x": reg.alpha_x, "z": reg.alpha_z}
        for name in selected_axes:
            axis = axis_map.get(name)
            alpha = alpha_map.get(name, 0.0)
            if axis is None or alpha <= 0.0 or arr.shape[axis] <= 1:
                continue
            diff = np.diff(arr, axis=axis)
            if reg.kind == "blocky":
                diff = diff / np.sqrt(diff**2 + blocky_eps**2)
            parts.append(np.sqrt(alpha) * diff.reshape(-1))

    if not parts:
        return np.array([], dtype=float)
    return np.concatenate(parts)


def _default_axes(ndim: int) -> tuple[str, ...]:
    if ndim <= 1:
        return ("z",)
    return ("x", "z")
