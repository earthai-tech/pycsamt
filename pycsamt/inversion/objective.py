# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Objective-function helpers for inversion backends."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .doc import _inversion_param_docs

__all__ = [
    "ErrorModel",
    "component_errors",
    "component_mask",
    "error_model_from_config",
    "relative_errors",
    "weighted_rms",
]


@dataclass
class ErrorModel:
    rho_relative: float = 0.05
    rho_absolute: float = 1e-12
    phase_absolute: float = 3.0
    phase_relative: float = 0.0
    tdem_relative: float = 0.05
    tdem_absolute: float = 1e-30
    impedance_relative: float = 0.05
    impedance_absolute: float = 1e-12
    min_error: float = 1e-12
    masks: dict[str, Any] = field(default_factory=dict)

    def errors(
        self,
        values: Any,
        *,
        component: str,
        explicit: Any = None,
        relative: bool = False,
    ) -> np.ndarray:
        """Return data errors for one observed component.

        Parameters
        ----------
        values : array-like
            Observed or predicted values used to scale relative error floors.
        component : str
            Component name. Common aliases include ``"rho"``, ``"rho_a"``,
            ``"phase"``, ``"tdem"``, ``"dbdt"``, and ``"impedance"``.
        explicit : array-like, optional
            Explicit absolute errors. When provided, component floor settings
            are bypassed except for the global ``min_error`` lower bound.
        relative : bool, default False
            If ``True``, return relative errors by dividing absolute errors by
            ``abs(values)`` with ``min_error`` protection.

        Returns
        -------
        ndarray
            Error array with the same shape as ``values``.

        Examples
        --------
        >>> from pycsamt.inversion.objective import ErrorModel
        >>> model = ErrorModel(rho_relative=0.1, rho_absolute=1.0)
        >>> model.errors([5.0, 100.0], component="rho").tolist()
        [1.0, 10.0]
        >>> model.errors([45.0], component="phase").tolist()
        [3.0]
        """
        arr = np.asarray(values, dtype=float)
        comp = _canonical_component(component)
        if explicit is not None:
            err = np.asarray(explicit, dtype=float)
            if relative:
                err = err / np.maximum(np.abs(arr), self.min_error)
            return np.maximum(err, self.min_error)

        if comp == "rho":
            err = np.maximum(
                np.abs(arr) * self.rho_relative, self.rho_absolute
            )
        elif comp == "phase":
            err = np.maximum(
                np.abs(arr) * self.phase_relative, self.phase_absolute
            )
        elif comp == "tdem":
            err = np.maximum(
                np.abs(arr) * self.tdem_relative, self.tdem_absolute
            )
        elif comp in {"z_real", "z_imag", "impedance"}:
            err = np.maximum(
                np.abs(arr) * self.impedance_relative, self.impedance_absolute
            )
        else:
            err = np.maximum(np.abs(arr) * self.rho_relative, self.min_error)
        if relative:
            err = err / np.maximum(np.abs(arr), self.min_error)
        return np.maximum(err, self.min_error)

    def mask(self, values: Any, *, component: str) -> np.ndarray:
        """Return a boolean inclusion mask for one component.

        Parameters
        ----------
        values : array-like
            Data array whose shape defines the returned mask shape.
        component : str
            Component name. Component-specific masks are checked first, then
            station-level masks under ``"station"`` or ``"station_mask"``.

        Returns
        -------
        ndarray of bool
            Boolean mask broadcast to the shape of ``values``.

        Examples
        --------
        >>> from pycsamt.inversion.objective import ErrorModel
        >>> model = ErrorModel(masks={"station": [True, False]})
        >>> model.mask([[1.0, 2.0], [3.0, 4.0]], component="rho").tolist()
        [[True, True], [False, False]]
        """
        arr = np.asarray(values)
        comp = _canonical_component(component)
        raw = self.masks.get(comp, self.masks.get(component, None))
        if raw is None:
            raw = self.masks.get(
                "station", self.masks.get("station_mask", None)
            )
        if raw is None:
            return np.ones(arr.shape, dtype=bool)
        mask = np.asarray(raw, dtype=bool)
        if mask.shape == arr.shape:
            return mask
        if mask.ndim == 1 and arr.ndim == 2 and mask.size == arr.shape[0]:
            return np.broadcast_to(mask[:, None], arr.shape)
        if mask.ndim == 1 and arr.ndim == 2 and mask.size == arr.shape[1]:
            return np.broadcast_to(mask[None, :], arr.shape)
        return np.broadcast_to(mask, arr.shape)


ErrorModel.__doc__ = f"""
Component-aware data-error settings for inversion backends.

``ErrorModel`` centralizes the floors used when packing objective functions for
the built-in, SimPEG, pyGIMLi, Occam2D, and ModEM adapters. It supports
component-specific relative/absolute floors and boolean masks used to exclude
stations, frequencies, time gates, or individual samples.

Parameters
----------
rho_relative : float, default 0.05
    Relative apparent-resistivity error floor.
rho_absolute : float, default 1e-12
    Absolute apparent-resistivity error floor in ohm metres.
phase_absolute : float, default 3.0
    Absolute phase error floor in degrees.
phase_relative : float, default 0.0
    Optional relative phase error floor.
tdem_relative : float, default 0.05
    Relative TDEM decay-value error floor.
tdem_absolute : float, default 1e-30
    Absolute TDEM decay-value error floor.
impedance_relative : float, default 0.05
    Relative impedance-component error floor.
impedance_absolute : float, default 1e-12
    Absolute impedance-component error floor.
min_error : float, default 1e-12
    Global lower bound applied to all returned errors.
masks : dict, optional
    Component or station masks. Supported keys include component names such as
    ``"rho"``, ``"phase"``, ``"tdem"``, and station-level keys
    ``"station"``/``"station_mask"``.

Notes
-----
Configuration helpers read overrides from ``backend_options`` and nested
``backend_options["error_model"]`` dictionaries. Explicit observed-data errors
can still be passed to :meth:`errors` or :func:`component_errors`.

{_inversion_param_docs.errors.error_model_examples}

{_inversion_param_docs.errors.error_model_references}
"""


def relative_errors(values, floor: float = 0.05) -> np.ndarray:
    """Return absolute errors from a relative floor.

    Parameters
    ----------
    values : array-like
        Data values used for relative scaling.
    floor : float, default 0.05
        Relative floor as a fraction of absolute value.

    Returns
    -------
    ndarray
        ``max(abs(values) * floor, 1e-12)``.

    Examples
    --------
    >>> from pycsamt.inversion.objective import relative_errors
    >>> relative_errors([10.0, 100.0], floor=0.1).tolist()
    [1.0, 10.0]
    """
    arr = np.asarray(values, dtype=float)
    return np.maximum(np.abs(arr) * float(floor), 1e-12)


def error_model_from_config(cfg: Any) -> ErrorModel:
    """Build an :class:`ErrorModel` from inversion config options.

    Parameters
    ----------
    cfg : object
        Inversion config-like object exposing ``error_floor``, ``phase_error``,
        and optional ``backend_options``. Overrides may be placed directly in
        ``backend_options`` or inside ``backend_options["error_model"]``.

    Returns
    -------
    ErrorModel
        Component-aware error model.

    Examples
    --------
    >>> from pycsamt.inversion.config import InversionConfig
    >>> from pycsamt.inversion.objective import error_model_from_config
    >>> cfg = InversionConfig(error_floor=0.1, phase_error=2.0)
    >>> model = error_model_from_config(cfg)
    >>> model.rho_relative
    0.1
    >>> model.phase_absolute
    2.0
    """
    opts = dict(getattr(cfg, "backend_options", {}) or {})
    nested = dict(opts.get("error_model", {}) or {})

    def pick(name: str, default: Any) -> Any:
        return nested.get(name, opts.get(name, default))

    return ErrorModel(
        rho_relative=float(
            pick("rho_relative", getattr(cfg, "error_floor", 0.05))
        ),
        rho_absolute=float(pick("rho_absolute", 1e-12)),
        phase_absolute=float(
            pick("phase_absolute", getattr(cfg, "phase_error", 3.0))
        ),
        phase_relative=float(pick("phase_relative", 0.0)),
        tdem_relative=float(
            pick("tdem_relative", getattr(cfg, "error_floor", 0.05))
        ),
        tdem_absolute=float(pick("tdem_absolute", 1e-30)),
        impedance_relative=float(
            pick("impedance_relative", getattr(cfg, "error_floor", 0.05))
        ),
        impedance_absolute=float(pick("impedance_absolute", 1e-12)),
        min_error=float(pick("min_error", 1e-12)),
        masks=dict(pick("masks", opts.get("masks", {})) or {}),
    )


def component_errors(
    values: Any,
    cfg: Any,
    *,
    component: str,
    explicit: Any = None,
    relative: bool = False,
) -> np.ndarray:
    """Return errors for one data component using config settings.

    Parameters
    ----------
    values : array-like
        Data values used for scaling.
    cfg : object
        Inversion config-like object accepted by :func:`error_model_from_config`.
    component : str
        Component name or alias.
    explicit : array-like, optional
        Explicit data errors.
    relative : bool, default False
        Return relative errors instead of absolute errors.

    Returns
    -------
    ndarray
        Component error array.

    Examples
    --------
    >>> from pycsamt.inversion.config import InversionConfig
    >>> from pycsamt.inversion.objective import component_errors
    >>> cfg = InversionConfig(error_floor=0.05, phase_error=2.0)
    >>> component_errors([100.0, 200.0], cfg, component="rho").tolist()
    [5.0, 10.0]
    >>> component_errors([45.0], cfg, component="phase").tolist()
    [2.0]
    """
    return error_model_from_config(cfg).errors(
        values,
        component=component,
        explicit=explicit,
        relative=relative,
    )


def component_mask(values: Any, cfg: Any, *, component: str) -> np.ndarray:
    """Return a boolean mask for one component using config mask settings.

    Parameters
    ----------
    values : array-like
        Data array whose shape defines the returned mask shape.
    cfg : object
        Inversion config-like object accepted by :func:`error_model_from_config`.
    component : str
        Component name or alias.

    Returns
    -------
    ndarray of bool
        Boolean inclusion mask.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.inversion.config import InversionConfig
    >>> from pycsamt.inversion.objective import component_mask
    >>> cfg = InversionConfig(
    ...     backend_options={"masks": {"station": [True, False]}}
    ... )
    >>> component_mask(np.ones((2, 2)), cfg, component="rho").tolist()
    [[True, True], [False, False]]
    """
    return error_model_from_config(cfg).mask(values, component=component)


def weighted_rms(observed, predicted, errors=None) -> float:
    """Return normalized weighted RMS misfit.

    Parameters
    ----------
    observed, predicted : array-like
        Observed and predicted data arrays. They must be broadcast-compatible
        through NumPy conversion.
    errors : array-like, optional
        Absolute data errors. If omitted, unit errors are used.

    Returns
    -------
    float
        Root-mean-square of ``(predicted - observed) / errors`` over finite,
        positive-error samples. Returns ``nan`` when no valid samples exist.

    Examples
    --------
    >>> from pycsamt.inversion.objective import weighted_rms
    >>> round(weighted_rms([10.0, 20.0], [11.0, 18.0], [1.0, 2.0]), 6)
    1.0

    References
    ----------
    .. [weighted-rms-1] Tarantola, A. (2005). *Inverse Problem Theory and Methods for Model
       Parameter Estimation*. SIAM.
    """
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
    return float(np.sqrt(np.mean(resid**2)))


def _canonical_component(component: str) -> str:
    comp = str(component).lower()
    if comp in {"rho_a", "rhoa", "apparent_resistivity", "resistivity"}:
        return "rho"
    if comp in {"phi", "phase_deg"}:
        return "phase"
    if comp in {"values", "dbdt", "dbz_dt", "voltage", "tem"}:
        return "tdem"
    if comp in {"z", "zxy", "zyx", "impedance"}:
        return "impedance"
    return comp
