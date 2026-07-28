# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Complex-response residual diagnostics for EM forward comparisons.

These diagnostics break the ``L_response`` penalty from
:mod:`pycsamt.ai.losses.response` down by station, frequency, and
component, matching the Inversion and Field rows of the validation
matrix in the AI-inversion plan: "residual maps by station/
frequency/component".  Inputs share the canonical impedance shape
``(station, frequency, component)`` used by
:class:`~pycsamt.forward.maxwell.contracts.ForwardResult` and
:class:`~pycsamt.ai.data.contracts.SurveyData`.

All functions operate on plain NumPy arrays so the module stays
importable without an optional deep-learning backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from ..losses.response import ResponseLossResult, response_residual_loss

if TYPE_CHECKING:
    from pycsamt.ai.data import SurveyData
    from pycsamt.forward.maxwell import ForwardResult

__all__ = [
    "ResponseResidualReport",
    "response_residual_report",
    "response_residual_report_from_contracts",
]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    """Return a read-only copy of *value* as an ndarray."""
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _residual_magnitude(
    predicted: Any,
    observed: Any,
    *,
    errors: Any | None,
    valid: Any | None,
    kind: str,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Return per-cell residual magnitude, its mask, and normalized
    flag, mirroring :mod:`pycsamt.ai.losses.response` cell by cell.
    """
    if kind not in ("l1", "l2"):
        raise ValueError("kind must be 'l1' or 'l2'.")
    pred = np.asarray(predicted, dtype=complex)
    obs = np.asarray(observed, dtype=complex)
    if pred.shape != obs.shape:
        raise ValueError("predicted and observed must share one shape.")
    if pred.ndim != 3:
        raise ValueError(
            "predicted and observed must have canonical shape "
            "(station, frequency, component)."
        )
    if pred.size == 0:
        raise ValueError("predicted and observed must not be empty.")
    mask = np.isfinite(pred.real) & np.isfinite(pred.imag)
    mask &= np.isfinite(obs.real) & np.isfinite(obs.imag)
    normalized = errors is not None
    error_array = None
    if normalized:
        error_array = np.asarray(errors, dtype=float)
        if error_array.shape != pred.shape:
            raise ValueError(
                "errors must have the same shape as predicted/observed."
            )
        mask &= np.isfinite(error_array) & (error_array > 0)
    if valid is not None:
        supplied = np.asarray(valid, dtype=bool)
        if supplied.shape != pred.shape:
            raise ValueError(
                "valid must have the same shape as predicted/observed."
            )
        mask &= supplied
    if normalized:
        safe_errors = np.where(mask, error_array, 1.0)
        residual = np.where(mask, (pred - obs) / safe_errors, 0.0)
    else:
        residual = np.where(mask, pred - obs, 0.0)
    if kind == "l1":
        per_cell = np.abs(residual)
    else:
        per_cell = np.square(residual.real) + np.square(residual.imag)
    return per_cell, mask, normalized


def _axis_mean(
    per_cell: np.ndarray, mask: np.ndarray, keep_axis: int
) -> np.ndarray:
    """Return the masked mean of ``per_cell``, keeping only
    ``keep_axis`` and averaging over the other two axes; ``nan``
    where an axis position has no valid cell.
    """
    reduce_axes = tuple(
        axis for axis in range(per_cell.ndim) if axis != keep_axis
    )
    weight_sum = mask.sum(axis=reduce_axes, dtype=float)
    total = np.where(mask, per_cell, 0.0).sum(axis=reduce_axes)
    safe_weight = np.where(weight_sum > 0, weight_sum, 1.0)
    result = total / safe_weight
    return np.where(weight_sum > 0, result, np.nan)


@dataclass(frozen=True)
class ResponseResidualReport:
    """Immutable per-axis complex-impedance residual diagnostics.

    Parameters
    ----------
    overall : ResponseLossResult
        Reduced ``L_response`` penalty over every included cell.
    by_station, by_frequency, by_component : ndarray
        Masked mean per-cell penalty aggregated over the other two
        axes, in the same units as ``overall.value`` (e.g. mean
        squared residual for ``kind="l2"``, not RMS). ``nan`` where
        an axis position has no valid cell.
    station_names : tuple of str or None
        Optional station labels, length ``shape[0]``.
    frequencies_hz : ndarray or None
        Optional frequency labels, length ``shape[1]``.
    components : tuple of str or None
        Optional component labels, length ``shape[2]``.
    shape : tuple of int
        ``(station, frequency, component)`` shape of the compared
        arrays.

    Examples
    --------
    >>> import numpy as np
    >>> pred = np.array([[[1 + 0j], [2 + 0j]], [[0j], [0j]]])
    >>> obs = np.array([[[1 + 0j], [0j]], [[0j], [3 + 0j]]])
    >>> report = response_residual_report(pred, obs)
    >>> report.shape
    (2, 2, 1)
    >>> report.by_station.tolist()
    [2.0, 4.5]
    """

    overall: ResponseLossResult
    by_station: np.ndarray
    by_frequency: np.ndarray
    by_component: np.ndarray
    station_names: tuple[str, ...] | None
    frequencies_hz: np.ndarray | None
    components: tuple[str, ...] | None
    shape: tuple[int, int, int]

    def __post_init__(self) -> None:
        if not isinstance(self.overall, ResponseLossResult):
            raise TypeError("overall must be a ResponseLossResult.")
        shape = tuple(int(size) for size in self.shape)
        if len(shape) != 3:
            raise ValueError(
                "shape must be (station, frequency, component)."
            )
        object.__setattr__(self, "shape", shape)
        axes = {
            "by_station": (self.by_station, shape[0]),
            "by_frequency": (self.by_frequency, shape[1]),
            "by_component": (self.by_component, shape[2]),
        }
        for name, (array, expected_len) in axes.items():
            array = np.asarray(array, dtype=float)
            if array.shape != (expected_len,):
                raise ValueError(
                    f"{name} must have shape ({expected_len},)."
                )
            object.__setattr__(self, name, _readonly(array))
        if (
            self.station_names is not None
            and len(self.station_names) != shape[0]
        ):
            raise ValueError("station_names must have length shape[0].")
        if self.station_names is not None:
            object.__setattr__(
                self,
                "station_names",
                tuple(str(name) for name in self.station_names),
            )
        if self.frequencies_hz is not None:
            frequencies = np.asarray(self.frequencies_hz, dtype=float)
            if frequencies.shape != (shape[1],):
                raise ValueError(
                    "frequencies_hz must have shape (shape[1],)."
                )
            object.__setattr__(
                self, "frequencies_hz", _readonly(frequencies)
            )
        if (
            self.components is not None
            and len(self.components) != shape[2]
        ):
            raise ValueError("components must have length shape[2].")
        if self.components is not None:
            object.__setattr__(
                self,
                "components",
                tuple(str(value) for value in self.components),
            )


def response_residual_report(
    predicted: Any,
    observed: Any,
    *,
    errors: Any | None = None,
    valid: Any | None = None,
    kind: str = "l2",
    station_names: Any | None = None,
    frequencies_hz: Any | None = None,
    components: Any | None = None,
) -> ResponseResidualReport:
    """Break complex-impedance residuals down by station/frequency/
    component.

    Parameters
    ----------
    predicted : array-like of complex
        Forward-simulated impedance, canonical shape
        ``(station, frequency, component)``.
    observed : array-like of complex
        Observed impedance with the same shape.
    errors : array-like or None, optional
        Positive absolute standard errors used to normalize
        residuals, as in
        :func:`~pycsamt.ai.losses.response.response_residual_loss`.
    valid : array-like of bool or None, optional
        Explicit observation mask.
    kind : {"l1", "l2"}, default="l2"
        Elementwise penalty applied to each residual magnitude.
    station_names, components : sequence of str or None, optional
        Optional axis labels attached to the returned report.
    frequencies_hz : array-like or None, optional
        Optional frequency labels attached to the returned report.

    Returns
    -------
    ResponseResidualReport
        Combined per-axis residual diagnostics.

    Examples
    --------
    >>> import numpy as np
    >>> pred = np.array([[[1 + 0j], [2 + 0j]], [[0j], [0j]]])
    >>> obs = np.array([[[1 + 0j], [0j]], [[0j], [3 + 0j]]])
    >>> report = response_residual_report(pred, obs)
    >>> report.by_frequency.tolist()
    [0.0, 6.5]
    """
    overall = response_residual_loss(
        predicted, observed, errors=errors, valid=valid, kind=kind
    )
    per_cell, mask, _ = _residual_magnitude(
        predicted, observed, errors=errors, valid=valid, kind=kind
    )
    return ResponseResidualReport(
        overall=overall,
        by_station=_axis_mean(per_cell, mask, 0),
        by_frequency=_axis_mean(per_cell, mask, 1),
        by_component=_axis_mean(per_cell, mask, 2),
        station_names=station_names,
        frequencies_hz=frequencies_hz,
        components=components,
        shape=per_cell.shape,
    )


def response_residual_report_from_contracts(
    forward: ForwardResult,
    observed: SurveyData,
    *,
    kind: str = "l2",
    use_errors: bool = True,
) -> ResponseResidualReport:
    """Build a :class:`ResponseResidualReport` directly from a
    ``ForwardResult``/``SurveyData`` pair.

    Requires exact station, component, and frequency alignment, as
    in :func:`~pycsamt.ai.losses.response.response_loss_from_contracts`.

    Parameters
    ----------
    forward : ForwardResult
        Predicted impedance from a Maxwell backend adapter.
    observed : SurveyData
        Observed survey impedance to compare against.
    kind : {"l1", "l2"}, default="l2"
        Elementwise penalty applied to each residual magnitude.
    use_errors : bool, default=True
        Normalize residuals by ``observed.impedance_error`` when it
        is available.

    Returns
    -------
    ResponseResidualReport
        Combined per-axis residual diagnostics, labeled with the
        survey's station names, frequencies, and components.
    """
    from ..losses.response import response_loss_from_contracts

    overall = response_loss_from_contracts(
        forward, observed, kind=kind, use_errors=use_errors
    )
    errors = observed.impedance_error if use_errors else None
    mask = forward.valid & observed.valid
    per_cell, combined_mask, _ = _residual_magnitude(
        forward.impedance_v_a,
        observed.impedance,
        errors=errors,
        valid=mask,
        kind=kind,
    )
    return ResponseResidualReport(
        overall=overall,
        by_station=_axis_mean(per_cell, combined_mask, 0),
        by_frequency=_axis_mean(per_cell, combined_mask, 1),
        by_component=_axis_mean(per_cell, combined_mask, 2),
        station_names=observed.station_names,
        frequencies_hz=observed.frequencies_hz,
        components=observed.components,
        shape=per_cell.shape,
    )
