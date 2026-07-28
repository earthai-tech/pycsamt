# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Electromagnetic response-consistency losses.

These losses implement the ``L_response`` term of the staged
inversion objective in the AI-inversion plan::

    L = w_m * L_model + lambda_x * L_grad_x + lambda_z * L_grad_z
        + lambda_tv * L_TV + lambda_d * L_response

Inputs are predicted and observed complex impedance arrays sharing
the canonical shape ``(station, frequency, component)`` used by
:class:`~pycsamt.forward.maxwell.contracts.ForwardResult` and
:class:`~pycsamt.ai.data.contracts.SurveyData`.  When positive
standard errors are supplied, residuals are normalized by them
before the elementwise penalty is applied, matching the usual
normalized-RMS EM data-misfit convention.

All functions operate on plain NumPy arrays so the module stays
importable without an optional deep-learning backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from ._common import _reduce

if TYPE_CHECKING:
    from pycsamt.ai.data import SurveyData
    from pycsamt.forward.maxwell import ForwardResult

__all__ = [
    "ResponseLossResult",
    "ResponseLoss",
    "response_residual_loss",
    "response_loss_from_contracts",
]

_KINDS = ("l1", "l2")
_REDUCTIONS = ("mean", "sum")


@dataclass(frozen=True)
class ResponseLossResult:
    """Immutable scalar result of a response-consistency loss.

    Parameters
    ----------
    value : float
        Reduced loss value. ``nan`` when no cell was included and
        ``reduction="mean"``.
    kind : {"l1", "l2"}
        Elementwise penalty applied to each residual magnitude.
    reduction : {"mean", "sum"}
        Reduction applied over valid cells.
    n_valid : int
        Number of impedance cells included after masking.
    weight_sum : float
        Number of included cells; kept for interface parity with
        other loss results in this package.
    normalized : bool
        Whether residuals were divided by a positive standard error
        before the penalty was applied.

    Examples
    --------
    >>> import numpy as np
    >>> pred = np.array([1 + 1j, 2 + 2j])
    >>> obs = np.array([1 + 1j, 0 + 0j])
    >>> response_residual_loss(pred, obs).value
    4.0
    """

    value: float
    kind: str
    reduction: str
    n_valid: int
    weight_sum: float
    normalized: bool

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise ValueError(f"kind must be one of {_KINDS}.")
        if self.reduction not in _REDUCTIONS:
            raise ValueError(f"reduction must be one of {_REDUCTIONS}.")
        if self.n_valid < 0:
            raise ValueError("n_valid must be non-negative.")
        object.__setattr__(self, "value", float(self.value))
        object.__setattr__(self, "weight_sum", float(self.weight_sum))
        object.__setattr__(self, "n_valid", int(self.n_valid))
        object.__setattr__(self, "normalized", bool(self.normalized))


def _response_residual_loss(
    predicted: Any,
    observed: Any,
    *,
    errors: Any | None,
    valid: Any | None,
    kind: str,
    reduction: str,
) -> ResponseLossResult:
    """Compute one masked, optionally error-normalized response loss."""
    if kind not in _KINDS:
        raise ValueError(f"kind must be one of {_KINDS}.")
    pred = np.asarray(predicted, dtype=complex)
    obs = np.asarray(observed, dtype=complex)
    if pred.shape != obs.shape:
        raise ValueError("predicted and observed must share one shape.")
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
    active_weights = mask.astype(float)
    value, weight_sum = _reduce(per_cell, active_weights, reduction)
    return ResponseLossResult(
        value=value,
        kind=kind,
        reduction=reduction,
        n_valid=int(np.count_nonzero(mask)),
        weight_sum=weight_sum,
        normalized=normalized,
    )


def response_residual_loss(
    predicted: Any,
    observed: Any,
    *,
    errors: Any | None = None,
    valid: Any | None = None,
    kind: str = "l2",
    reduction: str = "mean",
) -> ResponseLossResult:
    """Compare predicted and observed complex impedance responses.

    Parameters
    ----------
    predicted : array-like of complex
        Forward-simulated impedance. Any shape is accepted; the
        canonical layout is
        ``(station, frequency, component)``.
    observed : array-like of complex
        Observed impedance with the same shape.
    errors : array-like or None, optional
        Positive absolute standard errors, same shape as
        ``predicted``. When given, each residual is divided by its
        error before the elementwise penalty, matching the
        normalized-RMS EM data-misfit convention. Entries with a
        non-finite or non-positive error are excluded rather than
        raising.
    valid : array-like of bool or None, optional
        Explicit observation mask, combined with finite-value
        masking of ``predicted``, ``observed``, and ``errors``.
    kind : {"l1", "l2"}, default="l2"
        Elementwise penalty applied to each residual magnitude.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over included cells.

    Returns
    -------
    ResponseLossResult
        Reduced response-consistency penalty, ``L_response``.

    Examples
    --------
    >>> import numpy as np
    >>> pred = np.array([1 + 1j, 2 + 2j])
    >>> obs = np.array([1 + 1j, 0 + 0j])
    >>> response_residual_loss(pred, obs, kind="l2").value
    4.0
    >>> response_residual_loss(
    ...     pred, obs, errors=np.array([1.0, 2.0]), kind="l2"
    ... ).value
    1.0
    """
    return _response_residual_loss(
        predicted,
        observed,
        errors=errors,
        valid=valid,
        kind=kind,
        reduction=reduction,
    )


def response_loss_from_contracts(
    forward: ForwardResult,
    observed: SurveyData,
    *,
    kind: str = "l2",
    reduction: str = "mean",
    use_errors: bool = True,
) -> ResponseLossResult:
    """Compute ``L_response`` directly from canonical result/survey.

    Requires exact station, component, and frequency alignment
    between ``forward`` and ``observed`` rather than silently
    interpolating or reordering either axis, per the survey-matching
    principle in the AI-inversion plan.

    Parameters
    ----------
    forward : ForwardResult
        Predicted impedance from a Maxwell backend adapter.
    observed : SurveyData
        Observed survey impedance to compare against.
    kind : {"l1", "l2"}, default="l2"
        Elementwise penalty applied to each residual magnitude.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over included cells.
    use_errors : bool, default=True
        Normalize residuals by ``observed.impedance_error`` when it
        is available.

    Returns
    -------
    ResponseLossResult
        Reduced response-consistency penalty, ``L_response``.

    Raises
    ------
    TypeError
        If ``forward`` or ``observed`` has the wrong type.
    ValueError
        If station names, components, or frequencies are not
        identical and identically ordered on both inputs.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data import SurveyData
    >>> from pycsamt.forward.maxwell import (
    ...     ForwardResult, SolverDiagnostics
    ... )
    >>> z = np.array([[[1 + 1j]]])
    >>> observed = SurveyData(z, [10.0], ["S1"], ["zxy"], [[0, 0]])
    >>> diagnostics = SolverDiagnostics([[True]], [[1]], [[0.0]], 0.01)
    >>> forward = ForwardResult(
    ...     "a" * 64, [10.0], ["S1"], ["zxy"], z, None,
    ...     "demo", "1", diagnostics,
    ... )
    >>> response_loss_from_contracts(forward, observed).value
    0.0
    """
    from pycsamt.ai.data import SurveyData as _SurveyData
    from pycsamt.forward.maxwell import ForwardResult as _ForwardResult

    if not isinstance(forward, _ForwardResult):
        raise TypeError("forward must be a ForwardResult.")
    if not isinstance(observed, _SurveyData):
        raise TypeError("observed must be a SurveyData.")
    if forward.receiver_names != observed.station_names:
        raise ValueError(
            "forward.receiver_names and observed.station_names must "
            "match exactly; align stations before computing "
            "L_response."
        )
    if forward.components != observed.components:
        raise ValueError(
            "forward.components and observed.components must match "
            "exactly; align components before computing L_response."
        )
    same_frequencies = forward.frequencies_hz.shape == (
        observed.frequencies_hz.shape
    ) and np.array_equal(forward.frequencies_hz, observed.frequencies_hz)
    if not same_frequencies:
        raise ValueError(
            "forward.frequencies_hz and observed.frequencies_hz must "
            "match exactly; use a survey-matched frequency selector "
            "instead of silently interpolating."
        )
    mask = forward.valid & observed.valid
    errors = observed.impedance_error if use_errors else None
    return _response_residual_loss(
        forward.impedance_v_a,
        observed.impedance,
        errors=errors,
        valid=mask,
        kind=kind,
        reduction=reduction,
    )


@dataclass(frozen=True)
class ResponseLoss:
    """Configurable, callable response-consistency loss.

    Parameters
    ----------
    kind : {"l1", "l2"}, default="l2"
        Elementwise penalty applied to each residual magnitude.
    reduction : {"mean", "sum"}, default="mean"
        Reduction applied over included cells.

    Examples
    --------
    >>> import numpy as np
    >>> loss = ResponseLoss()
    >>> pred = np.array([1 + 1j, 2 + 2j])
    >>> obs = np.array([1 + 1j, 0 + 0j])
    >>> loss(pred, obs).value
    4.0
    """

    kind: str = "l2"
    reduction: str = "mean"

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise ValueError(f"kind must be one of {_KINDS}.")
        if self.reduction not in _REDUCTIONS:
            raise ValueError(f"reduction must be one of {_REDUCTIONS}.")

    def __call__(
        self,
        predicted: Any,
        observed: Any,
        *,
        errors: Any | None = None,
        valid: Any | None = None,
    ) -> ResponseLossResult:
        """Evaluate the configured loss on one pair of impedances."""
        return _response_residual_loss(
            predicted,
            observed,
            errors=errors,
            valid=valid,
            kind=self.kind,
            reduction=self.reduction,
        )
