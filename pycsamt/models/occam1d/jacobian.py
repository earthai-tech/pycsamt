# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Analytic and reference Jacobians for native Occam1D inversion.

Sensitivities are taken with respect to layer log10 resistivities. The
analytic implementation differentiates the same impedance recursion used by
:mod:`.forward`; central differences provide an independent verification
path and a documented fallback for future forward-response extensions.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType
from typing import Any

import numpy as np

from .base import Occam1DBase
from .data import Occam1DData
from .forward import Occam1DForwardModel, Occam1DForwardResponse
from .model import Occam1DModel

__all__ = [
    "Occam1DJacobian",
    "Occam1DJacobianCheck",
    "Occam1DJacobianMethod",
    "Occam1DJacobianResult",
]

_LN10 = math.log(10.0)
_RAD_TO_DEG = 180.0 / math.pi
_RHO_CODE = 103
_PHASE_CODE = 104
_THICK_LAYER_LIMIT = 50.0


class Occam1DJacobianMethod(str, Enum):
    """Available sensitivity algorithms."""

    ANALYTIC = "analytic"
    CENTRAL = "central"


def _readonly(values, dtype=float) -> np.ndarray:
    """Return a detached immutable NumPy array."""
    result = np.array(values, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class Occam1DJacobianResult:
    r"""Immutable response derivative and its numerical diagnostics.

    Parameters
    ----------
    matrix : ndarray of shape (n_data, n_layers)
        Derivative of solver-unit predictions with respect to layer
        :math:`\log_{10}` resistivity.
    prediction : ndarray of shape (n_data,)
        Baseline predictions aligned with matrix rows.
    parameters : ndarray of shape (n_layers,)
        Baseline layer :math:`\log_{10}` resistivities.
    frequency_index, type_code : ndarray of shape (n_data,)
        One-based frequency indices and native type codes identifying rows.
    method : {"analytic", "central"}
        Sensitivity algorithm.
    step : ndarray of shape (n_layers,)
        Finite-difference steps. Analytic results contain zeros.
    weighted : bool
        Whether rows and predictions were divided by observational standard
        errors. Raw predictions remain available through an unweighted call.
    metadata : mapping, optional
        Calculation provenance.
    """

    matrix: np.ndarray
    prediction: np.ndarray
    parameters: np.ndarray
    frequency_index: np.ndarray
    type_code: np.ndarray
    method: str
    step: np.ndarray
    weighted: bool = False
    metadata: Any = None

    def __post_init__(self) -> None:
        """Detach arrays and validate shape and finite-value invariants."""
        matrix = _readonly(self.matrix)
        prediction = _readonly(self.prediction)
        parameters = _readonly(self.parameters)
        frequency_index = _readonly(self.frequency_index, int)
        type_code = _readonly(self.type_code, int)
        step = _readonly(self.step)
        if matrix.ndim != 2:
            raise ValueError("Jacobian matrix must be two-dimensional.")
        n_data, n_layers = matrix.shape
        if n_data < 1 or n_layers < 2:
            raise ValueError(
                "Jacobian requires at least one datum and two layers."
            )
        if prediction.shape != (n_data,):
            raise ValueError("prediction must align with Jacobian rows.")
        if frequency_index.shape != (n_data,):
            raise ValueError("frequency_index must align with Jacobian rows.")
        if type_code.shape != (n_data,):
            raise ValueError("type_code must align with Jacobian rows.")
        if parameters.shape != (n_layers,) or step.shape != (n_layers,):
            raise ValueError("parameters and step must align with columns.")
        arrays = (matrix, prediction, parameters, step)
        if any(np.any(~np.isfinite(values)) for values in arrays):
            raise ValueError("Jacobian result arrays must be finite.")
        if np.any(frequency_index < 1):
            raise ValueError("frequency_index values must be one-based.")
        if np.any(~np.isin(type_code, [_RHO_CODE, _PHASE_CODE])):
            raise ValueError("Jacobian rows must use type code 103 or 104.")
        try:
            method = Occam1DJacobianMethod(self.method).value
        except ValueError:
            raise ValueError(
                "method must be 'analytic' or 'central'."
            ) from None
        if not isinstance(self.weighted, bool):
            raise TypeError("weighted must be a bool.")
        metadata = {} if self.metadata is None else dict(self.metadata)
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "prediction", prediction)
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "frequency_index", frequency_index)
        object.__setattr__(self, "type_code", type_code)
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "step", step)
        object.__setattr__(self, "metadata", MappingProxyType(metadata))

    @property
    def shape(self) -> tuple[int, int]:
        """Jacobian matrix shape ``(n_data, n_layers)``."""
        return self.matrix.shape

    @property
    def n_data(self) -> int:
        """Number of scalar response rows."""
        return int(self.matrix.shape[0])

    @property
    def n_layers(self) -> int:
        """Number of log10-resistivity parameters."""
        return int(self.matrix.shape[1])

    @property
    def rank(self) -> int:
        """Numerical matrix rank using NumPy's scale-aware threshold."""
        return int(np.linalg.matrix_rank(self.matrix))

    @property
    def singular_values(self) -> np.ndarray:
        """Read-only singular values in descending order."""
        values = np.linalg.svd(self.matrix, compute_uv=False)
        values.setflags(write=False)
        return values

    @property
    def condition_number(self) -> float:
        """Two-norm condition number, possibly infinity."""
        return float(np.linalg.cond(self.matrix))

    @property
    def column_norm(self) -> np.ndarray:
        """Euclidean sensitivity norm for every model parameter."""
        values = np.linalg.norm(self.matrix, axis=0)
        values.setflags(write=False)
        return values

    def diagnostics(self) -> dict[str, Any]:
        """Return JSON-safe sensitivity and conditioning diagnostics."""
        condition = self.condition_number
        return {
            "method": self.method,
            "weighted": self.weighted,
            "shape": list(self.shape),
            "rank": self.rank,
            "condition_number": (
                condition if math.isfinite(condition) else None
            ),
            "column_norm": self.column_norm.tolist(),
            "step": self.step.tolist(),
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class Occam1DJacobianCheck:
    """Comparison between analytic and central-difference derivatives."""

    analytic: Occam1DJacobianResult
    reference: Occam1DJacobianResult
    absolute_error: np.ndarray
    relative_error: np.ndarray
    atol: float
    rtol: float

    def __post_init__(self) -> None:
        """Validate aligned immutable error arrays and tolerances."""
        if self.analytic.shape != self.reference.shape:
            raise ValueError("Compared Jacobians must have equal shapes.")
        absolute = _readonly(self.absolute_error)
        relative = _readonly(self.relative_error)
        if absolute.shape != self.analytic.shape:
            raise ValueError("absolute_error must match Jacobian shape.")
        if relative.shape != self.analytic.shape:
            raise ValueError("relative_error must match Jacobian shape.")
        object.__setattr__(self, "absolute_error", absolute)
        object.__setattr__(self, "relative_error", relative)

    @property
    def passed(self) -> bool:
        """Whether every element satisfies combined tolerances."""
        threshold = self.atol + self.rtol * np.abs(self.reference.matrix)
        return bool(np.all(self.absolute_error <= threshold))

    @property
    def max_absolute_error(self) -> float:
        """Largest elementwise absolute discrepancy."""
        return float(np.max(self.absolute_error))

    @property
    def max_relative_error(self) -> float:
        """Largest scale-protected elementwise relative discrepancy."""
        return float(np.max(self.relative_error))

    def diagnostics(self) -> dict[str, Any]:
        """Return JSON-safe verification statistics."""
        return {
            "passed": self.passed,
            "shape": list(self.analytic.shape),
            "atol": self.atol,
            "rtol": self.rtol,
            "max_absolute_error": self.max_absolute_error,
            "max_relative_error": self.max_relative_error,
        }


class Occam1DJacobian(Occam1DBase):
    r"""Compute MT response derivatives with respect to log resistivity.

    Parameters
    ----------
    model : Occam1DModel
        Fixed layer geometry and default physical resistivity model.
    forward : Occam1DForwardModel, optional
        Forward engine. When omitted, a new engine is bound to ``model``.
        A supplied engine is not mutated.
    difference_step : float, default=1e-4
        Positive central-difference step in log10 resistivity. This is an
        absolute logarithmic perturbation, so ``1e-4`` changes physical
        resistivity by approximately 0.023 percent.
    verbose, logger, metadata, stream
        Shared :class:`Occam1DBase` options.

    Notes
    -----
    For complex surface impedance :math:`Z` and its parameter derivative,

    .. math::

       \frac{\partial\log_{10}\rho_a}{\partial m_i}
       =\frac{2}{\ln 10}\Re\left(\frac{1}{Z}
       \frac{\partial Z}{\partial m_i}\right),

    .. math::

       \frac{\partial\phi}{\partial m_i}
       =\frac{180}{\pi}\Im\left(\frac{1}{Z}
       \frac{\partial Z}{\partial m_i}\right).

    The analytic recursion propagates all columns simultaneously. It avoids
    one forward solve per perturbation and is therefore the production path.
    """

    def __init__(
        self,
        model,
        *,
        forward=None,
        difference_step=1e-4,
        **kwargs,
    ):
        super().__init__(**kwargs)
        if not isinstance(model, Occam1DModel):
            raise TypeError("model must be an Occam1DModel instance.")
        model.validate()
        if np.any(~np.isfinite(model.resistivity)):
            raise ValueError("model must have finite layer resistivities.")
        if forward is not None and not isinstance(
            forward, Occam1DForwardModel
        ):
            raise TypeError(
                "forward must be an Occam1DForwardModel instance or None."
            )
        self._depth = np.array(model.depth, dtype=float, copy=True)
        self._default_parameters = np.log10(
            np.array(model.resistivity, dtype=float, copy=True)
        )
        self.permeability = (
            forward.permeability if forward is not None else 4e-7 * math.pi
        )
        self.forward = Occam1DForwardModel(
            Occam1DModel(self._depth, 10.0**self._default_parameters),
            permeability=self.permeability,
            logger=self.logger,
        )
        self.difference_step = self._positive_real(
            "difference_step", difference_step
        )

    @staticmethod
    def _positive_real(name, value) -> float:
        """Return a finite strictly positive real scalar."""
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError(f"{name} must be a real number.")
        result = float(value)
        if not math.isfinite(result) or result <= 0:
            raise ValueError(f"{name} must be finite and strictly positive.")
        return result

    @property
    def n_layers(self) -> int:
        """Number of sensitivity parameters."""
        return int(self._depth.size)

    @property
    def default_parameters(self) -> np.ndarray:
        """Independent default log10-resistivity parameter vector."""
        return self._default_parameters.copy()

    def _parameters(self, values) -> np.ndarray:
        """Return a finite parameter vector or the configured default."""
        if values is None:
            return self._default_parameters.copy()
        try:
            result = np.array(values, dtype=float, copy=True)
        except (TypeError, ValueError) as error:
            raise TypeError(
                "parameters must contain real log10 resistivities."
            ) from error
        if result.shape != (self.n_layers,):
            raise ValueError(
                f"parameters must have shape ({self.n_layers},)."
            )
        if np.any(~np.isfinite(result)):
            raise ValueError("parameters must contain only finite values.")
        return result

    @staticmethod
    def _method(value) -> Occam1DJacobianMethod:
        """Return a supported sensitivity method."""
        if isinstance(value, Occam1DJacobianMethod):
            return value
        try:
            return Occam1DJacobianMethod(str(value).strip().lower())
        except ValueError:
            raise ValueError(
                "method must be 'analytic' or 'central'."
            ) from None

    def _analytic_impedance(
        self,
        frequency: np.ndarray,
        parameters: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Propagate impedance and all exact parameter derivatives."""
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            rho = np.power(10.0, parameters)
        if np.any(~np.isfinite(rho)) or np.any(rho <= 0):
            raise ValueError(
                "parameters exceed representable physical resistivity."
            )
        omega = 2.0 * math.pi * frequency
        iwmu = 1j * omega * self.permeability
        n_frequency = frequency.size
        n_layers = self.n_layers
        impedance = np.sqrt(iwmu * rho[-1])
        derivative = np.zeros(
            (n_frequency, n_layers),
            dtype=complex,
        )
        derivative[:, -1] = 0.5 * _LN10 * impedance
        thickness = np.diff(self._depth)

        for layer in range(n_layers - 2, -1, -1):
            intrinsic = np.sqrt(iwmu * rho[layer])
            propagation = np.sqrt(iwmu / rho[layer])
            argument = propagation * thickness[layer]
            thick = argument.real >= _THICK_LAYER_LIMIT
            transfer = np.ones(n_frequency, dtype=complex)
            transfer[~thick] = np.tanh(argument[~thick])

            d_intrinsic = np.zeros_like(derivative)
            d_intrinsic[:, layer] = 0.5 * _LN10 * intrinsic
            d_transfer = np.zeros_like(derivative)
            d_transfer[~thick, layer] = (
                (1.0 - transfer[~thick] ** 2)
                * thickness[layer]
                * (-0.5 * _LN10 * propagation[~thick])
            )

            numerator = impedance + intrinsic * transfer
            denominator = intrinsic + impedance * transfer
            scale = np.maximum(
                np.abs(intrinsic),
                np.abs(impedance * transfer),
            )
            tolerance = np.finfo(float).eps * np.maximum(
                scale,
                np.finfo(float).tiny,
            )
            if np.any(np.abs(denominator) <= tolerance):
                raise FloatingPointError(
                    f"Jacobian recursion is singular at layer {layer + 1}."
                )
            d_numerator = (
                derivative
                + d_intrinsic * transfer[:, None]
                + intrinsic[:, None] * d_transfer
            )
            d_denominator = (
                d_intrinsic
                + derivative * transfer[:, None]
                + impedance[:, None] * d_transfer
            )
            next_impedance = intrinsic * numerator / denominator
            derivative = (
                d_intrinsic * (numerator / denominator)[:, None]
                + intrinsic[:, None]
                * (
                    d_numerator * denominator[:, None]
                    - numerator[:, None] * d_denominator
                )
                / denominator[:, None] ** 2
            )
            impedance = next_impedance
        if np.any(~np.isfinite(derivative)):
            raise FloatingPointError(
                "Analytic recursion produced non-finite sensitivities."
            )
        return impedance, derivative

    @staticmethod
    def _response_derivatives(impedance, derivative):
        """Transform impedance derivatives to log-rho and phase units."""
        ratio = derivative / impedance[:, None]
        log_rho = (2.0 / _LN10) * np.real(ratio)
        phase = _RAD_TO_DEG * np.imag(ratio)
        return log_rho, phase

    @staticmethod
    def _row_layout(data: Occam1DData):
        """Return native row frequency indices and type codes."""
        frequency_index = []
        type_code = []
        resistivity_mask = data.resistivity_mask
        phase_mask = data.phase_mask
        for index in range(data.n_frequencies):
            if resistivity_mask[index]:
                frequency_index.append(index + 1)
                type_code.append(_RHO_CODE)
            if phase_mask[index]:
                frequency_index.append(index + 1)
                type_code.append(_PHASE_CODE)
        return (
            np.asarray(frequency_index, dtype=int),
            np.asarray(type_code, dtype=int),
        )

    @staticmethod
    def _align_rows(log_rho, phase, frequency_index, type_code):
        """Align frequency response arrays with native observation rows."""
        rows = np.empty((type_code.size, log_rho.shape[1]), dtype=float)
        for row, (index, code) in enumerate(zip(frequency_index, type_code)):
            source = log_rho if code == _RHO_CODE else phase
            rows[row] = source[index - 1]
        return rows

    @staticmethod
    def _prediction(response, frequency_index, type_code):
        """Align one forward response with native observation rows."""
        values = np.empty(type_code.size, dtype=float)
        log_rho = response.log10_resistivity
        for row, (index, code) in enumerate(zip(frequency_index, type_code)):
            values[row] = (
                log_rho[index - 1]
                if code == _RHO_CODE
                else response.phase[index - 1]
            )
        return values

    @staticmethod
    def solver_errors(data: Occam1DData) -> np.ndarray:
        """Return standard errors aligned in native solver units.

        Resistivity uncertainty is converted from relative physical error to
        log10 units using ``sigma_log10 = relative_error / ln(10)``. Phase
        uncertainty remains in degrees.
        """
        if not isinstance(data, Occam1DData):
            raise TypeError("data must be an Occam1DData instance.")
        errors = []
        resistivity_mask = data.resistivity_mask
        phase_mask = data.phase_mask
        for index in range(data.n_frequencies):
            if resistivity_mask[index]:
                errors.append(data.resistivity_error[index] / _LN10)
            if phase_mask[index]:
                errors.append(data.phase_error[index])
        result = np.asarray(errors, dtype=float)
        if result.shape != (data.n_data,):
            raise RuntimeError("Could not align Occam1D solver errors.")
        if np.any(~np.isfinite(result)) or np.any(result <= 0):
            raise ValueError("Solver errors must be finite and positive.")
        return result

    def _column_delta(self, data, parameters, step, frequency_index, type_code, column):
        """Return one central-difference Jacobian column."""
        upper = parameters.copy()
        lower = parameters.copy()
        upper[column] += step[column]
        lower[column] -= step[column]
        plus = self.forward.predict_data(
            data,
            log10_resistivity=upper,
        )
        minus = self.forward.predict_data(
            data,
            log10_resistivity=lower,
        )
        plus_values = self._prediction(
            plus,
            frequency_index,
            type_code,
        )
        minus_values = self._prediction(
            minus,
            frequency_index,
            type_code,
        )
        delta = plus_values - minus_values
        phase_rows = type_code == _PHASE_CODE
        delta[phase_rows] = (
            delta[phase_rows] + 90.0
        ) % 180.0 - 90.0
        return delta / (2.0 * step[column])

    def _central(
        self,
        data,
        parameters,
        step,
        n_jobs=1,
    ) -> tuple[np.ndarray, np.ndarray, Occam1DForwardResponse]:
        """Return central differences using independent forward solves.

        Every column perturbs one layer and re-runs the forward model twice;
        columns do not depend on each other, so ``n_jobs != 1`` evaluates
        them with a thread pool via :mod:`joblib` (soft dependency; falls
        back to a plain sequential loop when joblib is not installed).
        Threads, not processes, because each column is only a couple of
        cheap forward solves -- too small to amortize process/pickling
        overhead, and NumPy's own C calls release the GIL during that time.
        """
        baseline = self.forward.predict_data(
            data,
            log10_resistivity=parameters,
        )
        frequency_index, type_code = self._row_layout(data)
        columns = None
        if n_jobs != 1:
            try:
                from joblib import Parallel, delayed
            except ImportError:
                columns = None
            else:
                columns = Parallel(n_jobs=n_jobs, prefer="threads")(
                    delayed(self._column_delta)(
                        data, parameters, step, frequency_index, type_code, column
                    )
                    for column in range(self.n_layers)
                )
        if columns is None:
            columns = [
                self._column_delta(
                    data, parameters, step, frequency_index, type_code, column
                )
                for column in range(self.n_layers)
            ]
        matrix = np.column_stack(columns)
        return matrix, step, baseline

    def compute(
        self,
        data,
        *,
        parameters=None,
        method="analytic",
        difference_step=None,
        weighted=False,
        n_jobs=1,
    ) -> Occam1DJacobianResult:
        """Compute a response Jacobian aligned with observed data rows.

        Parameters
        ----------
        data : Occam1DData
            Observational layout controlling frequencies, mode, and which
            resistivity/phase rows are included.
        parameters : array-like of shape (n_layers,), optional
            Baseline log10 layer resistivities.
        method : {"analytic", "central"}, default="analytic"
            Exact recursive differentiation or finite-difference reference.
        difference_step : float or array-like, optional
            Positive log10 perturbation for central differences. A scalar is
            broadcast to every layer.
        weighted : bool, default=False
            Divide each matrix row and baseline prediction by its native
            observational standard error.
        n_jobs : int, default=1
            Worker threads for ``method="central"``, one column per
            perturbation (see :meth:`_central`). Ignored for the default
            ``"analytic"`` method, which already evaluates every column in
            a single simultaneous recursion and has nothing to parallelize.

        Returns
        -------
        Occam1DJacobianResult
            Immutable sensitivity matrix and diagnostics.
        """
        if not isinstance(data, Occam1DData):
            raise TypeError("data must be an Occam1DData instance.")
        data.validate()
        if not isinstance(weighted, bool):
            raise TypeError("weighted must be a bool.")
        parameters = self._parameters(parameters)
        selected_method = self._method(method)
        frequency_index, type_code = self._row_layout(data)

        if selected_method is Occam1DJacobianMethod.ANALYTIC:
            impedance, d_impedance = self._analytic_impedance(
                data.frequency,
                parameters,
            )
            log_rho, phase = self._response_derivatives(
                impedance,
                d_impedance,
            )
            matrix = self._align_rows(
                log_rho,
                phase,
                frequency_index,
                type_code,
            )
            baseline = self.forward.predict_data(
                data,
                log10_resistivity=parameters,
            )
            step = np.zeros(self.n_layers, dtype=float)
        else:
            step = self._steps(difference_step)
            matrix, step, baseline = self._central(
                data,
                parameters,
                step,
                n_jobs=n_jobs,
            )
        prediction = self._prediction(
            baseline,
            frequency_index,
            type_code,
        )
        if weighted:
            errors = self.solver_errors(data)
            matrix = matrix / errors[:, None]
            prediction = prediction / errors
        if np.any(~np.isfinite(matrix)):
            raise FloatingPointError(
                "Jacobian calculation produced non-finite derivatives."
            )
        return Occam1DJacobianResult(
            matrix=matrix,
            prediction=prediction,
            parameters=parameters,
            frequency_index=frequency_index,
            type_code=type_code,
            method=selected_method.value,
            step=step,
            weighted=weighted,
            metadata={
                "mode": data.mode,
                "station": data.station,
                "permeability_h_per_m": self.permeability,
                "parameterization": "log10_resistivity",
            },
        )

    def _steps(self, value) -> np.ndarray:
        """Return positive finite per-parameter difference steps."""
        if value is None:
            return np.full(self.n_layers, self.difference_step, dtype=float)
        if isinstance(value, bool):
            raise TypeError("difference_step must contain real numbers.")
        if isinstance(value, (int, float)):
            result = np.full(self.n_layers, float(value), dtype=float)
        else:
            try:
                result = np.asarray(value, dtype=float)
            except (TypeError, ValueError) as error:
                raise TypeError(
                    "difference_step must contain real numbers."
                ) from error
            if result.shape != (self.n_layers,):
                raise ValueError(
                    f"difference_step must have shape ({self.n_layers},)."
                )
            result = result.copy()
        if np.any(~np.isfinite(result)) or np.any(result <= 0):
            raise ValueError(
                "difference_step values must be finite and positive."
            )
        return result

    def check(
        self,
        data,
        *,
        parameters=None,
        difference_step=None,
        atol=1e-6,
        rtol=2e-4,
        n_jobs=1,
    ) -> Occam1DJacobianCheck:
        """Compare the analytic Jacobian with central differences.

        This method is intended for tests, solver diagnostics, and new
        response conventions. Inversion iterations should call the analytic
        method directly to avoid ``2 * n_layers`` extra forward evaluations.
        ``n_jobs`` speeds up the reference central-difference evaluation
        only; see :meth:`compute`.
        """
        atol = self._positive_real("atol", atol)
        rtol = self._positive_real("rtol", rtol)
        analytic = self.compute(
            data,
            parameters=parameters,
            method="analytic",
        )
        reference = self.compute(
            data,
            parameters=parameters,
            method="central",
            difference_step=difference_step,
            n_jobs=n_jobs,
        )
        absolute = np.abs(analytic.matrix - reference.matrix)
        scale = np.maximum(np.abs(reference.matrix), atol)
        relative = absolute / scale
        return Occam1DJacobianCheck(
            analytic=analytic,
            reference=reference,
            absolute_error=absolute,
            relative_error=relative,
            atol=atol,
            rtol=rtol,
        )

    def diagnostics(self) -> dict[str, Any]:
        """Extend lifecycle diagnostics with sensitivity configuration."""
        values = super().diagnostics()
        values.update(
            {
                "n_layers": self.n_layers,
                "difference_step": self.difference_step,
                "permeability_h_per_m": self.permeability,
                "parameterization": "log10_resistivity",
                "methods": [value.value for value in Occam1DJacobianMethod],
            }
        )
        return values
