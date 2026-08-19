# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Native nonlinear Occam inversion for one-dimensional MT soundings.

This module joins the forward, sensitivity, and regularization layers.  The
solver works in absolute base-10 log resistivity, evaluates every trial with
the nonlinear forward model, and selects the smoothest model that satisfies
the requested normalized RMS. Native input/output remains the responsibility
of the builder and runner; this module writes only explicit JSON restart
checkpoints.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from ...api.view.progress import get_progress_bar
from .base import Occam1DBase
from .config import Occam1DConfig
from .data import Occam1DData
from .forward import Occam1DForwardModel
from .jacobian import Occam1DJacobian
from .model import Occam1DModel
from .regularization import Occam1DRegularization, Occam1DSolverPolicy
from .startup import Occam1DStartup

__all__ = [
    "Occam1DConvergence",
    "Occam1DFailedIteration",
    "Occam1DInversion",
    "Occam1DInversionResult",
    "Occam1DIteration",
    "Occam1DRejectedCandidate",
    "Occam1DRejectionReason",
    "Occam1DRestart",
]


def _readonly(values, dtype=float) -> np.ndarray:
    """Return an immutable detached array."""
    result = np.array(values, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def _data_fingerprint(data: Occam1DData) -> str:
    """Return a stable digest of solver-relevant sounding content."""
    digest = hashlib.sha256()
    digest.update(data.mode.encode("utf8"))
    for values in (
        data.frequency,
        data.resistivity,
        data.phase,
        data.resistivity_error,
        data.phase_error,
    ):
        canonical = np.asarray(values, dtype="<f8")
        digest.update(canonical.tobytes())
    return digest.hexdigest()


class Occam1DConvergence(str, Enum):
    """Terminal state of a native inversion.

    ``TARGET`` means the normalized RMS reached the requested target.
    ``STAGNATED`` means successive accepted models made insufficient change.
    ``MAX_ITERATIONS`` means the iteration budget was exhausted, while
    ``CANCELLED`` records a cooperative callback cancellation.
    """

    TARGET = "target"
    STAGNATED = "stagnated"
    MAX_ITERATIONS = "max_iterations"
    CANCELLED = "cancelled"
    FAILED = "failed"


class Occam1DRejectionReason(str, Enum):
    """Machine-readable reason an inversion candidate was rejected."""

    LINEAR_SOLVE_FAILED = "linear_solve_failed"
    FORWARD_EVALUATION_FAILED = "forward_evaluation_failed"
    STEP_CUT_DISCARDED = "step_cut_discarded"
    NOT_SELECTED = "not_selected"
    RMS_REGRESSION = "rms_regression"


@dataclass(frozen=True)
class Occam1DRejectedCandidate:
    """Structured audit record for one rejected candidate or step cut."""

    iteration: int
    multiplier: float
    reason: Occam1DRejectionReason
    message: str
    step_scale: float | None = None
    rms: float | None = None
    roughness: float | None = None
    exception_type: str | None = None

    def __post_init__(self) -> None:
        """Validate rejection identity, metrics, and explanation."""
        if (
            isinstance(self.iteration, bool)
            or not isinstance(self.iteration, int)
            or self.iteration < 1
        ):
            raise ValueError("rejection iteration must be a positive integer.")
        if not math.isfinite(self.multiplier) or self.multiplier < 0:
            raise ValueError(
                "rejection multiplier must be finite and non-negative."
            )
        if not isinstance(self.reason, Occam1DRejectionReason):
            raise TypeError("reason must be Occam1DRejectionReason.")
        if not isinstance(self.message, str) or not self.message:
            raise ValueError("rejection message must be non-empty.")
        for name in ("step_scale", "rms", "roughness"):
            value = getattr(self, name)
            if value is not None and (
                not math.isfinite(value) or value < 0
            ):
                raise ValueError(
                    f"rejection {name} must be finite and non-negative."
                )
        if self.exception_type is not None and (
            not isinstance(self.exception_type, str)
            or not self.exception_type
        ):
            raise ValueError("exception_type must be non-empty or None.")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe rejection record."""
        return {
            "iteration": self.iteration,
            "multiplier": self.multiplier,
            "reason": self.reason.value,
            "message": self.message,
            "step_scale": self.step_scale,
            "rms": self.rms,
            "roughness": self.roughness,
            "exception_type": self.exception_type,
        }

    @classmethod
    def from_dict(cls, values):
        """Restore a rejection record from JSON-safe values."""
        return cls(
            iteration=values["iteration"],
            multiplier=values["multiplier"],
            reason=Occam1DRejectionReason(values["reason"]),
            message=values["message"],
            step_scale=values.get("step_scale"),
            rms=values.get("rms"),
            roughness=values.get("roughness"),
            exception_type=values.get("exception_type"),
        )


@dataclass(frozen=True)
class Occam1DFailedIteration:
    """Structured record of an iteration that accepted no model update."""

    iteration: int
    current_rms: float
    attempted_candidates: int
    rejected_candidates: int
    reason: str
    message: str
    recoverable: bool

    def __post_init__(self) -> None:
        """Validate failed-iteration counters and terminal context."""
        if (
            isinstance(self.iteration, bool)
            or not isinstance(self.iteration, int)
            or self.iteration < 1
        ):
            raise ValueError("failed iteration must be a positive integer.")
        if not math.isfinite(self.current_rms) or self.current_rms < 0:
            raise ValueError("current_rms must be finite and non-negative.")
        for name in ("attempted_candidates", "rejected_candidates"):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, int)
                or value < 0
            ):
                raise ValueError(f"{name} must be a non-negative integer.")
        if not isinstance(self.reason, str) or not self.reason:
            raise ValueError("failure reason must be non-empty.")
        if not isinstance(self.message, str) or not self.message:
            raise ValueError("failure message must be non-empty.")
        if not isinstance(self.recoverable, bool):
            raise TypeError("recoverable must be a bool.")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe failed-iteration record."""
        return {
            "iteration": self.iteration,
            "current_rms": self.current_rms,
            "attempted_candidates": self.attempted_candidates,
            "rejected_candidates": self.rejected_candidates,
            "reason": self.reason,
            "message": self.message,
            "recoverable": self.recoverable,
        }

    @classmethod
    def from_dict(cls, values):
        """Restore a failed-iteration record."""
        return cls(**values)


@dataclass(frozen=True)
class Occam1DIteration:
    r"""Immutable record of one accepted nonlinear model.

    Parameters are absolute :math:`\log_{10}` resistivities. ``rms`` is the
    square root of the mean squared error-normalized data residual.
    ``roughness`` includes configured smoothness and preference penalties.
    """

    number: int
    parameters: np.ndarray
    prediction: np.ndarray
    residual: np.ndarray
    rms: float
    roughness: float
    multiplier: float
    step_scale: float
    rank: int
    condition_number: float
    active_lower: np.ndarray | None = None
    active_upper: np.ndarray | None = None
    solver: str = "unknown"
    solve_attempts: int = 0
    damping: float = 0.0
    initial_condition_number: float | None = None

    def __post_init__(self) -> None:
        """Validate record invariants and freeze numerical arrays."""
        if (
            isinstance(self.number, bool)
            or not isinstance(self.number, int)
            or self.number < 0
        ):
            raise ValueError("number must be a non-negative integer.")
        parameters = _readonly(self.parameters)
        prediction = _readonly(self.prediction)
        residual = _readonly(self.residual)
        active_lower = (
            np.zeros(parameters.size, dtype=bool)
            if self.active_lower is None
            else _readonly(self.active_lower, dtype=bool)
        )
        active_upper = (
            np.zeros(parameters.size, dtype=bool)
            if self.active_upper is None
            else _readonly(self.active_upper, dtype=bool)
        )
        if parameters.ndim != 1 or prediction.ndim != 1:
            raise ValueError("parameters and prediction must be vectors.")
        if residual.shape != prediction.shape:
            raise ValueError("residual must align with prediction.")
        if active_lower.shape != parameters.shape:
            raise ValueError("active_lower must align with parameters.")
        if active_upper.shape != parameters.shape:
            raise ValueError("active_upper must align with parameters.")
        if any(np.any(~np.isfinite(value)) for value in (
            parameters,
            prediction,
            residual,
        )):
            raise ValueError("Iteration arrays must contain finite values.")
        for name in ("rms", "roughness", "multiplier", "step_scale"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and non-negative.")
        if isinstance(self.rank, bool) or not isinstance(self.rank, int):
            raise TypeError("rank must be an integer.")
        if self.rank < 0:
            raise ValueError("rank must be non-negative.")
        if math.isnan(self.condition_number) or self.condition_number < 0:
            raise ValueError(
                "condition_number must be non-negative and not NaN."
            )
        if not isinstance(self.solver, str) or not self.solver:
            raise ValueError("solver must be a non-empty string.")
        if self.solve_attempts < 0:
            raise ValueError("solve_attempts must be non-negative.")
        if not math.isfinite(self.damping) or self.damping < 0:
            raise ValueError("damping must be finite and non-negative.")
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "prediction", prediction)
        object.__setattr__(self, "residual", residual)
        object.__setattr__(self, "active_lower", active_lower)
        object.__setattr__(self, "active_upper", active_upper)

    @property
    def resistivity(self) -> np.ndarray:
        """Physical layer resistivity in ohm metres."""
        return np.power(10.0, self.parameters)

    @property
    def n_active_bounds(self) -> int:
        """Number of layer parameters at either configured bound."""
        return int(np.count_nonzero(self.active_lower | self.active_upper))

    def diagnostics(self) -> dict[str, Any]:
        """Return a JSON-safe iteration summary."""
        return {
            "iteration": self.number,
            "rms": self.rms,
            "roughness": self.roughness,
            "multiplier": self.multiplier,
            "step_scale": self.step_scale,
            "rank": self.rank,
            "condition_number": (
                self.condition_number
                if math.isfinite(self.condition_number)
                else None
            ),
            "n_active_bounds": self.n_active_bounds,
            "active_lower": self.active_lower.tolist(),
            "active_upper": self.active_upper.tolist(),
            "solver": self.solver,
            "solve_attempts": self.solve_attempts,
            "damping": self.damping,
            "initial_condition_number": (
                self.initial_condition_number
                if self.initial_condition_number is not None
                and math.isfinite(self.initial_condition_number)
                else None
            ),
        }

    def to_dict(self) -> dict[str, Any]:
        """Return the complete JSON-safe accepted-model record."""
        values = self.diagnostics()
        values.update(
            {
                "parameters": self.parameters.tolist(),
                "prediction": self.prediction.tolist(),
                "residual": self.residual.tolist(),
            }
        )
        return values

    @classmethod
    def from_dict(cls, values):
        """Restore an iteration from :meth:`to_dict` output."""
        if not isinstance(values, dict):
            raise TypeError("iteration payload must be a dictionary.")
        required = {
            "iteration",
            "parameters",
            "prediction",
            "residual",
            "rms",
            "roughness",
            "multiplier",
            "step_scale",
            "rank",
            "condition_number",
        }
        missing = sorted(required.difference(values))
        if missing:
            raise ValueError(
                "Iteration payload is missing: " + ", ".join(missing)
            )
        return cls(
            number=values["iteration"],
            parameters=values["parameters"],
            prediction=values["prediction"],
            residual=values["residual"],
            rms=values["rms"],
            roughness=values["roughness"],
            multiplier=values["multiplier"],
            step_scale=values["step_scale"],
            rank=values["rank"],
            condition_number=(
                float("inf")
                if values["condition_number"] is None
                else values["condition_number"]
            ),
            active_lower=values.get("active_lower"),
            active_upper=values.get("active_upper"),
            solver=values.get("solver", "unknown"),
            solve_attempts=values.get("solve_attempts", 0),
            damping=values.get("damping", 0.0),
            initial_condition_number=values.get(
                "initial_condition_number"
            ),
        )


@dataclass(frozen=True)
class Occam1DRestart:
    r"""Immutable, serializable checkpoint for continuing an inversion.

    Parameters
    ----------
    iterations : sequence of Occam1DIteration
        Complete contiguous accepted history beginning at iteration zero.
    depth : array-like of shape (n_layers,)
        Layer-top depths in metres used by the inversion.
    target_rms : float
        Target normalized RMS associated with the original run.
    data_fingerprint : str, optional
        SHA-256 digest of solver-relevant data. A checkpoint created through
        :meth:`Occam1DInversion.restart` always includes this digest.

    Notes
    -----
    The active model, iteration number, and multiplier are authoritative
    properties of the final accepted record. They are not stored twice, so a
    checkpoint cannot contain contradictory continuation state.
    """

    iterations: tuple[Occam1DIteration, ...]
    depth: np.ndarray
    target_rms: float
    data_fingerprint: str | None = None
    format_version: int = 1
    previous_convergence: Occam1DConvergence | None = None
    previous_message: str | None = None
    rejected_candidates: tuple[Occam1DRejectedCandidate, ...] = ()
    failed_iterations: tuple[Occam1DFailedIteration, ...] = ()

    def __post_init__(self) -> None:
        """Validate history, geometry, and checkpoint provenance."""
        history = tuple(self.iterations)
        if not history:
            raise ValueError("restart history cannot be empty.")
        if not all(isinstance(item, Occam1DIteration) for item in history):
            raise TypeError(
                "restart history must contain Occam1DIteration objects."
            )
        numbers = tuple(item.number for item in history)
        if numbers != tuple(range(len(history))):
            raise ValueError(
                "Restart iteration numbers must be contiguous from zero."
            )
        depth = _readonly(self.depth)
        if depth.shape != history[-1].parameters.shape:
            raise ValueError("restart depth must align with model parameters.")
        if np.any(~np.isfinite(depth)) or np.any(np.diff(depth) <= 0):
            raise ValueError(
                "restart depth must be finite and strictly increasing."
            )
        n_data = history[0].prediction.size
        if any(item.prediction.size != n_data for item in history):
            raise ValueError(
                "All restart predictions must have the same data count."
            )
        if not math.isfinite(self.target_rms) or self.target_rms <= 0:
            raise ValueError("restart target_rms must be finite and positive.")
        if self.data_fingerprint is not None:
            if (
                not isinstance(self.data_fingerprint, str)
                or len(self.data_fingerprint) != 64
            ):
                raise ValueError(
                    "data_fingerprint must be a SHA-256 hexadecimal digest."
                )
            try:
                bytes.fromhex(self.data_fingerprint)
            except ValueError:
                raise ValueError(
                    "data_fingerprint must be hexadecimal."
                ) from None
        if self.format_version != 1:
            raise ValueError("Unsupported Occam1D restart format version.")
        if self.previous_convergence is not None and not isinstance(
            self.previous_convergence,
            Occam1DConvergence,
        ):
            raise TypeError(
                "previous_convergence must be Occam1DConvergence or None."
            )
        if self.previous_message is not None and (
            not isinstance(self.previous_message, str)
            or not self.previous_message
        ):
            raise ValueError(
                "previous_message must be a non-empty string or None."
            )
        rejected = tuple(self.rejected_candidates)
        failures = tuple(self.failed_iterations)
        if not all(
            isinstance(item, Occam1DRejectedCandidate) for item in rejected
        ):
            raise TypeError("Invalid restart rejected-candidate history.")
        if not all(
            isinstance(item, Occam1DFailedIteration) for item in failures
        ):
            raise TypeError("Invalid restart failed-iteration history.")
        object.__setattr__(self, "iterations", history)
        object.__setattr__(self, "depth", depth)
        object.__setattr__(self, "rejected_candidates", rejected)
        object.__setattr__(self, "failed_iterations", failures)

    @property
    def iteration_number(self) -> int:
        """Absolute number of the latest accepted iteration."""
        return self.iterations[-1].number

    @property
    def parameters(self) -> np.ndarray:
        """Independent active log10-resistivity model."""
        return self.iterations[-1].parameters.copy()

    @property
    def resistivity(self) -> np.ndarray:
        """Independent active physical model in ohm metres."""
        return self.iterations[-1].resistivity

    @property
    def multiplier(self) -> float:
        """Lagrange multiplier selected at the latest iteration."""
        return self.iterations[-1].multiplier

    @property
    def rms_history(self) -> np.ndarray:
        """Normalized RMS history through the checkpoint."""
        return np.asarray([item.rms for item in self.iterations])

    @property
    def roughness_history(self) -> np.ndarray:
        """Structure-objective history through the checkpoint."""
        return np.asarray([item.roughness for item in self.iterations])

    def to_dict(self) -> dict[str, Any]:
        """Return a complete JSON-safe checkpoint payload."""
        return {
            "format": "pycsamt_occam1d_restart",
            "format_version": self.format_version,
            "target_rms": self.target_rms,
            "depth_m": self.depth.tolist(),
            "data_fingerprint": self.data_fingerprint,
            "previous_convergence": (
                None
                if self.previous_convergence is None
                else self.previous_convergence.value
            ),
            "previous_message": self.previous_message,
            "iterations": [item.to_dict() for item in self.iterations],
            "rejected_candidates": [
                item.to_dict() for item in self.rejected_candidates
            ],
            "failed_iterations": [
                item.to_dict() for item in self.failed_iterations
            ],
        }

    @classmethod
    def from_dict(cls, values):
        """Restore and validate a checkpoint dictionary."""
        if not isinstance(values, dict):
            raise TypeError("restart payload must be a dictionary.")
        if values.get("format") != "pycsamt_occam1d_restart":
            raise ValueError("Not a pyCSAMT Occam1D restart payload.")
        try:
            iterations = tuple(
                Occam1DIteration.from_dict(item)
                for item in values["iterations"]
            )
            depth = values["depth_m"]
            target = values["target_rms"]
        except KeyError as error:
            raise ValueError(
                f"Restart payload is missing {error.args[0]!r}."
            ) from None
        return cls(
            iterations=iterations,
            depth=depth,
            target_rms=target,
            data_fingerprint=values.get("data_fingerprint"),
            format_version=values.get("format_version", 1),
            previous_convergence=(
                None
                if values.get("previous_convergence") is None
                else Occam1DConvergence(values["previous_convergence"])
            ),
            previous_message=values.get("previous_message"),
            rejected_candidates=tuple(
                Occam1DRejectedCandidate.from_dict(item)
                for item in values.get("rejected_candidates", ())
            ),
            failed_iterations=tuple(
                Occam1DFailedIteration.from_dict(item)
                for item in values.get("failed_iterations", ())
            ),
        )

    def diagnostics(self) -> dict[str, Any]:
        """Return JSON-safe checkpoint state and convergence traces."""
        return {
            "iteration": self.iteration_number,
            "multiplier": self.multiplier,
            "target_rms": self.target_rms,
            "n_layers": int(self.depth.size),
            "n_data": int(self.iterations[-1].prediction.size),
            "previous_convergence": (
                None
                if self.previous_convergence is None
                else self.previous_convergence.value
            ),
            "previous_message": self.previous_message,
            "rms_history": self.rms_history.tolist(),
            "roughness_history": self.roughness_history.tolist(),
            "has_data_fingerprint": self.data_fingerprint is not None,
            "n_rejected_candidates": len(self.rejected_candidates),
            "n_failed_iterations": len(self.failed_iterations),
        }

    def write(self, path) -> Path:
        """Atomically write this checkpoint as UTF-8 JSON."""
        target = Path(path).expanduser().resolve()
        if target.exists() and target.is_dir():
            raise IsADirectoryError(f"Restart path is a directory: {target}")
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf8",
                newline="\n",
                prefix=f".{target.name}.",
                suffix=".tmp",
                dir=target.parent,
                delete=False,
            ) as stream:
                temporary = Path(stream.name)
                json.dump(
                    self.to_dict(),
                    stream,
                    indent=2,
                    sort_keys=True,
                    allow_nan=False,
                )
                stream.write("\n")
            temporary.replace(target)
        except Exception:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
            raise
        return target

    @classmethod
    def read(cls, path):
        """Read and validate a JSON checkpoint from disk."""
        source = Path(path).expanduser().resolve()
        if not source.is_file():
            raise FileNotFoundError(f"Restart file does not exist: {source}")
        try:
            values = json.loads(source.read_text(encoding="utf8"))
        except json.JSONDecodeError as error:
            raise ValueError(
                f"Malformed Occam1D restart JSON in {source}: {error}"
            ) from error
        return cls.from_dict(values)


@dataclass(frozen=True)
class Occam1DInversionResult:
    """Complete immutable in-memory outcome of one inversion."""

    iterations: tuple[Occam1DIteration, ...]
    convergence: Occam1DConvergence
    target_rms: float
    message: str
    rejected_candidates: tuple[Occam1DRejectedCandidate, ...] = ()
    failed_iterations: tuple[Occam1DFailedIteration, ...] = ()

    def __post_init__(self) -> None:
        """Validate history ordering and terminal metadata."""
        history = tuple(self.iterations)
        if not history:
            raise ValueError("iterations must include the initial model.")
        expected = tuple(range(len(history)))
        if tuple(item.number for item in history) != expected:
            raise ValueError("Iteration numbers must be contiguous from zero.")
        if not isinstance(self.convergence, Occam1DConvergence):
            raise TypeError("convergence must be Occam1DConvergence.")
        if not math.isfinite(self.target_rms) or self.target_rms <= 0:
            raise ValueError("target_rms must be finite and positive.")
        if not isinstance(self.message, str) or not self.message:
            raise ValueError("message must be a non-empty string.")
        rejected = tuple(self.rejected_candidates)
        failures = tuple(self.failed_iterations)
        if not all(
            isinstance(item, Occam1DRejectedCandidate) for item in rejected
        ):
            raise TypeError(
                "rejected_candidates must contain rejection records."
            )
        if not all(
            isinstance(item, Occam1DFailedIteration) for item in failures
        ):
            raise TypeError(
                "failed_iterations must contain failure records."
            )
        object.__setattr__(self, "iterations", history)
        object.__setattr__(self, "rejected_candidates", rejected)
        object.__setattr__(self, "failed_iterations", failures)

    @property
    def initial(self) -> Occam1DIteration:
        """Initial-model evaluation."""
        return self.iterations[0]

    @property
    def final(self) -> Occam1DIteration:
        """Last accepted model."""
        return self.iterations[-1]

    @property
    def converged(self) -> bool:
        """Whether the target RMS was reached."""
        return self.convergence is Occam1DConvergence.TARGET

    @property
    def n_iterations(self) -> int:
        """Number of nonlinear updates, excluding iteration zero."""
        return len(self.iterations) - 1

    @property
    def rms_history(self) -> np.ndarray:
        """Independent normalized-RMS convergence history."""
        return np.asarray([item.rms for item in self.iterations])

    def diagnostics(self) -> dict[str, Any]:
        """Return JSON-safe terminal and iteration diagnostics."""
        return {
            "convergence": self.convergence.value,
            "converged": self.converged,
            "message": self.message,
            "target_rms": self.target_rms,
            "n_iterations": self.n_iterations,
            "initial_rms": self.initial.rms,
            "final_rms": self.final.rms,
            "final_roughness": self.final.roughness,
            "n_rejected_candidates": len(self.rejected_candidates),
            "n_failed_iterations": len(self.failed_iterations),
            "history": [item.diagnostics() for item in self.iterations],
            "rejected_candidates": [
                item.to_dict() for item in self.rejected_candidates
            ],
            "failed_iterations": [
                item.to_dict() for item in self.failed_iterations
            ],
        }


class Occam1DInversion(Occam1DBase):
    r"""Solve one layered-earth MT sounding by nonlinear Occam inversion.

    Parameters
    ----------
    data : Occam1DData
        Observed apparent resistivity in ohm metres, phase in degrees, and
        positive observational errors.
    model : Occam1DModel
        Fixed layer-top depths in metres and starting resistivities.
    config : Occam1DConfig, optional
        Target RMS, iteration limit, roughness order, and native controls.
    startup : Occam1DStartup, optional
        Starting log10 parameters and solver controls. When supplied these
        controls take precedence over ``config``.
    regularization : Occam1DRegularization, optional
        Explicit structure objective. The default uses the model geometry
        and configured paper-compatible roughness order.
    multiplier_factors : sequence of float, optional
        Positive factors applied to the current Lagrange multiplier. Defaults
        to 13 logarithmically spaced values from ``1e-3`` through ``1e3``.
    log_resistivity_bounds : pair, optional
        Lower and upper bounds for absolute base-10 log resistivity. Each
        side may be a scalar, a vector of shape ``(n_layers,)``, or ``None``
        for an open side. For example, ``(0, 6)`` restricts physical
        resistivity to 1 through 1,000,000 ohm metres in every layer.
    lagrange_bounds : pair of float, optional
        Absolute inclusive minimum and maximum Lagrange multipliers. When
        omitted, ``config.lagrange_min`` and ``config.lagrange_max`` are
        authoritative. Products of the current multiplier and search factors
        are saturated to this interval and deduplicated before solving.
    solver_policy : Occam1DSolverPolicy, optional
        Explicit SVD fallback, conditioning threshold, and diagonal-damping
        retry policy. The validated default retries ill-conditioned systems
        three times and permits an independent SciPy LAPACK fallback.
    rms_tolerance : float, default=0.01
        Relative tolerance above target RMS for an acceptable model.
    model_tolerance : float, default=1e-4
        Stagnation threshold for maximum absolute log10 parameter change.
    rms_change_tolerance : float, default=1e-4
        Stagnation threshold for relative RMS improvement.
    verbose, logger, metadata, stream
        Shared pyCSAMT diagnostics and user-output controls.

    Notes
    -----
    Candidate linear systems are error weighted. Every candidate is then
    passed through the nonlinear forward model; linearized objective values
    are never reported as achieved data fits. Phase residuals are wrapped to
    the nearest 180-degree branch.
    """

    def __init__(
        self,
        data,
        model,
        *,
        config=None,
        startup=None,
        regularization=None,
        multiplier_factors=None,
        log_resistivity_bounds=None,
        lagrange_bounds=None,
        solver_policy=None,
        rms_tolerance=0.01,
        model_tolerance=1e-4,
        rms_change_tolerance=1e-4,
        **kwargs,
    ):
        super().__init__(**kwargs)
        if not isinstance(data, Occam1DData):
            raise TypeError("data must be an Occam1DData instance.")
        if not isinstance(model, Occam1DModel):
            raise TypeError("model must be an Occam1DModel instance.")
        data.validate()
        model.validate()
        self.data = data
        self.model = model
        self.config = config or Occam1DConfig(n_layers=model.n_layers)
        if not isinstance(self.config, Occam1DConfig):
            raise TypeError("config must be an Occam1DConfig instance.")
        self.config.validate()
        if self.config.n_layers != model.n_layers:
            raise ValueError("config and model layer counts do not match.")
        self.startup = startup or Occam1DStartup.from_model(
            model,
            self.config,
        )
        if not isinstance(self.startup, Occam1DStartup):
            raise TypeError("startup must be an Occam1DStartup instance.")
        if self.startup.n_parameters != model.n_layers:
            raise ValueError("startup and model layer counts do not match.")
        self.regularization = regularization or Occam1DRegularization(
            model.n_layers,
            order=self.startup.roughness_type,
            depth=model.depth,
            logger=self.logger,
        )
        if not isinstance(self.regularization, Occam1DRegularization):
            raise TypeError(
                "regularization must be Occam1DRegularization."
            )
        if self.regularization.n_layers != model.n_layers:
            raise ValueError(
                "regularization and model layer counts do not match."
            )
        self.multiplier_factors = self._factors(multiplier_factors)
        self.solver_policy = solver_policy or Occam1DSolverPolicy()
        if not isinstance(self.solver_policy, Occam1DSolverPolicy):
            raise TypeError(
                "solver_policy must be an Occam1DSolverPolicy or None."
            )
        self._lagrange_min, self._lagrange_max = self._multiplier_bounds(
            lagrange_bounds
        )
        if not (
            self._lagrange_min
            <= self.startup.lagrange_start
            <= self._lagrange_max
        ):
            raise ValueError(
                "startup.lagrange_start must lie inside lagrange_bounds."
            )
        self._lower_bound, self._upper_bound = (
            self.regularization._parameter_bounds(
                log_resistivity_bounds
            )
        )
        self._validate_initial_bounds(self.startup.parameters)
        self.rms_tolerance = self._nonnegative(
            "rms_tolerance",
            rms_tolerance,
        )
        self.model_tolerance = self._nonnegative(
            "model_tolerance",
            model_tolerance,
        )
        self.rms_change_tolerance = self._nonnegative(
            "rms_change_tolerance",
            rms_change_tolerance,
        )
        computational_model = model.with_resistivity(
            self.startup.physical_resistivity
        )
        self.forward = Occam1DForwardModel(
            computational_model,
            logger=self.logger,
        )
        self.jacobian = Occam1DJacobian(
            computational_model,
            forward=self.forward,
            logger=self.logger,
        )
        self._last_result = None

    @staticmethod
    def _nonnegative(name, value) -> float:
        """Return a finite non-negative scalar."""
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError(f"{name} must be a real number.")
        result = float(value)
        if not math.isfinite(result) or result < 0:
            raise ValueError(f"{name} must be finite and non-negative.")
        return result

    @staticmethod
    def _factors(values) -> np.ndarray:
        """Return a unique sorted positive multiplier-factor vector."""
        source = np.logspace(-3.0, 3.0, 13) if values is None else values
        try:
            result = np.asarray(source, dtype=float)
        except (TypeError, ValueError) as error:
            raise TypeError(
                "multiplier_factors must contain real numbers."
            ) from error
        if result.ndim != 1 or not result.size:
            raise ValueError("multiplier_factors must be a non-empty vector.")
        if np.any(~np.isfinite(result)) or np.any(result <= 0):
            raise ValueError(
                "multiplier_factors must be finite and strictly positive."
            )
        return np.unique(result)

    def _multiplier_bounds(self, values) -> tuple[float, float]:
        """Return the validated absolute multiplier interval."""
        source = (
            (self.config.lagrange_min, self.config.lagrange_max)
            if values is None
            else values
        )
        if not isinstance(source, (tuple, list)) or len(source) != 2:
            raise TypeError(
                "lagrange_bounds must be a (minimum, maximum) pair."
            )
        minimum = self._nonnegative("minimum Lagrange multiplier", source[0])
        maximum = self._nonnegative("maximum Lagrange multiplier", source[1])
        if maximum <= 0:
            raise ValueError(
                "maximum Lagrange multiplier must be strictly positive."
            )
        if minimum > maximum:
            raise ValueError(
                "minimum Lagrange multiplier cannot exceed maximum."
            )
        return minimum, maximum

    def _multipliers(self, center) -> np.ndarray:
        """Return bounded, sorted, and deduplicated trial multipliers."""
        with np.errstate(over="ignore", invalid="ignore"):
            products = center * self.multiplier_factors
        if np.any(~np.isfinite(products)):
            products = np.nan_to_num(
                products,
                nan=self._lagrange_max,
                posinf=self._lagrange_max,
                neginf=self._lagrange_min,
            )
        return np.unique(
            np.clip(
                products,
                self._lagrange_min,
                self._lagrange_max,
            )
        )

    def _validate_initial_bounds(self, parameters) -> None:
        """Require the supplied starting model to be feasible."""
        below = np.flatnonzero(parameters < self._lower_bound)
        above = np.flatnonzero(parameters > self._upper_bound)
        if below.size or above.size:
            pieces = []
            if below.size:
                pieces.append(f"below lower bound: {below.tolist()}")
            if above.size:
                pieces.append(f"above upper bound: {above.tolist()}")
            detail = "; ".join(pieces)
            raise ValueError(
                "Starting log10-resistivity parameters violate bounds "
                f"({detail}). Supply a feasible startup model."
            )

    @property
    def log_resistivity_bounds(self) -> tuple[np.ndarray, np.ndarray]:
        """Independent lower and upper log10-resistivity bounds."""
        return self._lower_bound.copy(), self._upper_bound.copy()

    @property
    def lagrange_bounds(self) -> tuple[float, float]:
        """Absolute inclusive multiplier search bounds."""
        return self._lagrange_min, self._lagrange_max

    def candidate_multipliers(self, center=None) -> np.ndarray:
        """Return trial multipliers after absolute bounding and deduplication.

        Parameters
        ----------
        center : float, optional
            Multiplier around which configured factors are applied. The
            startup multiplier is used when omitted.

        Returns
        -------
        ndarray
            Sorted unique multipliers inside :attr:`lagrange_bounds`.
        """
        value = (
            self.startup.lagrange_start
            if center is None
            else self._nonnegative("center multiplier", center)
        )
        return self._multipliers(value).copy()

    @property
    def last_result(self) -> Occam1DInversionResult | None:
        """Most recent completed result, if :meth:`run` has been called."""
        return self._last_result

    @property
    def controls(self):
        """Read-only effective native inversion controls."""
        return MappingProxyType(
            {
                "target_rms": self.startup.target_misfit,
                "max_iterations": self.startup.max_iterations,
                "stepsize_cut_count": self.startup.stepsize_cut_count,
                "roughness_type": self.startup.roughness_type,
                "lagrange_start": self.startup.lagrange_start,
                "lagrange_min": self._lagrange_min,
                "lagrange_max": self._lagrange_max,
            }
        )

    def _evaluate(
        self,
        parameters,
        number,
        multiplier,
        scale,
        rank,
        condition,
        active_lower=None,
        active_upper=None,
        solver="initial",
        solve_attempts=0,
        damping=0.0,
        initial_condition_number=None,
    ):
        """Evaluate one model with the exact nonlinear forward response."""
        prediction = self.forward.solver_vector(
            self.data,
            log10_resistivity=parameters,
        )
        observed = self.regularization.observed_solver_vector(self.data)
        errors = Occam1DJacobian.solver_errors(self.data)
        difference = prediction - observed
        phase_rows = []
        for index in range(self.data.n_frequencies):
            if self.data.resistivity_mask[index]:
                phase_rows.append(False)
            if self.data.phase_mask[index]:
                phase_rows.append(True)
        phase_rows = np.asarray(phase_rows, dtype=bool)
        difference[phase_rows] = (
            difference[phase_rows] + 90.0
        ) % 180.0 - 90.0
        residual = difference / errors
        rms = float(np.sqrt(np.mean(residual**2)))
        structure = self.regularization.evaluate(parameters)
        return Occam1DIteration(
            number=number,
            parameters=parameters,
            prediction=prediction,
            residual=residual,
            rms=rms,
            roughness=structure.objective,
            multiplier=multiplier,
            step_scale=scale,
            rank=rank,
            condition_number=condition,
            active_lower=active_lower,
            active_upper=active_upper,
            solver=solver,
            solve_attempts=solve_attempts,
            damping=damping,
            initial_condition_number=initial_condition_number,
        )

    def _trial(self, current, candidate, number, rejected):
        """Apply deterministic step cuts and return the best true trial."""
        direction = candidate.parameters - current.parameters
        trials = []
        for cut in range(self.startup.stepsize_cut_count + 1):
            scale = 0.5**cut
            parameters = current.parameters + scale * direction
            try:
                trial = self._evaluate(
                    parameters,
                    number,
                    candidate.multiplier,
                    scale,
                    candidate.rank,
                    candidate.condition_number,
                    candidate.active_lower,
                    candidate.active_upper,
                    candidate.solver,
                    candidate.solve_attempts,
                    candidate.damping,
                    candidate.initial_condition_number,
                )
            except (ArithmeticError, RuntimeError, ValueError) as error:
                rejected.append(
                    Occam1DRejectedCandidate(
                        iteration=number,
                        multiplier=candidate.multiplier,
                        step_scale=scale,
                        reason=(
                            Occam1DRejectionReason.FORWARD_EVALUATION_FAILED
                        ),
                        message=str(error) or "Forward evaluation failed.",
                        exception_type=type(error).__name__,
                    )
                )
                continue
            trials.append(trial)
            if trial.rms <= current.rms:
                break
        if not trials:
            return None
        selected = min(trials, key=lambda item: item.rms)
        for trial in trials:
            if trial is selected:
                continue
            rejected.append(
                Occam1DRejectedCandidate(
                    iteration=number,
                    multiplier=trial.multiplier,
                    step_scale=trial.step_scale,
                    rms=trial.rms,
                    roughness=trial.roughness,
                    reason=Occam1DRejectionReason.STEP_CUT_DISCARDED,
                    message=(
                        "Another step cut from the same linear candidate "
                        "had a lower nonlinear RMS."
                    ),
                )
            )
        return selected

    def _select(self, trials):
        """Select the smoothest target trial or best-fitting trial."""
        limit = self.startup.target_misfit * (1.0 + self.rms_tolerance)
        acceptable = [item for item in trials if item.rms <= limit]
        if acceptable:
            return min(
                acceptable,
                key=lambda item: (item.roughness, item.rms),
            )
        return min(trials, key=lambda item: (item.rms, item.roughness))

    def run(
        self,
        *,
        restart: Occam1DRestart | None = None,
        callback: Callable[[Occam1DIteration], Any] | None = None,
        cancel: Callable[[], bool] | None = None,
    ) -> Occam1DInversionResult:
        """Run the nonlinear inversion and return complete model history.

        Parameters
        ----------
        restart : Occam1DRestart, optional
            Previously accepted history to continue. Iteration numbering and
            the selected multiplier resume from its final record. The
            configured maximum iteration count remains an absolute total.
        callback : callable, optional
            Called after iteration zero and every accepted update. Returning
            values are ignored; callback exceptions propagate.
        cancel : callable returning bool, optional
            Cooperative cancellation check performed before every update.

        Returns
        -------
        Occam1DInversionResult
            Immutable history, terminal reason, and final model.

        Raises
        ------
        TypeError
            If callbacks are not callable.

        Notes
        -----
        Expected numerical candidate failures are returned through
        ``rejected_candidates`` and ``failed_iterations``. Programming errors
        and callback exceptions are not swallowed.
        """
        if callback is not None and not callable(callback):
            raise TypeError("callback must be callable or None.")
        if cancel is not None and not callable(cancel):
            raise TypeError("cancel must be callable or None.")
        if restart is None:
            initial = self._evaluate(
                self.startup.parameters,
                0,
                self.startup.lagrange_start,
                0.0,
                0,
                float("inf"),
            )
            history = [initial]
            multiplier = self.startup.lagrange_start
            start_number = 1
            if callback is not None:
                callback(initial)
        else:
            self._validate_restart(restart)
            history = list(restart.iterations)
            initial = history[0]
            multiplier = restart.multiplier
            start_number = restart.iteration_number + 1
        rejected = (
            [] if restart is None else list(restart.rejected_candidates)
        )
        failures = (
            [] if restart is None else list(restart.failed_iterations)
        )
        target = self.startup.target_misfit
        if history[-1].rms <= target * (1.0 + self.rms_tolerance):
            return self._finish(
                history,
                Occam1DConvergence.TARGET,
                "The active model satisfies the target RMS.",
                rejected,
                failures,
            )

        status = Occam1DConvergence.MAX_ITERATIONS
        message = "The maximum iteration count was reached."
        remaining = max(0, self.startup.max_iterations - start_number + 1)
        with get_progress_bar(
            total=remaining,
            desc=f"Occam1D {self.data.station}",
            unit="iteration",
            verbose=self.verbose,
        ) as progress:
            for number in range(
                start_number,
                self.startup.max_iterations + 1,
            ):
                if cancel is not None and cancel():
                    status = Occam1DConvergence.CANCELLED
                    message = "The inversion was cancelled by the caller."
                    break
                current = history[-1]
                rejected_start = len(rejected)
                multipliers = self._multipliers(multiplier)
                try:
                    sensitivity = self.jacobian.compute(
                        self.data,
                        parameters=current.parameters,
                        weighted=True,
                    )
                except (ArithmeticError, RuntimeError, ValueError) as error:
                    status = Occam1DConvergence.FAILED
                    message = "Jacobian evaluation failed."
                    failures.append(
                        Occam1DFailedIteration(
                            iteration=number,
                            current_rms=current.rms,
                            attempted_candidates=0,
                            rejected_candidates=0,
                            reason="jacobian_failed",
                            message=str(error) or message,
                            recoverable=False,
                        )
                    )
                    break
                trials = []
                for trial_multiplier in multipliers:
                    try:
                        candidate = self.regularization.solve_linearized(
                            self.data,
                            sensitivity,
                            current.parameters,
                            trial_multiplier,
                            bounds=(self._lower_bound, self._upper_bound),
                            policy=self.solver_policy,
                        )
                    except (
                        ArithmeticError,
                        RuntimeError,
                        ValueError,
                    ) as error:
                        rejected.append(
                            Occam1DRejectedCandidate(
                                iteration=number,
                                multiplier=float(trial_multiplier),
                                reason=(
                                    Occam1DRejectionReason.LINEAR_SOLVE_FAILED
                                ),
                                message=(
                                    str(error) or "Linear solve failed."
                                ),
                                exception_type=type(error).__name__,
                            )
                        )
                        continue
                    trial = self._trial(
                        current,
                        candidate,
                        number,
                        rejected,
                    )
                    if trial is not None:
                        trials.append(trial)
                if not trials:
                    status = Occam1DConvergence.FAILED
                    message = (
                        "No finite nonlinear model candidate could be "
                        "evaluated."
                    )
                    failures.append(
                        Occam1DFailedIteration(
                            iteration=number,
                            current_rms=current.rms,
                            attempted_candidates=int(multipliers.size),
                            rejected_candidates=(
                                len(rejected) - rejected_start
                            ),
                            reason="no_valid_candidate",
                            message=message,
                            recoverable=True,
                        )
                    )
                    break
                selected = self._select(trials)
                for trial in trials:
                    if trial is selected:
                        continue
                    rejected.append(
                        Occam1DRejectedCandidate(
                            iteration=number,
                            multiplier=trial.multiplier,
                            step_scale=trial.step_scale,
                            rms=trial.rms,
                            roughness=trial.roughness,
                            reason=Occam1DRejectionReason.NOT_SELECTED,
                            message=(
                                "A different valid candidate better "
                                "satisfied the Occam selection policy."
                            ),
                        )
                    )
                if selected.rms > current.rms + (
                    np.finfo(float).eps * max(1.0, current.rms)
                ):
                    status = Occam1DConvergence.STAGNATED
                    message = (
                        "No step-cut candidate improved the nonlinear "
                        "normalized RMS."
                    )
                    rejected.append(
                        Occam1DRejectedCandidate(
                            iteration=number,
                            multiplier=selected.multiplier,
                            step_scale=selected.step_scale,
                            rms=selected.rms,
                            roughness=selected.roughness,
                            reason=Occam1DRejectionReason.RMS_REGRESSION,
                            message=message,
                        )
                    )
                    failures.append(
                        Occam1DFailedIteration(
                            iteration=number,
                            current_rms=current.rms,
                            attempted_candidates=int(multipliers.size),
                            rejected_candidates=(
                                len(rejected) - rejected_start
                            ),
                            reason="no_improving_candidate",
                            message=message,
                            recoverable=True,
                        )
                    )
                    break
                history.append(selected)
                multiplier = selected.multiplier
                progress.update(
                    1,
                    metrics={
                        "rms": selected.rms,
                        "roughness": selected.roughness,
                        "mu": selected.multiplier,
                    },
                )
                self.logger.debug(
                    "Occam1D iteration %d: rms=%g roughness=%g mu=%g",
                    number,
                    selected.rms,
                    selected.roughness,
                    selected.multiplier,
                )
                if callback is not None:
                    callback(selected)
                if selected.rms <= target * (1.0 + self.rms_tolerance):
                    status = Occam1DConvergence.TARGET
                    message = "The requested normalized RMS was reached."
                    break
                model_change = float(
                    np.max(np.abs(selected.parameters - current.parameters))
                )
                rms_change = max(0.0, current.rms - selected.rms) / max(
                    current.rms,
                    np.finfo(float).tiny,
                )
                if (
                    model_change <= self.model_tolerance
                    and rms_change <= self.rms_change_tolerance
                ):
                    status = Occam1DConvergence.STAGNATED
                    message = (
                        "Model and normalized RMS changes fell below the "
                        "configured stagnation tolerances."
                    )
                    break
        return self._finish(
            history,
            status,
            message,
            rejected,
            failures,
        )

    def _validate_restart(self, restart) -> None:
        """Validate checkpoint compatibility with this inversion."""
        if not isinstance(restart, Occam1DRestart):
            raise TypeError("restart must be an Occam1DRestart or None.")
        if restart.depth.shape != self.model.depth.shape or not np.allclose(
            restart.depth,
            self.model.depth,
            rtol=0.0,
            atol=0.0,
        ):
            raise ValueError(
                "Restart layer geometry does not match the inversion model."
            )
        if restart.iterations[-1].prediction.size != self.data.n_data:
            raise ValueError(
                "Restart data count does not match the inversion sounding."
            )
        expected = _data_fingerprint(self.data)
        if (
            restart.data_fingerprint is not None
            and restart.data_fingerprint != expected
        ):
            raise ValueError(
                "Restart sounding fingerprint does not match inversion data."
            )
        tolerance = np.finfo(float).eps * max(1.0, restart.target_rms)
        if abs(restart.target_rms - self.startup.target_misfit) > tolerance:
            raise ValueError(
                "Restart target RMS does not match inversion controls."
            )
        self._validate_initial_bounds(restart.parameters)
        if not self._lagrange_min <= restart.multiplier <= self._lagrange_max:
            raise ValueError(
                "Restart multiplier lies outside lagrange_bounds."
            )

    def restart(self, result=None) -> Occam1DRestart:
        """Create a continuation checkpoint from a completed result.

        Parameters
        ----------
        result : Occam1DInversionResult, optional
            Result to checkpoint. The most recent result is used when
            omitted.
        """
        selected = self.last_result if result is None else result
        if selected is None:
            raise RuntimeError("No inversion result is available to restart.")
        if not isinstance(selected, Occam1DInversionResult):
            raise TypeError("result must be an Occam1DInversionResult.")
        return Occam1DRestart(
            iterations=selected.iterations,
            depth=self.model.depth,
            target_rms=selected.target_rms,
            data_fingerprint=_data_fingerprint(self.data),
            previous_convergence=selected.convergence,
            previous_message=selected.message,
            rejected_candidates=selected.rejected_candidates,
            failed_iterations=selected.failed_iterations,
        )

    def export_result(self, directory, result=None) -> dict[str, Path]:
        """Export native inversion text products and diagnostic metadata.

        Parameters
        ----------
        directory : path-like
            Destination ``model-text`` directory.
        result : Occam1DInversionResult, optional
            Result to export. The latest completed result is used by default.

        Returns
        -------
        dict
            Product names mapped to absolute output paths.
        """
        selected = self._selected_result(result)
        output = Path(directory).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=True)
        prefix = self._safe_name(self.data.station)
        paths = {
            "model": output / f"{prefix}_model.csv",
            "response": output / f"{prefix}_response.csv",
            "iterations": output / f"{prefix}_iterations.csv",
            "rejected": output / f"{prefix}_rejected_candidates.csv",
            "failures": output / f"{prefix}_failed_iterations.csv",
            "summary": output / f"{prefix}_summary.txt",
            "metadata": output / f"{prefix}_metadata.json",
        }
        self._write_model_csv(paths["model"], selected)
        self._write_response_csv(paths["response"], selected)
        self._write_iteration_csv(paths["iterations"], selected)
        self._write_records_csv(
            paths["rejected"],
            [item.to_dict() for item in selected.rejected_candidates],
        )
        self._write_records_csv(
            paths["failures"],
            [item.to_dict() for item in selected.failed_iterations],
        )
        paths["summary"].write_text(
            self.result_summary(selected),
            encoding="utf8",
        )
        metadata = {
            "schema": "pycsamt.occam1d.native-result/v1",
            "station": self.data.station,
            "mode": self.data.mode,
            "config": self.config.to_dict(),
            "inversion": selected.diagnostics(),
            "solver": self.diagnostics(),
        }
        paths["metadata"].write_text(
            json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False)
            + "\n",
            encoding="utf8",
        )
        return paths

    def save_main_images(
        self,
        directory,
        result=None,
        *,
        dpi=180,
        fmt="png",
        style=None,
        style_overrides=None,
    ) -> dict[str, Path]:
        """Save model, response, convergence, and summary figures.

        ``style`` accepts a named Occam1D preset or an explicit style object.
        ``style_overrides`` applies double-underscore paths such as
        ``observed__marker`` or ``target__visible`` to this export only.
        """
        import matplotlib.pyplot as plt

        from ...api.occam1d import resolve_occam1d_style

        selected = self._selected_result(result)
        if isinstance(dpi, bool) or not isinstance(dpi, int) or dpi < 1:
            raise ValueError("dpi must be a positive integer.")
        if not isinstance(fmt, str) or not fmt.strip(". "):
            raise ValueError("fmt must be a non-empty format name.")
        if style_overrides is None:
            style_overrides = {}
        if not isinstance(style_overrides, dict):
            raise TypeError("style_overrides must be a mapping or None.")
        plot_style = resolve_occam1d_style(style, **style_overrides)
        output = Path(directory).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=True)
        prefix = self._safe_name(self.data.station)
        extension = fmt.lower().strip(". ")
        factories = {
            "model": lambda: self._plot_native_model(selected, plot_style),
            "response": lambda: self._plot_native_response(
                selected, plot_style
            ),
            "convergence": lambda: self._plot_native_convergence(
                selected, plot_style
            ),
            "summary": lambda: self._plot_native_summary(
                selected, plot_style
            ),
        }
        paths = {}
        for kind, factory in factories.items():
            figure = factory()
            path = output / f"{prefix}_occam1d_{kind}.{extension}"
            try:
                figure.savefig(
                    path,
                    dpi=dpi,
                    bbox_inches="tight",
                    format=extension,
                )
            finally:
                plt.close(figure)
            paths[kind] = path
        return paths

    def _selected_result(self, result):
        """Return an explicit or latest native result."""
        selected = self.last_result if result is None else result
        if selected is None:
            raise RuntimeError("No native inversion result is available.")
        if not isinstance(selected, Occam1DInversionResult):
            raise TypeError("result must be an Occam1DInversionResult.")
        return selected

    @staticmethod
    def _safe_name(value):
        """Return a conservative product filename prefix."""
        text = "".join(
            character if character.isalnum() or character in "-_" else "_"
            for character in str(value).strip()
        )
        return text or "station"

    @staticmethod
    def _write_records_csv(path, records):
        """Write homogeneous dictionaries, including an empty header."""
        fields = sorted({key for record in records for key in record})
        with path.open("w", newline="", encoding="utf8") as stream:
            if not fields:
                stream.write("record\n")
                return
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(records)

    def _write_model_csv(self, path, result):
        """Write final physical model by layer."""
        with path.open("w", newline="", encoding="utf8") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                ["layer", "depth_top_m", "log10_resistivity", "resistivity"]
            )
            for index, (depth, parameter, resistivity) in enumerate(
                zip(
                    self.model.depth,
                    result.final.parameters,
                    result.final.resistivity,
                ),
                start=1,
            ):
                writer.writerow([index, depth, parameter, resistivity])

    def _write_response_csv(self, path, result):
        """Write observations, final predictions, errors, and residuals."""
        observed = self.regularization.observed_solver_vector(self.data)
        errors = Occam1DJacobian.solver_errors(self.data)
        rows = []
        cursor = 0
        for index, frequency in enumerate(self.data.frequency):
            for kind, present in (
                ("log10_resistivity", self.data.resistivity_mask[index]),
                ("phase_deg", self.data.phase_mask[index]),
            ):
                if present:
                    rows.append(
                        (
                            frequency,
                            kind,
                            observed[cursor],
                            result.final.prediction[cursor],
                            errors[cursor],
                            result.final.residual[cursor],
                        )
                    )
                    cursor += 1
        with path.open("w", newline="", encoding="utf8") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                [
                    "frequency_hz",
                    "quantity",
                    "observed",
                    "predicted",
                    "standard_error",
                    "weighted_residual",
                ]
            )
            writer.writerows(rows)

    @staticmethod
    def _write_iteration_csv(path, result):
        """Write accepted convergence and solver diagnostics."""
        fields = list(result.iterations[0].diagnostics())
        with path.open("w", newline="", encoding="utf8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(item.diagnostics() for item in result.iterations)

    def result_summary(self, result=None) -> str:
        """Return a concise report for a native inversion result."""
        selected = self._selected_result(result)
        return (
            "pyCSAMT native Occam1D inversion\n"
            f"  station       : {self.data.station}\n"
            f"  mode          : {self.data.mode}\n"
            f"  frequencies   : {self.data.n_frequencies}\n"
            f"  observations  : {self.data.n_data}\n"
            f"  layers        : {self.model.n_layers}\n"
            f"  status        : {selected.convergence.value}\n"
            f"  iterations    : {selected.n_iterations}\n"
            f"  initial RMS   : {selected.initial.rms:.8g}\n"
            f"  final RMS     : {selected.final.rms:.8g}\n"
            f"  target RMS    : {selected.target_rms:.8g}\n"
            f"  roughness     : {selected.final.roughness:.8g}\n"
            f"  multiplier    : {selected.final.multiplier:.8g}\n"
            f"  rejected      : {len(selected.rejected_candidates)}\n"
            f"  failed steps  : {len(selected.failed_iterations)}\n"
            f"  message       : {selected.message}\n"
        )

    def _response_arrays(self, result):
        """Return final predictions on the original frequency grid."""
        rho = np.full(self.data.n_frequencies, np.nan)
        phase = np.full(self.data.n_frequencies, np.nan)
        cursor = 0
        for index in range(self.data.n_frequencies):
            if self.data.resistivity_mask[index]:
                rho[index] = 10.0 ** result.final.prediction[cursor]
                cursor += 1
            if self.data.phase_mask[index]:
                phase[index] = result.final.prediction[cursor]
                cursor += 1
        return rho, phase

    def _plot_native_model(self, result, style):
        """Return the native final-model figure."""
        import matplotlib.pyplot as plt

        figure, axis = plt.subplots(figsize=(5.2, 6.4))
        if style.model.visible:
            axis.step(
                result.final.resistivity,
                self.model.depth,
                where="post",
                **style.model.kwargs(),
            )
        axis.set_xscale("log")
        axis.invert_yaxis()
        axis.grid(
            True,
            which="both",
            alpha=style.grid_alpha,
            linewidth=style.grid_linewidth,
        )
        axis.set_xlabel("Resistivity (ohm m)")
        axis.set_ylabel("Depth (m)")
        self._native_legend(axis, style, style.model_legend)
        axis.set_title(f"{self.data.station} — Occam1D model")
        figure.tight_layout()
        return figure

    def _plot_native_response(self, result, style):
        """Return observed and predicted sounding panels."""
        import matplotlib.pyplot as plt

        predicted_rho, predicted_phase = self._response_arrays(result)
        period = self.data.period
        figure, axes = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)
        if style.observed.visible:
            axes[0].loglog(
                period, self.data.resistivity, **style.observed.kwargs()
            )
        if style.predicted.visible:
            axes[0].loglog(
                period, predicted_rho, **style.predicted.kwargs()
            )
        axes[0].set_ylabel("Apparent resistivity (ohm m)")
        if style.observed.visible:
            axes[1].semilogx(
                period, self.data.phase, **style.observed.kwargs()
            )
        if style.phase_predicted.visible:
            axes[1].semilogx(
                period, predicted_phase, **style.phase_predicted.kwargs()
            )
        axes[1].set_ylabel("Phase (degrees)")
        axes[1].set_xlabel("Period (s)")
        for axis in axes:
            axis.grid(
                True,
                which="both",
                alpha=style.grid_alpha,
                linewidth=style.grid_linewidth,
            )
            self._native_legend(axis, style, style.response_legend)
        figure.suptitle(f"{self.data.station} — Occam1D response")
        figure.tight_layout()
        return figure

    def _plot_native_convergence(self, result, style):
        """Return normalized RMS and roughness history."""
        import matplotlib.pyplot as plt

        numbers = [item.number for item in result.iterations]
        rms = [item.rms for item in result.iterations]
        roughness = [item.roughness for item in result.iterations]
        figure, left = plt.subplots(figsize=(7.2, 4.8))
        right = left.twinx()
        if style.iteration.visible:
            left.plot(numbers, rms, **style.iteration.kwargs())
        if style.target.visible:
            left.axhline(
                result.target_rms,
                **style.target.kwargs(include_marker=False),
            )
        if style.roughness.visible:
            right.plot(numbers, roughness, **style.roughness.kwargs())
        left.set_xlabel("Iteration")
        left.set_ylabel("Normalized RMS", color=style.iteration.color)
        right.set_ylabel("Roughness", color=style.roughness.color)
        left.grid(
            True,
            alpha=style.grid_alpha,
            linewidth=style.grid_linewidth,
        )
        self._native_legend(left, style, style.convergence_legend)
        self._native_legend(right, style, style.roughness_legend)
        left.set_title(f"{self.data.station} — Occam1D convergence")
        figure.tight_layout()
        return figure

    def _plot_native_summary(self, result, style):
        """Return a compact model, response, and convergence dashboard."""
        import matplotlib.pyplot as plt

        predicted_rho, predicted_phase = self._response_arrays(result)
        period = self.data.period
        figure = plt.figure(figsize=(12.0, 7.2), constrained_layout=True)
        grid = figure.add_gridspec(2, 3, width_ratios=(0.85, 1.2, 1.1))
        model_axis = figure.add_subplot(grid[:, 0])
        rho_axis = figure.add_subplot(grid[0, 1])
        phase_axis = figure.add_subplot(grid[1, 1], sharex=rho_axis)
        convergence_axis = figure.add_subplot(grid[:, 2])
        if style.model.visible:
            model_axis.step(
                result.final.resistivity,
                self.model.depth,
                where="post",
                **style.model.kwargs(),
            )
        model_axis.set_xscale("log")
        model_axis.invert_yaxis()
        model_axis.set_xlabel("Resistivity (ohm m)")
        model_axis.set_ylabel("Depth (m)")
        if style.observed.visible:
            rho_axis.loglog(
                period, self.data.resistivity, **style.observed.kwargs()
            )
        if style.predicted.visible:
            rho_axis.loglog(
                period, predicted_rho, **style.predicted.kwargs()
            )
        rho_axis.set_ylabel("App. resistivity")
        if style.observed.visible:
            phase_axis.semilogx(
                period, self.data.phase, **style.observed.kwargs()
            )
        if style.phase_predicted.visible:
            phase_axis.semilogx(
                period, predicted_phase, **style.phase_predicted.kwargs()
            )
        phase_axis.set_ylabel("Phase (degrees)")
        phase_axis.set_xlabel("Period (s)")
        numbers = [item.number for item in result.iterations]
        if style.iteration.visible:
            convergence_axis.plot(
                numbers, result.rms_history, **style.iteration.kwargs()
            )
        if style.target.visible:
            convergence_axis.axhline(
                result.target_rms,
                **style.target.kwargs(include_marker=False),
            )
        convergence_axis.set_xlabel("Iteration")
        convergence_axis.set_ylabel("Normalized RMS")
        for axis in (model_axis, rho_axis, phase_axis, convergence_axis):
            axis.grid(
                True,
                which="both",
                alpha=style.grid_alpha,
                linewidth=style.grid_linewidth,
            )
        self._native_legend(model_axis, style, style.model_legend)
        self._native_legend(rho_axis, style, style.response_legend)
        self._native_legend(
            convergence_axis, style, style.convergence_legend
        )
        figure.suptitle(
            f"{self.data.station} - native Occam1D "
            f"({result.convergence.value}, RMS={result.final.rms:.3g})"
        )
        return figure

    @staticmethod
    def _native_legend(axis, style, visible):
        """Draw the configured legend only when labeled artists exist."""
        handles, _ = axis.get_legend_handles_labels()
        if visible and handles:
            axis.legend(
                frameon=style.legend_frame,
                fontsize=style.legend_fontsize,
                loc=style.legend_location,
            )

    def _finish(
        self,
        history,
        status,
        message,
        rejected=(),
        failures=(),
    ):
        """Create, retain, and log one terminal result."""
        result = Occam1DInversionResult(
            iterations=tuple(history),
            convergence=status,
            target_rms=self.startup.target_misfit,
            message=message,
            rejected_candidates=tuple(rejected),
            failed_iterations=tuple(failures),
        )
        self._last_result = result
        self.logger.info(
            "Occam1D inversion ended as %s after %d updates (RMS=%g).",
            status.value,
            result.n_iterations,
            result.final.rms,
        )
        return result

    def diagnostics(self) -> dict[str, Any]:
        """Extend lifecycle diagnostics with solver configuration."""
        values = super().diagnostics()
        values.update(
            {
                "station": self.data.station,
                "n_data": self.data.n_data,
                "n_layers": self.model.n_layers,
                "controls": dict(self.controls),
                "multiplier_factors": self.multiplier_factors.tolist(),
                "lagrange_min": self._lagrange_min,
                "lagrange_max": self._lagrange_max,
                "log_resistivity_lower": [
                    None if np.isneginf(value) else float(value)
                    for value in self._lower_bound
                ],
                "log_resistivity_upper": [
                    None if np.isposinf(value) else float(value)
                    for value in self._upper_bound
                ],
                "rms_tolerance": self.rms_tolerance,
                "model_tolerance": self.model_tolerance,
                "rms_change_tolerance": self.rms_change_tolerance,
                "solver_policy": {
                    "condition_limit": self.solver_policy.condition_limit,
                    "ill_condition_action": (
                        self.solver_policy.ill_condition_action
                    ),
                    "svd_failure_action": (
                        self.solver_policy.svd_failure_action
                    ),
                    "max_retries": self.solver_policy.max_retries,
                    "damping_start": self.solver_policy.damping_start,
                    "damping_growth": self.solver_policy.damping_growth,
                },
                "has_result": self.last_result is not None,
            }
        )
        return values
