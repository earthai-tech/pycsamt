# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Occam-specific smoothness operators and stable candidate solves.

The paper-compatible default penalizes unscaled first or second differences
of layer log10 resistivity. Optional depth scaling and reference penalties
are explicit extensions; they never alter the default Occam objective.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from .base import Occam1DBase
from .data import Occam1DData
from .jacobian import Occam1DJacobian, Occam1DJacobianResult
from .model import Occam1DModel

__all__ = [
    "Occam1DCandidate",
    "Occam1DLinearSolveError",
    "Occam1DRegularization",
    "Occam1DRegularizationResult",
    "Occam1DSolverPolicy",
]

_LN10 = math.log(10.0)


class Occam1DLinearSolveError(RuntimeError):
    """Raised after an explicit linear-solver policy is exhausted."""


@dataclass(frozen=True)
class Occam1DSolverPolicy:
    r"""Control linear fallback and ill-conditioning stabilization.

    Parameters
    ----------
    condition_limit : float, default=1e12
        Maximum accepted augmented-system condition number.
    ill_condition_action : {"accept", "reject", "damp"}, default="damp"
        Accept the candidate with diagnostics, reject it, or retry after
        adding :math:`\sqrt{\lambda}I` rows to the augmented system.
    svd_failure_action : {"raise", "fallback"}, default="fallback"
        Raise immediately or retry using SciPy's independent LAPACK path.
    max_retries : int, default=3
        Maximum diagonal-damping retries after the initial solve.
    damping_start : float, default=1e-10
        First positive diagonal stabilization objective weight.
    damping_growth : float, default=100
        Multiplicative damping increase between retries.
    """

    condition_limit: float = 1e12
    ill_condition_action: str = "damp"
    svd_failure_action: str = "fallback"
    max_retries: int = 3
    damping_start: float = 1e-10
    damping_growth: float = 100.0

    def __post_init__(self) -> None:
        """Validate all policy domains."""
        for name in ("condition_limit", "damping_start", "damping_growth"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError(f"{name} must be a real number.")
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive.")
        if self.damping_growth <= 1:
            raise ValueError("damping_growth must exceed one.")
        if (
            isinstance(self.max_retries, bool)
            or not isinstance(self.max_retries, int)
            or self.max_retries < 0
        ):
            raise ValueError("max_retries must be a non-negative integer.")
        if self.ill_condition_action not in {"accept", "reject", "damp"}:
            raise ValueError(
                "ill_condition_action must be 'accept', 'reject', or 'damp'."
            )
        if self.svd_failure_action not in {"raise", "fallback"}:
            raise ValueError(
                "svd_failure_action must be 'raise' or 'fallback'."
            )


def _readonly(values, dtype=float) -> np.ndarray:
    """Return a detached immutable NumPy array."""
    result = np.array(values, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class Occam1DRegularizationResult:
    r"""Immutable decomposition of one model-structure objective.

    ``smoothness`` is :math:`D m`; ``reference`` is
    :math:`P(m-m_{ref})`. The reported objective is their combined squared
    Euclidean norm.
    """

    parameters: np.ndarray
    smoothness: np.ndarray
    reference: np.ndarray
    roughness: float
    reference_penalty: float
    objective: float
    order: int
    scaling: str

    def __post_init__(self) -> None:
        """Validate immutable arrays and scalar decomposition."""
        parameters = _readonly(self.parameters)
        smoothness = _readonly(self.smoothness)
        reference = _readonly(self.reference)
        for name, value in (
            ("roughness", self.roughness),
            ("reference_penalty", self.reference_penalty),
            ("objective", self.objective),
        ):
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and non-negative.")
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "smoothness", smoothness)
        object.__setattr__(self, "reference", reference)

    @property
    def norm(self) -> float:
        """Square root of the total structure objective."""
        return math.sqrt(self.objective)

    def diagnostics(self) -> dict[str, Any]:
        """Return JSON-safe structure metrics."""
        return {
            "order": self.order,
            "scaling": self.scaling,
            "roughness": self.roughness,
            "reference_penalty": self.reference_penalty,
            "objective": self.objective,
            "norm": self.norm,
        }


@dataclass(frozen=True)
class Occam1DCandidate:
    r"""One regularized solution of a linearized Occam subproblem.

    Parameters are log10 layer resistivities. ``data_objective`` and
    ``structure_objective`` are squared residual norms, so

    .. math::

       \Phi = \Phi_d + \mu\Phi_m.
    """

    parameters: np.ndarray
    multiplier: float
    data_residual: np.ndarray
    structure_residual: np.ndarray
    data_objective: float
    structure_objective: float
    objective: float
    rank: int
    singular_values: np.ndarray
    condition_number: float
    solver: str = "svd_lstsq"
    active_lower: np.ndarray | None = None
    active_upper: np.ndarray | None = None
    solve_attempts: int = 1
    damping: float = 0.0
    initial_condition_number: float | None = None

    def __post_init__(self) -> None:
        """Detach arrays and validate candidate diagnostics."""
        parameters = _readonly(self.parameters)
        data_residual = _readonly(self.data_residual)
        structure_residual = _readonly(self.structure_residual)
        singular_values = _readonly(self.singular_values)
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
        if parameters.ndim != 1 or parameters.size < 2:
            raise ValueError("candidate parameters must be one-dimensional.")
        arrays = (parameters, data_residual, structure_residual)
        if any(np.any(~np.isfinite(value)) for value in arrays):
            raise ValueError("candidate arrays must contain finite values.")
        if active_lower.shape != parameters.shape:
            raise ValueError("active_lower must align with parameters.")
        if active_upper.shape != parameters.shape:
            raise ValueError("active_upper must align with parameters.")
        if not math.isfinite(self.multiplier) or self.multiplier < 0:
            raise ValueError("multiplier must be finite and non-negative.")
        for name, value in (
            ("data_objective", self.data_objective),
            ("structure_objective", self.structure_objective),
            ("objective", self.objective),
        ):
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and non-negative.")
        if self.solve_attempts < 1:
            raise ValueError("solve_attempts must be at least one.")
        if not math.isfinite(self.damping) or self.damping < 0:
            raise ValueError("damping must be finite and non-negative.")
        initial_condition = self.initial_condition_number
        if initial_condition is not None and (
            math.isnan(initial_condition) or initial_condition < 0
        ):
            raise ValueError(
                "initial_condition_number must be non-negative or None."
            )
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "data_residual", data_residual)
        object.__setattr__(self, "structure_residual", structure_residual)
        object.__setattr__(self, "singular_values", singular_values)
        object.__setattr__(self, "active_lower", active_lower)
        object.__setattr__(self, "active_upper", active_upper)

    @property
    def parameter_norm(self) -> float:
        """Euclidean norm of the absolute parameter solution."""
        return float(np.linalg.norm(self.parameters))

    @property
    def n_active_bounds(self) -> int:
        """Number of parameters lying at either configured bound."""
        return int(np.count_nonzero(self.active_lower | self.active_upper))

    def diagnostics(self) -> dict[str, Any]:
        """Return JSON-safe linear-solve diagnostics."""
        return {
            "multiplier": self.multiplier,
            "data_objective": self.data_objective,
            "structure_objective": self.structure_objective,
            "objective": self.objective,
            "rank": self.rank,
            "condition_number": (
                self.condition_number
                if math.isfinite(self.condition_number)
                else None
            ),
            "solver": self.solver,
            "solve_attempts": self.solve_attempts,
            "damping": self.damping,
            "initial_condition_number": (
                self.initial_condition_number
                if self.initial_condition_number is not None
                and math.isfinite(self.initial_condition_number)
                else None
            ),
            "n_active_bounds": self.n_active_bounds,
            "active_lower": self.active_lower.tolist(),
            "active_upper": self.active_upper.tolist(),
            "singular_values": self.singular_values.tolist(),
        }


class Occam1DRegularization(Occam1DBase):
    r"""Construct and apply an Occam1D quadratic structure penalty.

    Parameters
    ----------
    n_layers : int
        Number of log10-resistivity parameters, including the half-space.
    order : {1, 2}, default=1
        First- or second-difference roughness.
    depth : array-like of shape (n_layers,), optional
        Strictly increasing layer-top depths in metres. Required only for
        ``scaling="depth"``.
    scaling : {"index", "depth"}, default="index"
        ``index`` reproduces equations (2) of Constable et al. (1987).
        ``depth`` approximates integrated squared physical derivatives using
        irregular layer-top spacing.
    row_weights : array-like, optional
        Non-negative amplitude weights applied to roughness rows. Length is
        ``n_layers - order``. Values multiply residuals, and are therefore
        squared in the objective.
    reference : array-like of shape (n_layers,), optional
        Preferred log10-resistivity model.
    reference_weights : array-like of shape (n_layers,), optional
        Non-negative amplitude weights for ``parameters - reference``.
        Positive values require ``reference``. Zero-weight rows are omitted.

    Notes
    -----
    The default objective is exactly ``||D m||**2``. A large Lagrange
    multiplier in :meth:`solve` therefore emphasizes smoothness, matching the
    convention adopted for the native pyCSAMT inversion engine.
    """

    def __init__(
        self,
        n_layers,
        *,
        order=1,
        depth=None,
        scaling="index",
        row_weights=None,
        reference=None,
        reference_weights=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.n_layers = self._integer("n_layers", n_layers, minimum=2)
        self.order = self._integer("order", order, minimum=1)
        if self.order not in {1, 2}:
            raise ValueError("order must be 1 or 2.")
        if not isinstance(scaling, str):
            raise TypeError("scaling must be a string.")
        self.scaling = scaling.strip().lower()
        if self.scaling not in {"index", "depth"}:
            raise ValueError("scaling must be 'index' or 'depth'.")
        self._depth = self._depth_vector(depth)
        self._operator = self._build_operator()
        self._row_weights = self._weights(
            row_weights,
            self._operator.shape[0],
            "row_weights",
            default=1.0,
        )
        self._operator *= self._row_weights[:, None]
        self._reference = self._optional_vector(reference, "reference")
        self._reference_weights = self._weights(
            reference_weights,
            self.n_layers,
            "reference_weights",
            default=0.0,
        )
        if np.any(self._reference_weights > 0) and self._reference is None:
            raise ValueError(
                "reference is required when reference_weights are positive."
            )

    @classmethod
    def from_model(
        cls,
        model,
        *,
        order=1,
        scaling="index",
        use_native_preference=False,
        **kwargs,
    ):
        """Construct from model geometry and optional native preferences.

        ``use_native_preference=True`` interprets ``model.preference`` as
        log10 resistivity and ``model.weight`` as residual amplitude. This is
        opt-in because legacy files may use those columns differently.
        """
        if not isinstance(model, Occam1DModel):
            raise TypeError("model must be an Occam1DModel instance.")
        if not isinstance(use_native_preference, bool):
            raise TypeError("use_native_preference must be a bool.")
        model.validate()
        options = dict(kwargs)
        if use_native_preference:
            options["reference"] = model.preference
            options["reference_weights"] = model.weight
        return cls(
            model.n_layers,
            order=order,
            depth=model.depth,
            scaling=scaling,
            **options,
        )

    @staticmethod
    def _integer(name, value, *, minimum):
        """Return an integer at or above ``minimum``."""
        if not isinstance(value, int) or isinstance(value, bool):
            raise TypeError(f"{name} must be an integer.")
        if value < minimum:
            raise ValueError(f"{name} must be at least {minimum}.")
        return int(value)

    def _depth_vector(self, values):
        """Validate optional layer-top depths."""
        if values is None:
            if self.scaling == "depth":
                raise ValueError("depth is required for depth scaling.")
            return None
        result = self._vector(values, "depth", self.n_layers)
        if np.any(result < 0) or np.any(np.diff(result) <= 0):
            raise ValueError(
                "depth must be non-negative and increase strictly."
            )
        return result

    @staticmethod
    def _vector(values, name, size):
        """Return a finite vector of an exact size."""
        try:
            result = np.array(values, dtype=float, copy=True)
        except (TypeError, ValueError) as error:
            raise TypeError(f"{name} must contain real numbers.") from error
        if result.shape != (size,):
            raise ValueError(f"{name} must have shape ({size},).")
        if np.any(~np.isfinite(result)):
            raise ValueError(f"{name} must contain only finite values.")
        return result

    def _optional_vector(self, values, name):
        """Return an optional finite layer vector."""
        return (
            None
            if values is None
            else self._vector(values, name, self.n_layers)
        )

    def _weights(self, values, size, name, *, default):
        """Return finite non-negative residual amplitude weights."""
        if values is None:
            return np.full(size, default, dtype=float)
        result = self._vector(values, name, size)
        if np.any(result < 0):
            raise ValueError(f"{name} must be non-negative.")
        return result

    def _build_operator(self) -> np.ndarray:
        """Build the requested unweighted difference matrix."""
        count = self.n_layers - self.order
        operator = np.zeros((count, self.n_layers), dtype=float)
        if self.order == 1:
            rows = np.arange(count)
            operator[rows, rows] = -1.0
            operator[rows, rows + 1] = 1.0
            if self.scaling == "depth":
                spacing = np.diff(self._depth)
                operator /= np.sqrt(spacing)[:, None]
            return operator

        rows = np.arange(count)
        operator[rows, rows] = 1.0
        operator[rows, rows + 1] = -2.0
        operator[rows, rows + 2] = 1.0
        if self.scaling == "depth":
            left = np.diff(self._depth)[:-1]
            right = np.diff(self._depth)[1:]
            scale = np.sqrt(0.5 * (left + right))
            operator.fill(0.0)
            operator[rows, rows] = 2.0 / (left * (left + right))
            operator[rows, rows + 1] = -2.0 / (left * right)
            operator[rows, rows + 2] = 2.0 / (
                right * (left + right)
            )
            operator *= scale[:, None]
        return operator

    @property
    def operator(self) -> np.ndarray:
        """Independent smoothness matrix ``D``."""
        return self._operator.copy()

    @property
    def reference(self) -> np.ndarray | None:
        """Independent preferred log10-resistivity model."""
        return None if self._reference is None else self._reference.copy()

    @property
    def reference_weights(self) -> np.ndarray:
        """Independent reference residual amplitude weights."""
        return self._reference_weights.copy()

    @property
    def n_smoothness(self) -> int:
        """Number of roughness residuals."""
        return int(self._operator.shape[0])

    @property
    def n_reference(self) -> int:
        """Number of active preference residuals."""
        return int(np.count_nonzero(self._reference_weights > 0))

    def _parameters(self, values) -> np.ndarray:
        """Validate a layer log10-resistivity vector."""
        return self._vector(values, "parameters", self.n_layers)

    def system(self) -> tuple[np.ndarray, np.ndarray]:
        r"""Return structure matrix and right-hand side ``(R, q)``.

        The complete residual is ``R @ parameters - q``. Preference rows are
        appended only where their amplitude weight is positive.
        """
        matrix = self._operator.copy()
        rhs = np.zeros(self.n_smoothness, dtype=float)
        active = self._reference_weights > 0
        if np.any(active):
            preference = np.diag(self._reference_weights[active])
            columns = np.flatnonzero(active)
            expanded = np.zeros((columns.size, self.n_layers), dtype=float)
            expanded[:, columns] = preference
            matrix = np.vstack((matrix, expanded))
            rhs = np.r_[
                rhs,
                self._reference_weights[active] * self._reference[active],
            ]
        return matrix, rhs

    def evaluate(self, parameters) -> Occam1DRegularizationResult:
        """Evaluate roughness, preference penalty, and total objective."""
        parameters = self._parameters(parameters)
        smoothness = self._operator @ parameters
        active = self._reference_weights > 0
        reference = (
            self._reference_weights[active]
            * (parameters[active] - self._reference[active])
            if np.any(active)
            else np.empty(0, dtype=float)
        )
        roughness = float(smoothness @ smoothness)
        reference_penalty = float(reference @ reference)
        return Occam1DRegularizationResult(
            parameters=parameters,
            smoothness=smoothness,
            reference=reference,
            roughness=roughness,
            reference_penalty=reference_penalty,
            objective=roughness + reference_penalty,
            order=self.order,
            scaling=self.scaling,
        )

    @staticmethod
    def _linear_system(matrix, rhs, n_layers):
        """Validate a finite over- or under-determined data system."""
        try:
            design = np.array(matrix, dtype=float, copy=True)
            target = np.array(rhs, dtype=float, copy=True)
        except (TypeError, ValueError) as error:
            raise TypeError(
                "matrix and rhs must contain real values."
            ) from error
        if design.ndim != 2 or design.shape[1] != n_layers:
            raise ValueError(
                f"matrix must have shape (n_data, {n_layers})."
            )
        if not design.shape[0]:
            raise ValueError("matrix must contain at least one data row.")
        if target.shape != (design.shape[0],):
            raise ValueError("rhs must align with matrix rows.")
        if np.any(~np.isfinite(design)) or np.any(~np.isfinite(target)):
            raise ValueError("matrix and rhs must contain finite values.")
        return design, target

    def solve(
        self,
        matrix,
        rhs,
        multiplier,
        *,
        rcond=None,
        bounds=None,
        policy=None,
    ) -> Occam1DCandidate:
        r"""Solve one augmented least-squares Occam candidate.

        The solved objective is

        .. math::

           \|A m-b\|_2^2 + \mu\|R m-q\|_2^2.

        No normal-equation matrix or explicit inverse is formed. NumPy's SVD
        least-squares solver is used when ``bounds`` is omitted. With bounds,
        SciPy's trust-region reflective bounded least-squares solver minimizes
        the same augmented objective directly. The unconstrained answer is
        never merely clipped after solving.
        """
        design, target = self._linear_system(
            matrix,
            rhs,
            self.n_layers,
        )
        if isinstance(multiplier, bool) or not isinstance(
            multiplier, (int, float)
        ):
            raise TypeError("multiplier must be a non-negative real number.")
        multiplier = float(multiplier)
        if not math.isfinite(multiplier) or multiplier < 0:
            raise ValueError("multiplier must be finite and non-negative.")
        if rcond is not None:
            if isinstance(rcond, bool) or not isinstance(rcond, (int, float)):
                raise TypeError(
                    "rcond must be a positive real number or None."
                )
            if not math.isfinite(rcond) or rcond <= 0:
                raise ValueError("rcond must be finite and positive.")
        structure, structure_rhs = self.system()
        if multiplier > 0:
            root = math.sqrt(multiplier)
            augmented = np.vstack((design, root * structure))
            augmented_rhs = np.r_[target, root * structure_rhs]
        else:
            augmented = design
            augmented_rhs = target
        lower, upper = self._parameter_bounds(bounds)
        bounded = np.any(np.isfinite(lower)) or np.any(np.isfinite(upper))
        if policy is None:
            policy = Occam1DSolverPolicy()
        if not isinstance(policy, Occam1DSolverPolicy):
            raise TypeError("policy must be an Occam1DSolverPolicy or None.")
        working = augmented
        working_rhs = augmented_rhs
        damping = 0.0
        attempts = 0
        initial_condition = None
        while True:
            attempts += 1
            (
                parameters,
                rank,
                singular_values,
                condition,
                solver,
                active_lower,
                active_upper,
            ) = self._solve_augmented(
                working,
                working_rhs,
                lower,
                upper,
                bounded,
                rcond,
                policy,
            )
            if initial_condition is None:
                initial_condition = condition
            if condition <= policy.condition_limit:
                break
            if policy.ill_condition_action == "accept":
                break
            if policy.ill_condition_action == "reject":
                raise Occam1DLinearSolveError(
                    "Augmented system condition number "
                    f"{condition:.6g} exceeds policy limit "
                    f"{policy.condition_limit:.6g}."
                )
            if attempts > policy.max_retries:
                raise Occam1DLinearSolveError(
                    "Ill-conditioning persisted after "
                    f"{policy.max_retries} damping retries: "
                    f"condition={condition:.6g}, "
                    f"limit={policy.condition_limit:.6g}."
                )
            damping = policy.damping_start * (
                policy.damping_growth ** (attempts - 1)
            )
            working = np.vstack(
                (augmented, math.sqrt(damping) * np.eye(self.n_layers))
            )
            working_rhs = np.r_[augmented_rhs, np.zeros(self.n_layers)]
        if np.any(~np.isfinite(parameters)):
            raise FloatingPointError(
                "Regularized least-squares solve produced non-finite values."
            )
        data_residual = design @ parameters - target
        structure_residual = structure @ parameters - structure_rhs
        data_objective = float(data_residual @ data_residual)
        structure_objective = float(
            structure_residual @ structure_residual
        )
        objective = data_objective + multiplier * structure_objective
        return Occam1DCandidate(
            parameters=parameters,
            multiplier=multiplier,
            data_residual=data_residual,
            structure_residual=structure_residual,
            data_objective=data_objective,
            structure_objective=structure_objective,
            objective=objective,
            rank=int(rank),
            singular_values=singular_values,
            condition_number=condition,
            solver=solver,
            active_lower=active_lower,
            active_upper=active_upper,
            solve_attempts=attempts,
            damping=damping,
            initial_condition_number=initial_condition,
        )

    def _solve_augmented(
        self,
        matrix,
        rhs,
        lower,
        upper,
        bounded,
        rcond,
        policy,
    ):
        """Solve once and return complete conditioning diagnostics."""
        active_lower = np.zeros(self.n_layers, dtype=bool)
        active_upper = np.zeros(self.n_layers, dtype=bool)
        if bounded:
            from scipy.optimize import lsq_linear

            solution = lsq_linear(
                matrix,
                rhs,
                bounds=(lower, upper),
                method="trf",
                lsmr_tol="auto",
            )
            if not solution.success:
                raise Occam1DLinearSolveError(
                    "Bounded least-squares solve failed: "
                    f"{solution.message}"
                )
            parameters = solution.x
            active_lower = solution.active_mask == -1
            active_upper = solution.active_mask == 1
            solver = "scipy_lsq_linear_trf"
            singular_values = self._singular_values(matrix, policy)
            rank = self._rank(singular_values, matrix.shape)
        else:
            try:
                parameters, _, rank, singular_values = np.linalg.lstsq(
                    matrix,
                    rhs,
                    rcond=rcond,
                )
                solver = "svd_lstsq"
            except np.linalg.LinAlgError as error:
                if policy.svd_failure_action == "raise":
                    raise Occam1DLinearSolveError(
                        "NumPy SVD least-squares solve failed."
                    ) from error
                from scipy.linalg import lstsq

                try:
                    parameters, _, rank, _ = lstsq(
                        matrix,
                        rhs,
                        cond=rcond,
                        lapack_driver="gelsy",
                    )
                except Exception as fallback_error:
                    raise Occam1DLinearSolveError(
                        "Both NumPy SVD and SciPy pivoted-QR fallback "
                        "failed."
                    ) from fallback_error
                singular_values = self._singular_values(matrix, policy)
                solver = "scipy_gelsy_fallback"
        condition = self._condition(singular_values)
        return (
            parameters,
            int(rank),
            singular_values,
            condition,
            solver,
            active_lower,
            active_upper,
        )

    @staticmethod
    def _singular_values(matrix, policy):
        """Compute singular values with the configured independent fallback."""
        try:
            return np.linalg.svd(matrix, compute_uv=False)
        except np.linalg.LinAlgError as error:
            if policy.svd_failure_action == "raise":
                raise Occam1DLinearSolveError(
                    "SVD conditioning diagnostics failed."
                ) from error
            from scipy.linalg import svdvals

            try:
                return svdvals(matrix, check_finite=False)
            except Exception as fallback_error:
                raise Occam1DLinearSolveError(
                    "Primary and fallback SVD diagnostics failed."
                ) from fallback_error

    @staticmethod
    def _rank(singular_values, shape):
        """Return numerical matrix rank from precomputed singular values."""
        if not singular_values.size:
            return 0
        tolerance = (
            max(shape) * np.finfo(float).eps * singular_values[0]
        )
        return int(np.count_nonzero(singular_values > tolerance))

    @staticmethod
    def _condition(singular_values):
        """Return the two-norm condition number from singular values."""
        return (
            float(singular_values[0] / singular_values[-1])
            if singular_values.size
            and singular_values[-1] > np.finfo(float).tiny
            else float("inf")
        )

    def _parameter_bounds(self, bounds):
        """Return validated lower and upper log10 parameter bounds."""
        if bounds is None:
            return (
                np.full(self.n_layers, -np.inf),
                np.full(self.n_layers, np.inf),
            )
        if not isinstance(bounds, (tuple, list)) or len(bounds) != 2:
            raise TypeError("bounds must be a (lower, upper) pair or None.")
        output = []
        for name, values, default in (
            ("lower bounds", bounds[0], -np.inf),
            ("upper bounds", bounds[1], np.inf),
        ):
            if values is None:
                vector = np.full(self.n_layers, default)
            elif isinstance(values, bool):
                raise TypeError(f"{name} must contain real numbers.")
            elif isinstance(values, (int, float)):
                vector = np.full(self.n_layers, float(values))
            else:
                try:
                    vector = np.asarray(values, dtype=float)
                except (TypeError, ValueError) as error:
                    raise TypeError(
                        f"{name} must contain real numbers."
                    ) from error
                if vector.shape != (self.n_layers,):
                    raise ValueError(
                        f"{name} must have shape ({self.n_layers},)."
                    )
                vector = vector.copy()
            if np.any(np.isnan(vector)):
                raise ValueError(f"{name} cannot contain NaN.")
            output.append(vector)
        lower, upper = output
        if np.any(lower >= upper):
            raise ValueError(
                "Every lower parameter bound must be below its upper bound."
            )
        return lower, upper

    @staticmethod
    def observed_solver_vector(data: Occam1DData) -> np.ndarray:
        """Return observations in native row order and solver units."""
        if not isinstance(data, Occam1DData):
            raise TypeError("data must be an Occam1DData instance.")
        data.validate()
        values = []
        resistivity_mask = data.resistivity_mask
        phase_mask = data.phase_mask
        for index in range(data.n_frequencies):
            if resistivity_mask[index]:
                values.append(math.log10(data.resistivity[index]))
            if phase_mask[index]:
                values.append(data.phase[index])
        return np.asarray(values, dtype=float)

    def linearized_system(
        self,
        data,
        jacobian,
        current_parameters,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Build the absolute-model system for one nonlinear iteration.

        For current model ``m``, prediction ``F(m)``, and Jacobian ``J``, the
        right-hand side is ``d - F(m) + J @ m``. If the Jacobian result is
        weighted, ``d`` is divided by the same native standard errors.
        """
        if not isinstance(jacobian, Occam1DJacobianResult):
            raise TypeError(
                "jacobian must be an Occam1DJacobianResult instance."
            )
        if jacobian.n_layers != self.n_layers:
            raise ValueError(
                "Jacobian column count does not match regularization."
            )
        current = self._parameters(current_parameters)
        observed = self.observed_solver_vector(data)
        if observed.shape != (jacobian.n_data,):
            raise ValueError(
                "Observed data rows do not match the Jacobian result."
            )
        if jacobian.weighted:
            observed = observed / Occam1DJacobian.solver_errors(data)
        rhs = observed - jacobian.prediction + jacobian.matrix @ current
        return jacobian.matrix.copy(), rhs

    def solve_linearized(
        self,
        data,
        jacobian,
        current_parameters,
        multiplier,
        *,
        rcond=None,
        bounds=None,
        policy=None,
    ) -> Occam1DCandidate:
        """Build and solve one linearized absolute-model candidate."""
        matrix, rhs = self.linearized_system(
            data,
            jacobian,
            current_parameters,
        )
        return self.solve(
            matrix,
            rhs,
            multiplier,
            rcond=rcond,
            bounds=bounds,
            policy=policy,
        )

    def diagnostics(self) -> dict[str, Any]:
        """Extend lifecycle diagnostics with structure configuration."""
        values = super().diagnostics()
        values.update(
            {
                "n_layers": self.n_layers,
                "order": self.order,
                "scaling": self.scaling,
                "n_smoothness": self.n_smoothness,
                "n_reference": self.n_reference,
                "row_weights": self._row_weights.tolist(),
                "reference_weights": self._reference_weights.tolist(),
            }
        )
        return values
