# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Canonical analytic benchmarks for Maxwell backend validation.

The module defines executable benchmark cases rather than certifying any
backend by name. A backend becomes verified only when stored benchmark
outcomes identify its exact version, problem hashes, thresholds, and metrics.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

import numpy as np

from .backends import MaxwellBackend
from .contracts import ForwardResult, MaxwellMesh, MaxwellProblem, ReceiverSet

__all__ = [
    "BenchmarkThresholds",
    "BenchmarkMetrics",
    "BenchmarkOutcome",
    "MaxwellBenchmark",
    "BenchmarkReport",
    "half_space_impedance",
    "layered_earth_impedance",
    "half_space_benchmark",
    "layered_earth_benchmark",
    "run_benchmarks",
]

_MU0 = 4.0e-7 * np.pi
_COMPONENTS = ("zxx", "zxy", "zyx", "zyy")


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def _json_mapping(value: Mapping[str, Any]) -> Mapping[str, Any]:
    try:
        encoded = json.dumps(dict(value), sort_keys=True, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "metadata must contain finite JSON-compatible values."
        ) from exc
    return MappingProxyType(json.loads(encoded))


def _positive(value: float, name: str, *, allow_zero: bool = False) -> float:
    result = float(value)
    lower_ok = result >= 0 if allow_zero else result > 0
    if not np.isfinite(result) or not lower_ok:
        qualifier = "non-negative" if allow_zero else "positive"
        raise ValueError(f"{name} must be finite and {qualifier}.")
    return result


@dataclass(frozen=True)
class BenchmarkThresholds:
    """Define quantitative acceptance limits for one benchmark.

    Parameters
    ----------
    maximum_normalized_rms : float, default=0.05
        Maximum complex root-sum-square error normalized by the reference.
    maximum_amplitude_relative_error : float, default=0.05
        Maximum pointwise relative impedance-amplitude error.
    maximum_phase_error_deg : float, default=2.0
        Maximum absolute circular phase error in degrees.
    minimum_valid_fraction : float, default=1.0
        Minimum fraction of output values marked valid.
    require_convergence : bool, default=True
        Require every solve in backend diagnostics to converge.

    Examples
    --------
    >>> limits = BenchmarkThresholds(maximum_phase_error_deg=1)
    >>> limits.maximum_phase_error_deg
    1.0
    """

    maximum_normalized_rms: float = 0.05
    maximum_amplitude_relative_error: float = 0.05
    maximum_phase_error_deg: float = 2.0
    minimum_valid_fraction: float = 1.0
    require_convergence: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "maximum_normalized_rms",
            _positive(
                self.maximum_normalized_rms,
                "maximum_normalized_rms",
                allow_zero=True,
            ),
        )
        object.__setattr__(
            self,
            "maximum_amplitude_relative_error",
            _positive(
                self.maximum_amplitude_relative_error,
                "maximum_amplitude_relative_error",
                allow_zero=True,
            ),
        )
        object.__setattr__(
            self,
            "maximum_phase_error_deg",
            _positive(
                self.maximum_phase_error_deg,
                "maximum_phase_error_deg",
                allow_zero=True,
            ),
        )
        valid = float(self.minimum_valid_fraction)
        if not np.isfinite(valid) or valid < 0 or valid > 1:
            raise ValueError("minimum_valid_fraction must be in [0, 1].")
        object.__setattr__(self, "minimum_valid_fraction", valid)

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible acceptance limits.

        Returns
        -------
        dict
            Versioned threshold state.

        Examples
        --------
        >>> BenchmarkThresholds().to_dict()["schema_version"]
        1
        """
        return {
            "schema_version": 1,
            "maximum_normalized_rms": self.maximum_normalized_rms,
            "maximum_amplitude_relative_error": (self.maximum_amplitude_relative_error),
            "maximum_phase_error_deg": self.maximum_phase_error_deg,
            "minimum_valid_fraction": self.minimum_valid_fraction,
            "require_convergence": self.require_convergence,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> BenchmarkThresholds:
        """Restore validated benchmark thresholds.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        BenchmarkThresholds
            Restored limits.

        Examples
        --------
        >>> limits = BenchmarkThresholds(maximum_normalized_rms=0.1)
        >>> BenchmarkThresholds.from_dict(limits.to_dict()) == limits
        True
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported BenchmarkThresholds schema version.")
        values = dict(data)
        values.pop("schema_version")
        return cls(**values)


@dataclass(frozen=True)
class BenchmarkMetrics:
    """Store errors measured against one analytic reference.

    Parameters
    ----------
    normalized_rms : float
        Complex normalized root-sum-square error.
    maximum_amplitude_relative_error : float
        Worst pointwise relative amplitude error.
    maximum_phase_error_deg : float
        Worst circular phase difference.
    valid_fraction : float
        Fraction of output values marked valid and finite.
    converged : bool
        Whether all backend solves converged.

    Examples
    --------
    >>> BenchmarkMetrics(0.01, 0.02, 0.5, 1, True).converged
    True
    """

    normalized_rms: float
    maximum_amplitude_relative_error: float
    maximum_phase_error_deg: float
    valid_fraction: float
    converged: bool

    def __post_init__(self) -> None:
        names = (
            "normalized_rms",
            "maximum_amplitude_relative_error",
            "maximum_phase_error_deg",
        )
        for name in names:
            value = _positive(getattr(self, name), name, allow_zero=True)
            object.__setattr__(self, name, value)
        valid = float(self.valid_fraction)
        if not np.isfinite(valid) or valid < 0 or valid > 1:
            raise ValueError("valid_fraction must be in [0, 1].")
        object.__setattr__(self, "valid_fraction", valid)
        object.__setattr__(self, "converged", bool(self.converged))

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible benchmark metrics.

        Returns
        -------
        dict
            Error values, validity, and convergence status.

        Examples
        --------
        >>> BenchmarkMetrics(0, 0, 0, 1, True).to_dict()["valid_fraction"]
        1.0
        """
        return {
            "normalized_rms": self.normalized_rms,
            "maximum_amplitude_relative_error": (self.maximum_amplitude_relative_error),
            "maximum_phase_error_deg": self.maximum_phase_error_deg,
            "valid_fraction": self.valid_fraction,
            "converged": self.converged,
        }


@dataclass(frozen=True)
class BenchmarkOutcome:
    """Record the auditable outcome of one backend benchmark.

    Parameters
    ----------
    benchmark_name, benchmark_hash : str
        Stable case identity and full content digest.
    backend_name, backend_version : str
        Exact adapter identity from the result.
    passed : bool
        Whether every configured acceptance criterion passed.
    metrics : BenchmarkMetrics
        Quantitative comparison with the reference.
    failures : tuple of str
        Human-readable failed criteria.

    Examples
    --------
    >>> metrics = BenchmarkMetrics(0, 0, 0, 1, True)
    >>> outcome = BenchmarkOutcome(
    ...     "half-space", "0" * 64, "demo", "1", True, metrics
    ... )
    >>> outcome.passed
    True
    """

    benchmark_name: str
    benchmark_hash: str
    backend_name: str
    backend_version: str
    passed: bool
    metrics: BenchmarkMetrics
    failures: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        name = str(self.benchmark_name).strip()
        backend = str(self.backend_name).strip()
        version = str(self.backend_version).strip()
        digest = str(self.benchmark_hash).strip().lower()
        if not name or not backend or not version:
            raise ValueError("benchmark and backend identity cannot be empty.")
        if len(digest) != 64 or any(
            value not in "0123456789abcdef" for value in digest
        ):
            raise ValueError("benchmark_hash must be a SHA-256 digest.")
        if not isinstance(self.metrics, BenchmarkMetrics):
            raise TypeError("metrics must be BenchmarkMetrics.")
        failures = tuple(str(value).strip() for value in self.failures)
        if any(not value for value in failures):
            raise ValueError("failures cannot contain empty messages.")
        if bool(self.passed) == bool(failures):
            raise ValueError("passed must be true exactly when failures is empty.")
        object.__setattr__(self, "benchmark_name", name)
        object.__setattr__(self, "benchmark_hash", digest)
        object.__setattr__(self, "backend_name", backend)
        object.__setattr__(self, "backend_version", version)
        object.__setattr__(self, "passed", bool(self.passed))
        object.__setattr__(self, "failures", failures)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible benchmark outcome.

        Returns
        -------
        dict
            Case identity, backend identity, metrics, and failures.

        Examples
        --------
        >>> metrics = BenchmarkMetrics(0, 0, 0, 1, True)
        >>> value = BenchmarkOutcome(
        ...     "case", "0" * 64, "demo", "1", True, metrics
        ... )
        >>> value.to_dict()["passed"]
        True
        """
        return {
            "benchmark_name": self.benchmark_name,
            "benchmark_hash": self.benchmark_hash,
            "backend_name": self.backend_name,
            "backend_version": self.backend_version,
            "passed": self.passed,
            "metrics": self.metrics.to_dict(),
            "failures": list(self.failures),
        }


@dataclass(frozen=True)
class MaxwellBenchmark:
    """Define one immutable problem and its expected impedance.

    Parameters
    ----------
    name, description : str
        Stable identifier and scientific purpose.
    problem : MaxwellProblem
        Exact solver input.
    reference_impedance_v_a : complex ndarray
        Analytic reference with canonical problem output shape.
    thresholds : BenchmarkThresholds, optional
        Quantitative acceptance criteria.
    tags : sequence of str, optional
        Searchable labels such as ``analytic`` and ``half-space``.
    metadata : mapping, optional
        Finite JSON-compatible provenance.

    Examples
    --------
    Cases are normally built with :func:`half_space_benchmark` or
    :func:`layered_earth_benchmark`.
    """

    name: str
    description: str
    problem: MaxwellProblem
    reference_impedance_v_a: np.ndarray
    thresholds: BenchmarkThresholds = BenchmarkThresholds()
    tags: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(
        default_factory=lambda: MappingProxyType({})
    )

    def __post_init__(self) -> None:
        name = str(self.name).strip().lower().replace("_", "-")
        description = str(self.description).strip()
        if not name or not description:
            raise ValueError("name and description cannot be empty.")
        if not isinstance(self.problem, MaxwellProblem):
            raise TypeError("problem must be a MaxwellProblem.")
        if not isinstance(self.thresholds, BenchmarkThresholds):
            raise TypeError("thresholds must be BenchmarkThresholds.")
        reference = np.asarray(self.reference_impedance_v_a, dtype=complex)
        expected = (
            self.problem.receivers.count,
            len(self.problem.frequencies_hz),
            len(self.problem.components),
        )
        if reference.shape != expected:
            raise ValueError(f"reference_impedance_v_a must have shape {expected}.")
        if not np.all(np.isfinite(reference)):
            raise ValueError("reference_impedance_v_a must be finite.")
        if np.any(np.abs(reference) == 0):
            raise ValueError("reference impedance magnitudes must be nonzero.")
        tags = tuple(str(value).strip().lower() for value in self.tags)
        if len(set(tags)) != len(tags) or any(not value for value in tags):
            raise ValueError("tags must contain unique non-empty values.")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "description", description)
        object.__setattr__(self, "reference_impedance_v_a", _readonly(reference))
        object.__setattr__(self, "tags", tags)
        object.__setattr__(self, "metadata", _json_mapping(self.metadata))

    @property
    def benchmark_hash(self) -> str:
        """Return a deterministic digest of case inputs and thresholds.

        Returns
        -------
        str
            SHA-256 benchmark identity.

        Examples
        --------
        Benchmark hashes contain 64 hexadecimal characters.
        """
        digest = hashlib.sha256()
        reference = np.ascontiguousarray(
            self.reference_impedance_v_a,
            dtype="<c16",
        )
        digest.update(reference.tobytes())
        digest.update(
            json.dumps(
                self.provenance(),
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        )
        return digest.hexdigest()

    def evaluate(self, result: ForwardResult) -> BenchmarkOutcome:
        """Compare one canonical result with the analytic reference.

        Parameters
        ----------
        result : ForwardResult
            Backend result for this exact problem.

        Returns
        -------
        BenchmarkOutcome
            Metrics and every failed acceptance criterion.

        Raises
        ------
        ValueError
            If the result belongs to a different problem or output axes.

        Examples
        --------
        Exact references produce zero error when wrapped as backend results.
        """
        if not isinstance(result, ForwardResult):
            raise TypeError("result must be a ForwardResult.")
        result.validate_against(self.problem)
        finite = np.isfinite(result.impedance_v_a.real)
        finite &= np.isfinite(result.impedance_v_a.imag)
        valid = result.valid & finite
        valid_fraction = float(np.mean(valid))
        if np.any(valid):
            actual = result.impedance_v_a[valid]
            expected = self.reference_impedance_v_a[valid]
            difference = actual - expected
            denominator = float(np.sum(np.abs(expected) ** 2))
            normalized_rms = float(
                np.sqrt(np.sum(np.abs(difference) ** 2) / denominator)
            )
            amplitude_error = np.abs(np.abs(actual) - np.abs(expected))
            amplitude_error /= np.abs(expected)
            maximum_amplitude = float(np.max(amplitude_error))
            phase_error = np.angle(actual / expected, deg=True)
            maximum_phase = float(np.max(np.abs(phase_error)))
        else:
            unavailable = float(np.finfo(float).max)
            normalized_rms = unavailable
            maximum_amplitude = unavailable
            maximum_phase = unavailable
        metrics = BenchmarkMetrics(
            normalized_rms,
            maximum_amplitude,
            maximum_phase,
            valid_fraction,
            result.diagnostics.success,
        )
        failures = self._failures(metrics)
        return BenchmarkOutcome(
            self.name,
            self.benchmark_hash,
            result.backend_name,
            result.backend_version,
            not failures,
            metrics,
            failures,
        )

    def run(self, backend: MaxwellBackend) -> BenchmarkOutcome:
        """Execute and evaluate this case with a conforming backend.

        Parameters
        ----------
        backend : MaxwellBackend
            Backend compatible with the benchmark problem.

        Returns
        -------
        BenchmarkOutcome
            Auditable validation outcome.

        Examples
        --------
        Backend exceptions propagate so infrastructure failures cannot be
        mistaken for numerical benchmark failures.
        """
        if not isinstance(backend, MaxwellBackend):
            raise TypeError("backend must implement MaxwellBackend.")
        return self.evaluate(backend.solve(self.problem))

    def provenance(self) -> dict[str, Any]:
        """Return JSON-compatible benchmark provenance.

        Returns
        -------
        dict
            Case identity, problem hash, limits, tags, and metadata.

        Examples
        --------
        The full conductivity model remains identified by ``problem_hash``.
        """
        return {
            "schema_version": 1,
            "name": self.name,
            "description": self.description,
            "problem_hash": self.problem.problem_hash,
            "thresholds": self.thresholds.to_dict(),
            "tags": list(self.tags),
            "metadata": dict(self.metadata),
        }

    def _failures(self, metrics: BenchmarkMetrics) -> tuple[str, ...]:
        failures = []
        limits = self.thresholds
        if metrics.normalized_rms > limits.maximum_normalized_rms:
            failures.append(
                "normalized RMS "
                f"{metrics.normalized_rms:.6g} exceeds "
                f"{limits.maximum_normalized_rms:.6g}"
            )
        maximum_amplitude = metrics.maximum_amplitude_relative_error
        if maximum_amplitude > limits.maximum_amplitude_relative_error:
            failures.append(
                "maximum amplitude relative error "
                f"{maximum_amplitude:.6g} exceeds "
                f"{limits.maximum_amplitude_relative_error:.6g}"
            )
        if metrics.maximum_phase_error_deg > limits.maximum_phase_error_deg:
            failures.append(
                "maximum phase error "
                f"{metrics.maximum_phase_error_deg:.6g} deg exceeds "
                f"{limits.maximum_phase_error_deg:.6g} deg"
            )
        if metrics.valid_fraction < limits.minimum_valid_fraction:
            failures.append(
                f"valid fraction {metrics.valid_fraction:.6g} is below "
                f"{limits.minimum_valid_fraction:.6g}"
            )
        if limits.require_convergence and not metrics.converged:
            failures.append("one or more solver systems did not converge")
        return tuple(failures)


@dataclass(frozen=True)
class BenchmarkReport:
    """Aggregate ordered outcomes from one backend benchmark run.

    Parameters
    ----------
    outcomes : sequence of BenchmarkOutcome
        Non-empty ordered outcomes from one backend version.

    Examples
    --------
    >>> metrics = BenchmarkMetrics(0, 0, 0, 1, True)
    >>> outcome = BenchmarkOutcome(
    ...     "case", "0" * 64, "demo", "1", True, metrics
    ... )
    >>> BenchmarkReport((outcome,)).passed
    True
    """

    outcomes: tuple[BenchmarkOutcome, ...]

    def __post_init__(self) -> None:
        outcomes = tuple(self.outcomes)
        if not outcomes or any(
            not isinstance(value, BenchmarkOutcome) for value in outcomes
        ):
            raise ValueError("outcomes must contain BenchmarkOutcome values.")
        identities = {(value.backend_name, value.backend_version) for value in outcomes}
        if len(identities) != 1:
            raise ValueError("all outcomes must use one backend version.")
        names = [value.benchmark_name for value in outcomes]
        if len(set(names)) != len(names):
            raise ValueError("benchmark names must be unique in a report.")
        object.__setattr__(self, "outcomes", outcomes)

    @property
    def passed(self) -> bool:
        """Return whether every benchmark passed.

        Returns
        -------
        bool
            Aggregate acceptance status.

        Examples
        --------
        A report fails if any contained outcome fails.
        """
        return all(value.passed for value in self.outcomes)

    @property
    def pass_fraction(self) -> float:
        """Return the fraction of cases that passed.

        Returns
        -------
        float
            Passed case count divided by total cases.

        Examples
        --------
        The value lies in the closed interval ``[0, 1]``.
        """
        return sum(value.passed for value in self.outcomes) / len(self.outcomes)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible benchmark report.

        Returns
        -------
        dict
            Backend identity, summary, and ordered outcomes.

        Examples
        --------
        Serialized reports retain every individual failure message.
        """
        first = self.outcomes[0]
        return {
            "schema_version": 1,
            "backend_name": first.backend_name,
            "backend_version": first.backend_version,
            "passed": self.passed,
            "pass_fraction": self.pass_fraction,
            "outcomes": [value.to_dict() for value in self.outcomes],
        }


def half_space_impedance(
    resistivity_ohm_m: float,
    frequencies_hz: Any,
    *,
    time_dependence: str = "exp(+iwt)",
) -> np.ndarray:
    """Return analytic plane-wave impedance of a uniform half-space.

    Parameters
    ----------
    resistivity_ohm_m : float
        Positive half-space resistivity.
    frequencies_hz : array-like
        Positive frequencies.
    time_dependence : {"exp(+iwt)", "exp(-iwt)"}
        Complex phasor convention.

    Returns
    -------
    ndarray of complex
        ``sqrt(i*omega*mu0*rho)`` or its conjugate convention.

    Examples
    --------
    >>> round(float(np.angle(half_space_impedance(100, 1), deg=True)))
    45
    """
    resistivity = _positive(resistivity_ohm_m, "resistivity_ohm_m")
    frequencies = np.asarray(frequencies_hz, dtype=float)
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0):
        raise ValueError("frequencies_hz must be positive and finite.")
    if time_dependence not in {"exp(+iwt)", "exp(-iwt)"}:
        raise ValueError("unsupported time_dependence.")
    impedance = np.sqrt(1j * 2 * np.pi * frequencies * _MU0 * resistivity)
    if time_dependence == "exp(-iwt)":
        impedance = np.conjugate(impedance)
    return _readonly(impedance)


def layered_earth_impedance(
    resistivity_ohm_m: Sequence[float],
    thickness_m: Sequence[float],
    frequencies_hz: Any,
    *,
    time_dependence: str = "exp(+iwt)",
) -> np.ndarray:
    """Return analytic 1-D MT impedance by upward layer recursion.

    Parameters
    ----------
    resistivity_ohm_m : sequence of float
        Layer resistivities from surface to basal half-space.
    thickness_m : sequence of float
        Thickness of every layer except the basal half-space.
    frequencies_hz : array-like
        Positive evaluation frequencies.
    time_dependence : {"exp(+iwt)", "exp(-iwt)"}
        Complex phasor convention.

    Returns
    -------
    ndarray of complex
        Surface impedance in frequency order.

    Examples
    --------
    A single layer reduces exactly to the half-space expression:

    >>> np.allclose(
    ...     layered_earth_impedance([100], [], [1, 10]),
    ...     half_space_impedance(100, [1, 10]),
    ... )
    True
    """
    resistivity = np.asarray(resistivity_ohm_m, dtype=float)
    thickness = np.asarray(thickness_m, dtype=float)
    frequencies = np.asarray(frequencies_hz, dtype=float)
    if resistivity.ndim != 1 or len(resistivity) < 1:
        raise ValueError("resistivity_ohm_m must be a non-empty vector.")
    if thickness.shape != (len(resistivity) - 1,):
        raise ValueError("thickness_m must omit only the basal half-space.")
    if not np.all(np.isfinite(resistivity)) or np.any(resistivity <= 0):
        raise ValueError("resistivity_ohm_m must be positive and finite.")
    if not np.all(np.isfinite(thickness)) or np.any(thickness <= 0):
        raise ValueError("thickness_m must be positive and finite.")
    if frequencies.ndim != 1 or len(frequencies) < 1:
        raise ValueError("frequencies_hz must be a non-empty vector.")
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0):
        raise ValueError("frequencies_hz must be positive and finite.")
    if time_dependence not in {"exp(+iwt)", "exp(-iwt)"}:
        raise ValueError("unsupported time_dependence.")
    output = np.empty(len(frequencies), dtype=complex)
    for index, frequency in enumerate(frequencies):
        omega = 2 * np.pi * frequency
        wave_number = np.sqrt(1j * omega * _MU0 / resistivity[-1])
        impedance = 1j * omega * _MU0 / wave_number
        for layer in range(len(resistivity) - 2, -1, -1):
            wave_number = np.sqrt(1j * omega * _MU0 / resistivity[layer])
            intrinsic = 1j * omega * _MU0 / wave_number
            tangent = np.tanh(wave_number * thickness[layer])
            numerator = impedance + intrinsic * tangent
            denominator = intrinsic + impedance * tangent
            impedance = intrinsic * numerator / denominator
        output[index] = impedance
    if time_dependence == "exp(-iwt)":
        output = np.conjugate(output)
    return _readonly(output)


def half_space_benchmark(
    mesh: MaxwellMesh,
    receivers: ReceiverSet,
    frequencies_hz: Sequence[float],
    *,
    resistivity_ohm_m: float = 100.0,
    components: Sequence[str] = ("zxy", "zyx"),
    time_dependence: str = "exp(+iwt)",
    thresholds: BenchmarkThresholds | None = None,
) -> MaxwellBenchmark:
    """Build a uniform-earth analytic benchmark.

    Parameters
    ----------
    mesh, receivers : MaxwellMesh, ReceiverSet
        Solver geometry and observation locations.
    frequencies_hz : sequence of float
        Positive benchmark frequencies.
    resistivity_ohm_m : float, default=100
        Uniform earth resistivity.
    components : sequence of str, default=("zxy", "zyx")
        Requested components. Diagonal 3-D components are excluded because
        their analytic reference is zero and relative metrics are undefined.
    time_dependence : str, default="exp(+iwt)"
        Complex phasor convention.
    thresholds : BenchmarkThresholds or None, optional
        Acceptance limits.

    Returns
    -------
    MaxwellBenchmark
        Executable half-space case.

    Examples
    --------
    >>> mesh = MaxwellMesh([0,1,2], [0,1,2])
    >>> receivers = ReceiverSet([[0.5,0]], ["S"])
    >>> half_space_benchmark(mesh, receivers, [1]).name
    'half-space'
    """
    resistivity = _positive(resistivity_ohm_m, "resistivity_ohm_m")
    selected = _off_diagonal_components(components)
    problem = MaxwellProblem(
        mesh,
        np.full(mesh.shape, 1.0 / resistivity),
        frequencies_hz,
        receivers,
        selected,
        time_dependence=time_dependence,
        metadata={"benchmark": "half-space"},
    )
    base = half_space_impedance(
        resistivity,
        problem.frequencies_hz,
        time_dependence=time_dependence,
    )
    reference = _tensor_reference(problem, base)
    return MaxwellBenchmark(
        "half-space",
        "Uniform isotropic earth analytic MT limit.",
        problem,
        reference,
        BenchmarkThresholds() if thresholds is None else thresholds,
        ("analytic", "half-space", f"{mesh.dimension}d"),
        {"resistivity_ohm_m": resistivity},
    )


def layered_earth_benchmark(
    mesh: MaxwellMesh,
    receivers: ReceiverSet,
    frequencies_hz: Sequence[float],
    resistivity_ohm_m: Sequence[float],
    thickness_m: Sequence[float],
    *,
    components: Sequence[str] = ("zxy", "zyx"),
    time_dependence: str = "exp(+iwt)",
    thresholds: BenchmarkThresholds | None = None,
) -> MaxwellBenchmark:
    """Build a laterally uniform layered-earth benchmark.

    Parameters
    ----------
    mesh, receivers, frequencies_hz
        Solver geometry, observations, and frequencies.
    resistivity_ohm_m : sequence of float
        Layer resistivities ending with a basal half-space.
    thickness_m : sequence of float
        Finite-layer thicknesses. Every cumulative interface must coincide
        with a mesh z edge.
    components, time_dependence, thresholds
        Output components, phasor convention, and acceptance limits.

    Returns
    -------
    MaxwellBenchmark
        Executable analytic layered-earth case.

    Examples
    --------
    >>> mesh = MaxwellMesh([0,1,2], [0,1,2])
    >>> receivers = ReceiverSet([[0.5,0]], ["S"])
    >>> case = layered_earth_benchmark(
    ...     mesh, receivers, [1], [10, 100], [1]
    ... )
    >>> case.name
    'layered-earth'
    """
    resistivity = np.asarray(resistivity_ohm_m, dtype=float)
    thickness = np.asarray(thickness_m, dtype=float)
    base = layered_earth_impedance(
        resistivity,
        thickness,
        frequencies_hz,
        time_dependence=time_dependence,
    )
    if not np.isclose(mesh.z_edges_m[0], 0, atol=1e-12, rtol=0):
        raise ValueError("mesh z edges must begin at depth zero.")
    interfaces = np.cumsum(thickness)
    for interface in interfaces:
        tolerance = max(1e-10, abs(interface) * 1e-10)
        if not np.any(np.isclose(mesh.z_edges_m, interface, atol=tolerance)):
            raise ValueError(f"layer interface {interface:g} m is not a mesh z edge.")
    layer_index = np.searchsorted(
        interfaces,
        mesh.cell_centres_m["z"],
        side="right",
    )
    vertical = 1.0 / resistivity[layer_index]
    shape = (len(vertical),) + (1,) * (mesh.dimension - 1)
    conductivity = np.broadcast_to(vertical.reshape(shape), mesh.shape)
    selected = _off_diagonal_components(components)
    problem = MaxwellProblem(
        mesh,
        conductivity,
        frequencies_hz,
        receivers,
        selected,
        time_dependence=time_dependence,
        metadata={"benchmark": "layered-earth"},
    )
    reference = _tensor_reference(problem, base)
    return MaxwellBenchmark(
        "layered-earth",
        "Laterally uniform layered-earth analytic MT limit.",
        problem,
        reference,
        BenchmarkThresholds() if thresholds is None else thresholds,
        ("analytic", "layered-earth", f"{mesh.dimension}d"),
        {
            "resistivity_ohm_m": resistivity.tolist(),
            "thickness_m": thickness.tolist(),
        },
    )


def run_benchmarks(
    backend: MaxwellBackend,
    benchmarks: Sequence[MaxwellBenchmark],
) -> BenchmarkReport:
    """Run an ordered benchmark collection with one backend version.

    Parameters
    ----------
    backend : MaxwellBackend
        Backend under validation.
    benchmarks : sequence of MaxwellBenchmark
        Non-empty cases with unique names.

    Returns
    -------
    BenchmarkReport
        Aggregate and per-case outcomes.

    Examples
    --------
    Backend and numerical exceptions propagate rather than becoming false
    benchmark failures.
    """
    if not isinstance(backend, MaxwellBackend):
        raise TypeError("backend must implement MaxwellBackend.")
    cases = tuple(benchmarks)
    if not cases or any(not isinstance(value, MaxwellBenchmark) for value in cases):
        raise ValueError("benchmarks must contain MaxwellBenchmark values.")
    names = [value.name for value in cases]
    if len(set(names)) != len(names):
        raise ValueError("benchmark names must be unique.")
    return BenchmarkReport(tuple(value.run(backend) for value in cases))


def _off_diagonal_components(values: Sequence[str]) -> tuple[str, ...]:
    components = tuple(str(value).strip().lower() for value in values)
    if not components or len(set(components)) != len(components):
        raise ValueError("components must be non-empty and unique.")
    if any(value not in _COMPONENTS for value in components):
        raise ValueError("components contain an unknown impedance name.")
    if any(value in {"zxx", "zyy"} for value in components):
        raise ValueError("analytic relative metrics require zxy and/or zyx.")
    return components


def _tensor_reference(
    problem: MaxwellProblem,
    base: np.ndarray,
) -> np.ndarray:
    values = np.empty(
        (
            problem.receivers.count,
            len(problem.frequencies_hz),
            len(problem.components),
        ),
        dtype=complex,
    )
    for index, component in enumerate(problem.components):
        sign = 1.0 if component == "zxy" else -1.0
        values[:, :, index] = sign * base[None, :]
    return values
