# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Capability checks and lazy registration for Maxwell solver backends.

Backends are integrations, not deep-learning frameworks.  Each adapter must
declare its physical and numerical scope before it can receive a problem.
Factories are stored lazily so optional solver packages are imported only when
a caller explicitly creates that backend.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from threading import RLock
from types import MappingProxyType
from typing import (
    Any,
    Callable,
    Optional,
    Protocol,
    runtime_checkable,
)

import numpy as np

from .contracts import ForwardResult, MaxwellProblem

__all__ = [
    "BackendCapabilities",
    "CompatibilityReport",
    "MaxwellBackend",
    "BackendRegistration",
    "BackendRegistry",
    "backend_registry",
    "register_backend",
    "unregister_backend",
    "create_backend",
    "list_backends",
]

_COMPONENTS = ("zxx", "zxy", "zyx", "zyy")
_TIME_CONVENTIONS = ("exp(+iwt)", "exp(-iwt)")


def _identifier(value: str, label: str) -> str:
    result = str(value).strip().lower().replace("_", "-")
    if not result or any(
        not (character.isalnum() or character in "-.") for character in result
    ):
        raise ValueError(
            f"{label} must use letters, numbers, hyphens, or dots."
        )
    return result


def _positive_limit(value: int | None, label: str) -> int | None:
    if value is None:
        return None
    if (
        not isinstance(value, (int, np.integer))
        or isinstance(value, bool)
        or int(value) < 1
    ):
        raise ValueError(f"{label} must be a positive integer or None.")
    return int(value)


@dataclass(frozen=True)
class CompatibilityReport:
    """Describe whether a backend can solve a particular problem.

    Parameters
    ----------
    backend_name : str
        Normalized backend identifier.
    compatible : bool
        Whether all hard capability requirements are satisfied.
    errors, warnings : tuple of str, optional
        Hard incompatibilities and advisory concerns.

    Examples
    --------
    >>> report = CompatibilityReport("demo", False, ("3-D unsupported",))
    >>> report.require()
    Traceback (most recent call last):
    ...
    ValueError: backend 'demo' is incompatible: 3-D unsupported
    """

    backend_name: str
    compatible: bool
    errors: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        name = _identifier(self.backend_name, "backend_name")
        errors = tuple(str(value).strip() for value in self.errors)
        warnings = tuple(str(value).strip() for value in self.warnings)
        if any(not value for value in errors + warnings):
            raise ValueError("compatibility messages cannot be empty.")
        if bool(self.compatible) == bool(errors):
            raise ValueError(
                "compatible must be true exactly when errors is empty."
            )
        object.__setattr__(self, "backend_name", name)
        object.__setattr__(self, "compatible", bool(self.compatible))
        object.__setattr__(self, "errors", errors)
        object.__setattr__(self, "warnings", warnings)

    def require(self) -> None:
        """Raise a consolidated error when the report is incompatible.

        Raises
        ------
        ValueError
            If one or more hard incompatibilities were found.

        Examples
        --------
        >>> CompatibilityReport("demo", True).require()
        """
        if self.errors:
            raise ValueError(
                f"backend {self.backend_name!r} is incompatible: "
                + "; ".join(self.errors)
            )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible report.

        Returns
        -------
        dict
            Compatibility state and messages.

        Examples
        --------
        >>> CompatibilityReport("demo", True).to_dict()["compatible"]
        True
        """
        return {
            "backend_name": self.backend_name,
            "compatible": self.compatible,
            "errors": list(self.errors),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class BackendCapabilities:
    """Declare a Maxwell adapter's supported physical and numerical scope.

    Parameters
    ----------
    name, version : str
        Stable backend identifier and adapter/solver version.
    dimensions : tuple containing 2 and/or 3
        Spatial problem dimensions genuinely supported.
    components : tuple of str
        Impedance tensor components the adapter can return.
    time_conventions : tuple of str, default=("exp(+iwt)",)
        Phasor conventions accepted without conversion.
    supports_nonuniform_mesh : bool, default=True
        Whether variable cell widths are supported.
    supports_inactive_cells : bool, default=False
        Whether ``MaxwellProblem.active_cells`` is honored.
    supports_topography : bool, default=False
        Whether an inactive/air mask may describe non-flat terrain.
    supports_anisotropy : bool, default=False
        Reserved declaration for future tensor-conductivity contracts.
    maximum_cells, maximum_frequencies : int or None, optional
        Enforced adapter limits. None means no declared limit.
    verified_benchmarks : tuple of str, optional
        Stable benchmark identifiers passed by this adapter version.

    Examples
    --------
    >>> capability = BackendCapabilities("mt2d", "1.0", (2,), ("zxy", "zyx"))
    >>> capability.supports_dimension(2), capability.supports_component("zxy")
    (True, True)
    """

    name: str
    version: str
    dimensions: tuple[int, ...]
    components: tuple[str, ...]
    time_conventions: tuple[str, ...] = ("exp(+iwt)",)
    supports_nonuniform_mesh: bool = True
    supports_inactive_cells: bool = False
    supports_topography: bool = False
    supports_anisotropy: bool = False
    maximum_cells: int | None = None
    maximum_frequencies: int | None = None
    verified_benchmarks: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        name = _identifier(self.name, "name")
        version = str(self.version).strip()
        if not version:
            raise ValueError("version cannot be empty.")
        dimensions = tuple(int(value) for value in self.dimensions)
        if (
            not dimensions
            or len(set(dimensions)) != len(dimensions)
            or any(value not in (2, 3) for value in dimensions)
        ):
            raise ValueError(
                "dimensions must contain unique values drawn from (2, 3)."
            )
        components = tuple(
            str(value).strip().lower() for value in self.components
        )
        if (
            not components
            or len(set(components)) != len(components)
            or any(value not in _COMPONENTS for value in components)
        ):
            raise ValueError(
                f"components must contain unique values drawn from {_COMPONENTS}."
            )
        conventions = tuple(
            str(value).strip() for value in self.time_conventions
        )
        if (
            not conventions
            or len(set(conventions)) != len(conventions)
            or any(value not in _TIME_CONVENTIONS for value in conventions)
        ):
            raise ValueError(
                f"time_conventions must contain values drawn from {_TIME_CONVENTIONS}."
            )
        benchmarks = tuple(
            str(value).strip() for value in self.verified_benchmarks
        )
        if len(set(benchmarks)) != len(benchmarks) or any(
            not value for value in benchmarks
        ):
            raise ValueError(
                "verified_benchmarks must contain unique non-empty names."
            )
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "version", version)
        object.__setattr__(self, "dimensions", dimensions)
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "time_conventions", conventions)
        object.__setattr__(
            self,
            "maximum_cells",
            _positive_limit(self.maximum_cells, "maximum_cells"),
        )
        object.__setattr__(
            self,
            "maximum_frequencies",
            _positive_limit(self.maximum_frequencies, "maximum_frequencies"),
        )
        object.__setattr__(self, "verified_benchmarks", benchmarks)

    def supports_dimension(self, dimension: int) -> bool:
        """Return whether a spatial dimension is supported.

        Parameters
        ----------
        dimension : int
            Requested dimension.

        Returns
        -------
        bool
            Capability status.

        Examples
        --------
        >>> BackendCapabilities("b", "1", (2,), ("zxy",)).supports_dimension(3)
        False
        """
        return dimension in self.dimensions

    def supports_component(self, component: str) -> bool:
        """Return whether an impedance component is supported.

        Parameters
        ----------
        component : str
            Canonical tensor component name.

        Returns
        -------
        bool
            Capability status.

        Examples
        --------
        >>> BackendCapabilities("b", "1", (2,), ("zxy",)).supports_component(
        ...     "ZYX"
        ... )
        False
        """
        return str(component).strip().lower() in self.components

    def assess(self, problem: MaxwellProblem) -> CompatibilityReport:
        """Assess a problem without invoking the numerical backend.

        Parameters
        ----------
        problem : MaxwellProblem
            Validated problem contract.

        Returns
        -------
        CompatibilityReport
            All hard errors plus advisory validation warnings.

        Examples
        --------
        >>> from .contracts import MaxwellMesh, ReceiverSet
        >>> mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
        >>> problem = MaxwellProblem(
        ...     mesh, np.ones((2, 2)), [1], ReceiverSet([[0, 0]], ["S"])
        ... )
        >>> BackendCapabilities("b", "1", (2,), ("zxy", "zyx")).assess(
        ...     problem
        ... ).compatible
        True
        """
        if not isinstance(problem, MaxwellProblem):
            raise TypeError("problem must be a MaxwellProblem.")
        errors: list[str] = []
        warnings: list[str] = []
        if problem.mesh.dimension not in self.dimensions:
            errors.append(
                f"{problem.mesh.dimension}-D problems are unsupported"
            )
        missing = [
            value
            for value in problem.components
            if value not in self.components
        ]
        if missing:
            errors.append(f"unsupported impedance components: {missing}")
        if problem.time_dependence not in self.time_conventions:
            errors.append(
                f"time convention {problem.time_dependence!r} is unsupported"
            )
        widths = problem.mesh.cell_widths_m
        nonuniform = any(
            not np.allclose(value, value[0]) for value in widths.values()
        )
        if nonuniform and not self.supports_nonuniform_mesh:
            errors.append("nonuniform meshes are unsupported")
        inactive = ~problem.active_cells
        if np.any(inactive) and not self.supports_inactive_cells:
            errors.append("inactive cells are unsupported")
        if (
            np.any(inactive)
            and self.supports_inactive_cells
            and not self.supports_topography
        ):
            horizontal_axes = tuple(range(1, inactive.ndim))
            inactive_per_depth = np.all(inactive, axis=horizontal_axes)
            reconstructed = inactive_per_depth.reshape(
                (-1,) + (1,) * len(horizontal_axes)
            )
            if not np.array_equal(
                inactive, np.broadcast_to(reconstructed, inactive.shape)
            ):
                errors.append(
                    "laterally varying inactive cells require topography support"
                )
        cell_count = int(np.prod(problem.mesh.shape))
        if self.maximum_cells is not None and cell_count > self.maximum_cells:
            errors.append(
                f"cell count {cell_count} exceeds limit {self.maximum_cells}"
            )
        if (
            self.maximum_frequencies is not None
            and len(problem.frequencies_hz) > self.maximum_frequencies
        ):
            errors.append(
                f"frequency count {len(problem.frequencies_hz)} exceeds limit {self.maximum_frequencies}"
            )
        if not self.verified_benchmarks:
            warnings.append("backend declares no verified benchmarks")
        return CompatibilityReport(
            self.name, not errors, tuple(errors), tuple(warnings)
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible capability declaration.

        Returns
        -------
        dict
            Versioned capability state.

        Examples
        --------
        >>> BackendCapabilities("b", "1", (2,), ("zxy",)).to_dict()[
        ...     "dimensions"
        ... ]
        [2]
        """
        return {
            "schema_version": 1,
            "name": self.name,
            "version": self.version,
            "dimensions": list(self.dimensions),
            "components": list(self.components),
            "time_conventions": list(self.time_conventions),
            "supports_nonuniform_mesh": self.supports_nonuniform_mesh,
            "supports_inactive_cells": self.supports_inactive_cells,
            "supports_topography": self.supports_topography,
            "supports_anisotropy": self.supports_anisotropy,
            "maximum_cells": self.maximum_cells,
            "maximum_frequencies": self.maximum_frequencies,
            "verified_benchmarks": list(self.verified_benchmarks),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> BackendCapabilities:
        """Restore a validated capability declaration.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        BackendCapabilities
            Restored declaration.

        Examples
        --------
        >>> cap = BackendCapabilities("b", "1", (2,), ("zxy",))
        >>> BackendCapabilities.from_dict(cap.to_dict()).name
        'b'
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported BackendCapabilities schema version.")
        values = dict(data)
        values.pop("schema_version")
        for key in (
            "dimensions",
            "components",
            "time_conventions",
            "verified_benchmarks",
        ):
            values[key] = tuple(values.get(key, ()))
        return cls(**values)


@runtime_checkable
class MaxwellBackend(Protocol):
    """Runtime-checkable interface implemented by Maxwell adapters.

    Examples
    --------
    A conforming adapter exposes immutable capabilities and a solve method.

    >>> class Demo:
    ...     capabilities = BackendCapabilities(
    ...         "demo", "1", (2,), ("zxy", "zyx")
    ...     )
    ...
    ...     def solve(self, problem):
    ...         raise NotImplementedError
    >>> isinstance(Demo(), MaxwellBackend)
    True
    """

    @property
    def capabilities(self) -> BackendCapabilities:
        """Return the adapter's immutable capability declaration."""
        ...

    def solve(self, problem: MaxwellProblem) -> ForwardResult:
        """Solve a compatible problem and return canonical output."""
        ...


BackendFactory = Callable[..., MaxwellBackend]
AvailabilityProbe = Callable[[], tuple[bool, Optional[str]]]


@dataclass(frozen=True)
class BackendRegistration:
    """Store one lazy backend factory and its availability probe.

    Parameters
    ----------
    capabilities : BackendCapabilities
        Static capability declaration, available without importing the solver.
    factory : callable
        Factory returning a :class:`MaxwellBackend`.
    availability_probe : callable or None, optional
        Lightweight function returning ``(available, reason)``.

    Examples
    --------
    >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
    >>> registration = BackendRegistration(cap, lambda: object())
    >>> registration.availability()
    (True, None)
    """

    capabilities: BackendCapabilities
    factory: BackendFactory
    availability_probe: AvailabilityProbe | None = None

    def __post_init__(self) -> None:
        if not isinstance(
            self.capabilities, BackendCapabilities
        ) or not callable(self.factory):
            raise TypeError(
                "capabilities and a callable factory are required."
            )
        if self.availability_probe is not None and not callable(
            self.availability_probe
        ):
            raise TypeError("availability_probe must be callable or None.")

    def availability(self) -> tuple[bool, str | None]:
        """Return whether the optional backend can currently be created.

        Returns
        -------
        available : bool
            Probe status.
        reason : str or None
            Human-readable reason when unavailable.

        Examples
        --------
        >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
        >>> BackendRegistration(
        ...     cap, lambda: None, lambda: (False, "missing")
        ... ).availability()
        (False, 'missing')
        """
        if self.availability_probe is None:
            return True, None
        try:
            available, reason = self.availability_probe()
        except Exception as exc:
            return False, f"availability probe failed: {exc}"
        normalized_reason = (
            None if reason is None else str(reason).strip() or None
        )
        if not available and normalized_reason is None:
            normalized_reason = "backend reported itself unavailable"
        return bool(available), normalized_reason

    def create(self, **options: Any) -> MaxwellBackend:
        """Create and validate a backend instance.

        Parameters
        ----------
        **options
            Backend-specific constructor options.

        Returns
        -------
        MaxwellBackend
            Conforming adapter instance.

        Raises
        ------
        RuntimeError
            If unavailable or the factory violates its registration.

        Examples
        --------
        >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
        >>> class Demo:
        ...     capabilities = cap
        ...
        ...     def solve(self, problem):
        ...         raise NotImplementedError
        >>> BackendRegistration(cap, Demo).create().capabilities.name
        'demo'
        """
        available, reason = self.availability()
        if not available:
            raise RuntimeError(
                f"backend {self.capabilities.name!r} is unavailable: {reason}"
            )
        instance = self.factory(**options)
        if not isinstance(instance, MaxwellBackend):
            raise RuntimeError(
                "backend factory returned an object that violates MaxwellBackend."
            )
        if instance.capabilities != self.capabilities:
            raise RuntimeError(
                "backend instance capabilities differ from its registration."
            )
        return instance


class BackendRegistry:
    """Thread-safe registry of lazy Maxwell backend factories.

    Examples
    --------
    >>> registry = BackendRegistry()
    >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
    >>> registry.register(BackendRegistration(cap, lambda: None))
    >>> registry.names()
    ('demo',)
    """

    def __init__(self) -> None:
        self._registrations: dict[str, BackendRegistration] = {}
        self._lock = RLock()

    def register(
        self, registration: BackendRegistration, *, replace: bool = False
    ) -> None:
        """Register a backend under its declared capability name.

        Parameters
        ----------
        registration : BackendRegistration
            Lazy backend definition.
        replace : bool, default=False
            Explicitly replace an existing registration.

        Raises
        ------
        ValueError
            If the name already exists and replacement was not requested.

        Examples
        --------
        >>> registry = BackendRegistry()
        >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
        >>> registry.register(BackendRegistration(cap, lambda: None))
        """
        if not isinstance(registration, BackendRegistration):
            raise TypeError("registration must be a BackendRegistration.")
        name = registration.capabilities.name
        with self._lock:
            if name in self._registrations and not replace:
                raise ValueError(f"backend {name!r} is already registered.")
            self._registrations[name] = registration

    def unregister(self, name: str) -> BackendRegistration:
        """Remove and return a registration.

        Parameters
        ----------
        name : str
            Backend identifier.

        Returns
        -------
        BackendRegistration
            Removed registration.

        Examples
        --------
        >>> registry = BackendRegistry()
        >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
        >>> registry.register(BackendRegistration(cap, lambda: None))
        >>> registry.unregister("demo").capabilities.name
        'demo'
        """
        key = _identifier(name, "name")
        with self._lock:
            try:
                return self._registrations.pop(key)
            except KeyError as exc:
                raise KeyError(f"unknown Maxwell backend {key!r}.") from exc

    def get(self, name: str) -> BackendRegistration:
        """Return one registration without creating its backend.

        Parameters
        ----------
        name : str
            Backend identifier.

        Returns
        -------
        BackendRegistration
            Lazy registration.

        Examples
        --------
        >>> registry = BackendRegistry()
        >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
        >>> registry.register(BackendRegistration(cap, lambda: None))
        >>> registry.get("DEMO").capabilities.version
        '1'
        """
        key = _identifier(name, "name")
        with self._lock:
            try:
                return self._registrations[key]
            except KeyError as exc:
                raise KeyError(f"unknown Maxwell backend {key!r}.") from exc

    def names(self, *, available_only: bool = False) -> tuple[str, ...]:
        """Return sorted registered backend names.

        Parameters
        ----------
        available_only : bool, default=False
            Exclude registrations whose availability probe fails.

        Returns
        -------
        tuple of str
            Stable sorted names.

        Examples
        --------
        >>> BackendRegistry().names()
        ()
        """
        with self._lock:
            registrations = tuple(self._registrations.items())
        if available_only:
            registrations = tuple(
                (name, value)
                for name, value in registrations
                if value.availability()[0]
            )
        return tuple(sorted(name for name, _ in registrations))

    def describe(self) -> Mapping[str, Mapping[str, Any]]:
        """Return immutable availability and capability summaries.

        Returns
        -------
        mapping
            Backend names mapped to JSON-compatible summaries.

        Examples
        --------
        >>> BackendRegistry().describe() == {}
        True
        """
        with self._lock:
            registrations = tuple(self._registrations.items())
        descriptions = {}
        for name, registration in registrations:
            available, reason = registration.availability()
            descriptions[name] = MappingProxyType(
                {
                    "available": available,
                    "reason": reason,
                    "capabilities": registration.capabilities.to_dict(),
                }
            )
        return MappingProxyType(descriptions)

    def create(self, name: str, **options: Any) -> MaxwellBackend:
        """Create a named backend through its lazy factory.

        Parameters
        ----------
        name : str
            Registered backend identifier.
        **options
            Backend-specific constructor options.

        Returns
        -------
        MaxwellBackend
            Validated adapter.

        Examples
        --------
        See :meth:`BackendRegistration.create` for a complete adapter example.
        """
        return self.get(name).create(**options)


backend_registry = BackendRegistry()
"""Process-wide registry used by the convenience functions."""


def register_backend(
    registration: BackendRegistration, *, replace: bool = False
) -> None:
    """Register a lazy backend in the process-wide registry.

    Parameters
    ----------
    registration : BackendRegistration
        Backend definition.
    replace : bool, default=False
        Explicitly replace an existing name.

    Examples
    --------
    Prefer a private :class:`BackendRegistry` in isolated applications and
    tests; this function is intended for adapter package initialization.
    """
    backend_registry.register(registration, replace=replace)


def unregister_backend(name: str) -> BackendRegistration:
    """Remove a backend from the process-wide registry.

    Parameters
    ----------
    name : str
        Backend identifier.

    Returns
    -------
    BackendRegistration
        Removed definition.

    Examples
    --------
    This operation is primarily useful for plugin teardown and tests.
    """
    return backend_registry.unregister(name)


def create_backend(name: str, **options: Any) -> MaxwellBackend:
    """Create a registered Maxwell adapter lazily.

    Parameters
    ----------
    name : str
        Backend identifier.
    **options
        Backend-specific constructor options.

    Returns
    -------
    MaxwellBackend
        Validated backend instance.

    Examples
    --------
    Backend packages register their factories before this function is called.
    """
    return backend_registry.create(name, **options)


def list_backends(
    *, available_only: bool = False
) -> Mapping[str, Mapping[str, Any]]:
    """Describe registered Maxwell backends without creating them.

    Parameters
    ----------
    available_only : bool, default=False
        Exclude unavailable optional backends.

    Returns
    -------
    mapping
        Immutable capability and availability summaries.

    Examples
    --------
    >>> isinstance(list_backends(), Mapping)
    True
    """
    descriptions = backend_registry.describe()
    if not available_only:
        return descriptions
    return MappingProxyType(
        {
            name: value
            for name, value in descriptions.items()
            if value["available"]
        }
    )
