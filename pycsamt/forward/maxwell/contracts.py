# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Solver-neutral contracts for frequency-domain Maxwell simulations.

The canonical model order is ``(z, x)`` in 2-D and ``(z, y, x)`` in 3-D.
Canonical impedance output uses ``(station, frequency, component)``.  This
module contains no discretization, sparse solver, or optional backend import.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

__all__ = [
    "MaxwellMesh",
    "ReceiverSet",
    "MaxwellProblem",
    "SolverDiagnostics",
    "ForwardResult",
]

_COMPONENTS = ("zxx", "zxy", "zyx", "zyy")
_TIME_CONVENTIONS = ("exp(+iwt)", "exp(-iwt)")


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _axis_edges(value: Any, name: str) -> np.ndarray:
    edges = np.asarray(value, dtype=float)
    if edges.ndim != 1 or len(edges) < 3:
        raise ValueError(f"{name} must contain at least three cell edges.")
    if not np.all(np.isfinite(edges)) or not np.all(np.diff(edges) > 0):
        raise ValueError(f"{name} must be finite and strictly increasing.")
    return _readonly(edges)


def _json_mapping(value: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    try:
        encoded = json.dumps(dict(value), sort_keys=True, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{name} must contain finite JSON-compatible values."
        ) from exc
    return MappingProxyType(json.loads(encoded))


def _names(values: Sequence[str], expected: int, name: str) -> tuple[str, ...]:
    result = tuple(str(value).strip() for value in values)
    if len(result) != expected or any(not value for value in result):
        raise ValueError(f"{name} must contain {expected} non-empty entries.")
    if len(set(result)) != len(result):
        raise ValueError(f"{name} must be unique.")
    return result


@dataclass(frozen=True)
class MaxwellMesh:
    """Describe a rectilinear finite-volume or finite-element mesh.

    Parameters
    ----------
    x_edges_m, z_edges_m : array-like
        Strictly increasing cell-edge coordinates in metres. Depth ``z``
        increases downward.
    y_edges_m : array-like or None, optional
        Second horizontal axis. Omit for a 2-D mesh.
    crs : str or None, optional
        Coordinate reference system for horizontal coordinates.

    Examples
    --------
    >>> mesh = MaxwellMesh([0, 100, 250], [0, 50, 150])
    >>> mesh.shape, mesh.dimension
    ((2, 2), 2)
    >>> mesh.cell_centres_m["x"].tolist()
    [50.0, 175.0]
    """

    x_edges_m: np.ndarray
    z_edges_m: np.ndarray
    y_edges_m: np.ndarray | None = None
    crs: str | None = None

    def __post_init__(self) -> None:
        x = _axis_edges(self.x_edges_m, "x_edges_m")
        z = _axis_edges(self.z_edges_m, "z_edges_m")
        y = (
            None
            if self.y_edges_m is None
            else _axis_edges(self.y_edges_m, "y_edges_m")
        )
        crs = None if self.crs is None else str(self.crs).strip()
        if self.crs is not None and not crs:
            raise ValueError("crs cannot be empty.")
        object.__setattr__(self, "x_edges_m", x)
        object.__setattr__(self, "z_edges_m", z)
        object.__setattr__(self, "y_edges_m", y)
        object.__setattr__(self, "crs", crs)

    @property
    def dimension(self) -> int:
        """Return the spatial dimension.

        Returns
        -------
        {2, 3}
            Mesh dimension.

        Examples
        --------
        >>> MaxwellMesh([0, 1, 2], [0, 1, 2]).dimension
        2
        """
        return 2 if self.y_edges_m is None else 3

    @property
    def shape(self) -> tuple[int, ...]:
        """Return canonical cell-array shape.

        Returns
        -------
        tuple of int
            ``(nz, nx)`` or ``(nz, ny, nx)``.

        Examples
        --------
        >>> MaxwellMesh([0, 1, 2], [0, 1, 2, 3]).shape
        (3, 2)
        """
        nz, nx = len(self.z_edges_m) - 1, len(self.x_edges_m) - 1
        return (
            (nz, nx)
            if self.y_edges_m is None
            else (nz, len(self.y_edges_m) - 1, nx)
        )

    @property
    def cell_widths_m(self) -> Mapping[str, np.ndarray]:
        """Return read-only cell widths keyed by axis.

        Returns
        -------
        mapping
            Keys are ``x`` and ``z``, plus ``y`` for 3-D.

        Examples
        --------
        >>> mesh = MaxwellMesh([0, 2, 5], [0, 1, 3])
        >>> mesh.cell_widths_m["x"].tolist()
        [2.0, 3.0]
        """
        values = {
            "x": _readonly(np.diff(self.x_edges_m)),
            "z": _readonly(np.diff(self.z_edges_m)),
        }
        if self.y_edges_m is not None:
            values["y"] = _readonly(np.diff(self.y_edges_m))
        return MappingProxyType(values)

    @property
    def cell_centres_m(self) -> Mapping[str, np.ndarray]:
        """Return read-only cell centres keyed by axis.

        Returns
        -------
        mapping
            Centre coordinates in metres.

        Examples
        --------
        >>> mesh = MaxwellMesh([0, 2, 4], [0, 10, 20])
        >>> mesh.cell_centres_m["z"].tolist()
        [5.0, 15.0]
        """
        values = {
            "x": _readonly((self.x_edges_m[:-1] + self.x_edges_m[1:]) / 2),
            "z": _readonly((self.z_edges_m[:-1] + self.z_edges_m[1:]) / 2),
        }
        if self.y_edges_m is not None:
            values["y"] = _readonly(
                (self.y_edges_m[:-1] + self.y_edges_m[1:]) / 2
            )
        return MappingProxyType(values)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible mesh representation.

        Returns
        -------
        dict
            Versioned mesh state.

        Examples
        --------
        >>> MaxwellMesh([0, 1, 2], [0, 1, 2]).to_dict()["schema_version"]
        1
        """
        return {
            "schema_version": 1,
            "x_edges_m": self.x_edges_m.tolist(),
            "z_edges_m": self.z_edges_m.tolist(),
            "y_edges_m": None
            if self.y_edges_m is None
            else self.y_edges_m.tolist(),
            "crs": self.crs,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> MaxwellMesh:
        """Restore a validated mesh from serialized state.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        MaxwellMesh
            Restored immutable mesh.

        Examples
        --------
        >>> mesh = MaxwellMesh([0, 1, 2], [0, 2, 4])
        >>> MaxwellMesh.from_dict(mesh.to_dict()).shape
        (2, 2)
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported MaxwellMesh schema version.")
        return cls(
            data["x_edges_m"],
            data["z_edges_m"],
            data.get("y_edges_m"),
            data.get("crs"),
        )


@dataclass(frozen=True)
class ReceiverSet:
    """Define named receiver locations in mesh coordinates.

    Parameters
    ----------
    coordinates_m : array-like, shape (n, dimension)
        ``x,z`` locations for 2-D or ``x,y,z`` locations for 3-D.
    names : sequence of str
        Unique receiver or station identifiers.
    orientation_deg : float, default=0.0
        Clockwise rotation of receiver x/y axes from the model axes.

    Examples
    --------
    >>> receivers = ReceiverSet([[50, 0], [150, 0]], ["S00", "S01"])
    >>> receivers.dimension, receivers.count
    (2, 2)
    """

    coordinates_m: np.ndarray
    names: tuple[str, ...]
    orientation_deg: float = 0.0

    def __post_init__(self) -> None:
        coordinates = np.asarray(self.coordinates_m, dtype=float)
        if (
            coordinates.ndim != 2
            or coordinates.shape[0] < 1
            or coordinates.shape[1] not in (2, 3)
        ):
            raise ValueError("coordinates_m must have shape (n, 2) or (n, 3).")
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("coordinates_m must be finite.")
        names = _names(self.names, len(coordinates), "names")
        orientation = float(self.orientation_deg)
        if not np.isfinite(orientation):
            raise ValueError("orientation_deg must be finite.")
        object.__setattr__(self, "coordinates_m", _readonly(coordinates))
        object.__setattr__(self, "names", names)
        object.__setattr__(self, "orientation_deg", orientation % 360.0)

    @property
    def count(self) -> int:
        """Return the number of receivers.

        Returns
        -------
        int
            Receiver count.

        Examples
        --------
        >>> ReceiverSet([[0, 0]], ["S00"]).count
        1
        """
        return len(self.names)

    @property
    def dimension(self) -> int:
        """Return the coordinate dimension.

        Returns
        -------
        {2, 3}
            Number of coordinate columns.

        Examples
        --------
        >>> ReceiverSet([[0, 0, 0]], ["S00"]).dimension
        3
        """
        return self.coordinates_m.shape[1]

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible receiver state.

        Returns
        -------
        dict
            Versioned receiver coordinates and names.

        Examples
        --------
        >>> ReceiverSet([[0, 0]], ["S00"]).to_dict()["names"]
        ['S00']
        """
        return {
            "schema_version": 1,
            "coordinates_m": self.coordinates_m.tolist(),
            "names": list(self.names),
            "orientation_deg": self.orientation_deg,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ReceiverSet:
        """Restore receivers from serialized state.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        ReceiverSet
            Validated receiver collection.

        Examples
        --------
        >>> state = ReceiverSet([[0, 0]], ["S00"]).to_dict()
        >>> ReceiverSet.from_dict(state).names
        ('S00',)
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported ReceiverSet schema version.")
        return cls(
            data["coordinates_m"],
            tuple(data["names"]),
            data.get("orientation_deg", 0.0),
        )


@dataclass(frozen=True)
class MaxwellProblem:
    """Define one isotropic frequency-domain MT boundary-value problem.

    Parameters
    ----------
    mesh : MaxwellMesh
        Rectilinear simulation mesh.
    conductivity_s_m : array-like
        Positive isotropic conductivity in S/m, shaped like ``mesh``.
    frequencies_hz : array-like, shape (n_frequency,)
        Positive, unique frequencies. Input order is retained.
    receivers : ReceiverSet
        Observation locations with the same dimension as ``mesh``.
    components : sequence of {"zxx", "zxy", "zyx", "zyy"}
        Requested impedance components in explicit output order.
    active_cells : array-like of bool or None, optional
        Cells participating in the physical model. Air cells can remain active
        with a small conductivity or be marked inactive for capable backends.
    time_dependence : {"exp(+iwt)", "exp(-iwt)"}, default="exp(+iwt)"
        Complex phasor convention.
    magnetic_permeability_h_m : float, default=1.25663706212e-6
        Uniform scalar permeability in H/m.
    metadata : mapping, optional
        Finite JSON-compatible provenance; never interpreted by a backend.

    Examples
    --------
    >>> mesh = MaxwellMesh([0, 100, 200], [0, 50, 100])
    >>> receivers = ReceiverSet([[50, 0]], ["S00"])
    >>> problem = MaxwellProblem(
    ...     mesh, np.full(mesh.shape, 0.01), [10, 1], receivers
    ... )
    >>> problem.problem_hash == problem.problem_hash
    True
    """

    mesh: MaxwellMesh
    conductivity_s_m: np.ndarray
    frequencies_hz: np.ndarray
    receivers: ReceiverSet
    components: tuple[str, ...] = ("zxy", "zyx")
    active_cells: np.ndarray | None = None
    time_dependence: str = "exp(+iwt)"
    magnetic_permeability_h_m: float = 4.0e-7 * np.pi
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.mesh, MaxwellMesh) or not isinstance(
            self.receivers, ReceiverSet
        ):
            raise TypeError(
                "mesh and receivers must use Maxwell contract types."
            )
        if self.mesh.dimension != self.receivers.dimension:
            raise ValueError("mesh and receiver dimensions must match.")
        conductivity = np.asarray(self.conductivity_s_m, dtype=float)
        if (
            conductivity.shape != self.mesh.shape
            or not np.all(np.isfinite(conductivity))
            or np.any(conductivity <= 0)
        ):
            raise ValueError(
                f"conductivity_s_m must be positive, finite, and shaped {self.mesh.shape}."
            )
        frequencies = np.asarray(self.frequencies_hz, dtype=float)
        if (
            frequencies.ndim != 1
            or len(frequencies) < 1
            or not np.all(np.isfinite(frequencies))
            or np.any(frequencies <= 0)
        ):
            raise ValueError(
                "frequencies_hz must be a non-empty vector of positive finite values."
            )
        if len(np.unique(frequencies)) != len(frequencies):
            raise ValueError("frequencies_hz must be unique.")
        components = tuple(
            str(value).strip().lower() for value in self.components
        )
        if (
            not components
            or len(set(components)) != len(components)
            or any(value not in _COMPONENTS for value in components)
        ):
            raise ValueError(
                f"components must be unique values drawn from {_COMPONENTS}."
            )
        if self.mesh.dimension == 2 and any(
            value in {"zxx", "zyy"} for value in components
        ):
            raise ValueError(
                "2-D problems support only zxy and zyx impedance components."
            )
        active = (
            np.ones(self.mesh.shape, dtype=bool)
            if self.active_cells is None
            else np.asarray(self.active_cells, dtype=bool)
        )
        if active.shape != self.mesh.shape or not np.any(active):
            raise ValueError(
                f"active_cells must be shaped {self.mesh.shape} and contain an active cell."
            )
        if self.time_dependence not in _TIME_CONVENTIONS:
            raise ValueError(
                f"time_dependence must be one of {_TIME_CONVENTIONS}."
            )
        permeability = float(self.magnetic_permeability_h_m)
        if not np.isfinite(permeability) or permeability <= 0:
            raise ValueError(
                "magnetic_permeability_h_m must be positive and finite."
            )
        object.__setattr__(self, "conductivity_s_m", _readonly(conductivity))
        object.__setattr__(self, "frequencies_hz", _readonly(frequencies))
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "active_cells", _readonly(active, bool))
        object.__setattr__(self, "magnetic_permeability_h_m", permeability)
        object.__setattr__(
            self, "metadata", _json_mapping(self.metadata, "metadata")
        )

    @property
    def problem_hash(self) -> str:
        """Return a deterministic SHA-256 digest of all physical inputs.

        Returns
        -------
        str
            Digest suitable for cache keys.

        Examples
        --------
        >>> mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
        >>> p = MaxwellProblem(
        ...     mesh, np.ones(mesh.shape), [1], ReceiverSet([[0, 0]], ["S"])
        ... )
        >>> len(p.problem_hash)
        64
        """
        digest = hashlib.sha256()
        digest.update(
            np.ascontiguousarray(self.conductivity_s_m, dtype="<f8").tobytes()
        )
        digest.update(
            np.ascontiguousarray(self.frequencies_hz, dtype="<f8").tobytes()
        )
        digest.update(
            np.ascontiguousarray(self.active_cells, dtype=np.uint8).tobytes()
        )
        digest.update(
            json.dumps(
                self.provenance(), sort_keys=True, separators=(",", ":")
            ).encode("utf-8")
        )
        return digest.hexdigest()

    def provenance(self) -> dict[str, Any]:
        """Return JSON-compatible problem provenance excluding large arrays.

        Returns
        -------
        dict
            Mesh, receiver, convention, components, and metadata state.

        Examples
        --------
        >>> mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
        >>> p = MaxwellProblem(
        ...     mesh, np.ones(mesh.shape), [1], ReceiverSet([[0, 0]], ["S"])
        ... )
        >>> p.provenance()["components"]
        ['zxy', 'zyx']
        """
        return {
            "schema_version": 1,
            "mesh": self.mesh.to_dict(),
            "receivers": self.receivers.to_dict(),
            "components": list(self.components),
            "time_dependence": self.time_dependence,
            "magnetic_permeability_h_m": self.magnetic_permeability_h_m,
            "metadata": dict(self.metadata),
        }

    def to_npz(self, path: str | Path) -> Path:
        """Write a pickle-free problem archive.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination archive.

        Returns
        -------
        pathlib.Path
            Destination path.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> p = MaxwellProblem(
        ...     MaxwellMesh([0, 1, 2], [0, 1, 2]),
        ...     np.ones((2, 2)),
        ...     [1],
        ...     ReceiverSet([[0, 0]], ["S"]),
        ... )
        >>> with TemporaryDirectory() as d:
        ...     restored = MaxwellProblem.from_npz(p.to_npz(Path(d) / "p.npz"))
        >>> restored.problem_hash == p.problem_hash
        True
        """
        target = Path(path)
        np.savez_compressed(
            target,
            conductivity_s_m=self.conductivity_s_m,
            frequencies_hz=self.frequencies_hz,
            active_cells=self.active_cells,
            provenance_json=np.array(
                json.dumps(self.provenance(), sort_keys=True)
            ),
        )
        return target

    @classmethod
    def from_npz(cls, path: str | Path) -> MaxwellProblem:
        """Restore and validate a problem archive.

        Parameters
        ----------
        path : str or pathlib.Path
            Archive written by :meth:`to_npz`.

        Returns
        -------
        MaxwellProblem
            Restored problem.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> p = MaxwellProblem(
        ...     MaxwellMesh([0, 1, 2], [0, 1, 2]),
        ...     np.ones((2, 2)),
        ...     [1],
        ...     ReceiverSet([[0, 0]], ["S"]),
        ... )
        >>> with TemporaryDirectory() as d:
        ...     q = MaxwellProblem.from_npz(p.to_npz(Path(d) / "p.npz"))
        >>> np.array_equal(q.conductivity_s_m, p.conductivity_s_m)
        True
        """
        with np.load(Path(path), allow_pickle=False) as archive:
            state = json.loads(str(archive["provenance_json"].item()))
            if state.get("schema_version") != 1:
                raise ValueError("unsupported MaxwellProblem schema version.")
            return cls(
                MaxwellMesh.from_dict(state["mesh"]),
                archive["conductivity_s_m"],
                archive["frequencies_hz"],
                ReceiverSet.from_dict(state["receivers"]),
                tuple(state["components"]),
                archive["active_cells"],
                state["time_dependence"],
                state["magnetic_permeability_h_m"],
                state.get("metadata", {}),
            )


@dataclass(frozen=True)
class SolverDiagnostics:
    """Record convergence information for every frequency and source solve.

    Parameters
    ----------
    converged : array-like of bool, shape (n_frequency, n_source)
        Whether each linear solve met its tolerance.
    iterations : array-like of int, same shape
        Iteration count; zero is valid for direct solvers.
    relative_residual : array-like of float, same shape
        Final non-negative relative residual.
    runtime_s : float
        Total non-negative solver runtime in seconds.
    messages : sequence of str, optional
        Backend messages for failed or noteworthy solves.

    Examples
    --------
    >>> d = SolverDiagnostics(
    ...     [[True], [False]], [[4], [20]], [[1e-8], [1e-2]], 0.5
    ... )
    >>> d.success, d.maximum_relative_residual
    (False, 0.01)
    """

    converged: np.ndarray
    iterations: np.ndarray
    relative_residual: np.ndarray
    runtime_s: float
    messages: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        converged = np.asarray(self.converged, dtype=bool)
        iterations = np.asarray(self.iterations)
        residual = np.asarray(self.relative_residual, dtype=float)
        if converged.ndim != 2 or converged.size < 1:
            raise ValueError("converged must be a non-empty 2-D array.")
        if (
            iterations.shape != converged.shape
            or residual.shape != converged.shape
        ):
            raise ValueError("diagnostic arrays must have identical shapes.")
        if not np.issubdtype(iterations.dtype, np.integer) or np.any(
            iterations < 0
        ):
            raise ValueError("iterations must contain non-negative integers.")
        if not np.all(np.isfinite(residual)) or np.any(residual < 0):
            raise ValueError(
                "relative_residual must be finite and non-negative."
            )
        runtime = float(self.runtime_s)
        if not np.isfinite(runtime) or runtime < 0:
            raise ValueError("runtime_s must be finite and non-negative.")
        messages = tuple(str(value).strip() for value in self.messages)
        object.__setattr__(self, "converged", _readonly(converged, bool))
        object.__setattr__(self, "iterations", _readonly(iterations, np.int64))
        object.__setattr__(self, "relative_residual", _readonly(residual))
        object.__setattr__(self, "runtime_s", runtime)
        object.__setattr__(self, "messages", messages)

    @property
    def success(self) -> bool:
        """Return whether every solve converged.

        Returns
        -------
        bool
            True only when all convergence flags are true.

        Examples
        --------
        >>> SolverDiagnostics([[True]], [[1]], [[1e-9]], 0).success
        True
        """
        return bool(np.all(self.converged))

    @property
    def maximum_relative_residual(self) -> float:
        """Return the largest reported relative residual.

        Returns
        -------
        float
            Worst solve residual.

        Examples
        --------
        >>> SolverDiagnostics(
        ...     [[True]], [[1]], [[1e-7]], 0
        ... ).maximum_relative_residual
        1e-07
        """
        return float(np.max(self.relative_residual))

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible convergence diagnostics.

        Returns
        -------
        dict
            Versioned diagnostic state.

        Examples
        --------
        >>> SolverDiagnostics([[True]], [[2]], [[1e-8]], 0.1).to_dict()[
        ...     "runtime_s"
        ... ]
        0.1
        """
        return {
            "schema_version": 1,
            "converged": self.converged.tolist(),
            "iterations": self.iterations.tolist(),
            "relative_residual": self.relative_residual.tolist(),
            "runtime_s": self.runtime_s,
            "messages": list(self.messages),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> SolverDiagnostics:
        """Restore validated convergence diagnostics.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        SolverDiagnostics
            Restored diagnostic record.

        Examples
        --------
        >>> state = SolverDiagnostics([[True]], [[2]], [[1e-8]], 0.1).to_dict()
        >>> SolverDiagnostics.from_dict(state).success
        True
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported SolverDiagnostics schema version.")
        return cls(
            data["converged"],
            data["iterations"],
            data["relative_residual"],
            data["runtime_s"],
            tuple(data.get("messages", ())),
        )


@dataclass(frozen=True)
class ForwardResult:
    """Store canonical impedance predictions from a Maxwell backend.

    Parameters
    ----------
    problem_hash : str
        Hash of the exact :class:`MaxwellProblem` solved.
    frequencies_hz : array-like
        Frequency vector in problem order.
    receiver_names, components : sequence of str
        Explicit station and tensor-component axes.
    impedance_v_a : complex array, shape (station, frequency, component)
        Predicted SI impedance.
    valid : bool array or None, optional
        Validity mask with the same shape. Defaults to finite predictions.
    backend_name, backend_version : str
        Solver identity required for reproducibility.
    diagnostics : SolverDiagnostics
        Per-solve convergence record.
    metadata : mapping, optional
        Additional finite JSON-compatible backend provenance.

    Examples
    --------
    >>> d = SolverDiagnostics([[True]], [[3]], [[1e-9]], 0.01)
    >>> r = ForwardResult(
    ...     "a" * 64, [1], ["S"], ["zxy"], [[[1 + 2j]]], None, "demo", "1", d
    ... )
    >>> r.shape, r.success
    ((1, 1, 1), True)
    """

    problem_hash: str
    frequencies_hz: np.ndarray
    receiver_names: tuple[str, ...]
    components: tuple[str, ...]
    impedance_v_a: np.ndarray
    valid: np.ndarray | None
    backend_name: str
    backend_version: str
    diagnostics: SolverDiagnostics
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        problem_hash = str(self.problem_hash).lower()
        if len(problem_hash) != 64 or any(
            value not in "0123456789abcdef" for value in problem_hash
        ):
            raise ValueError(
                "problem_hash must be a 64-character hexadecimal SHA-256 digest."
            )
        frequencies = np.asarray(self.frequencies_hz, dtype=float)
        if (
            frequencies.ndim != 1
            or len(frequencies) < 1
            or not np.all(np.isfinite(frequencies))
            or np.any(frequencies <= 0)
            or len(np.unique(frequencies)) != len(frequencies)
        ):
            raise ValueError(
                "frequencies_hz must contain unique positive finite values."
            )
        impedance = np.asarray(self.impedance_v_a, dtype=complex)
        names = _names(
            self.receiver_names,
            impedance.shape[0] if impedance.ndim == 3 else 0,
            "receiver_names",
        )
        components = tuple(
            str(value).strip().lower() for value in self.components
        )
        expected = (len(names), len(frequencies), len(components))
        if (
            impedance.shape != expected
            or not components
            or len(set(components)) != len(components)
            or any(value not in _COMPONENTS for value in components)
        ):
            raise ValueError(
                f"impedance_v_a and axes must define shape {expected} with canonical components."
            )
        valid = (
            np.isfinite(impedance.real) & np.isfinite(impedance.imag)
            if self.valid is None
            else np.asarray(self.valid, dtype=bool)
        )
        if valid.shape != impedance.shape:
            raise ValueError("valid must have the impedance shape.")
        if np.any(
            valid
            & ~(np.isfinite(impedance.real) & np.isfinite(impedance.imag))
        ):
            raise ValueError("valid impedance entries must be finite.")
        backend_name, backend_version = (
            str(self.backend_name).strip(),
            str(self.backend_version).strip(),
        )
        if (
            not backend_name
            or not backend_version
            or not isinstance(self.diagnostics, SolverDiagnostics)
        ):
            raise ValueError(
                "backend identity and SolverDiagnostics are required."
            )
        if self.diagnostics.converged.shape[0] != len(frequencies):
            raise ValueError(
                "diagnostics first axis must match frequencies_hz."
            )
        object.__setattr__(self, "problem_hash", problem_hash)
        object.__setattr__(self, "frequencies_hz", _readonly(frequencies))
        object.__setattr__(self, "receiver_names", names)
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "impedance_v_a", _readonly(impedance))
        object.__setattr__(self, "valid", _readonly(valid, bool))
        object.__setattr__(self, "backend_name", backend_name)
        object.__setattr__(self, "backend_version", backend_version)
        object.__setattr__(
            self, "metadata", _json_mapping(self.metadata, "metadata")
        )

    @property
    def shape(self) -> tuple[int, int, int]:
        """Return canonical impedance shape.

        Returns
        -------
        tuple of int
            ``(station, frequency, component)``.

        Examples
        --------
        >>> d = SolverDiagnostics([[True]], [[0]], [[0]], 0)
        >>> ForwardResult(
        ...     "0" * 64, [1], ["S"], ["zxy"], [[[1j]]], None, "b", "1", d
        ... ).shape
        (1, 1, 1)
        """
        return self.impedance_v_a.shape

    @property
    def success(self) -> bool:
        """Return whether all solves converged and predictions are valid.

        Returns
        -------
        bool
            Combined numerical and observation validity status.

        Examples
        --------
        >>> d = SolverDiagnostics([[True]], [[0]], [[0]], 0)
        >>> ForwardResult(
        ...     "0" * 64, [1], ["S"], ["zxy"], [[[1j]]], None, "b", "1", d
        ... ).success
        True
        """
        return self.diagnostics.success and bool(np.all(self.valid))

    def validate_against(self, problem: MaxwellProblem) -> None:
        """Raise if this result does not exactly match a problem contract.

        Parameters
        ----------
        problem : MaxwellProblem
            Expected input problem.

        Raises
        ------
        ValueError
            If hash or any output axis differs.

        Examples
        --------
        >>> mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
        >>> p = MaxwellProblem(
        ...     mesh,
        ...     np.ones((2, 2)),
        ...     [1],
        ...     ReceiverSet([[0, 0]], ["S"]),
        ...     ("zxy",),
        ... )
        >>> d = SolverDiagnostics([[True]], [[0]], [[0]], 0)
        >>> ForwardResult(
        ...     p.problem_hash,
        ...     [1],
        ...     ["S"],
        ...     ["zxy"],
        ...     [[[1j]]],
        ...     None,
        ...     "b",
        ...     "1",
        ...     d,
        ... ).validate_against(p)
        """
        if self.problem_hash != problem.problem_hash:
            raise ValueError("result problem_hash does not match the problem.")
        if (
            not np.array_equal(self.frequencies_hz, problem.frequencies_hz)
            or self.receiver_names != problem.receivers.names
            or self.components != problem.components
        ):
            raise ValueError(
                "result frequency, receiver, or component axes do not match the problem."
            )

    def provenance(self) -> dict[str, Any]:
        """Return JSON-compatible solver and output-axis provenance.

        Returns
        -------
        dict
            Problem identity, axes, backend, diagnostics, and metadata.

        Examples
        --------
        >>> d = SolverDiagnostics([[True]], [[0]], [[0]], 0)
        >>> r = ForwardResult(
        ...     "0" * 64, [1], ["S"], ["zxy"], [[[1j]]], None, "b", "1", d
        ... )
        >>> r.provenance()["backend_name"]
        'b'
        """
        return {
            "schema_version": 1,
            "problem_hash": self.problem_hash,
            "receiver_names": list(self.receiver_names),
            "components": list(self.components),
            "backend_name": self.backend_name,
            "backend_version": self.backend_version,
            "diagnostics": self.diagnostics.to_dict(),
            "metadata": dict(self.metadata),
        }

    def to_npz(self, path: str | Path) -> Path:
        """Write a versioned, pickle-free result archive.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination archive.

        Returns
        -------
        pathlib.Path
            Destination path.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> d = SolverDiagnostics([[True]], [[0]], [[0]], 0)
        >>> r = ForwardResult(
        ...     "0" * 64, [1], ["S"], ["zxy"], [[[1j]]], None, "b", "1", d
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     restored = ForwardResult.from_npz(
        ...         r.to_npz(Path(directory) / "r.npz")
        ...     )
        >>> restored.backend_name, restored.shape
        ('b', (1, 1, 1))
        """
        target = Path(path)
        np.savez_compressed(
            target,
            frequencies_hz=self.frequencies_hz,
            impedance_v_a=self.impedance_v_a,
            valid=self.valid,
            provenance_json=np.array(
                json.dumps(self.provenance(), sort_keys=True)
            ),
        )
        return target

    @classmethod
    def from_npz(cls, path: str | Path) -> ForwardResult:
        """Restore and validate a result archive without pickle.

        Parameters
        ----------
        path : str or pathlib.Path
            Archive written by :meth:`to_npz`.

        Returns
        -------
        ForwardResult
            Restored canonical result.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> d = SolverDiagnostics([[True]], [[0]], [[0]], 0)
        >>> r = ForwardResult(
        ...     "0" * 64, [1], ["S"], ["zxy"], [[[1j]]], None, "b", "1", d
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     restored = ForwardResult.from_npz(
        ...         r.to_npz(Path(directory) / "r.npz")
        ...     )
        >>> np.array_equal(restored.impedance_v_a, r.impedance_v_a)
        True
        """
        with np.load(Path(path), allow_pickle=False) as archive:
            state = json.loads(str(archive["provenance_json"].item()))
            if state.get("schema_version") != 1:
                raise ValueError("unsupported ForwardResult schema version.")
            return cls(
                state["problem_hash"],
                archive["frequencies_hz"],
                tuple(state["receiver_names"]),
                tuple(state["components"]),
                archive["impedance_v_a"],
                archive["valid"],
                state["backend_name"],
                state["backend_version"],
                SolverDiagnostics.from_dict(state["diagnostics"]),
                state.get("metadata", {}),
            )
