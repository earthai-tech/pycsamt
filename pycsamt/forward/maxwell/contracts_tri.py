# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Solver-neutral contracts for unstructured triangular-mesh Maxwell physics.

Sibling of :mod:`pycsamt.forward.maxwell.contracts`, whose ``MaxwellMesh``/
``MaxwellProblem`` pair is rectilinear-only by design. :class:`TriMesh` and
:class:`TriProblem` are the same kind of solver-neutral, frozen, validated
contract, but describe an unstructured 2-D triangulation (nodes + triangle
connectivity) instead of tensor-grid edges. Per-triangle values replace the
``(nz, nx)`` cell array used by the rectilinear contracts.

:class:`~pycsamt.forward.maxwell.contracts.ReceiverSet`,
:class:`~pycsamt.forward.maxwell.contracts.SolverDiagnostics`, and
:class:`~pycsamt.forward.maxwell.contracts.ForwardResult` are reused
unchanged — none of them reference the mesh's internal shape, only
``problem_hash``/``frequencies_hz``/``receiver_names``/``components``, which
:class:`TriProblem` exposes identically to
:class:`~pycsamt.forward.maxwell.contracts.MaxwellProblem`.

This module contains no discretization, mesh-generation, or Triangle/PSLG
I/O. Building a :class:`TriMesh` from real survey geometry belongs to a
mesh-generation module (e.g. :mod:`pycsamt.models.mare2dem`); this module
only defines the validated data shape that generator hands off.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .contracts import (
    _COMPONENTS,
    _TIME_CONVENTIONS,
    ReceiverSet,
    _json_mapping,
    _readonly,
)

__all__ = ["TriMesh", "TriProblem"]


def _tri_nodes(value: Any, name: str) -> np.ndarray:
    nodes = np.asarray(value, dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] != 2 or nodes.shape[0] < 3:
        raise ValueError(
            f"{name} must have shape (n, 2) with at least three nodes."
        )
    if not np.all(np.isfinite(nodes)):
        raise ValueError(f"{name} must be finite.")
    return _readonly(nodes)


def _tri_connectivity(value: Any, name: str, n_nodes: int) -> np.ndarray:
    triangles = np.asarray(value, dtype=np.int64)
    if triangles.ndim != 2 or triangles.shape[1] != 3 or triangles.shape[0] < 1:
        raise ValueError(
            f"{name} must have shape (m, 3) with at least one triangle."
        )
    if np.any(triangles < 0) or np.any(triangles >= n_nodes):
        raise ValueError(f"{name} indices must reference nodes_m rows.")
    if (
        np.any(triangles[:, 0] == triangles[:, 1])
        or np.any(triangles[:, 1] == triangles[:, 2])
        or np.any(triangles[:, 0] == triangles[:, 2])
    ):
        raise ValueError(f"{name} triangles must reference three distinct nodes.")
    return _readonly(triangles, np.int64)


def _tri_segments(value: Any, name: str, n_nodes: int) -> np.ndarray:
    segments = np.asarray(value, dtype=np.int64)
    if segments.ndim != 2 or segments.shape[1] != 2 or segments.shape[0] < 1:
        raise ValueError(
            f"{name} must have shape (k, 2) with at least one segment."
        )
    if np.any(segments < 0) or np.any(segments >= n_nodes):
        raise ValueError(f"{name} indices must reference nodes_m rows.")
    if np.any(segments[:, 0] == segments[:, 1]):
        raise ValueError(f"{name} segments must reference two distinct nodes.")
    return _readonly(segments, np.int64)


def _signed_double_areas(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    p0 = nodes[triangles[:, 0]]
    p1 = nodes[triangles[:, 1]]
    p2 = nodes[triangles[:, 2]]
    return (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (
        p2[:, 0] - p0[:, 0]
    ) * (p1[:, 1] - p0[:, 1])


@dataclass(frozen=True)
class TriMesh:
    """Describe an unstructured 2-D triangular finite-element mesh.

    Parameters
    ----------
    nodes_m : array-like, shape (n_nodes, 2)
        Node ``(x, z)`` coordinates in metres. Depth ``z`` increases
        downward, matching :class:`~pycsamt.forward.maxwell.contracts.MaxwellMesh`.
    triangles : array-like of int, shape (n_triangles, 3)
        0-based node-index connectivity, one row per triangle.
    region_ids : array-like of int, shape (n_triangles,), optional
        Material/region label per triangle. Defaults to all zeros (a single
        region) when omitted.
    boundary_segments : array-like of int, shape (n_segments, 2), optional
        0-based node-index pairs marking PSLG boundary edges (region
        divides, survey extent). Purely descriptive; not required for a
        valid mesh.
    crs : str or None, optional
        Coordinate reference system for horizontal coordinates.

    Examples
    --------
    >>> mesh = TriMesh(
    ...     nodes_m=[[0, 0], [100, 0], [50, 50], [150, 60]],
    ...     triangles=[[0, 1, 2], [1, 3, 2]],
    ... )
    >>> mesh.dimension, mesh.n_nodes, mesh.n_triangles
    (2, 4, 2)
    >>> mesh.shape
    (2,)
    """

    nodes_m: np.ndarray
    triangles: np.ndarray
    region_ids: np.ndarray | None = None
    boundary_segments: np.ndarray | None = None
    crs: str | None = None

    def __post_init__(self) -> None:
        nodes = _tri_nodes(self.nodes_m, "nodes_m")
        triangles = _tri_connectivity(self.triangles, "triangles", len(nodes))
        areas = _signed_double_areas(nodes, triangles)
        if np.any(np.abs(areas) <= 0.0):
            raise ValueError(
                "triangles must have strictly positive area "
                "(degenerate or collinear triangle detected)."
            )
        region_ids = (
            np.zeros(len(triangles), dtype=np.int64)
            if self.region_ids is None
            else np.asarray(self.region_ids, dtype=np.int64)
        )
        if region_ids.shape != (len(triangles),):
            raise ValueError("region_ids must contain one entry per triangle.")
        boundary_segments = (
            None
            if self.boundary_segments is None
            else _tri_segments(
                self.boundary_segments, "boundary_segments", len(nodes)
            )
        )
        crs = None if self.crs is None else str(self.crs).strip()
        if self.crs is not None and not crs:
            raise ValueError("crs cannot be empty.")
        object.__setattr__(self, "nodes_m", nodes)
        object.__setattr__(self, "triangles", triangles)
        object.__setattr__(self, "region_ids", _readonly(region_ids, np.int64))
        object.__setattr__(self, "boundary_segments", boundary_segments)
        object.__setattr__(self, "crs", crs)

    @property
    def dimension(self) -> int:
        """Return the spatial dimension.

        Returns
        -------
        int
            Always 2 for this contract.

        Examples
        --------
        >>> TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]]).dimension
        2
        """
        return 2

    @property
    def n_nodes(self) -> int:
        """Return the number of mesh nodes.

        Examples
        --------
        >>> TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]]).n_nodes
        3
        """
        return int(self.nodes_m.shape[0])

    @property
    def n_triangles(self) -> int:
        """Return the number of mesh triangles.

        Examples
        --------
        >>> TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]]).n_triangles
        1
        """
        return int(self.triangles.shape[0])

    @property
    def shape(self) -> tuple[int]:
        """Return canonical per-triangle array shape.

        Returns
        -------
        tuple of int
            ``(n_triangles,)``, mirroring
            :attr:`~pycsamt.forward.maxwell.contracts.MaxwellMesh.shape`'s
            role for :class:`TriProblem` conductivity/active-cell arrays.

        Examples
        --------
        >>> TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]]).shape
        (1,)
        """
        return (self.n_triangles,)

    @property
    def triangle_centroids_m(self) -> np.ndarray:
        """Return the ``(x, z)`` centroid of every triangle.

        Returns
        -------
        ndarray, shape (n_triangles, 2)
            Mean of each triangle's three node coordinates.

        Examples
        --------
        >>> mesh = TriMesh([[0, 0], [3, 0], [0, 3]], [[0, 1, 2]])
        >>> mesh.triangle_centroids_m.tolist()
        [[1.0, 1.0]]
        """
        corners = self.nodes_m[self.triangles]
        return _readonly(corners.mean(axis=1))

    @property
    def triangle_areas_m2(self) -> np.ndarray:
        """Return the area of every triangle.

        Returns
        -------
        ndarray, shape (n_triangles,)
            Positive triangle areas in square metres.

        Examples
        --------
        >>> mesh = TriMesh([[0, 0], [4, 0], [0, 3]], [[0, 1, 2]])
        >>> mesh.triangle_areas_m2.tolist()
        [6.0]
        """
        areas = np.abs(_signed_double_areas(self.nodes_m, self.triangles)) / 2.0
        return _readonly(areas)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible mesh representation.

        Returns
        -------
        dict
            Versioned mesh state.

        Examples
        --------
        >>> TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]]).to_dict()[
        ...     "schema_version"
        ... ]
        1
        """
        return {
            "schema_version": 1,
            "nodes_m": self.nodes_m.tolist(),
            "triangles": self.triangles.tolist(),
            "region_ids": self.region_ids.tolist(),
            "boundary_segments": None
            if self.boundary_segments is None
            else self.boundary_segments.tolist(),
            "crs": self.crs,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> TriMesh:
        """Restore a validated mesh from serialized state.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        TriMesh
            Restored immutable mesh.

        Examples
        --------
        >>> mesh = TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]])
        >>> TriMesh.from_dict(mesh.to_dict()).n_triangles
        1
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported TriMesh schema version.")
        return cls(
            data["nodes_m"],
            data["triangles"],
            data.get("region_ids"),
            data.get("boundary_segments"),
            data.get("crs"),
        )


@dataclass(frozen=True)
class TriProblem:
    """Define one isotropic frequency-domain MT problem on a triangular mesh.

    Parameters
    ----------
    mesh : TriMesh
        Unstructured simulation mesh.
    conductivity_s_m : array-like, shape (n_triangles,)
        Positive isotropic conductivity in S/m, one value per triangle.
    frequencies_hz : array-like, shape (n_frequency,)
        Positive, unique frequencies. Input order is retained.
    receivers : ReceiverSet
        Observation locations with the same dimension as ``mesh``.
    components : sequence of {"zxy", "zyx"}
        Requested impedance components; only the 2-D pair is valid since
        :class:`TriMesh` is 2-D-only in this release.
    active_cells : array-like of bool or None, optional
        Triangles participating in the physical model.
    time_dependence : {"exp(+iwt)", "exp(-iwt)"}, default="exp(+iwt)"
        Complex phasor convention.
    magnetic_permeability_h_m : float, default=1.25663706212e-6
        Uniform scalar permeability in H/m.
    metadata : mapping, optional
        Finite JSON-compatible provenance; never interpreted by a backend.

    Examples
    --------
    >>> from pycsamt.forward.maxwell.contracts import ReceiverSet
    >>> mesh = TriMesh([[0, 0], [200, 0], [100, 100]], [[0, 1, 2]])
    >>> receivers = ReceiverSet([[100, 0]], ["S00"])
    >>> problem = TriProblem(mesh, [0.01], [10, 1], receivers)
    >>> problem.problem_hash == problem.problem_hash
    True
    """

    mesh: TriMesh
    conductivity_s_m: np.ndarray
    frequencies_hz: np.ndarray
    receivers: ReceiverSet
    components: tuple[str, ...] = ("zxy", "zyx")
    active_cells: np.ndarray | None = None
    time_dependence: str = "exp(+iwt)"
    magnetic_permeability_h_m: float = 4.0e-7 * np.pi
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.mesh, TriMesh) or not isinstance(
            self.receivers, ReceiverSet
        ):
            raise TypeError(
                "mesh must be a TriMesh and receivers a ReceiverSet."
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
        if any(value in {"zxx", "zyy"} for value in components):
            raise ValueError(
                "2-D triangular problems support only zxy and zyx impedance components."
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
        >>> from pycsamt.forward.maxwell.contracts import ReceiverSet
        >>> mesh = TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]])
        >>> p = TriProblem(mesh, [1], [1], ReceiverSet([[0, 0]], ["S"]))
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
        >>> from pycsamt.forward.maxwell.contracts import ReceiverSet
        >>> mesh = TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]])
        >>> p = TriProblem(mesh, [1], [1], ReceiverSet([[0, 0]], ["S"]))
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
        >>> from pycsamt.forward.maxwell.contracts import ReceiverSet
        >>> mesh = TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]])
        >>> p = TriProblem(mesh, [1], [1], ReceiverSet([[0, 0]], ["S"]))
        >>> with TemporaryDirectory() as d:
        ...     restored = TriProblem.from_npz(p.to_npz(Path(d) / "p.npz"))
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
    def from_npz(cls, path: str | Path) -> TriProblem:
        """Restore and validate a problem archive.

        Parameters
        ----------
        path : str or pathlib.Path
            Archive written by :meth:`to_npz`.

        Returns
        -------
        TriProblem
            Restored problem.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> from pycsamt.forward.maxwell.contracts import ReceiverSet
        >>> mesh = TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]])
        >>> p = TriProblem(mesh, [1], [1], ReceiverSet([[0, 0]], ["S"]))
        >>> with TemporaryDirectory() as d:
        ...     q = TriProblem.from_npz(p.to_npz(Path(d) / "p.npz"))
        >>> np.array_equal(q.conductivity_s_m, p.conductivity_s_m)
        True
        """
        with np.load(Path(path), allow_pickle=False) as archive:
            state = json.loads(str(archive["provenance_json"].item()))
            if state.get("schema_version") != 1:
                raise ValueError("unsupported TriProblem schema version.")
            return cls(
                TriMesh.from_dict(state["mesh"]),
                archive["conductivity_s_m"],
                archive["frequencies_hz"],
                ReceiverSet.from_dict(state["receivers"]),
                tuple(state["components"]),
                archive["active_cells"],
                state["time_dependence"],
                state["magnetic_permeability_h_m"],
                state.get("metadata", {}),
            )
