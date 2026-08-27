# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""In-memory schema for the pyCSAMT Common Subsurface Format (PCSF).

PCSF is the backend-neutral container every inversion result (Occam2D,
ModEM, MARE2DEM, DUHI) can be converted to, so :mod:`pycsamt.format.io`,
``app/mapview``, and the web 3-D view share one on-disk representation
instead of five ad hoc in-memory shapes. See
``PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md`` at the repository root for the
full design rationale.

This module defines the schema only — geometry and container
dataclasses, the frozen list of supported geometry kinds, and shape
validation. Reading/writing ``.pcsf`` files lives in
:mod:`pycsamt.format.io`; per-backend conversion lives in
:mod:`pycsamt.format.adapters` (Phases 2-4).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, ClassVar

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject

__all__ = [
    "PCSF_VERSION",
    "RESISTIVITY_UNIT",
    "GEOMETRY_KINDS",
    "RESISTIVITY_ENCODINGS",
    "TOPOGRAPHY_KINDS",
    "DERIVATION_METHODS",
    "Grid2DGeometry",
    "Grid3DGeometry",
    "UnstructuredMeshGeometry",
    "LineEntry",
    "DerivedVolume",
    "MultilineGeometry",
    "StationTable",
    "TopographyPerStation",
    "TopographyRaster",
    "PCSFModel",
]

PCSF_VERSION = "0.1.0"
RESISTIVITY_UNIT = "ohm.m"
GEOMETRY_KINDS = ("grid2d", "grid3d", "mesh_unstructured", "multiline")
RESISTIVITY_ENCODINGS = ("log10", "ln", "linear")
TOPOGRAPHY_KINDS = ("per_station", "raster")
DERIVATION_METHODS = ("linear_interp", "kriging", "idw")


def _arr(value: Any, *, dtype: Any = float) -> np.ndarray:
    return np.asarray(value, dtype=dtype)


def _opt_arr(value: Any, *, dtype: Any = float) -> np.ndarray | None:
    return None if value is None else np.asarray(value, dtype=dtype)


# ---------------------------------------------------------------------
# Geometry kinds
# ---------------------------------------------------------------------


@dataclass(repr=False)
class Grid2DGeometry(PyCSAMTObject):
    """Single-profile rectilinear geometry (Occam2D / DUHI-via-Occam2D).

    Parameters
    ----------
    x : ndarray (n_x,)
        Real station chainage, metres — never a solver's mesh-local
        frame (see the Occam2D coordinate-frame note in
        :meth:`pycsamt.interp._base.ResistivityModel.from_occam2d`).
    z : ndarray (n_z,)
        Depth cell centres, metres, positive downward.
    x_nodes, z_nodes : ndarray, optional
        Cell-edge coordinates, one longer than *x*/*z*.
    origin : ndarray (2,), optional
        Real-world offset when *x* is locally referenced.
    azimuth_deg : float, optional
        Profile bearing, for georeferencing back to the survey line.
    """

    kind: ClassVar[str] = "grid2d"

    x: np.ndarray
    z: np.ndarray
    x_nodes: np.ndarray | None = None
    z_nodes: np.ndarray | None = None
    origin: np.ndarray | None = None
    azimuth_deg: float | None = None

    def __post_init__(self) -> None:
        self.x = _arr(self.x)
        self.z = _arr(self.z)
        self.x_nodes = _opt_arr(self.x_nodes)
        self.z_nodes = _opt_arr(self.z_nodes)
        self.origin = _opt_arr(self.origin)

    def validate(self) -> None:
        if self.x.ndim != 1 or self.z.ndim != 1:
            raise ValueError("Grid2DGeometry.x and .z must be 1-D")
        if self.x_nodes is not None and self.x_nodes.shape != (
            self.x.shape[0] + 1,
        ):
            raise ValueError("x_nodes must have one more element than x")
        if self.z_nodes is not None and self.z_nodes.shape != (
            self.z.shape[0] + 1,
        ):
            raise ValueError("z_nodes must have one more element than z")
        if self.origin is not None and self.origin.shape != (2,):
            raise ValueError("origin must have shape (2,)")

    @property
    def resistivity_shape(self) -> tuple[int, int]:
        return (self.z.shape[0], self.x.shape[0])


@dataclass(repr=False)
class Grid3DGeometry(PyCSAMTObject):
    """Native 3-D tensor volume geometry (ModEM).

    Parameters
    ----------
    x, y, z : ndarray
        Cell-centre coordinates, metres.
    x_nodes, y_nodes, z_nodes : ndarray, optional
        Cell-edge coordinates.
    origin : ndarray (3,), optional
        Real-world grid origin.
    rotation_deg : float, default 0.0
        Grid rotation about the vertical axis.
    n_air : int, default 0
        Explicit air-layer count (unlike Occam2D's inferred count).
    """

    kind: ClassVar[str] = "grid3d"

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    x_nodes: np.ndarray | None = None
    y_nodes: np.ndarray | None = None
    z_nodes: np.ndarray | None = None
    origin: np.ndarray | None = None
    rotation_deg: float = 0.0
    n_air: int = 0

    def __post_init__(self) -> None:
        self.x = _arr(self.x)
        self.y = _arr(self.y)
        self.z = _arr(self.z)
        self.x_nodes = _opt_arr(self.x_nodes)
        self.y_nodes = _opt_arr(self.y_nodes)
        self.z_nodes = _opt_arr(self.z_nodes)
        self.origin = _opt_arr(self.origin)

    def validate(self) -> None:
        if self.x.ndim != 1 or self.y.ndim != 1 or self.z.ndim != 1:
            raise ValueError("Grid3DGeometry.x, .y, .z must be 1-D")
        if self.origin is not None and self.origin.shape != (3,):
            raise ValueError("origin must have shape (3,)")
        if self.n_air < 0:
            raise ValueError("n_air must be >= 0")

    @property
    def resistivity_shape(self) -> tuple[int, int, int]:
        # Deliberately (n_z, n_y, n_x) — matches ModEM's own native
        # order, to avoid introducing a second axis-order transpose bug.
        return (self.z.shape[0], self.y.shape[0], self.x.shape[0])


@dataclass(repr=False)
class UnstructuredMeshGeometry(PyCSAMTObject):
    """Native unstructured triangular mesh geometry (MARE2DEM).

    Preserves the mesh as-is (no forced regrid onto a tensor grid), so
    a MARE2DEM result keeps its real element resolution.

    Parameters
    ----------
    nodes : ndarray (n, 2) or (n, 3)
        Node coordinates, metres.
    connectivity : ndarray (m, 3), int
        Triangle node indices.
    region_ids : ndarray (m,), int
        Region id per triangle.
    plane : {"xz", "xy", "3d"}, default "xz"
        Physical plane the mesh lives in. MARE2DEM profiles are
        conventionally in ``(y, z)`` but stored generically as
        ``plane="xz"`` with *x* holding the profile's own along-line
        coordinate.
    """

    kind: ClassVar[str] = "mesh_unstructured"

    nodes: np.ndarray
    connectivity: np.ndarray
    region_ids: np.ndarray
    plane: str = "xz"

    def __post_init__(self) -> None:
        self.nodes = _arr(self.nodes)
        self.connectivity = _arr(self.connectivity, dtype=np.int64)
        self.region_ids = _arr(self.region_ids, dtype=np.int32)

    def validate(self) -> None:
        if self.nodes.ndim != 2 or self.nodes.shape[1] not in (2, 3):
            raise ValueError("nodes must have shape (n, 2) or (n, 3)")
        if self.connectivity.ndim != 2 or self.connectivity.shape[1] != 3:
            raise ValueError("connectivity must have shape (m, 3)")
        if self.region_ids.shape != (self.connectivity.shape[0],):
            raise ValueError(
                "region_ids must have one entry per triangle"
            )
        if self.plane not in ("xz", "xy", "3d"):
            raise ValueError(
                f"plane must be one of 'xz'/'xy'/'3d', got {self.plane!r}"
            )
        n_nodes = self.nodes.shape[0]
        if self.connectivity.size and (
            self.connectivity.min() < 0 or self.connectivity.max() >= n_nodes
        ):
            raise ValueError("connectivity references an out-of-range node")

    @property
    def n_regions(self) -> int:
        return int(np.unique(self.region_ids).size) if self.region_ids.size else 0


@dataclass(repr=False)
class LineEntry(PyCSAMTObject):
    """One profile within a :class:`MultilineGeometry`.

    Parameters
    ----------
    line_id : str
        Unique identifier for this line.
    geometry : Grid2DGeometry
        The line's own 2-D section geometry.
    resistivity : ndarray (n_z, n_x)
        Canonical linear ohm.m resistivity for this line.
    offset_y : float
        Cross-line position, metres.
    offset_kind : {"real", "synthetic"}, default "synthetic"
        Whether *offset_y* comes from real survey geometry or is a
        placeholder spacing for display only.
    azimuth_deg : float, optional
        Line bearing.
    """

    line_id: str
    geometry: Grid2DGeometry
    resistivity: np.ndarray
    offset_y: float = 0.0
    offset_kind: str = "synthetic"
    azimuth_deg: float | None = None

    def __post_init__(self) -> None:
        self.resistivity = _arr(self.resistivity)

    def validate(self) -> None:
        self.geometry.validate()
        if self.resistivity.shape != self.geometry.resistivity_shape:
            raise ValueError(
                f"line {self.line_id!r}: resistivity shape "
                f"{self.resistivity.shape} does not match geometry "
                f"{self.geometry.resistivity_shape}"
            )
        if self.offset_kind not in ("real", "synthetic"):
            raise ValueError("offset_kind must be 'real' or 'synthetic'")


@dataclass(repr=False)
class DerivedVolume(PyCSAMTObject):
    """Optional cached 3-D volume synthesized from stacked lines.

    Kept explicitly tagged as synthesized so a reader never mistakes a
    stack-interpolated volume for a native 3-D inversion (see
    ``derivation_method``/``synthesized`` in the design plan's §2).
    """

    grid: Grid3DGeometry
    resistivity: np.ndarray
    derivation_method: str = "linear_interp"
    derived_from: list[str] = field(default_factory=list)
    synthesized: bool = True

    def __post_init__(self) -> None:
        self.resistivity = _arr(self.resistivity)

    def validate(self) -> None:
        self.grid.validate()
        if self.resistivity.shape != self.grid.resistivity_shape:
            raise ValueError(
                f"derived_volume resistivity shape {self.resistivity.shape} "
                f"does not match grid {self.grid.resistivity_shape}"
            )
        if self.derivation_method not in DERIVATION_METHODS:
            raise ValueError(
                "derivation_method must be one of "
                f"{DERIVATION_METHODS}, got {self.derivation_method!r}"
            )


@dataclass(repr=False)
class MultilineGeometry(PyCSAMTObject):
    """A set of profiles plus real line geometry (fence/block views).

    Formalizes what ``pycsamt/app/web/callbacks/map3d.py`` currently
    reconstructs at render time from a stack of independent 2-D
    sections. Each line carries its own resistivity; the optional
    :attr:`derived_volume` is a documented, reproducible synthesis
    rather than a render-time-only side effect.
    """

    kind: ClassVar[str] = "multiline"

    lines: list[LineEntry] = field(default_factory=list)
    derived_volume: DerivedVolume | None = None

    def validate(self) -> None:
        if not self.lines:
            raise ValueError("MultilineGeometry needs at least one line")
        ids = [line.line_id for line in self.lines]
        if len(set(ids)) != len(ids):
            raise ValueError(f"duplicate line_id values: {ids}")
        for line in self.lines:
            line.validate()
        if self.derived_volume is not None:
            self.derived_volume.validate()


# ---------------------------------------------------------------------
# Container-level tables
# ---------------------------------------------------------------------


@dataclass(repr=False)
class StationTable(PyCSAMTObject):
    """Survey station positions, shared across geometry kinds.

    ``x``/``y``/``z`` are geometry-local (along-profile chainage for
    ``grid2d``, the model grid's own frame for ``grid3d``, whatever
    frame the caller supplied for ``mesh_unstructured``) — never
    assumed to be real-world geographic coordinates, per SPEC.md's own
    ``load_pcsf_lines`` convention. ``lon``/``lat``, when present, are
    the one explicit, unambiguous carrier of real-world position: WGS84
    decimal degrees, the same convention every other real-coordinate
    source in this codebase already uses (EDI headers, a ModEM ``.dat``
    file's ``GG_Lat``/``GG_Lon`` columns, :class:`pycsamt.map._core.StationRecord`).
    A single-line ``grid2d`` (or any other kind's) PCSF file that sets
    these needs no separate ``known_stations`` match to place its
    stations on a real basemap.
    """

    name: list[str] = field(default_factory=list)
    x: np.ndarray = field(default_factory=lambda: np.array([]))
    y: np.ndarray = field(default_factory=lambda: np.array([]))
    z: np.ndarray = field(default_factory=lambda: np.array([]))
    line_id: list[str] | None = None
    lon: np.ndarray | None = None
    lat: np.ndarray | None = None

    def __post_init__(self) -> None:
        self.x = _arr(self.x)
        self.y = _arr(self.y)
        self.z = _arr(self.z)
        self.lon = _opt_arr(self.lon)
        self.lat = _opt_arr(self.lat)

    def validate(self) -> None:
        n = len(self.name)
        for label, arr in (("x", self.x), ("y", self.y), ("z", self.z)):
            if arr.shape != (n,):
                raise ValueError(
                    f"StationTable.{label} must have shape ({n},) to "
                    f"match {n} station names, got {arr.shape}"
                )
        if self.line_id is not None and len(self.line_id) != n:
            raise ValueError("StationTable.line_id must match name length")
        for label, arr in (("lon", self.lon), ("lat", self.lat)):
            if arr is not None and arr.shape != (n,):
                raise ValueError(
                    f"StationTable.{label} must have shape ({n},) to "
                    f"match {n} station names, got {arr.shape}"
                )
        if (self.lon is None) != (self.lat is None):
            raise ValueError(
                "StationTable.lon and .lat must be set together (or "
                "both left unset), never only one"
            )


@dataclass(repr=False)
class TopographyPerStation(PyCSAMTObject):
    """Scalar-per-station topography (matches the existing convention
    in :mod:`pycsamt.map.topo`).
    """

    kind: ClassVar[str] = "per_station"

    station_id: list[str] = field(default_factory=list)
    elevation: np.ndarray = field(default_factory=lambda: np.array([]))

    def __post_init__(self) -> None:
        self.elevation = _arr(self.elevation)

    def validate(self) -> None:
        if self.elevation.shape != (len(self.station_id),):
            raise ValueError(
                "TopographyPerStation.elevation must match station_id length"
            )


@dataclass(repr=False)
class TopographyRaster(PyCSAMTObject):
    """Gridded-DEM topography — a regular elevation surface independent
    of any station table.

    Unlike :class:`TopographyPerStation`, this carries no station
    identifiers at all: it is a standalone terrain surface a consumer
    can sample at any coordinate, not a per-station lookup table. It
    introduces no GDAL/rasterio-class dependency — construction is via
    plain ``x``/``y``/``elevation`` arrays a caller has already
    obtained by whatever means it likes (see
    :func:`pycsamt.format.topography.topography_from_grid`); PCSF
    itself never parses a georeferenced raster file format.

    Parameters
    ----------
    x : ndarray (n_x,)
        Grid x-coordinates (or longitude), increasing.
    y : ndarray (n_y,)
        Grid y-coordinates (or latitude), increasing.
    elevation : ndarray (n_y, n_x)
        Elevation surface, metres, sampled on the ``(y, x)`` meshgrid
        implied by *x*/*y* — the same row-major convention
        ``numpy.meshgrid(x, y)`` (default ``indexing="xy"``) produces.
    """

    kind: ClassVar[str] = "raster"

    x: np.ndarray
    y: np.ndarray
    elevation: np.ndarray

    def __post_init__(self) -> None:
        self.x = _arr(self.x)
        self.y = _arr(self.y)
        self.elevation = _arr(self.elevation)

    def validate(self) -> None:
        if self.x.ndim != 1 or self.y.ndim != 1:
            raise ValueError("TopographyRaster.x and .y must be 1-D")
        expected = (self.y.shape[0], self.x.shape[0])
        if self.elevation.shape != expected:
            raise ValueError(
                f"TopographyRaster.elevation shape {self.elevation.shape} "
                f"does not match (len(y), len(x)) = {expected}"
            )


# ---------------------------------------------------------------------
# Root container
# ---------------------------------------------------------------------


@dataclass(repr=False)
class PCSFModel(PyCSAMTObject, MetadataMixin):
    """Backend-neutral inversion-result container (one PCSF file).

    Parameters
    ----------
    geometry : Grid2DGeometry | Grid3DGeometry | UnstructuredMeshGeometry | MultilineGeometry
        The model's geometry, discriminated by ``geometry.kind``.
    resistivity : ndarray, optional
        Canonical **linear ohm.m** resistivity. Required for
        ``grid2d``/``grid3d``/``mesh_unstructured`` geometries;
        must be ``None`` for ``multiline`` (each line carries its own
        resistivity — see :class:`LineEntry`).
    resistivity_native : ndarray, optional
        Passthrough of the source backend's own array, for provenance.
    resistivity_native_encoding : {"log10", "ln", "linear"}, optional
        Encoding of *resistivity_native*. Required whenever
        *resistivity_native* is given — never assumed.
    uncertainty, sensitivity : ndarray, optional
        Same shape as *resistivity*, when available from the source.
    resistivity_by_region : ndarray, optional
        Per-region resistivity table (``mesh_unstructured`` only),
        alongside the per-cell *resistivity* expanded from it.
    resistivity_by_node : ndarray, optional
        Per-node resistivity table (``mesh_unstructured`` only), shape
        ``(n_nodes,)`` -- the natural output shape of a graph-based
        model (e.g. a GCN) that predicts one value per mesh vertex
        rather than per cell. Kept alongside the per-cell *resistivity*
        expanded from it (see
        :func:`pycsamt.format.adapters.generic.mesh_to_pcsf`), the same
        provenance relationship *resistivity_by_region* has to its own
        per-cell expansion.
    stations : StationTable, optional
    topography : TopographyPerStation | TopographyRaster, optional
    survey : dict
        Free-form survey metadata. Adapters populate this from
        :mod:`pycsamt.metadata` objects (``SurveyMeta``, ``BBox``,
        ``ProvenanceMeta``) via their own dict conversion; PCSF itself
        does not require a specific metadata class here.
    history : dict of ndarray
        Optional per-iteration series (e.g. ``{"rms": ..., "lambda": ...}``)
        from an ``InversionHistory``-like source.
    source_backend : str, default "generic"
        ``"occam2d"`` | ``"modem3d"`` | ``"mare2dem"`` | ``"duhi"`` | ``"generic"``.
    crs : str, optional
        A pyproj-compatible CRS string.
    """

    geometry: (
        Grid2DGeometry
        | Grid3DGeometry
        | UnstructuredMeshGeometry
        | MultilineGeometry
    )
    resistivity: np.ndarray | None = None
    resistivity_native: np.ndarray | None = None
    resistivity_native_encoding: str | None = None
    uncertainty: np.ndarray | None = None
    sensitivity: np.ndarray | None = None
    resistivity_by_region: np.ndarray | None = None
    resistivity_by_node: np.ndarray | None = None
    stations: StationTable | None = None
    topography: TopographyPerStation | TopographyRaster | None = None
    survey: dict[str, Any] = field(default_factory=dict)
    history: dict[str, np.ndarray] = field(default_factory=dict)
    source_backend: str = "generic"
    created_by: str = ""
    created_at: str = ""
    crs: str | None = None
    description: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.resistivity = _opt_arr(self.resistivity)
        self.resistivity_native = _opt_arr(self.resistivity_native)
        self.uncertainty = _opt_arr(self.uncertainty)
        self.sensitivity = _opt_arr(self.sensitivity)
        self.resistivity_by_region = _opt_arr(self.resistivity_by_region)
        self.resistivity_by_node = _opt_arr(self.resistivity_by_node)
        self.history = {
            str(key): _arr(value) for key, value in dict(self.history).items()
        }

    @property
    def kind(self) -> str:
        return self.geometry.kind

    def validate(self) -> None:
        kind = self.geometry.kind
        if kind not in GEOMETRY_KINDS:
            raise ValueError(
                f"unknown geometry kind {kind!r}; expected one of "
                f"{GEOMETRY_KINDS}"
            )
        self.geometry.validate()

        if kind == "multiline":
            if self.resistivity is not None:
                raise ValueError(
                    "multiline geometry carries resistivity per line "
                    "(PCSFModel.resistivity must be None); set it on "
                    "each LineEntry instead"
                )
        else:
            if self.resistivity is None:
                raise ValueError(
                    f"resistivity is required for geometry kind {kind!r}"
                )
            if kind == "mesh_unstructured":
                # Per-cell (m,) or the region-collapsed (n_regions,) form.
                valid_shapes = {
                    (self.geometry.connectivity.shape[0],),
                    (self.geometry.n_regions,),
                }
                if self.resistivity.shape not in valid_shapes:
                    raise ValueError(
                        f"resistivity shape {self.resistivity.shape} does "
                        f"not match mesh cells or region count {valid_shapes}"
                    )
            else:
                expected = self.geometry.resistivity_shape
                if self.resistivity.shape != expected:
                    raise ValueError(
                        f"resistivity shape {self.resistivity.shape} does "
                        f"not match geometry {expected}"
                    )
            for label, arr in (
                ("uncertainty", self.uncertainty),
                ("sensitivity", self.sensitivity),
            ):
                if arr is not None and arr.shape != self.resistivity.shape:
                    raise ValueError(
                        f"{label} shape {arr.shape} does not match "
                        f"resistivity shape {self.resistivity.shape}"
                    )

        if self.resistivity_by_node is not None:
            if kind != "mesh_unstructured":
                raise ValueError(
                    "resistivity_by_node is only valid for geometry kind "
                    f"'mesh_unstructured', got {kind!r}"
                )
            n_nodes = self.geometry.nodes.shape[0]
            if self.resistivity_by_node.shape != (n_nodes,):
                raise ValueError(
                    "resistivity_by_node shape "
                    f"{self.resistivity_by_node.shape} does not match "
                    f"the mesh's node count ({n_nodes},)"
                )

        if self.resistivity_native is not None and (
            self.resistivity_native_encoding not in RESISTIVITY_ENCODINGS
        ):
            raise ValueError(
                "resistivity_native_encoding must be one of "
                f"{RESISTIVITY_ENCODINGS} when resistivity_native is set, "
                f"got {self.resistivity_native_encoding!r}"
            )

        if self.stations is not None:
            self.stations.validate()
        if self.topography is not None:
            self.topography.validate()
