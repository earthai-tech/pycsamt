# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""PCBH association and coordinate alignment for PCSF models."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

from ...api.property import PyCSAMTObject
from .jsonio import pcbh_from_dict, pcbh_to_dict
from .render import build_render_model
from .schema import PCBHDocument

if TYPE_CHECKING:
    from ..schema import PCSFModel

__all__ = [
    "PCBHReference",
    "PCBHAssociation",
    "AlignedTrajectoryPoint",
    "BoreholeBoundsResult",
    "PCBHAlignmentReport",
    "pcbh_document_checksum",
    "embed_pcbh",
    "reference_pcbh",
    "extract_pcbh",
    "align_pcbh_to_pcsf",
]


@dataclass(frozen=True, repr=False)
class PCBHReference(PyCSAMTObject):
    """External PCBH URI and integrity checksum."""

    uri: str
    sha256: str

    def validate(self) -> None:
        if not isinstance(self.uri, str) or not self.uri.strip():
            raise ValueError("PCBH reference URI must be non-empty")
        digest = self.sha256.lower()
        if len(digest) != 64 or any(
            c not in "0123456789abcdef" for c in digest
        ):
            raise ValueError(
                "PCBH reference sha256 must be 64 hexadecimal characters"
            )


@dataclass(frozen=True, repr=False)
class PCBHAssociation(PyCSAMTObject):
    """Embedded and/or externally referenced PCBH attachment."""

    embedded: PCBHDocument | None = None
    reference: PCBHReference | None = None

    def validate(self) -> None:
        if self.embedded is None and self.reference is None:
            raise ValueError(
                "PCBH association needs embedded data or a reference"
            )
        if self.embedded is not None:
            self.embedded.validate()
        if self.reference is not None:
            self.reference.validate()
        if self.embedded is not None and self.reference is not None:
            checksum = pcbh_document_checksum(self.embedded)
            if checksum != self.reference.sha256.lower():
                raise ValueError(
                    "embedded PCBH does not match reference checksum"
                )


@dataclass(frozen=True, repr=False)
class AlignedTrajectoryPoint(PyCSAMTObject):
    """A trajectory point expressed in PCSF model coordinates."""

    md: float
    x: float
    y: float
    z: float


@dataclass(frozen=True, repr=False)
class BoreholeBoundsResult(PyCSAMTObject):
    """Per-hole relationship with the model bounding box."""

    borehole_id: str
    relation: str
    points_inside: int
    points_total: int


@dataclass(frozen=True, repr=False)
class PCBHAlignmentReport(PyCSAMTObject):
    """Aligned paths plus compatibility and bounds diagnostics."""

    source_crs: str
    target_crs: str
    vertical_compatible: bool
    warnings: tuple[str, ...]
    model_bounds: tuple[float, float, float, float, float, float]
    trajectories: dict[str, tuple[AlignedTrajectoryPoint, ...]]
    bounds_results: tuple[BoreholeBoundsResult, ...]


def pcbh_document_checksum(document: PCBHDocument) -> str:
    """Return SHA-256 of canonical compact PCBH JSON."""
    payload = json.dumps(
        pcbh_to_dict(document),
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def embed_pcbh(
    model: PCSFModel,
    document: PCBHDocument,
    *,
    uri: str | None = None,
) -> PCSFModel:
    """Return a model copy with a validated portable PCBH attachment."""
    reference = None
    if uri is not None:
        reference = PCBHReference(uri, pcbh_document_checksum(document))
    association = PCBHAssociation(document, reference)
    association.validate()
    return replace(model, boreholes=association)


def reference_pcbh(model: PCSFModel, uri: str, sha256: str) -> PCSFModel:
    """Return a model copy with an external PCBH reference."""
    association = PCBHAssociation(reference=PCBHReference(uri, sha256))
    association.validate()
    return replace(model, boreholes=association)


def extract_pcbh(model: PCSFModel) -> PCBHDocument | None:
    """Return embedded PCBH data without resolving external resources."""
    association = model.boreholes
    return association.embedded if association is not None else None


def align_pcbh_to_pcsf(
    document: PCBHDocument,
    model: PCSFModel,
    *,
    vertical_offset: float | None = None,
    sampling_step_md: float = 10.0,
) -> PCBHAlignmentReport:
    """Transform PCBH paths into PCSF x/y/depth coordinates.

    ``vertical_offset`` is required when PCBH and PCSF vertical references
    cannot be established as identical. It is added to PCBH elevations before
    conversion to PCSF depth-positive-down coordinates.
    """
    model.validate()
    if not model.crs:
        raise ValueError("PCSF model CRS is required for PCBH alignment")
    render = build_render_model(document, sampling_step_md=sampling_step_md)
    bounds, origin, rotation = _model_frame(model)
    vertical_compatible = _vertical_crs_compatible(
        document.crs.vertical,
        model.crs,
    )
    warnings = []
    if not vertical_compatible and vertical_offset is None:
        raise ValueError(
            "PCBH and PCSF vertical datums are not demonstrably compatible; "
            "provide an explicit vertical_offset"
        )
    if vertical_offset is not None:
        warnings.append("explicit vertical_offset applied")
    transformer = _horizontal_transformer(document.crs.horizontal, model.crs)
    trajectories = {}
    results = []
    for hole in render.boreholes:
        points = tuple(
            _align_point(
                point,
                transformer=transformer,
                origin=origin,
                rotation_deg=rotation,
                vertical_offset=float(vertical_offset or 0.0),
            )
            for point in hole.centerline.points
        )
        trajectories[hole.borehole_id] = points
        inside = sum(_inside(point, bounds) for point in points)
        if inside == len(points):
            relation = "inside"
        elif inside == 0 and not _path_intersects_box(points, bounds):
            relation = "outside"
        else:
            relation = "intersects"
        results.append(
            BoreholeBoundsResult(
                hole.borehole_id,
                relation,
                inside,
                len(points),
            )
        )
    return PCBHAlignmentReport(
        document.crs.horizontal,
        model.crs,
        vertical_compatible and vertical_offset is None,
        tuple(warnings),
        bounds,
        trajectories,
        tuple(results),
    )


def _horizontal_transformer(source: str, target: str):
    try:
        from pyproj import CRS, Transformer
    except ImportError as error:
        if source != target:
            raise ImportError(
                "pyproj is required for PCBH CRS transformation"
            ) from error
        return None
    source_crs = CRS.from_user_input(source)
    target_crs = CRS.from_user_input(target)
    if source_crs == target_crs:
        return None
    return Transformer.from_crs(source_crs, target_crs, always_xy=True)


def _vertical_crs_compatible(source_vertical: str, target: str) -> bool:
    if source_vertical in {"", "unknown", None}:
        return False
    try:
        from pyproj import CRS

        source = CRS.from_user_input(source_vertical)
        target_crs = CRS.from_user_input(target)
    except (ImportError, ValueError):
        return False
    candidates = [target_crs, *target_crs.sub_crs_list]
    return any(item.is_vertical and item == source for item in candidates)


def _align_point(point, *, transformer, origin, rotation_deg, vertical_offset):
    east, north = point.x, point.y
    if transformer is not None:
        east, north = transformer.transform(east, north)
    dx, dy = east - origin[0], north - origin[1]
    angle = math.radians(-rotation_deg)
    x = dx * math.cos(angle) - dy * math.sin(angle)
    y = dx * math.sin(angle) + dy * math.cos(angle)
    depth = origin[2] - (point.z + vertical_offset)
    return AlignedTrajectoryPoint(point.md, x, y, depth)


def _axis_bounds(values, nodes=None):
    data = nodes if nodes is not None else values
    return float(min(data)), float(max(data))


def _model_frame(model):
    geometry = model.geometry
    if geometry.kind == "multiline":
        if geometry.derived_volume is None:
            raise ValueError("multiline alignment requires a derived volume")
        geometry = geometry.derived_volume.grid
    if geometry.kind != "grid3d":
        raise ValueError("PCBH block alignment requires grid3d geometry")
    xb = _axis_bounds(geometry.x, geometry.x_nodes)
    yb = _axis_bounds(geometry.y, geometry.y_nodes)
    zb = _axis_bounds(geometry.z, geometry.z_nodes)
    origin = (
        tuple(geometry.origin)
        if geometry.origin is not None
        else (0.0, 0.0, 0.0)
    )
    bounds = (xb[0], xb[1], yb[0], yb[1], zb[0], zb[1])
    return bounds, origin, geometry.rotation_deg


def _inside(point, bounds):
    return (
        bounds[0] <= point.x <= bounds[1]
        and bounds[2] <= point.y <= bounds[3]
        and bounds[4] <= point.z <= bounds[5]
    )


def _path_intersects_box(points, bounds):
    for left, right in zip(points, points[1:]):
        for step in range(1, 20):
            fraction = step / 20.0
            candidate = AlignedTrajectoryPoint(
                left.md + fraction * (right.md - left.md),
                left.x + fraction * (right.x - left.x),
                left.y + fraction * (right.y - left.y),
                left.z + fraction * (right.z - left.z),
            )
            if _inside(candidate, bounds):
                return True
    return False


def association_to_dict(association: PCBHAssociation) -> dict[str, Any]:
    """Serialize an association for the PCSF HDF5 adapter."""
    association.validate()
    return {
        "embedded": pcbh_to_dict(association.embedded)
        if association.embedded is not None
        else None,
        "reference": {
            "uri": association.reference.uri,
            "sha256": association.reference.sha256,
        }
        if association.reference is not None
        else None,
    }


def association_from_dict(payload: dict[str, Any]) -> PCBHAssociation:
    """Deserialize an association from the PCSF HDF5 adapter."""
    embedded = payload.get("embedded")
    reference = payload.get("reference")
    result = PCBHAssociation(
        embedded=pcbh_from_dict(embedded) if embedded is not None else None,
        reference=(
            PCBHReference(**reference) if reference is not None else None
        ),
    )
    result.validate()
    return result
