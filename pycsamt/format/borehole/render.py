# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Viewer-neutral render models for PCBH boreholes."""

from __future__ import annotations

import colorsys
import hashlib
import math
from dataclasses import dataclass
from typing import Any

from ...api.property import PyCSAMTObject
from .schema import PCBHBorehole, PCBHDocument, VocabularyEntry
from .trajectory import (
    DesurveyedTrajectory,
    TrajectoryPoint,
    desurvey,
    trajectory_checksum,
)

__all__ = [
    "DisplayRadiusPolicy",
    "CollarPrimitive",
    "PolylinePrimitive",
    "IntervalSegmentPrimitive",
    "ContactPrimitive",
    "StructureGlyphPrimitive",
    "MaterialBatch",
    "BoreholeRenderModel",
    "PCBHRenderModel",
    "BoreholeRenderBuilder",
    "build_render_model",
    "deterministic_color",
]

_BOUNDARY_EPS = 1e-9


@dataclass(frozen=True, repr=False)
class DisplayRadiusPolicy(PyCSAMTObject):
    """View-only borehole radius policy in document coordinate units."""

    mode: str = "auto"
    fixed_radius: float | None = None
    exaggeration: float = 1.0
    extent_fraction: float = 0.001
    minimum_radius: float = 0.05

    def __post_init__(self) -> None:
        if self.mode not in {"auto", "fixed", "exaggeration"}:
            raise ValueError(
                "radius mode must be auto, fixed, or exaggeration"
            )
        for name in ("exaggeration", "extent_fraction", "minimum_radius"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a finite number > 0")
            object.__setattr__(self, name, value)
        if self.mode == "fixed":
            if self.fixed_radius is None:
                raise ValueError("fixed mode requires fixed_radius")
            radius = float(self.fixed_radius)
            if not math.isfinite(radius) or radius <= 0:
                raise ValueError("fixed_radius must be a finite number > 0")
            object.__setattr__(self, "fixed_radius", radius)

    def resolve(self, diameter: float | None, model_extent: float) -> float:
        """Resolve display radius without modifying physical diameter."""
        automatic = max(
            float(model_extent) * self.extent_fraction,
            self.minimum_radius,
        )
        physical = float(diameter) / 2.0 if diameter is not None else automatic
        if self.mode == "fixed":
            return float(self.fixed_radius)
        if self.mode == "exaggeration":
            return max(physical * self.exaggeration, self.minimum_radius)
        return max(physical, automatic)


@dataclass(frozen=True, repr=False)
class CollarPrimitive(PyCSAMTObject):
    """Collar marker, label, selection state, and hover metadata."""

    borehole_id: str
    position: tuple[float, float, float]
    label: str
    selected: bool
    metadata: dict[str, Any]


@dataclass(frozen=True, repr=False)
class PolylinePrimitive(PyCSAMTObject):
    """One borehole centerline in document coordinates."""

    borehole_id: str
    points: tuple[TrajectoryPoint, ...]
    color: str
    selected: bool
    metadata: dict[str, Any]


@dataclass(frozen=True, repr=False)
class IntervalSegmentPrimitive(PyCSAMTObject):
    """Colored portion of a centerline for one logged interval."""

    segment_id: str
    borehole_id: str
    family: str
    from_md: float
    to_md: float
    points: tuple[TrajectoryPoint, ...]
    color: str
    display_radius: float
    selected: bool
    metadata: dict[str, Any]


@dataclass(frozen=True, repr=False)
class ContactPrimitive(PyCSAMTObject):
    """Contact-ring input located at an exact shared interval boundary."""

    contact_id: str
    borehole_id: str
    family: str
    md: float
    position: tuple[float, float, float]
    tangent: tuple[float, float, float]
    display_radius: float
    left_code: str | None
    right_code: str | None
    selected: bool
    metadata: dict[str, Any]


@dataclass(frozen=True, repr=False)
class StructureGlyphPrimitive(PyCSAMTObject):
    """Position and orientation inputs for one structural glyph."""

    glyph_id: str
    borehole_id: str
    kind: str
    md: float
    position: tuple[float, float, float]
    borehole_tangent: tuple[float, float, float]
    representation: str
    orientation: dict[str, float]
    display_radius: float
    selected: bool
    metadata: dict[str, Any]


@dataclass(frozen=True, repr=False)
class MaterialBatch(PyCSAMTObject):
    """Stable grouping of segments sharing one family and material color."""

    family: str
    color: str
    segment_ids: tuple[str, ...]


@dataclass(frozen=True, repr=False)
class BoreholeRenderModel(PyCSAMTObject):
    """Complete viewer-neutral primitives for one borehole."""

    borehole_id: str
    collar: CollarPrimitive
    centerline: PolylinePrimitive
    interval_segments: tuple[IntervalSegmentPrimitive, ...]
    contacts: tuple[ContactPrimitive, ...]
    structure_glyphs: tuple[StructureGlyphPrimitive, ...]
    source_checksum: str
    display_radius: float


@dataclass(frozen=True, repr=False)
class PCBHRenderModel(PyCSAMTObject):
    """Complete render contract for one PCBH document."""

    document_id: str
    crs_horizontal: str
    coordinate_unit: str
    boreholes: tuple[BoreholeRenderModel, ...]
    batches: tuple[MaterialBatch, ...]
    bounds: tuple[float, float, float, float, float, float]
    lod_tolerance: float
    sampling_step_md: float
    max_points_per_hole: int
    radius_policy: DisplayRadiusPolicy


class BoreholeRenderBuilder(PyCSAMTObject):
    """Build reusable render primitives while caching scientific paths."""

    def __init__(self) -> None:
        self._trajectory_cache: dict[str, DesurveyedTrajectory] = {}

    @property
    def cache_size(self) -> int:
        """Number of geometry checksums currently cached."""
        return len(self._trajectory_cache)

    def clear_cache(self) -> None:
        """Discard cached desurveyed centerlines."""
        self._trajectory_cache.clear()

    def build(
        self,
        document: PCBHDocument,
        *,
        family: str = "lithology",
        selected_ids: set[str] | None = None,
        radius_policy: DisplayRadiusPolicy | None = None,
        lod_tolerance: float = 0.0,
        sampling_step_md: float = 10.0,
        max_points_per_hole: int = 10_000,
    ) -> PCBHRenderModel:
        """Build a backend-independent render model for one PCBH document."""
        if not isinstance(document, PCBHDocument):
            raise TypeError("document must be a PCBHDocument")
        document.validate()
        if document.units.depth != document.crs.coordinate_unit:
            raise ValueError(
                "rendering requires depth and coordinate units to match; "
                "convert units before building geometry"
            )
        if not isinstance(family, str) or not family.strip():
            raise ValueError("family must be a non-empty string")
        tolerance = float(lod_tolerance)
        if not math.isfinite(tolerance) or tolerance < 0:
            raise ValueError("lod_tolerance must be a finite number >= 0")
        sampling_step = float(sampling_step_md)
        if not math.isfinite(sampling_step) or sampling_step <= 0:
            raise ValueError("sampling_step_md must be a finite number > 0")
        if (
            isinstance(max_points_per_hole, bool)
            or not isinstance(max_points_per_hole, int)
            or max_points_per_hole < 2
        ):
            raise ValueError("max_points_per_hole must be an integer >= 2")
        selected = set(selected_ids or ())
        policy = radius_policy or DisplayRadiusPolicy()
        if not isinstance(policy, DisplayRadiusPolicy):
            raise TypeError("radius_policy must be a DisplayRadiusPolicy")

        base_paths = {
            hole.id: self._trajectory(hole) for hole in document.boreholes
        }
        display_paths = {
            hole.id: base_paths[hole.id].split_at(
                _regular_depths(
                    hole.total_depth_md,
                    sampling_step,
                    max_points_per_hole,
                )
            )
            for hole in document.boreholes
        }
        bounds = _bounds(display_paths.values())
        extent = max(
            bounds[1] - bounds[0],
            bounds[3] - bounds[2],
            bounds[5] - bounds[4],
        )
        colors = _vocabulary_colors(document.lithologies)
        holes = tuple(
            self._build_hole(
                hole,
                display_paths[hole.id],
                family=family,
                selected=hole.id in selected,
                display_radius=policy.resolve(hole.diameter, extent),
                colors=colors,
                lod_tolerance=tolerance,
            )
            for hole in document.boreholes
        )
        return PCBHRenderModel(
            document_id=document.document_id,
            crs_horizontal=document.crs.horizontal,
            coordinate_unit=document.crs.coordinate_unit,
            boreholes=holes,
            batches=_batches(holes),
            bounds=bounds,
            lod_tolerance=tolerance,
            sampling_step_md=sampling_step,
            max_points_per_hole=max_points_per_hole,
            radius_policy=policy,
        )

    def _trajectory(self, hole: PCBHBorehole) -> DesurveyedTrajectory:
        checksum = trajectory_checksum(hole)
        if checksum not in self._trajectory_cache:
            self._trajectory_cache[checksum] = desurvey(hole)
        return self._trajectory_cache[checksum]

    def _build_hole(
        self,
        hole: PCBHBorehole,
        base: DesurveyedTrajectory,
        *,
        family: str,
        selected: bool,
        display_radius: float,
        colors: dict[str, str],
        lod_tolerance: float,
    ) -> BoreholeRenderModel:
        intervals = sorted(
            hole.interval_logs.get(family, []),
            key=lambda item: item.from_md,
        )
        boundaries = {
            depth
            for item in intervals
            for depth in (float(item.from_md), float(item.to_md))
        }
        split = base.split_at(boundaries)
        center_points = _simplify_protected(
            split.points, boundaries, lod_tolerance
        )
        segments = []
        for index, interval in enumerate(intervals):
            code = interval.code or interval.label or "unknown"
            color = colors.get(code) or colors.get(
                (interval.label or "").casefold()
            )
            color = color or deterministic_color(code)
            points = _points_between(
                split.points, interval.from_md, interval.to_md
            )
            points = _simplify(points, lod_tolerance)
            segment_id = f"{hole.id}:{family}:{index}"
            segments.append(
                IntervalSegmentPrimitive(
                    segment_id=segment_id,
                    borehole_id=hole.id,
                    family=family,
                    from_md=float(interval.from_md),
                    to_md=float(interval.to_md),
                    points=points,
                    color=color,
                    display_radius=display_radius,
                    selected=selected,
                    metadata={
                        "borehole_id": hole.id,
                        "family": family,
                        "from_md": interval.from_md,
                        "to_md": interval.to_md,
                        "code": interval.code,
                        "label": interval.label,
                        "resistivity_ohm_m": interval.resistivity_ohm_m,
                        "description": interval.description,
                        "properties": dict(interval.properties),
                    },
                )
            )
        contacts = _contacts(
            hole,
            split,
            intervals,
            family,
            display_radius,
            selected,
        )
        glyphs = tuple(
            _structure_glyph(
                hole,
                base,
                item,
                index,
                display_radius,
                selected,
            )
            for index, item in enumerate(hole.structures)
        )
        collar_metadata = {
            "borehole_id": hole.id,
            "name": hole.name,
            "kind": hole.kind,
            "status": hole.status,
            "total_depth_md": hole.total_depth_md,
            "collar_elevation": hole.collar.z,
        }
        return BoreholeRenderModel(
            borehole_id=hole.id,
            collar=CollarPrimitive(
                borehole_id=hole.id,
                position=(hole.collar.x, hole.collar.y, hole.collar.z),
                label=hole.name,
                selected=selected,
                metadata=collar_metadata,
            ),
            centerline=PolylinePrimitive(
                borehole_id=hole.id,
                points=center_points,
                color="#111827",
                selected=selected,
                metadata=collar_metadata,
            ),
            interval_segments=tuple(segments),
            contacts=contacts,
            structure_glyphs=glyphs,
            source_checksum=base.source_checksum,
            display_radius=display_radius,
        )


def build_render_model(
    document: PCBHDocument,
    *,
    family: str = "lithology",
    selected_ids: set[str] | None = None,
    radius_policy: DisplayRadiusPolicy | None = None,
    lod_tolerance: float = 0.0,
    sampling_step_md: float = 10.0,
    max_points_per_hole: int = 10_000,
) -> PCBHRenderModel:
    """Build a render model with a short-lived builder.

    Parameters mirror :meth:`BoreholeRenderBuilder.build`. Distances are in
    the document coordinate/depth unit and no CRS transformation is applied.
    """
    return BoreholeRenderBuilder().build(
        document,
        family=family,
        selected_ids=selected_ids,
        radius_policy=radius_policy,
        lod_tolerance=lod_tolerance,
        sampling_step_md=sampling_step_md,
        max_points_per_hole=max_points_per_hole,
    )


def deterministic_color(value: str) -> str:
    """Return a process-independent, readable color for a vocabulary key."""
    digest = hashlib.sha256(str(value).casefold().encode("utf-8")).digest()
    hue = int.from_bytes(digest[:2], "big") / 65535.0
    red, green, blue = colorsys.hsv_to_rgb(hue, 0.58, 0.78)
    components = (round(red * 255), round(green * 255), round(blue * 255))
    return "#{:02X}{:02X}{:02X}".format(*components)


def _vocabulary_colors(entries: list[VocabularyEntry]) -> dict[str, str]:
    colors = {}
    for entry in entries:
        color = entry.color or deterministic_color(entry.code)
        colors[entry.code] = color
        colors[entry.name.casefold()] = color
    return colors


def _bounds(paths):
    points = [point for path in paths for point in path.points]
    return (
        min(point.x for point in points),
        max(point.x for point in points),
        min(point.y for point in points),
        max(point.y for point in points),
        min(point.z for point in points),
        max(point.z for point in points),
    )


def _regular_depths(
    total_depth: float, step: float, max_points: int
) -> tuple[float, ...]:
    effective_step = max(step, float(total_depth) / (max_points - 1))
    count = int(float(total_depth) // effective_step)
    return tuple(
        index * effective_step
        for index in range(1, count + 1)
        if index * effective_step < total_depth
    )


def _points_between(points, start, end):
    return tuple(
        point
        for point in points
        if start - _BOUNDARY_EPS <= point.md <= end + _BOUNDARY_EPS
    )


def _tangent(point: TrajectoryPoint) -> tuple[float, float, float]:
    azimuth = math.radians(point.azimuth_deg)
    inclination = math.radians(point.inclination_deg)
    horizontal = math.sin(inclination)
    return (
        horizontal * math.sin(azimuth),
        horizontal * math.cos(azimuth),
        -math.cos(inclination),
    )


def _contacts(hole, path, intervals, family, radius, selected):
    contacts = []
    for index in range(len(intervals) - 1):
        left, right = intervals[index], intervals[index + 1]
        if not math.isclose(left.to_md, right.from_md, abs_tol=_BOUNDARY_EPS):
            continue
        point = path.at_md(left.to_md)
        contacts.append(
            ContactPrimitive(
                contact_id=f"{hole.id}:{family}:contact:{index}",
                borehole_id=hole.id,
                family=family,
                md=point.md,
                position=(point.x, point.y, point.z),
                tangent=_tangent(point),
                display_radius=radius * 1.15,
                left_code=left.code,
                right_code=right.code,
                selected=selected,
                metadata={
                    "borehole_id": hole.id,
                    "md": point.md,
                    "tvd": point.tvd,
                    "elevation": point.z,
                    "left_label": left.label,
                    "right_label": right.label,
                },
            )
        )
    return tuple(contacts)


def _structure_glyph(hole, path, item, index, radius, selected):
    md = (
        float(item.at_md)
        if item.at_md is not None
        else (float(item.from_md) + float(item.to_md)) / 2.0
    )
    point = path.at_md(md)
    orientation = {
        name: float(value)
        for name in (
            "strike_deg",
            "dip_deg",
            "dip_direction_deg",
            "trend_deg",
            "plunge_deg",
            "alpha_deg",
            "beta_deg",
        )
        if (value := getattr(item, name)) is not None
    }
    return StructureGlyphPrimitive(
        glyph_id=f"{hole.id}:structure:{index}",
        borehole_id=hole.id,
        kind=item.kind,
        md=md,
        position=(point.x, point.y, point.z),
        borehole_tangent=_tangent(point),
        representation=item.orientation_representation,
        orientation=orientation,
        display_radius=radius * 2.0,
        selected=selected,
        metadata={
            "borehole_id": hole.id,
            "kind": item.kind,
            "md": md,
            "tvd": point.tvd,
            "elevation": point.z,
            "from_md": item.from_md,
            "to_md": item.to_md,
            "aperture_m": item.aperture_m,
            "fill": item.fill,
            "confidence": item.confidence,
        },
    )


def _batches(holes):
    grouped: dict[tuple[str, str], list[str]] = {}
    for hole in holes:
        for segment in hole.interval_segments:
            grouped.setdefault((segment.family, segment.color), []).append(
                segment.segment_id
            )
    return tuple(
        MaterialBatch(family, color, tuple(segment_ids))
        for (family, color), segment_ids in sorted(grouped.items())
    )


def _simplify_protected(points, protected, tolerance):
    if tolerance <= 0 or len(points) <= 2:
        return tuple(points)
    indexes = [0]
    indexes.extend(
        index
        for index, point in enumerate(points[1:-1], start=1)
        if any(
            math.isclose(point.md, md, abs_tol=_BOUNDARY_EPS)
            for md in protected
        )
    )
    indexes.append(len(points) - 1)
    output = []
    for start, end in zip(indexes, indexes[1:]):
        simplified = _simplify(tuple(points[start : end + 1]), tolerance)
        output.extend(simplified if not output else simplified[1:])
    return tuple(output)


def _simplify(points, tolerance):
    points = tuple(points)
    if tolerance <= 0 or len(points) <= 2:
        return points
    start, end = points[0], points[-1]
    maximum, split = -1.0, 0
    for index, point in enumerate(points[1:-1], start=1):
        distance = _point_line_distance(point, start, end)
        if distance > maximum:
            maximum, split = distance, index
    if maximum <= tolerance:
        return (start, end)
    left = _simplify(points[: split + 1], tolerance)
    right = _simplify(points[split:], tolerance)
    return left[:-1] + right


def _point_line_distance(point, start, end):
    vector = (end.x - start.x, end.y - start.y, end.z - start.z)
    offset = (point.x - start.x, point.y - start.y, point.z - start.z)
    denominator = sum(value * value for value in vector)
    if denominator == 0:
        return math.sqrt(sum(value * value for value in offset))
    fraction = max(
        0.0,
        min(1.0, sum(a * b for a, b in zip(offset, vector)) / denominator),
    )
    residual = tuple(a - fraction * b for a, b in zip(offset, vector))
    return math.sqrt(sum(value * value for value in residual))
