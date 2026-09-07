# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Place PCBH boreholes in the 3-D map's local scene coordinates.

The fence / block / depth-slice / iso-surface builders in
:mod:`pycsamt.map.volume` do not work in real easting/northing. They
work in a synthetic *profile space*: ``x`` is along-strike distance,
``y`` is cross-strike line offset, and ``z`` is elevation with depth
running downward from a ``z = 0`` (or draped-topography) datum.

A PCBH document, by contrast, carries absolute collar coordinates in its
own CRS. Dropping those numbers straight into the scene puts the hole
kilometres away from the section. This module bridges the two: it
projects each collar into the survey's :class:`~pycsamt.map.geometry.
SurveyFrame`, walks the desurveyed trajectory, and returns per-hole
polylines and coloured interval segments already in scene coordinates,
plus a note on how the hole sits relative to the imaged volume.

Pure numpy + the PCBH render model — no Plotly / Dash / Qt.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..format.borehole import (
    DisplayRadiusPolicy,
    PCBHDocument,
    build_render_model,
)
from .geometry import SurveyFrame

__all__ = [
    "AlignedSegment",
    "AlignedHole",
    "SceneAlignment",
    "surface_from_sections",
    "align_boreholes_to_scene",
]

_ScenePoint = tuple[float, float, float]


@dataclass(frozen=True, repr=False)
class AlignedSegment:
    """One logged interval as a scene-space polyline."""

    borehole_id: str
    family: str
    from_md: float
    to_md: float
    color: str
    points: tuple[_ScenePoint, ...]
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, repr=False)
class AlignedHole:
    """One borehole placed in the 3-D map scene."""

    borehole_id: str
    name: str
    collar_scene: _ScenePoint
    centerline: tuple[_ScenePoint, ...]
    segments: tuple[AlignedSegment, ...]
    display_radius: float
    relation: str
    warnings: tuple[str, ...] = ()
    placed: bool = True


@dataclass(frozen=True, repr=False)
class SceneAlignment:
    """Result of aligning a whole document to one scene."""

    holes: tuple[AlignedHole, ...]
    family: str
    warnings: tuple[str, ...] = ()

    @property
    def placed(self) -> tuple[AlignedHole, ...]:
        return tuple(hole for hole in self.holes if hole.placed)


def surface_from_sections(
    sections: Sequence[tuple[Sequence[float], Sequence[float]]],
) -> Callable[[float], float] | None:
    """Build a ``u -> elevation`` interpolator from section terrain lines.

    ``sections`` is a list of ``(u_values, elevation_values)`` pairs, one
    per rendered line, in the same local metres as
    :class:`~pycsamt.map.geometry.SurveyFrame`. Returns ``None`` when no
    finite terrain is available (the caller then falls back to a flat
    ``z = 0`` datum).
    """
    us: list[float] = []
    zs: list[float] = []
    for u_values, z_values in sections:
        for u, z in zip(u_values, z_values):
            if np.isfinite(u) and np.isfinite(z):
                us.append(float(u))
                zs.append(float(z))
    if len(us) < 2:
        return None
    order = np.argsort(us)
    grid_u = np.asarray(us, dtype=float)[order]
    grid_z = np.asarray(zs, dtype=float)[order]
    grid_u, unique = np.unique(grid_u, return_index=True)
    grid_z = grid_z[unique]

    def surface(u: float) -> float:
        return float(np.interp(u, grid_u, grid_z, grid_z[0], grid_z[-1]))

    return surface


def align_boreholes_to_scene(
    document: PCBHDocument,
    frame: SurveyFrame,
    *,
    family: str = "lithology",
    azimuth_deg: float = 0.0,
    surface: Callable[[float], float] | None = None,
    datum: str = "surface",
    scene_bounds: (
        tuple[float, float, float, float, float, float] | None
    ) = None,
    depth_range: tuple[float, float] | None = None,
    radius_policy: DisplayRadiusPolicy | None = None,
    offset_shift: float = 0.0,
    offset_scale: float = 1.0,
    lean_deg: float = 0.0,
    lean_azimuth_deg: float = 0.0,
    lean_only_vertical: bool = True,
    sampling_step_md: float = 10.0,
    selected_ids: set[str] | None = None,
) -> SceneAlignment:
    """Align every borehole in *document* to a 3-D map scene.

    Parameters
    ----------
    document : PCBHDocument
    frame : SurveyFrame
        The survey's local frame (see :func:`pycsamt.map.geometry.
        survey_frame`), fitted from the same stations the scene uses.
    family : str, default='lithology'
        Interval-log family to colour the trajectory by.
    azimuth_deg : float, default=0
        Scene azimuth, matching ``VolumeMapOptions.azimuth`` — the hole's
        cross-strike offset is rotated by it exactly as line panels are.
    surface : callable, optional
        ``u -> elevation`` terrain function (see
        :func:`surface_from_sections`). Required for ``datum='surface'``.
    datum : {'surface', 'collar_z', 'zero'}, default='surface'
        Depth reference. ``surface`` drapes the collar onto the terrain
        line, ``collar_z`` trusts the PCBH collar elevation, ``zero`` puts
        the collar at ``z = 0`` (matching a no-topography scene).
    scene_bounds : tuple, optional
        ``(umin, umax, vmin, vmax, zmin, zmax)`` of the imaged volume, for
        the per-hole ``relation`` classification.
    depth_range : tuple, optional
        ``(lo, hi)`` metres-below-datum clip applied to the trajectory,
        matching the scene's depth filter.
    radius_policy : DisplayRadiusPolicy, optional
        View-only tube radius policy (physical diameter is never changed).
    offset_shift : float, default=0
        Cross-strike origin shift, in local metres, applied to every
        borehole ``v`` before ``offset_scale``. The
        :mod:`pycsamt.map.volume` builders normalise the line panels so
        the front-most line sits at ``v = 0`` (see
        :func:`pycsamt.map.geometry.normalize_offsets`); pass that same
        shift here or the holes land in front of / behind their line.
    offset_scale : float, default=1
        Cross-strike exaggeration matching ``VolumeMapOptions.line_spacing``
        — the scene stretches the line-to-line axis by this factor, so the
        hole's offset must stretch with it.
    lean_deg : float, default=0
        Apparent lean from vertical, in degrees, applied to holes that
        have no surveyed deviation (see ``lean_only_vertical``). ``0``
        keeps every hole plumb.
    lean_azimuth_deg : float, default=0
        Compass bearing (0 = grid north, 90 = east) the leaned hole tips
        toward.
    lean_only_vertical : bool, default=True
        When true, only holes whose desurveyed trajectory stays within
        ~1 m of the collar are leaned; real deviation surveys pass
        through untouched. Set false to force every hole to ``lean_deg``.
    sampling_step_md : float, default=10
        Trajectory sampling step in measured depth.
    selected_ids : set of str, optional
        Highlighted borehole ids (passed through to the render model).

    Returns
    -------
    SceneAlignment
    """
    if not isinstance(document, PCBHDocument):
        raise TypeError("document must be a PCBHDocument")
    if not isinstance(frame, SurveyFrame):
        raise TypeError("frame must be a SurveyFrame")
    if datum not in {"surface", "collar_z", "zero"}:
        raise ValueError("datum must be 'surface', 'collar_z', or 'zero'")

    model = build_render_model(
        document,
        family=family,
        selected_ids=selected_ids,
        radius_policy=radius_policy or DisplayRadiusPolicy(),
        sampling_step_md=sampling_step_md,
    )
    az = math.radians(float(azimuth_deg))
    sin_az, cos_az = math.sin(az), math.cos(az)
    lean = math.radians(max(0.0, float(lean_deg)))
    bearing = math.radians(float(lean_azimuth_deg))
    # scene (u, v) components of a unit horizontal step along the bearing
    dir_e, dir_n = math.sin(bearing), math.cos(bearing)
    lean_du = math.sin(lean) * (
        dir_e * frame.strike[0] + dir_n * frame.strike[1]
    )
    lean_dv = math.sin(lean) * (
        dir_e * frame.perp[0] + dir_n * frame.perp[1]
    )
    lean_dz = math.cos(lean)
    lonlat = _collar_lonlat(document)
    collars = {hole.id: hole.collar for hole in document.boreholes}

    holes: list[AlignedHole] = []
    doc_warnings: list[str] = []
    for render_hole in model.boreholes:
        hid = render_hole.borehole_id
        collar = collars[hid]
        latlon = lonlat.get(hid)
        warnings: list[str] = []
        if latlon is None:
            holes.append(
                AlignedHole(
                    borehole_id=hid,
                    name=render_hole.collar.label,
                    collar_scene=(math.nan, math.nan, math.nan),
                    centerline=(),
                    segments=(),
                    display_radius=render_hole.display_radius,
                    relation="unplaced",
                    warnings=(
                        "no WGS84 collar and CRS could not be reprojected; "
                        "set collar longitude/latitude or a projected CRS",
                    ),
                    placed=False,
                )
            )
            continue
        lat, lon = latlon
        cu, cv = frame.project(lat, lon)

        if datum == "surface" and surface is not None:
            z_top = surface(cu)
        elif datum == "collar_z":
            z_top = float(collar.z)
        else:
            z_top = 0.0
            if datum == "surface":
                warnings.append("no terrain available; collar placed at z=0")

        lean_active = False
        if lean > 1e-9:
            if lean_only_vertical:
                max_h = max(
                    (
                        math.hypot(
                            float(p.x) - float(collar.x),
                            float(p.y) - float(collar.y),
                        )
                        for p in render_hole.centerline.points
                    ),
                    default=0.0,
                )
                lean_active = max_h <= 1.0
            else:
                lean_active = True
            if lean_active:
                warnings.append(
                    f"shown leaning {float(lean_deg):.0f}° from vertical "
                    f"toward {float(lean_azimuth_deg):.0f}° (display only)"
                )

        place = _Placement(
            frame=frame,
            collar_east=float(collar.x),
            collar_north=float(collar.y),
            cu=cu,
            cv=cv,
            z_top=z_top,
            sin_az=sin_az,
            cos_az=cos_az,
            offset_shift=float(offset_shift),
            offset_scale=float(offset_scale),
            lean_active=lean_active,
            lean_du=lean_du,
            lean_dv=lean_dv,
            lean_dz=lean_dz,
        )
        centerline = tuple(
            place.to_scene(p)
            for p in render_hole.centerline.points
            if _within_depth(place.effective_tvd(p), depth_range)
        )
        segments = tuple(
            AlignedSegment(
                borehole_id=hid,
                family=family,
                from_md=seg.from_md,
                to_md=seg.to_md,
                color=seg.color,
                points=tuple(
                    place.to_scene(p)
                    for p in seg.points
                    if _within_depth(place.effective_tvd(p), depth_range)
                ),
                metadata=dict(seg.metadata),
            )
            for seg in render_hole.interval_segments
        )
        segments = tuple(seg for seg in segments if len(seg.points) >= 2)

        v_scene = (cv - float(offset_shift)) * float(offset_scale)
        collar_scene = (cu + v_scene * sin_az, v_scene * cos_az, z_top)
        relation = _relation(collar_scene, centerline, scene_bounds)
        if relation == "outside":
            warnings.append(
                "collar and trajectory fall outside the imaged volume"
            )
        holes.append(
            AlignedHole(
                borehole_id=hid,
                name=render_hole.collar.label,
                collar_scene=collar_scene,
                centerline=centerline,
                segments=segments,
                display_radius=render_hole.display_radius,
                relation=relation,
                warnings=tuple(warnings),
            )
        )

    return SceneAlignment(
        holes=tuple(holes), family=family, warnings=tuple(doc_warnings)
    )


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _Placement:
    """Per-hole scene transform for one desurveyed trajectory point."""

    frame: SurveyFrame
    collar_east: float
    collar_north: float
    cu: float
    cv: float
    z_top: float
    sin_az: float
    cos_az: float
    offset_shift: float = 0.0
    offset_scale: float = 1.0
    lean_active: bool = False
    lean_du: float = 0.0
    lean_dv: float = 0.0
    lean_dz: float = 1.0

    def effective_tvd(self, point: Any) -> float:
        """True-vertical depth for the depth-window filter.

        A leaned hole's vertical reach is ``md * cos(lean)``, shorter than
        the plumb-hole ``tvd`` the render model carries.
        """
        if self.lean_active:
            return float(point.md) * self.lean_dz
        return float(point.tvd)

    def to_scene(self, point: Any) -> _ScenePoint:
        if self.lean_active:
            md = float(point.md)
            u = self.cu + md * self.lean_du
            v_raw = self.cv + md * self.lean_dv
            dz = md * self.lean_dz
        else:
            d_east = float(point.x) - self.collar_east
            d_north = float(point.y) - self.collar_north
            du = (
                d_east * self.frame.strike[0]
                + d_north * self.frame.strike[1]
            )
            dv = d_east * self.frame.perp[0] + d_north * self.frame.perp[1]
            u = self.cu + du
            v_raw = self.cv + dv
            dz = float(point.tvd)
        v = (v_raw - self.offset_shift) * self.offset_scale
        return (
            u + v * self.sin_az,
            v * self.cos_az,
            self.z_top - dz,
        )


def _within_depth(tvd: float, depth_range: tuple[float, float] | None) -> bool:
    if depth_range is None:
        return True
    lo, hi = depth_range
    return lo - 1e-6 <= float(tvd) <= hi + 1e-6


def _relation(
    collar_scene: _ScenePoint,
    centerline: Sequence[_ScenePoint],
    scene_bounds: tuple[float, float, float, float, float, float] | None,
) -> str:
    if scene_bounds is None:
        return "unknown"
    umin, umax, vmin, vmax, zmin, zmax = scene_bounds
    points = [collar_scene, *centerline]
    inside = 0
    for x, y, z in points:
        if not (
            np.isfinite(x) and np.isfinite(y) and np.isfinite(z)
        ):
            continue
        if umin <= x <= umax and vmin <= y <= vmax and zmin <= z <= zmax:
            inside += 1
    if inside == 0:
        return "outside"
    if inside == len([p for p in points if all(np.isfinite(c) for c in p)]):
        return "inside"
    return "edge"


def _collar_lonlat(
    document: PCBHDocument,
) -> dict[str, tuple[float, float]]:
    """Return ``{borehole_id: (lat, lon)}`` for every locatable hole."""
    result: dict[str, tuple[float, float]] = {}
    transformer = None
    need_transform = False
    for hole in document.boreholes:
        collar = hole.collar
        if collar.latitude is not None and collar.longitude is not None:
            result[hole.id] = (
                float(collar.latitude),
                float(collar.longitude),
            )
        else:
            need_transform = True

    if not need_transform:
        return result

    horizontal = str(document.crs.horizontal or "")
    if horizontal.upper().startswith("LOCAL:") or not horizontal:
        return result
    try:
        from pyproj import Transformer

        transformer = Transformer.from_crs(
            horizontal, "EPSG:4326", always_xy=True
        )
    except Exception:  # noqa: BLE001 - pyproj missing / unknown CRS
        return result

    for hole in document.boreholes:
        if hole.id in result:
            continue
        collar = hole.collar
        try:
            lon, lat = transformer.transform(
                float(collar.x), float(collar.y)
            )
        except Exception:  # noqa: BLE001
            continue
        if np.isfinite(lon) and np.isfinite(lat):
            result[hole.id] = (float(lat), float(lon))
    return result
