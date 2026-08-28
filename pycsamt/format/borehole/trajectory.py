# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Deterministic 3-D trajectory derivation for PCBH boreholes.

Coordinates follow the PCBH contract: X is easting, Y is northing, Z is
elevation positive upward, and TVD is positive downward from the collar.
Survey inclination is measured from vertical down.  The implementation uses
minimum curvature and interpolates attitudes on the unit sphere, which also
makes north-crossing azimuths unambiguous.
"""

from __future__ import annotations

import hashlib
import json
import math
from bisect import bisect_right
from collections.abc import Iterable
from dataclasses import dataclass

from ...api.property import PyCSAMTObject
from .schema import PCBHBorehole, SurveyStation

__all__ = [
    "TrajectoryPoint",
    "DesurveyedTrajectory",
    "desurvey",
    "trajectory_checksum",
]

_ANGLE_EPS = 1e-12
_MD_EPS = 1e-10


@dataclass(frozen=True, repr=False)
class TrajectoryPoint(PyCSAMTObject):
    """One derived point and local borehole attitude."""

    md: float
    x: float
    y: float
    z: float
    tvd: float
    azimuth_deg: float
    inclination_deg: float

    def __post_init__(self) -> None:
        names = (
            "md",
            "x",
            "y",
            "z",
            "tvd",
            "azimuth_deg",
            "inclination_deg",
        )
        for name in names:
            value = getattr(self, name)
            if isinstance(value, bool):
                raise ValueError(
                    "trajectory point values must be finite numbers"
                )
            try:
                number = float(value)
            except (TypeError, ValueError) as error:
                raise ValueError(
                    "trajectory point values must be finite numbers"
                ) from error
            if not math.isfinite(number):
                raise ValueError(
                    "trajectory point values must be finite numbers"
                )
            object.__setattr__(self, name, number)
        if self.md < 0:
            raise ValueError("trajectory point md must be >= 0")
        if not 0 <= self.azimuth_deg < 360:
            raise ValueError(
                "trajectory point azimuth must be within [0, 360)"
            )
        if not 0 <= self.inclination_deg <= 180:
            raise ValueError(
                "trajectory point inclination must be within [0, 180]"
            )


@dataclass(frozen=True, repr=False)
class DesurveyedTrajectory(PyCSAMTObject):
    """Immutable centerline generated from one PCBH borehole."""

    points: tuple[TrajectoryPoint, ...]
    source_checksum: str
    method: str = "minimum_curvature"

    def __post_init__(self) -> None:
        if not self.points:
            raise ValueError(
                "a desurveyed trajectory needs at least one point"
            )
        if any(
            current.md <= previous.md
            for previous, current in zip(self.points, self.points[1:])
        ):
            raise ValueError("trajectory point MDs must increase strictly")
        if self.method != "minimum_curvature":
            raise ValueError("unsupported desurvey method")
        if (
            not isinstance(self.source_checksum, str)
            or not self.source_checksum
        ):
            raise ValueError("source_checksum must be a non-empty string")

    def at_md(self, md: float) -> TrajectoryPoint:
        """Return an exact or minimum-curvature-interpolated point at *md*."""
        target = _finite_md(md)
        first, last = self.points[0], self.points[-1]
        if target < first.md - _MD_EPS or target > last.md + _MD_EPS:
            raise ValueError(
                f"md must be within [{first.md}, {last.md}], got {target}"
            )
        if math.isclose(target, first.md, abs_tol=_MD_EPS):
            return first
        if math.isclose(target, last.md, abs_tol=_MD_EPS):
            return last

        depths = [point.md for point in self.points]
        upper = bisect_right(depths, target)
        start, end = self.points[upper - 1], self.points[upper]
        if math.isclose(target, start.md, abs_tol=_MD_EPS):
            return start
        fraction = (target - start.md) / (end.md - start.md)
        direction = _slerp(
            _direction(start.azimuth_deg, start.inclination_deg),
            _direction(end.azimuth_deg, end.inclination_deg),
            fraction,
        )
        azimuth, inclination = _attitude(direction)
        east, north, down = _minimum_curvature_displacement(
            target - start.md,
            _direction(start.azimuth_deg, start.inclination_deg),
            direction,
        )
        return TrajectoryPoint(
            md=target,
            x=start.x + east,
            y=start.y + north,
            z=start.z - down,
            tvd=start.tvd + down,
            azimuth_deg=azimuth,
            inclination_deg=inclination,
        )

    def split_at(
        self, measured_depths: Iterable[float]
    ) -> DesurveyedTrajectory:
        """Return a centerline containing every requested interval boundary."""
        depths = {point.md for point in self.points}
        for md in measured_depths:
            depths.add(_finite_md(md))
        points = tuple(self.at_md(md) for md in sorted(depths))
        return DesurveyedTrajectory(
            points=points,
            source_checksum=self.source_checksum,
            method=self.method,
        )


def desurvey(
    borehole: PCBHBorehole,
    *,
    boundaries: Iterable[float] = (),
) -> DesurveyedTrajectory:
    """Generate a 3-D centerline for *borehole*.

    Survey stations after total depth are rejected by schema validation.  If
    the first survey station is below the collar, its attitude is extended
    back to MD 0.  Likewise, the final attitude is extended to total depth.
    """
    if not isinstance(borehole, PCBHBorehole):
        raise TypeError("borehole must be a PCBHBorehole")
    borehole.validate()
    total_depth = _finite_md(borehole.total_depth_md)
    if total_depth <= 0:
        raise ValueError("total_depth_md must be greater than zero")

    collar = borehole.collar
    if borehole.trajectory.method == "vertical":
        stations = (
            SurveyStation(0.0, 0.0, 0.0),
            SurveyStation(total_depth, 0.0, 0.0),
        )
    else:
        stations = _stations_with_endpoints(
            borehole.trajectory.stations, total_depth
        )

    first = stations[0]
    points = [
        TrajectoryPoint(
            md=0.0,
            x=float(collar.x),
            y=float(collar.y),
            z=float(collar.z),
            tvd=0.0,
            azimuth_deg=float(first.azimuth_deg) % 360.0,
            inclination_deg=float(first.inclination_deg),
        )
    ]
    for start_station, end_station in zip(stations, stations[1:]):
        start = points[-1]
        end_direction = _direction(
            end_station.azimuth_deg, end_station.inclination_deg
        )
        east, north, down = _minimum_curvature_displacement(
            end_station.md - start_station.md,
            _direction(
                start_station.azimuth_deg, start_station.inclination_deg
            ),
            end_direction,
        )
        azimuth, inclination = _attitude(end_direction)
        points.append(
            TrajectoryPoint(
                md=float(end_station.md),
                x=start.x + east,
                y=start.y + north,
                z=start.z - down,
                tvd=start.tvd + down,
                azimuth_deg=azimuth,
                inclination_deg=inclination,
            )
        )

    result = DesurveyedTrajectory(
        points=tuple(points),
        source_checksum=trajectory_checksum(borehole),
    )
    return result.split_at(boundaries) if boundaries else result


def trajectory_checksum(borehole: PCBHBorehole) -> str:
    """Return a stable SHA-256 key for trajectory-defining inputs."""
    if not isinstance(borehole, PCBHBorehole):
        raise TypeError("borehole must be a PCBHBorehole")
    payload = {
        "collar": [borehole.collar.x, borehole.collar.y, borehole.collar.z],
        "total_depth_md": borehole.total_depth_md,
        "trajectory": {
            "method": borehole.trajectory.method,
            "north_reference": borehole.trajectory.north_reference,
            "desurvey_method": borehole.trajectory.desurvey_method,
            "stations": [
                [station.md, station.azimuth_deg, station.inclination_deg]
                for station in borehole.trajectory.stations
            ],
        },
    }
    encoded = json.dumps(
        payload,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _stations_with_endpoints(
    source: list[SurveyStation], total_depth: float
) -> tuple[SurveyStation, ...]:
    stations = [
        SurveyStation(
            float(station.md),
            float(station.azimuth_deg),
            float(station.inclination_deg),
        )
        for station in source
    ]
    if stations[0].md > 0:
        stations.insert(
            0,
            SurveyStation(
                0.0,
                stations[0].azimuth_deg,
                stations[0].inclination_deg,
            ),
        )
    if stations[-1].md < total_depth:
        stations.append(
            SurveyStation(
                total_depth,
                stations[-1].azimuth_deg,
                stations[-1].inclination_deg,
            )
        )
    return tuple(stations)


def _finite_md(value: float) -> float:
    if isinstance(value, bool):
        raise TypeError("md must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError("md must be a finite number") from error
    if not math.isfinite(result):
        raise ValueError("md must be a finite number")
    return result


def _direction(
    azimuth_deg: float, inclination_deg: float
) -> tuple[float, float, float]:
    azimuth = math.radians(float(azimuth_deg))
    inclination = math.radians(float(inclination_deg))
    horizontal = math.sin(inclination)
    return (
        horizontal * math.sin(azimuth),
        horizontal * math.cos(azimuth),
        math.cos(inclination),
    )


def _attitude(direction: tuple[float, float, float]) -> tuple[float, float]:
    east, north, down = direction
    inclination = math.degrees(math.acos(max(-1.0, min(1.0, down))))
    if math.hypot(east, north) <= _ANGLE_EPS:
        azimuth = 0.0
    else:
        azimuth = math.degrees(math.atan2(east, north)) % 360.0
        if math.isclose(azimuth, 360.0, abs_tol=1e-12):
            azimuth = 0.0
    return azimuth, inclination


def _dot(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
) -> float:
    return sum(a * b for a, b in zip(left, right))


def _slerp(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    fraction: float,
) -> tuple[float, float, float]:
    cosine = max(-1.0, min(1.0, _dot(start, end)))
    angle = math.acos(cosine)
    if angle <= _ANGLE_EPS:
        return start
    sine = math.sin(angle)
    if abs(sine) <= _ANGLE_EPS:
        raise ValueError("a 180-degree survey dogleg has no unique path")
    start_weight = math.sin((1.0 - fraction) * angle) / sine
    end_weight = math.sin(fraction * angle) / sine
    return tuple(start_weight * a + end_weight * b for a, b in zip(start, end))


def _minimum_curvature_displacement(
    delta_md: float,
    start: tuple[float, float, float],
    end: tuple[float, float, float],
) -> tuple[float, float, float]:
    cosine = max(-1.0, min(1.0, _dot(start, end)))
    dogleg = math.acos(cosine)
    if math.pi - dogleg <= _ANGLE_EPS:
        raise ValueError("a 180-degree survey dogleg has no unique path")
    if dogleg <= 1e-6:
        squared = dogleg * dogleg
        ratio_factor = 1.0 + squared / 12.0 + squared * squared / 120.0
    else:
        ratio_factor = 2.0 * math.tan(dogleg / 2.0) / dogleg
    scale = float(delta_md) * ratio_factor / 2.0
    return tuple(scale * (a + b) for a, b in zip(start, end))
