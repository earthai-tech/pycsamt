# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""In-memory schema and semantic validation for PCBH 0.1.

PCBH represents one or more boreholes independently of an inversion model.
This module contains no JSON, CSV, LAS, trajectory-desurvey, CRS-transform,
PCSF, or rendering code. Those responsibilities are intentionally deferred
to later implementation phases.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import datetime
from math import isfinite
from typing import Any

from ...api.property import MetadataMixin, PyCSAMTObject

__all__ = [
    "PCBH_VERSION",
    "SEVERITIES",
    "BOREHOLE_KINDS",
    "BOREHOLE_STATUSES",
    "DATA_NATURES",
    "NORTH_REFERENCES",
    "AZIMUTH_DIRECTIONS",
    "ValidationIssue",
    "PCBHValidationError",
    "CoordinateReferenceSystem",
    "UnitSystem",
    "Collar",
    "VocabularyEntry",
    "SurveyStation",
    "Trajectory",
    "LogInterval",
    "StructureObservation",
    "PCBHBorehole",
    "PCBHDocument",
]

PCBH_VERSION = "0.1.0"
SEVERITIES = ("error", "warning", "info")
BOREHOLE_KINDS = (
    "water",
    "mining_exploration",
    "mining_production",
    "geotechnical",
    "environmental",
    "petroleum",
    "geothermal",
    "scientific",
    "monitoring",
    "unknown",
)
BOREHOLE_STATUSES = (
    "planned",
    "drilling",
    "completed",
    "suspended",
    "abandoned",
    "decommissioned",
    "unknown",
)
DATA_NATURES = ("observed", "interpreted", "derived", "unknown")
NORTH_REFERENCES = ("true", "grid", "magnetic", "unknown")
AZIMUTH_DIRECTIONS = ("clockwise",)
TRAJECTORY_METHODS = ("vertical", "survey")
DESURVEY_METHODS = ("minimum_curvature",)
ORIENTATION_REPRESENTATIONS = (
    "none",
    "global_plane",
    "global_line",
    "core_alpha_beta",
)


def _join_path(base: str, field_name: str) -> str:
    return f"{base}.{field_name}" if base else field_name


def _is_namespaced(value: str) -> bool:
    left, sep, right = value.partition(":")
    return bool(sep and left and right)


def _is_allowed_or_namespaced(value: str, allowed: tuple[str, ...]) -> bool:
    return value in allowed or _is_namespaced(value)


def _finite_number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and isfinite(float(value))
    )


@dataclass(frozen=True, repr=False)
class ValidationIssue(PyCSAMTObject):
    """One machine-readable PCBH validation finding."""

    severity: str
    code: str
    message: str
    path: str = ""
    source: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.severity not in SEVERITIES:
            raise ValueError(f"severity must be one of {SEVERITIES}")
        if not str(self.code).strip():
            raise ValueError("validation issue code must be non-empty")
        if not str(self.message).strip():
            raise ValueError("validation issue message must be non-empty")


class PCBHValidationError(ValueError):
    """Raised when error-severity issues remain after validation."""

    def __init__(self, issues: Iterable[ValidationIssue]):
        self.issues = tuple(
            issue for issue in issues if issue.severity == "error"
        )
        preview = "; ".join(
            f"{issue.path or '$'}: {issue.message}"
            for issue in self.issues[:5]
        )
        if len(self.issues) > 5:
            preview += f"; ... ({len(self.issues) - 5} more)"
        super().__init__(preview or "PCBH validation failed")


class _Validatable:
    def collect_issues(self, path: str = "") -> list[ValidationIssue]:
        raise NotImplementedError

    def validate(self) -> None:
        issues = self.collect_issues()
        if any(issue.severity == "error" for issue in issues):
            raise PCBHValidationError(issues)


def _issue(
    code: str, message: str, path: str, severity: str = "error"
) -> ValidationIssue:
    return ValidationIssue(
        severity=severity, code=code, message=message, path=path
    )


def _required_text(value: Any, path: str, code: str) -> list[ValidationIssue]:
    if not isinstance(value, str) or not value.strip():
        return [_issue(code, "must be a non-empty string", path)]
    return []


def _optional_finite(
    value: Any, path: str, code: str
) -> list[ValidationIssue]:
    if value is not None and not _finite_number(value):
        return [_issue(code, "must be a finite number or null", path)]
    return []


@dataclass(repr=False)
class CoordinateReferenceSystem(_Validatable, PyCSAMTObject):
    """Document-wide horizontal and vertical coordinate references."""

    horizontal: str
    vertical: str = "unknown"
    axis_order: str = "xy"
    coordinate_unit: str = "m"

    def collect_issues(self, path: str = "crs") -> list[ValidationIssue]:
        issues = _required_text(
            self.horizontal, _join_path(path, "horizontal"), "crs.horizontal"
        )
        if self.axis_order != "xy":
            issues.append(
                _issue(
                    "crs.axis_order",
                    "must be 'xy' in PCBH 0.1",
                    _join_path(path, "axis_order"),
                )
            )
        issues += _required_text(
            self.vertical, _join_path(path, "vertical"), "crs.vertical"
        )
        issues += _required_text(
            self.coordinate_unit,
            _join_path(path, "coordinate_unit"),
            "crs.coordinate_unit",
        )
        return issues


@dataclass(repr=False)
class UnitSystem(_Validatable, PyCSAMTObject):
    """Units shared by all boreholes in one PCBH 0.1 document."""

    depth: str = "m"
    diameter: str = "m"
    angle: str = "deg"
    resistivity: str = "ohm.m"

    def collect_issues(self, path: str = "units") -> list[ValidationIssue]:
        issues: list[ValidationIssue] = []
        for name in ("depth", "diameter", "angle", "resistivity"):
            issues += _required_text(
                getattr(self, name), _join_path(path, name), f"units.{name}"
            )
        return issues


@dataclass(repr=False)
class Collar(_Validatable, PyCSAMTObject):
    """Borehole reference/origin position in the document CRS."""

    x: float
    y: float
    z: float
    longitude: float | None = None
    latitude: float | None = None
    position_uncertainty_m: float | None = None

    def collect_issues(self, path: str = "collar") -> list[ValidationIssue]:
        issues: list[ValidationIssue] = []
        for name in ("x", "y", "z"):
            if not _finite_number(getattr(self, name)):
                issues.append(
                    _issue(
                        f"collar.{name}",
                        "must be a finite number",
                        _join_path(path, name),
                    )
                )
        issues += _optional_finite(
            self.longitude, _join_path(path, "longitude"), "collar.longitude"
        )
        issues += _optional_finite(
            self.latitude, _join_path(path, "latitude"), "collar.latitude"
        )
        if (
            self.longitude is not None
            and _finite_number(self.longitude)
            and not -180 <= float(self.longitude) <= 180
        ):
            issues.append(
                _issue(
                    "collar.longitude_range",
                    "must be within [-180, 180]",
                    _join_path(path, "longitude"),
                )
            )
        if (
            self.latitude is not None
            and _finite_number(self.latitude)
            and not -90 <= float(self.latitude) <= 90
        ):
            issues.append(
                _issue(
                    "collar.latitude_range",
                    "must be within [-90, 90]",
                    _join_path(path, "latitude"),
                )
            )
        if (self.longitude is None) != (self.latitude is None):
            issues.append(
                _issue(
                    "collar.lonlat_pair",
                    "longitude and latitude must be set together",
                    path,
                )
            )
        issues += _optional_finite(
            self.position_uncertainty_m,
            _join_path(path, "position_uncertainty_m"),
            "collar.uncertainty",
        )
        if (
            _finite_number(self.position_uncertainty_m)
            and float(self.position_uncertainty_m) < 0
        ):
            issues.append(
                _issue(
                    "collar.uncertainty_range",
                    "must be >= 0",
                    _join_path(path, "position_uncertainty_m"),
                )
            )
        return issues


@dataclass(repr=False)
class VocabularyEntry(_Validatable, PyCSAMTObject):
    """Self-contained geological or domain vocabulary entry."""

    code: str
    name: str
    color: str | None = None
    description: str = ""
    external_ids: dict[str, str] = field(default_factory=dict)
    properties: dict[str, Any] = field(default_factory=dict)

    def collect_issues(
        self, path: str = "vocabulary"
    ) -> list[ValidationIssue]:
        issues = _required_text(
            self.code, _join_path(path, "code"), "vocabulary.code"
        )
        issues += _required_text(
            self.name, _join_path(path, "name"), "vocabulary.name"
        )
        if self.color is not None:
            color = str(self.color)
            if (
                len(color) != 7
                or not color.startswith("#")
                or any(ch not in "0123456789abcdefABCDEF" for ch in color[1:])
            ):
                issues.append(
                    _issue(
                        "vocabulary.color",
                        "must be a #RRGGBB color",
                        _join_path(path, "color"),
                    )
                )
        if not isinstance(self.external_ids, dict) or not all(
            isinstance(k, str) and isinstance(v, str)
            for k, v in self.external_ids.items()
        ):
            issues.append(
                _issue(
                    "vocabulary.external_ids",
                    "must map strings to strings",
                    _join_path(path, "external_ids"),
                )
            )
        return issues


@dataclass(repr=False)
class SurveyStation(_Validatable, PyCSAMTObject):
    """One measured-depth orientation station."""

    md: float
    azimuth_deg: float
    inclination_deg: float

    def collect_issues(
        self, path: str = "survey_station"
    ) -> list[ValidationIssue]:
        issues: list[ValidationIssue] = []
        for name in ("md", "azimuth_deg", "inclination_deg"):
            if not _finite_number(getattr(self, name)):
                issues.append(
                    _issue(
                        f"trajectory.{name}",
                        "must be a finite number",
                        _join_path(path, name),
                    )
                )
        if _finite_number(self.md) and float(self.md) < 0:
            issues.append(
                _issue(
                    "trajectory.md_range",
                    "must be >= 0",
                    _join_path(path, "md"),
                )
            )
        if (
            _finite_number(self.azimuth_deg)
            and not 0 <= float(self.azimuth_deg) < 360
        ):
            issues.append(
                _issue(
                    "trajectory.azimuth_range",
                    "must be within [0, 360)",
                    _join_path(path, "azimuth_deg"),
                )
            )
        if (
            _finite_number(self.inclination_deg)
            and not 0 <= float(self.inclination_deg) <= 180
        ):
            issues.append(
                _issue(
                    "trajectory.inclination_range",
                    "must be within [0, 180]",
                    _join_path(path, "inclination_deg"),
                )
            )
        if (
            _finite_number(self.inclination_deg)
            and float(self.inclination_deg) > 90
        ):
            issues.append(
                _issue(
                    "trajectory.reentry",
                    "inclination above 90 degrees represents an "
                    "upward/re-entry segment",
                    _join_path(path, "inclination_deg"),
                    "warning",
                )
            )
        return issues


@dataclass(repr=False)
class Trajectory(_Validatable, PyCSAMTObject):
    """Vertical shorthand or measured survey for a borehole path."""

    method: str = "vertical"
    north_reference: str = "unknown"
    desurvey_method: str = "minimum_curvature"
    stations: list[SurveyStation] = field(default_factory=list)

    def collect_issues(
        self, path: str = "trajectory"
    ) -> list[ValidationIssue]:
        issues: list[ValidationIssue] = []
        if self.method not in TRAJECTORY_METHODS:
            issues.append(
                _issue(
                    "trajectory.method",
                    f"must be one of {TRAJECTORY_METHODS}",
                    _join_path(path, "method"),
                )
            )
        if self.north_reference not in NORTH_REFERENCES:
            issues.append(
                _issue(
                    "trajectory.north_reference",
                    f"must be one of {NORTH_REFERENCES}",
                    _join_path(path, "north_reference"),
                )
            )
        if self.desurvey_method not in DESURVEY_METHODS:
            issues.append(
                _issue(
                    "trajectory.desurvey_method",
                    f"must be one of {DESURVEY_METHODS}",
                    _join_path(path, "desurvey_method"),
                )
            )
        if self.method == "vertical" and self.stations:
            issues.append(
                _issue(
                    "trajectory.vertical_stations",
                    "vertical trajectory must not contain survey stations",
                    _join_path(path, "stations"),
                )
            )
        if self.method == "survey" and len(self.stations) < 2:
            issues.append(
                _issue(
                    "trajectory.station_count",
                    "survey trajectory needs at least two stations",
                    _join_path(path, "stations"),
                )
            )
        previous: float | None = None
        for index, station in enumerate(self.stations):
            station_path = f"{path}.stations[{index}]"
            if not isinstance(station, SurveyStation):
                issues.append(
                    _issue(
                        "trajectory.station_type",
                        "must be a SurveyStation",
                        station_path,
                    )
                )
                continue
            issues.extend(station.collect_issues(station_path))
            if _finite_number(station.md):
                md = float(station.md)
                if previous is not None and md <= previous:
                    issues.append(
                        _issue(
                            "trajectory.md_order",
                            "station measured depths must increase strictly",
                            _join_path(station_path, "md"),
                        )
                    )
                previous = md
        if (
            self.method == "survey"
            and self.stations
            and _finite_number(self.stations[0].md)
            and float(self.stations[0].md) != 0
        ):
            issues.append(
                _issue(
                    "trajectory.missing_collar_station",
                    "first survey station should be at md 0",
                    f"{path}.stations[0].md",
                    "warning",
                )
            )
        if self.north_reference == "magnetic":
            issues.append(
                _issue(
                    "trajectory.magnetic_provenance",
                    "magnetic azimuth requires declination/date/source "
                    "provenance before conversion",
                    _join_path(path, "north_reference"),
                    "warning",
                )
            )
        return issues


@dataclass(repr=False)
class LogInterval(_Validatable, PyCSAMTObject):
    """One measured-depth interval in a named log family."""

    from_md: float
    to_md: float
    code: str | None = None
    label: str | None = None
    description: str = ""
    resistivity_ohm_m: float | None = None
    data_nature: str = "unknown"
    confidence: float | None = None
    properties: dict[str, Any] = field(default_factory=dict)

    def collect_issues(self, path: str = "interval") -> list[ValidationIssue]:
        issues: list[ValidationIssue] = []
        for name in ("from_md", "to_md"):
            if not _finite_number(getattr(self, name)):
                issues.append(
                    _issue(
                        f"interval.{name}",
                        "must be a finite number",
                        _join_path(path, name),
                    )
                )
        if _finite_number(self.from_md) and float(self.from_md) < 0:
            issues.append(
                _issue(
                    "interval.from_range",
                    "must be >= 0",
                    _join_path(path, "from_md"),
                )
            )
        if (
            _finite_number(self.from_md)
            and _finite_number(self.to_md)
            and float(self.to_md) <= float(self.from_md)
        ):
            issues.append(
                _issue(
                    "interval.order",
                    "to_md must be greater than from_md",
                    path,
                )
            )
        if not (isinstance(self.code, str) and self.code.strip()) and not (
            isinstance(self.label, str) and self.label.strip()
        ):
            issues.append(
                _issue(
                    "interval.identity",
                    "at least one of code or label must be non-empty",
                    path,
                )
            )
        if self.data_nature not in DATA_NATURES:
            issues.append(
                _issue(
                    "interval.data_nature",
                    f"must be one of {DATA_NATURES}",
                    _join_path(path, "data_nature"),
                )
            )
        issues += _optional_finite(
            self.resistivity_ohm_m,
            _join_path(path, "resistivity_ohm_m"),
            "interval.resistivity",
        )
        if (
            _finite_number(self.resistivity_ohm_m)
            and float(self.resistivity_ohm_m) <= 0
        ):
            issues.append(
                _issue(
                    "interval.resistivity_range",
                    "must be > 0",
                    _join_path(path, "resistivity_ohm_m"),
                )
            )
        issues += _optional_finite(
            self.confidence,
            _join_path(path, "confidence"),
            "interval.confidence",
        )
        if (
            _finite_number(self.confidence)
            and not 0 <= float(self.confidence) <= 1
        ):
            issues.append(
                _issue(
                    "interval.confidence_range",
                    "must be within [0, 1]",
                    _join_path(path, "confidence"),
                )
            )
        return issues


@dataclass(repr=False)
class StructureObservation(_Validatable, PyCSAMTObject):
    """Point or interval structural observation along a borehole."""

    kind: str
    at_md: float | None = None
    from_md: float | None = None
    to_md: float | None = None
    orientation_representation: str = "none"
    strike_deg: float | None = None
    dip_deg: float | None = None
    dip_direction_deg: float | None = None
    trend_deg: float | None = None
    plunge_deg: float | None = None
    alpha_deg: float | None = None
    beta_deg: float | None = None
    aperture_m: float | None = None
    fill: str | None = None
    data_nature: str = "unknown"
    confidence: float | None = None

    def collect_issues(self, path: str = "structure") -> list[ValidationIssue]:
        issues = _required_text(
            self.kind, _join_path(path, "kind"), "structure.kind"
        )
        point = self.at_md is not None
        interval = self.from_md is not None or self.to_md is not None
        if point == interval:
            issues.append(
                _issue(
                    "structure.location",
                    "set exactly one of at_md or a from_md/to_md interval",
                    path,
                )
            )
        if point:
            issues += _optional_finite(
                self.at_md, _join_path(path, "at_md"), "structure.at_md"
            )
            if _finite_number(self.at_md) and float(self.at_md) < 0:
                issues.append(
                    _issue(
                        "structure.at_md_range",
                        "must be >= 0",
                        _join_path(path, "at_md"),
                    )
                )
        if interval:
            if self.from_md is None or self.to_md is None:
                issues.append(
                    _issue(
                        "structure.interval_pair",
                        "from_md and to_md must be set together",
                        path,
                    )
                )
            elif (
                not _finite_number(self.from_md)
                or not _finite_number(self.to_md)
                or float(self.from_md) < 0
                or float(self.to_md) <= float(self.from_md)
            ):
                issues.append(
                    _issue(
                        "structure.interval_range",
                        "requires finite 0 <= from_md < to_md",
                        path,
                    )
                )
        if self.orientation_representation not in ORIENTATION_REPRESENTATIONS:
            issues.append(
                _issue(
                    "structure.orientation_representation",
                    f"must be one of {ORIENTATION_REPRESENTATIONS}",
                    _join_path(path, "orientation_representation"),
                )
            )
        required_by_representation = {
            "global_plane": ("dip_deg", "dip_direction_deg"),
            "global_line": ("trend_deg", "plunge_deg"),
            "core_alpha_beta": ("alpha_deg", "beta_deg"),
        }
        for name in required_by_representation.get(
            self.orientation_representation, ()
        ):
            if not _finite_number(getattr(self, name)):
                issues.append(
                    _issue(
                        "structure.orientation_required",
                        f"{name} is required for "
                        f"{self.orientation_representation}",
                        _join_path(path, name),
                    )
                )
        for name in (
            "strike_deg",
            "dip_direction_deg",
            "trend_deg",
            "alpha_deg",
            "beta_deg",
        ):
            value = getattr(self, name)
            if value is not None and (
                not _finite_number(value) or not 0 <= float(value) < 360
            ):
                issues.append(
                    _issue(
                        "structure.bearing_range",
                        "must be within [0, 360)",
                        _join_path(path, name),
                    )
                )
        for name in ("dip_deg", "plunge_deg"):
            value = getattr(self, name)
            if value is not None and (
                not _finite_number(value) or not 0 <= float(value) <= 90
            ):
                issues.append(
                    _issue(
                        "structure.angle_range",
                        "must be within [0, 90]",
                        _join_path(path, name),
                    )
                )
        issues += _optional_finite(
            self.aperture_m,
            _join_path(path, "aperture_m"),
            "structure.aperture",
        )
        if _finite_number(self.aperture_m) and float(self.aperture_m) < 0:
            issues.append(
                _issue(
                    "structure.aperture_range",
                    "must be >= 0",
                    _join_path(path, "aperture_m"),
                )
            )
        if self.data_nature not in DATA_NATURES:
            issues.append(
                _issue(
                    "structure.data_nature",
                    f"must be one of {DATA_NATURES}",
                    _join_path(path, "data_nature"),
                )
            )
        issues += _optional_finite(
            self.confidence,
            _join_path(path, "confidence"),
            "structure.confidence",
        )
        if (
            _finite_number(self.confidence)
            and not 0 <= float(self.confidence) <= 1
        ):
            issues.append(
                _issue(
                    "structure.confidence_range",
                    "must be within [0, 1]",
                    _join_path(path, "confidence"),
                )
            )
        return issues


@dataclass(repr=False)
class PCBHBorehole(_Validatable, PyCSAMTObject, MetadataMixin):
    """One borehole and its trajectory, logs, and structures."""

    id: str
    name: str
    kind: str
    status: str
    collar: Collar
    total_depth_md: float
    trajectory: Trajectory = field(default_factory=Trajectory)
    diameter: float | None = None
    interval_logs: dict[str, list[LogInterval]] = field(default_factory=dict)
    structures: list[StructureObservation] = field(default_factory=list)
    aliases: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    extensions: dict[str, Any] = field(default_factory=dict)

    def collect_issues(self, path: str = "borehole") -> list[ValidationIssue]:
        issues = _required_text(self.id, _join_path(path, "id"), "borehole.id")
        issues += _required_text(
            self.name, _join_path(path, "name"), "borehole.name"
        )
        if not _is_allowed_or_namespaced(self.kind, BOREHOLE_KINDS):
            issues.append(
                _issue(
                    "borehole.kind",
                    "must be a standard or namespaced kind; standard "
                    f"values are {BOREHOLE_KINDS}",
                    _join_path(path, "kind"),
                )
            )
        if not _is_allowed_or_namespaced(self.status, BOREHOLE_STATUSES):
            issues.append(
                _issue(
                    "borehole.status",
                    "must be a standard or namespaced status; standard "
                    f"values are {BOREHOLE_STATUSES}",
                    _join_path(path, "status"),
                )
            )
        if (
            not _finite_number(self.total_depth_md)
            or float(self.total_depth_md) <= 0
        ):
            issues.append(
                _issue(
                    "borehole.total_depth",
                    "must be a finite number > 0",
                    _join_path(path, "total_depth_md"),
                )
            )
        issues += _optional_finite(
            self.diameter, _join_path(path, "diameter"), "borehole.diameter"
        )
        if _finite_number(self.diameter) and float(self.diameter) <= 0:
            issues.append(
                _issue(
                    "borehole.diameter_range",
                    "must be > 0",
                    _join_path(path, "diameter"),
                )
            )
        if not isinstance(self.collar, Collar):
            issues.append(
                _issue(
                    "borehole.collar_type",
                    "must be a Collar",
                    _join_path(path, "collar"),
                )
            )
        else:
            issues.extend(
                self.collar.collect_issues(_join_path(path, "collar"))
            )
        if not isinstance(self.trajectory, Trajectory):
            issues.append(
                _issue(
                    "borehole.trajectory_type",
                    "must be a Trajectory",
                    _join_path(path, "trajectory"),
                )
            )
        else:
            issues.extend(
                self.trajectory.collect_issues(_join_path(path, "trajectory"))
            )
        max_depth = (
            float(self.total_depth_md)
            if _finite_number(self.total_depth_md)
            else None
        )
        if isinstance(self.trajectory, Trajectory):
            for index, station in enumerate(self.trajectory.stations):
                if (
                    max_depth is not None
                    and isinstance(station, SurveyStation)
                    and _finite_number(station.md)
                    and float(station.md) > max_depth
                ):
                    issues.append(
                        _issue(
                            "borehole.station_beyond_td",
                            "survey station exceeds total_depth_md",
                            f"{path}.trajectory.stations[{index}].md",
                        )
                    )
        known_families = {
            "lithology",
            "formation",
            "weathering",
            "alteration",
            "mineralization",
            "oxidation",
            "hydrostratigraphy",
            "geotechnical",
            "interpretation",
        }
        for family, intervals in self.interval_logs.items():
            family_path = f"{path}.interval_logs.{family}"
            if family not in known_families and not _is_namespaced(family):
                issues.append(
                    _issue(
                        "interval.family",
                        "custom log family must be namespaced",
                        family_path,
                    )
                )
            valid_intervals: list[tuple[int, LogInterval]] = []
            for index, interval in enumerate(intervals):
                interval_path = f"{family_path}[{index}]"
                if not isinstance(interval, LogInterval):
                    issues.append(
                        _issue(
                            "interval.type",
                            "must be a LogInterval",
                            interval_path,
                        )
                    )
                    continue
                issues.extend(interval.collect_issues(interval_path))
                valid_intervals.append((index, interval))
                if (
                    max_depth is not None
                    and _finite_number(interval.to_md)
                    and float(interval.to_md) > max_depth
                ):
                    issues.append(
                        _issue(
                            "borehole.interval_beyond_td",
                            "interval exceeds total_depth_md",
                            _join_path(interval_path, "to_md"),
                        )
                    )
                if interval.code and family == "lithology":
                    pass
            ordered = sorted(
                valid_intervals,
                key=lambda item: (
                    float(item[1].from_md)
                    if _finite_number(item[1].from_md)
                    else float("inf")
                ),
            )
            previous_to: float | None = None
            for index, interval in ordered:
                if (
                    previous_to is not None
                    and _finite_number(interval.from_md)
                    and float(interval.from_md) < previous_to
                ):
                    issues.append(
                        _issue(
                            "interval.overlap",
                            "overlaps another interval in exclusive family "
                            f"{family!r}",
                            f"{family_path}[{index}]",
                        )
                    )
                if _finite_number(interval.to_md):
                    previous_to = max(
                        previous_to or float("-inf"), float(interval.to_md)
                    )
        for index, structure in enumerate(self.structures):
            structure_path = f"{path}.structures[{index}]"
            if not isinstance(structure, StructureObservation):
                issues.append(
                    _issue(
                        "structure.type",
                        "must be a StructureObservation",
                        structure_path,
                    )
                )
                continue
            issues.extend(structure.collect_issues(structure_path))
            depth = (
                structure.at_md
                if structure.at_md is not None
                else structure.to_md
            )
            if (
                max_depth is not None
                and _finite_number(depth)
                and float(depth) > max_depth
            ):
                issues.append(
                    _issue(
                        "borehole.structure_beyond_td",
                        "structure exceeds total_depth_md",
                        structure_path,
                    )
                )
        for key in self.extensions:
            if not isinstance(key, str) or not _is_namespaced(key):
                issues.append(
                    _issue(
                        "extensions.namespace",
                        "extension keys must be namespaced",
                        f"{path}.extensions.{key}",
                    )
                )
        return issues


@dataclass(repr=False)
class PCBHDocument(_Validatable, PyCSAMTObject, MetadataMixin):
    """Root PCBH document containing one or more boreholes."""

    document_id: str
    created_at: str
    created_by: str
    crs: CoordinateReferenceSystem
    boreholes: list[PCBHBorehole]
    units: UnitSystem = field(default_factory=UnitSystem)
    title: str = ""
    description: str = ""
    lithologies: list[VocabularyEntry] = field(default_factory=list)
    formations: list[VocabularyEntry] = field(default_factory=list)
    pcbh_version: str = PCBH_VERSION
    azimuth_direction: str = "clockwise"
    inclination_reference: str = "vertical_down"
    depth_reference: str = "collar"
    z_positive: str = "up"
    metadata: dict[str, Any] = field(default_factory=dict)
    extensions: dict[str, Any] = field(default_factory=dict)

    def collect_issues(self, path: str = "$") -> list[ValidationIssue]:
        issues: list[ValidationIssue] = []
        try:
            version_parts = tuple(
                int(part) for part in self.pcbh_version.split(".")
            )
        except (AttributeError, ValueError):
            version_parts = ()
        if len(version_parts) != 3 or version_parts[0] != 0:
            issues.append(
                _issue(
                    "document.version",
                    f"must be compatible with PCBH {PCBH_VERSION!r}",
                    f"{path}.pcbh_version",
                )
            )
        issues += _required_text(
            self.document_id, f"{path}.document_id", "document.id"
        )
        issues += _required_text(
            self.created_at, f"{path}.created_at", "document.created_at"
        )
        if isinstance(self.created_at, str) and self.created_at.strip():
            try:
                datetime.fromisoformat(
                    self.created_at.strip().replace("Z", "+00:00")
                )
            except ValueError:
                issues.append(
                    _issue(
                        "document.created_at_format",
                        "must be an ISO-8601 timestamp",
                        f"{path}.created_at",
                    )
                )
        issues += _required_text(
            self.created_by, f"{path}.created_by", "document.created_by"
        )
        if self.azimuth_direction not in AZIMUTH_DIRECTIONS:
            issues.append(
                _issue(
                    "document.azimuth_direction",
                    "must be 'clockwise'",
                    f"{path}.azimuth_direction",
                )
            )
        if self.inclination_reference != "vertical_down":
            issues.append(
                _issue(
                    "document.inclination_reference",
                    "must be 'vertical_down'",
                    f"{path}.inclination_reference",
                )
            )
        if self.depth_reference != "collar":
            issues.append(
                _issue(
                    "document.depth_reference",
                    "must be 'collar'",
                    f"{path}.depth_reference",
                )
            )
        if self.z_positive != "up":
            issues.append(
                _issue(
                    "document.z_positive", "must be 'up'", f"{path}.z_positive"
                )
            )
        if not isinstance(self.crs, CoordinateReferenceSystem):
            issues.append(
                _issue(
                    "document.crs_type",
                    "must be a CoordinateReferenceSystem",
                    f"{path}.crs",
                )
            )
        else:
            issues.extend(self.crs.collect_issues(f"{path}.crs"))
        if not isinstance(self.units, UnitSystem):
            issues.append(
                _issue(
                    "document.units_type",
                    "must be a UnitSystem",
                    f"{path}.units",
                )
            )
        else:
            issues.extend(self.units.collect_issues(f"{path}.units"))
        if not self.boreholes:
            issues.append(
                _issue(
                    "document.boreholes_empty",
                    "must contain at least one borehole",
                    f"{path}.boreholes",
                )
            )
        seen_ids: set[str] = set()
        seen_folded: dict[str, str] = {}
        for index, borehole in enumerate(self.boreholes):
            borehole_path = f"{path}.boreholes[{index}]"
            if not isinstance(borehole, PCBHBorehole):
                issues.append(
                    _issue(
                        "borehole.type",
                        "must be a PCBHBorehole",
                        borehole_path,
                    )
                )
                continue
            issues.extend(borehole.collect_issues(borehole_path))
            if borehole.id in seen_ids:
                issues.append(
                    _issue(
                        "borehole.duplicate_id",
                        f"duplicate borehole id {borehole.id!r}",
                        _join_path(borehole_path, "id"),
                    )
                )
            seen_ids.add(borehole.id)
            folded = borehole.id.casefold()
            if folded in seen_folded and seen_folded[folded] != borehole.id:
                issues.append(
                    _issue(
                        "borehole.case_collision",
                        "id differs only by case from "
                        f"{seen_folded[folded]!r}",
                        _join_path(borehole_path, "id"),
                        "warning",
                    )
                )
            seen_folded[folded] = borehole.id
        for vocabulary_name, entries in (
            ("lithologies", self.lithologies),
            ("formations", self.formations),
        ):
            codes: set[str] = set()
            for index, entry in enumerate(entries):
                entry_path = f"{path}.{vocabulary_name}[{index}]"
                if not isinstance(entry, VocabularyEntry):
                    issues.append(
                        _issue(
                            "vocabulary.type",
                            "must be a VocabularyEntry",
                            entry_path,
                        )
                    )
                    continue
                issues.extend(entry.collect_issues(entry_path))
                if entry.code in codes:
                    issues.append(
                        _issue(
                            "vocabulary.duplicate_code",
                            f"duplicate code {entry.code!r}",
                            _join_path(entry_path, "code"),
                        )
                    )
                codes.add(entry.code)
        lithology_codes = {
            entry.code
            for entry in self.lithologies
            if isinstance(entry, VocabularyEntry)
        }
        formation_codes = {
            entry.code
            for entry in self.formations
            if isinstance(entry, VocabularyEntry)
        }
        for bh_index, borehole in enumerate(self.boreholes):
            if not isinstance(borehole, PCBHBorehole):
                continue
            for family, known_codes in (
                ("lithology", lithology_codes),
                ("formation", formation_codes),
            ):
                for interval_index, interval in enumerate(
                    borehole.interval_logs.get(family, [])
                ):
                    if (
                        isinstance(interval, LogInterval)
                        and interval.code
                        and interval.code not in known_codes
                    ):
                        issues.append(
                            _issue(
                                "interval.unresolved_code",
                                f"code {interval.code!r} is absent from the "
                                f"{family} vocabulary",
                                f"{path}.boreholes[{bh_index}].interval_logs.{family}[{interval_index}].code",
                            )
                        )
        for key in self.extensions:
            if not isinstance(key, str) or not _is_namespaced(key):
                issues.append(
                    _issue(
                        "extensions.namespace",
                        "extension keys must be namespaced",
                        f"{path}.extensions.{key}",
                    )
                )
        return issues
