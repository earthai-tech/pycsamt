# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.metadata.survey
=======================

Campaign-level survey descriptor.

:class:`SurveyMeta` captures the essential context of an MT/AMT/CSAMT/TEM
measurement campaign — who, what, when, where — as a first-class Python
object.  It serialises to / deserialises from JSON and YAML, can be built
directly from a :class:`~pycsamt.site.base.Sites` collection, and can push
its fields into an EDI ``>HEAD`` section.

Quick start
-----------
::

    from pycsamt.metadata.survey import SurveyMeta, BBox

    meta = SurveyMeta(
        name="WILLY_2023",
        project="Copper-gold exploration",
        operator="EarthAI-Tech",
        method="AMT",
        bbox=BBox(lat_min=27.8, lat_max=28.9, lon_min=101.5, lon_max=103.2),
    )

    # build from a Sites collection
    meta2 = SurveyMeta.from_sites(sites, name="WILLY_2023")

    # round-trip to JSON
    meta.to_json("survey.json")
    meta3 = SurveyMeta.from_json("survey.json")
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

__all__ = ["BBox", "SurveyMeta"]

_VALID_METHODS = frozenset({"MT", "AMT", "CSAMT", "CSEM", "TEM", "BBMT", "LAMT", "LMT"})


# ---------------------------------------------------------------------------
# BBox
# ---------------------------------------------------------------------------


@dataclass
class BBox:
    """Geographic bounding box in decimal degrees (WGS-84 default).

    Parameters
    ----------
    lat_min, lat_max : float
        Latitude bounds in decimal degrees (south → north).
    lon_min, lon_max : float
        Longitude bounds in decimal degrees (west → east).

    Examples
    --------
    ::

        bbox = BBox(27.8, 28.9, 101.5, 103.2)
        assert 28.0 in bbox  # latitude membership test
        print(bbox.centre)  # (28.35, 102.35)
        print(bbox.area_deg2)  # 1.1 × 1.7 ≈ 1.87 deg²
    """

    lat_min: float
    lat_max: float
    lon_min: float
    lon_max: float

    def __post_init__(self) -> None:
        if self.lat_min > self.lat_max:
            raise ValueError("lat_min must be ≤ lat_max")
        if self.lon_min > self.lon_max:
            raise ValueError("lon_min must be ≤ lon_max")

    @property
    def centre(self) -> tuple[float, float]:
        """Return (lat, lon) centre of the box."""
        return (
            (self.lat_min + self.lat_max) / 2.0,
            (self.lon_min + self.lon_max) / 2.0,
        )

    @property
    def area_deg2(self) -> float:
        """Return the area of the box in squared degrees."""
        return (self.lat_max - self.lat_min) * (self.lon_max - self.lon_min)

    @property
    def span_lat(self) -> float:
        return self.lat_max - self.lat_min

    @property
    def span_lon(self) -> float:
        return self.lon_max - self.lon_min

    def contains(self, lat: float, lon: float) -> bool:
        """Return True when *(lat, lon)* lies inside or on the boundary."""
        return (
            self.lat_min <= lat <= self.lat_max and self.lon_min <= lon <= self.lon_max
        )

    def __contains__(self, latlon: Any) -> bool:
        """Support ``(lat, lon) in bbox`` and single latitude check."""
        if isinstance(latlon, (int, float)):
            return self.lat_min <= latlon <= self.lat_max
        try:
            lat, lon = latlon
            return self.contains(lat, lon)
        except (TypeError, ValueError):
            return False

    def to_dict(self) -> dict[str, float]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> BBox:
        return cls(**{k: float(v) for k, v in d.items()})

    @classmethod
    def from_coords(
        cls,
        lats: list[float],
        lons: list[float],
    ) -> BBox:
        """Build a tight bounding box from lists of coordinates."""
        if not lats or not lons:
            raise ValueError("lats and lons must be non-empty")
        return cls(
            lat_min=min(lats),
            lat_max=max(lats),
            lon_min=min(lons),
            lon_max=max(lons),
        )

    def __repr__(self) -> str:
        return (
            f"BBox(lat=[{self.lat_min:.4f}, {self.lat_max:.4f}]"
            f"  lon=[{self.lon_min:.4f}, {self.lon_max:.4f}])"
        )


# ---------------------------------------------------------------------------
# SurveyMeta
# ---------------------------------------------------------------------------


@dataclass
class SurveyMeta:
    """Campaign-level descriptor for an MT/AMT/CSAMT/TEM survey.

    Parameters
    ----------
    name : str
        Short survey identifier (e.g. ``"WILLY_L18_2023"``).
    project : str, optional
        Parent project name (e.g. ``"Copper-gold exploration"``).
    operator : str, optional
        Acquisition company or institution.
    method : str, default ``"MT"``
        EM method: ``"MT"``, ``"AMT"``, ``"CSAMT"``, ``"TEM"``,
        ``"CSEM"``, ``"BBMT"``, ``"LAMT"``, or ``"LMT"``.
    bbox : BBox, optional
        Geographic bounding box of all station locations.
    date_start : date, optional
        First day of field acquisition.
    date_end : date, optional
        Last day of field acquisition.
    crs : str, default ``"WGS84"``
        Coordinate reference system name.
    n_stations : int, optional
        Total number of stations.
    notes : str
        Free-form annotation.
    extra : dict
        Arbitrary additional fields preserved on round-trips.

    Examples
    --------
    ::

        meta = SurveyMeta(
            name="WILLY_L18",
            project="Phase-III drill targeting",
            operator="EarthAI-Tech",
            method="AMT",
            bbox=BBox(27.8, 28.9, 101.5, 103.2),
            n_stations=128,
        )
        print(meta.duration_days)  # None (no dates set)

        # build from a Sites collection:
        meta2 = SurveyMeta.from_sites(sites, name="survey_A")

        # JSON round-trip:
        meta.to_json("survey_meta.json")
        meta3 = SurveyMeta.from_json("survey_meta.json")
    """

    name: str
    project: str | None = None
    operator: str | None = None
    method: str = "MT"
    bbox: BBox | None = None
    date_start: date | None = None
    date_end: date | None = None
    crs: str = "WGS84"
    n_stations: int | None = None
    notes: str = ""
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        m = self.method.upper()
        if m not in _VALID_METHODS:
            raise ValueError(
                f"method {self.method!r} is not one of " f"{sorted(_VALID_METHODS)}."
            )
        self.method = m
        if self.date_start and self.date_end:
            if self.date_start > self.date_end:
                raise ValueError("date_start must be ≤ date_end")

    # ------------------------------------------------------------------
    # Computed properties
    # ------------------------------------------------------------------

    @property
    def duration_days(self) -> int | None:
        """Number of acquisition days, or None when dates are incomplete."""
        if self.date_start and self.date_end:
            return (self.date_end - self.date_start).days
        return None

    @property
    def is_complete(self) -> bool:
        """Return True when all key descriptors are populated."""
        return all(
            [
                self.name,
                self.method,
                self.bbox is not None,
                self.n_stations is not None,
            ]
        )

    # ------------------------------------------------------------------
    # Integration helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_sites(
        cls,
        sites: Any,
        name: str = "survey",
        method: str = "MT",
        **kwargs: Any,
    ) -> SurveyMeta:
        """Build a :class:`SurveyMeta` from a Sites collection.

        Parameters
        ----------
        sites :
            :class:`~pycsamt.site.base.Sites` or any iterable of Site-like
            objects exposing ``.coords`` (lat, lon, elev).
        name : str
            Survey name to assign.
        method : str
            EM method label.
        **kwargs :
            Additional fields forwarded to the constructor.
        """
        lats, lons = [], []
        n = 0
        for s in sites:
            try:
                lat, lon, _ = s.coords
                if lat is not None and lon is not None:
                    lats.append(float(lat))
                    lons.append(float(lon))
            except (AttributeError, TypeError, ValueError):
                pass
            n += 1

        bbox = BBox.from_coords(lats, lons) if lats and lons else None
        return cls(
            name=name,
            method=method,
            bbox=bbox,
            n_stations=n,
            **kwargs,
        )

    def update_edi_head(self, head: Any) -> None:
        """Push survey fields into an EDI :class:`~pycsamt.seg.heads.Head`.

        Only non-None fields are written so existing values are not
        overwritten unless a replacement is explicitly provided.

        Parameters
        ----------
        head : Head
            An EDI :class:`~pycsamt.seg.heads.Head` instance.
        """
        if self.project is not None:
            try:
                head.project = self.project
            except AttributeError:
                pass
        if self.operator is not None:
            try:
                head.acqby = self.operator
            except AttributeError:
                pass

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """Return all fields as a plain JSON-serialisable dict."""
        d: dict[str, Any] = {
            "name": self.name,
            "project": self.project,
            "operator": self.operator,
            "method": self.method,
            "bbox": self.bbox.to_dict() if self.bbox else None,
            "date_start": self.date_start.isoformat() if self.date_start else None,
            "date_end": self.date_end.isoformat() if self.date_end else None,
            "crs": self.crs,
            "n_stations": self.n_stations,
            "notes": self.notes,
        }
        d.update(self.extra)
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> SurveyMeta:
        """Reconstruct from a plain dict (as produced by :meth:`to_dict`)."""
        d = dict(d)
        bbox_raw = d.pop("bbox", None)
        bbox = BBox.from_dict(bbox_raw) if bbox_raw else None

        def _parse_date(v: Any) -> date | None:
            if v is None:
                return None
            if isinstance(v, date):
                return v
            return date.fromisoformat(str(v))

        known = {
            "name",
            "project",
            "operator",
            "method",
            "date_start",
            "date_end",
            "crs",
            "n_stations",
            "notes",
        }
        extra = {k: v for k, v in d.items() if k not in known}
        return cls(
            name=d.get("name", ""),
            project=d.get("project"),
            operator=d.get("operator"),
            method=d.get("method", "MT"),
            bbox=bbox,
            date_start=_parse_date(d.get("date_start")),
            date_end=_parse_date(d.get("date_end")),
            crs=d.get("crs", "WGS84"),
            n_stations=d.get("n_stations"),
            notes=d.get("notes", ""),
            extra=extra,
        )

    def to_json(self, path: str | Path | None = None) -> str:
        """Serialise to a JSON string (optionally write to *path*)."""
        text = json.dumps(self.to_dict(), indent=2, default=str)
        if path is not None:
            Path(path).write_text(text, encoding="utf-8")
        return text

    @classmethod
    def from_json(cls, path: str | Path) -> SurveyMeta:
        """Load from a JSON file."""
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_dict(data)

    def to_yaml(self, path: str | Path | None = None) -> str:
        """Serialise to a YAML string (optionally write to *path*)."""
        try:
            import yaml  # noqa: PLC0415
        except ImportError as exc:
            raise ImportError("PyYAML is required for YAML serialisation.") from exc
        text = yaml.safe_dump(
            self.to_dict(), default_flow_style=False, allow_unicode=True
        )
        if path is not None:
            Path(path).write_text(text, encoding="utf-8")
        return text

    @classmethod
    def from_yaml(cls, path: str | Path) -> SurveyMeta:
        """Load from a YAML file."""
        try:
            import yaml  # noqa: PLC0415
        except ImportError as exc:
            raise ImportError("PyYAML is required for YAML deserialisation.") from exc
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        return cls.from_dict(data)

    def summary(self) -> str:
        """Return a compact one-paragraph summary string."""
        parts = [f"Survey: {self.name!r}  method={self.method}"]
        if self.project:
            parts.append(f"  Project : {self.project}")
        if self.operator:
            parts.append(f"  Operator: {self.operator}")
        if self.bbox:
            c = self.bbox.centre
            parts.append(
                f"  Area    : lat [{self.bbox.lat_min:.2f}, {self.bbox.lat_max:.2f}]  "
                f"lon [{self.bbox.lon_min:.2f}, {self.bbox.lon_max:.2f}]  "
                f"centre ({c[0]:.3f}°N, {c[1]:.3f}°E)"
            )
        if self.n_stations is not None:
            parts.append(f"  Stations: {self.n_stations}")
        if self.date_start:
            parts.append(
                f"  Dates   : {self.date_start} → {self.date_end or '?'}"
                + (f"  ({self.duration_days} days)" if self.duration_days else "")
            )
        if self.notes:
            parts.append(f"  Notes   : {self.notes}")
        return "\n".join(parts)

    def __repr__(self) -> str:
        n = self.n_stations
        return (
            f"SurveyMeta(name={self.name!r}  method={self.method}"
            + (f"  n={n}" if n is not None else "")
            + ")"
        )
