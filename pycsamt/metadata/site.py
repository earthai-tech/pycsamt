# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Site-level metadata for electromagnetic transfer functions."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from math import isfinite
from typing import Any

from ..api.property import PyCSAMTObject

__all__ = ["LocationMeta", "SiteMeta"]


@dataclass(repr=False)
class LocationMeta(PyCSAMTObject):
    """Geographic location and optional magnetic declination metadata."""

    latitude: float | None = None
    longitude: float | None = None
    elevation: float | None = None
    datum: str | None = "WGS84"
    elevation_units: str = "meters"
    declination: float | None = None
    declination_epoch: float | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        if self.latitude is not None:
            lat = float(self.latitude)
            if not isfinite(lat) or not -90.0 <= lat <= 90.0:
                raise ValueError(
                    "latitude must be finite and within [-90, 90]"
                )
            self.latitude = lat
        if self.longitude is not None:
            lon = float(self.longitude)
            if not isfinite(lon) or not -180.0 <= lon <= 180.0:
                raise ValueError(
                    "longitude must be finite and within [-180, 180]"
                )
            self.longitude = lon
        for attr in ("elevation", "declination", "declination_epoch"):
            value = getattr(self, attr)
            if value is None:
                continue
            numeric = float(value)
            if not isfinite(numeric):
                raise ValueError(f"{attr} must be finite")
            setattr(self, attr, numeric)
        if self.datum is not None:
            datum = str(self.datum).strip()
            self.datum = datum or None
        units = str(self.elevation_units).strip()
        if not units:
            raise ValueError("elevation_units must be non-empty")
        self.elevation_units = units
        self.extra = dict(self.extra or {})

    @property
    def has_horizontal_coordinates(self) -> bool:
        return self.latitude is not None and self.longitude is not None


@dataclass(repr=False)
class SiteMeta(PyCSAMTObject):
    """Site identity, acquisition interval, and geographic location."""

    project: str | None = None
    survey: str | None = None
    year_collected: int | None = None
    country: str | None = None
    site_id: str | int | None = None
    name: str | None = None
    location: LocationMeta | None = None
    acquired_by: str | None = None
    start: datetime | str | None = None
    end: datetime | str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        for attr in ("project", "survey", "country", "name", "acquired_by"):
            value = getattr(self, attr)
            if value is not None:
                text = str(value).strip()
                setattr(self, attr, text or None)
        if self.site_id is not None and isinstance(self.site_id, str):
            value = self.site_id.strip()
            self.site_id = value or None
        if self.year_collected is not None:
            year = int(self.year_collected)
            if year < 0:
                raise ValueError("year_collected must be non-negative")
            self.year_collected = year
        if self.location is not None and not isinstance(
            self.location, LocationMeta
        ):
            raise TypeError("location must be a LocationMeta")

        for attr in ("start", "end"):
            value = getattr(self, attr)
            if value is not None and not isinstance(value, (datetime, str)):
                raise TypeError(f"{attr} must be datetime, str, or None")
            if isinstance(value, str):
                text = value.strip()
                setattr(self, attr, text or None)
        if isinstance(self.start, datetime) and isinstance(self.end, datetime):
            if self.start > self.end:
                raise ValueError("site start must not be after site end")
        self.extra = dict(self.extra or {})

    @property
    def preferred_name(self) -> str | None:
        """Return human site name, falling back to the site identifier."""
        if self.name:
            return self.name
        if self.site_id is None:
            return None
        return str(self.site_id)
