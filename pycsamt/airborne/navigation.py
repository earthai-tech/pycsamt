# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Navigation and attitude containers for airborne EM surveys."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..core.base import CoreObject
from ..metadata import BBox

__all__ = ["NavigationTrack"]


def _normalize_sample_ids(values: Any) -> tuple[str, ...]:
    if values is None:
        raise ValueError("sample_ids must be provided")
    out = tuple(str(value).strip() for value in values)
    if not out:
        raise ValueError("sample_ids must be non-empty")
    if any(not value for value in out):
        raise ValueError("sample_ids must not contain empty identifiers")
    if len(out) != len(set(out)):
        raise ValueError("sample_ids must be unique")
    return out


def _normalize_numeric_vector(
    value: Any | None,
    *,
    name: str,
    size: int,
) -> np.ndarray | None:
    if value is None:
        return None
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1-D array")
    if arr.size != size:
        raise ValueError(
            f"{name} length must match sample_ids: {arr.size} != {size}"
        )
    if np.any(np.isinf(arr)):
        raise ValueError(f"{name} must not contain infinite values")
    return arr


def _normalize_object_vector(
    value: Any | None,
    *,
    name: str,
    size: int,
) -> tuple[Any, ...] | None:
    if value is None:
        return None
    out = tuple(value)
    if len(out) != size:
        raise ValueError(
            f"{name} length must match sample_ids: {len(out)} != {size}"
        )
    return out


@dataclass(repr=False)
class NavigationTrack(CoreObject):
    """Sample-aligned navigation and attitude information for one line.

    The class intentionally contains only concepts that are broadly common
    across airborne EM systems.  Technology-specific fields should be kept in
    :attr:`attrs` until genuine delivery data justify a stable public field.

    Parameters
    ----------
    sample_ids : sequence of str-like
        Ordered identifiers defining the common sample axis.
    latitude, longitude : array-like, optional
        Geographic coordinates in decimal degrees. They must be supplied
        together when used. NaN values are allowed for individual missing
        samples, while finite values are range-checked.
    easting, northing : array-like, optional
        Projected coordinates. They must be supplied together when used.
    terrain_elevation, platform_elevation : array-like, optional
        Elevations on a common vertical datum, normally metres.
    clearance : array-like, optional
        Explicit receiver/platform clearance above ground. If absent,
        :attr:`clearance_values` derives it from platform minus terrain where
        both elevations are available.
    heading, pitch, roll : array-like, optional
        Attitude channels in source-system angular convention. They are
        retained without imposing a proprietary sign convention.
    timestamps : sequence, optional
        Sample-aligned acquisition times. Values are preserved as supplied.
    crs : str, optional
        Coordinate reference system label for projected coordinates.
    datum : str, default ``"WGS84"``
        Geographic datum label.
    attrs : dict
        Extension point for system-specific navigation metadata.
    """

    sample_ids: Any
    latitude: Any | None = None
    longitude: Any | None = None
    easting: Any | None = None
    northing: Any | None = None
    terrain_elevation: Any | None = None
    platform_elevation: Any | None = None
    clearance: Any | None = None
    heading: Any | None = None
    pitch: Any | None = None
    roll: Any | None = None
    timestamps: Any | None = None
    crs: str | None = None
    datum: str | None = "WGS84"
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.sample_ids = _normalize_sample_ids(self.sample_ids)
        n = len(self.sample_ids)
        for name in (
            "latitude",
            "longitude",
            "easting",
            "northing",
            "terrain_elevation",
            "platform_elevation",
            "clearance",
            "heading",
            "pitch",
            "roll",
        ):
            setattr(
                self,
                name,
                _normalize_numeric_vector(
                    getattr(self, name),
                    name=name,
                    size=n,
                ),
            )
        self.timestamps = _normalize_object_vector(
            self.timestamps,
            name="timestamps",
            size=n,
        )

        if (self.latitude is None) != (self.longitude is None):
            raise ValueError(
                "latitude and longitude must be supplied together"
            )
        if (self.easting is None) != (self.northing is None):
            raise ValueError("easting and northing must be supplied together")

        if self.latitude is not None:
            finite = np.isfinite(self.latitude)
            if np.any((self.latitude[finite] < -90.0) | (
                self.latitude[finite] > 90.0
            )):
                raise ValueError(
                    "finite latitude values must be within [-90, 90]"
                )
            finite = np.isfinite(self.longitude)
            if np.any((self.longitude[finite] < -180.0) | (
                self.longitude[finite] > 180.0
            )):
                raise ValueError(
                    "finite longitude values must be within [-180, 180]"
                )

        for attr in ("crs", "datum"):
            value = getattr(self, attr)
            if value is not None:
                text = str(value).strip()
                setattr(self, attr, text or None)
        self.attrs = dict(self.attrs or {})

    @property
    def n_samples(self) -> int:
        """Number of navigation samples."""
        return len(self.sample_ids)

    @property
    def has_geographic_coordinates(self) -> bool:
        """Whether latitude/longitude arrays are present."""
        return self.latitude is not None and self.longitude is not None

    @property
    def has_projected_coordinates(self) -> bool:
        """Whether easting/northing arrays are present."""
        return self.easting is not None and self.northing is not None

    @property
    def clearance_values(self) -> np.ndarray | None:
        """Return explicit or safely derived clearance values.

        Explicit clearance always takes precedence.  When it is unavailable,
        the difference ``platform_elevation - terrain_elevation`` is returned
        without changing the stored metadata.
        """
        if self.clearance is not None:
            return np.array(self.clearance, copy=True)
        if self.platform_elevation is None or self.terrain_elevation is None:
            return None
        return np.asarray(self.platform_elevation) - np.asarray(
            self.terrain_elevation
        )

    @property
    def bbox(self) -> BBox | None:
        """Tight geographic bounding box over finite coordinates."""
        if not self.has_geographic_coordinates:
            return None
        mask = np.isfinite(self.latitude) & np.isfinite(self.longitude)
        if not np.any(mask):
            return None
        return BBox.from_coords(
            self.latitude[mask].tolist(),
            self.longitude[mask].tolist(),
        )

    def index_of(self, sample_id: str) -> int:
        """Return the navigation index for *sample_id*."""
        key = str(sample_id).strip()
        try:
            return self.sample_ids.index(key)
        except ValueError as exc:
            raise KeyError(f"unknown airborne sample: {sample_id!r}") from exc
