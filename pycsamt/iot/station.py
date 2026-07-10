"""Field-station configuration for IoT-enabled AMT/CSAMT acquisition.

A :class:`DeviceConfig` describes an IoT *node* (its transport, channels,
and role). A :class:`StationConfig` describes the *ground location* that
one or more nodes occupy: coordinates, profile/line membership, electrode
geometry, and sensor orientation. Separating the two mirrors real field
practice, where a recorder may be redeployed across several stations.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import (
    Any,
)

import pandas as pd

from ..api.property import MetadataMixin, PyCSAMTObject
from ..api.view import maybe_wrap_frame
from . import _common as _c

__all__ = [
    "StationConfig",
    "station_table",
]


@dataclass
class StationConfig(PyCSAMTObject, MetadataMixin):
    """Geospatial and acquisition metadata for one field station.

    Parameters
    ----------
    station_id : str
        Stable station identifier (e.g. ``"S01"``). Used as the join key
        across telemetry, QC, and provenance records.
    lat, lon : float, optional
        Geographic coordinates in decimal degrees. Latitude is validated
        to ``[-90, 90]``; longitude accepts both ``[-180, 180]`` and
        ``[0, 360]`` conventions.
    elevation : float, optional
        Ground elevation in metres.
    profile : str, optional
        Survey line/profile label the station belongs to (e.g. ``"L1"``).
    position_m : float, optional
        Chainage (distance along the profile) in metres. Enables ordering
        of stations into a 2D section.
    channels : list of str, optional
        Acquisition channels recorded at the station (e.g.
        ``["ex", "ey", "hx", "hy"]``).
    dipole_length_m : float, optional
        Electric dipole length in metres, used for E-field scaling and
        contact-resistance checks.
    ex_azimuth_deg, ey_azimuth_deg : float, optional
        Orientation of the Ex/Ey electric dipoles, degrees clockwise from
        geographic north.
    device_ids : list of str, optional
        Identifiers of the IoT nodes occupying this station.
    operator : str, optional
        Field operator responsible for the occupation.
    notes : str
        Free-form field notes.
    metadata : dict
        Free-form structured metadata.
    """

    station_id: str
    lat: float | None = None
    lon: float | None = None
    elevation: float | None = None
    profile: str | None = None
    position_m: float | None = None
    channels: list[str] = field(default_factory=list)
    dipole_length_m: float | None = None
    ex_azimuth_deg: float | None = None
    ey_azimuth_deg: float | None = None
    device_ids: list[str] = field(default_factory=list)
    operator: str | None = None
    notes: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise station fields."""
        self.station_id = _c.as_nonempty_str(self.station_id, "station_id")
        self.lat = _c.as_latitude(self.lat)
        self.lon = _c.as_longitude(self.lon)
        self.elevation = _c.as_elevation(self.elevation)
        self.profile = _c.as_optional_str(self.profile, "profile")
        self.position_m = _c.as_optional_finite_float(
            self.position_m, "position_m"
        )
        self.channels = _c.as_channel_list(self.channels)
        self.dipole_length_m = _c.as_optional_positive(
            self.dipole_length_m, "dipole_length_m"
        )
        self.ex_azimuth_deg = _c.as_optional_finite_float(
            self.ex_azimuth_deg, "ex_azimuth_deg"
        )
        self.ey_azimuth_deg = _c.as_optional_finite_float(
            self.ey_azimuth_deg, "ey_azimuth_deg"
        )
        self.device_ids = [
            _c.as_nonempty_str(d, "device_id")
            for d in _c.unique_preserving(list(self.device_ids or []))
        ]
        self.operator = _c.as_optional_str(self.operator, "operator")
        self.notes = str(self.notes or "")
        if not isinstance(self.metadata, dict):
            self.metadata = dict(self.metadata or {})

    @property
    def has_location(self) -> bool:
        """Return whether latitude and longitude are both known."""
        return self.lat is not None and self.lon is not None

    @property
    def coords(self) -> tuple[float, float, float] | None:
        """Return ``(lat, lon, elevation)`` if a location is set."""
        if not self.has_location:
            return None
        return (self.lat, self.lon, float(self.elevation or 0.0))

    @property
    def n_channels(self) -> int:
        """Number of declared acquisition channels."""
        return len(self.channels)

    def attach_device(self, device_id: str) -> StationConfig:
        """Associate an IoT node with this station (idempotent)."""
        key = _c.as_nonempty_str(device_id, "device_id")
        if key not in self.device_ids:
            self.device_ids.append(key)
        return self

    def as_dict(self) -> dict[str, Any]:
        """Return a serialisable station dictionary."""
        return dict(
            station_id=self.station_id,
            lat=self.lat,
            lon=self.lon,
            elevation=self.elevation,
            profile=self.profile,
            position_m=self.position_m,
            channels=list(self.channels),
            dipole_length_m=self.dipole_length_m,
            ex_azimuth_deg=self.ex_azimuth_deg,
            ey_azimuth_deg=self.ey_azimuth_deg,
            device_ids=list(self.device_ids),
            operator=self.operator,
            notes=self.notes,
            metadata=dict(self.metadata),
        )

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> StationConfig:
        """Build a station config from a mapping, ignoring unknown keys."""
        allowed = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in dict(data).items() if k in allowed})


def station_table(
    stations: StationConfig | Iterable[StationConfig],
    *,
    api: bool | None = None,
) -> Any:
    """Return one row per station describing its acquisition metadata."""
    items = (
        [stations] if isinstance(stations, StationConfig)
        else list(stations)
    )
    rows: list[dict[str, Any]] = []
    for station in items:
        station.validate()
        row = station.as_dict()
        row["channels"] = ";".join(row["channels"])
        row["device_ids"] = ";".join(row["device_ids"])
        row["n_channels"] = station.n_channels
        row["has_location"] = station.has_location
        rows.append(row)
    df = pd.DataFrame.from_records(rows)
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_station_table",
        kind="iot.station",
        source=items,
        description="IoT field stations and acquisition metadata.",
    )
