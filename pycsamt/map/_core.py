# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Core data extraction helpers for :mod:`pycsamt.map`."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from pycsamt.emtools._core import ensure_sites


@dataclass
class StationRecord:
    """A normalized station row used by map renderers."""

    id: str
    latitude: float | None = None
    longitude: float | None = None
    elevation: float | None = None
    line: str | None = None
    index: int = 0
    source: Any = None

    def __post_init__(self) -> None:
        """Normalize station metadata."""
        self.id = str(self.id)
        self.latitude = _finite_or_none(self.latitude)
        self.longitude = _finite_or_none(self.longitude)
        self.elevation = _finite_or_none(self.elevation)
        if self.line not in (None, ""):
            self.line = str(self.line)


@dataclass
class ProfileLine:
    """A named sequence of stations in one profile."""

    name: str
    stations: tuple[StationRecord, ...] = ()


@dataclass
class MapData:
    """Normalized survey data shared by map renderers."""

    sites: Any
    stations: tuple[StationRecord, ...] = ()
    profiles: tuple[ProfileLine, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize station/profile containers."""
        self.stations = tuple(self.stations or ())
        if self.profiles:
            self.profiles = tuple(self.profiles)
        else:
            self.profiles = _build_profiles(self.stations)
        self.metadata.setdefault(
            "n_stations",
            len(self.stations),
        )
        self.metadata.setdefault(
            "n_profiles",
            len(self.profiles),
        )

    @property
    def station_ids(self) -> tuple[str, ...]:
        """Return normalized station IDs."""
        return tuple(station.id for station in self.stations)

    @property
    def lines(self) -> tuple[str, ...]:
        """Return available profile/line names."""
        return tuple(
            profile.name for profile in self.profiles
        )

    @property
    def has_geo(self) -> bool:
        """Whether stations have finite map coordinates."""
        if not self.stations:
            return False
        return all(
            station.latitude is not None
            and station.longitude is not None
            for station in self.stations
        )

    def iter_edis(self) -> tuple[Any, ...]:
        """Return EDI-like source objects."""
        return _as_edi_tuple(self.sites)


TENSOR_COMPONENTS: dict[str, tuple[int, int]] = {
    "xx": (0, 0),
    "xy": (0, 1),
    "yx": (1, 0),
    "yy": (1, 1),
}

COMPONENT_ALIASES: dict[str, str] = {
    "zxx": "xx",
    "zxy": "xy",
    "zyx": "yx",
    "zyy": "yy",
    "rho_xy": "xy",
    "rho_yx": "yx",
    "phase_xy": "xy",
    "phase_yx": "yx",
    "determinant": "det",
    "det": "det",
    "avg": "avg",
    "average": "avg",
    "mean": "avg",
}


@dataclass(frozen=True)
class ComponentSpec:
    """Parsed impedance component request."""

    name: str
    indices: tuple[int, int] | None = None
    mode: str = "tensor"


@dataclass(frozen=True)
class FrequencySelection:
    """Nearest-frequency selection metadata."""

    requested: float | None
    actual: float
    index: int
    delta: float
    relative_delta: float
    within_tolerance: bool = True


@dataclass(frozen=True)
class FrequencyValue:
    """Station value and selected-frequency metadata."""

    station: str
    value: float
    selection: FrequencySelection


def ensure_map_data(
    data: Any,
    *,
    recursive: bool = True,
    line_map: dict[str, Iterable[str]] | None = None,
    verbose: int = 0,
) -> MapData:
    """Return normalized map data.

    Parameters
    ----------
    data :
        EDI path, directory, iterable, or ``Sites``.
    recursive :
        Passed to :func:`pycsamt.emtools._core.ensure_sites`.
    line_map :
        Optional ``line -> station names`` mapping used when
        line metadata is not embedded in EDI objects.
    verbose :
        Verbosity passed to :func:`ensure_sites`.
    """
    if isinstance(data, MapData):
        return data

    sites = ensure_sites(
        data,
        recursive=recursive,
        verbose=verbose,
    )
    edis = _as_edi_tuple(sites)
    stations = tuple(
        _station_record(edi, index=i, line_map=line_map)
        for i, edi in enumerate(edis)
    )
    profiles = _build_profiles(stations)
    metadata = {
        "n_stations": len(stations),
        "n_profiles": len(profiles),
    }
    return MapData(
        sites=sites,
        stations=stations,
        profiles=profiles,
        metadata=metadata,
    )


def station_records(
    data: Any,
    **kwargs: Any,
) -> tuple[StationRecord, ...]:
    """Return normalized station records for *data*."""
    return ensure_map_data(data, **kwargs).stations


def profile_lines(
    data: Any,
    **kwargs: Any,
) -> tuple[ProfileLine, ...]:
    """Return normalized profile lines for *data*."""
    return ensure_map_data(data, **kwargs).profiles


def station_dataframe(data: MapData) -> pd.DataFrame:
    """Return station records as a DataFrame."""
    rows = [
        {
            "ID": s.id,
            "Latitude": s.latitude,
            "Longitude": s.longitude,
            "Elevation": s.elevation,
            "Line": s.line or "line",
            "Index": s.index,
        }
        for s in data.stations
    ]
    return pd.DataFrame(rows)


def normalize_component(component: str | None) -> str:
    """Return the canonical map component name."""
    raw = str(component or "xy").strip().lower()
    raw = raw.replace("-", "_")
    name = COMPONENT_ALIASES.get(raw, raw)
    if name in TENSOR_COMPONENTS or name in {"det", "avg"}:
        return name
    choices = sorted([*TENSOR_COMPONENTS, "avg", "det"])
    msg = (
        f"Unknown impedance component {component!r}. "
        f"Expected one of {choices}."
    )
    raise ValueError(msg)


def component_spec(component: str | None) -> ComponentSpec:
    """Return a parsed impedance component specification."""
    name = normalize_component(component)
    if name in TENSOR_COMPONENTS:
        return ComponentSpec(
            name=name,
            indices=TENSOR_COMPONENTS[name],
        )
    return ComponentSpec(name=name, mode=name)


def component_index(component: str) -> tuple[int, int]:
    """Return tensor indices for an impedance component."""
    spec = component_spec(component)
    if spec.indices is None:
        msg = (
            f"Component {component!r} is derived and has no "
            "single tensor index."
        )
        raise ValueError(msg)
    return spec.indices


def component_values(
    arr: Any,
    component: str | None,
    *,
    quantity: str = "rho",
) -> np.ndarray:
    """Extract a 1-D value series for one component."""
    values = np.asarray(arr, dtype=float)
    if values.ndim < 3 or values.shape[-2:] != (2, 2):
        msg = (
            "Expected impedance values with shape "
            "(nfreq, 2, 2)."
        )
        raise ValueError(msg)
    spec = component_spec(component)
    if spec.indices is not None:
        r, c = spec.indices
        return values[:, r, c]
    xy = values[:, 0, 1]
    yx = values[:, 1, 0]
    if spec.name == "det":
        if quantity.lower() in {"phase", "phi"}:
            return np.nanmean(
                np.column_stack([xy, yx]),
                axis=1,
            )
        return np.sqrt(np.abs(xy * yx))
    return np.nanmean(
        np.column_stack([xy, yx]),
        axis=1,
    )


def frequency_axis(data: MapData) -> np.ndarray:
    """Return the first finite frequency axis found."""
    for edi in data.iter_edis():
        z_obj = getattr(edi, "Z", None)
        if z_obj is None:
            continue
        freq = np.asarray(
            getattr(z_obj, "freq", []),
            dtype=float,
        )
        freq = freq[np.isfinite(freq) & (freq > 0)]
        if freq.size:
            return freq
    return np.array([], dtype=float)


def select_frequency(
    frequencies: Any,
    requested: float | None = None,
    *,
    tolerance: float | None = None,
) -> FrequencySelection | None:
    """Select the closest finite positive frequency."""
    freq = np.asarray(frequencies, dtype=float)
    good = np.isfinite(freq) & (freq > 0)
    if not good.any():
        return None
    valid_idx = np.flatnonzero(good)
    valid = freq[good]
    if requested is None:
        local_idx = 0
        req = None
        delta = 0.0
        rel = 0.0
    else:
        req = float(requested)
        local_idx = int(np.nanargmin(np.abs(valid - req)))
        delta = float(abs(valid[local_idx] - req))
        rel = delta / abs(req) if req else delta
    actual = float(valid[local_idx])
    within = True
    if tolerance is not None and requested is not None:
        within = delta <= float(tolerance)
    return FrequencySelection(
        requested=req,
        actual=actual,
        index=int(valid_idx[local_idx]),
        delta=delta,
        relative_delta=float(rel),
        within_tolerance=within,
    )


def station_distance_km(data: MapData) -> np.ndarray:
    """Return approximate station distance in km."""
    stations = list(data.stations)
    if not stations:
        return np.array([], dtype=float)
    lat = np.array(
        [
            s.latitude if s.latitude is not None else np.nan
            for s in stations
        ],
        dtype=float,
    )
    lon = np.array(
        [
            s.longitude if s.longitude is not None else np.nan
            for s in stations
        ],
        dtype=float,
    )
    if (
        not np.isfinite(lat).all()
        or not np.isfinite(lon).all()
    ):
        return np.arange(len(stations), dtype=float)
    x = lon * 111.320 * np.cos(np.deg2rad(np.nanmean(lat)))
    y = lat * 110.574
    dist = np.zeros(len(stations), dtype=float)
    for i in range(1, len(stations)):
        dist[i] = dist[i - 1] + float(np.hypot(
            x[i] - x[i - 1],
            y[i] - y[i - 1],
        ))
    return dist


def value_at_frequency(
    data: MapData,
    *,
    frequency: float,
    quantity: str = "rho",
    component: str = "xy",
    tolerance: float | None = None,
) -> dict[str, float]:
    """Return station values at the closest frequency."""
    details = value_at_frequency_details(
        data,
        frequency=frequency,
        quantity=quantity,
        component=component,
        tolerance=tolerance,
    )
    return {
        station: item.value
        for station, item in details.items()
    }


def value_at_frequency_details(
    data: MapData,
    *,
    frequency: float,
    quantity: str = "rho",
    component: str = "xy",
    tolerance: float | None = None,
) -> dict[str, FrequencyValue]:
    """Return values with selected-frequency metadata."""
    out: dict[str, FrequencyValue] = {}
    for edi in data.iter_edis():
        sid = _station_id_from_edi(edi)
        z_obj = getattr(edi, "Z", None)
        if z_obj is None:
            continue
        freq = np.asarray(
            getattr(z_obj, "freq", []),
            dtype=float,
        )
        selection = select_frequency(
            freq,
            frequency,
            tolerance=tolerance,
        )
        if (
            selection is None
            or not selection.within_tolerance
        ):
            continue
        q_name = quantity.lower()
        if q_name in {"phase", "phi"}:
            arr = getattr(z_obj, "phase", None)
        else:
            arr = getattr(z_obj, "resistivity", None)
        if arr is None:
            continue
        values = component_values(
            arr,
            component,
            quantity=q_name,
        )
        if selection.index >= values.size:
            continue
        value = float(values[selection.index])
        if np.isfinite(value):
            out[sid] = FrequencyValue(
                station=sid,
                value=value,
                selection=selection,
            )
    return out


def skin_depth_at_frequency(
    data: MapData,
    *,
    frequency: float,
    component: str = "xy",
    tolerance: float | None = None,
) -> dict[str, float]:
    """Return skin depth at the closest frequency."""
    details = value_at_frequency_details(
        data,
        frequency=frequency,
        quantity="rho",
        component=component,
        tolerance=tolerance,
    )
    out: dict[str, float] = {}
    for sid, item in details.items():
        freq = item.selection.actual
        value = item.value
        if value > 0 and freq > 0:
            depth = 503.0 * np.sqrt(value / freq)
            out[sid] = float(depth)
    return out


def pseudosection_table(
    data: MapData,
    *,
    quantity: str = "rho",
    component: str = "xy",
) -> pd.DataFrame:
    """Return long-form profile data for pseudosections."""
    rows: list[dict[str, float | str]] = []
    distances = station_distance_km(data)
    dist_by_id = {
        station.id: float(distances[i])
        for i, station in enumerate(data.stations)
    }
    for edi in data.iter_edis():
        sid = _station_id_from_edi(edi)
        z_obj = getattr(edi, "Z", None)
        if z_obj is None:
            continue
        freq = np.asarray(
            getattr(z_obj, "freq", []),
            dtype=float,
        )
        if freq.size == 0:
            continue
        arr = getattr(
            z_obj,
            "resistivity" if quantity == "rho" else "phase",
            None,
        )
        if arr is None:
            continue
        vals = component_values(
            arr,
            component,
            quantity=quantity,
        )
        periods = np.where(freq > 0, 1.0 / freq, np.nan)
        for period, value in zip(periods, vals):
            if (
                not np.isfinite(period)
                or not np.isfinite(value)
            ):
                continue
            if quantity == "rho" and value <= 0:
                continue
            rows.append(
                {
                    "station": sid,
                    "distance": dist_by_id.get(sid, np.nan),
                    "period": float(period),
                    "value": float(value),
                }
            )
    return pd.DataFrame(rows)


def _station_record(
    edi: Any,
    *,
    index: int,
    line_map: dict[str, Iterable[str]] | None,
) -> StationRecord:
    station_id = _station_id_from_edi(edi) or f"S{index:03d}"
    return StationRecord(
        id=station_id,
        latitude=_first_float(
            edi, "latitude", "lat", "Latitude", "LAT"
        ),
        longitude=_first_float(
            edi,
            "longitude",
            "lon",
            "long",
            "Longitude",
            "LON",
        ),
        elevation=_first_float(
            edi, "elevation", "elev", "Elevation", "ELEV"
        ),
        line=_resolve_line(edi, station_id, line_map),
        index=index,
        source=edi,
    )


def _first_float(obj: Any, *names: str) -> float | None:
    for name in names:
        if not hasattr(obj, name):
            continue
        value = getattr(obj, name)
        try:
            value_f = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(value_f):
            return value_f
    return None


def _finite_or_none(value: Any) -> float | None:
    if value is None:
        return None
    try:
        value_f = float(value)
    except (TypeError, ValueError):
        return None
    if np.isfinite(value_f):
        return value_f
    return None


def _station_id_from_edi(edi: Any) -> str:
    for name in ("station", "id", "name"):
        value = getattr(edi, name, None)
        if value not in (None, ""):
            return str(value)
    return ""


def _as_edi_tuple(sites: Any) -> tuple[Any, ...]:
    if sites is None:
        return ()
    if hasattr(sites, "as_list"):
        return tuple(sites.as_list())
    if isinstance(sites, Iterable) and not isinstance(
        sites,
        (str, bytes),
    ):
        return tuple(sites)
    return (sites,)


def _resolve_line(
    edi: Any,
    station_id: str,
    line_map: dict[str, Iterable[str]] | None,
) -> str | None:
    names = (
        "line",
        "profile",
        "survey_line",
        "Line",
        "Profile",
    )
    for name in names:
        if hasattr(edi, name):
            value = getattr(edi, name)
            if value not in (None, ""):
                return str(value)
    if line_map:
        for line, stations in line_map.items():
            if station_id in {str(s) for s in stations}:
                return str(line)
    return None


def _build_profiles(
    stations: tuple[StationRecord, ...],
) -> tuple[ProfileLine, ...]:
    grouped: dict[str, list[StationRecord]] = {}
    for station in stations:
        line = station.line or "line"
        grouped.setdefault(line, []).append(station)
    return tuple(
        ProfileLine(name=name, stations=tuple(items))
        for name, items in grouped.items()
    )
