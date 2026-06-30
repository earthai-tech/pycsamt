# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Elevation / topography sources for :mod:`pycsamt.map`.

Station elevations carried in the EDI metadata are the default terrain
source for 3-D draping.  This module lets callers override them from a
file (CSV / HDF5 / NPZ) or fetch them online from station coordinates,
so :class:`~pycsamt.map.MapView` and the platform can drape real
terrain rather than only the values embedded in the EDIs.
"""

from __future__ import annotations

import base64
import io
from dataclasses import replace
from typing import Any

import numpy as np

from ._core import MapData

__all__ = [
    "apply_elevations",
    "fetch_elevations",
    "parse_elevation_file",
]

_ID_KEYS = ("station", "id", "name", "sta", "station_names")
_ELEV_KEYS = ("elevation", "elev", "z", "alt", "altitude", "height")


def apply_elevations(
    data: MapData,
    elev_map: dict[str, float],
) -> MapData:
    """Return a copy of *data* with station elevations overridden.

    Parameters
    ----------
    data :
        Source survey data.
    elev_map :
        Mapping of ``station_id -> elevation (m)``. Stations not present
        keep their existing elevation.
    """
    if not elev_map:
        return data
    lookup = {str(k): v for k, v in elev_map.items()}
    stations = tuple(
        replace(s, elevation=float(lookup[s.id]))
        if s.id in lookup and _finite(lookup[s.id])
        else s
        for s in data.stations
    )
    return MapData(
        sites=data.sites,
        stations=stations,
        profiles=(),
        metadata=dict(data.metadata),
    )


def fetch_elevations(
    data: MapData,
    *,
    api_name: str = "open_meteo",
) -> dict[str, float]:
    """Fetch station elevations online from their coordinates.

    Uses :func:`pycsamt.gis.utils.get_elevation_from_api` (needs an
    internet connection and ``requests``).

    Returns
    -------
    dict
        ``station_id -> elevation (m)`` for stations with valid
        coordinates and a finite fetched value.
    """
    from pycsamt.gis.utils import get_elevation_from_api

    ids: list[str] = []
    lats: list[float] = []
    lons: list[float] = []
    for s in data.stations:
        if s.latitude is not None and s.longitude is not None:
            ids.append(s.id)
            lats.append(float(s.latitude))
            lons.append(float(s.longitude))
    if not ids:
        return {}
    elevs = np.atleast_1d(
        np.asarray(
            get_elevation_from_api(lats, lons, api_name=api_name),
            dtype=float,
        )
    )
    return {
        sid: float(e)
        for sid, e in zip(ids, elevs)
        if np.isfinite(e)
    }


def parse_elevation_file(
    content: str | bytes,
    filename: str,
) -> dict[str, float]:
    """Parse an uploaded elevation file into ``{station_id: elev}``.

    Supports CSV, HDF5 (``.h5``/``.hdf5``) and NPZ.  The file must hold
    a station-id column/array (one of ``station``, ``id``, ``name``,
    ``sta``) and an elevation column/array (``elevation``, ``elev``,
    ``z``, ``alt``, ``altitude``, ``height``). Returns ``{}`` on any
    parse failure rather than raising.
    """
    raw = _decode(content)
    if raw is None:
        return {}
    name = (filename or "").lower()
    try:
        if name.endswith(".csv"):
            return _parse_csv(raw)
        if name.endswith((".h5", ".hdf5")):
            return _parse_h5(raw)
        if name.endswith(".npz"):
            return _parse_npz(raw)
    except Exception:  # noqa: BLE001 - best-effort parser
        return {}
    return {}


# ── internals ──────────────────────────────────────────


def _decode(content: str | bytes) -> bytes | None:
    if isinstance(content, bytes):
        return content
    if not isinstance(content, str):
        return None
    if "," in content and content.split(",", 1)[0].startswith("data:"):
        content = content.split(",", 1)[1]
    try:
        return base64.b64decode(content)
    except (ValueError, TypeError):
        return None


def _pick(keys: Any, candidates: tuple[str, ...]) -> str | None:
    for key in keys:
        if str(key).strip().lower() in candidates:
            return key
    return None


def _parse_csv(raw: bytes) -> dict[str, float]:
    import pandas as pd

    df = pd.read_csv(io.BytesIO(raw), dtype=str)
    df.columns = [c.strip() for c in df.columns]
    id_col = _pick(df.columns, _ID_KEYS)
    elev_col = _pick(df.columns, _ELEV_KEYS)
    if id_col is None or elev_col is None:
        return {}
    out: dict[str, float] = {}
    for _, row in df.iterrows():
        name = str(row[id_col]).strip()
        val = _to_float(row[elev_col])
        if name and val is not None:
            out[name] = val
    return out


def _parse_h5(raw: bytes) -> dict[str, float]:
    import h5py

    with h5py.File(io.BytesIO(raw), "r") as f:
        id_key = _pick(list(f), _ID_KEYS)
        elev_key = _pick(list(f), _ELEV_KEYS)
        if id_key is None or elev_key is None:
            return {}
        names = [_as_name(n) for n in f[id_key][:]]
        elevs = np.asarray(f[elev_key][:], dtype=float)
        return _zip_clean(names, elevs)


def _parse_npz(raw: bytes) -> dict[str, float]:
    data = np.load(io.BytesIO(raw), allow_pickle=True)
    id_key = _pick(list(data), _ID_KEYS)
    elev_key = _pick(list(data), _ELEV_KEYS)
    if id_key is None or elev_key is None:
        return {}
    names = [_as_name(n) for n in data[id_key]]
    elevs = np.asarray(data[elev_key], dtype=float)
    return _zip_clean(names, elevs)


def _zip_clean(names, elevs) -> dict[str, float]:
    out: dict[str, float] = {}
    for name, e in zip(names, elevs):
        if name and np.isfinite(e):
            out[name] = float(e)
    return out


def _as_name(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode().strip()
    return str(value).strip()


def _to_float(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def _finite(value: Any) -> bool:
    try:
        return np.isfinite(float(value))
    except (TypeError, ValueError):
        return False
