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
from pathlib import Path
from typing import Any

import numpy as np

from ._core import MapData, normalize_station_id

__all__ = [
    "apply_elevations",
    "export_elevations",
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
    # Normalized fallback so minor id-formatting differences between
    # the elevation source and the survey's own station ids don't
    # block an otherwise-obvious match.
    norm_lookup = {normalize_station_id(k): v for k, v in elev_map.items()}

    def _resolve(station_id: str) -> float | None:
        if station_id in lookup and _finite(lookup[station_id]):
            return float(lookup[station_id])
        key = normalize_station_id(station_id)
        if key in norm_lookup and _finite(norm_lookup[key]):
            return float(norm_lookup[key])
        return None

    stations = tuple(
        replace(s, elevation=resolved)
        if (resolved := _resolve(s.id)) is not None
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
    return {sid: float(e) for sid, e in zip(ids, elevs) if np.isfinite(e)}


def export_elevations(
    data: MapData,
    path: str | Path,
    *,
    fmt: str | None = None,
) -> Path:
    """Export station id/elevation/coordinates to CSV or HDF5.

    EDI files carry real, field-surveyed elevation; a ModEM (or
    Occam2D/MARE2DEM) inversion result does not. Exporting a survey's
    EDI-derived topography here produces a small, portable lookup
    table — ``station`` + ``elevation`` (plus ``latitude``/
    ``longitude``/``line`` for provenance) — that :func:`parse_elevation_file`
    already knows how to read back in. That means it can be applied
    to an inversion-sourced view later via the "Upload file"
    elevation source, in a *different* session that never reloads
    the original EDIs, matching stations by id (with the same
    normalized fallback used by :func:`apply_elevations`).

    Parameters
    ----------
    data : MapData
        Survey to export (typically EDI-sourced, where ``elevation``
        is already populated).
    path : path-like
        Destination file. ``fmt`` defaults to the file's own suffix
        (``.csv`` / ``.h5`` / ``.hdf5``); pass it explicitly to
        override.
    fmt : {"csv", "h5", "hdf5"}, optional

    Returns
    -------
    pathlib.Path
        The path written to.

    Raises
    ------
    ValueError
        If no station in *data* has a known elevation, or *fmt* is
        unsupported.
    """
    path = Path(path)
    resolved_fmt = (fmt or path.suffix.lstrip(".") or "csv").lower()

    rows = [
        (s.id, float(s.elevation), s.latitude, s.longitude, s.line or "")
        for s in data.stations
        if s.elevation is not None
    ]
    if not rows:
        msg = (
            "No stations with a known elevation to export — drape "
            "topography (station EDIs, upload, or fetch online) first."
        )
        raise ValueError(msg)

    ids, elevs, lats, lons, lines = (list(col) for col in zip(*rows))
    path.parent.mkdir(parents=True, exist_ok=True)

    if resolved_fmt == "csv":
        import pandas as pd

        pd.DataFrame(
            {
                "station": ids,
                "elevation": elevs,
                "latitude": lats,
                "longitude": lons,
                "line": lines,
            }
        ).to_csv(path, index=False)
    elif resolved_fmt in ("h5", "hdf5"):
        import h5py

        str_dtype = h5py.string_dtype()
        with h5py.File(path, "w") as fh:
            fh.create_dataset(
                "station", data=np.asarray(ids, dtype=object), dtype=str_dtype
            )
            fh.create_dataset(
                "elevation", data=np.asarray(elevs, dtype=float)
            )
            fh.create_dataset(
                "latitude",
                data=np.asarray(
                    [v if v is not None else np.nan for v in lats],
                    dtype=float,
                ),
            )
            fh.create_dataset(
                "longitude",
                data=np.asarray(
                    [v if v is not None else np.nan for v in lons],
                    dtype=float,
                ),
            )
            fh.create_dataset(
                "line", data=np.asarray(lines, dtype=object), dtype=str_dtype
            )
    else:
        msg = (
            f"Unsupported export format: {resolved_fmt!r} (use 'csv' or 'h5')"
        )
        raise ValueError(msg)
    return path


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
