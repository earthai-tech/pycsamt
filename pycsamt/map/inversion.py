# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Import geophysical inversion results into :class:`~pycsamt.map.MapData`.

A ModEM run inverts one 3-D resistivity volume for the whole survey —
individual lines are not separate model files, they are paths of
stations *through* that one volume. :func:`load_modem_lines` slices
that volume into one 2-D vertical curtain per survey line (see
:mod:`pycsamt.models.modem.section`) and returns a :class:`MapData`
carrying real station coordinates plus the precomputed sections, so
the existing 3-D fence/depth builders in :mod:`pycsamt.map.volume`
render inversion-sourced lines exactly like EDI-sourced ones — same
real-geometry line placement, same UI.

Geo-referencing station coordinates
------------------------------------
Three sources are tried, in order, for each ModEM station:

1. ``known_stations`` — a previously-loaded EDI :class:`MapData`'s
   stations, matched by station id. Recommended: works for any
   inversion backend and lets a matched station's line/elevation
   override the (less complete) values recoverable from the ModEM
   files alone.
2. The ModEM ``.dat`` file's own ``GG_Lat``/``GG_Lon`` columns
   (see :attr:`pycsamt.models.modem.data.ModEmData.site_lonlat`).
3. Neither — the station keeps ``latitude=longitude=None`` and the
   3-D builder falls back to synthetic index-based line spacing for
   that line (see :mod:`pycsamt.map.geometry`).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np

from ._core import MapData, StationRecord, normalize_station_id

__all__ = ["group_modem_stations", "load_modem_lines"]


def group_modem_stations(
    station_names: Iterable[str],
    *,
    known_stations: Iterable[StationRecord] | None = None,
) -> dict[str, list[str]]:
    """Group ModEM station names into survey lines.

    Prefers matching each name against *known_stations* (e.g.
    previously loaded EDI stations) and using their ``line`` tag.
    Falls back to parsing the line token out of the station-name
    convention ``{survey}-{line}-{station}{suffix}`` (e.g.
    ``23-18-001A`` -> line ``18``) — a heuristic, used only when a
    station has no match in *known_stations*.
    """
    known_by_id = _index_known_stations(known_stations)
    groups: dict[str, list[str]] = {}
    for name in station_names:
        match = _lookup(known_by_id, name)
        line = str(match.line) if match is not None and match.line else _line_token(name)
        groups.setdefault(line, []).append(str(name))
    return groups


def _line_token(name: str) -> str:
    """Best-effort line id from a ``{survey}-{line}-{station}`` name."""
    parts = str(name).split("-")
    if len(parts) >= 3:
        return parts[1]
    if len(parts) == 2:
        return parts[0]
    return "line"


def _index_known_stations(
    known_stations: Iterable[StationRecord] | None,
) -> dict[str, StationRecord]:
    index: dict[str, StationRecord] = {}
    for station in known_stations or ():
        index[str(station.id)] = station
        index[str(station.id).strip().lower()] = station
        # Normalized fallback tier — matches even when the two
        # sources format ids slightly differently (dashes vs
        # underscores vs spaces, mixed case, ...).
        index.setdefault(normalize_station_id(station.id), station)
    return index


def _lookup(
    known_by_id: dict[str, StationRecord],
    name: str,
) -> StationRecord | None:
    exact = known_by_id.get(str(name)) or known_by_id.get(
        str(name).strip().lower()
    )
    if exact is not None:
        return exact
    return known_by_id.get(normalize_station_id(name))


def load_modem_lines(
    folder: str | Path,
    *,
    known_stations: Iterable[StationRecord] | None = None,
    fetch_elevation: bool = True,
    verbose: int = 0,
) -> MapData:
    """Load a ModEM 3-D inversion result folder as a multi-line MapData.

    Parameters
    ----------
    folder : path-like
        A ModEM output directory. The matching final-iteration
        ``.rho``/``.dat`` pair is auto-detected by
        :class:`pycsamt.models.modem.results.InversionResult`.
    known_stations : iterable of StationRecord, optional
        Previously-loaded EDI stations (e.g. ``existing_map_data.stations``)
        used to geo-reference and group ModEM stations by real
        coordinates/line name — see the module docstring.
    fetch_elevation : bool, default True
        ModEM output carries no real elevation (unlike EDI, where
        it's already in the file header, so "Drape topography" just
        works). When ``True``, any station still missing an
        elevation after ``known_stations`` matching gets a
        best-effort online lookup (Open-Meteo) so topography isn't
        silently flat by default. Failures (offline, API error, …)
        are swallowed — elevation simply stays unset, same as
        passing ``False``.
    verbose : int, default 0
        Verbosity forwarded to the ModEM readers.

    Returns
    -------
    MapData
        ``sites=None`` (no EDI backing); ``stations`` carries one
        :class:`StationRecord` per ModEM station with coordinates
        resolved where possible; ``metadata["sections"]`` carries the
        precomputed per-line ``(x, z, rho)`` curtains consumed
        directly by :mod:`pycsamt.map.volume`.
    """
    from pycsamt.models.modem.results import InversionResult
    from pycsamt.models.modem.section import station_curtain

    result = InversionResult(str(folder), verbose=verbose)
    model = result.model_final
    data = result.data_obs
    if model is None or data is None or not data.site_names:
        msg = f"No usable ModEM model/data found under {folder!r}"
        raise ValueError(msg)

    known_by_id = _index_known_stations(known_stations)
    groups = group_modem_stations(data.site_names, known_stations=known_stations)

    stations: list[StationRecord] = []
    sections: dict[str, dict[str, Any]] = {}
    for idx, (line, names) in enumerate(groups.items()):
        curtain = station_curtain(model, data, names)
        if not curtain.station_names:
            continue
        elevs: list[float] = []
        for name in curtain.station_names:
            match = _lookup(known_by_id, name)
            lon, lat = _resolve_lonlat(name, data, match)
            elev = match.elevation if match is not None else None
            stations.append(
                StationRecord(
                    id=name,
                    latitude=lat,
                    longitude=lon,
                    elevation=elev,
                    line=str(line),
                    index=idx,
                )
            )
            elevs.append(elev if elev is not None else np.nan)
        sections[str(line)] = {
            "z": curtain.z,
            "rho": curtain.rho,
            "stations": np.array(curtain.station_names, dtype=object),
            "elev": np.array(elevs, dtype=float),
        }

    metadata = {
        "source": "modem",
        "workdir": str(folder),
        "sections": sections,
        "rms": result.final_rms,
    }
    data = MapData(sites=None, stations=tuple(stations), metadata=metadata)
    if fetch_elevation:
        data = _try_fetch_elevations(data)
    return data


def _try_fetch_elevations(data: MapData) -> MapData:
    """Best-effort online elevation fetch for stations missing one.

    Only touches stations that don't already have an elevation (e.g.
    from ``known_stations`` matching), so a real, previously-resolved
    value is never overwritten by a coarser online lookup.
    """
    from .topo import apply_elevations, fetch_elevations

    missing = tuple(s for s in data.stations if s.elevation is None)
    if not missing:
        return data
    probe = MapData(sites=None, stations=missing)
    try:
        elev_map = fetch_elevations(probe)
    except Exception:  # noqa: BLE001 - never let a network hiccup break the import
        return data
    if not elev_map:
        return data
    return apply_elevations(data, elev_map)


def _resolve_lonlat(
    name: str,
    data: Any,
    known_match: StationRecord | None,
) -> tuple[float | None, float | None]:
    if (
        known_match is not None
        and known_match.longitude is not None
        and known_match.latitude is not None
    ):
        return known_match.longitude, known_match.latitude
    lonlat = data.site_lonlat.get(name)
    if lonlat is not None:
        return float(lonlat[0]), float(lonlat[1])
    return None, None
