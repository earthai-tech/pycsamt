# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Attach real per-station spatial coordinates (and elevation) to a PCSF
adapter's :class:`~pycsamt.format.schema.StationTable` from an external
*topo source*.

Every PCSF adapter (:mod:`pycsamt.format.adapters`) already accepts a
``station_elevations``/``station_lonlat``-style mapping the caller has
to build by hand. This module is the "smart" layer on top of that: a
single ``topo=`` argument that can be

- a ``.bln`` (Golden Software Surfer blanking) file — a bare ``x, y``
  point list, no station identity, matched to stations *positionally*
  (in survey order) since the format itself carries none;
- a ``.stn`` (Zonge station-location) file, parsed by the same
  battle-tested low-level reader :mod:`pycsamt.zonge` already uses
  (:func:`pycsamt.zonge.utils.read_stn`) — never re-implemented here;
- a ``.csv`` (or ``.txt``) file with a header pycsamt can recognise
  (``station``/``lat``/``lon``/``easting``/``northing``/``elevation``
  and common synonyms);
- an already-geo-located ``Sites``/``MapData``/iterable-of-
  ``StationRecord``-like object (duck-typed: anything exposing a
  station id plus ``latitude``/``longitude``/``elevation``) — so a
  survey's own EDI collection can supply topography directly, with no
  intermediate file at all;
- a plain ``{station_name: (lon, lat)}`` / ``{station_name: (lon, lat,
  elevation)}`` mapping — the same shape ``station_lonlat`` already
  accepted, now handled by the same code path;
- for a multiline build, a ``{line_id: <any of the above>}`` mapping,
  or a sequence of one source per line in line order — one ``.bln``/
  ``.stn`` file per surveyed line is a common real-world layout.

Design choices, stated explicitly because they are the crux of "smart"
attribution rather than an implementation detail:

1. **A named source (``.stn``, ``.csv`` with a station column, a
   plain dict, or a Sites/MapData object) matches by station id** —
   exact, then a normalized fallback (:func:`pycsamt.map._core.normalize_station_id`,
   the same one every other cross-source station match in this
   codebase already uses). Unmatched stations on either side are
   reported, never silently dropped or fabricated.
2. **A name-less source (a bare ``.bln``, or a ``.csv`` without a
   station column) is positional** and therefore *requires* an exact
   station-count match — the whole reason to "first detect the number
   of stations passed and compare to the inversion stations" before
   attributing anything. A mismatch raises by default
   (``on_mismatch="raise"``); ``on_mismatch="warn"`` instead issues a
   :class:`UserWarning` and attributes only the overlapping prefix, a
   deliberate best-effort escape hatch rather than the default.
3. **When ``topo`` is given, it takes precedence over any other
   lon/lat source** (``station_lonlat``, a backend's own native
   coordinates) for every station it successfully attributes — with a
   :class:`UserWarning` when both were supplied, so the override is
   never silent. Elevation merges rather than fully overriding: a
   station the topo source has no elevation for keeps whatever
   ``station_elevations`` already gave it.
4. **Passing no ``topo`` behaves exactly as before this module
   existed** — every adapter's existing ``station_elevations``/
   ``station_lonlat`` parameters are unaffected when ``topo`` is
   ``None``.
5. **Projected coordinates (easting/northing) need `epsg` or
   `utm_zone` to become lon/lat** — conversion is delegated entirely
   to :func:`pycsamt.gis.utils.to_ll` (the project's one existing
   UTM/EPSG conversion utility; ``pyproj`` stays an optional
   dependency, imported lazily only when a conversion is actually
   requested, matching PCSF's own "no new mandatory dependency"
   principle in ``SPEC.md`` S2). A source already carrying lat/lon
   columns, or a ``.bln`` explicitly marked ``latlon=True``, needs
   neither.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "TopoTable",
    "TopoAttribution",
    "read_topo_file",
    "topo_from_sites",
    "attribute_topo",
    "resolve_topo",
]

_NAME_ALIASES = ("station", "name", "site", "id", "sta", "dot")
_LAT_ALIASES = ("lat", "latitude")
_LON_ALIASES = ("lon", "long", "longitude")
_EAST_ALIASES = ("easting", "east", "gride", "x", "e")
_NORTH_ALIASES = ("northing", "north", "gridn", "y", "n")
_ELEV_ALIASES = ("elev", "elevation", "alt", "altitude", "z", "h")


def _find_col(columns: Sequence[str], aliases: Sequence[str]) -> str | None:
    lower = {str(c).strip().lower(): c for c in columns}
    for alias in aliases:
        if alias in lower:
            return lower[alias]
    # Substring fallback (handles e.g. "Elev(m)") -- only for aliases
    # long enough that a false-positive substring hit is implausible;
    # a single-letter alias like "n"/"e"/"x"/"y" must match the whole
    # column name exactly (already checked above) or not at all, or it
    # would match almost anything ("note" contains "n").
    for col in columns:
        cl = str(col).strip().lower()
        for alias in aliases:
            if len(alias) >= 3 and alias in cl:
                return col
    return None


@dataclass
class TopoTable:
    """A parsed topo source, before it is matched to any station names.

    Attributes
    ----------
    names : list of str, optional
        Station identity carried by the source itself. ``None`` for a
        positional-only source (e.g. a bare ``.bln``).
    lon, lat : ndarray, shape (n,)
        WGS84 decimal degrees — already converted from
        easting/northing if the source needed that.
    elevation : ndarray, shape (n,), optional
        Metres, ``nan`` where genuinely unknown.
    source : str
        Human-readable provenance (file path, or a short description
        for an in-memory source), surfaced in error/warning messages.
    """

    lon: np.ndarray
    lat: np.ndarray
    elevation: np.ndarray | None = None
    names: list[str] | None = None
    source: str = "<unknown>"

    def __post_init__(self) -> None:
        self.lon = np.asarray(self.lon, dtype=float)
        self.lat = np.asarray(self.lat, dtype=float)
        if self.elevation is not None:
            self.elevation = np.asarray(self.elevation, dtype=float)
        if self.names is not None and len(self.names) != self.lon.shape[0]:
            raise ValueError(
                f"{self.source}: {len(self.names)} name(s) but "
                f"{self.lon.shape[0]} coordinate(s)"
            )

    @property
    def n(self) -> int:
        return int(self.lon.shape[0])


@dataclass
class TopoAttribution:
    """Per-station real coordinates resolved from a topo source, ready
    to populate :attr:`~pycsamt.format.schema.StationTable.lon`/``lat``
    (and merge into elevation)."""

    lon: dict[str, float] = field(default_factory=dict)
    lat: dict[str, float] = field(default_factory=dict)
    elevation: dict[str, float] = field(default_factory=dict)
    matched: list[str] = field(default_factory=list)
    unmatched_stations: list[str] = field(default_factory=list)
    source: str = "none"

    def __bool__(self) -> bool:
        return bool(self.matched)


# ---------------------------------------------------------------------
# Coordinate conversion (lazy pyproj/GDAL via pycsamt.gis.utils)
# ---------------------------------------------------------------------


def _convert_xy(
    x: np.ndarray, y: np.ndarray, *, epsg: int | None, utm_zone: Any | None, source: str
) -> tuple[np.ndarray, np.ndarray]:
    if epsg is None and utm_zone is None:
        raise ValueError(
            f"{source}: gives projected easting/northing but neither "
            "epsg nor utm_zone was provided -- pass one so it can be "
            "converted to lon/lat (see pycsamt.gis.utils.to_ll)."
        )
    from ..gis.utils import to_ll

    lat, lon = to_ll(np.asarray(x, dtype=float), np.asarray(y, dtype=float), zone=utm_zone, epsg=epsg)
    return np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)


# ---------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------


def _read_bln(
    path: Path, *, epsg: int | None, utm_zone: Any | None, latlon: bool
) -> TopoTable:
    """Golden Software Surfer ``.bln``: a ``n_points[,flag]`` header
    line, then *n_points* rows of ``x, y`` (an optional 3rd numeric
    column is read as elevation). The format carries no station
    identity -- points are positional, in file order."""
    raw_lines = path.read_text(encoding="utf-8").splitlines()
    lines = [
        ln.strip()
        for ln in raw_lines
        if ln.strip() and not ln.strip().startswith(("#", "//"))
    ]
    if not lines:
        raise ValueError(f"{path}: empty .bln file")
    header = [p.strip() for p in lines[0].split(",") if p.strip() != ""]
    try:
        n_points = int(float(header[0]))
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{path}: malformed .bln header {lines[0]!r}") from exc

    body = lines[1 : 1 + n_points]
    if len(body) != n_points:
        raise ValueError(
            f"{path}: header declares {n_points} point(s), found {len(body)}"
        )
    rows = [
        [float(v) for v in ln.replace("\t", ",").split(",") if v.strip() != ""]
        for ln in body
    ]
    arr = np.asarray(rows, dtype=float)
    x, y = arr[:, 0], arr[:, 1]
    elevation = arr[:, 2] if arr.shape[1] >= 3 else None
    if latlon:
        lon, lat = x, y
    else:
        lon, lat = _convert_xy(x, y, epsg=epsg, utm_zone=utm_zone, source=str(path))
    return TopoTable(lon=lon, lat=lat, elevation=elevation, names=None, source=str(path))


def _dataframe_to_topotable(
    df: Any, source: str, *, epsg: int | None, utm_zone: Any | None
) -> TopoTable:
    columns = list(df.columns)
    name_col = _find_col(columns, _NAME_ALIASES)
    lat_col = _find_col(columns, _LAT_ALIASES)
    lon_col = _find_col(columns, _LON_ALIASES)
    east_col = _find_col(columns, _EAST_ALIASES)
    north_col = _find_col(columns, _NORTH_ALIASES)
    elev_col = _find_col(columns, _ELEV_ALIASES)

    if lat_col and lon_col:
        lon = df[lon_col].to_numpy(dtype=float)
        lat = df[lat_col].to_numpy(dtype=float)
    elif east_col and north_col:
        lon, lat = _convert_xy(
            df[east_col].to_numpy(dtype=float),
            df[north_col].to_numpy(dtype=float),
            epsg=epsg,
            utm_zone=utm_zone,
            source=str(source),
        )
    else:
        raise ValueError(
            f"{source}: could not find lat/lon or easting/northing "
            f"columns among {columns}"
        )
    elevation = df[elev_col].to_numpy(dtype=float) if elev_col else None
    names = [str(v) for v in df[name_col].tolist()] if name_col else None
    return TopoTable(lon=lon, lat=lat, elevation=elevation, names=names, source=str(source))


def read_topo_file(
    path: str | Path,
    *,
    epsg: int | None = None,
    utm_zone: Any | None = None,
    latlon: bool = False,
) -> TopoTable:
    """Parse a ``.bln``/``.csv``/``.stn`` topo file into a
    :class:`TopoTable`.

    Parameters
    ----------
    path : path-like
        A ``.bln``, ``.csv``/``.txt``, or ``.stn`` file.
    epsg : int, optional
        EPSG code of the source's projected CRS, when it stores
        easting/northing rather than lat/lon. Takes precedence over
        *utm_zone* when both are given (matches
        :func:`pycsamt.gis.utils.to_ll`'s own precedence).
    utm_zone : optional
        UTM zone designator (e.g. ``"32N"``), an alternative to *epsg*.
    latlon : bool, default False
        ``.bln`` only: set ``True`` when the file's own ``x, y``
        columns are already ``lon, lat`` (a ``.bln`` carries no CRS
        metadata to detect this from). Ignored for ``.csv``/``.stn``,
        which are only treated as already-geographic when their own
        header says ``lat``/``lon``.

    Raises
    ------
    ValueError
        Unrecognised extension, a malformed ``.bln`` header/body, no
        recognisable coordinate columns in a ``.csv``/``.stn`` file,
        or projected coordinates with neither *epsg* nor *utm_zone*.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".bln":
        return _read_bln(path, epsg=epsg, utm_zone=utm_zone, latlon=latlon)
    if suffix == ".stn":
        from ..zonge.utils import read_stn

        df = read_stn(path)
        return _dataframe_to_topotable(df, str(path), epsg=epsg, utm_zone=utm_zone)
    if suffix in (".csv", ".txt"):
        import pandas as pd

        df = pd.read_csv(path, comment="#", skip_blank_lines=True)
        return _dataframe_to_topotable(df, str(path), epsg=epsg, utm_zone=utm_zone)
    raise ValueError(
        f"{path}: unrecognised topo file extension {suffix!r} -- expected "
        "one of '.bln', '.csv', '.txt', '.stn'."
    )


# ---------------------------------------------------------------------
# Sites / MapData extraction
# ---------------------------------------------------------------------


def topo_from_sites(source: Any) -> TopoTable:
    """Extract lon/lat/elevation from an already-geo-located
    ``Sites``/``MapData``/iterable-of-station-record object.

    Duck-typed on purpose: works with anything iterable whose items
    (or whose ``.stations`` attribute's items, for a ``MapData``-like
    container) expose a station identifier (``id``/``name``/
    ``station``) plus ``longitude``/``latitude`` and, optionally,
    ``elevation`` -- :class:`pycsamt.map._core.StationRecord` and
    similar objects all qualify without any adapter code.

    Raises
    ------
    ValueError
        If no station in *source* carries a usable id + lon/lat pair.
    """
    stations = getattr(source, "stations", source)
    names: list[str] = []
    lons: list[float] = []
    lats: list[float] = []
    elevs: list[float] = []
    for st in stations:
        name = getattr(st, "id", None) or getattr(st, "name", None) or getattr(st, "station", None)
        lon = getattr(st, "longitude", None)
        lat = getattr(st, "latitude", None)
        if name is None or lon is None or lat is None:
            continue
        elev = getattr(st, "elevation", None)
        names.append(str(name))
        lons.append(float(lon))
        lats.append(float(lat))
        elevs.append(float(elev) if elev is not None else np.nan)
    if not names:
        raise ValueError(
            "topo_from_sites: no geo-located station (id + longitude + "
            "latitude) found in the given source"
        )
    return TopoTable(
        lon=np.array(lons),
        lat=np.array(lats),
        elevation=np.array(elevs),
        names=names,
        source=f"{type(source).__name__} ({len(names)} geo-located station(s))",
    )


# ---------------------------------------------------------------------
# Coercion: any accepted "topo" argument -> TopoTable
# ---------------------------------------------------------------------


def _mapping_to_topotable(mapping: Mapping[str, Any], source: str) -> TopoTable:
    names = list(mapping.keys())
    values = list(mapping.values())
    lon = np.array([float(v[0]) for v in values], dtype=float)
    lat = np.array([float(v[1]) for v in values], dtype=float)
    elevation = np.array(
        [float(v[2]) if len(v) > 2 else np.nan for v in values], dtype=float
    )
    return TopoTable(lon=lon, lat=lat, elevation=elevation, names=names, source=source)


def _coerce_topo_table(
    source: Any, *, epsg: int | None, utm_zone: Any | None, latlon: bool
) -> TopoTable:
    if isinstance(source, TopoTable):
        return source
    if isinstance(source, (str, Path)):
        return read_topo_file(source, epsg=epsg, utm_zone=utm_zone, latlon=latlon)
    if isinstance(source, Mapping):
        return _mapping_to_topotable(source, source="in-memory {name: (lon, lat[, elev])} mapping")
    return topo_from_sites(source)


# ---------------------------------------------------------------------
# Attribution
# ---------------------------------------------------------------------


def attribute_topo(
    topo: TopoTable,
    station_names: Sequence[str],
    *,
    on_mismatch: str = "raise",
) -> TopoAttribution:
    """Match a parsed :class:`TopoTable` onto *station_names*.

    Name-based when ``topo.names`` is populated (exact match, then a
    normalized fallback); positional (in order, requiring an exact
    count match) otherwise. See the module docstring for the full
    rationale.

    Parameters
    ----------
    on_mismatch : {"raise", "warn"}, default "raise"
        Only consulted for a positional (name-less) source. ``"raise"``
        rejects a station-count mismatch outright; ``"warn"`` issues a
        :class:`UserWarning` and attributes only the overlapping
        prefix (``min(topo.n, len(station_names))`` points, in order).
    """
    if on_mismatch not in ("raise", "warn"):
        raise ValueError(f"on_mismatch must be 'raise' or 'warn', got {on_mismatch!r}")

    station_names = [str(n) for n in station_names]
    lon_by: dict[str, float] = {}
    lat_by: dict[str, float] = {}
    elev_by: dict[str, float] = {}
    matched: list[str] = []

    if topo.names is not None:
        from ..map._core import normalize_station_id

        index: dict[str, int] = {}
        for i, name in enumerate(topo.names):
            index.setdefault(name, i)
            index.setdefault(name.strip().lower(), i)
            index.setdefault(normalize_station_id(name), i)
        for name in station_names:
            i = (
                index.get(name)
                if index.get(name) is not None
                else index.get(name.strip().lower(), index.get(normalize_station_id(name)))
            )
            if i is None:
                continue
            lon_by[name] = float(topo.lon[i])
            lat_by[name] = float(topo.lat[i])
            if topo.elevation is not None and np.isfinite(topo.elevation[i]):
                elev_by[name] = float(topo.elevation[i])
            matched.append(name)
    else:
        if topo.n != len(station_names):
            msg = (
                f"{topo.source}: has {topo.n} point(s) but "
                f"{len(station_names)} station(s) were expected -- a "
                "name-less topo source must carry exactly one point per "
                "station, in survey order."
            )
            if on_mismatch == "raise":
                raise ValueError(msg)
            warnings.warn(
                msg + " Attributing only the overlapping prefix (on_mismatch='warn').",
                UserWarning,
                stacklevel=2,
            )
        n = min(topo.n, len(station_names))
        for i in range(n):
            name = station_names[i]
            lon_by[name] = float(topo.lon[i])
            lat_by[name] = float(topo.lat[i])
            if topo.elevation is not None and np.isfinite(topo.elevation[i]):
                elev_by[name] = float(topo.elevation[i])
            matched.append(name)

    unmatched = [n for n in station_names if n not in lon_by]
    return TopoAttribution(
        lon=lon_by,
        lat=lat_by,
        elevation=elev_by,
        matched=matched,
        unmatched_stations=unmatched,
        source=topo.source,
    )


def _merge_attributions(parts: Sequence[TopoAttribution]) -> TopoAttribution:
    merged = TopoAttribution(source="; ".join(p.source for p in parts) or "none")
    for part in parts:
        merged.lon.update(part.lon)
        merged.lat.update(part.lat)
        merged.elevation.update(part.elevation)
        merged.matched.extend(part.matched)
        merged.unmatched_stations.extend(part.unmatched_stations)
    return merged


def resolve_topo(
    topo: Any,
    station_names: Sequence[str] | Mapping[str, Sequence[str]],
    *,
    epsg: int | None = None,
    utm_zone: Any | None = None,
    latlon: bool = False,
    on_mismatch: str = "raise",
) -> TopoAttribution:
    """Resolve any accepted ``topo=`` argument into a
    :class:`TopoAttribution` against *station_names*.

    Parameters
    ----------
    topo : None, path-like, TopoTable, Sites/MapData-like, mapping, or sequence
        - ``None`` -- returns an empty attribution (every adapter's
          existing ``station_elevations``/``station_lonlat`` behaviour
          is unaffected).
        - a ``.bln``/``.csv``/``.stn`` path, or an already-parsed
          :class:`TopoTable`.
        - a ``Sites``/``MapData``/iterable-of-station-record object
          (see :func:`topo_from_sites`).
        - a plain ``{station_name: (lon, lat[, elevation])}`` mapping.
        - when *station_names* is a mapping (multiline, ``{line_id:
          [names, ...]}``): a ``{line_id: <any of the above>}``
          mapping, or a sequence with exactly one source per line, in
          the same order as *station_names*'s own keys. A single
          non-mapping, non-per-line-sequence source is instead matched
          by name across *all* lines' stations combined -- the natural
          choice for one combined ``.stn``/``.csv``/Sites source
          covering a whole multiline survey.
    station_names : sequence of str, or mapping of str to sequence of str
        The inversion's own station names (flat), or ``{line_id:
        names}`` for a multiline build.
    epsg, utm_zone, latlon, on_mismatch
        Forwarded to :func:`read_topo_file`/:func:`attribute_topo` for
        every file-based source encountered.

    Returns
    -------
    TopoAttribution
        Empty (all fields blank, ``source="none"``) when *topo* is
        ``None``.
    """
    if isinstance(station_names, Mapping):
        lines = list(station_names.items())
        if isinstance(topo, Mapping) and not isinstance(topo, TopoTable):
            parts = []
            for line_id, names in lines:
                src = topo.get(line_id)
                if src is None:
                    parts.append(
                        TopoAttribution(unmatched_stations=list(names), source="none")
                    )
                    continue
                table = _coerce_topo_table(src, epsg=epsg, utm_zone=utm_zone, latlon=latlon)
                parts.append(attribute_topo(table, names, on_mismatch=on_mismatch))
            return _merge_attributions(parts)

        if (
            isinstance(topo, Sequence)
            and not isinstance(topo, (str, Path))
            and len(lines) > 1
        ):
            if len(topo) != len(lines):
                raise ValueError(
                    f"got {len(topo)} topo source(s) but {len(lines)} "
                    f"line(s) were expected ({[lid for lid, _ in lines]}); "
                    "pass one topo source per line, in the same order, "
                    "or a single named source (.stn/.csv/dict/Sites) "
                    "matched by station name across every line."
                )
            parts = [
                attribute_topo(
                    _coerce_topo_table(src, epsg=epsg, utm_zone=utm_zone, latlon=latlon),
                    names,
                    on_mismatch=on_mismatch,
                )
                for (_, names), src in zip(lines, topo)
            ]
            return _merge_attributions(parts)

        all_names = [n for _, names in lines for n in names]
        if topo is None:
            return TopoAttribution(unmatched_stations=all_names, source="none")
        table = _coerce_topo_table(topo, epsg=epsg, utm_zone=utm_zone, latlon=latlon)
        return attribute_topo(table, all_names, on_mismatch=on_mismatch)

    if topo is None:
        return TopoAttribution(unmatched_stations=[str(n) for n in station_names], source="none")
    table = _coerce_topo_table(topo, epsg=epsg, utm_zone=utm_zone, latlon=latlon)
    return attribute_topo(table, station_names, on_mismatch=on_mismatch)
