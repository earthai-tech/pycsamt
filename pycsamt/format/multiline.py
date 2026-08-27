# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Multiline PCSF builder/reader — Phase 5 of the PCSF format plan.

Formalizes what ``pycsamt/app/web/callbacks/map3d.py`` currently
reconstructs at render time from a stack of independent 2-D sections:
``_line_real_offsets``/``_assemble_3d_grid`` derive real-vs-synthetic
line offsets and resample each line onto a common (x, z) grid, purely
in memory, from whatever cached ``InversionResult`` or pseudo-section
data the current session happens to hold. Nothing is ever persisted.

This module is the inverse pair:

- :func:`build_multiline_pcsf` turns a *profiles* dict — the same
  ``{line_name: {"x", "z", "rho", "sta_x", "sta_names", ...}}`` shape
  ``map3d.py``'s own ``_profiles_from_pseudo``/
  ``_profiles_from_inversion_result`` already produce — into a
  :class:`~pycsamt.format.schema.PCSFModel` with a ``multiline``
  geometry, one real (unresampled) :class:`~pycsamt.format.schema.LineEntry`
  per line.
- :func:`multiline_pcsf_to_profiles` reconstructs that same profiles
  dict shape from a loaded file, so it is a drop-in alternate data
  source for ``map3d.py``'s existing renderers
  (``_build_fence_fig``/``_build_block_fig``/``_assemble_3d_grid``) —
  no rendering code needs to change to consume a persisted file.
- :func:`line_offsets_from_stations` and
  :func:`stack_lines_to_common_grid` are the same real-offset /
  common-grid-resampling algorithms ``map3d.py`` already implements
  privately as ``_line_real_offsets``/``_assemble_3d_grid``, promoted
  here so both the live rendering path and the persisted-file path
  are provably the same computation, not two copies that can drift.
"""

from __future__ import annotations

import warnings
from typing import Any, Mapping, Sequence

import numpy as np

from ..map.geometry import normalize_offsets, resolve_offset, survey_uv
from .schema import (
    DerivedVolume,
    Grid2DGeometry,
    Grid3DGeometry,
    LineEntry,
    MultilineGeometry,
    PCSFModel,
    StationTable,
    TopographyPerStation,
)
from .topo_source import resolve_topo

__all__ = [
    "line_offsets_from_stations",
    "stack_lines_to_common_grid",
    "build_multiline_pcsf",
    "multiline_pcsf_to_profiles",
]


def _survey_to_dict(survey: Any | Mapping[str, Any] | None) -> dict[str, Any]:
    if survey is None:
        return {}
    to_dict = getattr(survey, "to_dict", None)
    if callable(to_dict):
        return dict(to_dict())
    return dict(survey)


def line_offsets_from_stations(
    profiles: Mapping[str, Mapping[str, Any]],
) -> dict[str, float] | None:
    """Cross-strike offset (m) for each line, from real station lat/lon.

    Same algorithm as ``map3d.py``'s private ``_line_real_offsets``
    (built on :func:`pycsamt.map.geometry.survey_uv`), so both the
    live-cache rendering path and this persisted-file path place lines
    identically. Requires every line to carry ``sta_lat``/``sta_lon``/
    ``sta_names`` of equal length; returns ``None`` otherwise so
    callers fall back to a synthetic index-based stack via
    :func:`~pycsamt.map.geometry.resolve_offset`.
    """
    ids: list[str] = []
    lats: list[float] = []
    lons: list[float] = []
    lines: list[str] = []
    per_line_ids: dict[str, list[str]] = {}
    for name, p in profiles.items():
        lat = p.get("sta_lat") or []
        lon = p.get("sta_lon") or []
        names = p.get("sta_names") or []
        if not lat or not lon or len(lat) != len(lon) or len(lat) != len(names):
            return None
        line_ids = [f"{name}::{n}" for n in names]
        per_line_ids[name] = line_ids
        ids.extend(line_ids)
        lats.extend(lat)
        lons.extend(lon)
        lines.extend([name] * len(line_ids))
    if len(ids) < 2:
        return None
    uv = survey_uv(ids, lats, lons, lines)
    if not uv:
        return None
    raw: dict[str, float] = {}
    for name, line_ids in per_line_ids.items():
        vs = [uv[i][1] for i in line_ids if i in uv]
        if not vs:
            return None
        raw[name] = float(np.median(vs))
    return normalize_offsets(raw)


def _resample_line_to_grid(
    x_i: np.ndarray,
    z_i: np.ndarray,
    rho_i_raw: np.ndarray,
    x_arr: np.ndarray,
    z_arr: np.ndarray,
) -> np.ndarray:
    """Resample one line's rho onto ``(x_arr, z_arr)``, returned (n_x, n_z).

    Mirrors ``_assemble_3d_grid``'s per-line branch exactly, including
    its two-step shape check (a cheap transpose first, a full
    ``RegularGridInterpolator`` resample only when that alone does not
    land on ``(n_x, n_z)``) so a cached :class:`DerivedVolume` is
    numerically identical to what live rendering would compute.
    """
    n_x, n_z = len(x_arr), len(z_arr)
    rho_i = np.asarray(rho_i_raw)
    if rho_i.shape == (n_z, n_x):
        rho_i = rho_i.T
    if rho_i.shape != (n_x, n_z):
        from scipy.interpolate import RegularGridInterpolator

        try:
            interp = RegularGridInterpolator(
                (np.asarray(z_i, float), np.asarray(x_i, float)),
                np.asarray(rho_i_raw),
                bounds_error=False,
                fill_value=np.nan,
            )
            xi = np.array(
                np.meshgrid(z_arr, x_arr, indexing="ij")
            ).T.reshape(-1, 2)
            rho_i = interp(xi).reshape(n_x, n_z)
        except Exception:
            rho_i = np.full((n_x, n_z), np.nanmean(rho_i))
    return rho_i


def stack_lines_to_common_grid(
    profiles: Mapping[str, Mapping[str, Any]],
    *,
    line_spacing: float = 1.0,
    fallback_unit: float = 1000.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Resample every line onto the first line's own (x, z) grid.

    Same algorithm as ``map3d.py``'s private ``_assemble_3d_grid``:
    the first profile's grid is the reference; other lines are
    resampled onto it via :func:`_resample_line_to_grid`.

    Returns
    -------
    x_arr, z_arr : ndarray
        The reference line's own coordinates.
    y_arr : ndarray, shape (n_lines,)
        Per-line cross-strike offset, from
        :func:`line_offsets_from_stations` when available, otherwise a
        synthetic ``idx * spacing * fallback_unit`` stack (see
        :func:`pycsamt.map.geometry.resolve_offset`).
    rho_3d : ndarray, shape (n_lines, n_x, n_z)
    """
    line_names = list(profiles.keys())
    if not line_names:
        raise ValueError("profiles must have at least one line")
    real_offsets = line_offsets_from_stations(profiles)
    ref = profiles[line_names[0]]
    x_arr = np.asarray(ref["x"], dtype=float)
    z_arr = np.asarray(ref["z"], dtype=float)
    n_x, n_z = len(x_arr), len(z_arr)

    rho_3d = np.zeros((len(line_names), n_x, n_z))
    for i, name in enumerate(line_names):
        p = profiles[name]
        rho_3d[i] = _resample_line_to_grid(
            p["x"], p["z"], p["rho"], x_arr, z_arr
        )

    y_arr = np.array(
        [
            resolve_offset(name, i, real_offsets, fallback_unit, line_spacing)
            for i, name in enumerate(line_names)
        ],
        dtype=float,
    )
    return x_arr, z_arr, y_arr, rho_3d


def build_multiline_pcsf(
    profiles: Mapping[str, Mapping[str, Any]],
    *,
    line_spacing: float = 1.0,
    fallback_unit: float = 1000.0,
    cache_derived_volume: bool = True,
    topo: Any = None,
    epsg: int | None = None,
    utm_zone: Any | None = None,
    latlon: bool = False,
    on_mismatch: str = "raise",
    survey: Any | Mapping[str, Any] | None = None,
    source_backend: str = "generic",
    created_by: str = "",
    crs: str | None = None,
    description: str = "",
) -> PCSFModel:
    r"""Build a ``multiline`` :class:`PCSFModel` from a profiles dict.

    Parameters
    ----------
    profiles : mapping of str to mapping
        ``{line_name: {"x": (n_x,), "z": (n_z,), "rho": (n_z, n_x), ...}}``,
        the exact shape ``map3d.py``'s ``_profiles_from_pseudo``/
        ``_profiles_from_inversion_result`` already produce. ``rho``
        must already be linear ohm.m (call
        ``_rho_log_to_ohm_m``-equivalent conversion first, matching
        every other PCSF adapter's canonical-linear convention).
        Optional per-line keys ``sta_x``, ``sta_names``, ``sta_elev``,
        ``sta_lat``, ``sta_lon`` populate :attr:`PCSFModel.stations`/
        :attr:`PCSFModel.topography` when present. ``sta_lat``/
        ``sta_lon`` do double duty: :func:`line_offsets_from_stations`
        uses them (via :func:`pycsamt.map.geometry.survey_uv`) to
        compute each line's *real* cross-strike offset when every line
        carries them, falling back to a synthetic index-based stack
        otherwise (see ``offset_kind`` below) -- and the same values
        are also persisted into :attr:`PCSFModel.stations`'
        ``lon``/``lat``, so a loaded multiline file is
        self-sufficiently geo-referenced too, not just correctly
        spaced.
    line_spacing, fallback_unit : float
        Forwarded to :func:`stack_lines_to_common_grid` for the
        optional cached :class:`DerivedVolume` (real per-line
        geometry itself never depends on these — only the synthetic
        offset fallback does).
    cache_derived_volume : bool, default True
        When ``True`` and there are at least two lines, also cache a
        :class:`DerivedVolume` (each line resampled onto a common
        grid) so a large multiline file does not need to re-resample
        on every render. Set ``False`` to keep the file smaller when
        that convenience volume is not needed.
    topo : optional
        A "smart" real-coordinate source resolved via
        :func:`pycsamt.format.topo_source.resolve_topo` -- see
        :func:`pycsamt.format.adapters.occam2d.occam2d_to_pcsf`'s
        identical parameter for the full description of accepted
        source types. Populates each line's own ``sta_lat``/``sta_lon``
        *before* the real-offset computation above runs, so passing
        *topo* is enough to get both a real cross-strike offset per
        line and a self-georeferenced file -- no separate offset step
        is needed. Accepts either a single source matched by station
        name across every line combined (a `.stn`/`.csv`/Sites source
        covering the whole survey), or a ``{line_name: source}``
        mapping / one-source-per-line sequence (in ``profiles``'s own
        key order) for e.g. one ``.bln`` file per surveyed line. Takes
        precedence over any ``sta_lat``/``sta_lon`` already present in
        *profiles* for every station it resolves (with a
        :class:`UserWarning` if both were supplied). If *topo* only
        partially covers a line's stations, the resulting ``nan``
        entries make :func:`line_offsets_from_stations` fall back to a
        *synthetic* offset for every line (mixing a real and a
        synthetic stack would look inconsistent -- see that function's
        own all-or-nothing behaviour) rather than silently using a
        partially-real one; the stations themselves still keep
        whatever real lon/lat *topo* did resolve. ``None`` (the
        default) leaves this function's behaviour exactly as it was
        before *topo* existed.
    epsg, utm_zone, latlon, on_mismatch
        Forwarded to :func:`pycsamt.format.topo_source.resolve_topo`;
        see *occam2d_to_pcsf*'s identical parameters.
    survey, source_backend, created_by, crs, description :
        Passed straight through to :class:`PCSFModel`.

    Raises
    ------
    ValueError
        If *profiles* is empty, or any line is missing ``x``/``z``/``rho``.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.format.multiline import build_multiline_pcsf
    >>> profiles = {
    ...     "L1": {"x": np.array([0.0, 100.0]), "z": np.array([10.0, 50.0]),
    ...            "rho": np.array([[100.0, 110.0], [50.0, 55.0]])},
    ...     "L2": {"x": np.array([0.0, 100.0]), "z": np.array([10.0, 50.0]),
    ...            "rho": np.array([[200.0, 210.0], [90.0, 95.0]])},
    ... }
    >>> model = build_multiline_pcsf(profiles, cache_derived_volume=False)
    >>> model.kind
    'multiline'
    >>> [line.line_id for line in model.geometry.lines]
    ['L1', 'L2']
    """
    if not profiles:
        raise ValueError("profiles must have at least one line")
    line_names = list(profiles.keys())
    for name in line_names:
        missing = {"x", "z", "rho"} - set(profiles[name])
        if missing:
            raise ValueError(f"line {name!r} is missing required key(s) {missing}")

    if topo is not None:
        if any(
            profiles[name].get("sta_lat") or profiles[name].get("sta_lon")
            for name in line_names
        ):
            warnings.warn(
                "build_multiline_pcsf: 'topo' takes precedence over any "
                "'sta_lat'/'sta_lon' already present in profiles for "
                "every station it resolves.",
                UserWarning,
                stacklevel=2,
            )
        station_names_by_line = {
            name: list(profiles[name].get("sta_names") or []) for name in line_names
        }
        attr = resolve_topo(
            topo,
            station_names_by_line,
            epsg=epsg,
            utm_zone=utm_zone,
            latlon=latlon,
            on_mismatch=on_mismatch,
        )
        patched: dict[str, Mapping[str, Any]] = {}
        for name in line_names:
            p = dict(profiles[name])
            sta_names = list(p.get("sta_names") or [])
            if sta_names:
                old_lat = list(p.get("sta_lat") or [])
                old_lon = list(p.get("sta_lon") or [])
                old_elev = list(p.get("sta_elev") or [])
                p["sta_lat"] = [
                    attr.lat.get(n, old_lat[i] if i < len(old_lat) else float("nan"))
                    for i, n in enumerate(sta_names)
                ]
                p["sta_lon"] = [
                    attr.lon.get(n, old_lon[i] if i < len(old_lon) else float("nan"))
                    for i, n in enumerate(sta_names)
                ]
                p["sta_elev"] = [
                    attr.elevation.get(n, old_elev[i] if i < len(old_elev) else float("nan"))
                    for i, n in enumerate(sta_names)
                ]
            patched[name] = p
        profiles = patched

    real_offsets = line_offsets_from_stations(profiles)
    offset_kind = "real" if real_offsets is not None else "synthetic"

    lines: list[LineEntry] = []
    station_names: list[str] = []
    station_x: list[float] = []
    station_y: list[float] = []
    station_z: list[float] = []
    station_line_id: list[str] = []
    station_lon: list[float] = []
    station_lat: list[float] = []
    topo_ids: list[str] = []
    topo_elev: list[float] = []

    for i, name in enumerate(line_names):
        p = profiles[name]
        geometry = Grid2DGeometry(x=p["x"], z=p["z"])
        offset_y = resolve_offset(
            name, i, real_offsets, fallback_unit, line_spacing
        )
        lines.append(
            LineEntry(
                line_id=str(name),
                geometry=geometry,
                resistivity=p["rho"],
                offset_y=offset_y,
                offset_kind=offset_kind,
            )
        )

        sta_names = list(p.get("sta_names") or [])
        sta_x = list(p.get("sta_x") or [])
        sta_elev = list(p.get("sta_elev") or [])
        sta_lat = list(p.get("sta_lat") or [])
        sta_lon = list(p.get("sta_lon") or [])
        if sta_names and len(sta_x) == len(sta_names):
            station_names.extend(sta_names)
            station_x.extend(float(v) for v in sta_x)
            station_y.extend([offset_y] * len(sta_names))
            if len(sta_elev) == len(sta_names):
                station_z.extend(float(v) for v in sta_elev)
            else:
                station_z.extend([float("nan")] * len(sta_names))
            station_line_id.extend([str(name)] * len(sta_names))
            if len(sta_elev) == len(sta_names):
                topo_ids.extend(sta_names)
                topo_elev.extend(float(v) for v in sta_elev)
            # Same real lat/lon already used above to compute this
            # line's cross-strike offset_y (via line_offsets_from_stations
            # / survey_uv) -- persisted into StationTable too, not just
            # consumed transiently for the offset math, so a loaded
            # multiline PCSF file stays self-sufficiently geo-referenced
            # (multiline_pcsf_to_profiles reads it straight back out).
            if len(sta_lat) == len(sta_names) and len(sta_lon) == len(sta_names):
                station_lat.extend(float(v) for v in sta_lat)
                station_lon.extend(float(v) for v in sta_lon)
            else:
                station_lat.extend([float("nan")] * len(sta_names))
                station_lon.extend([float("nan")] * len(sta_names))

    lon_arr = np.asarray(station_lon, dtype=float)
    lat_arr = np.asarray(station_lat, dtype=float)
    has_lonlat = lon_arr.size and not np.all(np.isnan(lon_arr))
    stations = (
        StationTable(
            name=station_names,
            x=np.asarray(station_x, dtype=float),
            y=np.asarray(station_y, dtype=float),
            z=np.asarray(station_z, dtype=float),
            line_id=station_line_id,
            lon=lon_arr if has_lonlat else None,
            lat=lat_arr if has_lonlat else None,
        )
        if station_names
        else None
    )
    topography = (
        TopographyPerStation(
            station_id=topo_ids, elevation=np.asarray(topo_elev, dtype=float)
        )
        if topo_ids
        else None
    )

    derived_volume = None
    if cache_derived_volume and len(line_names) >= 2:
        x_arr, z_arr, y_arr, rho_3d = stack_lines_to_common_grid(
            profiles, line_spacing=line_spacing, fallback_unit=fallback_unit
        )
        grid = Grid3DGeometry(x=x_arr, y=y_arr, z=z_arr)
        # rho_3d is (n_lines, n_x, n_z); Grid3DGeometry's canonical
        # order is (n_z, n_y, n_x) (matches ModEM's own convention).
        resistivity = np.transpose(rho_3d, (2, 0, 1))
        derived_volume = DerivedVolume(
            grid=grid,
            resistivity=resistivity,
            derivation_method="linear_interp",
            derived_from=list(line_names),
            synthesized=True,
        )

    geometry = MultilineGeometry(lines=lines, derived_volume=derived_volume)
    return PCSFModel(
        geometry=geometry,
        stations=stations,
        topography=topography,
        survey=_survey_to_dict(survey),
        source_backend=source_backend,
        created_by=created_by,
        crs=crs,
        description=description,
    )


def multiline_pcsf_to_profiles(model: PCSFModel) -> dict[str, dict[str, Any]]:
    """Reconstruct a profiles dict from a ``multiline`` :class:`PCSFModel`.

    The exact inverse of :func:`build_multiline_pcsf`'s per-line
    conversion — the result is a drop-in alternate data source for
    ``map3d.py``'s existing ``_build_fence_fig``/``_assemble_3d_grid``/
    ``_build_block_fig``, which only need the ``{"x", "z", "rho", ...}``
    shape, not any particular origin.

    Parameters
    ----------
    model : PCSFModel
        Must have ``geometry.kind == "multiline"``.

    Returns
    -------
    dict
        ``{line_id: {"x", "z", "rho", "sta_x", "sta_names", "sta_elev",
        "sta_lat", "sta_lon"}}``. Station keys are populated from
        :attr:`PCSFModel.stations` when present (matched to each line
        via ``StationTable.line_id``); ``sta_lat``/``sta_lon`` come
        from ``StationTable.lat``/``.lon`` when the file has them
        (see :func:`build_multiline_pcsf`), else stay empty lists --
        an older file written before those fields existed round-trips
        the same way it always did.

    Raises
    ------
    ValueError
        If *model* is not a ``multiline`` geometry.
    """
    if model.kind != "multiline":
        raise ValueError(
            f"multiline_pcsf_to_profiles needs a multiline PCSFModel, "
            f"got kind={model.kind!r}"
        )

    stations_by_line: dict[str, list[int]] = {}
    if model.stations is not None and model.stations.line_id is not None:
        for idx, line_id in enumerate(model.stations.line_id):
            stations_by_line.setdefault(line_id, []).append(idx)

    profiles: dict[str, dict[str, Any]] = {}
    for line in model.geometry.lines:
        entry: dict[str, Any] = {
            "x": np.asarray(line.geometry.x),
            "z": np.asarray(line.geometry.z),
            "rho": np.asarray(line.resistivity),
        }
        idxs = stations_by_line.get(line.line_id, [])
        if idxs and model.stations is not None:
            entry["sta_x"] = [float(model.stations.x[i]) for i in idxs]
            entry["sta_names"] = [str(model.stations.name[i]) for i in idxs]
            entry["sta_elev"] = [float(model.stations.z[i]) for i in idxs]
            if model.stations.lat is not None and model.stations.lon is not None:
                entry["sta_lat"] = [float(model.stations.lat[i]) for i in idxs]
                entry["sta_lon"] = [float(model.stations.lon[i]) for i in idxs]
            else:
                entry["sta_lat"] = []
                entry["sta_lon"] = []
        profiles[line.line_id] = entry
    return profiles
