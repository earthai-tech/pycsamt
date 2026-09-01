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

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np

from ._core import (
    MapData,
    StationRecord,
    normalize_station_id,
)

__all__ = ["group_modem_stations", "load_modem_lines", "load_pcsf_lines"]

#: Cell resistivity (ohm.m) above which a ``grid3d`` PCSF cell is
#: treated as ModEM's above-topography "air" / padding fill and dropped
#: from a sliced curtain. A defensive net for ``.pcsf`` files written
#: before :func:`pycsamt.format.adapters.modem3d.modem3d_to_pcsf` began
#: masking that fill itself (mirrors its ``DEFAULT_AIR_THRESHOLD_OHM_M``);
#: harmless on newer files, whose air cells are already ``nan``.
_AIR_FILL_THRESHOLD = 1e8

#: A ModEM 3-D grid grades geometrically to tens of km for the boundary
#: conditions; those deep cells carry no interpretable resolution and,
#: left in a sliced curtain, dominate the vertical extent of every 3-D
#: view. Trim trailing cells once a cell is this many times thicker than
#: the shallowest earth cell (keeping at least ``_MIN_EARTH_CELLS``).
_BC_PADDING_GROWTH = 15.0
_MIN_EARTH_CELLS = 8


def _earth_cell_cut(z_nodes_earth: np.ndarray) -> int:
    """Index of the first boundary-condition padding cell in an earth
    depth-node array (``0`` at the earth top). Cells at and beyond it are
    dropped from a sliced curtain."""
    dz = np.diff(np.asarray(z_nodes_earth, dtype=float))
    if dz.size <= _MIN_EARTH_CELLS:
        return dz.size
    too_thick = dz > _BC_PADDING_GROWTH * dz[0]
    if not too_thick.any():
        return dz.size
    return max(_MIN_EARTH_CELLS, int(np.argmax(too_thick)))


def _model_top_datum(model: Any, elevs: np.ndarray) -> float | None:
    """Elevation (m a.s.l.) of a ``grid3d`` model's own ``z = 0`` plane.

    ModEM's model z-axis is measured downward from the *model top*, which
    for a run that represents topography through a high-resistivity fill
    (rather than explicit air layers) sits well above the real ground --
    so a curtain sliced straight off ``geometry.z`` is referenced to the
    model top, not each station's own surface, and a "200 m" depth slice
    is really ~(200 - (model_top - station_elev)) m below ground.

    :func:`pycsamt.format.adapters.modem3d.modem3d_to_pcsf` records the
    datum in ``metadata["station_z"]["datum_masl"]`` when the source
    ``.dat`` file names one (``... m a.s.l.``). Absent that, the highest
    station elevation is a safe floor: the model must extend at least to
    the highest ground. ``None`` when there is no topography to reference
    at all (flat/zero station z).
    """
    meta = (getattr(model, "metadata", None) or {}).get("station_z") or {}
    datum = meta.get("datum_masl")
    if datum is not None:
        return float(datum)
    finite = elevs[np.isfinite(elevs)]
    if finite.size == 0 or not np.any(finite != 0.0):
        return None
    return float(np.max(finite))


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
        line = (
            str(match.line)
            if match is not None and match.line
            else _line_token(name)
        )
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


def _grid3d_offset(
    node_extent: tuple[float, float] | None,
    centers: np.ndarray,
    coords: np.ndarray,
) -> float:
    """Registration offset between a PCSF ``grid3d`` file's own grid
    frame and its station coordinate frame.

    Same symmetric-padding heuristic
    :func:`pycsamt.models.modem.section.station_curtain` uses for a
    *live* ModEM folder (there, ``_grid_offset``), reimplemented here
    from a PCSF file's plain arrays instead of a live
    ``ModEmModel3D``/``ModEmData`` pair: ModEM builds its horizontal
    grid with symmetric padding around the station footprint, so the
    grid's own extent midpoint always coincides with the midpoint of
    the station coordinate range -- a closed-form offset that needs no
    parsed padding/core-zone boundaries either way. Grid rotation is
    not applied here, matching ``station_curtain``'s own scope -- both
    are exact for an axis-aligned grid (``rotation_deg == 0``, the
    common case) and best-effort otherwise.
    """
    finite = coords[np.isfinite(coords)]
    if finite.size == 0 or centers.size == 0:
        return 0.0
    data_mid = (float(finite.min()) + float(finite.max())) / 2.0
    if node_extent is not None:
        grid_mid = (node_extent[0] + node_extent[1]) / 2.0
    else:
        grid_mid = (float(centers[0]) + float(centers[-1])) / 2.0
    return grid_mid - data_mid


def _nearest_cell(centers: np.ndarray, value: float) -> int:
    idx = int(np.searchsorted(centers, value))
    if idx <= 0:
        return 0
    if idx >= centers.size:
        return centers.size - 1
    return (
        idx
        if abs(centers[idx] - value) < abs(centers[idx - 1] - value)
        else idx - 1
    )


def _grid3d_sections(
    model: Any,
    known_stations: Iterable[StationRecord] | None,
) -> tuple[list[StationRecord], dict[str, dict[str, Any]]]:
    """Slice one 2-D vertical curtain per survey line from a PCSF
    ``grid3d`` volume.

    The same nearest-cell sampling
    :func:`pycsamt.models.modem.section.station_curtain` already does
    for a *live* ModEM folder — a single ModEM 3-D run inverts one
    volume for a whole (possibly multi-line) survey, so a curtain has
    to be sliced out per line, not read off a per-line file. Lines are
    grouped the same way :func:`load_modem_lines` groups a live
    folder's stations (:func:`group_modem_stations`, by
    ``known_stations`` line tags first, else a station-name-prefix
    heuristic) whenever the file carries no explicit
    ``stations.line_id`` --
    :func:`pycsamt.format.adapters.modem3d.modem3d_to_pcsf` does not
    currently set one.
    """
    geo = model.geometry
    st = model.stations
    known_by_id = _index_known_stations(known_stations)

    if st.line_id:
        groups: dict[str, list[str]] = {}
        for name, line in zip(st.name, st.line_id):
            groups.setdefault(str(line), []).append(str(name))
    else:
        groups = group_modem_stations(st.name, known_stations=known_stations)

    name_to_idx = {str(name): i for i, name in enumerate(st.name)}
    x_nodes = geo.x_nodes if geo.x_nodes is not None else geo.x
    y_nodes = geo.y_nodes if geo.y_nodes is not None else geo.y
    offset_x = _grid3d_offset(
        (float(x_nodes[0]), float(x_nodes[-1])), geo.x, st.x
    )
    offset_y = _grid3d_offset(
        (float(y_nodes[0]), float(y_nodes[-1])), geo.y, st.y
    )

    n_air = int(getattr(geo, "n_air", 0) or 0)
    # Re-zero the earth's depth axis so 0 is the *top of the earth
    # domain*, not the model top (matches
    # :func:`pycsamt.models.modem.section.station_curtain`); the two
    # differ by the air-layer thickness for an ``n_air > 0`` model.
    cut = None
    if geo.z_nodes is not None and len(geo.z_nodes) >= n_air + 2:
        z_nodes_earth = np.asarray(geo.z_nodes[n_air:], dtype=float)
        z_nodes_earth = z_nodes_earth - z_nodes_earth[0]
        cut = _earth_cell_cut(z_nodes_earth)
        z_earth = (z_nodes_earth[:-1] + z_nodes_earth[1:]) / 2.0
    else:
        z_earth = np.asarray(geo.z[n_air:], dtype=float)
        z_earth = z_earth - float(z_earth[0]) if z_earth.size else z_earth
    rho_earth = np.asarray(model.resistivity[n_air:, :, :], dtype=float)
    if cut is not None and cut < z_earth.size:
        # Drop the deep boundary-condition padding (no interpretable
        # resolution; only bloats every 3-D view's vertical extent).
        z_earth = z_earth[:cut]
        rho_earth = rho_earth[:cut, :, :]
    # Defensive: drop any residual above-topography air fill (very high
    # ohm.m) a pre-mask .pcsf file still carries -- newer files already
    # have nan here, so this is a no-op for them.
    rho_earth = np.where(rho_earth > _AIR_FILL_THRESHOLD, np.nan, rho_earth)

    # An ``n_air == 0`` ModEM model represents topography through a
    # high-resistivity fill, not explicit air layers, so its ``z = 0`` is
    # the flat model top and a station's real ground surface sits
    # ``datum - elev`` below it (see :func:`_model_top_datum`).
    # Re-reference every column to that station's own surface so the
    # returned ``z`` axis is a true depth-below-surface -- otherwise a
    # depth slice reads ~100 m shallower than its label and warps with
    # topography. An ``n_air > 0`` model already carries topography in its
    # own air layers, so the earth-top-relative axis above is already
    # correct and no per-station shift is applied.
    datum = None
    if n_air == 0:
        datum = _model_top_datum(
            model, np.asarray(getattr(st, "z", []), dtype=float)
        )

    stations: list[StationRecord] = []
    sections: dict[str, dict[str, Any]] = {}
    for line_id, names in groups.items():
        line_names: list[str] = []
        elevs: list[float] = []
        columns: list[np.ndarray] = []
        for name in names:
            i = name_to_idx.get(str(name))
            if i is None:
                continue
            match = _lookup(known_by_id, name)
            ix = _nearest_cell(geo.x, float(st.x[i]) + offset_x)
            iy = _nearest_cell(geo.y, float(st.y[i]) + offset_y)
            elev = float(st.z[i])
            if not np.isfinite(elev) and match is not None:
                elev = (
                    match.elevation
                    if match.elevation is not None
                    else np.nan
                )
            lon, lat = _resolve_pcsf_lonlat(st, i, match)
            stations.append(
                StationRecord(
                    id=str(name),
                    latitude=lat,
                    longitude=lon,
                    elevation=elev if np.isfinite(elev) else None,
                    line=str(line_id),
                    index=len(stations),
                )
            )
            line_names.append(str(name))
            elevs.append(elev)
            columns.append(
                _rezero_column(rho_earth[:, iy, ix], z_earth, datum, elev)
            )
        if not line_names:
            continue
        sections[str(line_id)] = {
            "z": z_earth,
            "rho": np.column_stack(columns),
            "stations": np.array(line_names, dtype=object),
            "elev": np.array(elevs, dtype=float),
        }
    return stations, sections


def _rezero_column(
    col: np.ndarray,
    z_earth: np.ndarray,
    datum: float | None,
    elev: float,
) -> np.ndarray:
    """Shift one resistivity column from model-top-referenced depth to
    depth below *elev*'s own ground surface, resampled onto *z_earth*.

    Interpolation is in ``log10`` (resistivity spans decades). Cells the
    shift would pull from above the model top, or from deeper than the
    model reaches, come back ``nan``. A no-op when there is no datum or
    no real elevation.
    """
    col = np.asarray(col, dtype=float)
    if datum is None or not np.isfinite(elev):
        return col
    surf_below_top = float(datum) - float(elev)
    if abs(surf_below_top) < 1e-6:
        return col
    # A cell at ``z_earth[k]`` below the earth top sits
    # ``z_earth[k] - surf_below_top`` below *this* station's surface.
    src_depth = z_earth - surf_below_top
    good = np.isfinite(col) & (col > 0) & np.isfinite(src_depth)
    if good.sum() < 2:
        return np.full(z_earth.shape, np.nan)
    log_out = np.interp(
        z_earth,
        src_depth[good],
        np.log10(col[good]),
        left=np.nan,
        right=np.nan,
    )
    return np.power(10.0, log_out)


_MESH_DEFAULT_N_Z = 60


def _mesh_unstructured_sections(
    model: Any,
    known_stations: Iterable[StationRecord] | None,
    n_z: int,
) -> tuple[list[StationRecord], dict[str, dict[str, Any]]]:
    """Slice one per-station curtain from a PCSF ``mesh_unstructured``
    (MARE2DEM) mesh via point-location on its real triangulation.

    Unlike ``grid2d``/``multiline``/``grid3d``, a mesh has no
    rectilinear index to look up a station's nearest column by — a
    triangle's own connectivity is the only structure there is, so a
    curtain needs an actual point-in-triangle query at each
    ``(x_station, z_sample)`` location. :class:`matplotlib.tri.Triangulation`
    (already a hard pycsamt dependency, no new install) builds that
    query structure directly from the mesh's own real ``nodes``/
    ``connectivity`` — no re-triangulation, no approximation of the
    mesh's own geometry.

    ``z_sample`` itself has no natural resolution to reuse (a mesh
    carries no separate z axis the way a rectilinear grid does), so
    *n_z* evenly spaced samples across the mesh's own node z-range are
    used — the same kind of pragmatic, documented default the rest of
    this format already makes rather than inventing an unbounded
    resolution choice (cf. PCSM's own row-width cap).

    A query point outside the mesh (above the mesh's shallowest
    triangle at that x, or beyond its lateral/depth extent) returns
    ``nan`` for that depth, never a fabricated value.
    """
    from matplotlib.tri import Triangulation

    geo = model.geometry
    st = model.stations
    known_by_id = _index_known_stations(known_stations)

    connectivity = np.asarray(geo.connectivity, dtype=np.int64)
    nodes = np.asarray(geo.nodes, dtype=float)
    resistivity = np.asarray(model.resistivity, dtype=float)
    if resistivity.shape != (connectivity.shape[0],):
        raise ValueError(
            "load_pcsf_lines needs per-triangle mesh_unstructured "
            f"resistivity (shape ({connectivity.shape[0]},)), got "
            f"shape {resistivity.shape} -- expand resistivity_by_region "
            "onto each triangle's region id first (every "
            "pycsamt.format.adapters writer already does this)."
        )

    triangulation = Triangulation(
        nodes[:, 0], nodes[:, 1], triangles=connectivity
    )
    trifinder = triangulation.get_trifinder()
    z_samples = np.linspace(
        float(nodes[:, 1].min()), float(nodes[:, 1].max()), n_z
    )

    if st.line_id:
        groups: dict[str, list[str]] = {}
        for name, line in zip(st.name, st.line_id):
            groups.setdefault(str(line), []).append(str(name))
    else:
        groups = {"line1": [str(name) for name in st.name]}
    name_to_idx = {str(name): i for i, name in enumerate(st.name)}

    stations: list[StationRecord] = []
    sections: dict[str, dict[str, Any]] = {}
    for line_id, names in groups.items():
        line_names: list[str] = []
        elevs: list[float] = []
        columns: list[np.ndarray] = []
        for name in names:
            i = name_to_idx.get(str(name))
            if i is None:
                continue
            match = _lookup(known_by_id, name)
            x_station = float(st.x[i])
            tri_idx = trifinder(np.full(n_z, x_station), z_samples)
            column = np.full(n_z, np.nan)
            inside = tri_idx >= 0
            column[inside] = resistivity[tri_idx[inside]]
            elev = float(st.z[i])
            if not np.isfinite(elev) and match is not None:
                elev = (
                    match.elevation
                    if match.elevation is not None
                    else np.nan
                )
            lon, lat = _resolve_pcsf_lonlat(st, i, match)
            stations.append(
                StationRecord(
                    id=str(name),
                    latitude=lat,
                    longitude=lon,
                    elevation=elev if np.isfinite(elev) else None,
                    line=str(line_id),
                    index=len(stations),
                )
            )
            line_names.append(str(name))
            elevs.append(elev)
            columns.append(column)
        if not line_names:
            continue
        sections[str(line_id)] = {
            "z": z_samples,
            "rho": np.column_stack(columns),
            "stations": np.array(line_names, dtype=object),
            "elev": np.array(elevs, dtype=float),
        }
    return stations, sections


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
    groups = group_modem_stations(
        data.site_names, known_stations=known_stations
    )

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


def load_pcsf_lines(
    path: str | Path,
    *,
    known_stations: Iterable[StationRecord] | None = None,
    fetch_elevation: bool = True,
    verbose: int = 0,
    mesh_z_samples: int = _MESH_DEFAULT_N_Z,
) -> MapData:
    """Load a backend-neutral ``.pcsf``/``.pcsm``/``.pcsm.gz`` file as a
    multi-line MapData.

    Phase 7 of the PCSF format plan: the same route
    :func:`load_modem_lines` already established, but working from any
    backend's converted ``grid2d``/``multiline`` PCSF file instead of a
    live ModEM folder — the view layer here carries no
    Occam2D/ModEM/MARE2DEM-specific logic, only PCSF's own schema.
    Both encodings of the format are accepted transparently — see
    :func:`pycsamt.format.read_pcsf_or_pcsm` — since PCSM is a
    lossless ASCII projection of the exact same in-memory model.

    ``pycsamt.map.volume``'s 3-D builders need one real ``(x, z, rho)``
    curtain *per named station*. For ``grid2d``/``multiline``, a PCSF
    file's own dense mesh column at each real station's along-profile
    position is extracted (the nearest-column convention
    :meth:`pycsamt.interp._base.ResistivityModel.column_nearest` also
    uses) — the full mesh resolution between stations is not carried
    into ``MapData``, the same simplification :func:`load_modem_lines`
    already makes via :func:`pycsamt.models.modem.section.station_curtain`.
    For ``grid3d`` (native ModEM 3-D), the same nearest-cell sampling
    ``station_curtain`` uses for a *live* ModEM folder is applied here
    to the PCSF file's own volume instead — see :func:`_grid3d_sections`.
    For ``mesh_unstructured`` (MARE2DEM), a station's column instead
    comes from point-location on the mesh's real triangulation (which
    triangle contains a given ``(x, z)`` query point) — a genuinely
    different algorithm from the other three kinds' index lookups,
    since a triangular mesh carries no rectilinear index to begin with
    — see :func:`_mesh_unstructured_sections`.

    Parameters
    ----------
    path : path-like
        A ``.pcsf``, ``.pcsm``, or ``.pcsm.gz`` file with
        ``geometry.kind`` ``"grid2d"``, ``"multiline"``, ``"grid3d"``,
        or ``"mesh_unstructured"``, and a populated
        :attr:`PCSFModel.stations` table.
    known_stations : iterable of StationRecord, optional
        Previously-loaded EDI stations, used the same way
        :func:`load_modem_lines` uses them: when a station matches one
        here, its ``line``/elevation/lon/lat take priority. Real
        geo-referencing does not strictly require this, though --
        ``StationTable``'s own ``lon``/``lat`` (populated by
        ``occam2d_to_pcsf(station_lonlat=...)``,
        ``modem3d_to_pcsf``'s ``GG_Lat``/``GG_Lon`` passthrough, or
        ``build_multiline_pcsf``'s ``sta_lat``/``sta_lon``), when
        present, are used as the fallback source -- only an
        unmatched station in a file with neither falls back to
        synthetic line spacing.
    fetch_elevation : bool, default True
        Best-effort online elevation lookup for any station still
        missing one after ``known_stations`` matching (see
        :func:`load_modem_lines`).
    verbose : int, default 0
        Unused today; accepted for signature parity with
        :func:`load_modem_lines`.
    mesh_z_samples : int, default 60
        ``mesh_unstructured`` only: number of evenly spaced depth
        samples across the mesh's own node z-range to query per
        station. A triangular mesh carries no separate z axis to reuse
        the way a rectilinear grid does, so this is a pragmatic,
        overridable resolution choice, not a property of the file
        itself. Ignored for every other geometry kind.

    Returns
    -------
    MapData
        ``sites=None``; ``metadata["sections"]`` carries one precomputed
        curtain per line, consumed directly by :mod:`pycsamt.map.volume`.

    Raises
    ------
    ValueError
        If the file has no station table, or (``mesh_unstructured``
        only) its resistivity is the compact per-region form rather
        than already expanded onto each triangle.
    """
    from pycsamt.format import read_pcsf_or_pcsm

    model = read_pcsf_or_pcsm(path)
    if model.stations is None or not list(model.stations.name):
        raise ValueError(
            f"{path}: PCSF file has no station table -- MapView needs "
            "named stations to build 2-D curtains."
        )

    if model.kind == "grid3d":
        stations, sections = _grid3d_sections(model, known_stations)
    elif model.kind == "mesh_unstructured":
        stations, sections = _mesh_unstructured_sections(
            model, known_stations, mesh_z_samples
        )
    else:
        if model.kind == "grid2d":
            line_specs = [
                (
                    "line1",
                    model.geometry.x,
                    model.geometry.z,
                    model.resistivity,
                )
            ]
        else:
            line_specs = [
                (
                    line.line_id,
                    line.geometry.x,
                    line.geometry.z,
                    line.resistivity,
                )
                for line in model.geometry.lines
            ]

        st = model.stations
        station_by_line: dict[str, list[int]] = {}
        for i in range(len(st.name)):
            line_id = st.line_id[i] if st.line_id else "line1"
            station_by_line.setdefault(str(line_id), []).append(i)

        known_by_id = _index_known_stations(known_stations)
        stations = []
        sections = {}
        for line_id, x_geo, z_geo, rho in line_specs:
            idxs = station_by_line.get(str(line_id), [])
            if not idxs:
                continue
            names: list[str] = []
            elevs: list[float] = []
            columns: list[np.ndarray] = []
            for i in idxs:
                name = str(st.name[i])
                match = _lookup(known_by_id, name)
                col_idx = int(np.argmin(np.abs(x_geo - st.x[i])))
                elev = float(st.z[i])
                if not np.isfinite(elev) and match is not None:
                    elev = (
                        match.elevation
                        if match.elevation is not None
                        else np.nan
                    )
                lon, lat = _resolve_pcsf_lonlat(st, i, match)
                stations.append(
                    StationRecord(
                        id=name,
                        latitude=lat,
                        longitude=lon,
                        elevation=elev if np.isfinite(elev) else None,
                        line=str(line_id),
                        index=len(stations),
                    )
                )
                names.append(name)
                elevs.append(elev)
                columns.append(rho[:, col_idx])
            sections[str(line_id)] = {
                "z": np.asarray(z_geo, dtype=float),
                "rho": np.column_stack(columns),
                "stations": np.array(names, dtype=object),
                "elev": np.array(elevs, dtype=float),
            }

    metadata = {
        "source": "pcsf",
        "path": str(path),
        "sections": sections,
        "rms": (model.metadata or {}).get("final_rms"),
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


def _resolve_pcsf_lonlat(
    st: Any,
    i: int,
    known_match: StationRecord | None,
) -> tuple[float | None, float | None]:
    """Same priority :func:`_resolve_lonlat` uses for a live ModEM
    folder's ``site_lonlat`` -- ``known_stations`` first (lets an
    already-geo-located EDI survey override/enrich), else the PCSF
    file's own ``StationTable.lon``/``.lat`` (populated by
    ``occam2d_to_pcsf``'s ``station_lonlat``,
    ``modem3d_to_pcsf``'s own ``GG_Lat``/``GG_Lon`` passthrough, or
    ``build_multiline_pcsf``'s ``sta_lat``/``sta_lon``), else unknown.
    """
    if (
        known_match is not None
        and known_match.longitude is not None
        and known_match.latitude is not None
    ):
        return known_match.longitude, known_match.latitude
    if st.lon is not None and st.lat is not None:
        lon, lat = float(st.lon[i]), float(st.lat[i])
        if np.isfinite(lon) and np.isfinite(lat):
            return lon, lat
    return None, None
