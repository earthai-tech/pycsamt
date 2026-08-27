# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Occam2D -> PCSF adapter (Phase 2 of the PCSF format plan).

Converts a completed :class:`pycsamt.models.occam2d.results.InversionResult`
into a backend-neutral :class:`~pycsamt.format.schema.PCSFModel` with a
``grid2d`` geometry. Built directly on
:meth:`pycsamt.interp._base.ResistivityModel.from_occam2d`, which already
recovers real station chainage from Occam2D's mesh-local coordinate frame
(see that method's docstring for the padding-symmetry argument) — this
module reuses that correction rather than re-deriving it, and only adds
the cell-edge (node) coordinates and history/station bookkeeping PCSF
needs on top.

DUHI (:mod:`pycsamt.ai.inversion`) has no separate adapter: its output
only becomes a real, final resistivity model once folded back into an
Occam2D run (see ``mapping2d.map_ai_grid_to_occam``), so a DUHI-produced
:class:`InversionResult` converts through this same function.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any, Mapping, Sequence

import numpy as np

from ...interp._base import ResistivityModel
from ..schema import Grid2DGeometry, PCSFModel, StationTable, TopographyPerStation
from ..topo_source import resolve_topo

if TYPE_CHECKING:
    from pycsamt.models.occam2d.results import InversionResult

__all__ = ["occam2d_to_pcsf"]


def _station_elevation_array(
    names: Sequence[str], station_elevations: Mapping[str, float] | None
) -> np.ndarray:
    if not station_elevations:
        return np.full(len(names), np.nan)
    return np.array(
        [float(station_elevations.get(name, np.nan)) for name in names]
    )


def _station_lonlat_arrays(
    names: Sequence[str],
    station_lonlat: Mapping[str, tuple[float, float]] | None,
) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
    if not station_lonlat:
        return None, None
    lonlat = np.array(
        [station_lonlat.get(name, (np.nan, np.nan)) for name in names],
        dtype=float,
    )
    if np.all(np.isnan(lonlat)):
        return None, None
    return lonlat[:, 0], lonlat[:, 1]


def _survey_to_dict(survey: Any | Mapping[str, Any] | None) -> dict[str, Any]:
    if survey is None:
        return {}
    to_dict = getattr(survey, "to_dict", None)
    if callable(to_dict):
        return dict(to_dict())
    return dict(survey)


def occam2d_to_pcsf(
    result: InversionResult,
    *,
    station_elevations: Mapping[str, float] | None = None,
    station_lonlat: Mapping[str, tuple[float, float]] | None = None,
    topo: Any = None,
    epsg: int | None = None,
    utm_zone: Any | None = None,
    latlon: bool = False,
    on_mismatch: str = "raise",
    survey: Any | Mapping[str, Any] | None = None,
    origin: np.ndarray | Sequence[float] | None = None,
    azimuth_deg: float | None = None,
    created_by: str = "",
    crs: str | None = None,
    description: str = "",
) -> PCSFModel:
    r"""Convert an Occam2D :class:`InversionResult` to a :class:`PCSFModel`.

    Parameters
    ----------
    result : InversionResult
        A loaded Occam2D post-inversion result (``workdir`` scanned,
        ``rho_2d`` and ``mesh`` populated).
    station_elevations : mapping of str to float, optional
        ``station_name -> elevation (m)``, e.g. from
        :func:`pycsamt.map.topo.fetch_elevations` or a survey's own
        EDI-derived topography. Occam2D itself carries no elevation, so
        this is the only way :attr:`PCSFModel.topography` gets
        populated. Stations without a known elevation are recorded as
        ``nan`` rather than a fabricated flat value.
    station_lonlat : mapping of str to (lon, lat), optional
        ``station_name -> (longitude, latitude)``, WGS84 decimal
        degrees -- e.g. from the same EDI headers a
        :class:`pycsamt.metadata.SurveyMeta` or
        :class:`pycsamt.map._core.StationRecord` would carry. Occam2D
        itself has no real-world coordinate concept (only mesh-local
        chainage, recovered into :attr:`Grid2DGeometry.x`/
        :class:`StationTable`'s own ``x``), so this is the only way a
        *single* Occam2D line's PCSF file becomes self-sufficiently
        geo-referenced -- without it, placing this line on a real
        basemap needs a separate ``known_stations`` match at
        :func:`pycsamt.map.MapView.from_pcsf` load time instead.
        Superseded per-station by *topo* when both are given (see
        below); a station *topo* has no data for still falls back to
        this mapping.
    topo : path-like, TopoTable, Sites/MapData-like, or mapping, optional
        A "smart" real-coordinate source resolved via
        :func:`pycsamt.format.topo_source.resolve_topo` — a
        ``.bln``/``.csv``/``.stn`` topo file, an already-geo-located
        ``Sites``/``MapData`` object (e.g.
        ``pycsamt.map.load_lines(edi_folder)``), or a plain
        ``{station_name: (lon, lat[, elevation])}`` mapping. A
        name-less source (a bare ``.bln``, or a ``.csv`` without a
        station column) is matched *positionally*, in
        :attr:`ResistivityModel.station_names`'s own along-profile
        order, and therefore requires exactly one point per station
        (see *on_mismatch*). When given, *topo*'s own lon/lat/
        elevation take precedence over *station_lonlat*/
        *station_elevations* for every station it resolves — with a
        :class:`UserWarning` if both were supplied, so the override is
        never silent; a station *topo* has no data for keeps whatever
        *station_lonlat*/*station_elevations* already gave it. Passing
        ``None`` (the default) leaves this adapter's behaviour exactly
        as it was before *topo* existed.
    epsg, utm_zone : optional
        Forwarded to :func:`pycsamt.format.topo_source.resolve_topo`
        for a *topo* file storing projected easting/northing rather
        than lon/lat (``.stn`` files always do; a ``.csv``/``.bln``
        does when its columns/units are projected). Ignored unless
        *topo* is given and actually needs conversion.
    latlon : bool, default False
        ``.bln`` *topo* files only: set ``True`` when the file's own
        ``x, y`` columns are already ``lon, lat`` (a ``.bln`` carries
        no CRS metadata to detect this from).
    on_mismatch : {"raise", "warn"}, default "raise"
        How a *topo* station-count mismatch is handled for a
        name-less (positional) source — see
        :func:`pycsamt.format.topo_source.attribute_topo`.
    survey : SurveyMeta or mapping, optional
        Survey-level metadata (e.g. a :class:`pycsamt.metadata.SurveyMeta`).
        Stored via its own ``to_dict()`` when available, otherwise
        copied as a plain mapping.
    origin : ndarray or sequence of float, optional
        Real-world ``(x, y)`` offset for the profile, when known.
    azimuth_deg : float, optional
        Profile bearing, when known.
    created_by, crs, description : str, optional
        Passed straight through to :class:`PCSFModel`.

    Returns
    -------
    PCSFModel
        ``geometry.kind == "grid2d"``, canonical linear-ohm.m
        resistivity in :attr:`PCSFModel.resistivity`, the original
        log10 grid preserved in :attr:`PCSFModel.resistivity_native`,
        and iteration history (RMS, roughness, Lagrange multiplier,
        step size) from :attr:`InversionResult.log` when available.

    Raises
    ------
    ValueError
        If *result* has no ``rho_2d``/``mesh`` (mirrors
        :meth:`ResistivityModel.from_occam2d`'s own check).

    Examples
    --------
    >>> from pycsamt.models.occam2d.results import InversionResult
    >>> from pycsamt.format.adapters.occam2d import occam2d_to_pcsf
    >>> from pycsamt.format import write_pcsf
    >>> result = InversionResult(workdir="data/occam2D")  # doctest: +SKIP
    >>> model = occam2d_to_pcsf(result)  # doctest: +SKIP
    >>> write_pcsf(model, "occam2d_run.pcsf")  # doctest: +SKIP
    """
    if result.rho_2d is None or result.mesh is None:
        raise ValueError(
            "InversionResult has no rho_2d or mesh — ensure the workdir "
            "contains mesh, model, and iter files."
        )

    rm = ResistivityModel.from_occam2d(result)
    mesh = result.mesh

    # Cell centres are already corrected for the mesh-local -> real
    # station-chainage shift (see ResistivityModel.from_occam2d's
    # docstring). Recover that same shift for the node (cell-edge)
    # coordinates, which the neutral ResistivityModel does not carry.
    x_centers_raw = (mesh.x_nodes[:-1] + mesh.x_nodes[1:]) / 2.0
    shift = (
        float(rm.x_centers[0] - x_centers_raw[0])
        if rm.x_centers.size and x_centers_raw.size
        else 0.0
    )
    x_nodes = mesh.x_nodes + shift
    z_nodes = mesh.z_nodes.copy()

    geometry = Grid2DGeometry(
        x=rm.x_centers,
        z=rm.z_centers,
        x_nodes=x_nodes,
        z_nodes=z_nodes,
        origin=None if origin is None else np.asarray(origin, dtype=float),
        azimuth_deg=azimuth_deg,
    )

    resistivity = 10.0**rm.rho_2d

    stations = None
    topography = None
    if rm.station_names:
        elevation = _station_elevation_array(
            rm.station_names, station_elevations
        )
        lon, lat = _station_lonlat_arrays(rm.station_names, station_lonlat)

        if topo is not None:
            if station_lonlat or station_elevations:
                warnings.warn(
                    "occam2d_to_pcsf: 'topo' takes precedence over "
                    "'station_lonlat'/'station_elevations' for every "
                    "station it resolves.",
                    UserWarning,
                    stacklevel=2,
                )
            attr = resolve_topo(
                topo,
                list(rm.station_names),
                epsg=epsg,
                utm_zone=utm_zone,
                latlon=latlon,
                on_mismatch=on_mismatch,
            )
            lon = np.array(
                [attr.lon.get(n, lon[i] if lon is not None else np.nan) for i, n in enumerate(rm.station_names)]
            )
            lat = np.array(
                [attr.lat.get(n, lat[i] if lat is not None else np.nan) for i, n in enumerate(rm.station_names)]
            )
            if np.all(np.isnan(lon)):
                lon = lat = None
            elevation = np.array(
                [attr.elevation.get(n, elevation[i]) for i, n in enumerate(rm.station_names)]
            )

        stations = StationTable(
            name=list(rm.station_names),
            x=rm.station_x,
            y=np.zeros_like(rm.station_x),
            z=elevation,
            lon=lon,
            lat=lat,
        )
        if station_elevations or (topo is not None and not np.all(np.isnan(elevation))):
            topography = TopographyPerStation(
                station_id=list(rm.station_names), elevation=elevation
            )

    history: dict[str, np.ndarray] = {}
    if result.log is not None:
        log = result.log
        history = {
            "iteration": np.asarray(log.iterations, dtype=float),
            "rms": np.asarray(log.rms, dtype=float),
            "roughness": np.asarray(log.roughness, dtype=float),
            "lagrange": np.asarray(log.lagrange, dtype=float),
            "stepsize": np.asarray(log.stepsize, dtype=float),
        }

    metadata: dict[str, Any] = {
        "workdir": str(result.workdir),
        "final_rms": rm.rms,
        "n_iterations": result.n_iterations,
    }
    if result.log is not None:
        metadata["converged"] = bool(result.log.converged)

    return PCSFModel(
        geometry=geometry,
        resistivity=resistivity,
        resistivity_native=rm.rho_2d,
        resistivity_native_encoding="log10",
        stations=stations,
        topography=topography,
        survey=_survey_to_dict(survey),
        history=history,
        source_backend="occam2d",
        created_by=created_by,
        crs=crs,
        description=description,
        metadata=metadata,
    )
