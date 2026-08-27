# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEM 3-D -> PCSF adapter (Phase 3 of the PCSF format plan).

Converts a completed :class:`pycsamt.models.modem.results.InversionResult`
(3-D mode) into a backend-neutral :class:`~pycsamt.format.schema.PCSFModel`
with a ``grid3d`` geometry — the first genuinely new persisted 3-D volume
artifact in the project (see ``PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md``,
§1: nothing else in the codebase currently writes one to disk; the web
3-D view only ever *synthesizes* a volume at render time from stacked
2-D sections).

Unlike Occam2D, ModEM's model already carries its own real-world grid
centre and rotation when a genuine ModEM run wrote the file — see the
``origin``/``rotation_deg`` attributes added to
:class:`pycsamt.models.modem.model3d.ModEmModel3D` alongside this
adapter, which previously parsed and then silently discarded that
trailing line.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any, Mapping

import numpy as np

from ..schema import Grid3DGeometry, PCSFModel, StationTable, TopographyPerStation
from ..topo_source import resolve_topo

if TYPE_CHECKING:
    from pycsamt.models.modem.data import ModEmData
    from pycsamt.models.modem.model3d import ModEmModel3D
    from pycsamt.models.modem.results import InversionResult

__all__ = ["modem3d_to_pcsf"]


def _survey_to_dict(survey: Any | Mapping[str, Any] | None) -> dict[str, Any]:
    if survey is None:
        return {}
    to_dict = getattr(survey, "to_dict", None)
    if callable(to_dict):
        return dict(to_dict())
    return dict(survey)


def _stations_from_modem_data(
    data: ModEmData,
    station_elevations: Mapping[str, float] | None,
    *,
    topo: Any = None,
    epsg: int | None = None,
    utm_zone: Any | None = None,
    latlon: bool = False,
    on_mismatch: str = "raise",
) -> tuple[StationTable, dict[str, float]]:
    """Returns the built table plus the *effective* elevation-override
    mapping (``station_elevations`` merged with ``topo``'s own
    elevation, ``topo`` taking priority) -- the caller uses that merged
    mapping to decide which stations get a :class:`TopographyPerStation`
    entry, since ModEM's own z is a real recorded value (typically a
    flat ``0.0`` placeholder), not distinguishable from "known" by
    itself."""
    names = list(data.site_names)
    coords = data.site_coords
    xyz = np.array(
        [coords.get(name, (np.nan, np.nan, np.nan)) for name in names],
        dtype=float,
    )
    z = xyz[:, 2] if xyz.size else np.array([])
    elevation_overrides: dict[str, float] = dict(station_elevations or {})
    if elevation_overrides:
        # ModEM's own z is commonly a flat 0.0 placeholder (no real
        # topography in the .dat file); an explicit override replaces
        # it for matched stations only, real unmatched values (0.0 is
        # a real recorded value here, not "unknown") are left as-is.
        z = z.copy()
        for i, name in enumerate(names):
            if name in elevation_overrides:
                z[i] = float(elevation_overrides[name])

    # A ModEM .dat file's own GG_Lat/GG_Lon columns, when present, are
    # the only real-world spatial reference this backend carries --
    # populate StationTable.lon/lat from them directly, the same way
    # pycsamt.map.inversion.load_modem_lines already does for a *live*
    # folder import (see _resolve_lonlat there); no separate
    # known_stations match should be required just to place a grid3d
    # PCSF file's stations on a real basemap.
    lon = lat = None
    if data.site_lonlat:
        lonlat = np.array(
            [data.site_lonlat.get(name, (np.nan, np.nan)) for name in names],
            dtype=float,
        )
        if not np.all(np.isnan(lonlat)):
            lon, lat = lonlat[:, 0], lonlat[:, 1]

    if topo is not None:
        if data.site_lonlat or station_elevations:
            warnings.warn(
                "modem3d_to_pcsf: 'topo' takes precedence over the "
                "'.dat' file's own GG_Lat/GG_Lon and over "
                "'station_elevations' for every station it resolves.",
                UserWarning,
                stacklevel=3,
            )
        attr = resolve_topo(
            topo, names, epsg=epsg, utm_zone=utm_zone, latlon=latlon, on_mismatch=on_mismatch
        )
        lon = np.array(
            [attr.lon.get(n, lon[i] if lon is not None else np.nan) for i, n in enumerate(names)]
        )
        lat = np.array(
            [attr.lat.get(n, lat[i] if lat is not None else np.nan) for i, n in enumerate(names)]
        )
        if np.all(np.isnan(lon)):
            lon = lat = None
        z = np.array([attr.elevation.get(n, z[i]) for i, n in enumerate(names)])
        elevation_overrides.update(attr.elevation)

    station_table = StationTable(
        name=names,
        x=xyz[:, 0] if xyz.size else np.array([]),
        y=xyz[:, 1] if xyz.size else np.array([]),
        z=z,
        lon=lon,
        lat=lat,
    )
    return station_table, elevation_overrides


def modem3d_to_pcsf(
    result: InversionResult,
    *,
    model: ModEmModel3D | None = None,
    data: ModEmData | None = None,
    station_elevations: Mapping[str, float] | None = None,
    topo: Any = None,
    epsg: int | None = None,
    utm_zone: Any | None = None,
    latlon: bool = False,
    on_mismatch: str = "raise",
    survey: Any | Mapping[str, Any] | None = None,
    created_by: str = "",
    crs: str | None = None,
    description: str = "",
) -> PCSFModel:
    r"""Convert a ModEM 3-D :class:`InversionResult` to a :class:`PCSFModel`.

    Parameters
    ----------
    result : InversionResult
        A loaded ModEM working directory (``result.mode == "3d"``).
    model : ModEmModel3D, optional
        Model to convert. Defaults to ``result.model_final``, falling
        back to ``result.model_initial`` when no final model was
        parsed.
    data : ModEmData, optional
        Source of station coordinates. Defaults to
        ``result.data_obs``, falling back to ``result.data_pred``.
        ``None`` when neither is available (``PCSFModel.stations``
        stays ``None`` rather than fabricating positions).
    station_elevations : mapping of str to float, optional
        ``station_name -> elevation (m)``, overriding ModEM's own
        station z (commonly a flat ``0.0`` placeholder — a real ModEM
        ``.dat`` file carries no topography). Matched stations also
        populate :attr:`PCSFModel.topography`; unmatched stations keep
        their real recorded z as-is (``0.0`` is a genuine value here,
        not "unknown", unlike Occam2D's equivalent parameter).
        Superseded per-station by *topo* when both are given.
    topo : path-like, TopoTable, Sites/MapData-like, or mapping, optional
        A "smart" real-coordinate source resolved via
        :func:`pycsamt.format.topo_source.resolve_topo` -- see
        :func:`pycsamt.format.adapters.occam2d.occam2d_to_pcsf`'s
        identical parameter for the full description. When given, it
        takes precedence over both *station_elevations* and the
        ``.dat`` file's own ``GG_Lat``/``GG_Lon`` for every station it
        resolves (with a :class:`UserWarning` if either was also
        supplied); a station it has no data for keeps its existing
        value. ``None`` (the default) leaves this adapter's behaviour
        unaffected.
    epsg, utm_zone, latlon, on_mismatch
        Forwarded to :func:`pycsamt.format.topo_source.resolve_topo`;
        see *occam2d_to_pcsf*'s identical parameters.
    survey : SurveyMeta or mapping, optional
        Survey-level metadata, stored the same way as in
        :func:`pycsamt.format.adapters.occam2d.occam2d_to_pcsf`.
    created_by, crs, description : str, optional
        Passed straight through to :class:`PCSFModel`.

    Returns
    -------
    PCSFModel
        ``geometry.kind == "grid3d"``, canonical linear-ohm.m
        resistivity in :attr:`PCSFModel.resistivity`
        (``model.rho_linear``), the original natural-log grid
        preserved in :attr:`PCSFModel.resistivity_native`, and
        iteration history (RMS, objective, model norm, Lagrange
        multiplier, step-size scaling ``alpha``) from
        :attr:`InversionResult.log` when available.
        :attr:`PCSFModel.stations`'s ``lon``/``lat`` are populated
        from *data*'s own ``GG_Lat``/``GG_Lon`` columns
        (:attr:`ModEmData.site_lonlat`) when present, so the file is
        self-sufficiently geo-referenced without needing a separate
        ``known_stations`` match at load time.

    Raises
    ------
    ValueError
        If *result* is not a 3-D ModEM result, or no model
        (explicit or resolved from *result*) is available.

    Examples
    --------
    >>> from pycsamt.models.modem.results import InversionResult
    >>> from pycsamt.format.adapters.modem3d import modem3d_to_pcsf
    >>> from pycsamt.format import write_pcsf
    >>> result = InversionResult("modem_run")  # doctest: +SKIP
    >>> model = modem3d_to_pcsf(result)  # doctest: +SKIP
    >>> write_pcsf(model, "modem3d_run.pcsf")  # doctest: +SKIP
    """
    if result.mode != "3d":
        raise ValueError(
            f"modem3d_to_pcsf needs a 3-D ModEM result, got mode={result.mode!r}"
        )

    resolved_model = model or result.model_final or result.model_initial
    if resolved_model is None:
        raise ValueError(
            "InversionResult has no model_final/model_initial and no "
            "explicit model was passed — ensure the workdir contains a "
            "readable ModEM 3-D model file (.ws or Modular_NLCG_NNN.rho)."
        )

    x_nodes = resolved_model.x_nodes
    y_nodes = resolved_model.y_nodes
    z_nodes = resolved_model.z_nodes
    x_c = (x_nodes[:-1] + x_nodes[1:]) / 2.0
    y_c = (y_nodes[:-1] + y_nodes[1:]) / 2.0
    z_c = (z_nodes[:-1] + z_nodes[1:]) / 2.0

    # ModEmModel3D.origin/.rotation default to zeros/0.0 (same
    # convention as read_mackie3d) rather than None, so there is no
    # reliable way to tell "really at the origin" from "unknown" —
    # passed through as-is either way.
    geometry = Grid3DGeometry(
        x=x_c,
        y=y_c,
        z=z_c,
        x_nodes=x_nodes,
        y_nodes=y_nodes,
        z_nodes=z_nodes,
        origin=np.asarray(resolved_model.origin, dtype=float),
        rotation_deg=float(resolved_model.rotation),
        n_air=int(resolved_model.n_air),
    )

    resolved_data = data or result.data_obs or result.data_pred
    stations = None
    elevation_overrides: dict[str, float] = {}
    if resolved_data is not None:
        stations, elevation_overrides = _stations_from_modem_data(
            resolved_data,
            station_elevations,
            topo=topo,
            epsg=epsg,
            utm_zone=utm_zone,
            latlon=latlon,
            on_mismatch=on_mismatch,
        )
    topography = None
    if stations is not None and elevation_overrides:
        matched = [
            (name, elevation_overrides[name])
            for name in stations.name
            if name in elevation_overrides
        ]
        if matched:
            topography = TopographyPerStation(
                station_id=[name for name, _ in matched],
                elevation=np.array([elev for _, elev in matched], dtype=float),
            )

    history: dict[str, np.ndarray] = {}
    if result.log is not None:
        log = result.log
        history = {
            "iteration": np.asarray(log.iterations, dtype=float),
            "rms": np.asarray(log.rms, dtype=float),
            "objective": np.asarray(log.objective, dtype=float),
            "model_norm": np.asarray(log.model_norm, dtype=float),
            "lagrange": np.asarray(log.lagrange, dtype=float),
            "alpha": np.asarray(log.alpha, dtype=float),
        }

    metadata: dict[str, Any] = {
        "workdir": str(result.workdir),
        "mode": result.mode,
        "final_rms": float(result.final_rms),
        "n_iter": int(result.n_iter),
        "model_keys": list(result.models),
    }

    return PCSFModel(
        geometry=geometry,
        resistivity=resolved_model.rho_linear,
        resistivity_native=resolved_model.rho_loge,
        resistivity_native_encoding="ln",
        stations=stations,
        topography=topography,
        survey=_survey_to_dict(survey),
        history=history,
        source_backend="modem3d",
        created_by=created_by,
        crs=crs,
        description=description,
        metadata=metadata,
    )
