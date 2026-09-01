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

#: Default cell resistivity (ohm.m) above which a ModEM 3-D cell is
#: treated as the model's above-topography "air" / padding fill rather
#: than resolved earth. A genuine rock cell is at most a few times
#: 1e4 ohm.m even for a dry, unweathered resistor; ModEM's air fill is
#: 1e8-1e13 ohm.m, so any reasonable threshold in that gap works. Used
#: only when ``modem3d_to_pcsf(air_threshold_ohm_m=...)`` is left at its
#: default -- pass ``None`` there for an exact, unmasked passthrough.
DEFAULT_AIR_THRESHOLD_OHM_M = 1e8


def _mask_air_fill(
    rho_linear: np.ndarray, threshold: float
) -> tuple[np.ndarray, int]:
    """Return ``(masked_rho, n_air)`` for a ModEM 3-D linear-ohm.m volume.

    ``masked_rho`` is a copy of *rho_linear* with every cell exceeding
    *threshold* set to ``nan`` -- ModEM leaves the model cells above the
    real ground surface filled with air-like resistivity (1e8-1e13
    ohm.m) even when the run declares no dedicated air layers, and those
    cells otherwise dominate any auto colour scale and crush the real
    earth into a featureless band (see the module docstring / the PCSF
    format plan).

    ``n_air`` counts leading *complete* air layers -- consecutive top
    z-layers in which **every** cell is above the threshold. Uneven
    topography usually means no layer is 100% air (some stations sit
    higher than others), so this is commonly ``0`` even when many
    individual cells were masked; the per-cell ``nan`` mask is the
    operative fix either way. ``rho_linear`` axis order is ModEM's own
    ``(nz, ny, nx)``.
    """
    rho = np.asarray(rho_linear, dtype=float)
    air = rho > threshold  # nan > x is False -> genuine gaps stay nan
    masked = np.where(air, np.nan, rho)
    n_air = 0
    if rho.ndim == 3 and rho.shape[0]:
        for is_air_layer in air.all(axis=(1, 2)):
            if not is_air_layer:
                break
            n_air += 1
    return masked, n_air


#: Phrases in a ModEM ``.dat`` file's own ``#`` comment lines that mark
#: its ``Z(m)`` station column as a *positive-downward depth* below a
#: datum rather than a positive-up elevation (ModEM's model z-axis
#: points down, so a populated station Z is a depth by the format's own
#: convention -- but pycsamt's own writer, and some third-party ones,
#: put a plain a.s.l. elevation there instead, so the file has to say).
_DEPTH_DOWN_MARKERS = (
    "depth below",
    "positive down",
    "positive-down",
    "z(m) is depth",
    "z is depth",
    "z down",
    "z-down",
)


def _comment_says_depth_down(comment: str | None) -> bool:
    text = (comment or "").lower()
    return any(marker in text for marker in _DEPTH_DOWN_MARKERS)


def _parse_masl_datum(comment: str | None) -> float | None:
    """Pull an ``a.s.l.`` datum elevation out of a ``.dat`` comment.

    Matches e.g. ``top = 224 m a.s.l.`` / ``datum 1035m asl`` -- the
    reference the positive-down ``Z`` column is measured from, so a real
    metres-a.s.l. station elevation can be recovered as ``datum - Z``.
    Returns ``None`` when no such hint is present (a *relative*
    elevation is used instead).
    """
    import re

    m = re.search(
        r"(-?\d+(?:\.\d+)?)\s*m\s*a\.?\s*s\.?\s*l", (comment or "").lower()
    )
    return float(m.group(1)) if m else None


def _modem_z_to_elevation(
    z: np.ndarray, comment: str | None, convention: str
) -> tuple[np.ndarray, dict[str, Any]]:
    """Return ``(elevation, info)`` -- ModEM station ``Z`` as a positive-up
    elevation.

    *convention* is ``"elevation"`` (trust the column as-is),
    ``"depth_down"`` (force ``datum - Z``), or ``"auto"`` (flip only
    when the ``.dat`` comment says the column is a positive-down depth,
    otherwise leave it -- the safe default, since a pycsamt-written
    ModEM file already stores a plain elevation there). A flat all-zero
    column is always left untouched (it is ModEM's "no topography"
    placeholder).
    """
    z = np.asarray(z, dtype=float)
    info: dict[str, Any] = {"convention": convention, "flipped": False}
    if z.size == 0 or not np.any(np.isfinite(z) & (z != 0.0)):
        return z, info
    if convention == "elevation":
        return z, info
    depth_down = convention == "depth_down" or (
        convention == "auto" and _comment_says_depth_down(comment)
    )
    if not depth_down:
        return z, info
    datum = _parse_masl_datum(comment)
    if datum is not None:
        elevation = datum - z
        info["datum_masl"] = datum
    else:
        # No a.s.l. reference in the file -- keep the *shape* (deepest
        # station at 0, higher ground positive) without inventing an
        # absolute datum.
        elevation = float(np.nanmax(z)) - z
        info["datum_masl"] = None
    info["flipped"] = True
    return elevation, info


def _warn_if_lonlat_quantized(lonlat: np.ndarray) -> None:
    """Warn when the ``.dat`` file's ``GG_Lat``/``GG_Lon`` are truncated
    to 3 decimals -- ModEM's fixed ``f9.3`` output format, ~100 m at mid
    latitudes. Common when the data file is a ModEM ``-R`` rewrite echo
    rather than the real input; on a survey with ~100 m station spacing
    it makes profile lines zig-zag and collapses nearby stations onto
    one point.
    """
    vals = np.asarray(lonlat, dtype=float).ravel()
    vals = vals[np.isfinite(vals)]
    if vals.size < 6:  # < 3 stations
        return
    if np.all(np.abs(vals - np.round(vals, 3)) < 1e-9):
        warnings.warn(
            "modem3d_to_pcsf: the '.dat' file's GG_Lat/GG_Lon are "
            "quantized to 3 decimal places (~100 m) -- likely a ModEM "
            "'-R' rewrite echo, not the real input. Profile geometry "
            "will look ragged and nearby stations may collapse onto one "
            "point. Pass a full-precision 'data=' file, or 'topo=' / "
            "'known_stations' with real coordinates.",
            UserWarning,
            stacklevel=3,
        )


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
    station_z_convention: str = "auto",
    topo: Any = None,
    epsg: int | None = None,
    utm_zone: Any | None = None,
    latlon: bool = False,
    on_mismatch: str = "raise",
) -> tuple[StationTable, dict[str, float], dict[str, Any]]:
    """Returns the built table, the *effective* elevation-override
    mapping (``station_elevations`` merged with ``topo``'s own
    elevation, ``topo`` taking priority -- the caller uses it to decide
    which stations get a :class:`TopographyPerStation` entry, since
    ModEM's own z is a real recorded value, typically a flat ``0.0``
    placeholder, not distinguishable from "known" by itself), and a
    small ``station_z`` info dict describing whether the ModEM ``Z``
    column was flipped from a positive-down depth to an elevation."""
    names = list(data.site_names)
    coords = data.site_coords
    xyz = np.array(
        [coords.get(name, (np.nan, np.nan, np.nan)) for name in names],
        dtype=float,
    )
    z = xyz[:, 2] if xyz.size else np.array([])
    # ModEM's model z-axis points *down*, so a populated station Z is a
    # depth below a datum, not an elevation -- MapView's topography
    # drape (and TopographyPerStation) need a positive-up elevation, or
    # hills render as hollows. Convert here, before any override merge.
    z, station_z_info = _modem_z_to_elevation(
        z, getattr(data, "comment", None), station_z_convention
    )
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
            _warn_if_lonlat_quantized(lonlat[np.isfinite(lonlat).all(axis=1)])

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
    return station_table, elevation_overrides, station_z_info


def modem3d_to_pcsf(
    result: InversionResult,
    *,
    model: ModEmModel3D | None = None,
    data: ModEmData | None = None,
    station_elevations: Mapping[str, float] | None = None,
    air_threshold_ohm_m: float | None = DEFAULT_AIR_THRESHOLD_OHM_M,
    station_z_convention: str = "auto",
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
    air_threshold_ohm_m : float or None, default 1e8
        Cell resistivity (ohm.m) above which a cell is treated as
        ModEM's above-topography "air" / padding fill rather than
        resolved earth, and set to ``nan`` in the canonical
        :attr:`PCSFModel.resistivity`. ModEM leaves the model cells
        above the real ground surface filled with air-like resistivity
        (1e8-1e13 ohm.m) even when the run declares **no** dedicated air
        layers (``n_air == 0``), and those cells otherwise hijack any
        auto colour scale (e.g. MapView's fence view), crushing the real
        earth into a featureless band.
        :attr:`PCSFModel.resistivity_native` keeps the source model's
        raw natural-log array **unchanged** for provenance. Leading
        *complete* air layers detected this way are also written to
        :attr:`Grid3DGeometry.n_air` (when the source model reported
        ``0``); ``metadata["air_mask"]`` records what was done. Pass
        ``None`` for an exact, unmasked passthrough of every cell.
    station_z_convention : {"auto", "elevation", "depth_down"}, default "auto"
        How to read a ModEM ``.dat`` file's own ``Z(m)`` station
        column. ModEM's model z-axis points **downward**, so a
        populated station ``Z`` is a *depth below a datum*, not an
        elevation -- copied verbatim into :attr:`StationTable.z` it
        makes MapView's topography drape render hills as hollows.
        ``"depth_down"`` forces the flip to an elevation
        (``datum - Z``; the datum is read from an ``... m a.s.l.`` hint
        in the file's ``#`` comments when present, else the deepest
        station is taken as ``0`` and the result is a *relative*
        elevation). ``"elevation"`` trusts the column as-is (what
        pycsamt's own ModEM writer produces). ``"auto"`` (default)
        flips only when the ``.dat`` comment says the column is a
        positive-down depth; a flat all-zero column is always left
        untouched. ``metadata["station_z"]`` records the outcome.
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
        (``model.rho_linear``, with above-topography air fill set to
        ``nan`` -- see *air_threshold_ohm_m*), the original natural-log
        grid preserved unchanged in
        :attr:`PCSFModel.resistivity_native`, and
        iteration history (RMS, objective, model norm, Lagrange
        multiplier, step-size scaling ``alpha``) from
        :attr:`InversionResult.log` when available.
        :attr:`PCSFModel.stations`'s ``lon``/``lat`` are populated
        from *data*'s own ``GG_Lat``/``GG_Lon`` columns
        (:attr:`ModEmData.site_lonlat`) when present, so the file is
        self-sufficiently geo-referenced without needing a separate
        ``known_stations`` match at load time.
        :attr:`StationTable.z` is a positive-up elevation -- ModEM's
        positive-down ``Z`` column is flipped when needed (see
        *station_z_convention*), so MapView's topography drape shows
        hills as hills.

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

    source_n_air = int(resolved_model.n_air)
    resistivity = resolved_model.rho_linear
    n_air = source_n_air
    air_mask_meta: dict[str, Any] | None = None
    if air_threshold_ohm_m is not None:
        masked, detected_n_air = _mask_air_fill(
            resolved_model.rho_linear, float(air_threshold_ohm_m)
        )
        newly_nan = np.isnan(masked) & ~np.isnan(resistivity)
        n_masked = int(np.count_nonzero(newly_nan))
        if n_masked:
            resistivity = masked
            # Only *raise* n_air from a source that reported none -- never
            # override a real, larger air-layer count a proper ModEM run
            # already declared.
            if source_n_air == 0 and detected_n_air > 0:
                n_air = detected_n_air
            air_mask_meta = {
                "threshold_ohm_m": float(air_threshold_ohm_m),
                "n_cells_masked": n_masked,
                "n_air_detected": detected_n_air,
            }

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
        n_air=n_air,
    )

    resolved_data = data or result.data_obs or result.data_pred
    stations = None
    elevation_overrides: dict[str, float] = {}
    station_z_info: dict[str, Any] = {}
    if resolved_data is not None:
        stations, elevation_overrides, station_z_info = (
            _stations_from_modem_data(
                resolved_data,
                station_elevations,
                station_z_convention=station_z_convention,
                topo=topo,
                epsg=epsg,
                utm_zone=utm_zone,
                latlon=latlon,
                on_mismatch=on_mismatch,
            )
        )
        if station_z_info.get("flipped"):
            datum = station_z_info.get("datum_masl")
            warnings.warn(
                "modem3d_to_pcsf: the ModEM '.dat' Z column reads as a "
                "positive-down depth; flipped it to a positive-up "
                + (
                    f"elevation (datum {datum} m a.s.l.)."
                    if datum is not None
                    else "relative elevation (deepest station = 0; no "
                    "a.s.l. datum in the file)."
                )
                + " Pass station_z_convention='elevation' to keep it "
                "verbatim.",
                UserWarning,
                stacklevel=2,
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
    if air_mask_meta is not None:
        metadata["air_mask"] = air_mask_meta
    if station_z_info:
        metadata["station_z"] = station_z_info

    return PCSFModel(
        geometry=geometry,
        resistivity=resistivity,
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
