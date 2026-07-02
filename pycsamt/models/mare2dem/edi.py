# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Build MARE2DEM MT inputs directly from EDI data.

Bridges the pycsamt data model (``Sites`` / EDI collections / EDI
directories) to the MARE2DEM ``.emdata`` writer: each station's impedance
tensor becomes a :class:`~pycsamt.models.mare2dem.zmm.ZMMStation`
(TE = Zxy, TM = Zyx) which then flows through the same
profile-projection / error-floor / data-block pipeline as the ZMM path
(:func:`~pycsamt.models.mare2dem.zmm.make_mt_data_from_stations`).

Apparent resistivity follows the survey convention for field-unit
impedances, :math:`\\rho_a = 0.2\\,|Z|^2 / f`, consistent with the rest
of the package (see :mod:`pycsamt.emtools`).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .iotools.emdata import EMDataFile
from .zmm import ZMMStation, make_mt_data_from_stations

__all__ = ["stations_from_edi", "make_mt_data_from_edi"]


def _latlon(ed: Any) -> tuple[float, float]:
    """Return ``(lat, lon)`` for an EDI-like object, NaN when absent.

    Uses the canonical header accessor
    :func:`pycsamt.site.utils.get_coords` first (the coordinates live in
    the EDI HEAD section, not as attributes), then falls back to probing
    ``lat``/``latitude`` attributes on the object and its ``.edi``
    wrapper."""
    for obj in (ed, getattr(ed, "edi", None)):
        if obj is None:
            continue
        try:
            from ...site.utils import get_coords
            c = get_coords(obj)
            if c is not None:
                return float(c.lat), float(c.lon)
        except Exception:  # noqa: BLE001
            pass
        lat = getattr(obj, "lat", None)
        if lat is None:
            lat = getattr(obj, "latitude", None)
        lon = getattr(obj, "lon", None)
        if lon is None:
            lon = getattr(obj, "longitude", None)
        if lat is not None and lon is not None:
            try:
                return float(lat), float(lon)
            except (TypeError, ValueError):
                continue
    return float("nan"), float("nan")


def stations_from_edi(
    source: Any,
    *,
    default_rel_error: float = 0.05,
) -> list[ZMMStation]:
    """Convert an EDI source into :class:`ZMMStation` objects.

    Parameters
    ----------
    source : path-like, Sites, or EDI collection
        Anything accepted by :func:`pycsamt.emtools._core.ensure_sites`.
    default_rel_error : float, default 0.05
        Relative impedance error assumed when the EDIs carry no error
        block (5 %).  Error floors applied later can only raise it.

    Returns
    -------
    list of ZMMStation
        One station per EDI with valid impedance, coordinates, and a
        frequency table matching the first station's.  TE is Zxy, TM is
        Zyx; apparent resistivity uses the field-units convention
        ρ_a = 0.2 |Z|² / f and TM phase is wrapped by +180° into the
        first quadrant (same conventions as :func:`~.zmm.read_zmm`).

    Raises
    ------
    ValueError
        When no station yields usable impedance + coordinates.
    """
    from ...emtools._core import (
        ensure_sites, _iter_items, _name, _get_z_block, _unwrap,
    )

    sites = ensure_sites(source, verbose=0)
    out: list[ZMMStation] = []
    skipped: list[str] = []
    n_ref: int | None = None

    for i, ed in enumerate(_iter_items(sites)):
        name = _name(ed, i)
        got = _get_z_block(ed, with_errors=True)
        Z_obj, z, fr = got[0], got[1], got[2]
        ze = got[3] if len(got) > 3 else None
        if z is None or fr is None:
            skipped.append(name)
            continue
        try:
            raw = _unwrap(ed)
        except Exception:  # noqa: BLE001
            raw = ed
        lat, lon = _latlon(raw)
        if not (np.isfinite(lat) and np.isfinite(lon)):
            lat, lon = _latlon(ed)
        if not (np.isfinite(lat) and np.isfinite(lon)):
            skipped.append(name)
            continue
        # All stations must share the first station's frequency table:
        # the .emdata DATA block indexes one common frequency list.
        if n_ref is None:
            n_ref = len(fr)
        elif len(fr) != n_ref:
            skipped.append(name)
            continue

        # sort by ascending period (EDI frequency tables are usually
        # descending in frequency, i.e. already ascending in period)
        per = 1.0 / np.where(fr == 0, np.nan, fr)
        order = np.argsort(per)
        per = per[order]
        zxy = z[order, 0, 1]
        zyx = z[order, 1, 0]
        fr_s = fr[order]

        with np.errstate(divide="ignore", invalid="ignore"):
            apres_te = 0.2 * np.abs(zxy) ** 2 / fr_s
            apres_tm = 0.2 * np.abs(zyx) ** 2 / fr_s
            phase_te = np.degrees(np.angle(zxy))
            phase_tm = np.degrees(np.angle(zyx)) + 180.0  # TM phase wrap

            if ze is not None:
                ze_s = np.asarray(ze)[order]
                rel_te = np.abs(ze_s[:, 0, 1]) / np.maximum(
                    np.abs(zxy), 1e-30
                )
                rel_tm = np.abs(ze_s[:, 1, 0]) / np.maximum(
                    np.abs(zyx), 1e-30
                )
            else:
                rel_te = np.full(per.shape, default_rel_error)
                rel_tm = np.full(per.shape, default_rel_error)
            # guard degenerate stored errors (zero / NaN)
            rel_te = np.where(
                np.isfinite(rel_te) & (rel_te > 0),
                rel_te, default_rel_error,
            )
            rel_tm = np.where(
                np.isfinite(rel_tm) & (rel_tm > 0),
                rel_tm, default_rel_error,
            )

        st = ZMMStation(
            name=name,
            latitude=lat,
            longitude=lon,
            periods=per,
            apres_te=apres_te,
            phase_te=phase_te,
            # σ_ρ/ρ = 2 σ_Z/|Z| ;  σ_φ = σ_Z/|Z| rad
            apres_te_se=2.0 * rel_te * apres_te,
            phase_te_se=np.degrees(rel_te),
            apres_tm=apres_tm,
            phase_tm=phase_tm,
            apres_tm_se=2.0 * rel_tm * apres_tm,
            phase_tm_se=np.degrees(rel_tm),
        )
        out.append(st)

    if not out:
        raise ValueError(
            "No EDI station provided usable impedance and coordinates "
            f"(skipped: {skipped[:8]}{'…' if len(skipped) > 8 else ''})."
        )
    return out


def make_mt_data_from_edi(
    source: Any,
    out_file: str | Path,
    *,
    output_modes: str = "all",
    error_floor_te: float = 0.05,
    error_floor_tm: float = 0.05,
    error_floor_tipper: float = 0.0,
    default_rel_error: float = 0.05,
    **kwargs: Any,
) -> EMDataFile:
    """Create a MARE2DEM MT ``.emdata`` file from EDI data.

    Parameters
    ----------
    source : path-like, Sites, or EDI collection
        Anything accepted by :func:`pycsamt.emtools._core.ensure_sites`.
    out_file : path-like
        Output ``.emdata`` filename.
    output_modes, error_floor_te, error_floor_tm, error_floor_tipper
        See :func:`~.zmm.make_mt_data_from_zmm`.
    default_rel_error : float, default 0.05
        Assumed relative impedance error when the EDIs carry none.
    **kwargs :
        Remaining keyword arguments of
        :func:`~.zmm.make_mt_data_from_stations` (``omit_periods``,
        ``line_orientation``, ``utm0``, ``utm_zone``, ``topo``,
        ``rx_z_offset``, ``declination``).

    Returns
    -------
    EMDataFile
        Constructed data file (also written to *out_file*).

    Examples
    --------
    >>> from pycsamt.models.mare2dem.edi import make_mt_data_from_edi
    >>> em = make_mt_data_from_edi(
    ...     "data/AMT/WILLY_DATA/L22PLT",
    ...     "run/mare2dem.emdata",
    ...     error_floor_te=0.05, error_floor_tm=0.05,
    ... )
    >>> em.n_mt_receivers
    25
    """
    stations = stations_from_edi(
        source, default_rel_error=default_rel_error,
    )
    return make_mt_data_from_stations(
        stations, out_file,
        output_modes=output_modes,
        error_floor_te=error_floor_te,
        error_floor_tm=error_floor_tm,
        error_floor_tipper=error_floor_tipper,
        **kwargs,
    )
