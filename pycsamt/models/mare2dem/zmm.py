# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ZMM impedance file reader and MARE2DEM MT data file builder.

Port of ``m2d_makeMTDataFromZmm.m`` and its helpers.

The ``.zmm`` format is the standard ASCII output of Gary Egbert's
EMTF (``occam2d``) processing code.  Each file contains the MT
impedance tensor and tipper for one station at a set of periods.

:func:`read_zmm` reads a single ``.zmm`` file.
:func:`make_mt_data_from_zmm` reads a list of ``.zmm`` files,
rotates impedances to the survey-line coordinate system, applies
error floors, projects stations onto the 2-D profile, and writes a
MARE2DEM ``.emdata`` data file.

ZMM file format (abbreviated)::

    Station  SiteName
    Latitude   42.0
    Longitude  -70.0
    Declination  0.0
    NumFreq  20
    Period   1.00000   ZXX.r  ZXX.i  ...
    ...
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from .geom.line_orientation import (
    get_line_orientation,
    project_onto_line,
)
from .geom.topo import parse_topo
from .geom.utm import lonlat_to_utm
from .iotools.emdata import (
    EMDataFile,
    MTConfig,
    UTMOrigin,
    write_emdata,
)

__all__ = [
    "ZMMStation",
    "read_zmm",
    "make_mt_data_from_zmm",
    "make_mt_data_from_stations",
]

# ---------------------------------------------------------------------------
# Data class for one ZMM station
# ---------------------------------------------------------------------------


@dataclass
class ZMMStation:
    """Impedance tensor data from one ``.zmm`` file.

    Attributes
    ----------
    name : str
        Station name.
    latitude : float
        Geographic latitude in decimal degrees.
    longitude : float
        Geographic longitude in decimal degrees.
    declination : float
        Geomagnetic declination to apply (degrees).
    periods : numpy.ndarray
        Period samples in seconds.
    apres_te : numpy.ndarray
        TE-mode apparent resistivity (Ω·m).
    phase_te : numpy.ndarray
        TE-mode phase (degrees).
    apres_te_se : numpy.ndarray
        TE apparent resistivity standard error.
    phase_te_se : numpy.ndarray
        TE phase standard error (degrees).
    apres_tm : numpy.ndarray
        TM-mode apparent resistivity.
    phase_tm : numpy.ndarray
        TM-mode phase.
    apres_tm_se : numpy.ndarray
    phase_tm_se : numpy.ndarray
    tipper_zy : numpy.ndarray or None
        Complex TE tipper values.
    tipper_zy_se : numpy.ndarray or None
        Tipper standard errors.
    """

    name: str = ""
    latitude: float = 0.0
    longitude: float = 0.0
    declination: float = 0.0
    periods: np.ndarray = field(default_factory=lambda: np.empty(0))
    apres_te: np.ndarray = field(default_factory=lambda: np.empty(0))
    phase_te: np.ndarray = field(default_factory=lambda: np.empty(0))
    apres_te_se: np.ndarray = field(default_factory=lambda: np.empty(0))
    phase_te_se: np.ndarray = field(default_factory=lambda: np.empty(0))
    apres_tm: np.ndarray = field(default_factory=lambda: np.empty(0))
    phase_tm: np.ndarray = field(default_factory=lambda: np.empty(0))
    apres_tm_se: np.ndarray = field(default_factory=lambda: np.empty(0))
    phase_tm_se: np.ndarray = field(default_factory=lambda: np.empty(0))
    tipper_zy: np.ndarray | None = None
    tipper_zy_se: np.ndarray | None = None

    # projected profile coordinates (filled by make_mt_data_from_zmm)
    x_profile: float = 0.0
    y_profile: float = 0.0
    z_profile: float = 0.0


# ---------------------------------------------------------------------------
# ZMM reader
# ---------------------------------------------------------------------------


def read_zmm(path: str | Path) -> ZMMStation:
    """Read one EMTF ``.zmm`` impedance file.

    Parameters
    ----------
    path : path-like
        File to read.

    Returns
    -------
    ZMMStation
        Parsed station data.

    Raises
    ------
    FileNotFoundError
        When *path* does not exist.

    Notes
    -----
    The reader parses the ZMM 2-column-per-component format:

        Period  ZXX.r  ZXX.i  ZXY.r  ZXY.i  ZYX.r  ZYX.i  ZYY.r  ZYY.i \\
                TZX.r  TZX.i  TZY.r  TZY.i  Coh_XY  Coh_YX

    Apparent resistivities and phases are derived from the complex
    impedance components ZXY (TE) and ZYX (TM) using:

        ρ_a = |Z|² / (ω · μ₀)
        φ = atan2(Z.imag, Z.real)

    where ``ω = 2π / T`` and ``μ₀ = 4π × 10⁻⁷`` H/m.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"ZMM file not found: {path}")

    st = ZMMStation()
    st.name = path.stem

    lines = path.read_text(errors="replace").splitlines()
    i = 0
    periods_list: list[float] = []
    zxy_list: list[complex] = []
    zyx_list: list[complex] = []
    tzy_list: list[complex] = []
    coh_xy_list: list[float] = []
    coh_yx_list: list[float] = []

    while i < len(lines):
        raw = lines[i].strip()
        i += 1
        low = raw.lower()

        if low.startswith("station") or low.startswith("site"):
            parts = raw.split()
            if len(parts) >= 2:
                st.name = parts[1]

        elif low.startswith("latitude"):
            try:
                st.latitude = float(raw.split()[1])
            except (IndexError, ValueError):
                pass

        elif low.startswith("longitude"):
            try:
                st.longitude = float(raw.split()[1])
            except (IndexError, ValueError):
                pass

        elif low.startswith("declination"):
            try:
                st.declination = float(raw.split()[1])
            except (IndexError, ValueError):
                pass

        elif low.startswith("numfreq") or low.startswith("nfreq"):
            try:
                int(raw.split()[1])
            except (IndexError, ValueError):
                pass

        elif low.startswith("period"):
            # Data line: Period  Zxx.r Zxx.i  Zxy.r Zxy.i  Zyx.r Zyx.i  Zyy.r Zyy.i  [Tzx.r Tzx.i Tzy.r Tzy.i ...]
            nums = [float(v) for v in raw.split() if _is_number(v)]
            if len(nums) >= 9:
                T = nums[0]
                periods_list.append(T)
                zxy_r, zxy_i = nums[3], nums[4]
                zyx_r, zyx_i = nums[5], nums[6]
                zxy_list.append(complex(zxy_r, zxy_i))
                zyx_list.append(complex(zyx_r, zyx_i))

                if len(nums) >= 13:
                    tzy_r, tzy_i = nums[11], nums[12]
                    tzy_list.append(complex(tzy_r, tzy_i))
                else:
                    tzy_list.append(complex(0.0, 0.0))

                if len(nums) >= 15:
                    coh_xy_list.append(nums[13])
                    coh_yx_list.append(nums[14])
                else:
                    coh_xy_list.append(1.0)
                    coh_yx_list.append(1.0)

    if not periods_list:
        return st

    MU0 = 4e-7 * math.pi
    periods = np.array(periods_list)
    omega = 2.0 * math.pi / periods

    zxy = np.array(zxy_list)
    zyx = np.array(zyx_list)
    tzy = np.array(tzy_list) if tzy_list else None

    # apparent resistivity (Ω·m) and phase (degrees)
    st.periods = periods

    with np.errstate(divide="ignore", invalid="ignore"):
        st.apres_te = np.abs(zxy) ** 2 / (omega * MU0)
        st.phase_te = np.degrees(np.angle(zxy))

        st.apres_tm = np.abs(zyx) ** 2 / (omega * MU0)
        st.phase_tm = np.degrees(np.angle(zyx)) + 180.0  # TM phase wrap

        # Standard errors: use coherence as proxy (σ_Z ≈ (1-coh)*|Z|)
        coh_xy = np.clip(np.array(coh_xy_list), 0, 1)
        coh_yx = np.clip(np.array(coh_yx_list), 0, 1)

        se_z_xy = (1.0 - coh_xy) * np.abs(zxy)
        se_z_yx = (1.0 - coh_yx) * np.abs(zyx)

        # propagate to apres: σ_ρ = 2|Z|σ_Z/(ω·μ₀)
        st.apres_te_se = 2 * np.abs(zxy) * se_z_xy / (omega * MU0)
        st.apres_tm_se = 2 * np.abs(zyx) * se_z_yx / (omega * MU0)

        # propagate to phase: σ_φ = σ_Z / |Z| (radians) → degrees
        st.phase_te_se = np.degrees(se_z_xy / np.maximum(np.abs(zxy), 1e-30))
        st.phase_tm_se = np.degrees(se_z_yx / np.maximum(np.abs(zyx), 1e-30))

    if tzy is not None and np.any(tzy != 0):
        st.tipper_zy = tzy
        st.tipper_zy_se = (
            np.abs(tzy) * 0.1
        )  # 10% as default; no direct estimate

    return st


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


# ---------------------------------------------------------------------------
# Error-floor application
# ---------------------------------------------------------------------------


def _apply_error_floor(
    stations: list[ZMMStation],
    err_floor_te: float,
    err_floor_tm: float,
    err_floor_tipper: float,
) -> list[ZMMStation]:
    for s in stations:
        if err_floor_te > 0:
            floor = s.apres_te * err_floor_te
            s.apres_te_se = np.maximum(s.apres_te_se, floor)
            s.phase_te_se = np.maximum(
                s.phase_te_se, err_floor_te * 180.0 / math.pi / 2
            )
        if err_floor_tm > 0:
            floor = s.apres_tm * err_floor_tm
            s.apres_tm_se = np.maximum(s.apres_tm_se, floor)
            s.phase_tm_se = np.maximum(
                s.phase_tm_se, err_floor_tm * 180.0 / math.pi / 2
            )
        if err_floor_tipper > 0 and s.tipper_zy_se is not None:
            s.tipper_zy_se = np.maximum(s.tipper_zy_se, err_floor_tipper)
    return stations


# ---------------------------------------------------------------------------
# Main public function
# ---------------------------------------------------------------------------


def make_mt_data_from_zmm(
    zmm_files: list[str | Path],
    out_file: str | Path,
    *,
    output_modes: str = "all",
    error_floor_te: float = 0.0,
    error_floor_tm: float = 0.0,
    error_floor_tipper: float = 0.0,
    omit_periods: np.ndarray | None = None,
    line_orientation: float | None = None,
    declination: float = 0.0,
    utm0: tuple[float, float] | None = None,
    utm_zone: str = "",
    topo: Any = 0.0,
    rx_z_offset: float = -0.1,
) -> EMDataFile:
    """Create a MARE2DEM MT data file from a list of ``.zmm`` files.

    Port of ``m2d_makeMTDataFromZmm.m``.

    Parameters
    ----------
    zmm_files : list of path-like
        ``.zmm`` impedance files in profile order.
    out_file : path-like
        Output ``.emdata`` filename.
    output_modes : str, default "all"
        Comma-separated data types to include:
        ``"TE"``, ``"TM"``, ``"tipper"``, ``"TE+tipper"``,
        ``"all impedance"``, or ``"all"`` (TE + TM + tipper).
    error_floor_te : float
        Relative error floor for TE (e.g. 0.05 = 5 %).
    error_floor_tm : float
        Relative error floor for TM.
    error_floor_tipper : float
        Absolute tipper error floor.
    omit_periods : array-like, shape (n, 2) or None
        Period bands to omit: each row is ``[T_min, T_max]`` in seconds.
    line_orientation : float or None
        Profile orientation in degrees clockwise from North.  Auto-fit
        from station positions when ``None``.
    declination : float, default 0.0
        Geomagnetic declination correction in degrees.
    utm0 : (float, float) or None
        Profile UTM origin ``(northing, easting)`` in metres.
    utm_zone : str, default ""
        UTM zone string, e.g. ``"19N"``.  Auto-detected when empty.
    topo : float or array-like
        Topography for receiver placement (see :func:`parse_topo`).
    rx_z_offset : float, default -0.1
        Vertical offset of receivers relative to topography (m).
        Negative → above (marine), positive → below (land).

    Returns
    -------
    EMDataFile
        Constructed data file (also written to *out_file*).

    Examples
    --------
    >>> from pycsamt.models.mare2dem.zmm import make_mt_data_from_zmm
    >>> em = make_mt_data_from_zmm(
    ...     ["S001.zmm", "S002.zmm"],
    ...     "line1_mt.emdata",
    ...     output_modes="TE+tipper",
    ...     error_floor_te=0.05,
    ... )
    >>> em.n_mt_receivers
    2
    """
    # ---- read ----
    stations = [read_zmm(f) for f in zmm_files]
    return make_mt_data_from_stations(
        stations,
        out_file,
        output_modes=output_modes,
        error_floor_te=error_floor_te,
        error_floor_tm=error_floor_tm,
        error_floor_tipper=error_floor_tipper,
        omit_periods=omit_periods,
        line_orientation=line_orientation,
        declination=declination,
        utm0=utm0,
        utm_zone=utm_zone,
        topo=topo,
        rx_z_offset=rx_z_offset,
    )


def make_mt_data_from_stations(
    stations: list[ZMMStation],
    out_file: str | Path,
    *,
    output_modes: str = "all",
    error_floor_te: float = 0.0,
    error_floor_tm: float = 0.0,
    error_floor_tipper: float = 0.0,
    omit_periods: np.ndarray | None = None,
    line_orientation: float | None = None,
    declination: float = 0.0,
    utm0: tuple[float, float] | None = None,
    utm_zone: str = "",
    topo: Any = 0.0,
    rx_z_offset: float = -0.1,
) -> EMDataFile:
    """Create a MARE2DEM MT data file from prepared :class:`ZMMStation` objects.

    Backend shared by :func:`make_mt_data_from_zmm` (stations parsed from
    ``.zmm`` files) and
    :func:`pycsamt.models.mare2dem.edi.make_mt_data_from_edi` (stations
    converted from EDI impedances).  See :func:`make_mt_data_from_zmm` for
    the parameter documentation; the only difference is that *stations*
    are supplied directly.
    """
    # ---- UTM zone ----
    if not utm_zone:
        lons = [s.longitude for s in stations]
        lats = [s.latitude for s in stations]
        _, _, zone_int, sh = lonlat_to_utm(np.array(lons), np.array(lats))
        hemi = "S" if sh else "N"
        utm_zone = f"{int(zone_int)}{hemi}"
    else:
        zone_int = int(re.match(r"(\d+)", utm_zone).group(1))
        sh = "S" in utm_zone.upper()

    # ---- convert to UTM ----
    for s in stations:
        e, n, _, _ = lonlat_to_utm(
            np.array([s.longitude]),
            np.array([s.latitude]),
            zone=zone_int,
            south_hemi=sh,
        )
        s._easting = float(e[0])
        s._northing = float(n[0])

    northings = np.array([s._northing for s in stations])
    eastings = np.array([s._easting for s in stations])

    # ---- line orientation ----
    if line_orientation is None:
        line_orientation = get_line_orientation(northings, eastings)

    # ---- rotation angle and UTM0 ----
    # sort stations by position along profile
    if line_orientation >= 315 or line_orientation <= 45:
        sort_idx = np.argsort(northings)
    elif 135 < line_orientation <= 225:
        sort_idx = np.argsort(northings)[::-1]
    elif 45 < line_orientation <= 135:
        sort_idx = np.argsort(eastings)
    else:
        sort_idx = np.argsort(eastings)[::-1]

    stations = [stations[i] for i in sort_idx]
    northings = northings[sort_idx]
    eastings = eastings[sort_idx]

    if utm0 is None:
        utm0 = (float(northings[0]), float(eastings[0]))

    # project onto profile
    x_arr, y_arr = project_onto_line(
        northings, eastings, utm0[0], utm0[1], line_orientation
    )
    for i, s in enumerate(stations):
        s.x_profile = float(x_arr[i])
        s.y_profile = float(y_arr[i])
        s.z_profile = 0.0

    # ---- omit periods ----
    if omit_periods is not None:
        op = np.asarray(omit_periods, dtype=float).reshape(-1, 2)
        ref_periods = stations[0].periods if stations else np.array([])
        keep = np.ones(len(ref_periods), dtype=bool)
        for row in op:
            keep &= ~((ref_periods >= row[0]) & (ref_periods <= row[1]))
        for s in stations:
            for attr in (
                "periods",
                "apres_te",
                "phase_te",
                "apres_te_se",
                "phase_te_se",
                "apres_tm",
                "phase_tm",
                "apres_tm_se",
                "phase_tm_se",
            ):
                arr = getattr(s, attr)
                if len(arr) == len(keep):
                    setattr(s, attr, arr[keep])
            if s.tipper_zy is not None and len(s.tipper_zy) == len(keep):
                s.tipper_zy = s.tipper_zy[keep]
                s.tipper_zy_se = (
                    s.tipper_zy_se[keep]
                    if s.tipper_zy_se is not None
                    else None
                )

    # ---- apply error floors ----
    stations = _apply_error_floor(
        stations, error_floor_te, error_floor_tm, error_floor_tipper
    )

    # ---- place receivers relative to topography ----
    y_rx = np.array([s.y_profile for s in stations])
    z_topo, slope, _ = parse_topo(topo, y_rx)
    z_rx = z_topo + rx_z_offset
    beta = slope.copy()

    # ---- MT config ----
    if not stations:
        em = EMDataFile()
        write_emdata(em, out_file)
        return em

    periods = stations[0].periods
    freqs = 1.0 / periods
    n_rx = len(stations)
    rx_arr = np.zeros((n_rx, 8))
    for i, s in enumerate(stations):
        rx_arr[i, :] = [
            s.x_profile,
            s.y_profile,
            z_rx[i],
            0.0,
            0.0,
            beta[i],
            0.0,
            0.0,
        ]
    rx_names = [s.name for s in stations]

    mt_cfg = MTConfig(
        frequencies=freqs,
        receivers=rx_arr,
        receiver_name=rx_names,
    )

    # ---- DATA block ----
    modes = output_modes.lower().strip()
    rows: list[list[float]] = []

    for ifreq in range(1, len(freqs) + 1):
        periods[ifreq - 1]
        for irx, s in enumerate(stations, start=1):
            if ifreq > len(s.periods):
                continue
            fi = ifreq - 1

            if modes in ("all", "all impedance", "te", "te+tipper"):
                if s.apres_te_se[fi] > 0 and not np.isnan(s.apres_te[fi]):
                    rows.append(
                        [
                            123,
                            ifreq,
                            irx,
                            irx,
                            math.log10(max(s.apres_te[fi], 1e-30)),
                            s.apres_te_se[fi] / s.apres_te[fi] * 0.4343,
                        ]
                    )
                    rows.append(
                        [
                            104,
                            ifreq,
                            irx,
                            irx,
                            s.phase_te[fi],
                            s.phase_te_se[fi],
                        ]
                    )

            if modes in ("all", "all impedance", "tm"):
                if s.apres_tm_se[fi] > 0 and not np.isnan(s.apres_tm[fi]):
                    rows.append(
                        [
                            125,
                            ifreq,
                            irx,
                            irx,
                            math.log10(max(s.apres_tm[fi], 1e-30)),
                            s.apres_tm_se[fi] / s.apres_tm[fi] * 0.4343,
                        ]
                    )
                    rows.append(
                        [
                            106,
                            ifreq,
                            irx,
                            irx,
                            s.phase_tm[fi],
                            s.phase_tm_se[fi],
                        ]
                    )

            if modes in ("all", "tipper", "te+tipper"):
                if (
                    s.tipper_zy is not None
                    and fi < len(s.tipper_zy)
                    and s.tipper_zy_se is not None
                ):
                    rows.append(
                        [
                            133,
                            ifreq,
                            irx,
                            irx,
                            s.tipper_zy[fi].real,
                            s.tipper_zy_se[fi],
                        ]
                    )
                    rows.append(
                        [
                            134,
                            ifreq,
                            irx,
                            irx,
                            s.tipper_zy[fi].imag,
                            s.tipper_zy_se[fi],
                        ]
                    )

    data = (
        np.array(rows, dtype=float) if rows else np.empty((0, 6), dtype=float)
    )

    # ---- UTM block ----
    theta = line_orientation - 90.0
    utm = UTMOrigin(
        grid=int(zone_int),
        hemi="S" if sh else "N",
        north0=utm0[0],
        east0=utm0[1],
        theta=theta,
    )

    em = EMDataFile(
        format="EMData_2.3",
        comment=f"Created with pycsamt make_mt_data_from_stations on {datetime.now().strftime('%Y-%m-%d')}",
        utm=utm,
        mt=mt_cfg,
        data=data,
    )
    write_emdata(em, out_file)
    return em
