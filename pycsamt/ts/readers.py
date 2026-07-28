# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Readers for magnetotelluric field time-series formats.

The entry point is :func:`read_ts`, which sniffs the file
format and dispatches to a concrete reader.  All readers
return a :class:`pycsamt.ts.TSData`.

Supported formats
-----------------
``lims``
    SAMTEX / GSC LiMS ASCII produced by ``mp2ts`` /
    ``tssplice`` (e.g. ``data/MT/TS/kap103as.ts``): ``#``
    comment block, ``>KEY :value`` INFO block terminated by
    ``>INFO_END``, then free-format rows of ``NCHAN`` floats.
    Missing samples equal the ``MIS_DATA`` sentinel.
``emslab``
    EMSLAB Lincoln-Line ASCII (e.g.
    ``data/MT/TS/emsl01.asc``): hour records made of a
    16-character header ``<STATION><yymmddhhmm>`` followed by
    60 lines in ``15I5`` fixed width — 900 integers = 180
    samples of the multiplexed channels Hx, Hy, Hz, Ex, Ey in
    units of 0.1 nT / 0.1 mV/km; ``-9999`` marks missing.
``edi``
    A SEG EDI file carrying ``>=TSERIESSECT`` / ``>TSERIES``
    blocks (delegates to :mod:`pycsamt.seg.time_series`).
``ascii``
    Generic whitespace table, one column per channel
    (requires ``dt`` and optionally ``chan``).

References
----------
.. [1] SEG EDI MT/EMAP standard (1987).  MTNet.
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
.. [2] Jones, A.G. et al. (2009). Area selection for diamonds
       using magnetotellurics. *Lithos*, 112S, 83-92 (SAMTEX).
.. [3] EMSLAB Group (1988). The EMSLAB electromagnetic
       sounding experiment. *EOS Trans. AGU*, 69.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

from ..exceptions import FileHandlingError
from ..log.logger import get_logger
from ._base import CHANNEL_ORDER, TSData

logger = get_logger(__name__)

__all__ = [
    "read_ts",
    "read_lims",
    "read_emslab",
    "read_edi_tseries",
    "read_ascii",
    "sniff_format",
]

# EMSLAB hour-record header: 6-char station + yymmddhhmm
_EMSLAB_HDR = re.compile(r"^([A-Z][A-Z0-9]{2,5})(\d{10})\s*$")
#: EMSLAB multiplex order (README.emslab)
_EMSLAB_CHAN: tuple[str, ...] = ("HX", "HY", "HZ", "EX", "EY")
_EMSLAB_UNITS = {
    "HX": "nT",
    "HY": "nT",
    "HZ": "nT",
    "EX": "mV/km",
    "EY": "mV/km",
}


# --------------------------------------------------------------- sniffing
def sniff_format(path: str | Path) -> str:
    """
    Guess the time-series file format.

    Returns one of ``"lims"``, ``"emslab"``, ``"edi"`` or
    ``"ascii"`` (fallback).
    """
    p = Path(path)
    head: list[str] = []
    with p.open("r", encoding="utf-8-sig", errors="replace") as f:
        for _ in range(400):
            ln = f.readline()
            if not ln:
                break
            head.append(ln.rstrip("\n"))

    stripped = [ln.strip() for ln in head if ln.strip()]
    if not stripped:
        raise FileHandlingError(f"{p} is empty.")

    first = stripped[0].upper()
    if first.startswith(">HEAD"):
        return "edi"
    if _EMSLAB_HDR.match(stripped[0]):
        return "emslab"
    for ln in stripped:
        u = ln.upper()
        if u.startswith(">INFO_START") or u.startswith(">NCHAN"):
            return "lims"
    return "ascii"


# ------------------------------------------------------------ LiMS reader
def _lims_header(
    lines: list[str],
) -> tuple[dict[str, str], list[str], int]:
    """
    Split a LiMS file head into ``(info, comments, data_start)``.

    ``info`` maps upper-case INFO keys to raw string values,
    ``comments`` keeps the ``#`` block, and ``data_start`` is
    the line index of the first sample row.
    """
    info: dict[str, str] = {}
    comments: list[str] = []
    data_start = 0
    for i, raw in enumerate(lines):
        s = raw.strip()
        if not s:
            continue
        if s.startswith("#"):
            comments.append(s)
            continue
        if s.startswith(">"):
            body = s.lstrip(">").strip()
            if ":" in body:
                key, val = body.split(":", 1)
                info[key.strip().upper()] = val.strip()
            if body.upper().startswith("INFO_END"):
                data_start = i + 1
                break
            continue
        # a bare numeric row before INFO_END should not
        # happen; treat it as the start of data anyway.
        data_start = i
        break
    else:  # no break — header only
        data_start = len(lines)
    return info, comments, data_start


def read_lims(
    path: str | Path,
    *,
    verbose: int = 0,
) -> TSData:
    """
    Read a SAMTEX / LiMS ``mp2ts`` / ``tssplice`` ASCII file.

    Parameters
    ----------
    path : str or Path
        File to read (e.g. ``kap103as.ts``).
    verbose : int, default 0
        Verbosity level.

    Returns
    -------
    TSData
        Channels as declared by ``CHAN_n`` (typically HX, HY,
        HZ, EX, EY), ``MIS_DATA`` sentinel converted to
        ``NaN``, and site metadata populated.
    """
    p = Path(path)
    if not p.is_file():
        raise FileHandlingError(f"{p} is not a file.")

    # header is small: read the first chunk line by line
    head_lines: list[str] = []
    with p.open("r", encoding="utf-8-sig", errors="replace") as f:
        for _ in range(2000):
            ln = f.readline()
            if not ln:
                break
            head_lines.append(ln.rstrip("\n"))

    info, comments, data_start = _lims_header(head_lines)
    if not info:
        raise FileHandlingError(f"{p} has no LiMS '>KEY :value' INFO block.")

    nchan = int(float(info.get("NCHAN", 0) or 0))
    if nchan <= 0:
        raise FileHandlingError(f"{p}: NCHAN missing or invalid in INFO block.")

    chan: list[str] = []
    azim: dict[str, float] = {}
    units: dict[str, str] = {}
    gain: dict[str, float] = {}
    baseline: dict[str, float] = {}
    sensor: dict[str, str] = {}
    for k in range(1, nchan + 1):
        cid = (info.get(f"CHAN_{k}", f"CH{k}") or f"CH{k}").upper()
        chan.append(cid)
        if f"AZIM_{k}" in info:
            try:
                azim[cid] = float(info[f"AZIM_{k}"])
            except ValueError:
                pass
        if f"UNITS_{k}" in info:
            units[cid] = info[f"UNITS_{k}"]
        if f"GAIN_{k}" in info:
            try:
                gain[cid] = float(info[f"GAIN_{k}"])
            except ValueError:
                pass
        if f"BASELINE_{k}" in info:
            try:
                baseline[cid] = float(info[f"BASELINE_{k}"])
            except ValueError:
                pass
        if f"SENSOR_{k}" in info:
            sensor[cid] = info[f"SENSOR_{k}"]

    def _f(key: str) -> float | None:
        v = info.get(key)
        if v in (None, ""):
            return None
        try:
            return float(v)
        except ValueError:
            return None

    dt = _f("DELTA_T")
    missing = _f("MIS_DATA")

    # electric dipole lengths live in the '#' comments
    dipole: dict[str, float] = {}
    for c in comments:
        m = re.search(
            r"#\s*(E[xy])\s+line\s+length\s*\(m\)\s*:\s*" r"([0-9.+-Ee]+)",
            c,
        )
        if m:
            try:
                dipole[m.group(1).upper()] = float(m.group(2))
            except ValueError:
                pass

    # ---- bulk-load the numeric payload (fast C parser).
    # sep='\s+' maps to the pandas C-engine whitespace mode.
    try:
        df = pd.read_csv(
            p,
            sep=r"\s+",
            skiprows=data_start,
            header=None,
            comment="#",
            dtype=np.float64,
        )
        values = df.to_numpy()
    except Exception as exc:  # pragma: no cover - fallback
        logger.warning(
            "pandas parse failed (%s); using numpy loader",
            exc,
        )
        values = np.loadtxt(str(p), skiprows=data_start, comments="#")

    if values.ndim == 1:
        values = values[None, :]
    if values.shape[1] != nchan:
        raise FileHandlingError(
            f"{p}: expected {nchan} columns, found {values.shape[1]}."
        )

    if missing is not None:
        bad = np.abs(values) >= 0.99 * abs(missing)
        if bad.any():
            values = values.copy()
            values[bad] = np.nan

    data = {cid: values[:, j] for j, cid in enumerate(chan)}

    ts = TSData(
        data=data,
        dt=dt,
        name=info.get("WINDOW") or info.get("STATION"),
        verbose=verbose,
        station=info.get("STATION"),
        lat=_f("LATITUDE"),
        lon=_f("LONGITUDE"),
        elev=_f("ELEVATION"),
        declination=_f("DECLIN"),
        coordsys=info.get("COORD_SYS"),
        start=info.get("STARTTIME"),
        stop=info.get("ENDTIME"),
        azim=azim,
        units=units,
        gain=gain,
        baseline=baseline,
        dipole=dipole,
        sensor=sensor,
        missing=missing,
        source=str(p),
        format="lims",
        instrument=info.get("INSTRUMENT"),
        comments=comments,
        info=info,
    )
    if verbose:
        logger.info(
            "LiMS %s: %d chan x %d samples, dt=%s s",
            ts.station,
            ts.n_chan,
            ts.n_samples,
            ts.dt,
        )
    return ts


# --------------------------------------------------------- EMSLAB reader
def _emslab_time(stamp: str) -> datetime:
    """``yymmddhhmm`` → :class:`datetime` (1950–2049 pivot)."""
    yy = int(stamp[:2])
    year = 2000 + yy if yy < 50 else 1900 + yy
    return datetime(
        year,
        int(stamp[2:4]),
        int(stamp[4:6]),
        int(stamp[6:8]),
        int(stamp[8:10]),
    )


def read_emslab(
    path: str | Path,
    *,
    dt: float = 20.0,
    scale: float = 0.1,
    chan: Sequence[str] = _EMSLAB_CHAN,
    declination: float | None = None,
    verbose: int = 0,
) -> TSData:
    """
    Read an EMSLAB Lincoln-Line multiplexed ASCII file.

    Parameters
    ----------
    path : str or Path
        File to read (e.g. ``emsl01.asc``).
    dt : float, default 20.0
        Sampling interval in seconds.
    scale : float, default 0.1
        Multiplier applied to the stored integers (data are
        recorded in 0.1 nT and 0.1 mV/km).
    chan : sequence of str
        Multiplex order (default Hx, Hy, Hz, Ex, Ey).
    declination : float, optional
        Declination to record in metadata.  The EMSLAB README
        suggests ``-19.5`` degrees for all Lincoln-Line sites
        (data are in geomagnetic coordinates).
    verbose : int, default 0
        Verbosity level.

    Returns
    -------
    TSData
        Continuous record; hours missing from the file are
        padded with ``NaN`` so sample timing stays uniform.
    """
    p = Path(path)
    if not p.is_file():
        raise FileHandlingError(f"{p} is not a file.")

    nch = len(chan)
    station: str | None = None
    hours: list[tuple[datetime, list[str]]] = []
    current: list[str] | None = None

    with p.open("r", encoding="utf-8-sig", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            m = _EMSLAB_HDR.match(line.strip())
            if m:
                try:
                    tstamp = _emslab_time(m.group(2))
                except ValueError:
                    tstamp = None  # not a real timestamp
                if tstamp is not None:
                    station = station or m.group(1)
                    current = []
                    hours.append((tstamp, current))
                    continue
            if current is None:
                # ignore any preamble before the first header
                continue
            if line.strip():
                current.append(line)

    if not hours:
        raise FileHandlingError(f"{p}: no EMSLAB hour records found.")

    # fixed-width decode of one hour block -> (nsamp, nch)
    def _decode(lines: list[str]) -> np.ndarray:
        vals: list[float] = []
        for ln in lines:
            # 15 fields of width 5 (last line may be short)
            n = min(len(ln), 75)
            for i in range(0, n, 5):
                tok = ln[i : i + 5].strip()
                if tok:
                    vals.append(float(tok))
        arr = np.asarray(vals, dtype=float)
        usable = (arr.size // nch) * nch
        arr = arr[:usable].reshape(-1, nch)
        return arr

    t0 = hours[0][0]
    blocks: list[np.ndarray] = []
    prev_end: datetime | None = None
    samples_per_hour = int(round(3600.0 / dt))

    for tstamp, lines in hours:
        if prev_end is not None and tstamp > prev_end:
            # pad fully missing hours to keep timing uniform
            gap_h = int(round((tstamp - prev_end).total_seconds() / 3600.0))
            if gap_h > 0:
                blocks.append(
                    np.full(
                        (gap_h * samples_per_hour, nch),
                        np.nan,
                    )
                )
        block = _decode(lines)
        if block.shape[0] < samples_per_hour:
            pad = np.full(
                (samples_per_hour - block.shape[0], nch),
                np.nan,
            )
            block = np.vstack([block, pad])
        blocks.append(block[:samples_per_hour])
        prev_end = tstamp + timedelta(hours=1)

    values = np.vstack(blocks)
    values = np.where(values == -9999.0, np.nan, values)
    values = values * float(scale)

    data = {str(cid).upper(): values[:, j] for j, cid in enumerate(chan)}
    # first sample is one interval AFTER the hour mark
    start = t0 + timedelta(seconds=dt)
    stop = start + timedelta(seconds=dt * (values.shape[0] - 1))

    ts = TSData(
        data=data,
        dt=dt,
        name=station,
        verbose=verbose,
        station=station,
        declination=declination,
        coordsys="GEOMAGNETIC NORTH",
        start=start.strftime("%Y-%m-%d %H:%M:%S"),
        stop=stop.strftime("%Y-%m-%d %H:%M:%S"),
        units={str(c).upper(): _EMSLAB_UNITS.get(str(c).upper(), "") for c in chan},
        missing=-9999.0,
        source=str(p),
        format="emslab",
    )
    if verbose:
        logger.info(
            "EMSLAB %s: %d chan x %d samples, dt=%s s",
            station,
            ts.n_chan,
            ts.n_samples,
            ts.dt,
        )
    return ts


# ------------------------------------------------------------ EDI reader
def read_edi_tseries(
    path: str | Path,
    *,
    verbose: int = 0,
) -> TSData:
    """
    Read ``>TSERIES`` blocks embedded in a SEG EDI file.

    Thin bridge over :class:`pycsamt.seg.time_series.TSect`
    and :class:`~pycsamt.seg.time_series.TSIO`.
    """
    from ..seg.time_series import TSIO, TimeSeries, TSect

    sect = TSect.from_file(str(path))
    io = TSIO.from_file(str(path), start_line=sect.start_data_lines_num)
    seg_ts = TimeSeries.from_io(sect, io)

    dt = None
    data: dict[str, np.ndarray] = {}
    for cid in seg_ts.channels():
        data[cid] = np.asarray(seg_ts.get(cid), float)
        dt = seg_ts.dt_map.get(cid, dt)

    return TSData(
        data=data,
        dt=dt,
        name=getattr(sect, "sectid", None),
        verbose=verbose,
        station=getattr(sect, "sectid", None),
        source=str(path),
        format="edi",
    )


# -------------------------------------------------------- generic ASCII
def read_ascii(
    path: str | Path,
    *,
    dt: float,
    chan: Sequence[str] | None = None,
    missing: float | None = None,
    skiprows: int = 0,
    verbose: int = 0,
) -> TSData:
    """
    Read a plain whitespace table (one column per channel).

    Parameters
    ----------
    path : str or Path
        Text file with numeric columns.
    dt : float
        Sampling interval in seconds (required).
    chan : sequence of str, optional
        Column names.  Defaults to the first ``n`` entries of
        :data:`pycsamt.ts.CHANNEL_ORDER`.
    missing : float, optional
        Sentinel converted to ``NaN``.
    skiprows : int, default 0
        Header rows to skip.
    """
    p = Path(path)
    df = pd.read_csv(
        p,
        sep=r"\s+",
        skiprows=skiprows,
        header=None,
        comment="#",
        dtype=np.float64,
    )
    values = df.to_numpy()
    ncol = values.shape[1]
    order = (
        [str(c).upper() for c in chan]
        if chan is not None
        else list(CHANNEL_ORDER[:ncol])
    )
    if len(order) != ncol:
        raise FileHandlingError(
            f"{p}: {ncol} columns but {len(order)} channel names given."
        )
    if missing is not None:
        values = np.where(
            np.abs(values) >= 0.99 * abs(missing),
            np.nan,
            values,
        )
    data = {cid: values[:, j] for j, cid in enumerate(order)}
    return TSData(
        data=data,
        dt=dt,
        verbose=verbose,
        missing=missing,
        source=str(p),
        format="ascii",
    )


# ---------------------------------------------------------------- entry
def read_ts(
    path: str | Path,
    format: str | None = None,
    **kws,
) -> TSData:
    """
    Read a magnetotelluric field time-series file.

    Parameters
    ----------
    path : str or Path
        Input file.
    format : {'lims', 'emslab', 'edi', 'ascii'}, optional
        Force a reader.  When omitted the format is sniffed
        with :func:`sniff_format`.
    **kws
        Forwarded to the concrete reader (e.g. ``dt=`` and
        ``chan=`` for ``ascii``).

    Returns
    -------
    TSData
        Parsed multichannel time series.

    Examples
    --------
    >>> from pycsamt.ts import read_ts
    >>> ts = read_ts("data/MT/TS/kap103as.ts/kap103as.ts")
    >>> ts.channels()
    ['HX', 'HY', 'HZ', 'EX', 'EY']

    See Also
    --------
    read_lims, read_emslab, read_edi_tseries, read_ascii
    """
    fmt = (format or sniff_format(path)).strip().lower()
    if fmt == "lims":
        return read_lims(path, **kws)
    if fmt == "emslab":
        return read_emslab(path, **kws)
    if fmt == "edi":
        return read_edi_tseries(path, **kws)
    if fmt == "ascii":
        return read_ascii(path, **kws)
    raise FileHandlingError(
        f"Unknown time-series format {fmt!r}. Expected "
        "'lims', 'emslab', 'edi' or 'ascii'."
    )
