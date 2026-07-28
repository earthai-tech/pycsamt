# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Time series → impedance → EDI conversion.

High-level orchestration on top of
:mod:`pycsamt.ts.readers` (field-format IO) and
:mod:`pycsamt.ts.process` (spectral estimation):

* :func:`ts_to_z` — recover the impedance tensor (and tipper)
  from a record, by single-site least squares or remote
  reference.
* :func:`ts_to_edi` — write a complete SEG EDI file
  (``>HEAD``, ``>INFO``, ``>=DEFINEMEAS``, ``>=MTSECT``,
  ``>FREQ``, impedance/rho/phase and tipper blocks) whose
  header sections are synthesized from the time-series
  metadata.

Examples
--------
One-liner, LiMS field file to EDI::

    from pycsamt.ts import ts_to_edi

    out = ts_to_edi(
        "data/MT/TS/kap103as.ts/kap103as.ts",
        savepath="edi_from_ts",
        nfft=10240,
    )

Remote reference with a second simultaneous record::

    from pycsamt.ts import read_ts, ts_to_edi

    local = read_ts("kap103as.ts")
    remote = read_ts("kap112as.ts")
    out = ts_to_edi(local, remote=remote, estimator="rr")

See Also
--------
pycsamt.seg.spectra.Spectra.to_Z
    The least-squares estimator reused here.
pycsamt.seg.edi.EDIFile
    Writer that assembles the final file.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from ..exceptions import EdIDataError
from ..log.logger import get_logger
from ._base import TSData
from .process import ts_to_spectra

logger = get_logger(__name__)

__all__ = ["ts_to_z", "ts_to_edifile", "ts_to_edi"]

#: SEG-recommended channel index (10*site + C)
_CHAN_C = {"HX": 1, "HY": 2, "HZ": 3, "EX": 4, "EY": 5}


def _as_tsdata(
    source: TSData | str | Path,
    verbose: int = 0,
    reader_kws: dict | None = None,
) -> TSData:
    if isinstance(source, TSData):
        return source
    from .readers import read_ts

    return read_ts(source, verbose=verbose, **(reader_kws or {}))


# ------------------------------------------------- remote reference (RR)
def _rr_estimate(
    sp,
    *,
    ridge: float | None = None,
    estimate_error: bool = True,
):
    """
    Remote-reference impedance ``Z = S_ER inv(S_HR)`` (and
    tipper ``T = S_zR inv(S_HR)``) from a spectra stack whose
    channel list includes ``RHX``/``RHY``.

    Error bars reuse the least-squares Wishart formula with
    the local blocks — a mild approximation for RR.
    """
    from ..seg.ops import (
        effective_dof_from_meta,
        z_error_from_blocks,
    )
    from ..z.tipper import Tipper
    from ..z.z import Z

    kinds = [str(c).upper() for c in sp.chan_ids]
    idx = {c: k for k, c in enumerate(kinds)}
    need = ("EX", "EY", "HX", "HY", "RHX", "RHY")
    missing = [c for c in need if c not in idx]
    if missing:
        raise EdIDataError(
            "Remote-reference estimation needs channels " f"{need}; missing {missing}."
        )
    e = [idx["EX"], idx["EY"]]
    h = [idx["HX"], idx["HY"]]
    r = [idx["RHX"], idx["RHY"]]
    hz = idx.get("HZ")

    S = np.asarray(sp._S)
    nf = S.shape[0]
    z_arr = np.empty((nf, 2, 2), complex)
    tip_arr = np.empty((nf, 1, 2), complex) if hz is not None else None
    z_err = np.full((nf, 2, 2), np.nan, float) if estimate_error else None

    dof = effective_dof_from_meta(
        segnum=getattr(sp, "segnum", None),
        avgt=getattr(sp, "avgt", None),
        bw=getattr(sp, "bw", None),
    )

    I2 = np.eye(2)
    for k in range(nf):
        Shr = S[k][np.ix_(h, r)]
        if ridge:
            Shr = Shr + float(ridge) * I2
        try:
            inv_Shr = np.linalg.inv(Shr)
        except np.linalg.LinAlgError as exc:
            raise EdIDataError(f"S_HR singular at band {k}: {exc}") from exc
        z_arr[k] = S[k][np.ix_(e, r)] @ inv_Shr
        if tip_arr is not None:
            tip_arr[k, 0, :] = S[k][np.ix_([hz], r)] @ inv_Shr
        if z_err is not None and dof is not None:
            M_k = float(np.atleast_1d(dof)[k])
            z_err[k] = z_error_from_blocks(
                S_EE=S[k][np.ix_(e, e)],
                S_EH=S[k][np.ix_(e, h)],
                S_HH=S[k][np.ix_(h, h)],
                M=M_k,
                ridge=ridge,
            )

    if z_err is not None and not (np.all(np.isfinite(z_err)) and np.all(z_err >= 0)):
        z_err = None

    freq = np.asarray(sp.freq, float)
    nm = getattr(sp, "name", None)
    z_obj = (
        Z(z_array=z_arr, freq=freq, name=nm, z_err_array=z_err)
        if z_err is not None
        else Z(z_array=z_arr, freq=freq, name=nm)
    )
    tip_obj = None
    if tip_arr is not None:
        tip_obj = Tipper()
        tip_obj._freq = freq.copy()
        tip_obj._tipper = tip_arr
        tip_obj.compute_amp_phase()
        tip_obj.compute_mag_direction()
    return z_obj, tip_obj


# ---------------------------------------------------------------- ts→Z
def ts_to_z(
    source: TSData | str | Path,
    *,
    estimator: str = "ls",
    ridge: float | None = None,
    estimate_error: bool = True,
    reader_kws: dict | None = None,
    verbose: int = 0,
    **spectra_kws,
):
    """
    Estimate impedance (and tipper) from a field record.

    Parameters
    ----------
    source : TSData or path
        Record or field file (loaded with
        :func:`pycsamt.ts.readers.read_ts`).
    estimator : {'ls', 'rr'}, default 'ls'
        ``'ls'``: single-site least squares
        ``Z = S_EH inv(S_HH)`` via
        :meth:`pycsamt.seg.spectra.Spectra.to_Z`.
        ``'rr'``: remote reference ``Z = S_ER inv(S_HR)``
        (requires ``remote=`` in ``spectra_kws``).
    ridge : float, optional
        Tikhonov stabilization of the inverted block.
    estimate_error : bool, default True
        Attach 1-sigma errors (DoF from the band averaging).
    reader_kws : dict, optional
        Extra options for
        :func:`pycsamt.ts.readers.read_ts` when ``source``
        is a path.
    verbose : int, default 0
        Verbosity.
    **spectra_kws
        Forwarded to :func:`pycsamt.ts.process.ts_to_spectra`
        (``nfft``, ``overlap``, ``per_decade``, ``fmin``,
        ``fmax``, ``freqs``, ``remote``, ``robust``, ...).

    Returns
    -------
    (Z, Tipper or None, Spectra)
        Impedance object, tipper when HZ exists, and the
        spectra stack used for the estimate.
    """
    ts = _as_tsdata(source, verbose=verbose, reader_kws=reader_kws)
    sp = ts_to_spectra(ts, verbose=verbose, **spectra_kws)
    est = str(estimator or "ls").strip().lower()
    if est == "rr":
        z_obj, tip = _rr_estimate(
            sp,
            ridge=ridge,
            estimate_error=estimate_error,
        )
    elif est == "ls":
        z_obj, tip = sp.to_Z(
            ridge=ridge,
            estimate_error=estimate_error,
        )
    else:
        raise EdIDataError(f"Unknown estimator {estimator!r}; use 'ls' or 'rr'.")
    return z_obj, tip, sp


# ----------------------------------------------------- header synthesis
def _edi_date(stamp: str | None) -> str | None:
    """``YYYY-MM-DD hh:mm:ss`` → EDI ``MM/DD/YY`` date."""
    if not stamp:
        return None
    s = str(stamp).strip()
    try:
        y, m, d = s[:10].split("-")
        return f"{int(m):02d}/{int(d):02d}/{y[2:4]}"
    except Exception:
        return s


def _meas_ids(ts: TSData, site_num: int = 1) -> dict[str, str]:
    """SEG-style measurement ids ``(10*site + C).run``."""
    ids: dict[str, str] = {}
    extra = 6
    for cid in ts.channels():
        c = _CHAN_C.get(cid)
        if c is None:
            c, extra = extra, extra + 1
        ids[cid] = f"{10 * site_num + c:04.0f}.001"
    return ids


def _make_head(ts: TSData, station: str):
    from ..seg.heads import Head

    kw = dict(
        dataid=station,
        acqdate=_edi_date(ts.start),
        enddate=_edi_date(ts.stop),
        fileby="pycsamt.ts",
        empty=1.0e32,
    )
    if ts.lat is not None:
        kw["lat"] = float(ts.lat)
    if ts.lon is not None:
        kw["long"] = float(ts.lon)
    if ts.elev is not None:
        kw["elev"] = float(ts.elev)
    if ts.declination is not None:
        kw["declination"] = float(ts.declination)
    if ts.coordsys:
        kw["coordsys"] = str(ts.coordsys).title()
    instrument = ts.meta.get("instrument")
    if instrument:
        kw["acqby"] = f"LiMS-{instrument}"
    return Head(**{k: v for k, v in kw.items() if v is not None})


def _make_info(ts: TSData, sp, estimator: str):
    from ..seg.heads import Info

    m = getattr(sp, "meta", {}) or {}
    freq = np.asarray(sp.freq, float)
    lines = [
        "Impedance estimated from field time series by pycsamt.ts",
        f"source file      : {m.get('source')}",
        f"station          : {ts.station}",
        f"sampling dt      : {ts.dt} s",
        f"samples          : {ts.n_samples}",
        f"channels         : {','.join(ts.channels())}",
        f"segments         : {m.get('n_segments')} "
        f"(nfft={m.get('nfft')}, "
        f"overlap={m.get('overlap')})",
        f"bands            : {freq.size} in "
        f"[{freq.min():.4E}, {freq.max():.4E}] Hz",
        "stacking         : Huber-weighted section averages",
        f"estimator        : {estimator.upper()} "
        + ("(Z = S_ER inv(S_HR))" if estimator == "rr" else "(Z = S_EH inv(S_HH))"),
    ]
    if ts.start or ts.stop:
        lines.insert(3, f"acquisition      : {ts.start} -> {ts.stop}")
    return Info(info_text=lines, maxinfo=999)


def _make_definemeas(ts: TSData, ids: dict[str, str]):
    from ..seg.meas import (
        DefineMeas,
        Emeasurement,
        Hmeasurement,
    )

    dm = DefineMeas(
        maxchan=ts.n_chan,
        maxrun=999,
        maxmeas=99,
        units="M",
        reftype="CART",
    )
    if ts.lat is not None:
        dm.reflat = float(ts.lat)
    if ts.lon is not None:
        dm.reflong = float(ts.lon)
    if ts.elev is not None:
        dm.refelev = float(ts.elev)

    for cid in ts.channels():
        az = float(ts.azim.get(cid, 0.0) or 0.0)
        if cid.startswith(("H", "RH")):
            dm.hmeas.append(
                Hmeasurement(
                    id=ids[cid],
                    chtype=cid,
                    x=0.0,
                    y=0.0,
                    azm=az,
                    acqchan=f"CH{_CHAN_C.get(cid, 0)}",
                    sensor=ts.sensor.get(cid),
                    gain=ts.gain.get(cid),
                )
            )
        else:  # electric dipole
            L = float(ts.dipole.get(cid, 0.0) or 0.0)
            rad = math.radians(az)
            dx = 0.5 * L * math.cos(rad)
            dy = 0.5 * L * math.sin(rad)
            dm.emeas.append(
                Emeasurement(
                    id=ids[cid],
                    chtype=cid,
                    x=round(-dx, 2),
                    y=round(-dy, 2),
                    x2=round(dx, 2),
                    y2=round(dy, 2),
                    acqchan=f"CH{_CHAN_C.get(cid, 0)}",
                    gain=ts.gain.get(cid),
                )
            )
    return dm


def _make_mtsect(station: str, nfreq: int, ids: dict[str, str]):
    from ..seg.mtemap import MTEMAP

    kw = {c.lower(): ids[c] for c in ("HX", "HY", "HZ", "EX", "EY") if c in ids}
    for c in ("RHX", "RHY"):
        if c in ids:
            kw["r" + c[-1].lower()] = ids[c]
    return MTEMAP(sectid=str(station), nfreq=int(nfreq), **kw)


# --------------------------------------------------------------- ts→EDI
def ts_to_edifile(
    source: TSData | str | Path,
    *,
    station: str | None = None,
    estimator: str = "ls",
    ridge: float | None = None,
    estimate_error: bool = True,
    include_spectra: bool = False,
    include_tseries: bool = False,
    tseries_max_samples: int = 20000,
    reader_kws: dict | None = None,
    verbose: int = 0,
    **spectra_kws,
):
    """
    Build an in-memory impedance :class:`EDIFile` from a
    field time series (no disk write).

    This is the assembly core of :func:`ts_to_edi`: it runs
    the spectral estimation, recovers ``Z`` (and the tipper),
    synthesizes the EDI header sections from the time-series
    metadata, and returns the populated
    :class:`~pycsamt.seg.edi.EDIFile`.  Call
    :meth:`~pycsamt.seg.edi.EDIFile.write` yourself — or use
    :func:`ts_to_edi` — to persist it.

    Parameters
    ----------
    source, station, estimator, ridge, estimate_error, \
include_spectra, include_tseries, tseries_max_samples, \
reader_kws, verbose, **spectra_kws
        Same meaning as in :func:`ts_to_edi`.

    Returns
    -------
    EDIFile
        Fully populated container (``head``, ``info``,
        ``definemeas``, ``mtsect`` sections, ``Z`` and, when
        available, ``Tip``), not yet written to disk.

    Examples
    --------
    >>> from pycsamt.ts import ts_to_edifile
    >>> ed = ts_to_edifile("kap103as.ts", nfft=10240)
    >>> ed.Z.n_freq > 0
    True
    >>> out = ed.write(new_edifn="kap103_ts.edi",
    ...                savepath="edi_out")

    See Also
    --------
    ts_to_edi
        Build **and** write in one call.
    """
    from ..seg.edi import EDIFile

    ts = _as_tsdata(source, verbose=verbose, reader_kws=reader_kws)
    sid = str(station or ts.station or ts.name or "SITE")

    z_obj, tip, sp = ts_to_z(
        ts,
        estimator=estimator,
        ridge=ridge,
        estimate_error=estimate_error,
        verbose=verbose,
        **spectra_kws,
    )
    z_obj.name = sid

    ids = _meas_ids(ts)
    ed = EDIFile(verbose=verbose)
    ed.add_section("head", _make_head(ts, sid))
    ed.add_section("info", _make_info(ts, sp, estimator))
    ed.add_section("definemeas", _make_definemeas(ts, ids))
    ed.add_section(
        "mtsect",
        _make_mtsect(sid, z_obj.n_freq, ids),
    )
    ed.Z = z_obj
    if tip is not None:
        ed.Tip = tip

    if include_spectra:
        ed.add_section("spectra", sp)
    if include_tseries:
        ed.add_section(
            "timeseries",
            ts.to_seg_timeseries(max_samples=tseries_max_samples),
        )
    return ed


def ts_to_edi(
    source: TSData | str | Path,
    out: str | None = None,
    *,
    savepath: str | Path | None = None,
    station: str | None = None,
    estimator: str = "ls",
    ridge: float | None = None,
    estimate_error: bool = True,
    include_spectra: bool = False,
    include_tseries: bool = False,
    tseries_max_samples: int = 20000,
    reader_kws: dict | None = None,
    verbose: int = 0,
    **spectra_kws,
) -> str:
    """
    Convert a field time series to a complete impedance EDI.

    The output file carries ``>HEAD``, ``>INFO``,
    ``>=DEFINEMEAS`` (with ``>HMEAS``/``>EMEAS`` built from
    the recorded azimuths and dipole lengths), ``>=MTSECT``,
    ``>FREQ``, the impedance blocks with 1-sigma variances,
    apparent resistivity / phase blocks, and — when a
    vertical magnetic channel exists — the tipper blocks.

    Parameters
    ----------
    source : TSData or path
        Field record (``lims``, ``emslab``, ``edi`` or
        ``ascii`` formats; see
        :func:`pycsamt.ts.readers.read_ts`).
    out : str, optional
        Output file name (``.edi`` appended when missing).
        Defaults to ``"<station>_ts.edi"``.
    savepath : str or Path, optional
        Output directory (created if needed).  Defaults to
        ``./edi_out`` (the :meth:`EDIFile.write` default).
    station : str, optional
        Override the station id used for ``DATAID`` and
        ``SECTID``.
    estimator : {'ls', 'rr'}, default 'ls'
        See :func:`ts_to_z`.
    ridge : float, optional
        Stabilization forwarded to the estimator.
    estimate_error : bool, default True
        Attach impedance variances to the EDI.
    include_spectra : bool, default False
        Also embed the estimated cross-spectra as
        ``>=SPECTRASECT`` / ``>SPECTRA`` blocks.
    include_tseries : bool, default False
        Also embed (a decimated head of) the raw series as
        ``>=TSERIESSECT`` / ``>TSERIES`` blocks.
    tseries_max_samples : int, default 20000
        Per-channel cap used when ``include_tseries``.
    reader_kws : dict, optional
        Extra options for
        :func:`pycsamt.ts.readers.read_ts` when ``source``
        is a path (e.g. ``{"declination": -19.5}``).
    verbose : int, default 0
        Verbosity.
    **spectra_kws
        Estimation options forwarded to
        :func:`pycsamt.ts.process.ts_to_spectra`
        (``nfft``, ``overlap``, ``per_decade``, ``fmin``,
        ``fmax``, ``freqs``, ``remote``, ``robust``, ...).

    Returns
    -------
    str
        Path of the written EDI file.

    Examples
    --------
    >>> from pycsamt.ts import ts_to_edi
    >>> out = ts_to_edi(
    ...     "data/MT/TS/kap103as.ts/kap103as.ts",
    ...     savepath="edi_from_ts",
    ...     nfft=10240,
    ... )
    >>> from pycsamt.seg.edi import EDIFile
    >>> ed = EDIFile(out)  # round-trips as regular EDI
    >>> ed.Z.n_freq > 0
    True
    """
    ed = ts_to_edifile(
        source,
        station=station,
        estimator=estimator,
        ridge=ridge,
        estimate_error=estimate_error,
        include_spectra=include_spectra,
        include_tseries=include_tseries,
        tseries_max_samples=tseries_max_samples,
        reader_kws=reader_kws,
        verbose=verbose,
        **spectra_kws,
    )
    sid = ed.station or "SITE"
    name = out or f"{sid}_ts.edi"
    path = ed.write(new_edifn=name, savepath=savepath)
    if verbose:
        logger.info("EDI written: %s", path)
    return path
