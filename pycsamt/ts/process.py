# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Spectral estimation for MT field time series.

This module turns a :class:`pycsamt.ts.TSData` record into
band-averaged cross-power spectral matrices — the quantity
stored in EDI ``>SPECTRA`` blocks and consumed by
:meth:`pycsamt.seg.spectra.Spectra.to_Z` to recover the
impedance tensor and tipper.

Processing chain
----------------
1. Small interior gaps are interpolated; segments still
   containing missing samples are dropped.
2. The record is cut into ``nfft``-sample segments with
   fractional ``overlap``; each segment is linearly detrended
   and tapered (Hann window).
3. Segment FFTs are turned into per-segment cross-power
   matrices, averaged over the Fourier bins of each target
   frequency band (log-spaced bands, Welch scaling).
4. Segments are stacked with optional Huber-style robust
   weights computed from the misfit of the predicted electric
   field (``E - Z H`` residual power), which strongly
   downweights noise bursts and non-plane-wave intervals.

The result is a genuine :class:`pycsamt.seg.spectra.Spectra`
object, so everything downstream (``to_Z``, ``to_edi``,
plotting in :mod:`pycsamt.emtools.spectra`) applies untouched.

References
----------
.. [1] Bendat, J. S., & Piersol, A. G. (2011). *Random Data:
       Analysis and Measurement Procedures*. Wiley.
.. [2] Chave, A. D., & Jones, A. G. (2012). *The
       Magnetotelluric Method: Theory and Practice*. CUP.
.. [3] Chave, A. D., & Thomson, D. J. (2004). Bounded
       influence magnetotelluric response function
       estimation. *GJI*, 157(3), 988-1006.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..exceptions import EdIDataError
from ..log.logger import get_logger
from ._base import CHANNEL_ORDER, TSData

logger = get_logger(__name__)

__all__ = [
    "preprocess",
    "target_frequencies",
    "cross_spectra",
    "ts_to_spectra",
]


# ------------------------------------------------------------ preprocess
def _fill_small_gaps(
    x: np.ndarray, max_gap: int
) -> np.ndarray:
    """
    Linearly interpolate interior ``NaN`` runs of length
    ``<= max_gap`` samples.  Longer runs (and leading/trailing
    gaps) are left as ``NaN`` so the segmentation step can
    reject them.
    """
    x = np.asarray(x, float).copy()
    isnan = np.isnan(x)
    if not isnan.any() or max_gap <= 0:
        return x
    n = x.size
    idx = np.flatnonzero(isnan)
    # split indices into contiguous runs
    splits = np.flatnonzero(np.diff(idx) > 1) + 1
    for run in np.split(idx, splits):
        if run.size == 0 or run.size > max_gap:
            continue
        lo, hi = run[0] - 1, run[-1] + 1
        if lo < 0 or hi >= n:
            continue  # edge gap: leave missing
        x[run] = np.interp(
            run.astype(float),
            [float(lo), float(hi)],
            [x[lo], x[hi]],
        )
    return x


def preprocess(
    ts: TSData,
    *,
    fill_gap: int = 5,
    remove_mean: bool = True,
) -> TSData:
    """
    Light global conditioning of a :class:`TSData`.

    Parameters
    ----------
    ts : TSData
        Input record.
    fill_gap : int, default 5
        Maximum interior gap (in samples) repaired by linear
        interpolation.  Longer gaps stay ``NaN`` and the
        affected segments are dropped later.
    remove_mean : bool, default True
        Subtract the (finite-sample) mean of every channel.
        Per-segment linear detrending is applied later in
        :func:`cross_spectra` regardless.

    Returns
    -------
    TSData
        A new conditioned record (input left untouched).
    """
    out = ts.copy_meta()
    for cid in ts.channels():
        x = _fill_small_gaps(ts.get(cid), int(fill_gap))
        if remove_mean:
            mu = np.nanmean(x)
            if np.isfinite(mu):
                x = x - mu
        out.add_channel(cid, x)
    return out


# ------------------------------------------------------- frequency grid
def target_frequencies(
    dt: float,
    nfft: int,
    *,
    per_decade: int = 7,
    fmin: Optional[float] = None,
    fmax: Optional[float] = None,
) -> np.ndarray:
    """
    Log-spaced target frequencies for band averaging.

    Parameters
    ----------
    dt : float
        Sampling interval (s).
    nfft : int
        FFT segment length (samples).
    per_decade : int, default 7
        Number of estimation frequencies per decade.
    fmin : float, optional
        Lowest target frequency.  Defaults to ``4 / (nfft*dt)``
        (at least four Fourier bins below the first target).
    fmax : float, optional
        Highest target frequency.  Defaults to a quarter of
        the sampling rate (half the Nyquist frequency), a
        conservative anti-alias margin.

    Returns
    -------
    ndarray
        Frequencies in Hz, **descending** (the normalized
        pycsamt order).
    """
    fs = 1.0 / float(dt)
    lo = float(fmin) if fmin else 4.0 / (nfft * float(dt))
    hi = float(fmax) if fmax else fs / 4.0
    if not (0.0 < lo < hi <= fs / 2.0):
        raise EdIDataError(
            f"Invalid frequency range [{lo}, {hi}] for "
            f"fs={fs} Hz."
        )
    ndec = np.log10(hi / lo)
    n = max(2, int(round(ndec * per_decade)) + 1)
    freqs = np.geomspace(hi, lo, n)  # descending
    return freqs


# --------------------------------------------------------- cross spectra
def _segment_starts(
    n: int, nfft: int, step: int
) -> np.ndarray:
    if n < nfft:
        return np.empty(0, dtype=int)
    return np.arange(0, n - nfft + 1, step, dtype=int)


def _band_edges(freqs_desc: np.ndarray) -> np.ndarray:
    """
    Geometric band edges around each (descending) target
    frequency: ``edges[j], edges[j+1]`` bracket ``freqs[j]``.
    """
    f = np.asarray(freqs_desc, float)
    mids = np.sqrt(f[:-1] * f[1:])
    ratio = (
        np.sqrt(f[0] / f[1]) if f.size > 1 else np.sqrt(2.0)
    )
    hi = f[0] * ratio
    lo = f[-1] / ratio
    return np.concatenate([[hi], mids, [lo]])


def cross_spectra(
    ts: TSData,
    *,
    nfft: Optional[int] = None,
    overlap: float = 0.5,
    freqs: Optional[Sequence[float]] = None,
    per_decade: int = 7,
    fmin: Optional[float] = None,
    fmax: Optional[float] = None,
    remote: Optional[TSData] = None,
    robust: Optional[str] = "huber",
    n_iter: int = 2,
    huber_k: float = 1.5,
    min_segments: int = 2,
    fill_gap: int = 5,
    verbose: int = 0,
) -> Dict[str, np.ndarray]:
    """
    Band-averaged cross-power spectral matrices of a record.

    Parameters
    ----------
    ts : TSData
        Local site record (requires ``dt``).  Channels are
        used in canonical order HX, HY, HZ, EX, EY (whatever
        subset exists), followed by any extra channels.
    nfft : int, optional
        Segment length in samples.  Default: largest power of
        two ``<= n_samples / 32`` (≈ 60+ segments at 50 %
        overlap), capped at 65536 and floored at 256.
    overlap : float, default 0.5
        Fractional segment overlap in ``[0, 0.9]``.
    freqs : sequence of float, optional
        Explicit target frequencies (Hz).  When omitted a
        log-spaced grid is built with
        :func:`target_frequencies`.
    per_decade, fmin, fmax
        Grid parameters forwarded to
        :func:`target_frequencies`.
    remote : TSData, optional
        Simultaneous remote-reference record.  Its horizontal
        magnetic channels are appended as ``RHX``/``RHY``
        (same ``dt`` required; the overlapping head of both
        records is used).
    robust : {'huber', None}, default 'huber'
        Segment-stacking scheme.  ``'huber'`` iteratively
        downweights segments whose electric field is poorly
        predicted by the current impedance estimate;
        ``None`` is a plain average.
    n_iter : int, default 2
        Robust reweighting iterations.
    huber_k : float, default 1.5
        Huber tuning constant (in units of the median segment
        residual power).
    min_segments : int, default 2
        Bands averaged over fewer valid segments are dropped.
    fill_gap : int, default 5
        Forwarded to :func:`preprocess` when the record still
        contains ``NaN``.
    verbose : int, default 0
        Verbosity level.

    Returns
    -------
    dict
        Keys: ``freq`` (nf, descending), ``S`` (nf, nc, nc)
        Hermitian cross-power matrices (Welch scaling,
        units²/Hz), ``chan_ids`` (nc,), ``segnum`` (nf,)
        effective number of independent averages, ``bw``
        (nf,) band width in Hz, ``avgt`` (nf,) averaged time
        in seconds, ``n_segments`` (scalar), ``weights``
        (n_segments,) final stacking weights.

    Notes
    -----
    Per segment ``s`` and band ``j`` the raw estimate is

    .. math::

       S_{sj} = \\frac{2\\,\\Delta t}{N\\,\\bar{w^2}}\\;
                \\langle X_s(f_b) X_s(f_b)^H
                \\rangle_{b \\in j}

    and the stacked matrix is the (weighted) mean over
    segments.  The Huber weights are computed from the
    residual power ``tr(S_EE - Z S_HE - (Z S_HE)^H +
    Z S_HH Z^H)`` summed over bands, so one weight applies
    per segment (a noisy burst pollutes all its bands).
    """
    if ts.dt is None:
        raise EdIDataError(
            "TSData.dt is not set; cannot process."
        )
    dt = float(ts.dt)

    # ---- channel order: canonical first, extras after
    order = [c for c in CHANNEL_ORDER if c in ts]
    order += [c for c in ts.channels() if c not in order]

    M, order = ts.matrix(order)

    # ---- optional remote reference channels
    if remote is not None:
        if remote.dt is None or float(remote.dt) != dt:
            raise EdIDataError(
                "remote record must share the same dt."
            )
        n_common = min(
            M.shape[0], remote.n_samples
        )
        M = M[:n_common]
        cols = []
        for src in ("HX", "HY"):
            if src in remote:
                cols.append(
                    remote.get(src)[:n_common]
                )
                order.append("RH" + src[-1])
        if cols:
            M = np.column_stack([M] + cols)

    n, nc = M.shape

    # ---- fill small gaps (idempotent if already done)
    if np.isnan(M).any():
        for j in range(nc):
            M[:, j] = _fill_small_gaps(
                M[:, j], int(fill_gap)
            )

    # ---- segmentation: default targets ~60+ segments at
    # 50% overlap, a robust averaging/low-frequency tradeoff
    if nfft is None:
        nfft = int(
            2 ** np.floor(np.log2(max(n // 32, 256)))
        )
        nfft = int(min(nfft, 65536))
    nfft = int(nfft)
    if nfft > n:
        raise EdIDataError(
            f"nfft={nfft} exceeds record length {n}."
        )
    overlap = float(np.clip(overlap, 0.0, 0.9))
    step = max(1, int(round(nfft * (1.0 - overlap))))
    starts = _segment_starts(n, nfft, step)

    # ---- frequency grid and band -> bin mapping
    fbin = np.fft.rfftfreq(nfft, dt)
    if freqs is None:
        fgrid = target_frequencies(
            dt,
            nfft,
            per_decade=per_decade,
            fmin=fmin,
            fmax=fmax,
        )
    else:
        fgrid = np.sort(
            np.asarray(list(freqs), float)
        )[::-1]
    edges = _band_edges(fgrid)

    band_bins: List[np.ndarray] = []
    keep: List[int] = []
    for j in range(fgrid.size):
        hi, lo = edges[j], edges[j + 1]
        sel = np.flatnonzero(
            (fbin > 0.0) & (fbin > lo) & (fbin <= hi)
        )
        if sel.size:
            band_bins.append(sel)
            keep.append(j)
    if not keep:
        raise EdIDataError(
            "No Fourier bins fall inside the requested "
            "bands; increase nfft or widen the range."
        )
    fgrid = fgrid[keep]
    nf = fgrid.size

    # ---- window and Welch scale
    w = np.hanning(nfft)
    scale = 2.0 * dt / np.sum(w * w)

    # ---- per-segment band cross-powers
    seg_S: List[np.ndarray] = []
    used_starts: List[int] = []
    tc = np.arange(nfft, dtype=float)
    tc -= tc.mean()
    tvar = float(tc @ tc)

    for s0 in starts:
        seg = M[s0 : s0 + nfft, :]
        if not np.isfinite(seg).all():
            continue
        # linear detrend (vectorized least squares)
        mu = seg.mean(axis=0)
        slope = ((seg - mu).T @ tc) / tvar
        seg = seg - mu - np.outer(tc, slope)
        X = np.fft.rfft(seg * w[:, None], axis=0)

        # S[i, k] = <x_i x_k^*>  (Hermitian, matches the
        # E = Z H convention used by Spectra.to_Z)
        Sj = np.empty((nf, nc, nc), dtype=complex)
        for j, sel in enumerate(band_bins):
            xb = X[sel, :]  # (nb, nc)
            Sj[j] = scale * (
                xb.T @ np.conj(xb)
            ) / float(sel.size)
        seg_S.append(Sj)
        used_starts.append(int(s0))

    n_seg = len(seg_S)
    if n_seg < max(1, int(min_segments)):
        raise EdIDataError(
            f"Only {n_seg} valid segments (min "
            f"{min_segments}); record too short or too "
            "gappy for nfft="
            f"{nfft}."
        )
    A = np.asarray(seg_S)  # (ns, nf, nc, nc)

    # ---- stacking (plain or Huber-reweighted)
    weights = np.ones(n_seg, float)
    S = A.mean(axis=0)

    idx = {c: k for k, c in enumerate(order)}
    can_z = all(
        c in idx for c in ("HX", "HY", "EX", "EY")
    )
    if robust and str(robust).lower() == "huber" and can_z:
        h = [idx["HX"], idx["HY"]]
        e = [idx["EX"], idx["EY"]]
        for _ in range(max(1, int(n_iter))):
            # current Z per band from the stacked matrices
            Z = np.empty((nf, 2, 2), complex)
            ok = np.ones(nf, bool)
            for j in range(nf):
                Shh = S[j][np.ix_(h, h)]
                Seh = S[j][np.ix_(e, h)]
                try:
                    Z[j] = Seh @ np.linalg.inv(Shh)
                except np.linalg.LinAlgError:
                    ok[j] = False
                    Z[j] = 0.0
            # residual power per segment (summed over bands)
            r = np.zeros(n_seg, float)
            for j in np.flatnonzero(ok):
                See = A[:, j][:, e][:, :, e]
                She = A[:, j][:, h][:, :, e]
                Shh = A[:, j][:, h][:, :, h]
                Zj = Z[j]
                pred = Zj @ She  # (ns, 2, 2) E x E cross
                res = (
                    See
                    - pred
                    - np.conj(
                        np.swapaxes(pred, -1, -2)
                    )
                    + Zj @ Shh @ np.conj(Zj.T)
                )
                tr = np.trace(
                    res, axis1=-2, axis2=-1
                ).real
                # normalize by band E power -> comparable
                pe = np.trace(
                    See, axis1=-2, axis2=-1
                ).real
                good = pe > 0
                r[good] += np.clip(
                    tr[good] / pe[good], 0.0, None
                )
            med = float(np.median(r[r > 0])) if (
                r > 0
            ).any() else 0.0
            if med <= 0:
                break
            weights = np.minimum(
                1.0, huber_k * med / np.maximum(r, 1e-30)
            )
            wsum = float(weights.sum())
            if wsum <= 0:
                weights = np.ones(n_seg, float)
                break
            S = np.tensordot(
                weights, A, axes=(0, 0)
            ) / wsum

    # enforce Hermitian symmetry (numerical hygiene)
    S = 0.5 * (S + np.conj(np.swapaxes(S, -1, -2)))

    # ---- effective DoF per band
    neff = (
        float(weights.sum()) ** 2
        / float((weights**2).sum())
    )
    nb = np.array([b.size for b in band_bins], float)
    segnum = np.maximum(
        1, np.rint(neff * nb)
    ).astype(int)
    bw = np.array(
        [
            fbin[b].max() - fbin[b].min()
            + (fbin[1] - fbin[0])
            for b in band_bins
        ],
        float,
    )
    avgt = np.full(nf, n_seg * nfft * dt, float)

    if verbose:
        logger.info(
            "cross_spectra: %d segments (nfft=%d, "
            "overlap=%.0f%%), %d bands "
            "[%.3e, %.3e] Hz, neff=%.1f",
            n_seg, nfft, 100 * overlap, nf,
            fgrid.min(), fgrid.max(), neff,
        )

    return {
        "freq": fgrid,
        "S": S,
        "chan_ids": list(order),
        "segnum": segnum,
        "bw": bw,
        "avgt": avgt,
        "n_segments": n_seg,
        "weights": weights,
        "nfft": nfft,
        "overlap": overlap,
    }


# ------------------------------------------------------------- spectra
def ts_to_spectra(
    ts: Union[TSData, str],
    *,
    name: Optional[str] = None,
    preprocess_first: bool = True,
    fill_gap: int = 5,
    reader_kws: Optional[dict] = None,
    verbose: int = 0,
    **kws,
):
    """
    Estimate a :class:`pycsamt.seg.spectra.Spectra` from a
    field time series.

    Parameters
    ----------
    ts : TSData or str
        Record or path (paths are loaded with
        :func:`pycsamt.ts.readers.read_ts`).
    name : str, optional
        Display name (defaults to the station id).
    preprocess_first : bool, default True
        Apply :func:`preprocess` before estimation.
    fill_gap : int, default 5
        Gap-repair length forwarded to :func:`preprocess`.
    reader_kws : dict, optional
        Extra options for
        :func:`pycsamt.ts.readers.read_ts` when ``ts`` is a
        path (e.g. ``{"declination": -19.5}`` for EMSLAB).
    verbose : int, default 0
        Verbosity level.
    **kws
        Estimation options forwarded to
        :func:`cross_spectra` (``nfft``, ``overlap``,
        ``freqs``, ``per_decade``, ``fmin``, ``fmax``,
        ``remote``, ``robust``, ...).

    Returns
    -------
    pycsamt.seg.spectra.Spectra
        Populated spectra container (``_freq``, ``_S``,
        ``chan_ids``, ``bw``, ``avgt``, ``segnum``) ready for
        :meth:`~pycsamt.seg.spectra.Spectra.to_Z` /
        :meth:`~pycsamt.seg.spectra.Spectra.to_edi`.

    Examples
    --------
    >>> from pycsamt.ts import read_ts, ts_to_spectra
    >>> rec = read_ts("data/MT/TS/kap103as.ts/kap103as.ts")
    >>> sp = ts_to_spectra(rec, nfft=8192)
    >>> Zhat, tip = sp.to_Z(estimate_error=True)

    See Also
    --------
    cross_spectra
        The estimation engine (returns raw arrays).
    pycsamt.seg.spectra.Spectra.to_Z
        Impedance/tipper recovery from the result.
    """
    from ..seg.spectra import Spectra

    if isinstance(ts, (str,)) or hasattr(ts, "__fspath__"):
        from .readers import read_ts as _read

        ts = _read(
            ts, verbose=verbose, **(reader_kws or {})
        )

    rec = (
        preprocess(ts, fill_gap=fill_gap)
        if preprocess_first
        else ts
    )
    est = cross_spectra(
        rec, fill_gap=fill_gap, verbose=verbose, **kws
    )

    sp = Spectra(
        name=name or ts.station or ts.name,
        verbose=verbose,
    )
    sp._freq = np.asarray(est["freq"], float)
    sp._S = np.asarray(est["S"], complex)
    sp.chan_ids = list(est["chan_ids"])

    nf = sp._freq.size
    sp.bw = np.asarray(est["bw"], float)
    sp.avgt = np.asarray(est["avgt"], float)
    sp.avgf = np.full(nf, np.nan, float)
    sp.rotspec = np.zeros(nf, float)
    sp.segnum = np.asarray(est["segnum"], int)
    sp.band = [
        f"{f:.4E}Hz" for f in np.asarray(est["freq"])
    ]
    # processing provenance for >INFO synthesis downstream
    sp.meta = {
        "n_segments": est["n_segments"],
        "nfft": est["nfft"],
        "overlap": est["overlap"],
        "source": ts.meta.get("source"),
        "dt": ts.dt,
    }
    return sp
