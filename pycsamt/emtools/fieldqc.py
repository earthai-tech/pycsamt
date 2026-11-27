# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Field-data quality control for noisy MT/AMT surveys.

Tools distilled from processing a noisy AMT survey whose EDIs had
station-to-station polarity reversals, a duplicated station, a
systematically corrupted frequency band and strong static shifts.
Every check works on plain arrays; the ``*_table`` wrappers accept a
``Sites`` object or a sequence of ``EDIFile``.

Conventions
-----------
Impedances ``z`` have shape ``(n_freq, 2, 2)`` (complex, field units
or SI), errors ``z_err`` the same shape (standard deviation of
:math:`|Z_{ij}|`), frequencies in Hz. The sign convention is
:math:`e^{+i\omega t}`, for which a 1-D or 2-D earth gives

.. math::

    \phi_{xy} = \arg Z_{xy} \in [0, 90^\circ], \qquad
    \phi_{yx} = \arg Z_{yx} \in [-180^\circ, -90^\circ].

Polarity
--------
In a reference band :math:`B` the mean phase direction of a mode is
:math:`\bar\phi = \arg \sum_{k \in B} e^{i\phi_k}` with resultant
length :math:`R = |\sum_{k\in B} e^{i\phi_k}| / |B|` (1 = perfectly
consistent phase). A mode whose :math:`\bar\phi` is closer to the
opposite quadrant than to the expected one is reversed. A channel
sign error maps onto the tensor as :math:`Z' = S_E Z S_H` with
:math:`S_E, S_H = \mathrm{diag}(\pm 1, \pm 1)`: a reversed electric
channel flips a *row*, a reversed magnetic channel a *column*. The
off-diagonal phases alone cannot tell the two apart; choose with
``assume`` in :func:`correct_polarity`.

Point flags
-----------
``err``: :math:`\delta|Z| / |Z| >` ``err_max``. ``quadrant``:
polarity-normalised phase outside :math:`[-t, 90 + t]`. ``spike``:
Hampel test along frequency on :math:`\log_{10}\rho_a` and
:math:`\phi`. ``neighbour``: :math:`|\phi - \mathrm{median}(\phi_{nb})|`
over the nearest stations exceeds a tolerance (phase is free of
galvanic static shift, so a large departure is local noise or 3-D).

Static shift
------------
For each station and mode, the shift is the median over a band of
:math:`\log_{10}\rho - \mathrm{median}_{nb}\log_{10}\rho` using only
QC-good points; estimates beyond ``max_abs`` decades are returned but
flagged as not applicable (they reflect bad data or neighbours rather
than galvanic distortion).
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

__all__ = [
    "circular_mean_deg",
    "detect_polarity",
    "correct_polarity",
    "polarity_table",
    "find_duplicate_stations",
    "hampel_mask",
    "point_qc_flags",
    "point_qc_table",
    "detect_band_break",
    "estimate_static_shift_neighbours",
]

_MODES = {"xy": (0, 1), "yx": (1, 0)}
_TARGET = {"xy": 45.0, "yx": -135.0}


# ---------------------------------------------------------------- utils
def _wrap180(a):
    return -((-np.asarray(a, float) + 180.0) % 360.0 - 180.0)


def circular_mean_deg(angles: np.ndarray) -> tuple[float, float]:
    """Mean direction (deg) and resultant length of angles (deg).

    NaN values are ignored. Returns ``(nan, 0.0)`` for no data.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.emtools.fieldqc import circular_mean_deg
    >>> m, r = circular_mean_deg(np.array([170.0, -170.0]))
    >>> round(abs(m), 6), round(r, 3)
    (180.0, 0.985)
    """
    a = np.asarray(angles, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan"), 0.0
    u = np.mean(np.exp(1j * np.radians(a)))
    return float(np.degrees(np.angle(u))), float(np.abs(u))


def _station_arrays(sites: Any):
    """(names, freq list, z list, z_err list) from Sites/EDIFiles."""
    from pathlib import Path

    from pycsamt.emtools._core import _iter_items

    names, fs, zs, es = [], [], [], []
    for i, it in enumerate(_iter_items(sites)):
        ed = getattr(it, "edi", it)
        zo = getattr(ed, "Z", None)
        if zo is None or getattr(zo, "z", None) is None:
            continue
        path = getattr(ed, "path", None)
        name = (Path(path).stem if path is not None
                else getattr(ed, "station", None) or f"S{i:03d}")
        names.append(str(name))
        fs.append(np.asarray(zo.freq, float))
        zs.append(np.asarray(zo.z, complex))
        e = getattr(zo, "z_err", None)
        es.append(np.full(zo.z.shape, np.nan) if e is None
                  else np.asarray(e, float))
    return names, fs, zs, es


# ------------------------------------------------------------ polarity
def detect_polarity(
    z: np.ndarray,
    freq: np.ndarray,
    *,
    band: tuple[float, float] = (3.0, 30.0),
    r_min: float = 0.6,
) -> dict[str, Any]:
    """Detect reversed impedance modes of one station.

    Parameters
    ----------
    z : ndarray of complex, shape (n_freq, 2, 2)
        Impedance tensor.
    freq : ndarray, shape (n_freq,)
        Frequencies (Hz).
    band : tuple of float, default (3.0, 30.0)
        Reference band (Hz) where the signal is coherent.
    r_min : float, default 0.6
        Minimum resultant length for the decision to be *trusted*.

    Returns
    -------
    dict
        For each mode ``m`` in ``("xy", "yx")``: ``flip_m`` (bool),
        ``phase_m`` (band mean phase, deg), ``R_m``, ``trusted_m``
        (``R >= r_min``) and ``anomalous_m`` (mean phase more than
        45 deg from both expected quadrants, i.e. in Q2/Q4).

    Raises
    ------
    ValueError
        If the band holds no frequency.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.emtools.fieldqc import detect_polarity
    >>> f = np.logspace(0, 2, 9)
    >>> z = np.zeros((9, 2, 2), complex)
    >>> z[:, 0, 1] = -(1 + 1j)          # reversed xy (Q3)
    >>> z[:, 1, 0] = -(1 + 1j)          # normal yx (Q3)
    >>> d = detect_polarity(z, f, band=(1, 100))
    >>> d["flip_xy"], d["flip_yx"]
    (True, False)
    """
    f = np.asarray(freq, float)
    sel = (f >= band[0]) & (f <= band[1])
    if not sel.any():
        raise ValueError(f"no frequency inside band {band}.")
    out: dict[str, Any] = {}
    for m, (i, j) in _MODES.items():
        ph = np.degrees(np.angle(np.asarray(z)[sel, i, j]))
        mean, r = circular_mean_deg(ph)
        if not np.isfinite(mean):
            out.update({f"flip_{m}": False, f"phase_{m}": np.nan,
                        f"R_{m}": 0.0, f"trusted_{m}": False,
                        f"anomalous_{m}": False})
            continue
        d_keep = abs(_wrap180(mean - _TARGET[m]))
        d_flip = abs(_wrap180(mean + 180.0 - _TARGET[m]))
        out.update({
            f"flip_{m}": bool(d_flip < d_keep),
            f"phase_{m}": mean,
            f"R_{m}": r,
            f"trusted_{m}": bool(r >= r_min),
            f"anomalous_{m}": bool(min(d_keep, d_flip) > 45.0),
        })
    return out


def correct_polarity(
    z: np.ndarray,
    *,
    flip_xy: bool,
    flip_yx: bool,
    assume: str = "E",
) -> np.ndarray:
    """Undo mode reversals by flipping channel signs.

    Parameters
    ----------
    z : ndarray of complex, shape (n_freq, 2, 2)
        Impedance tensor (not modified).
    flip_xy, flip_yx : bool
        Modes found reversed (see :func:`detect_polarity`).
    assume : {"E", "H"}, default "E"
        Channel blamed for a reversal. ``"E"``: electric, flips the
        row (``xy`` -> row x: Zxx, Zxy; ``yx`` -> row y). ``"H"``:
        magnetic, flips the column (``xy`` -> column y: Zxy, Zyy;
        ``yx`` -> column x). Both give the same off-diagonals; they
        differ for the diagonal terms and the phase-tensor skew
        sign.

    Returns
    -------
    ndarray
        Corrected copy. Errors are unchanged by sign flips.

    Raises
    ------
    ValueError
        If *assume* is not ``"E"`` or ``"H"``.
    """
    if assume not in ("E", "H"):
        raise ValueError("assume must be 'E' or 'H'.")
    out = np.array(z, dtype=complex, copy=True)
    if assume == "E":
        if flip_xy:
            out[:, 0, :] *= -1
        if flip_yx:
            out[:, 1, :] *= -1
    else:
        if flip_xy:
            out[:, :, 1] *= -1
        if flip_yx:
            out[:, :, 0] *= -1
    return out


def polarity_table(
    sites: Any,
    *,
    band: tuple[float, float] = (3.0, 30.0),
    r_min: float = 0.6,
) -> pd.DataFrame:
    """Run :func:`detect_polarity` on every station.

    Returns
    -------
    pandas.DataFrame
        One row per station (``station`` = file stem or name).
    """
    names, fs, zs, _ = _station_arrays(sites)
    rows = []
    for n, f, z in zip(names, fs, zs):
        r = detect_polarity(z, f, band=band, r_min=r_min)
        rows.append({"station": n, **r})
    return pd.DataFrame(rows)


# ---------------------------------------------------------- duplicates
def find_duplicate_stations(
    sites: Any,
    *,
    rtol: float = 1e-9,
) -> list[list[str]]:
    """Groups of stations carrying identical impedance data.

    Two stations are duplicates when their frequencies and all four
    impedance elements agree within *rtol* (NaN equal to NaN), a
    sign of a file copied under a second station name.

    Returns
    -------
    list of list of str
        Groups of two or more station names; empty when none.

    Examples
    --------
    >>> import numpy as np
    >>> from types import SimpleNamespace as NS
    >>> from pycsamt.emtools.fieldqc import find_duplicate_stations
    >>> f = np.array([10.0, 1.0]); z = np.ones((2, 2, 2), complex)
    >>> mk = lambda s, zz: NS(station=s, Z=NS(freq=f, z=zz))
    >>> find_duplicate_stations([mk("A", z), mk("B", z),
    ...                          mk("C", 2 * z)])
    [['A', 'B']]
    """
    names, fs, zs, _ = _station_arrays(sites)
    groups: list[list[int]] = []
    used = set()
    for a in range(len(names)):
        if a in used:
            continue
        grp = [a]
        for b in range(a + 1, len(names)):
            if b in used or fs[a].shape != fs[b].shape \
                    or zs[a].shape != zs[b].shape:
                continue
            if np.allclose(fs[a], fs[b], rtol=rtol) and np.allclose(
                    zs[a], zs[b], rtol=rtol, atol=0, equal_nan=True):
                grp.append(b)
        if len(grp) > 1:
            used.update(grp)
            groups.append(grp)
    return [[names[i] for i in g] for g in groups]


# -------------------------------------------------------- point flags
def hampel_mask(
    x: np.ndarray,
    *,
    half: int = 2,
    nsig: float = 3.0,
    floor: float = 0.0,
) -> np.ndarray:
    """Outlier mask of a Hampel filter along a 1-D series.

    A finite point is an outlier when it departs from the median of
    its window (``half`` points each side, excluding itself, NaN
    ignored) by more than ``max(nsig * 1.4826 * MAD, floor)``.

    Returns
    -------
    ndarray of bool
        True for outliers; False for NaN or under-populated windows
        (fewer than two finite neighbours).
    """
    x = np.asarray(x, float)
    bad = np.zeros(x.shape, bool)
    for i in range(x.size):
        if not np.isfinite(x[i]):
            continue
        w = np.r_[x[max(0, i - half): i], x[i + 1: i + half + 1]]
        w = w[np.isfinite(w)]
        if w.size < 2:
            continue
        med = np.median(w)
        mad = 1.4826 * np.median(np.abs(w - med))
        bad[i] = abs(x[i] - med) > max(nsig * mad, floor)
    return bad


def point_qc_flags(
    z: np.ndarray,
    z_err: np.ndarray,
    freq: np.ndarray,
    *,
    err_max: float = 0.3,
    quadrant_tol: float = 5.0,
    hampel_half: int = 2,
    hampel_nsig: float = 3.0,
    floor_log_rho: float = 0.2,
    floor_phase: float = 10.0,
) -> pd.DataFrame:
    """Per-frequency flags of one station (polarity-normalised z).

    Parameters
    ----------
    z, z_err : ndarray, shape (n_freq, 2, 2)
        Impedance and its error, **after** polarity correction.
    freq : ndarray, shape (n_freq,)
        Frequencies (Hz).
    err_max : float, default 0.3
        Maximum relative error :math:`\\delta|Z|/|Z|`.
    quadrant_tol : float, default 5.0
        Tolerance (deg) outside :math:`[0, 90]`.
    hampel_half, hampel_nsig : int, float
        Hampel window half-width and threshold.
    floor_log_rho, floor_phase : float
        Minimum Hampel departure in decades and degrees.

    Returns
    -------
    pandas.DataFrame
        Columns ``freq``, ``log_rho_m`` (``0.2|Z|^2/f``, field units),
        ``phase_m`` (deg, yx shifted by +180 to [0, 90]),
        ``relerr_m`` and boolean ``err_m``, ``quadrant_m``,
        ``spike_m`` for ``m`` in ``xy, yx``.
    """
    f = np.asarray(freq, float)
    df = pd.DataFrame({"freq": f})
    for m, (i, j) in _MODES.items():
        zz = np.asarray(z)[:, i, j]
        with np.errstate(all="ignore"):
            lr = np.log10(0.2 / f * np.abs(zz) ** 2)
            rel = np.asarray(z_err)[:, i, j] / np.abs(zz)
        ph = np.degrees(np.angle(zz))
        if m == "yx":
            ph = _wrap180(ph + 180.0)
        df[f"log_rho_{m}"] = lr
        df[f"phase_{m}"] = ph
        df[f"relerr_{m}"] = rel
        df[f"err_{m}"] = ~(rel <= err_max)
        df[f"quadrant_{m}"] = (ph < -quadrant_tol) | (
            ph > 90.0 + quadrant_tol)
        df[f"spike_{m}"] = (
            hampel_mask(lr, half=hampel_half, nsig=hampel_nsig,
                        floor=floor_log_rho)
            | hampel_mask(ph, half=hampel_half, nsig=hampel_nsig,
                          floor=floor_phase))
    return df


def point_qc_table(
    sites: Any,
    *,
    band: tuple[float, float] = (3.0, 30.0),
    assume: str = "E",
    neighbours: int = 2,
    neighbour_tol: float = 20.0,
    **flag_kws,
) -> pd.DataFrame:
    """Polarity-normalised point flags for a whole profile.

    Stations are taken in the given order (profile order). Each is
    polarity-corrected (:func:`detect_polarity`,
    :func:`correct_polarity`), flagged with :func:`point_qc_flags`,
    then compared with up to *neighbours* stations on each side at
    the same frequency (frequencies matched within 0.02 %).

    Parameters
    ----------
    sites : Sites or sequence of EDIFile
        Stations in profile order.
    band : tuple of float, default (3.0, 30.0)
        Reference band for :func:`detect_polarity`.
    assume : {"E", "H"}, default "E"
        Passed to :func:`correct_polarity`.
    neighbours : int, default 2
        Stations on each side for the neighbour test.
    neighbour_tol : float, default 20.0
        Maximum phase departure from the neighbour median (deg).
    **flag_kws
        Passed to :func:`point_qc_flags`.

    Returns
    -------
    pandas.DataFrame
        Long table, one row per station and frequency, with the
        point flags, ``neighbour_m`` and an overall ``bad_m``
        (any flag).
    """
    names, fs, zs, es = _station_arrays(sites)
    frames = []
    for k, (n, f, z, e) in enumerate(zip(names, fs, zs, es)):
        pol = detect_polarity(z, f, band=band)
        zc = correct_polarity(z, flip_xy=pol["flip_xy"],
                              flip_yx=pol["flip_yx"], assume=assume)
        d = point_qc_flags(zc, e, f, **flag_kws)
        d.insert(0, "station", n)
        d.insert(1, "order", k)
        frames.append(d)
    tab = pd.concat(frames, ignore_index=True)
    # frequencies are matched across stations on log10(f) rounded to
    # 4 decimals (0.02 % relative tolerance)
    tab["fkey"] = np.round(np.log10(tab.freq), 4)
    for m in _MODES:
        piv = tab.pivot_table(index="order", columns="fkey",
                              values=f"phase_{m}", aggfunc="first")
        med = np.full(piv.shape, np.nan)
        arr = piv.to_numpy()
        for r in range(arr.shape[0]):
            nb = [q for q in range(max(0, r - neighbours),
                                   min(arr.shape[0],
                                       r + neighbours + 1)) if q != r]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                med[r] = np.nanmedian(arr[nb], axis=0) if nb else np.nan
        lut = pd.DataFrame(med, index=piv.index, columns=piv.columns)
        nbv = np.array([lut.at[o, fk] for o, fk in
                        zip(tab.order, tab.fkey)])
        tab[f"neighbour_{m}"] = np.abs(tab[f"phase_{m}"] - nbv) \
            > neighbour_tol
        tab[f"bad_{m}"] = tab[[f"err_{m}", f"quadrant_{m}",
                               f"spike_{m}", f"neighbour_{m}"]].any(
            axis=1)
    return tab.drop(columns="fkey")


def detect_band_break(
    table: pd.DataFrame,
    *,
    mode: str = "xy",
    min_jump: float = 0.3,
) -> pd.DataFrame:
    """Frequencies where the survey-wide rejection rate jumps.

    A systematic problem confined to part of the spectrum (a bad
    calibration or band join, a strong coherent noise band) shows as
    a step in the fraction of rejected points between two adjacent
    frequencies across *all* stations.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of :func:`point_qc_table`.
    mode : {"xy", "yx"}, default "xy"
    min_jump : float, default 0.3
        Minimum change of rejected fraction between neighbours.

    Returns
    -------
    pandas.DataFrame
        ``f_low``, ``f_high``, ``bad_low``, ``bad_high`` and
        ``jump`` (high minus low) for each step, largest first.
    """
    rate = table.groupby("freq")[f"bad_{mode}"].mean().sort_index()
    f, r = rate.index.to_numpy(), rate.to_numpy(float)
    rows = [{"f_low": f[k], "f_high": f[k + 1], "bad_low": r[k],
             "bad_high": r[k + 1], "jump": r[k + 1] - r[k]}
            for k in range(f.size - 1)
            if abs(r[k + 1] - r[k]) >= min_jump]
    out = pd.DataFrame(rows, columns=["f_low", "f_high", "bad_low",
                                      "bad_high", "jump"])
    return out.reindex(out.jump.abs().sort_values(
        ascending=False).index).reset_index(drop=True)


# ------------------------------------------------------- static shift
def estimate_static_shift_neighbours(
    table: pd.DataFrame,
    *,
    band: tuple[float, float] = (10.0, 30.0),
    neighbours: int = 3,
    min_points: int = 3,
    min_abs: float = 0.1,
    max_abs: float = 1.0,
    stations: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Static shift of each mode relative to neighbouring stations.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of :func:`point_qc_table` (profile order).
    band : tuple of float, default (10.0, 30.0)
        Frequencies (Hz) used for the estimate.
    neighbours : int, default 3
        Stations taken on each side as reference.
    min_points : int, default 3
        Minimum QC-good frequencies for an estimate to apply.
    min_abs, max_abs : float
        An estimate is applicable when
        ``min_abs < |shift| <= max_abs`` decades.
    stations : sequence of str, optional
        Restrict reference *and* target stations (e.g. to good
        grades); others get NaN.

    Returns
    -------
    pandas.DataFrame
        Per station: ``shift_m`` (decades, rho_obs / rho_ref),
        ``n_m`` and ``apply_m``. Divide :math:`Z` rows by
        :math:`\\sqrt{10^{shift}}` to correct (row x for xy, row y
        for yx).
    """
    sel = table[(table.freq >= band[0]) & (table.freq <= band[1])]
    order = (table[["station", "order"]].drop_duplicates()
             .sort_values("order").station.tolist())
    if stations is not None:
        keep = set(stations)
        order = [s for s in order if s in keep]
    out = []
    for idx, st in enumerate(order):
        nb = order[max(0, idx - neighbours): idx] + \
            order[idx + 1: idx + 1 + neighbours]
        row = {"station": st}
        for m in _MODES:
            good = sel[~sel[f"bad_{m}"]]
            me = good[good.station == st].set_index("freq")
            ref = (good[good.station.isin(nb)]
                   .groupby("freq")[f"log_rho_{m}"]
                   .agg(["median", "count"]))
            ref = ref[ref["count"] >= 2]
            common = me.index.intersection(ref.index)
            d = (me.loc[common, f"log_rho_{m}"]
                 - ref.loc[common, "median"]).to_numpy(float)
            s = float(np.median(d)) if d.size else np.nan
            row[f"shift_{m}"] = s
            row[f"n_{m}"] = int(d.size)
            row[f"apply_{m}"] = bool(
                d.size >= min_points and np.isfinite(s)
                and min_abs < abs(s) <= max_abs)
        out.append(row)
    return pd.DataFrame(out)
