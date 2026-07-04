from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..api.station import PYCSAMT_STATION_RENDERING
from ..api.labels import LOG10_PERIOD_LABEL
from ._core import (
    ensure_sites,
    _apply_each,
    _axes_list,
    _iter_items,
    _get_z_block,
    _get_t_block,
    _name,
    _station_positions,
)


@dataclass
class EMAPFilterResult:
    """Container returned by confidence-gated EMAP filtering."""

    sites: Any
    report: pd.DataFrame
    decisions: pd.DataFrame
    method: str
    confidence_method: str
    ci_hi: float
    ci_lo: float

    @property
    def n_preserved(self) -> int:
        """Number of station-frequency rows left unchanged."""
        if self.decisions.empty:
            return 0
        return int((self.decisions["action"] == "preserved").sum())

    @property
    def n_blended(self) -> int:
        """Number of station-frequency rows partially blended."""
        if self.decisions.empty:
            return 0
        return int((self.decisions["action"] == "blended").sum())

    @property
    def n_filtered(self) -> int:
        """Number of station-frequency rows fully filtered."""
        if self.decisions.empty:
            return 0
        return int((self.decisions["action"] == "filtered").sum())

    def summary(self) -> str:
        """Return a compact text summary."""
        return (
            "EMAPFilterResult("
            f"method={self.method!r}, confidence={self.confidence_method!r}, "
            f"preserved={self.n_preserved}, blended={self.n_blended}, "
            f"filtered={self.n_filtered})"
        )

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()

# ----------------------------- SNR table -------------------------------- #

def snr_table(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    rows = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        out = _get_z_block(ed, with_errors=True)
        if len(out) == 4:
            Z, z, fr, ze = out
        else:
            Z, z, fr, ze = (out + (None,))[:4]
        if Z is None:
            continue
        a = np.sqrt(np.nanmean(np.abs(z) ** 2, axis=(1, 2)))
        if isinstance(ze, np.ndarray) and ze.shape == z.shape:
            e = np.sqrt(np.nanmean(np.abs(ze) ** 2, axis=(1, 2)))
        else:
            e = np.full_like(a, np.nan)
        snr = a / (e + 1e-12)
        for f, s in zip(fr, snr):
            rows.append(dict(station=st, freq=float(f), snr=float(s)))
    return pd.DataFrame.from_records(rows)


# --------------------------- power-line notching ------------------------- #

def _harm_mask(fr: np.ndarray, mains: float, n_harm: int,
               tol_hz: float) -> np.ndarray:
    kk = np.arange(1, int(n_harm) + 1, dtype=float)
    fH = kk * float(mains)
    m = np.zeros(fr.size, dtype=bool)
    for fh in fH:
        m |= np.abs(fr - fh) <= float(tol_hz)
    return m


def _interp_rows(
    y: np.ndarray,
    good: np.ndarray,
) -> np.ndarray:
    # y: (n, ...) complex; interp along axis 0
    if y.ndim == 1:
        y = y[:, None]
    x = np.arange(y.shape[0])
    yi = y.copy()
    g = good & np.all(np.isfinite(y), axis=tuple(range(1, y.ndim)))
    if g.sum() < 2:
        yi[~g] = np.nan
        return yi.squeeze()
    xi = x[g]
    for idx in np.argwhere(~g).ravel():
        # nearest two good neighbors (linear)
        j = np.searchsorted(xi, idx)
        j0 = max(0, min(j - 1, xi.size - 1))
        j1 = max(0, min(j, xi.size - 1))
        a, b = xi[j0], xi[j1]
        if a == b:
            yi[idx] = y[a]
        else:
            t = (idx - a) / (b - a)
            yi[idx] = (1 - t) * y[a] + t * y[b]
    return yi.squeeze()


def notch_powerline(
    sites: Any,
    *,
    mains_hz: float = 50.0,   # 50 or 60
    n_harm: int = 30,
    tol_hz: float = 0.08,     # Hz window around each harmonic
    mode: str = "interp",     # mask|interp
    also: str = "both",       # z|tipper|both
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        mH = _harm_mask(fr, mains_hz, n_harm, tol_hz)
        z2 = z.copy()
        if mode == "mask":
            z2[mH] = np.nan
        else:
            z2[mH] = np.nan  # mark, then fill
            good = ~np.isnan(z2).any(axis=(1, 2))
            z2 = _interp_rows(z2, good)
        Z.z = z2
        if also in ("tipper", "both"):
            T, t, ft = _get_t_block(ed)
            if T is not None and t is not None:
                mt = _harm_mask(ft, mains_hz, n_harm, tol_hz)
                t2 = t.copy()
                if mode == "mask":
                    t2[mt] = np.nan
                else:
                    t2[mt] = np.nan
                    good = ~np.isnan(t2).any(axis=1)
                    t2 = _interp_rows(t2, good)
                T.tipper = t2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# ------------------------- log-frequency smoothing ----------------------- #

def _smooth1d(y: np.ndarray, win: int, kind: str) -> np.ndarray:
    if win <= 1:
        return y
    w = np.ones(int(win), dtype=float)
    if kind == "tri":
        w = np.convolve(w, w, mode="full")
    w = w / np.sum(w)
    if y.ndim == 1:
        yy = np.convolve(y, w, mode="same")
        return yy
    # apply per column
    out = np.empty_like(y)
    for j in range(y.shape[1]):
        out[:, j] = np.convolve(y[:, j], w, mode="same")
    return out


def smooth_logfreq(
    sites: Any,
    *,
    win: int = 5,
    kind: str = "tri",         # box|tri
    also: str = "both",        # z|tipper|both
    gate_snr: Optional[float] = None,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    ST = snr_table(S) if gate_snr is not None else None

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        y = z.reshape(z.shape[0], -1)
        ok = np.isfinite(y).all(axis=1)
        if gate_snr is not None and not ST.empty:
            st = _name(ed, 0)
            sdf = ST[ST["station"] == st]
            if not sdf.empty:
                idx = np.searchsorted(sdf["freq"].to_numpy(), fr)
                idx = np.clip(idx, 0, len(sdf) - 1)
                ok &= (sdf["snr"].to_numpy()[idx] >= gate_snr)
        y2 = y.copy()
        y2[ok] = _smooth1d(y[ok], win, kind)
        Z.z = y2.reshape(z.shape)
        if also in ("tipper", "both"):
            T, t, ft = _get_t_block(ed)
            if T is not None and t is not None:
                yt = t
                okt = np.isfinite(yt).all(axis=1)
                T.tipper = _smooth1d(yt, win, kind) if okt.any() else yt
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# ----------------------- rho/phase trend smoothing ----------------------- #

_COMPONENT_INDEX: Dict[str, Tuple[int, int]] = {
    "xx": (0, 0),
    "xy": (0, 1),
    "yx": (1, 0),
    "yy": (1, 1),
}


def _resolve_components(components: str | Sequence[str]) -> Tuple[str, ...]:
    if isinstance(components, str):
        key = components.strip().lower().replace("-", "_")
        aliases = {
            "offdiag": ("xy", "yx"),
            "off_diagonal": ("xy", "yx"),
            "offdiagonal": ("xy", "yx"),
            "diagonal": ("xx", "yy"),
            "diag": ("xx", "yy"),
            "all": ("xx", "xy", "yx", "yy"),
            "both": ("xy", "yx"),
        }
        if key in aliases:
            return aliases[key]
        if key in _COMPONENT_INDEX:
            return (key,)
        parts = [p.strip().lower() for p in key.replace(",", " ").split()]
    else:
        parts = [str(p).strip().lower() for p in components]
    out: List[str] = []
    for part in parts:
        if part not in _COMPONENT_INDEX:
            raise ValueError(
                "components must be one of 'offdiag', 'diagonal', 'all', "
                "'xx', 'xy', 'yx', 'yy', or a sequence of those values; "
                f"got {part!r}."
            )
        if part not in out:
            out.append(part)
    if not out:
        raise ValueError("components cannot be empty.")
    return tuple(out)


def _robust_polyfit_predict(
    x: np.ndarray,
    y: np.ndarray,
    *,
    degree: int,
    min_points: int,
    robust: bool,
    robust_iters: int,
) -> np.ndarray:
    out = np.asarray(y, dtype=float).copy()
    good = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(good) < max(2, int(min_points)):
        return out
    xv = np.asarray(x[good], dtype=float)
    yv = np.asarray(y[good], dtype=float)
    deg = int(max(0, min(int(degree), xv.size - 1)))
    if deg == 0:
        out[good] = np.nanmedian(yv)
        return out

    weights = np.ones_like(yv, dtype=float)
    coeff = None
    n_iter = max(1, int(robust_iters) if robust else 1)
    for _ in range(n_iter):
        try:
            coeff = np.polyfit(xv, yv, deg=deg, w=weights)
        except (np.linalg.LinAlgError, ValueError, TypeError):
            return out
        pred = np.polyval(coeff, xv)
        if not robust:
            break
        resid = yv - pred
        scale = 1.4826 * np.nanmedian(np.abs(resid - np.nanmedian(resid)))
        if not np.isfinite(scale) or scale <= 0:
            break
        u = resid / (4.685 * scale)
        weights = (1.0 - u**2) ** 2
        weights[np.abs(u) >= 1.0] = 0.0
        if np.count_nonzero(weights > 0) <= deg:
            break
    if coeff is not None:
        out[good] = np.polyval(coeff, xv)
    return out


def _smooth_phase_deg(
    x: np.ndarray,
    phase_deg: np.ndarray,
    *,
    degree: int,
    min_points: int,
    robust: bool,
    robust_iters: int,
) -> np.ndarray:
    phase_deg = np.asarray(phase_deg, dtype=float)
    out = phase_deg.copy()
    good = np.isfinite(x) & np.isfinite(phase_deg)
    if np.count_nonzero(good) < max(2, int(min_points)):
        return out
    idx = np.flatnonzero(good)
    order = idx[np.argsort(x[idx])]
    unwrapped = np.unwrap(np.deg2rad(phase_deg[order]))
    smoothed = _robust_polyfit_predict(
        x[order],
        unwrapped,
        degree=degree,
        min_points=min_points,
        robust=robust,
        robust_iters=robust_iters,
    )
    out[order] = np.rad2deg(smoothed)
    return out


def smooth_rho_phase(
    sites: Any,
    *,
    components: str | Sequence[str] = "offdiag",
    degree: int = 3,
    min_points: Optional[int] = None,
    smooth_rho: bool = True,
    smooth_phase: bool = True,
    robust: bool = True,
    robust_iters: int = 3,
    blend: float = 1.0,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    r"""
    Smooth apparent resistivity and phase trends, then rebuild ``Z``.

    The function operates station-by-station along the frequency axis.  It
    fits polynomial trends versus :math:`\log_{10}(f)` to
    :math:`\log_{10}(\rho_a)` and to unwrapped impedance phase, then writes
    the corresponding complex impedance back into each selected tensor
    component.  This keeps apparent resistivity and phase physically coupled
    through the same complex ``Z`` tensor instead of only smoothing display
    arrays.

    Parameters
    ----------
    sites : object
        Any input accepted by :func:`ensure_sites`.
    components : {"offdiag", "diagonal", "all", "xx", "xy", "yx", "yy"} \
or sequence, default "offdiag"
        Tensor components to smooth.  The default targets ``xy`` and ``yx``
        because they are the usual MT/CSAMT apparent-resistivity and phase
        components used for interpretation and 2-D preparation.
    degree : int, default 3
        Polynomial degree for the log-frequency trend.  It is automatically
        reduced when a station has too few valid frequencies.
    min_points : int or None, default None
        Minimum number of finite points required per component.  If ``None``,
        uses ``degree + 2``.
    smooth_rho, smooth_phase : bool, default True
        Select whether the impedance amplitude, phase angle, or both are
        replaced by the fitted trend.
    robust : bool, default True
        Use a Tukey-style iteratively reweighted polynomial fit to reduce the
        influence of isolated spikes.
    robust_iters : int, default 3
        Maximum robust reweighting iterations.
    blend : float, default 1.0
        Blend between original and smoothed curves.  ``1`` fully applies the
        trend; ``0.5`` applies half of the correction.
    inplace : bool, default False
        If ``True``, mutate the normalized input sites.  Otherwise work on a
        best-effort copy of the underlying EDI objects and return new
        ``Sites``.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`ensure_sites` / ``to_edis``.

    Returns
    -------
    pycsamt.site.base.Sites
        Sites containing the smoothed impedance tensors.

    Notes
    -----
    Apparent resistivity is smoothed in logarithmic space because
    :math:`\rho_a` commonly spans orders of magnitude.  Phase is unwrapped
    before fitting so that crossings near :math:`\pm 180^\circ` do not create
    artificial jumps.
    """
    if not smooth_rho and not smooth_phase:
        return ensure_sites(
            sites, recursive=recursive, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )

    comps = _resolve_components(components)
    blend = float(np.clip(blend, 0.0, 1.0))
    degree = int(max(0, degree))
    min_points = int(min_points if min_points is not None else degree + 2)

    if inplace:
        S = ensure_sites(
            sites, recursive=recursive, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )
    else:
        from ..site.base import to_edis

        edis = to_edis(
            sites,
            copy=True,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
        S = ensure_sites(
            edis, recursive=False, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )

    for ed in _iter_items(S):
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        z2 = np.asarray(z, dtype=np.complex128).copy()
        fr = np.asarray(fr, dtype=float).ravel()
        n = min(z2.shape[0], fr.size)
        if n == 0:
            continue
        z2 = z2[:n].copy()
        fr = fr[:n]
        valid_freq = np.isfinite(fr) & (fr > 0.0)
        if np.count_nonzero(valid_freq) < max(2, min_points):
            continue
        x = np.full(fr.shape, np.nan, dtype=float)
        x[valid_freq] = np.log10(fr[valid_freq])

        for comp in comps:
            i, j = _COMPONENT_INDEX[comp]
            zij = z2[:, i, j]
            amp = np.abs(zij)
            phase = np.rad2deg(np.angle(zij))

            amp_new = amp.copy()
            if smooth_rho:
                rho = np.full_like(fr, np.nan, dtype=float)
                ok_rho = valid_freq & np.isfinite(amp) & (amp > 0.0)
                rho[ok_rho] = (amp[ok_rho] ** 2) / (5.0 * fr[ok_rho])
                log_rho = np.full_like(rho, np.nan, dtype=float)
                ok_log = ok_rho & (rho > 0.0)
                log_rho[ok_log] = np.log10(rho[ok_log])
                log_rho_fit = _robust_polyfit_predict(
                    x,
                    log_rho,
                    degree=degree,
                    min_points=min_points,
                    robust=robust,
                    robust_iters=robust_iters,
                )
                if blend < 1.0:
                    log_rho_fit = (
                        (1.0 - blend) * log_rho + blend * log_rho_fit
                    )
                rho_fit = 10.0 ** log_rho_fit
                ok_fit = ok_rho & np.isfinite(rho_fit) & (rho_fit > 0.0)
                amp_new[ok_fit] = np.sqrt(5.0 * fr[ok_fit] * rho_fit[ok_fit])

            phase_new = phase.copy()
            if smooth_phase:
                phase_fit = _smooth_phase_deg(
                    x,
                    phase,
                    degree=degree,
                    min_points=min_points,
                    robust=robust,
                    robust_iters=robust_iters,
                )
                if blend < 1.0:
                    phase_fit = (1.0 - blend) * phase + blend * phase_fit
                ok_phi = valid_freq & np.isfinite(phase_fit)
                phase_new[ok_phi] = phase_fit[ok_phi]

            ok = valid_freq & np.isfinite(amp_new) & np.isfinite(phase_new)
            z2[ok, i, j] = amp_new[ok] * np.exp(1j * np.deg2rad(phase_new[ok]))

        try:
            z_full = np.asarray(getattr(Z, "z"), dtype=np.complex128).copy()
            z_full[:n] = z2
            Z.z = z_full
        except Exception:
            Z.z = z2

    return S


# --------------------- group-trend shrinkage (robust) -------------------- #

def _build_group_trend(
    sites: Any,
    groups: Dict[str, List[str]],
) -> Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    # trend[grp][station_ref] = (fr_union, z_med_on_union)
    trend: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]] = {}
    # union frequency per group
    for g, sts in groups.items():
        G: List[float] = []
        pool = []
        for i, ed in enumerate(_iter_items(sites)):
            st = _name(ed, i)
            if st not in sts:
                continue
            Z, z, fr = _get_z_block(ed)
            if Z is None:
                continue
            G.append(fr)
            pool.append((st, fr, z))
        if not pool:
            continue
        Gu = np.unique(np.concatenate(G))
        Zm = []
        for k in range(Gu.size):
            vals = []
            for st, fr, z in pool:
                j = np.searchsorted(fr, Gu[k])
                j = np.clip(j, 0, fr.size - 1)
                vals.append(z[j])
            vals = np.asarray(vals)
            Zm.append(np.nanmedian(vals, axis=0))
        trend[g] = {"_union": (Gu, np.asarray(Zm))}
    return trend


def shrink_to_group_trend(
    sites: Any,
    *,
    groups: Optional[Dict[str, List[str]]] = None,
    group_key: Optional[str] = None,
    lam: float = 0.25,          # 0..1; 0=no change, 1=trend
    gate_harm: bool = True,
    mains_hz: float = 50.0,
    n_harm: int = 30,
    tol_hz: float = 0.08,
    also: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if groups is None:
        # auto: single group with all stations
        names = [_name(ed, i) for i, ed in enumerate(_iter_items(S))]
        groups = {"ALL": names}
    Tref = _build_group_trend(S, groups)

    def _one(Si):
        ed = next(_iter_items(Si))
        st = _name(ed, 0)
        # find group containing st
        g = None
        for k, v in groups.items():
            if st in v:
                g = k; break
        if g is None or g not in Tref:
            return Si
        (Gu, Zu) = Tref[g]["_union"]
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        idx = np.searchsorted(Gu, fr)
        idx = np.clip(idx, 0, Gu.size - 1)
        tr = Zu[idx]
        if gate_harm:
            mH = _harm_mask(fr, mains_hz, n_harm, tol_hz)
        else:
            mH = np.ones(fr.size, dtype=bool)
        z2 = z.copy()
        z2[mH] = (1.0 - lam) * z[mH] + lam * tr[mH]
        Z.z = z2
        if also in ("tipper", "both"):
            T, t, ft = _get_t_block(ed)
            if T is not None and t is not None:
                # simple shrink on tipper to its group median
                tt = t.copy()
                T.tipper = (1.0 - lam) * t + lam * tt
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# ------------------------------ pipeline --------------------------------- #

def remove_noise_pipeline(
    sites: Any,
    *,
    mains_hz: float = 50.0,
    n_harm: int = 30,
    tol_hz: float = 0.08,
    notch_mode: str = "interp",
    smooth_win: int = 5,
    smooth_kind: str = "tri",
    gate_snr: Optional[float] = 2.5,
    group_shrink: bool = False,
    shrink_lam: float = 0.25,
    groups: Optional[Dict[str, List[str]]] = None,
    also: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = notch_powerline(
        sites,
        mains_hz=mains_hz,
        n_harm=n_harm,
        tol_hz=tol_hz,
        mode=notch_mode,
        also=also,
        inplace=False,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    S = smooth_logfreq(
        S,
        win=smooth_win,
        kind=smooth_kind,
        also=also,
        gate_snr=gate_snr,
        inplace=False,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    if group_shrink:
        S = shrink_to_group_trend(
            S,
            groups=groups,
            lam=shrink_lam,
            mains_hz=mains_hz,
            n_harm=n_harm,
            tol_hz=tol_hz,
            also=also,
            inplace=False,
            recursive=False,
            on_dup=on_dup,
            strict=False,
            verbose=verbose,
        )
    if inplace:
        # copy back into input sites
        _ = _apply_each(
            sites, lambda Si: Si, inplace=True, verbose=verbose
        )
        return sites
    return S

# ------------------------- 1) Hampel outlier filter --------------------- #

def _hampel_1d(y: np.ndarray, k: int, nsig: float) -> np.ndarray:
    if y.size == 0 or k <= 0:
        return y
    y2 = y.copy()
    n = y.size
    for i in range(n):
        lo = max(0, i - k); hi = min(n, i + k + 1)
        win = y[lo:hi]
        med = np.nanmedian(win)
        mad = np.nanmedian(np.abs(win - med)) + 1e-12
        if np.abs(y[i] - med) > nsig * 1.4826 * mad:
            y2[i] = med
    return y2


def hampel_filter_freq(
    sites: Any,
    *,
    win: int = 3,
    nsig: float = 3.0,
    on: str = "both",          # z|tipper|both
    domain: str = "reim",      # reim|magphase
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is not None:
            Y = z.reshape(z.shape[0], -1)
            if domain == "magphase":
                m = np.abs(Y); p = np.angle(Y)
                m2 = np.vstack([
                    _hampel_1d(m[:, j], win, nsig)
                    for j in range(m.shape[1])
                ]).T
                Y2 = m2 * np.exp(1j * p)
            else:
                R = np.real(Y); I = np.imag(Y)
                R2 = np.vstack([
                    _hampel_1d(R[:, j], win, nsig)
                    for j in range(R.shape[1])
                ]).T
                I2 = np.vstack([
                    _hampel_1d(I[:, j], win, nsig)
                    for j in range(I.shape[1])
                ]).T
                Y2 = R2 + 1j * I2
            Z.z = Y2.reshape(z.shape)
        if on in ("tipper", "both"):
            from ._core import _get_t_block
            T, t, ft = _get_t_block(ed)
            if T is not None:
                r = np.real(t); im = np.imag(t)
                r2 = np.vstack([
                    _hampel_1d(r[:, j], win, nsig)
                    for j in range(r.shape[1])
                ]).T
                i2 = np.vstack([
                    _hampel_1d(im[:, j], win, nsig)
                    for j in range(im.shape[1])
                ]).T
                T.tipper = r2 + 1j * i2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# -------- 2) Spatial median smooth (across stations at fixed freq) ------ #

def spatial_median_filter(
    sites: Any,
    *,
    half_window: int = 2,
    lam: float = 0.25,          # 0..1 shrink toward local median
    on: str = "z",              # z|tipper|both
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    items = list(_iter_items(S))
    n = len(items)

    def _nbrs(i: int) -> List[int]:
        lo = max(0, i - half_window)
        hi = min(n - 1, i + half_window)
        ids = list(range(lo, hi + 1))
        if i in ids:
            ids.remove(i)
        return ids

    # prefetch Z-blocks
    ZB = []
    for ed in items:
        ZB.append(_get_z_block(ed))

    def _one(Si):
        ed = next(_iter_items(Si))
        i = items.index(ed)
        Z, z, fr = ZB[i]
        if Z is not None:
            z2 = z.copy()
            for k in range(z.shape[0]):
                pool = []
                for j in _nbrs(i):
                    Zj, zj, frj = ZB[j]
                    if Zj is None:
                        continue
                    jj = np.clip(
                        np.searchsorted(frj, fr[k]), 0, frj.size - 1
                    )
                    pool.append(zj[jj])
                if pool:
                    med = np.nanmedian(np.asarray(pool), axis=0)
                    z2[k] = (1.0 - lam) * z[k] + lam * med
            Z.z = z2
        if on in ("tipper", "both"):
            T, t, ft = _get_t_block(ed)
            if T is not None:
                # simple neighbor median on tipper too
                t2 = t.copy()
                for k in range(t.shape[0]):
                    pool = []
                    for j in _nbrs(i):
                        Tj, tj, ftj = None, None, None
                        try:
                            Tj, tj, ftj = _get_t_block(items[j])
                        except Exception:
                            pass
                        if Tj is None:
                            continue
                        jj = np.clip(
                            np.searchsorted(ftj, ft[k]),
                            0, ftj.size - 1
                        )
                        pool.append(tj[jj])
                    if pool:
                        med = np.nanmedian(np.asarray(pool), axis=0)
                        t2[k] = (1.0 - lam) * t[k] + lam * med
                T.tipper = t2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --------- 3) Low-rank denoise (profile RPCA-ish on |Z_xy|,|Z_yx|) ------ #

def _svd_rank_k(M: np.ndarray, r: int) -> np.ndarray:
    try:
        U, s, Vh = np.linalg.svd(M, full_matrices=False)
    except np.linalg.LinAlgError:
        U, s, Vh = np.linalg.svd(
            M + 1e-12 * np.random.randn(*M.shape),
            full_matrices=False,
        )
    r = int(max(1, min(r, min(M.shape))))
    return (U[:, :r] * s[:r]) @ Vh[:r, :]


def rpca_offdiag_denoise(
    sites: Any,
    *,
    rank: int = 2,
    keep_phase: bool = True,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    def _sta_name(ed) -> str:
        """Stable station key that survives the copies made by
        ``_apply_each(inplace=False)`` (object identity does not)."""
        for attr in ("station", "name", "sta"):
            v = getattr(ed, attr, None)
            if v:
                return str(v)
        site = getattr(ed, "site", None)
        if site is not None:
            v = (
                getattr(site, "name", None)
                or getattr(site, "station", None)
            )
            if v:
                return str(v)
        return repr(ed)

    items = list(_iter_items(S))
    # build matrix M(station, freq) from log |Z_xy| median(|Z_xy|,|Z_yx|)
    ST, FR, MAG, NAMES = [], [], [], []
    for i, ed in enumerate(items):
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        m = np.nanmedian(
            np.stack([np.abs(z[:, 0, 1]), np.abs(z[:, 1, 0])], axis=1),
            axis=1,
        )
        ST.append(i); FR.append(fr); MAG.append(np.log10(m + 1e-24))
        NAMES.append(_sta_name(ed))
    if not ST:
        return S
    row_of = {name: k for k, name in enumerate(NAMES)}
    G = np.unique(np.concatenate(FR))
    M = np.full((len(ST), G.size), np.nan, dtype=float)
    for row, (fr, lg) in enumerate(zip(FR, MAG)):
        idx = np.searchsorted(G, fr)
        idx = np.clip(idx, 0, G.size - 1)
        M[row, idx] = lg
    # simple fill of gaps by nearest valid along freq
    for i in range(M.shape[0]):
        r = M[i]
        g = np.isfinite(r)
        if g.sum() >= 2:
            xi = np.where(g)[0]
            for j in np.where(~g)[0]:
                k = np.searchsorted(xi, j)
                k0 = max(0, min(k - 1, xi.size - 1))
                k1 = max(0, min(k, xi.size - 1))
                a, b = xi[k0], xi[k1]
                t = 0.0 if a == b else (j - a) / (b - a)
                r[j] = (1 - t) * r[a] + t * r[b]
            M[i] = r
    # Rows with fewer than two finite samples keep NaNs after the gap
    # fill above — patch them with the global median so the SVD is
    # well-posed (an all-NaN matrix means nothing usable: return as-is).
    if not np.isfinite(M).any():
        return S
    M = np.where(np.isfinite(M), M, np.nanmedian(M))
    L = _svd_rank_k(M, rank)
    # push back magnitudes; keep phase if requested
    def _one(Si):
        ed = next(_iter_items(Si))
        i = row_of.get(_sta_name(ed))
        if i is None:
            return Si
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        idx = np.searchsorted(G, fr)
        idx = np.clip(idx, 0, G.size - 1)
        mag_clean = 10 ** L[i, idx]
        # per off-diagonal component scale
        for (a, b) in [(0, 1), (1, 0)]:
            comp = z[:, a, b]
            if keep_phase:
                ph = np.angle(comp)
                Z.z[:, a, b] = mag_clean * np.exp(1j * ph)
            else:
                # scale original by ratio on magnitude
                r = mag_clean / (np.abs(comp) + 1e-24)
                Z.z[:, a, b] = comp * r
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# -------- 4) Off-diagonal consistency (sym/anti-sym blend) -------------- #

def enforce_offdiag_consistency(
    sites: Any,
    *,
    mode: str = "anti",        # anti|sym
    lam: float = 0.5,          # blend toward constraint
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        xy = z[:, 0, 1]; yx = z[:, 1, 0]
        if mode == "sym":
            target_xy = 0.5 * (xy + yx)
            target_yx = target_xy
        else:
            target_xy = 0.5 * (xy - yx)
            target_yx = -target_xy
        z[:, 0, 1] = (1.0 - lam) * xy + lam * target_xy
        z[:, 1, 0] = (1.0 - lam) * yx + lam * target_yx
        Z.z = z
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# ------------- 5) Incoherent-frequency vote mask (SNR-based) ------------ #

def mask_incoherent_freqs(
    sites: Any,
    *,
    snr_thresh: float = 2.5,
    min_frac: float = 0.4,     # keep if ≥ this fraction pass
    also: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    T = snr_table(S)
    if T.empty:
        return S
    # union freq grid
    G = np.unique(T["freq"].to_numpy(dtype=float))
    # for each freq, fraction of stations with SNR ≥ thresh
    frac = np.zeros(G.size, dtype=float)
    for k, f in enumerate(G):
        m = T["freq"] == f
        sn = T.loc[m, "snr"].to_numpy(dtype=float)
        frac[k] = np.nanmean(sn >= snr_thresh)
    keepG = frac >= float(min_frac)

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is not None:
            idx = np.searchsorted(G, fr)
            idx = np.clip(idx, 0, G.size - 1)
            keep = keepG[idx]
            z2 = z.copy(); z2[~keep] = np.nan
            Z.z = z2
        if also in ("tipper", "both"):
            Tt, t, ft = _get_t_block(ed)
            if Tt is not None:
                idx = np.searchsorted(G, ft)
                idx = np.clip(idx, 0, G.size - 1)
                keep = keepG[idx]
                t2 = t.copy(); t2[~keep] = np.nan
                Tt.tipper = t2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def drop_freqs_manual(
    sites: Any,
    *,
    drop_freqs: Sequence[float] = (),
    tol_rel: float = 0.005,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    """Drop Z (and tipper) rows at user-specified frequencies.

    Parameters
    ----------
    sites : object
        Any input accepted by :func:`ensure_sites`.
    drop_freqs : sequence of float
        Frequencies (Hz) to remove.  Each value is matched within
        ``tol_rel`` relative tolerance: ``|f - f_drop| / f_drop < tol_rel``.
    tol_rel : float, default 0.005
        Relative frequency tolerance for matching (≈0.5%).
    inplace : bool, default False
        If ``True`` mutate; otherwise work on a copy.

    Returns
    -------
    Sites
        Sites with the specified frequency rows removed from Z, Z errors,
        apparent resistivity/phase, and tipper arrays where present.
    """
    drop_arr = np.asarray([float(f) for f in (drop_freqs or ())], dtype=float)
    drop_arr = drop_arr[np.isfinite(drop_arr) & (drop_arr > 0.0)]

    if drop_arr.size == 0:
        return ensure_sites(
            sites, recursive=recursive, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )

    if inplace:
        S = ensure_sites(
            sites, recursive=recursive, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )
    else:
        from ..site.base import to_edis
        edis = to_edis(
            sites, copy=True,
            recursive=recursive, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )
        S = ensure_sites(edis, recursive=False, on_dup=on_dup,
                         strict=strict, verbose=verbose)

    def _mask_for(fr: np.ndarray) -> np.ndarray:
        fr = np.asarray(fr, dtype=float)
        mask = np.zeros(fr.size, dtype=bool)
        for fd in drop_arr:
            mask |= np.abs(fr - fd) / (fd + 1e-30) < float(tol_rel)
        return mask

    def _assign(obj, public_name: str, private_name: str, value) -> None:
        assigned = False
        try:
            setattr(obj, public_name, value)
            assigned = True
        except Exception:
            pass
        if private_name and (not assigned or hasattr(obj, private_name)):
            try:
                setattr(obj, private_name, value)
            except Exception:
                pass

    def _slice_first_axis(arr, keep):
        if arr is None:
            return None
        try:
            a = np.asarray(arr)
        except Exception:
            return None
        if a.ndim < 1 or a.shape[0] != keep.size:
            return None
        return a[keep].copy()

    def _trim_z_block(Z, z, fr, keep):
        _assign(Z, "freq", "_freq", np.asarray(fr, dtype=float)[keep].copy())
        _assign(Z, "z", "_z", np.asarray(z, dtype=np.complex128)[keep].copy())
        for public_name, private_name in (
            ("z_err", "_z_err"),
            ("rho", "_rho"),
            ("resistivity", "_rho"),
            ("phase", "_phase"),
            ("rho_err", "_rho_err"),
            ("resistivity_err", "_rho_err"),
            ("phase_err", "_phase_err"),
        ):
            value = getattr(Z, public_name, None)
            sliced = _slice_first_axis(value, keep)
            if sliced is not None:
                _assign(Z, public_name, private_name, sliced)
        try:
            rot = getattr(Z, "rotation_angle", None)
            sliced = _slice_first_axis(rot, keep)
            if sliced is not None:
                Z.rotation_angle = sliced
        except Exception:
            pass
        try:
            Z.compute_resistivity_phase()
        except Exception:
            pass

    def _trim_t_block(Tt, t, ft, keep):
        _assign(Tt, "freq", "_freq", np.asarray(ft, dtype=float)[keep].copy())
        _assign(Tt, "tipper", "_tipper", np.asarray(t)[keep].copy())
        for public_name, private_name in (
            ("tipper_err", "_tipper_err"),
            ("amplitude", "_amplitude"),
            ("phase", "_phase"),
            ("amplitude_err", "_amplitude_err"),
            ("phase_err", "_phase_err"),
        ):
            value = getattr(Tt, public_name, None)
            sliced = _slice_first_axis(value, keep)
            if sliced is not None:
                _assign(Tt, public_name, private_name, sliced)
        try:
            Tt.compute_amp_phase()
        except Exception:
            pass

    for ed in _iter_items(S):
        Z, z, fr = _get_z_block(ed)
        if Z is not None:
            mask = _mask_for(fr)
            if mask.any():
                keep = ~mask
                if keep.any():
                    _trim_z_block(Z, z, fr, keep)
        Tt, t, ft = _get_t_block(ed)
        if Tt is not None:
            mask_t = _mask_for(ft)
            if mask_t.any():
                keep_t = ~mask_t
                if keep_t.any():
                    _trim_t_block(Tt, t, ft, keep_t)

    return S


# --------- 6) Static-shift removal — Torres-Verdín & Bostick (1992) ----- #

def _hanning_weights(dx: np.ndarray, w_H: float) -> np.ndarray:
    """Hanning (von Hann) spatial filter weights (Torres-Verdín & Bostick 1992).

    h(x) = (1 + cos(2π x / W_H)) / W_H  for |x| ≤ W_H/2, else 0.
    W_H is the full window width [m].
    """
    half = w_H / 2.0
    return np.where(
        np.abs(dx) <= half,
        np.maximum((1.0 + np.cos(2.0 * np.pi * dx / w_H)) / w_H, 0.0),
        0.0,
    )


def correct_static_shift(
    sites: Any,
    *,
    window_m: float = 1500.0,
    spacing_m: float = 200.0,
    comp: str = "det",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    """
    Remove static shift via Hanning adaptive moving-average (AMA) spatial filter.

    Implements the Torres-Verdín & Bostick (1992) approach used in
    Kouadio et al. (2024):

    1. For each frequency, build the spatial log(ρ_a) profile across stations.
    2. Apply a Hanning low-pass spatial filter with full-width ``window_m``.
    3. The static-shift correction factor at station *i* is
       ``C_i = sqrt(ρ_smooth_i / ρ_obs_i)`` (in log space:
       ``log C = 0.5 (log ρ_smooth − log ρ_obs)``).
    4. Update every Z component: ``Z_corrected = Z × C``.

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
        Input sites.
    window_m : float
        Full Hanning window width [m] (Torres-Verdín W_H).  Stations
        further than ``window_m/2`` contribute zero weight.
    spacing_m : float
        Fallback station spacing [m] used when EDI metadata carries no
        coordinate information.
    comp : {"det", "xy", "yx"}
        Apparent-resistivity component used to estimate the static shift.
        ``"det"`` uses the arithmetic mean of |Z_xy|² and |Z_yx|².
    inplace : bool
        If *True*, modify the input sites in place.  If *False* (default),
        return a new Sites object with corrected Z tensors.
    recursive, on_dup, strict, verbose
        Passed to :func:`ensure_sites`.

    Returns
    -------
    Sites
        Sites with static-shift-corrected Z tensors (when ``inplace=False``).
    """
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )

    items = list(_iter_items(S))
    if not items:
        return S

    names     = [_name(ed, i) for i, ed in enumerate(items)]
    positions = _station_positions(items, spacing_m)
    N         = len(items)

    # --- collect ρ_a on native frequency grids ---------------------------
    rho_data: List[Tuple[Optional[np.ndarray], Optional[np.ndarray]]] = []
    all_freqs: List[np.ndarray] = []

    for ed in items:
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            rho_data.append((None, None))
            continue
        fr_safe = np.maximum(fr, 1e-24)
        if comp == "xy":
            mag2 = np.abs(z[:, 0, 1]) ** 2
        elif comp == "yx":
            mag2 = np.abs(z[:, 1, 0]) ** 2
        else:   # det
            mag2 = 0.5 * (np.abs(z[:, 0, 1]) ** 2 + np.abs(z[:, 1, 0]) ** 2)
        # 0.2*|Z|^2/f, for Z in practical units (mV/km per nT) — matches
        # pycsamt.emtools.csumt's convention. The correction factor below
        # is a same-frequency ratio, so this was numerically safe even
        # before this fix (any constant scale factor cancels), but the
        # correct formula is used regardless for consistency.
        rho = 0.2 * mag2 / fr_safe
        lr  = np.where(rho > 0.0, np.log(rho), np.nan)
        rho_data.append((fr, lr))
        all_freqs.append(fr)

    if not all_freqs:
        return S

    G = np.unique(np.concatenate(all_freqs))
    F = G.size

    # --- map each station's log(ρ_a) onto the union frequency grid --------
    log_rho_mat = np.full((N, F), np.nan, dtype=float)
    for i, (fr_i, lr_i) in enumerate(rho_data):
        if fr_i is None:
            continue
        idx = np.searchsorted(G, fr_i)
        idx = np.clip(idx, 0, F - 1)
        log_rho_mat[i, idx] = lr_i

    # --- vectorised Hanning spatial smooth --------------------------------
    dx_mat = positions[:, None] - positions[None, :]   # [N, N]
    w_mat  = _hanning_weights(dx_mat, window_m)         # [N, N]

    valid = np.isfinite(log_rho_mat).astype(float)     # [N, F]
    # broadcast: w_mat[N,N,1] × valid[1,N,F] × log_rho[1,N,F]
    w3    = w_mat[:, :, None]
    lr3   = log_rho_mat[None, :, :]
    v3    = valid[None, :, :]
    num   = np.sum(w3 * v3 * lr3, axis=1)              # [N, F]
    denom = np.sum(w3 * v3,       axis=1)               # [N, F]
    log_rho_smooth = np.where(denom > 1e-30, num / denom, log_rho_mat)

    # --- log correction factor C = sqrt(ρ_smooth / ρ_obs) ----------------
    log_C = np.where(
        np.isfinite(log_rho_mat) & np.isfinite(log_rho_smooth),
        0.5 * (log_rho_smooth - log_rho_mat),
        0.0,
    )
    corr_map = {names[i]: (G, np.exp(log_C[i])) for i in range(N)}

    # --- apply per station -----------------------------------------------
    def _one(Si):
        ed = next(_iter_items(Si))
        st = _name(ed, 0)
        Z, z, fr = _get_z_block(ed)
        if Z is None or st not in corr_map:
            return Si
        G_c, C_arr = corr_map[st]
        idx = np.searchsorted(G_c, fr)
        idx = np.clip(idx, 0, G_c.size - 1)
        C   = C_arr[idx]
        Z.z = z * C[:, None, None]
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def _emap_component_indices(component: str) -> Tuple[Tuple[int, int], ...]:
    """Return tensor component indices selected by an EMAP filter."""
    component = str(component).lower()
    mapping = {
        "xx": ((0, 0),),
        "xy": ((0, 1),),
        "yx": ((1, 0),),
        "yy": ((1, 1),),
        "offdiag": ((0, 1), (1, 0)),
        "all": ((0, 0), (0, 1), (1, 0), (1, 1)),
    }
    if component not in mapping:
        msg = (
            "component must be one of 'all', 'offdiag', "
            "'xx', 'xy', 'yx', or 'yy'."
        )
        raise ValueError(msg)
    return mapping[component]


def _nearest_frequency_row(
    z: np.ndarray,
    fr: np.ndarray,
    freq: float,
    *,
    rtol: float,
) -> Optional[np.ndarray]:
    """Return the nearest tensor row when the frequency matches."""
    if fr.size == 0 or not np.isfinite(freq):
        return None
    idx = int(np.nanargmin(np.abs(fr - freq)))
    if not np.isclose(fr[idx], freq, rtol=rtol, atol=1e-12):
        return None
    return z[idx]


def _boxcar_spatial_filter(
    sites: Any,
    *,
    window: int,
    trim: bool,
    component: str,
    frequency_rtol: float,
    inplace: bool,
    recursive: bool,
    on_dup: str,
    strict: bool,
    verbose: int,
) -> Any:
    """Apply a count-based EMAP moving average along station order."""
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    items = list(_iter_items(S))
    if not items:
        return S
    window = max(1, int(window))
    if window % 2 == 0:
        window += 1
    half = window // 2
    comps = _emap_component_indices(component)
    blocks = []
    for ed in items:
        blocks.append(_get_z_block(ed))
    name_to_index = {
        _name(ed, i): i
        for i, ed in enumerate(items)
    }

    def _neighbors(index: int) -> List[int]:
        lo = max(0, index - half)
        hi = min(len(items), index + half + 1)
        return list(range(lo, hi))

    def _filtered_row(index: int, freq: float) -> Optional[np.ndarray]:
        rows = []
        for j in _neighbors(index):
            Zj, zj, frj = blocks[j]
            if Zj is None:
                continue
            row = _nearest_frequency_row(
                zj,
                frj,
                freq,
                rtol=frequency_rtol,
            )
            if row is not None:
                rows.append(row)
        if not rows:
            return None
        arr = np.asarray(rows, dtype=complex)
        out = np.nanmean(arr, axis=0)
        if trim and arr.shape[0] >= 5:
            for a, b in comps:
                vals = arr[:, a, b]
                finite = np.isfinite(vals)
                if np.count_nonzero(finite) < 3:
                    continue
                keep_vals = vals[finite]
                order = np.argsort(np.abs(keep_vals))
                trimmed = keep_vals[order][1:-1]
                if trimmed.size:
                    out[a, b] = np.nanmean(trimmed)
        return out

    def _one(Si):
        ed = next(_iter_items(Si))
        index = name_to_index.get(_name(ed, 0), 0)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        z2 = z.copy()
        for k, freq in enumerate(fr):
            row = _filtered_row(index, float(freq))
            if row is None:
                continue
            for a, b in comps:
                z2[k, a, b] = row[a, b]
        Z.z = z2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def fixed_length_moving_average(
    sites: Any,
    *,
    window: int = 5,
    component: str = "all",
    frequency_rtol: float = 1e-6,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    """Apply a fixed-length EMAP moving average along a profile.

    The filter replaces each selected impedance component by the local
    arithmetic mean of neighboring stations at the same frequency.  It is a
    v2 functional equivalent of the classic FLMA idea and preserves the
    input frequency grids.
    """
    return _boxcar_spatial_filter(
        sites,
        window=window,
        trim=False,
        component=component,
        frequency_rtol=frequency_rtol,
        inplace=inplace,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )


def trimmed_moving_average(
    sites: Any,
    *,
    window: int = 5,
    component: str = "all",
    frequency_rtol: float = 1e-6,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    """Apply a trimmed EMAP moving average along a profile.

    The filter is similar to :func:`fixed_length_moving_average`, but when a
    full enough window is available it removes the smallest and largest
    magnitudes before averaging.  This makes the profile smoothing less
    sensitive to isolated bad stations.
    """
    return _boxcar_spatial_filter(
        sites,
        window=window,
        trim=True,
        component=component,
        frequency_rtol=frequency_rtol,
        inplace=inplace,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )


def apply_emap_filter(
    sites: Any,
    *,
    method: str = "ama",
    window: int = 5,
    window_m: float = 1500.0,
    spacing_m: float = 200.0,
    component: str = "all",
    comp: str = "det",
    frequency_rtol: float = 1e-6,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    """Apply one EMAP-style spatial filter to MT/AMT sites.

    ``method='ama'`` delegates to :func:`correct_static_shift`, the existing
    Hanning adaptive moving-average correction.  ``'flma'`` and ``'tma'``
    apply count-based spatial smoothing along station order.
    """
    method = str(method).lower()
    if method in {"ama", "adaptive"}:
        return correct_static_shift(
            sites,
            window_m=window_m,
            spacing_m=spacing_m,
            comp=comp,
            inplace=inplace,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
    if method in {"flma", "fixed"}:
        return fixed_length_moving_average(
            sites,
            window=window,
            component=component,
            frequency_rtol=frequency_rtol,
            inplace=inplace,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
    if method in {"tma", "trimmed"}:
        return trimmed_moving_average(
            sites,
            window=window,
            component=component,
            frequency_rtol=frequency_rtol,
            inplace=inplace,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
    msg = "method must be one of 'ama', 'flma', or 'tma'."
    raise ValueError(msg)


def _confidence_vector_for_station(
    table: pd.DataFrame,
    station: str,
    fr: np.ndarray,
) -> np.ndarray:
    """Return confidence values aligned to one station frequency grid."""
    out = np.full(fr.size, np.nan, dtype=float)
    if table.empty:
        return out
    sub = table[table["station"].astype(str) == str(station)]
    if sub.empty:
        return out
    ftab = sub["frequency_hz"].to_numpy(dtype=float)
    ctab = sub["confidence"].to_numpy(dtype=float)
    for i, freq in enumerate(fr):
        if not np.isfinite(freq) or ftab.size == 0:
            continue
        idx = int(np.nanargmin(np.abs(ftab - freq)))
        if np.isclose(ftab[idx], freq, rtol=1e-6, atol=1e-12):
            out[i] = ctab[idx]
    return out


def _emap_decision_report(decisions: pd.DataFrame) -> pd.DataFrame:
    """Aggregate confidence-gated EMAP decisions by station."""
    if decisions.empty:
        return pd.DataFrame(
            columns=[
                "station",
                "n_freq",
                "n_preserved",
                "n_blended",
                "n_filtered",
                "mean_blend_weight",
                "median_confidence",
                "median_delta_log10_abs_z",
            ]
        )
    grouped = decisions.groupby("station", sort=False)
    rows = []
    for station, sub in grouped:
        rows.append(
            dict(
                station=station,
                n_freq=int(len(sub)),
                n_preserved=int((sub["action"] == "preserved").sum()),
                n_blended=int((sub["action"] == "blended").sum()),
                n_filtered=int((sub["action"] == "filtered").sum()),
                mean_blend_weight=float(
                    np.nanmean(sub["blend_weight"].to_numpy(dtype=float)),
                ),
                median_confidence=float(
                    np.nanmedian(sub["confidence"].to_numpy(dtype=float)),
                ),
                median_delta_log10_abs_z=float(
                    np.nanmedian(
                        sub["delta_log10_abs_z"].to_numpy(dtype=float),
                    ),
                ),
            )
        )
    return pd.DataFrame.from_records(rows)


def confidence_gated_emap_filter(
    sites: Any,
    *,
    before_sites: Optional[Any] = None,
    method: str = "flma",
    confidence_method: str = "composite",
    component: str = "xy",
    ci_hi: float = 0.90,
    ci_lo: float = 0.50,
    weights: Optional[Dict[str, float]] = None,
    blend_power: float = 1.0,
    window: int = 5,
    window_m: float = 1500.0,
    spacing_m: float = 200.0,
    comp: str = "det",
    frequency_rtol: float = 1e-6,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> EMAPFilterResult:
    """Apply EMAP filtering only as strongly as confidence requires.

    Rows with confidence greater than or equal to ``ci_hi`` are preserved.
    Rows below ``ci_lo`` are fully replaced by the EMAP-filtered estimate.
    Rows between the two limits are linearly blended, with optional
    ``blend_power`` shaping.
    """
    from .qc import frequency_confidence_table

    baseline = before_sites if before_sites is not None else sites
    before = ensure_sites(
        baseline,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    filtered = apply_emap_filter(
        sites,
        method=method,
        window=window,
        window_m=window_m,
        spacing_m=spacing_m,
        component=component,
        comp=comp,
        frequency_rtol=frequency_rtol,
        inplace=False,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    conf_table = frequency_confidence_table(
        before,
        method=confidence_method,
        weights=weights,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        spacing_m=spacing_m,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    before_map = {
        _name(ed, i): ed
        for i, ed in enumerate(_iter_items(before))
    }
    comps = _emap_component_indices(component)
    decisions = []
    for i, edf in enumerate(_iter_items(filtered)):
        station = _name(edf, i)
        ed0 = before_map.get(station)
        if ed0 is None:
            continue
        Z0, z0, fr0 = _get_z_block(ed0)
        Zf, zf, frf = _get_z_block(edf)
        if Z0 is None or Zf is None:
            continue
        conf = _confidence_vector_for_station(conf_table, station, fr0)
        znew = zf.copy()
        for j, freq in enumerate(frf):
            if fr0.size == 0:
                continue
            idx0 = int(np.nanargmin(np.abs(fr0 - freq)))
            if not np.isclose(fr0[idx0], freq, rtol=frequency_rtol, atol=1e-12):
                continue
            ci = conf[idx0]
            if not np.isfinite(ci):
                alpha = 0.0
            else:
                alpha = (float(ci_hi) - ci) / max(
                    float(ci_hi) - float(ci_lo),
                    1e-12,
                )
                alpha = float(np.clip(alpha, 0.0, 1.0))
                alpha = alpha ** max(float(blend_power), 1e-12)
            action = (
                "preserved"
                if alpha <= 1e-9
                else "filtered"
                if alpha >= 1.0 - 1e-9
                else "blended"
            )
            before_mag = []
            after_mag = []
            for a, b in comps:
                original = z0[idx0, a, b]
                filtered_value = zf[j, a, b]
                znew[j, a, b] = (
                    (1.0 - alpha) * original
                    + alpha * filtered_value
                )
                before_mag.append(abs(original))
                after_mag.append(abs(znew[j, a, b]))
            before_ref = float(np.nanmedian(before_mag))
            after_ref = float(np.nanmedian(after_mag))
            decisions.append(
                dict(
                    station=station,
                    frequency_hz=float(freq),
                    period_s=float(1.0 / freq) if freq else np.nan,
                    log10_period=(
                        float(np.log10(1.0 / freq))
                        if freq > 0
                        else np.nan
                    ),
                    confidence=float(ci) if np.isfinite(ci) else np.nan,
                    blend_weight=float(alpha),
                    action=action,
                    delta_log10_abs_z=float(
                        np.log10(after_ref + 1e-24)
                        - np.log10(before_ref + 1e-24)
                    ),
                )
            )
        Zf.z = znew
    decision_table = pd.DataFrame.from_records(decisions)
    return EMAPFilterResult(
        sites=filtered,
        report=_emap_decision_report(decision_table),
        decisions=decision_table,
        method=str(method),
        confidence_method=str(confidence_method),
        ci_hi=float(ci_hi),
        ci_lo=float(ci_lo),
    )


# -------------------- NOISE REMOVAL QC PLOTS ---------------------------- #

def _offdiag_logmag(z: np.ndarray) -> np.ndarray:
    m = np.nanmedian(
        np.stack([np.abs(z[:, 0, 1]), np.abs(z[:, 1, 0])], axis=1),
        axis=1,
    )
    return np.log10(np.maximum(m, 1e-24))


def _union_offdiag_matrix(
    sites: Any,
) -> Tuple[List[str], np.ndarray, np.ndarray]:
    S = ensure_sites(sites, recursive=False, strict=False)
    sts, Gs, Ms = [], [], []
    for i, ed in enumerate(_iter_items(S)):
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        sts.append(_name(ed, i))
        Gs.append(fr)
        Ms.append(_offdiag_logmag(z))
    if not sts:
        return [], np.array([]), np.empty((0, 0))
    G = np.unique(np.concatenate(Gs))
    M = np.full((len(sts), G.size), np.nan, dtype=float)
    for row, (fr, lm) in enumerate(zip(Gs, Ms)):
        idx = np.searchsorted(G, fr)
        idx = np.clip(idx, 0, G.size - 1)
        M[row, idx] = lm
        # fill small gaps by nearest neighbors
        r = M[row]
        g = np.isfinite(r)
        if g.sum() >= 2:
            xi = np.where(g)[0]
            for j in np.where(~g)[0]:
                k = np.searchsorted(xi, j)
                k0 = max(0, min(k - 1, xi.size - 1))
                k1 = max(0, min(k, xi.size - 1))
                a, b = xi[k0], xi[k1]
                t = 0.0 if a == b else (j - a) / (b - a)
                r[j] = (1 - t) * r[a] + t * r[b]
            M[row] = r
    return sts, G, M


def _union_component_logmag_matrix(
    sites: Any,
    component: str,
) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """Return station, union-frequency, log10 selected-component matrix."""
    S = ensure_sites(sites, recursive=False, strict=False)
    a, b = _component_index(component)
    sts, Gs, Ms = [], [], []
    for i, ed in enumerate(_iter_items(S)):
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        sts.append(_name(ed, i))
        Gs.append(fr)
        Ms.append(np.log10(np.abs(z[:, a, b]) + 1e-24))
    if not sts:
        return [], np.array([]), np.empty((0, 0))
    G = np.unique(np.concatenate(Gs))
    M = np.full((len(sts), G.size), np.nan, dtype=float)
    for row, (fr, lm) in enumerate(zip(Gs, Ms)):
        idx = np.searchsorted(G, fr)
        idx = np.clip(idx, 0, G.size - 1)
        M[row, idx] = lm
        r = M[row]
        g = np.isfinite(r)
        if g.sum() >= 2:
            xi = np.where(g)[0]
            r[~g] = np.interp(np.where(~g)[0], xi, r[g])
            M[row] = r
    return sts, G, M


def _denoise_sites(
    sites: Any,
    method: str,
    **kws: Any,
):
    m = method.lower()
    if m == "pipeline":
        return remove_noise_pipeline(sites, **kws)
    if m == "notch":
        return notch_powerline(sites, **kws)
    if m == "smooth":
        return smooth_logfreq(sites, **kws)
    if m == "rpca":
        return rpca_offdiag_denoise(sites, **kws)
    if m == "spatial":
        return spatial_median_filter(sites, **kws)
    if m == "hampel":
        return hampel_filter_freq(sites, **kws)
    if m in {"emap", "ama", "flma", "tma"}:
        emap_method = kws.pop("emap_method", "ama" if m == "emap" else m)
        return apply_emap_filter(sites, method=emap_method, **kws)
    raise ValueError(f"unknown method: {method}")


def _station_style_top(
    ax: plt.Axes,
    labels: Sequence[str],
    *,
    station_label_step: Optional[int],
    station_preset: str,
    station_style: Optional[Any],
) -> None:
    """Apply top station rendering for EMTools profile plots."""
    import copy
    from pycsamt.api.station import PYCSAMT_STATION_RENDERING

    style = station_style or PYCSAMT_STATION_RENDERING.style_for(
        station_preset,
    )
    style = copy.copy(style)
    style.side = "top"
    style.max_labels = max(int(style.max_labels), len(labels))
    style.every = 1 if station_label_step is None else int(
        station_label_step,
    )
    x = np.arange(len(labels), dtype=float)
    style.apply(ax, x, labels, xlim=(-0.5, len(labels) - 0.5))


def _component_index(component: str) -> Tuple[int, int]:
    """Return one tensor component index."""
    component = str(component).lower()
    mapping = {
        "xx": (0, 0),
        "xy": (0, 1),
        "yx": (1, 0),
        "yy": (1, 1),
    }
    if component not in mapping:
        msg = "component must be one of 'xx', 'xy', 'yx', or 'yy'."
        raise ValueError(msg)
    return mapping[component]


def _paired_z_items(before_sites: Any,
                    after_sites: Any) -> List[Tuple[str, Any, Any]]:
    """Pair before/after Z-bearing items by station name."""
    before = ensure_sites(before_sites, recursive=False, strict=False)
    after = ensure_sites(after_sites, recursive=False, strict=False)
    after_map = {
        _name(ed, i): ed
        for i, ed in enumerate(_iter_items(after))
    }
    pairs = []
    for i, ed0 in enumerate(_iter_items(before)):
        station = _name(ed0, i)
        ed1 = after_map.get(station)
        if ed1 is None:
            continue
        if _get_z_block(ed0)[0] is None or _get_z_block(ed1)[0] is None:
            continue
        pairs.append((station, ed0, ed1))
    return pairs


def _nearest_component_value(ed: Any,
                             freq_hz: Optional[float],
                             component: str) -> Tuple[float, complex]:
    """Return nearest frequency and complex tensor component magnitude."""
    _, z, fr = _get_z_block(ed)
    a, b = _component_index(component)
    if freq_hz is None:
        row = int(np.nanargmin(fr))
    else:
        row = int(np.nanargmin(np.abs(fr - float(freq_hz))))
    return float(fr[row]), complex(z[row, a, b])


def emap_filter_report(
    before_sites: Any,
    after_sites: Any,
    *,
    component: str = "xy",
    period_s: Optional[float] = None,
    frequency_hz: Optional[float] = None,
) -> pd.DataFrame:
    """Summarize station-level changes after an EMAP-style filter.

    The report compares the selected impedance component before and after a
    filter.  When ``period_s`` or ``frequency_hz`` is provided, the table also
    includes a reference-row before/after profile value for each station.
    """
    rows = []
    ref_freq = frequency_hz
    if ref_freq is None and period_s is not None:
        ref_freq = 1.0 / max(float(period_s), 1e-24)
    a, b = _component_index(component)
    for station, ed0, ed1 in _paired_z_items(before_sites, after_sites):
        _, z0, fr0 = _get_z_block(ed0)
        _, z1, fr1 = _get_z_block(ed1)
        vals = []
        for k, freq in enumerate(fr0):
            if fr1.size == 0:
                continue
            idx = int(np.nanargmin(np.abs(fr1 - freq)))
            if not np.isclose(fr1[idx], freq, rtol=1e-6, atol=1e-12):
                continue
            before = z0[k, a, b]
            after = z1[idx, a, b]
            if np.isfinite(before) and np.isfinite(after):
                vals.append(
                    np.log10(np.abs(after) + 1e-24)
                    - np.log10(np.abs(before) + 1e-24)
                )
        vals_arr = np.asarray(vals, dtype=float)
        row = dict(
            station=station,
            component=component,
            n_matched_freq=int(vals_arr.size),
            median_delta_log10_abs_z=(
                float(np.nanmedian(vals_arr))
                if vals_arr.size
                else np.nan
            ),
            rms_delta_log10_abs_z=(
                float(np.sqrt(np.nanmean(vals_arr ** 2)))
                if vals_arr.size
                else np.nan
            ),
        )
        if ref_freq is not None:
            f0, v0 = _nearest_component_value(ed0, ref_freq, component)
            f1, v1 = _nearest_component_value(ed1, ref_freq, component)
            row.update(
                reference_frequency_hz=float(f0),
                reference_frequency_after_hz=float(f1),
                before_log10_abs_z=float(np.log10(np.abs(v0) + 1e-24)),
                after_log10_abs_z=float(np.log10(np.abs(v1) + 1e-24)),
                reference_delta_log10_abs_z=float(
                    np.log10(np.abs(v1) + 1e-24)
                    - np.log10(np.abs(v0) + 1e-24)
                ),
            )
        rows.append(row)
    return pd.DataFrame.from_records(rows)


def plot_emap_filter_profile(
    before_sites: Any,
    after_sites: Optional[Any] = None,
    *,
    method: str = "flma",
    component: str = "xy",
    period_s: Optional[float] = None,
    frequency_hz: Optional[float] = None,
    window: int = 5,
    window_m: float = 1500.0,
    spacing_m: float = 200.0,
    comp: str = "det",
    figsize: Tuple[float, float] = (9.5, 4.0),
    station_label_step: Optional[int] = 1,
    station_preset: str = "pseudosection",
    station_style: Optional[Any] = None,
    ax: Optional[plt.Axes] = None,
    **filter_kws: Any,
) -> plt.Axes:
    """Plot a before/after EMAP filter station profile."""
    from pycsamt.api.style import PYCSAMT_STYLE

    if after_sites is None:
        after_sites = apply_emap_filter(
            before_sites,
            method=method,
            window=window,
            window_m=window_m,
            spacing_m=spacing_m,
            component=component,
            comp=comp,
            inplace=False,
            **filter_kws,
        )
    ref_freq = frequency_hz
    if ref_freq is None and period_s is not None:
        ref_freq = 1.0 / max(float(period_s), 1e-24)
    pairs = _paired_z_items(before_sites, after_sites)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if not pairs:
        ax.text(0.5, 0.5, "no paired stations", ha="center", va="center")
        return ax
    labels, before_vals, after_vals = [], [], []
    used_freqs = []
    for station, ed0, ed1 in pairs:
        f0, v0 = _nearest_component_value(ed0, ref_freq, component)
        f1, v1 = _nearest_component_value(ed1, f0, component)
        labels.append(station)
        before_vals.append(float(np.log10(np.abs(v0) + 1e-24)))
        after_vals.append(float(np.log10(np.abs(v1) + 1e-24)))
        used_freqs.append(f0 if np.isfinite(f0) else f1)
    x = np.arange(len(labels))
    ax.plot(
        x,
        before_vals,
        **PYCSAMT_STYLE.correction.before.plot_kwargs(),
    )
    ax.plot(
        x,
        after_vals,
        **PYCSAMT_STYLE.correction.after.plot_kwargs(),
    )
    _station_style_top(
        ax,
        labels,
        station_label_step=station_label_step,
        station_preset=station_preset,
        station_style=station_style,
    )
    freq_label = np.nanmedian(np.asarray(used_freqs, dtype=float))
    ax.set_ylabel(rf"$\log_{{10}}|Z_{{{component.upper()}}}|$")
    ax.set_title(
        f"{method.upper()} profile at {freq_label:.4g} Hz",
        fontsize=10,
    )
    ax.grid(True, ls=":", alpha=0.35)
    ax.legend(fontsize=8)
    return ax


def plot_emap_filter_psection(
    before_sites: Any,
    after_sites: Optional[Any] = None,
    *,
    method: str = "flma",
    component: str = "xy",
    window: int = 5,
    window_m: float = 1500.0,
    spacing_m: float = 200.0,
    comp: str = "det",
    cmap: str = "RdYlBu_r",
    delta_cmap: str = "RdBu_r",
    clim: Optional[Tuple[float, float]] = None,
    clim_pct: Tuple[float, float] = (2.0, 98.0),
    delta_vlim: Optional[float] = None,
    delta_vlim_pct: float = 95.0,
    axes=None,
    figsize: Tuple[float, float] = (11.0, 8.2),
    station_label_step: Optional[int] = 1,
    station_preset: str = "pseudosection",
    station_style: Optional[Any] = None,
    **filter_kws: Any,
) -> plt.Figure:
    """Plot before/after/delta pseudo-sections for an EMAP filter."""
    if after_sites is None:
        after_sites = apply_emap_filter(
            before_sites,
            method=method,
            window=window,
            window_m=window_m,
            spacing_m=spacing_m,
            component=component,
            comp=comp,
            inplace=False,
            **filter_kws,
        )
    st0, G0, M0 = _union_component_logmag_matrix(before_sites, component)
    st1, G1, M1 = _union_component_logmag_matrix(after_sites, component)
    axes_given = _axes_list(axes, 3) if axes is not None else None
    if axes_given is None:
        fig, axes_arr = plt.subplots(
            3,
            1,
            figsize=figsize,
            sharex=True,
            gridspec_kw={"hspace": 0.24},
        )
        axes_arr = np.asarray(axes_arr, dtype=object).ravel()
    else:
        axes_arr = np.asarray(axes_given, dtype=object)
        fig = axes_arr[0].figure
    if not st0 or not st1:
        axes_arr[1].text(0.5, 0.5, "no paired data", ha="center", va="center")
        return fig

    labels = [station for station in st0 if station in set(st1)]
    I0 = [st0.index(station) for station in labels]
    I1 = [st1.index(station) for station in labels]
    G = np.unique(np.concatenate([G0, G1]))

    def _resample(M: np.ndarray, Gs: np.ndarray) -> np.ndarray:
        out = np.full((M.shape[0], G.size), np.nan, dtype=float)
        for i in range(M.shape[0]):
            idx = np.searchsorted(G, Gs)
            idx = np.clip(idx, 0, G.size - 1)
            out[i, idx] = M[i]
            row = out[i]
            ok = np.isfinite(row)
            if ok.sum() >= 2:
                xi = np.where(ok)[0]
                row[~ok] = np.interp(np.where(~ok)[0], xi, row[ok])
        return out

    before = _resample(M0[I0], G0)
    after = _resample(M1[I1], G1)
    delta = after - before
    yvals = np.log10(1.0 / np.maximum(G, 1e-24))
    order = np.argsort(yvals)
    y_sorted = yvals[order]

    if clim is None:
        values = np.concatenate([before.ravel(), after.ravel()])
        values = values[np.isfinite(values)]
        if values.size:
            vmin, vmax = np.percentile(values, clim_pct)
        else:
            vmin, vmax = 0.0, 1.0
    else:
        vmin, vmax = float(clim[0]), float(clim[1])
    if delta_vlim is None:
        dvals = np.abs(delta[np.isfinite(delta)])
        delta_vlim = (
            float(np.nanpercentile(dvals, delta_vlim_pct))
            if dvals.size
            else 0.25
        )
        delta_vlim = max(delta_vlim, 1e-6)

    extent = (-0.5, len(labels) - 0.5, y_sorted.min(), y_sorted.max())
    panels = [
        (before[:, order].T, axes_arr[0], cmap, vmin, vmax, "Before"),
        (after[:, order].T, axes_arr[1], cmap, vmin, vmax, "After"),
        (
            delta[:, order].T,
            axes_arr[2],
            delta_cmap,
            -float(delta_vlim),
            float(delta_vlim),
            r"$\Delta$ after-before",
        ),
    ]
    images = []
    for data, ax, cm, lo, hi, title in panels:
        im = ax.imshow(
            data,
            aspect="auto",
            origin="lower",
            interpolation="nearest",
            cmap=cm,
            vmin=lo,
            vmax=hi,
            extent=extent,
        )
        images.append(im)
        ax.set_ylabel(LOG10_PERIOD_LABEL)
        ax.set_title(title, fontsize=9)
        ax.grid(False)
        if not ax.yaxis_inverted():
            ax.invert_yaxis()

    _station_style_top(
        axes_arr[0],
        labels,
        station_label_step=station_label_step,
        station_preset=station_preset,
        station_style=station_style,
    )
    for ax in axes_arr[1:]:
        ax.tick_params(
            axis="x",
            which="both",
            top=False,
            bottom=False,
            labeltop=False,
            labelbottom=False,
        )
    cb_main = fig.colorbar(
        images[0],
        ax=[axes_arr[0], axes_arr[1]],
        fraction=0.018,
        pad=0.012,
        aspect=35,
    )
    cb_main.set_label(rf"$\log_{{10}}|Z_{{{component.upper()}}}|$")
    cb_delta = fig.colorbar(
        images[2],
        ax=axes_arr[2],
        fraction=0.018,
        pad=0.012,
        aspect=35,
    )
    cb_delta.set_label(r"$\Delta\log_{10}|Z|$")
    return fig


# 1) Δ log|Z_off| pseudosection (station × log-period) ------------------- #

def nr_qc_delta_offdiag_psection(
    sites: Any,
    *,
    method: str = "pipeline",
    vlim: Optional[float] = None,
    figsize: Tuple[float, float] = (9.0, 4.8),
    ax: Optional[plt.Axes] = None,
    **denoise: Any,
) -> plt.Axes:
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _denoise_sites(S0, method, inplace=False, **denoise)
    st0, G0, M0 = _union_offdiag_matrix(S0)
    st1, G1, M1 = _union_offdiag_matrix(S1)
    if not st0 or not st1:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    # align stations intersection, and G union
    labs = sorted(list(set(st0) & set(st1)))
    I0 = [st0.index(s) for s in labs]
    I1 = [st1.index(s) for s in labs]
    G = np.unique(np.concatenate([G0, G1]))
    def _resample(M, Gs):
        out = np.full((M.shape[0], G.size), np.nan)
        for i in range(M.shape[0]):
            idx = np.searchsorted(G, Gs)
            idx = np.clip(idx, 0, G.size - 1)
            out[i, idx] = M[i]
            # fill as above
            r = out[i]; g = np.isfinite(r)
            if g.sum() >= 2:
                xi = np.where(g)[0]
                for j in np.where(~g)[0]:
                    k = np.searchsorted(xi, j)
                    k0 = max(0, min(k - 1, xi.size - 1))
                    k1 = max(0, min(k, xi.size - 1))
                    a, b = xi[k0], xi[k1]
                    t = 0 if a == b else (j - a) / (b - a)
                    r[j] = (1 - t) * r[a] + t * r[b]
        return out
    A = _resample(M0[I0], G0)
    B = _resample(M1[I1], G1)
    D = (B - A).T  # (freq, station)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    v = D[np.isfinite(D)]
    if vlim is None and v.size:
        vlim = float(max(0.1, np.nanpercentile(np.abs(v), 95)))
    lp = np.log10(1.0 / G)
    order = np.argsort(lp)
    lp = lp[order]
    D = D[order]
    im = ax.imshow(
        D,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-(vlim or 0.5),
        vmax=(vlim or 0.5),
    )
    ax.set_ylabel(LOG10_PERIOD_LABEL)
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(len(labs)),
        labs,
        preset="pseudosection",
        xlim=(-0.5, len(labs) - 0.5),
    )
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    yt = np.linspace(0, len(lp) - 1, num=min(8, len(lp)))
    yv = np.linspace(lp.min(), lp.max(), num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("Δ log10 |Z_off| (after−before)")
    return ax


# 2) SNR gain profile (per-station dB improvement) ----------------------- #

def nr_qc_snr_gain_profile(
    sites: Any,
    *,
    method: str = "pipeline",
    pband: Optional[Tuple[float, float]] = None,
    figsize: Tuple[float, float] = (8.6, 3.6),
    ax: Optional[plt.Axes] = None,
    **denoise: Any,
) -> plt.Axes:
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _denoise_sites(S0, method, inplace=False, **denoise)
    T0 = snr_table(S0)
    T1 = snr_table(S1)
    if T0.empty or T1.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no snr", ha="center", va="center")
        return ax
    def _band(v, f):
        if pband is None:
            return np.nanmedian(v)
        lo, hi = pband
        m = (1.0 / f >= lo) & (1.0 / f <= hi)
        return np.nanmedian(v[m]) if np.any(m) else np.nan
    labs, gain = [], []
    for st in sorted(set(T0["station"]) & set(T1["station"])):
        s0 = T0[T0["station"] == st]
        s1 = T1[T1["station"] == st]
        if s0.empty or s1.empty:
            continue
        f0 = s0["freq"].to_numpy(dtype=float)
        f1 = s1["freq"].to_numpy(dtype=float)
        sn0 = s0["snr"].to_numpy(dtype=float)
        sn1 = s1["snr"].to_numpy(dtype=float)
        g0 = _band(sn0, f0)
        g1 = _band(sn1, f1)
        if np.isfinite(g0) and np.isfinite(g1) and g0 > 0:
            labs.append(st)
            gain.append(20.0 * np.log10(g1 / g0))
    if not labs:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no gain", ha="center", va="center")
        return ax
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(labs))
    ax.axhline(0.0, color="0.7", lw=1.0)
    ax.bar(x, gain, width=0.8)
    ax.set_ylabel("SNR gain (dB)")
    ax.set_xlabel("Station")
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=90)
    return ax


# 3) Harmonic waterfall (Δ dB at mains harmonics) ------------------------ #

def nr_qc_harmonic_waterfall(
    sites: Any,
    *,
    method: str = "notch",
    mains_hz: float = 50.0,
    n_harm: int = 30,
    tol_hz: float = 0.08,
    figsize: Tuple[float, float] = (9.0, 4.6),
    ax: Optional[plt.Axes] = None,
    **denoise: Any,
) -> plt.Axes:
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _denoise_sites(
        S0, method, inplace=False,
        mains_hz=mains_hz, n_harm=n_harm, tol_hz=tol_hz,
        **denoise,
    )
    sts = []
    for i, (S, tag) in enumerate([(S0, "b"), (S1, "a")]):
        rows = []
        labs = []
        for j, ed in enumerate(_iter_items(S)):
            Z, z, fr = _get_z_block(ed)
            if Z is None:
                continue
            lm = _offdiag_logmag(z)
            kk = np.arange(1, int(n_harm) + 1, dtype=float)
            hv = []
            for k in kk:
                m = np.abs(fr - k * mains_hz) <= float(tol_hz)
                hv.append(
                    float(np.nanmedian(lm[m])) if np.any(m) else np.nan
                )
            rows.append(hv)
            if i == 0:
                labs.append(_name(ed, j))
        if i == 0:
            sts = labs
            Hb = np.array(rows, float)
        else:
            Ha = np.array(rows, float)
    if not sts:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    D = (Hb - Ha).T  # reduction (dB-like log units)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        D,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap="viridis",
    )
    ax.set_ylabel("Harmonic index k (k·mains)")
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(len(sts)),
        sts,
        preset="pseudosection",
        xlim=(-0.5, len(sts) - 0.5),
    )
    yt = np.linspace(0, D.shape[0] - 1, num=min(8, D.shape[0]))
    yv = np.linspace(1, D.shape[0], num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{int(v)}" for v in yv])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("Δ log10 |Z_off| (reduction)")
    return ax


# 4) Station curves (off-diag mag with mains guides) --------------------- #

def nr_qc_station_offdiag_curves(
    sites: Any,
    *,
    method: str = "pipeline",
    station: Optional[str] = None,
    mains_hz: float = 50.0,
    n_harm: int = 12,
    tol_hz: float = 0.08,
    figsize: Tuple[float, float] = (8.0, 4.2),
    ax: Optional[plt.Axes] = None,
    **denoise: Any,
) -> plt.Axes:
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _denoise_sites(S0, method, inplace=False, **denoise)
    # pick station
    fm = {}
    for i, ed in enumerate(_iter_items(S0)):
        fm[_name(ed, i)] = ed
    if station is None and fm:
        station = sorted(fm.keys())[0]
    ed0 = fm.get(station, None)
    if ed0 is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "station not found",
                ha="center", va="center")
        return ax
    # match from S1
    for ed in _iter_items(S1):
        if _name(ed, 0) == station:
            ed1 = ed; break
    else:
        ed1 = None
    if ed1 is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no after", ha="center", va="center")
        return ax
    Z0, z0, f0 = _get_z_block(ed0)
    Z1, z1, f1 = _get_z_block(ed1)
    if Z0 is None or Z1 is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no Z", ha="center", va="center")
        return ax
    m0 = 10 ** _offdiag_logmag(z0)
    m1 = 10 ** _offdiag_logmag(z1)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.set_xscale("log")
    ax.plot(1.0 / f0, m0, "o-", lw=1.4, label="before")
    ax.plot(1.0 / f1, m1, "s-", lw=1.4, label="after")
    # mains guides
    kk = np.arange(1, int(n_harm) + 1, dtype=float)
    for k in kk:
        f = k * mains_hz
        ax.axvspan(
            1.0 / (f + tol_hz),
            1.0 / max(1e-12, (f - tol_hz)),
            alpha=0.08, color="r",
        )
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("|Z_off|")
    ax.set_title(str(station))
    ax.grid(True, alpha=0.25, which="both")
    ax.legend()
    return ax
