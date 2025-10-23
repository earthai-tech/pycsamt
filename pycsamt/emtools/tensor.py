
from __future__ import annotations 

from typing import Optional, Tuple, List, Any
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import colors as mcolors


from pycsamt.site.base import Sites
from pycsamt.site import edit as _edit
from pycsamt.site import compute as _compute
from pycsamt.z import utils as zutils

from ._core import (
    ensure_sites,
    _iter_items,
    _name,
    _get_z_block,
    _apply_each
)

_BACKWARD_SINCE = "2.0.0"
_BACKWARD_REMOVE = "2.17.0"

#
def rotate_to_strike(
    sites,
    *,
    method: str = "swift",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ang = 0.0
        try:
            r = _compute.strike_estimate(Si, method=method)
            if isinstance(r, (int, float)):
                ang = float(r)
            elif isinstance(r, dict):
                ang = float(r.get("angle", 0.0))
            else:
                ang = float(getattr(r, "angle", 0.0))
        except Exception:
            ang = 0.0
        return _edit.rotate(
            Si,
            angle=ang,
            inplace=inplace,
        )

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 2) rotate by angle or map ------------------------------------------- #

def rotate(
    sites,
    angle: float,
    *,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    return _edit.rotate(S, angle=angle, inplace=inplace)


def rotate_by_map(
    sites,
    angle_by_station: dict,
    *,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        name = getattr(ed, "station", None) or getattr(ed, "name", None)
        ang = float(angle_by_station.get(name, 0.0))
        return _edit.rotate(Si, angle=ang, inplace=inplace)

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 3) enforce off-diagonal antisymmetry -------------------------------- #

def antisymmetrize(
    sites,
    *,
    how: str = "rms",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        Z = getattr(ed, "Z", None) or getattr(ed, "z", None)
        if Z is None:
            return Si
        z = getattr(Z, "z", None)
        ze = getattr(Z, "z_err", None)
        if z is None:
            return Si
        z2, ze2 = zutils.enforce_offdiag_antisymmetry(
            z, ze, how=how
        )
        try:
            Z.z = z2
            if ze is not None and ze2 is not None:
                Z.z_err = ze2
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 4) invert impedance -------------------------------------------------- #

def invert(
    sites,
    *,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        Z = getattr(ed, "Z", None) or getattr(ed, "z", None)
        if Z is None:
            return Si
        z = getattr(Z, "z", None)
        ze = getattr(Z, "z_err", None)
        if z is None:
            return Si
        z2, ze2 = zutils.invert_z(z, ze)
        try:
            Z.z = z2
            if ze is not None and ze2 is not None:
                Z.z_err = ze2
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 5) sensor-orientation correction ------------------------------------ #

def orient_from_sensors(
    sites,
    ex: float,
    ey: float,
    bx: float,
    by: float,
    *,
    degrees: bool = True,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        Z = getattr(ed, "Z", None) or getattr(ed, "z", None)
        if Z is None:
            return Si
        z = getattr(Z, "z", None)
        ze = getattr(Z, "z_err", None)
        if z is None:
            return Si
        z2, ze2 = zutils.correct_for_sensor_orientation(
            z,
            ex=ex, ey=ey, bx=bx, by=by,
            degrees=degrees,
            z_err=ze,
        )
        try:
            Z.z = z2
            if ze is not None and ze2 is not None:
                Z.z_err = ze2
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 6) sigma-clip outliers in Z ----------------------------------------- #

def sigma_clip_z(
    sites,
    *,
    sigma: float = 3.0,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        Z = getattr(ed, "Z", None) or getattr(ed, "z", None)
        if Z is None:
            return Si
        z = getattr(Z, "z", None)
        if z is None:
            return Si
        mask = zutils.sigma_clip_mask(z, sigma=sigma)
        z2 = z.copy()
        z2[~mask] = np.nan
        try:
            Z.z = z2
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- 7) balance off-diagonal magnitudes ---------------------------------- #

def balance_offdiag(
    sites,
    *,
    mode: str = "avgabs",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Sites:

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Sites) -> Sites:
        ed = next(_iter_items(Si))
        Z = getattr(ed, "Z", None) or getattr(ed, "z", None)
        if Z is None:
            return Si
        z = getattr(Z, "z", None)
        if z is None:
            return Si
        zx = z[:, 0, 1]
        zy = z[:, 1, 0]
        if mode == "avgabs":
            m = 0.5 * (np.abs(zx) + np.abs(zy)) + 1e-12
            sx = m / (np.abs(zx) + 1e-12)
            sy = m / (np.abs(zy) + 1e-12)
            z2 = z.copy()
            z2[:, 0, 1] = zx * sx
            z2[:, 1, 0] = zy * sy
        else:
            z2 = z
        try:
            Z.z = z2
        except Exception:
            pass
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)




def _get_z(ed: Any) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    Zobj = getattr(ed, "Z", None) or getattr(ed, "z", None)
    if Zobj is None:
        return None, None
    z = getattr(Zobj, "z", None)
    fr = getattr(Zobj, "freq", None)
    if isinstance(z, np.ndarray) and isinstance(fr, np.ndarray):
        return z, fr
    if isinstance(Zobj, np.ndarray) and Zobj.ndim == 3:
        fr2 = getattr(ed, "freq", None)
        if isinstance(fr2, np.ndarray):
            return Zobj, fr2
    return None, None


def _period(fr: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        p = 1.0 / fr
    return p


def _angles_deg(a: np.ndarray, b: np.ndarray,
                c: np.ndarray, d: np.ndarray) -> Tuple[
                    np.ndarray, np.ndarray]:
    num_b = b + c
    den_b = a - d
    beta = 0.5 * np.degrees(np.arctan2(num_b, den_b))
    num_a = -(b - c)
    den_a = a + d
    alpha = 0.5 * np.degrees(np.arctan2(num_a, den_a))
    return alpha, beta


def _phi_from_z(z: np.ndarray) -> np.ndarray:
    # z: (n, 2, 2) complex
    X = z.real
    Y = z.imag
    n = z.shape[0]
    phi = np.empty((n, 2, 2), dtype=float)
    for i in range(n):
        Xi = X[i]
        Yi = Y[i]
        try:
            invX = np.linalg.inv(Xi)
        except np.linalg.LinAlgError:
            invX = np.linalg.pinv(Xi)
        phi[i] = invX @ Yi
    return phi


def _svd_axes(phi: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # returns s (n,2) and angle of major axis (deg), shape (n,)
    n = phi.shape[0]
    s = np.empty((n, 2), dtype=float)
    ang = np.empty(n, dtype=float)
    for i in range(n):
        U, S, VT = np.linalg.svd(phi[i])
        s[i] = S
        vx = U[0, 0]
        vy = U[1, 0]
        ang[i] = np.degrees(np.arctan2(vy, vx))
    return s, ang


def _ellipticity(s: np.ndarray) -> np.ndarray:
    num = s[:, 0] - s[:, 1]
    den = s[:, 0] + s[:, 1] + 1e-12
    return num / den


def _df_phase_tensor_single(
    st: str,
    fr: np.ndarray,
    phi: np.ndarray,
    s: np.ndarray,
    theta: np.ndarray,
    alpha: np.ndarray,
    beta: np.ndarray,
) -> pd.DataFrame:
    p = _period(fr)
    e = _ellipticity(s)
    rec = dict(
        station=[st] * len(fr),
        freq=fr,
        period=p,
        s1=s[:, 0],
        s2=s[:, 1],
        theta=theta,
        alpha=alpha,
        beta=beta,
        skew=beta,
        ellipt=e,
    )
    return pd.DataFrame(rec)


# ------------------------ public computation ----------------------------- #

def build_phase_tensor_table(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    dfs: List[pd.DataFrame] = []
    for i, ed in enumerate(_iter_items(S)):
        nm = _name(ed, i)
        z, fr = _get_z(ed)
        if z is None or fr is None:
            continue
        phi = _phi_from_z(z)
        a = phi[:, 0, 0]
        b = phi[:, 0, 1]
        c = phi[:, 1, 0]
        d = phi[:, 1, 1]
        alpha, beta = _angles_deg(a, b, c, d)
        s, theta = _svd_axes(phi)
        df = _df_phase_tensor_single(
            nm, fr, phi, s, theta, alpha, beta
        )
        dfs.append(df)
    if not dfs:
        return pd.DataFrame(
            columns=[
                "station", "freq", "period", "s1", "s2",
                "theta", "alpha", "beta", "skew", "ellipt",
            ]
        )
    return pd.concat(dfs, ignore_index=True)


# ------------------------------ plotting --------------------------------- #

def plot_phase_tensor_psection(
    sites: Any,
    *,
    axis_y: str = "logperiod",
    scale: float = 0.12,
    c_by: str = "skew",
    figsize: Tuple[float, float] = (9.0, 5.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    st_list = list(df["station"].unique())
    st_map = {s: i for i, s in enumerate(st_list)}
    x = df["station"].map(st_map).to_numpy()
    if axis_y == "logperiod":
        y = np.log10(df["period"].to_numpy())
        ylab = "LogPeriod (s)"
    else:
        y = df["period"].to_numpy()
        ylab = "Period (s)"
    w = scale * df["s1"].to_numpy()
    h = scale * df["s2"].to_numpy()
    ang = df["theta"].to_numpy()
    cval = df.get(c_by, df["skew"]).to_numpy()
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    vmin = np.nanpercentile(cval, 5)
    vmax = np.nanpercentile(cval, 95)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("coolwarm")
    for xi, yi, wi, hi, ai, ci in zip(
        x, y, w, h, ang, cval
    ):
        if not np.all(np.isfinite([xi, yi, wi, hi, ai, ci])):
            continue
        e = Ellipse(
            (xi, yi),
            width=wi,
            height=hi,
            angle=ai,
            fill=True,
            lw=0.5,
        )
        e.set_facecolor(cmap(norm(ci)))
        e.set_edgecolor("k")
        ax.add_patch(e)
    ax.set_xlim(-0.5, len(st_list) - 0.5)
    ymin = np.nanmin(y) if np.isfinite(y).any() else 0.0
    ymax = np.nanmax(y) if np.isfinite(y).any() else 1.0
    pad = 0.02 * (ymax - ymin + 1e-9)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_xlabel("Station")
    ax.set_ylabel(ylab)
    ax.set_xticks(np.arange(len(st_list)))
    ax.set_xticklabels(st_list, rotation=90)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = plt.colorbar(sm, ax=ax)
    cb.set_label(c_by)
    return ax


def plot_phase_tensor_skewmap(
    sites: Any,
    *,
    axis_y: str = "logperiod",
    agg: str = "median",
    figsize: Tuple[float, float] = (9.0, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    if axis_y == "logperiod":
        df = df.copy()
        df["y"] = np.log10(df["period"].to_numpy())
        ylab = "LogPeriod (s)"
    else:
        df = df.copy()
        df["y"] = df["period"].to_numpy()
        ylab = "Period (s)"
    piv = df.pivot_table(
        index="y",
        columns="station",
        values="skew",
        aggfunc=agg,
    )
    piv = piv.sort_index()
    Z = piv.to_numpy(dtype=float).T
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel(ylab)
    ax.set_xticks(np.arange(len(piv.columns)))
    ax.set_xticklabels(piv.columns, rotation=90)
    yt = np.linspace(0, Z.shape[0] - 1, num=min(8, Z.shape[0]))
    yvals = np.linspace(piv.index.min(), piv.index.max(),
                        num=min(8, len(piv.index)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.3g}" for v in yvals])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("skew")
    return ax

# --- more plots ----------------------------------------------------------- #

def _coords_of_sites(sites) -> dict:
    m = {}
    for i, ed in enumerate(_iter_items(sites)):
        nm = _name(ed, i)
        lat = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
        lon = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
        try:
            if lat is None or lon is None:
                continue
            m[nm] = (float(lat), float(lon))
        except Exception:
            continue
    return m


def plot_theta_vs_period(
    sites: Any,
    *,
    figsize: Tuple[float, float] = (8.0, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df.empty:
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    for st, sdf in df.groupby("station"):
        x = sdf["period"].to_numpy()
        y = sdf["theta"].to_numpy()
        ax.plot(x, y, ".", ms=3, label=st)
    ax.set_xscale("log")
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("theta (deg)")
    ax.legend(ncols=2, fontsize=8)
    return ax


def plot_ellipticity_psection(
    sites: Any,
    *,
    figsize: Tuple[float, float] = (8.5, 4.0),
    agg: str = "median",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df.empty:
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    df = df.copy()
    df["logp"] = np.log10(df["period"].to_numpy())
    piv = df.pivot_table(
        index="logp",
        columns="station",
        values="ellipt",
        aggfunc=agg,
    )
    piv = piv.sort_index()
    Z = piv.to_numpy(dtype=float).T
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(Z.shape[0]))
    ax.set_xticklabels(piv.columns, rotation=90)
    yt = np.linspace(0, Z.shape[1] - 1, num=min(8, Z.shape[1]))
    yvals = np.linspace(piv.index.min(), piv.index.max(),
                        num=min(8, len(piv.index)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yvals])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("ellipticity")
    return ax


def plot_dimensionality_psection(
    sites: Any,
    *,
    skew_th: float = 3.0,
    ellipt_th: float = 0.2,
    figsize: Tuple[float, float] = (8.5, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df.empty:
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    df = df.copy()
    df["logp"] = np.log10(df["period"].to_numpy())
    # 0=1D, 1=2D, 2=3D (simple rule)
    lab = np.zeros(len(df), dtype=int)
    lab[(np.abs(df["skew"]) <= skew_th) &
        (df["ellipt"].abs() <= ellipt_th)] = 0
    lab[(np.abs(df["skew"]) <= skew_th) &
        (df["ellipt"].abs() > ellipt_th)] = 1
    lab[(np.abs(df["skew"]) > skew_th)] = 2
    df["dim"] = lab
    piv = df.pivot_table(
        index="logp",
        columns="station",
        values="dim",
        aggfunc="median",
    )
    piv = piv.sort_index()
    Z = piv.to_numpy(dtype=float).T
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(Z.shape[0]))
    ax.set_xticklabels(piv.columns, rotation=90)
    yt = np.linspace(0, Z.shape[1] - 1, num=min(8, Z.shape[1]))
    yvals = np.linspace(piv.index.min(), piv.index.max(),
                        num=min(8, len(piv.index)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yvals])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("dim (0=1D,1=2D,2=3D)")
    return ax


def plot_phase_tensor_rose(
    sites: Any,
    *,
    band: Tuple[float, float] = (1.0, 10.0),
    bins: int = 24,
    figsize: Tuple[float, float] = (5.0, 5.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(
                subplot_kw=dict(polar=True),
                figsize=figsize,
            )
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    lo, hi = band
    sel = (df["period"] >= lo) & (df["period"] <= hi)
    th = np.radians(df.loc[sel, "theta"].to_numpy())
    th = (th + 2 * np.pi) % (2 * np.pi)
    if ax is None:
        _, ax = plt.subplots(
            subplot_kw=dict(polar=True),
            figsize=figsize,
        )
    bins = int(max(8, bins))
    edges = np.linspace(0, 2 * np.pi, bins + 1)
    hist, _ = np.histogram(th, bins=edges)
    ang = 0.5 * (edges[1:] + edges[:-1])
    ax.bar(ang, hist, width=edges[1] - edges[0])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_title(f"theta rose [{lo},{hi}] s")
    return ax

def plot_phase_tensor_map(
    sites: Any,
    *,
    period: float = 10.0,
    scale: float = 0.15,
    c_by: str = "skew",
    figsize: Tuple[float, float] = (8.0, 6.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    df = build_phase_tensor_table(
        S,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    if c_by in ("rho_det", "phi_det"):
        df = _attach_det_metrics(df, S)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df.empty:
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    # coords
    coords = {}
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        lat = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
        lon = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
        try:
            if lat is None or lon is None:
                continue
            coords[st] = (float(lat), float(lon))
        except Exception:
            continue
    if not coords:
        ax.text(0.5, 0.5, "no coords",
                ha="center", va="center")
        return ax
    # one row per station nearest to target period
    rows = []
    for st, sdf in df.groupby("station"):
        if st not in coords:
            continue
        p = sdf["period"].to_numpy()
        i = int(np.nanargmin(np.abs(p - period)))
        rows.append(sdf.iloc[i])
    if not rows:
        ax.text(0.5, 0.5, "no match @period",
                ha="center", va="center")
        return ax
    vv = np.array([r.get(c_by, np.nan) for r in rows], float)
    v0 = np.nanpercentile(vv, 5)
    v1 = np.nanpercentile(vv, 95)
    norm = plt.Normalize(vmin=v0, vmax=v1)
    cmap = plt.get_cmap("coolwarm")
    for row in rows:
        st = row["station"]
        lat, lon = coords[st]
        s1 = float(row["s1"])
        s2 = float(row["s2"])
        th = float(row["theta"])
        w = scale * s1
        h = scale * s2
        e = Ellipse(
            (lon, lat),
            width=w,
            height=h,
            angle=th,
            fill=True,
            lw=0.5,
        )
        e.set_edgecolor("k")
        e.set_facecolor(cmap(norm(float(row.get(c_by, np.nan)))))
        ax.add_patch(e)
        ax.plot(lon, lat, "k.", ms=2)
    ax.set_xlabel("Lon")
    ax.set_ylabel("Lat")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = plt.colorbar(sm, ax=ax)
    cb.set_label(c_by)
    return ax



def phase_tensor_legend(
    *,
    size: float = 1.0,
    figsize: Tuple[float, float] = (2.5, 2.5),
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    e = Ellipse((0.0, 0.0), width=size, height=size,
                angle=0.0, fill=False, lw=1.0)
    ax.add_patch(e)
    ax.plot([0, 0], [0, size * 0.6], "-", lw=1.0)
    ax.set_xlim(-size, size)
    ax.set_ylim(-size, size)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    return ax



# -------------------- small utilities (local) --------------------------- #

def _rho_det_from_z(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    zx = z[:, 0, 1]
    zy = z[:, 1, 0]
    rx = 0.2 * (np.abs(zx) ** 2) / (fr + 1e-24)
    ry = 0.2 * (np.abs(zy) ** 2) / (fr + 1e-24)
    rdet = np.sqrt(rx * ry)
    return rdet


def _phi_det_from_z(z: np.ndarray) -> np.ndarray:
    d = z[:, 0, 0] * z[:, 1, 1] - z[:, 0, 1] * z[:, 1, 0]
    return np.degrees(np.angle(d))


def _attach_det_metrics(
    df: pd.DataFrame,
    S,
) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df["rho_det"] = np.nan
    df["phi_det"] = np.nan
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        rho = _rho_det_from_z(z, fr)
        phi = _phi_det_from_z(z)
        per = 1.0 / fr
        m = df["station"] == st
        if not m.any():
            continue
        p = df.loc[m, "period"].to_numpy()
        idx = np.searchsorted(per, p)
        idx = np.clip(idx, 0, len(per) - 1)
        df.loc[m, "rho_det"] = rho[idx]
        df.loc[m, "phi_det"] = phi[idx]
    return df


def _color_from_hsv(h: np.ndarray, s: np.ndarray,
                    v: np.ndarray) -> np.ndarray:
    hsv = np.stack([h, s, v], axis=-1)
    rgb = mcolors.hsv_to_rgb(hsv)
    return rgb

# ---------------- 2) dimensionality grid (1D/2D/3D) --------------------- #
def plot_dimensionality_grid(
    sites: Any,
    *,
    skew_th: float = 3.0,
    ellipt_th: float = 0.2,
    figsize: Tuple[float, float] = (8.5, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    df = df.copy()
    df["logp"] = np.log10(df["period"].to_numpy())
    lab = np.zeros(len(df), dtype=int)
    a = np.abs(df["skew"].to_numpy())
    e = np.abs(df["ellipt"].to_numpy())
    lab[:] = 2
    lab[(a <= skew_th) & (e <= ellipt_th)] = 0
    lab[(a <= skew_th) & (e > ellipt_th)] = 1
    df["dim"] = lab
    piv = df.pivot_table(
        index="logp",
        columns="station",
        values="dim",
        aggfunc="median",
    )
    piv = piv.sort_index()
    Z = piv.to_numpy(dtype=float).T
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(len(piv.columns)))
    ax.set_xticklabels(piv.columns, rotation=90)
    yt = np.linspace(0, Z.shape[1] - 1, num=min(8, Z.shape[1]))
    yv = np.linspace(piv.index.min(), piv.index.max(),
                     num=min(8, len(piv.index)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("dim (0=1D,1=2D,2=3D)")
    return ax


# ---------------- 3) theta stability stripe (HSV image) ----------------- #

def plot_theta_stability_stripe(
    sites: Any,
    *,
    win: int = 5,
    figsize: Tuple[float, float] = (9.0, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    df = df.copy()
    df["logp"] = np.log10(df["period"].to_numpy())
    sts = list(df["station"].unique())
    X = []
    H = []
    for st in sts:
        s = df[df["station"] == st].sort_values("logp")
        th = s["theta"].to_numpy()
        lp = s["logp"].to_numpy()
        # hue from theta in [0,180) → [0,1)
        h = ((th % 180.0) / 180.0)
        # variance in sliding window
        k = max(3, int(win))
        if len(h) < k:
            v = np.full_like(h, np.nan)
        else:
            v = np.convolve(
                (h - np.nanmean(h)) ** 2,
                np.ones(k) / k,
                mode="same",
            )
        # saturation from 1 - norm(var)
        v0 = np.nanpercentile(v[np.isfinite(v)], 5) if np.isfinite(
            v
        ).any() else 0.0
        v1 = np.nanpercentile(v[np.isfinite(v)], 95) if np.isfinite(
            v
        ).any() else 1.0
        s_sat = 1.0 - np.clip((v - v0) / (v1 - v0 + 1e-12), 0, 1)
        H.append(np.vstack([h, s_sat, np.ones_like(h)]))
        X.append(lp)
    # align to common y-grid for imshow
    yall = np.unique(np.concatenate(X))
    Himg = np.zeros((len(yall), len(sts), 3))
    for j, (lp, hs) in enumerate(zip(X, H)):
        i = np.searchsorted(yall, lp)
        i = np.clip(i, 0, len(yall) - 1)
        h, s, v = hs
        rgb = _color_from_hsv(h, s, v)
        for r, k in enumerate(i):
            Himg[k, j, :] = rgb[r]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.imshow(
        Himg,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(len(sts)))
    ax.set_xticklabels(sts, rotation=90)
    yt = np.linspace(0, len(yall) - 1, num=min(8, len(yall)))
    yv = np.linspace(yall.min(), yall.max(),
                     num=min(8, len(yall)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    return ax


# ---------------- 4) skew–ellipticity density (hexbin + contours) ------- #

def plot_skew_ellipt_density(
    sites: Any,
    *,
    band: Optional[Tuple[float, float]] = None,
    gridsize: int = 40,
    figsize: Tuple[float, float] = (6.5, 5.5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return ax
    a = np.abs(df["beta"].to_numpy())
    e = np.abs(df["ellipt"].to_numpy())
    if band is not None:
        lo, hi = band
        m = (df["period"] >= lo) & (df["period"] <= hi)
        a = a[m]
        e = e[m]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.hexbin(
        a, e, gridsize=gridsize, mincnt=1, linewidths=0.0
    )
    ax.set_xlabel("|beta| (deg)")
    ax.set_ylabel("|ellipticity|")
    # contours from 2D hist
    nx = ny = gridsize
    H, xedges, yedges = np.histogram2d(a, e, bins=[nx, ny])
    Xc = 0.5 * (xedges[1:] + xedges[:-1])
    Yc = 0.5 * (yedges[1:] + yedges[:-1])
    lev = np.linspace(H.max() * 0.2, H.max() * 0.9, 4)
    ax.contour(Xc, Yc, H.T, levels=lev, colors="k", linewidths=0.6)
    return ax


# ---------------- 5) theta rose per band (small multiples) -------------- #
def plot_theta_rose_grid(
    sites: Any,
    *,
    n_bands: int = 6,
    figsize: Tuple[float, float] = (10.0, 4.5),
    bins: int = 24,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, polar=True)
        ax.text(0.5, 0.5, "no phase tensor",
                ha="center", va="center")
        return fig
    p = df["period"].to_numpy()
    lo = float(np.nanmin(p))
    hi = float(np.nanmax(p))
    lo = max(lo, 1e-6)
    edges = np.logspace(np.log10(lo), np.log10(hi), n_bands + 1)
    bins = int(max(8, bins))
    fig = plt.figure(figsize=figsize)
    for i in range(n_bands):
        ax = fig.add_subplot(1, n_bands, i + 1, polar=True)
        m = (p >= edges[i]) & (p < edges[i + 1])
        th = np.radians(df.loc[m, "theta"].to_numpy())
        th = (th + 2 * np.pi) % (2 * np.pi)
        edges_a = np.linspace(0, 2 * np.pi, bins + 1)
        h, _ = np.histogram(th, bins=edges_a)
        ang = 0.5 * (edges_a[1:] + edges_a[:-1])
        ax.bar(ang, h, width=edges_a[1] - edges_a[0])
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_title(f"[{edges[i]:.2g},{edges[i+1]:.2g}]s",
                     fontsize=8)
    fig.tight_layout()
    return fig

