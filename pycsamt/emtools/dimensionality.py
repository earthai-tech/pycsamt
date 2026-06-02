# pycsamt/emtools/dimensionality.py
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import numpy as np
import pandas as pd

from ..site import edit as _edit
from ._core import (
    ensure_sites,
    _apply_each,
    _iter_items,
    _name,
    _get_t_block, 
    _get_z_block, 
)
from .tensor import build_phase_tensor_table

# -------------------------- local helpers ------------------------------- #

def _det_phase_from_z(z: np.ndarray) -> np.ndarray:
    detz = z[:, 0, 0] * z[:, 1, 1] - z[:, 0, 1] * z[:, 1, 0]
    return np.degrees(np.angle(detz))


def _rho_det_from_z(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    zx = z[:, 0, 1]
    zy = z[:, 1, 0]
    rx = 0.2 * (np.abs(zx) ** 2) / (fr + 1e-24)
    ry = 0.2 * (np.abs(zy) ** 2) / (fr + 1e-24)
    rdet = np.sqrt(rx * ry)
    return rdet


def _tip_amp(t: Optional[np.ndarray]) -> Optional[np.ndarray]:
    if t is None:
        return None
    return np.sqrt(np.abs(t[:, 0]) ** 2 + np.abs(t[:, 1]) ** 2)


def phase_features_table(
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
    pt = build_phase_tensor_table(
        S,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    if pt.empty:
        cols = [
            "station", "freq", "period", "beta_abs",
            "ellipt_abs", "logrho_det", "phi_det", "tip_amp",
        ]
        return pd.DataFrame(columns=cols)
    rows: List[Dict[str, float]] = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        T, t, _ = _get_t_block(ed)
        rho = _rho_det_from_z(z, fr)
        lgr = np.log10(np.maximum(rho, 1e-12))
        ph = _det_phase_from_z(z)
        ta = _tip_amp(t) if T is not None else None
        mask = pt["station"] == st
        pdf = pt.loc[mask]
        if pdf.empty:
            continue
        # align by nearest period
        p_pt = pdf["period"].to_numpy()
        p_z = 1.0 / fr
        idx = np.searchsorted(p_pt, p_z)
        idx = np.clip(idx, 0, len(p_pt) - 1)
        beta_abs = np.abs(pdf["beta"].to_numpy()[idx])
        ellipt = np.abs(pdf["ellipt"].to_numpy()[idx])
        for j in range(len(fr)):
            rows.append(
                dict(
                    station=st,
                    freq=float(fr[j]),
                    period=float(p_z[j]),
                    beta_abs=float(beta_abs[j]),
                    ellipt_abs=float(ellipt[j]),
                    logrho_det=float(lgr[j]),
                    phi_det=float(ph[j]),
                    tip_amp=float(ta[j]) if ta is not None else np.nan,
                )
            )
    df = pd.DataFrame.from_records(rows)
    return df

def classify_dimensionality(
    sites: Any,
    *,
    skew_th: float = 3.0,
    ellipt_th: float = 0.2,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    df = phase_features_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        return df
    lab = np.full(len(df), 2, dtype=int)
    ok2 = (df["beta_abs"] <= skew_th)
    lab[ok2 & (df["ellipt_abs"] <= ellipt_th)] = 0
    lab[ok2 & (df["ellipt_abs"] > ellipt_th)] = 1
    out = df.copy()
    out["dim"] = lab
    # 0=1D, 1=2D, 2=3D
    return out

# ---------------------- site-level masking/projection -------------------- #

def mask_by_dimensionality(
    sites: Any,
    *,
    keep: Tuple[int, ...] = (0, 1),
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    tbl = classify_dimensionality(
        S,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    if tbl.empty:
        return S

    def _one(Si):
        ed = next(_iter_items(Si))
        st = getattr(ed, "station", None) or getattr(ed, "name", None)
        Z, z, fr = _get_z_block(ed)
        T, t, ft = _get_t_block(ed)
        if Z is None:
            return Si
        sdf = tbl[tbl["station"] == st]
        if sdf.empty:
            return Si
        # map per freq via nearest period
        per = 1.0 / fr
        p_ref = sdf["period"].to_numpy()
        idx = np.searchsorted(p_ref, per)
        idx = np.clip(idx, 0, len(p_ref) - 1)
        dsel = sdf["dim"].to_numpy()[idx]
        mkeep = np.isin(dsel, keep)
        z2 = z.copy()
        z2[~mkeep] = np.nan
        Z.z = z2
        if T is not None and t is not None:
            t2 = t.copy()
            t2[~mkeep] = np.nan
            T.tipper = t2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def project_to_2d(
    sites: Any,
    *,
    strike: Optional[float] = None,
    method: str = "swift",
    antisym: bool = True,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # rotate
    if strike is None:
        S = _edit.rotate_to_strike(  # type: ignore
            S,
            method=method,
            inplace=inplace,
        )
    else:
        S = _edit.rotate(  # type: ignore
            S,
            angle=float(strike),
            inplace=inplace,
        )
    # antisymmetrize off-diagonals
    if antisym:
        from .tensor import antisymmetrize
        S = antisymmetrize(
            S,
            how="rms",
            inplace=True if inplace else False,
        )
    return S


# --------------------- dictionary learning (MOD + ISTA) ------------------ #

def _standardize(X: np.ndarray):
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0) + 1e-12
    Z = (X - mu) / sd
    Z[np.isnan(Z)] = 0.0
    return Z, mu, sd


def _soft(x: np.ndarray, t: float) -> np.ndarray:
    return np.sign(x) * np.maximum(np.abs(x) - t, 0.0)


def _ista(D: np.ndarray, x: np.ndarray, lam: float,
          n_iter: int) -> np.ndarray:
    # min 0.5||x - D a||^2 + lam||a||
    a = np.zeros(D.shape[1], dtype=float)
    smax = np.linalg.svd(D, compute_uv=False)[0]
    L = (smax ** 2) + 1e-12
    t = 1.0 / L
    for _ in range(n_iter):
        r = x - D @ a
        a = _soft(a + t * (D.T @ r), lam * t)
    return a


def _mod_update(X: np.ndarray, A: np.ndarray) -> np.ndarray:
    # D = X A^T (A A^T)^-1 ; normalize atoms
    At = A.T
    G = A @ At
    eps = 1e-8 * np.eye(G.shape[0])
    D = (X @ At) @ np.linalg.pinv(G + eps)
    # normalize columns
    n = np.linalg.norm(D, axis=0) + 1e-12
    D = D / n
    return D


def _feature_matrix(df: pd.DataFrame) -> Tuple[np.ndarray, List[str]]:
    cols = ["beta_abs", "ellipt_abs", "logrho_det", "tip_amp"]
    X = df[cols].to_numpy(dtype=float)
    X[np.isnan(X)] = 0.0
    return X, cols


def learn_dim_dictionary(
    sites: Any,
    *,
    n_atoms: int = 6,
    lam: float = 0.05,
    n_iter: int = 40,
    code_iter: int = 50,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Dict[str, Any]:
    df = phase_features_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        return dict(
            D=None, A=None, mu=None, sd=None,
            feat=[], meta=dict(samples=0),
        )
    X, feats = _feature_matrix(df)
    Z, mu, sd = _standardize(X)
    n, f = Z.shape
    k = int(max(2, min(n_atoms, max(2, f * 2))))
    rng = np.random.default_rng(1234)
    D = rng.normal(size=(f, k)).astype(float)
    # normalize initial atoms
    D = D / (np.linalg.norm(D, axis=0) + 1e-12)
    A = np.zeros((k, n), dtype=float)
    for it in range(n_iter):
        # code step
        for i in range(n):
            A[:, i] = _ista(D, Z[i], lam, code_iter)
        # dict step
        D = _mod_update(Z.T, A).T
    meta = dict(
        samples=n,
        stations=df["station"].tolist(),
        period=df["period"].to_numpy(),
        feats=feats,
    )
    return dict(D=D, A=A, mu=mu, sd=sd, feat=feats, meta=meta)


def _auto_label_atoms(D: np.ndarray, feats: List[str]) -> np.ndarray:
    # simple rule on atom means (beta vs ellipticity)
    jB = feats.index("beta_abs")
    jE = feats.index("ellipt_abs")
    b = np.abs(D[jB, :])
    e = np.abs(D[jE, :])
    lab = np.full(D.shape[1], 2, dtype=int)
    lab[(b <= 0.35) & (e <= 0.15)] = 0
    lab[(b <= 0.35) & (e > 0.15)] = 1
    return lab


def encode_dimensionality(
    sites: Any,
    model: Dict[str, Any],
    *,
    lam: float = 0.05,
    code_iter: int = 50,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    D = model.get("D", None)
    mu = model.get("mu", None)
    sd = model.get("sd", None)
    feats = model.get("feat", [])
    if D is None or mu is None or sd is None:
        return pd.DataFrame()
    df = phase_features_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        return df
    X, _ = _feature_matrix(df)
    Z = (X - mu) / (sd + 1e-12)
    Z[np.isnan(Z)] = 0.0
    k = D.shape[1]
    codes = np.zeros((len(df), k), dtype=float)
    for i in range(len(df)):
        codes[i] = _ista(D, Z[i], lam, code_iter)
    atoms_lab = _auto_label_atoms(D, feats)
    pred = atoms_lab[np.argmax(np.abs(codes), axis=1)]
    out = df.copy()
    for j in range(k):
        out[f"a{j}"] = codes[:, j]
    out["dim_pred"] = pred
    return out


def mask_by_dictionary(
    sites: Any,
    model: Dict[str, Any],
    *,
    keep: Tuple[int, ...] = (0, 1),
    lam: float = 0.05,
    code_iter: int = 50,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    tbl = encode_dimensionality(
        S,
        model,
        lam=lam,
        code_iter=code_iter,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    if tbl.empty:
        return S

    def _one(Si):
        ed = next(_iter_items(Si))
        st = getattr(ed, "station", None) or getattr(ed, "name", None)
        Z, z, fr = _get_z_block(ed)
        T, t, ft = _get_t_block(ed)
        if Z is None:
            return Si
        sdf = tbl[tbl["station"] == st]
        if sdf.empty:
            return Si
        per = 1.0 / fr
        p_ref = sdf["period"].to_numpy()
        idx = np.searchsorted(p_ref, per)
        idx = np.clip(idx, 0, len(p_ref) - 1)
        dsel = sdf["dim_pred"].to_numpy()[idx]
        mkeep = np.isin(dsel, keep)
        z2 = z.copy()
        z2[~mkeep] = np.nan
        Z.z = z2
        if T is not None and t is not None:
            t2 = t.copy()
            t2[~mkeep] = np.nan
            T.tipper = t2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)

# ---- novelty dimensionality plots (max 3) ------------------------------- #


def _dim_table_with_conf(
    sites: Any,
    *,
    skew_th: float,
    ellipt_th: float,
    recursive: bool,
    on_dup: str,
    strict: bool,
    verbose: int,
) -> pd.DataFrame:
    df = classify_dimensionality(
        sites,
        skew_th=skew_th,
        ellipt_th=ellipt_th,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        return df
    b = df["beta_abs"].to_numpy(dtype=float)
    e = df["ellipt_abs"].to_numpy(dtype=float)
    d = df["dim"].to_numpy(dtype=int)
    c = np.zeros_like(b, dtype=float)
    # raw margins (bigger → more confident)
    m0 = np.minimum(skew_th - b, ellipt_th - e)
    m1 = np.minimum(skew_th - b, e - ellipt_th)
    m2 = b - skew_th
    c[d == 0] = np.maximum(m0[d == 0], 0.0)
    c[d == 1] = np.maximum(m1[d == 1], 0.0)
    c[d == 2] = np.maximum(m2[d == 2], 0.0)
    # robust 0..1 normalization
    v = c[np.isfinite(c)]
    v0 = np.nanpercentile(v, 5) if v.size else 0.0
    v1 = np.nanpercentile(v, 95) if v.size else 1.0
    cn = (c - v0) / (v1 - v0 + 1e-12)
    cn = np.clip(cn, 0.0, 1.0)
    out = df.copy()
    out["conf"] = cn
    return out


def plot_atom_psection(
    sites: Any,
    model: dict[str, Any],
    *,
    energy: str = "l2",  # l2|l1|max
    figsize: tuple[float, float] = (9.0, 4.8),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    df = encode_dimensionality(
        S,
        model,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no codes", ha="center", va="center")
        return ax
    # gather code columns
    acols = [c for c in df.columns if c.startswith("a")]
    if not acols:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no atom cols", ha="center", va="center")
        return ax
    A = df[acols].to_numpy(dtype=float)
    P = df["period"].to_numpy(dtype=float)
    st = df["station"].astype(str).to_numpy()
    # dominant atom id and energy scalar
    Am = np.abs(A)
    aid = np.argmax(Am, axis=1)
    if energy == "l1":
        eng = np.sum(Am, axis=1)
    elif energy == "max":
        eng = np.max(Am, axis=1)
    else:
        eng = np.sqrt(np.sum(Am ** 2, axis=1))
    v = eng[np.isfinite(eng)]
    v0 = np.nanpercentile(v, 5) if v.size else 0.0
    v1 = np.nanpercentile(v, 95) if v.size else 1.0
    alp = np.clip((eng - v0) / (v1 - v0 + 1e-12), 0.0, 1.0)
    # station index and shared log-period grid
    sts = list(pd.unique(st))
    sidx = {s: i for i, s in enumerate(sts)}
    x = np.array([sidx[s] for s in st], dtype=int)
    lp = np.log10(np.maximum(P, 1e-9))
    ygrid = np.unique(lp)
    ny, nx = ygrid.size, len(sts)
    # color palette per atom
    k = A.shape[1]
    base = plt.get_cmap("tab20").colors
    rep = int(np.ceil(k / len(base)))
    cols = (list(base) * rep)[:k]
    # build RGBA image
    img = np.zeros((ny, nx, 4), dtype=float)
    for xi, yi, ai, al in zip(x, lp, aid, alp):
        yi = int(np.clip(np.searchsorted(ygrid, yi), 0, ny - 1))
        r, g, b = cols[int(ai)]
        # keep strongest alpha per cell
        if al >= img[yi, xi, 3]:
            img[yi, xi, :] = (r, g, b, al)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.imshow(
        img,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(nx))
    ax.set_xticklabels(sts, rotation=90)
    yt = np.linspace(0, ny - 1, num=min(8, ny))
    yv = np.linspace(ygrid.min(), ygrid.max(), num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    # legend: colored squares for a few atoms
    from matplotlib.lines import Line2D
    max_lab = min(k, 10)
    h = [
        Line2D([0], [0], marker="s", ls="", color=cols[j],
               label=f"a{j}")
        for j in range(max_lab)
    ]
    ax.legend(handles=h, ncols=min(max_lab, 10),
              fontsize=7, loc="upper right")
    return ax


def plot_dim_confidence_grid(
    sites: Any,
    *,
    skew_th: float = 3.0,
    ellipt_th: float = 0.2,
    figsize: tuple[float, float] = (8.8, 4.2),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    df = _dim_table_with_conf(
        S,
        skew_th=skew_th,
        ellipt_th=ellipt_th,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    df = df.copy()
    df["logp"] = np.log10(df["period"].to_numpy())
    piv_d = df.pivot_table(
        index="logp",
        columns="station",
        values="dim",
        aggfunc="median",
    ).sort_index()
    piv_c = df.pivot_table(
        index="logp",
        columns="station",
        values="conf",
        aggfunc="median",
    ).reindex(index=piv_d.index, columns=piv_d.columns)
    # piv_d / piv_c have shape (n_logp, n_stations); no transpose so that
    # imshow maps x → stations, y → log-periods.
    D = piv_d.to_numpy(dtype=float)
    C = piv_c.to_numpy(dtype=float)
    # map classes to colors; alpha from confidence
    pal = {
        0: (0.20, 0.60, 0.20),
        1: (0.20, 0.30, 0.85),
        2: (0.85, 0.25, 0.20),
    }
    rgb = np.zeros((D.shape[0], D.shape[1], 3))
    a = np.zeros((D.shape[0], D.shape[1]))
    for k, col in pal.items():
        m = (np.round(D) == k)
        for j in range(3):
            rgb[:, :, j][m] = col[j]
    a[:, :] = np.nan_to_num(C, nan=0.0)
    rgba = np.dstack([rgb, a])
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.imshow(
        rgba,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(D.shape[1]))          # shape[1] = n_stations
    ax.set_xticklabels(piv_d.columns, rotation=90)
    yt = np.linspace(0, D.shape[0] - 1, num=min(8, D.shape[0]))  # shape[0] = n_logp
    yv = np.linspace(piv_d.index.min(), piv_d.index.max(),
                     num=min(8, len(piv_d.index)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    # legend
    handles = [
        Line2D([0], [0], marker="s", ls="",
               color=pal[0], label="1D"),
        Line2D([0], [0], marker="s", ls="",
               color=pal[1], label="2D"),
        Line2D([0], [0], marker="s", ls="",
               color=pal[2], label="3D"),
    ]
    ax.legend(handles=handles, ncols=3, fontsize=8, loc="upper right")
    return ax


def plot_dim_occupancy_area(
    sites: Any,
    *,
    skew_th: float = 3.0,
    ellipt_th: float = 0.2,
    n_bands: int = 24,
    figsize: tuple[float, float] = (8.8, 3.6),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    df = _dim_table_with_conf(
        sites,
        skew_th=skew_th,
        ellipt_th=ellipt_th,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    p = df["period"].to_numpy(dtype=float)
    lo = max(float(np.nanmin(p)), 1e-6)
    hi = float(np.nanmax(p))
    edges = np.logspace(np.log10(lo), np.log10(hi), n_bands + 1)
    centers = np.sqrt(edges[:-1] * edges[1:])
    frac = np.zeros((3, n_bands), dtype=float)
    for i in range(n_bands):
        m = (p >= edges[i]) & (p < edges[i + 1])
        if not np.any(m):
            continue
        c = np.bincount(
            df.loc[m, "dim"].to_numpy(dtype=int),
            minlength=3,
        ).astype(float)
        frac[:, i] = c / (np.sum(c) + 1e-12)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    x = centers
    ax.set_xscale("log")
    y0 = frac[0]
    y1 = y0 + frac[1]
    y2 = y1 + frac[2]
    ax.fill_between(x, 0.0, y0, color=(0.20, 0.60, 0.20))
    ax.fill_between(x, y0, y1, color=(0.20, 0.30, 0.85))
    ax.fill_between(x, y1, y2, color=(0.85, 0.25, 0.20))
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("fraction")
    ax.grid(True, alpha=0.2, which="both")
    return ax


def plot_dim_map(
    sites: Any,
    *,
    period: float = 10.0,
    skew_th: float = 3.0,
    ellipt_th: float = 0.2,
    figsize: tuple[float, float] = (8.0, 6.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    df = _dim_table_with_conf(
        S,
        skew_th=skew_th,
        ellipt_th=ellipt_th,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
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
        ax.text(0.5, 0.5, "no coords", ha="center", va="center")
        return ax
    # nearest row to target period per station
    rows = []
    for st, sdf in df.groupby("station"):
        if st not in coords:
            continue
        p = sdf["period"].to_numpy(dtype=float)
        i = int(np.nanargmin(np.abs(p - period)))
        rows.append(sdf.iloc[i])
    if not rows:
        ax.text(0.5, 0.5, "no match @period",
                ha="center", va="center")
        return ax
    pal = {
        0: ((0.20, 0.60, 0.20), "o"),
        1: ((0.20, 0.30, 0.85), "s"),
        2: ((0.85, 0.25, 0.20), "^"),
    }
    for row in rows:
        st = row["station"]
        lat, lon = coords[st]
        d = int(row["dim"])
        c, m = pal.get(d, ((0.5, 0.5, 0.5), "o"))
        s = 30.0 + 120.0 * float(row["conf"])
        ax.scatter([lon], [lat], s=s, c=[c], marker=m)
    ax.set_xlabel("Lon")
    ax.set_ylabel("Lat")
    # legend
    handles = [
        Line2D([0], [0], marker="o", ls="", color=pal[0][0],
               label="1D"),
        Line2D([0], [0], marker="s", ls="", color=pal[1][0],
               label="2D"),
        Line2D([0], [0], marker="^", ls="", color=pal[2][0],
               label="3D"),
    ]
    ax.legend(handles=handles, ncols=3, fontsize=8, loc="best")
    return ax
