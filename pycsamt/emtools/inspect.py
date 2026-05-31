from __future__ import annotations 


# pycsamt/emtools/inspect.py
from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ._core import ensure_sites


# ------------------------- small helpers --------------------------------- #

def _is_df(x: Any) -> bool:
    return isinstance(x, pd.DataFrame)


def _name(ed: Any, idx: int) -> str:
    for k in ("station", "name", "site", "id"):
        v = getattr(ed, k, None)
        if isinstance(v, str) and v:
            return v
    return f"site_{idx}"


def _coords(ed: Any) -> Tuple[Optional[float], Optional[float]]:
    lat = getattr(ed, "lat", None)
    lon = getattr(ed, "lon", None)
    if lat is None:
        lat = getattr(ed, "latitude", None)
    if lon is None:
        lon = getattr(ed, "longitude", None)
    try:
        lat = float(lat) if lat is not None else None
        lon = float(lon) if lon is not None else None
    except Exception:
        lat, lon = None, None
    return lat, lon


def _has(ed: Any, sect: str) -> bool:
    f = getattr(ed, "has_section", None)
    if callable(f):
        try:
            return bool(f(sect))
        except Exception:
            pass
    # duck-typing fallback
    return hasattr(ed, sect) or hasattr(ed, sect.capitalize())


def _get_freq(ed: Any) -> Optional[np.ndarray]:
    # try Z first
    for a in ("Z", "z", "ResPhase", "resphase", "Tipper", "tipper"):
        obj = getattr(ed, a, None)
        if obj is None:
            continue
        fr = getattr(obj, "freq", None)
        if fr is not None:
            try:
                ar = np.asarray(fr, dtype=float)
                if ar.ndim == 1 and ar.size > 0:
                    return ar
            except Exception:
                pass
    # last resort
    fr = getattr(ed, "freq", None)
    if fr is not None:
        try:
            ar = np.asarray(fr, dtype=float)
            if ar.ndim == 1 and ar.size > 0:
                return ar
        except Exception:
            pass
    return None


def _iter_items(sites: Any) -> Iterable[Any]:
    # Sites is iterable in most builds
    try:
        for it in sites:
            yield it
        return
    except Exception:
        pass
    # fallback: dict-like
    items = getattr(sites, "items", None)
    if isinstance(items, dict):
        for _, it in items.items():
            yield it


def _period(fr: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        p = 1.0 / fr
    return p


def _union_freq(freqs: List[np.ndarray]) -> np.ndarray:
    if not freqs:
        return np.array([], dtype=float)
    allf = np.unique(np.concatenate(freqs))
    allf = allf[np.isfinite(allf) & (allf > 0.0)]
    return allf


def _presence_vec(fr: np.ndarray, grid: np.ndarray) -> np.ndarray:
    if fr.size == 0 or grid.size == 0:
        return np.zeros(grid.size, dtype=bool)
    idx = np.searchsorted(grid, fr)
    idx = np.clip(idx, 0, grid.size - 1)
    # mark nearest only when close
    tol = 1e-10
    near = np.abs(grid[idx] - fr) <= tol
    m = np.zeros(grid.size, dtype=bool)
    m[idx[near]] = True
    return m


def _df_from_kv(
    rows: List[Dict[str, Any]],
    cols: List[str],
) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame.from_records(rows, columns=cols)


# --------------------------- public API ---------------------------------- #

def sites_summary(
    sites: Any,
    *,
    fields: Tuple[str, ...] = (
        "station",
        "n_freq",
        "has_tipper",
        "period_min",
        "period_max",
        "lat",
        "lon",
    ),
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
    rows: List[Dict[str, Any]] = []
    for i, ed in enumerate(_iter_items(S)):
        fr = _get_freq(ed)
        nm = _name(ed, i)
        lat, lon = _coords(ed)
        has_tip = _has(ed, "tipper")
        nf = int(fr.size) if isinstance(fr, np.ndarray) else 0
        pmin = None
        pmax = None
        if isinstance(fr, np.ndarray) and fr.size > 0:
            pe = _period(fr)
            pmin = float(np.nanmin(pe))
            pmax = float(np.nanmax(pe))
        rec = dict(
            station=nm,
            n_freq=nf,
            has_tipper=bool(has_tip),
            period_min=pmin,
            period_max=pmax,
            lat=lat,
            lon=lon,
        )
        rows.append(rec)
    df = _df_from_kv(rows, list(fields))
    return df[list(fields)]


def list_missing_sections(
    sites: Any,
    *,
    require: Tuple[str, ...] = ("mt", "tipper"),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Dict[str, List[str]]:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    out: Dict[str, List[str]] = {}
    for i, ed in enumerate(_iter_items(S)):
        nm = _name(ed, i)
        miss: List[str] = []
        for sec in require:
            if not _has(ed, sec):
                miss.append(sec)
        if miss:
            out[nm] = miss
    return out


def frequency_coverage(
    sites: Any,
    *,
    mode: str = "per-site",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    freqs: Dict[str, np.ndarray] = {}
    for i, ed in enumerate(_iter_items(S)):
        nm = _name(ed, i)
        fr = _get_freq(ed)
        if isinstance(fr, np.ndarray):
            freqs[nm] = fr
    if mode == "per-site":
        return freqs
    grid = _union_freq(list(freqs.values()))
    if mode == "union":
        return grid
    if mode == "intersection":
        if not freqs:
            return np.array([], dtype=float)
        # numeric set-intersection with tol
        g = None
        for fr in freqs.values():
            if g is None:
                g = fr
            else:
                a = np.intersect1d(g, fr)
                g = a
        return g if g is not None else np.array([], dtype=float)
    raise ValueError("mode must be per-site|union|intersection")


def plot_coverage(
    sites: Any,
    *,
    axis: str = "period",
    show_mask: bool = True,
    figsize: Tuple[float, float] = (7.0, 4.0),
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
    data = frequency_coverage(
        S,
        mode="per-site",
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    labs = list(data.keys())
    arr = list(data.values())
    grid = _union_freq(arr)
    if grid.size == 0:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    M = np.vstack([_presence_vec(fr, grid) for fr in arr])
    y = _period(grid) if axis == "period" else grid
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        M.T,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("site")
    ax.set_ylabel(axis)
    ax.set_xticks(np.arange(len(labs)))
    ax.set_xticklabels(labs, rotation=90)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("presence")
    return ax

# ------------------------ impedance curves -------------------------------- #

def _df_resphase(
    sites: Any,
    *,
    kind: str = "resphase",
    recursive: bool = False,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Optional[pd.DataFrame]:
    tdf = None
    get_df = getattr(sites, "to_dataframe", None)
    if callable(get_df):
        try:
            tdf = get_df(kind=kind)
        except TypeError:
            try:
                tdf = get_df(kind)
            except Exception:
                tdf = None
    return tdf if _is_df(tdf) else None


def plot_rhoa_phi(
    sites: Any,
    *,
    components: Tuple[str, ...] = ("xy", "yx"),
    axis: str = "period",
    errorbar: bool = True,
    figsize: Tuple[float, float] = (7.5, 6.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax_r: Optional[plt.Axes] = None,
    ax_p: Optional[plt.Axes] = None,
) -> Tuple[plt.Axes, plt.Axes]:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    df = _df_resphase(S, kind="resphase")
    if df is None or df.empty:
        # Fallback: try z dataframe and convert later if needed
        df = _df_resphase(S, kind="z")
    if df is None or df.empty:
        if ax_r is None or ax_p is None:
            fig, (ax_r, ax_p) = plt.subplots(
                2, 1, figsize=figsize, sharex=True
            )
        ax_r.text(0.5, 0.5, "no data", ha="center", va="center")
        ax_p.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax_r, ax_p
    if ax_r is None or ax_p is None:
        fig, (ax_r, ax_p) = plt.subplots(
            2, 1, figsize=figsize, sharex=True
        )
    # Expected columns are package-defined; be defensive
    # We accept generic patterns: station, freq, period,
    # rho_xy, rho_yx, phi_xy, phi_yx, and maybe errors.
    xkey = "period" if axis == "period" else "freq"
    if xkey not in df.columns and "freq" in df.columns:
        df["period"] = 1.0 / df["freq"].replace(0, np.nan)
        xkey = "period" if axis == "period" else "freq"
    stations = sorted(df.get("station", pd.Series()).unique())
    if not stations:
        # try without station split
        stations = ["__all__"]
        df = df.copy()
        df["station"] = stations[0]
    for st in stations:
        sdf = df[df["station"] == st]
        for comp in components:
            rk = f"rho_{comp}"
            pk = f"phi_{comp}"
            rke = f"{rk}_err"
            pke = f"{pk}_err"
            if rk in sdf.columns and pk in sdf.columns:
                x = sdf[xkey].to_numpy()
                r = sdf[rk].to_numpy()
                p = sdf[pk].to_numpy()
                if errorbar and rke in sdf.columns:
                    re = sdf[rke].to_numpy()
                    ax_r.errorbar(x, r, yerr=re, label=f"{st}:{comp}")
                else:
                    ax_r.plot(x, r, ".", label=f"{st}:{comp}")
                if errorbar and pke in sdf.columns:
                    pe = sdf[pke].to_numpy()
                    ax_p.errorbar(x, p, yerr=pe, label=f"{st}:{comp}")
                else:
                    ax_p.plot(x, p, ".", label=f"{st}:{comp}")
    ax_r.set_xscale("log")
    ax_r.set_yscale("log")
    ax_r.set_ylabel("rho_a")
    ax_p.set_xscale("log")
    ax_p.set_ylabel("phi (deg)")
    ax_p.set_xlabel(axis)
    ax_r.legend(ncols=2, fontsize=8)
    ax_p.legend(ncols=2, fontsize=8)
    return ax_r, ax_p


# ----------------------------- tipper ------------------------------------- #

def plot_tipper_components(
    sites: Any,
    *,
    kind: Tuple[str, ...] = ("real", "imag"),
    axis: str = "period",
    figsize: Tuple[float, float] = (7.5, 4.5),
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
    df = _df_resphase(S, kind="tipper")
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if df is None or df.empty:
        ax.text(0.5, 0.5, "no tipper", ha="center", va="center")
        return ax
    xkey = "period" if axis == "period" else "freq"
    if xkey not in df.columns and "freq" in df.columns:
        df["period"] = 1.0 / df["freq"].replace(0, np.nan)
        xkey = "period" if axis == "period" else "freq"
    stations = sorted(df.get("station", pd.Series()).unique())
    if not stations:
        stations = ["__all__"]
        df = df.copy()
        df["station"] = stations[0]
    comps = ("tx", "ty")
    for st in stations:
        sdf = df[df["station"] == st]
        for c in comps:
            for k in kind:
                col = f"{c}_{k}"
                if col in sdf.columns:
                    ax.plot(
                        sdf[xkey].to_numpy(),
                        sdf[col].to_numpy(),
                        ".",
                        label=f"{st}:{col}",
                    )
    ax.set_xscale("log")
    ax.set_xlabel(axis)
    ax.set_ylabel("tipper")
    ax.legend(ncols=2, fontsize=8)
    return ax


# --------------------------- pseudosection -------------------------------- #

def pseudosection(
    sites: Any,
    *,
    quantity: str = "rho_xy",
    axis_x: str = "station",
    axis_y: str = "period",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize: Tuple[float, float] = (7.5, 4.5),
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
    df = _df_resphase(S, kind="resphase")
    if df is None or df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no resphase", ha="center", va="center")
        return ax
    # ensure axes keys
    if axis_y not in df.columns:
        if "freq" in df.columns and axis_y == "period":
            df = df.copy()
            df["period"] = 1.0 / df["freq"].replace(0, np.nan)
        else:
            raise ValueError("axis_y not available")
    if axis_x not in df.columns:
        raise ValueError("axis_x not available")
    # pivot
    piv = df.pivot_table(
        index=axis_y,
        columns=axis_x,
        values=quantity,
        aggfunc="median",
    )
    piv = piv.sort_index()
    # X = np.arange(piv.shape[1])
    Y = piv.index.to_numpy()
    Z = piv.to_numpy(dtype=float)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        Z.T,
        aspect="auto",
        origin="lower",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_xlabel(axis_x)
    ax.set_ylabel(axis_y)
    ax.set_yticks(
        np.linspace(0, len(Y) - 1, num=min(8, len(Y)))
    )
    yt = np.linspace(0, len(Y) - 1, num=min(8, len(Y)))
    ax.set_yticklabels([f"{Y[int(i)]:.3g}" for i in yt])
    ax.set_xticks(np.arange(len(piv.columns)))
    ax.set_xticklabels(piv.columns, rotation=90)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(quantity)
    return ax
