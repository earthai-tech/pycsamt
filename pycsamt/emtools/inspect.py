from __future__ import annotations 


# pycsamt/emtools/inspect.py
from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from ..api.style import PYCSAMT_STYLE
from ..api.view import maybe_wrap_frame
from ._core import ensure_sites, _axes_list, _get_z_block, _get_t_block

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
    api: bool | None = None,
) -> Any:
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
    df = df[list(fields)]

    return maybe_wrap_frame(
        df,
        api=api,
        name="sites_summary",
        kind="edi.summary",
        source=sites,
        description=(
            "Per-site EDI frequency, tipper, and coordinate summary."
        ),
    )


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
    cb = ax.get_figure().colorbar(im, ax=ax)
    cb.set_label("presence")
    return ax

# ------------------------ impedance curves -------------------------------- #

_RESPHASE_RENAME = {
    # Site.to_dataframe() uses z-prefix (rho_zxy); plot functions expect rho_xy
    "rho_zxx": "rho_xx", "phi_zxx": "phi_xx",
    "rho_zxy": "rho_xy", "phi_zxy": "phi_xy",
    "rho_zyx": "rho_yx", "phi_zyx": "phi_yx",
    "rho_zyy": "rho_yy", "phi_zyy": "phi_yy",
}


def _df_resphase(
    sites: Any,
    *,
    kind: str = "resphase",
    recursive: bool = False,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Optional[pd.DataFrame]:
    # Try the object's own to_dataframe first (works for a single Site)
    get_df = getattr(sites, "to_dataframe", None)
    if callable(get_df):
        try:
            tdf = get_df(kind=kind)
        except TypeError:
            try:
                tdf = get_df(kind)
            except Exception:
                tdf = None
        else:
            pass
        normed = _normalise_resphase_df(tdf, station=_name(sites, 0))
        if normed is not None:
            return normed

    # Fallback: sites is a collection (Sites is iterable, yields Site objects)
    dfs: List[pd.DataFrame] = []
    try:
        for i, item in enumerate(_iter_items(sites)):
            item_get = getattr(item, "to_dataframe", None)
            if not callable(item_get):
                continue
            try:
                idf = item_get(kind=kind)
            except Exception:
                continue
            normed = _normalise_resphase_df(idf, station=_name(item, i))
            if normed is not None:
                dfs.append(normed)
    except Exception:
        pass

    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return None


def _unwrap_frame(obj: Any) -> Optional[pd.DataFrame]:
    """Return a plain pd.DataFrame from obj, unwrapping APIFrame if needed."""
    if isinstance(obj, pd.DataFrame):
        return obj
    # APIFrame stores the real DataFrame in .df
    inner = getattr(obj, "df", None)
    if isinstance(inner, pd.DataFrame):
        return inner
    return None


def _normalise_resphase_df(raw: Any, station: str) -> Optional[pd.DataFrame]:
    """Add station/freq columns and normalise column names from to_dataframe output."""
    df = _unwrap_frame(raw)
    if df is None:
        return None
    df = df.copy()
    # Flatten freq index produced by Site.to_dataframe (index.name == 'f')
    if df.index.name in ("f", "freq") and "freq" not in df.columns:
        df = df.reset_index()
        df = df.rename(columns={"f": "freq"})
    if "station" not in df.columns:
        df["station"] = station
    # Rename rho_zxy → rho_xy etc.
    df = df.rename(columns=_RESPHASE_RENAME)
    return df


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
    from ._core import _iter_items, _name, _get_t_block

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    _STYLES = {
        ("tx", "real"): dict(color="#1f77b4", ls="-",  marker="o", ms=3, lw=0.9),
        ("tx", "imag"): dict(color="#1f77b4", ls="--", marker="o", ms=3, lw=0.9),
        ("ty", "real"): dict(color="#d62728", ls="-",  marker="s", ms=3, lw=0.9),
        ("ty", "imag"): dict(color="#d62728", ls="--", marker="s", ms=3, lw=0.9),
    }

    n_plotted = 0
    for i, ed in enumerate(_iter_items(S)):
        nm  = _name(ed, i)
        out = _get_t_block(ed)
        _, t, fr = out[0], out[1], out[2]
        if t is None or fr is None:
            continue
        t  = np.asarray(t, complex)   # (n_freq, 2) — Tx, Ty
        fr = np.asarray(fr, float)
        if t.ndim == 3 and t.shape[1] == 1 and t.shape[2] == 2:
            t = t[:, 0, :]
        if t.ndim != 2 or t.shape[1] < 2:
            continue

        x  = (1.0 / np.where(fr == 0, np.nan, fr)
               if axis == "period" else fr)

        for ci, comp in enumerate(("tx", "ty")):
            col_t = t[:, ci]
            for k in kind:
                vals = col_t.real if k == "real" else col_t.imag
                sty  = _STYLES.get((comp, k), {})
                ax.plot(x, vals, label=f"{nm}:{comp}_{k}", **sty)
                n_plotted += 1

    ax.set_xscale("log")
    ax.set_xlabel("Period (s)" if axis == "period" else "Frequency (Hz)", fontsize=9)
    ax.set_ylabel("Tipper", fontsize=9)
    ax.axhline(0, color="k", lw=0.5, ls=":")
    ax.tick_params(labelsize=8)
    if n_plotted:
        ax.legend(ncols=2, fontsize=7)
    else:
        ax.text(0.5, 0.5, "no tipper data", ha="center", va="center",
                transform=ax.transAxes, color="0.5")
    return ax


# --------------------------- pseudosection -------------------------------- #

def pseudosection(
    sites: Any,
    *,
    quantity: str = "rho_xy",
    axis_x: str = "station",
    axis_y: str = "period",
    period_range: Optional[Tuple[float, float]] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize: Tuple[float, float] = (7.5, 4.5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
    topo: Optional[bool] = None,
    dark: bool = True,
) -> plt.Axes:
    """Draw a period-vs-station pseudosection.

    Parameters
    ----------
    sites : Sites or compatible
        Station data.
    quantity : str
        Column name to plot (e.g. ``"rho_xy"``, ``"phi_xy"``).
    topo : bool or None
        Override the global :data:`~pycsamt.topo.PYCSAMT_TOPO` setting.
        ``None`` (default) reads the global singleton.
    dark : bool
        Use dark-palette styling for the topo strip.
    """
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
    # optional period-range filter (T in seconds; applied before pivot)
    if period_range is not None and axis_y == "period":
        T_min, T_max = float(period_range[0]), float(period_range[1])
        df = df[(df["period"] >= T_min) & (df["period"] <= T_max)]
        if df.empty:
            if ax is None:
                _, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "No data in selected period range",
                    ha="center", va="center", fontsize=9)
            return ax
    # pivot
    piv = df.pivot_table(
        index=axis_y,
        columns=axis_x,
        values=quantity,
        aggfunc="median",
    )
    # Sort periods ascending; shallow (short period) will go to the top.
    piv = piv.sort_index()
    Y = piv.index.to_numpy()          # period values, ascending
    Z = piv.to_numpy(dtype=float)     # shape (n_period, n_station)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    # Use Z directly (NOT Z.T) so that:
    #   x-axis = station index  (cols 0..n_station-1)  ← matches set_xticks
    #   y-axis = period index   (rows 0..n_period-1)
    # origin="upper" → row 0 (smallest period, shallowest) at the TOP,
    # which is the standard MT pseudosection convention.
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="upper",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )

    # ── X-axis: station labels at the TOP (MT convention) ────────────
    station_labels = list(piv.columns)
    ax.set_xticks(np.arange(len(station_labels)))
    ax.set_xticklabels(station_labels, rotation=90, fontsize=7, ha="center")
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel(axis_x, fontsize=9, labelpad=4)

    # ── Y-axis: period values ─────────────────────────────────────────
    n_yticks = min(8, len(Y))
    yt = np.linspace(0, len(Y) - 1, num=n_yticks)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{Y[int(i)]:.3g}" for i in yt], fontsize=8)
    ax.set_ylabel(axis_y, fontsize=9)

    # ── Colorbar ──────────────────────────────────────────────────────
    cb = ax.get_figure().colorbar(im, ax=ax)
    cb.set_label(quantity, fontsize=8)

    # ── Topography strip (optional) ───────────────────────────────────
    _topo_enabled = topo if topo is not None else _topo_cfg().enabled
    if _topo_enabled:
        try:
            from pycsamt.topo.extract import (
                extract_elevation, extract_chainage, extract_station_names,
            )
            from pycsamt.topo.overlay import draw_topo_strip
            chain_km = extract_chainage(S)
            elev_m   = extract_elevation(S)
            names    = extract_station_names(S)
            # Reorder to match pivot column order (station labels)
            name_to_idx = {n: i for i, n in enumerate(names)}
            ordered_idx = [name_to_idx.get(str(s), i)
                           for i, s in enumerate(station_labels)]
            chain_ord = chain_km[ordered_idx] if len(chain_km) else chain_km
            elev_ord  = elev_m[ordered_idx]   if len(elev_m)   else elev_m
            strip_names = [str(s) for s in station_labels]
            draw_topo_strip(
                ax.get_figure(), ax,
                chain_ord, elev_ord, strip_names,
                dark=dark,
            )
        except Exception:
            pass   # topo strip is optional; never break the main plot

    return ax


def _topo_cfg():
    """Return global TopoConfig singleton (lazy import to avoid circular deps)."""
    from pycsamt.topo.config import PYCSAMT_TOPO
    return PYCSAMT_TOPO


# ──────────────────────────────────────────────────────────────────────────────
# plot_station_response — full 4-component impedance tensor + tipper view
# ──────────────────────────────────────────────────────────────────────────────

_COMP_IDX: dict = {"xx": (0, 0), "xy": (0, 1), "yx": (1, 0), "yy": (1, 1)}
_TIP_PANELS: list = [
    ("Re($T_x$)",  "tx", "real"),
    ("Im($T_x$)",  "tx", "imag"),
    ("Re($T_y$)",  "ty", "real"),
    ("Im($T_y$)",  "ty", "imag"),
]


def _rho_phase_err(rho, z, z_err):
    """Return (rho_err, phase_err_deg) from impedance absolute errors."""
    z_abs = np.abs(z)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(z_abs > 0, np.abs(z_err) / z_abs, 0.0)
    rho_err    = 2.0 * rho * rel
    phase_err  = np.degrees(rel)
    return rho_err, phase_err


def _rms_log(obs, mod):
    """Root-mean-square misfit in log10(rho) space."""
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.log10(obs / np.where(mod > 0, mod, np.nan))
    r = r[np.isfinite(r)]
    return float(np.sqrt(np.mean(r ** 2))) if r.size else np.nan


def _pick_station(S, station):
    """Return the first (or named) Site object from S."""
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        if station is None or st == station:
            return st, ed
    raise RuntimeError(
        f"Station {station!r} not found." if station
        else "No sites in collection."
    )


def plot_station_response(
    sites: Any,
    *,
    station: Optional[str] = None,
    sites_model: Optional[Any] = None,
    components: Tuple[str, ...] = ("xx", "xy", "yx", "yy"),
    period_range: Optional[Tuple[float, float]] = None,
    rho_lim: Optional[Tuple[float, float]] = None,
    phase_lim: Optional[Tuple[float, float]] = None,
    tipper_lim: Tuple[float, float] = (-0.5, 0.5),
    show_tipper: bool = True,
    show_error_bars: bool = True,
    show_rms: bool = True,
    title: str = "",
    axes=None,
    figsize: Optional[Tuple[float, float]] = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Full four-component impedance tensor + tipper response for one station.

    Renders a 3-row × N-column figure (N = len(*components*)) following the
    :class:`~pycsamt.api.style.MTComponentStyle` colour scheme:

    * **Row 0** — apparent resistivity ρa (Ω·m) vs period, log–log.
    * **Row 1** — phase φ (°) vs period, semilog-x.
    * **Row 2** — tipper magnitude: Re(Tx), Im(Tx), Re(Ty), Im(Ty)
      vs period, semilog-x.  Hidden when *show_tipper* is ``False`` or when
      no tipper data are found.

    Each column corresponds to one impedance tensor component (Zxx, Zxy,
    Zyx, Zyy).  The tipper row always uses four fixed sub-panels regardless
    of which Z components are selected.

    An optional *sites_model* argument overlays a second (model/forward)
    dataset on the same axes as dotted lines.  When both observed and model
    data are present, a per-component RMS is computed in log10(ρa) space and
    appended to each column header.

    Parameters
    ----------
    sites : any
        Observed EDI data — path, :class:`~pycsamt.core.base.SitesCollection`,
        or anything accepted by :func:`~pycsamt.emtools._core.ensure_sites`.
    station : str or None
        Name of the station to plot.  ``None`` picks the first available.
    sites_model : any or None
        Optional forward-model or inversion-response EDI data (same API as
        *sites*).  When provided, overlay dashed lines and show RMS.
    components : tuple of {"xx","xy","yx","yy"}, default all four
        Z-tensor components to display.  Order determines column order.
    period_range : (T_min, T_max) or None
        Clip period axis to this window (seconds).
    rho_lim : (vmin, vmax) or None
        ρa y-axis limits.  ``None`` → matplotlib auto.
    phase_lim : (lo, hi) or None
        Phase y-axis limits in degrees.  ``None`` → auto.
    tipper_lim : (lo, hi), default (-0.5, 0.5)
        Tipper y-axis limits.
    show_tipper : bool, default True
        Whether to add the tipper row.
    show_error_bars : bool, default True
        Draw error bars on observed data.
    show_rms : bool, default True
        Append per-component RMS to column titles when *sites_model* is set.
    title : str
        Override the auto-derived figure title.
    figsize : (float, float) or None
        Figure size.  Auto-computed when ``None``.
    recursive, on_dup, strict, verbose
        Passed to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    Observed data only:

    >>> from pycsamt.emtools import plot_station_response
    >>> fig = plot_station_response("path/to/edis/", station="S07")

    With model overlay:

    >>> fig = plot_station_response(
    ...     obs_edis, station="HBH03_IMP",
    ...     sites_model=model_edis,
    ...     period_range=(1e-4, 1.0),
    ... )
    """

    mt = PYCSAMT_STYLE.mt

    # ── 1. load observed site ─────────────────────────────────────────────
    S_obs = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                         strict=strict, verbose=verbose)
    st_name, ed_obs = _pick_station(S_obs, station)

    Z_obs, z_obs, fr_obs = _get_z_block(ed_obs)
    rho_obs   = getattr(ed_obs, "rho",   None)
    phase_obs = getattr(ed_obs, "phase", None)
    z_err_obs = getattr(Z_obs,  "z_err", None) if Z_obs else None
    T_obs, t_obs, _ = _get_t_block(ed_obs)
    tip_err_obs = getattr(T_obs, "tipper_err", None) if T_obs else None

    axes_given = _axes_list(axes, 1) if axes is not None else None
    if rho_obs is None or phase_obs is None or z_obs is None:
        if axes_given is None:
            fig, ax = plt.subplots(figsize=figsize or (10, 6))
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no impedance data", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    per_obs = 1.0 / np.where(fr_obs == 0, np.nan, fr_obs)

    # ── 2. period mask ────────────────────────────────────────────────────
    if period_range is not None:
        lo, hi = float(period_range[0]), float(period_range[1])
        mask = np.isfinite(per_obs) & (per_obs >= lo) & (per_obs <= hi)
    else:
        mask = np.isfinite(per_obs)

    per = per_obs[mask]
    rho   = rho_obs[mask]
    phase = phase_obs[mask]
    z_m   = z_obs[mask]
    z_e   = z_err_obs[mask] if z_err_obs is not None else None
    tip   = t_obs[mask] if t_obs is not None else None
    tip_e = tip_err_obs[mask] if tip_err_obs is not None else None

    rho_err, phase_err = (
        _rho_phase_err(rho, z_m, z_e) if z_e is not None
        else (np.zeros_like(rho), np.zeros_like(phase))
    )

    # ── 3. optional model site ────────────────────────────────────────────
    rho_mod = phase_mod = per_mod = tip_mod = None
    if sites_model is not None:
        S_mod = ensure_sites(sites_model, recursive=recursive, on_dup=on_dup,
                             strict=False, verbose=0)
        try:
            _, ed_mod = _pick_station(S_mod, st_name)
            Z_mod, z_mod_arr, fr_mod = _get_z_block(ed_mod)
            rho_mod_raw   = getattr(ed_mod, "rho",   None)
            phase_mod_raw = getattr(ed_mod, "phase", None)
            if rho_mod_raw is not None and fr_mod is not None:
                per_m = 1.0 / np.where(fr_mod == 0, np.nan, fr_mod)
                msk_m = np.isfinite(per_m)
                if period_range is not None:
                    msk_m &= (per_m >= lo) & (per_m <= hi)
                per_mod   = per_m[msk_m]
                rho_mod   = rho_mod_raw[msk_m]
                phase_mod = phase_mod_raw[msk_m]
            T_mod, t_mod_arr, _ = _get_t_block(ed_mod)
            if t_mod_arr is not None:
                tip_mod = t_mod_arr[msk_m] if period_range else t_mod_arr
        except Exception:
            pass

    # ── 4. figure layout ──────────────────────────────────────────────────
    n_comp = len(components)
    has_tipper = show_tipper and tip is not None

    n_rows = 3 if has_tipper else 2
    h_ratios = ([2.2, 1.6, 1.0] if has_tipper else [2.2, 1.6])

    if figsize is None:
        figsize = (3.4 * max(n_comp, 4) + 0.6, 4.5 + (1.8 if has_tipper else 0))

    n_required_axes = (2 * n_comp) + (4 if has_tipper else 0)
    axes_given = _axes_list(axes, n_required_axes) if axes is not None else None
    if axes_given is None:
        fig = plt.figure(figsize=figsize, constrained_layout=True)
        gs = gridspec.GridSpec(
            n_rows, max(n_comp, 4) if has_tipper else n_comp,
            figure=fig,
            height_ratios=h_ratios,
            hspace=0.08,
            wspace=0.28,
        )

        # Build axes: [rho row, phase row] share x across columns
        ax_rho   = [fig.add_subplot(gs[0, c]) for c in range(n_comp)]
        ax_phase = [fig.add_subplot(gs[1, c], sharex=ax_rho[c]) for c in range(n_comp)]
        if has_tipper:
            ax_tip = [fig.add_subplot(gs[2, c]) for c in range(4)]
        else:
            ax_tip = []
    else:
        fig = axes_given[0].figure
        ax_rho = axes_given[:n_comp]
        ax_phase = axes_given[n_comp:2 * n_comp]
        ax_tip = axes_given[2 * n_comp:2 * n_comp + 4] if has_tipper else []

    # ── 5. per-component RMS store ────────────────────────────────────────
    rms_dict: dict = {}

    # ── 6. draw ρa and φ panels ───────────────────────────────────────────
    for ci, comp in enumerate(components):
        if comp.lower() not in _COMP_IDX:
            continue
        ri, ci_mat = _COMP_IDX[comp.lower()]
        cs = mt.component(comp.lower())
        ekw = cs.errorbar_kwargs()

        ax_r = ax_rho[ci]
        ax_p = ax_phase[ci]

        # ── observed ────────────────────────────────────────────────
        rho_c   = rho[:, ri, ci_mat]
        phi_c   = phase[:, ri, ci_mat]
        rerr_c  = rho_err[:, ri, ci_mat]   if show_error_bars else None
        perr_c  = phase_err[:, ri, ci_mat] if show_error_bars else None

        # mask non-positive ρa
        valid = rho_c > 0
        xp = per[valid]; rr = rho_c[valid]; pp = phi_c[valid]
        if rerr_c is not None:
            re = rerr_c[valid]
            yerr_r = re if np.any(re > 0) else None
        else:
            yerr_r = None
        if perr_c is not None:
            pe = perr_c[valid]
            yerr_p = pe if np.any(pe > 0) else None
        else:
            yerr_p = None

        obs_label = f"$Z_{{\\rm {comp.upper()}}}$"
        ax_r.errorbar(xp, rr, yerr=yerr_r, **{**ekw, "label": obs_label})
        ax_p.errorbar(xp, pp, yerr=yerr_p, **{**ekw, "label": obs_label})

        # ── model overlay ────────────────────────────────────────────
        if rho_mod is not None and per_mod is not None:
            rho_mc  = rho_mod[:, ri, ci_mat]
            phi_mc  = phase_mod[:, ri, ci_mat]
            valid_m = rho_mc > 0
            # compute RMS
            rms_val = _rms_log(rr, np.interp(xp, per_mod[valid_m], rho_mc[valid_m],
                                               left=np.nan, right=np.nan))
            if np.isfinite(rms_val):
                rms_dict[comp] = rms_val
            mod_lbl = f"$Z^m_{{\\rm {comp.upper()}}}$"
            if show_rms and np.isfinite(rms_val):
                mod_lbl += f"  rms={rms_val:.2f}"
            mk = cs.plot_kwargs(ls="--", lw=1.6, alpha=0.85, marker="",
                                label=mod_lbl)
            ax_r.plot(per_mod[valid_m], rho_mc[valid_m], **mk)
            ax_p.plot(per_mod, phi_mc, **{**mk, "label": ""})

        # ── axes styling ─────────────────────────────────────────────
        ax_r.set_xscale("log")
        ax_r.set_yscale("log")
        ax_r.grid(True, which="both", alpha=0.2, lw=0.5)
        ax_r.grid(True, which="major", alpha=0.4, lw=0.7)
        if rho_lim:
            ax_r.set_ylim(*rho_lim)
        plt.setp(ax_r.get_xticklabels(), visible=False)
        ax_r.tick_params(labelsize=7)
        ax_r.legend(fontsize=7, framealpha=0.85, loc="best",
                    handlelength=1.5, borderpad=0.4)

        ax_p.set_xscale("log")
        ax_p.grid(True, which="both", alpha=0.2, lw=0.5)
        ax_p.grid(True, which="major", alpha=0.4, lw=0.7)
        if phase_lim:
            ax_p.set_ylim(*phase_lim)
        ax_p.tick_params(labelsize=7)
        ax_p.set_xlabel("Period (s)", fontsize=8)

        # column header
        rms_str = (f"  rms={rms_dict[comp]:.2f}"
                   if (show_rms and rms_dict.get(comp) is not None) else "")
        ax_r.set_title(
            f"$Z_{{\\rm {comp.upper()}}}${rms_str}",
            fontsize=9, pad=6, fontweight="bold",
        )

        # shared y labels on leftmost column only
        if ci == 0:
            ax_r.set_ylabel(r"App. Res. (Ω·m)", fontsize=8)
            ax_p.set_ylabel("Phase (deg)", fontsize=8)

    # ── 7. tipper panels ──────────────────────────────────────────────────
    if has_tipper and len(ax_tip) == 4:
        # Tx=t[:,0], Ty=t[:,1]
        tip_series = {
            "tx": tip[:, 0] if tip is not None else None,
            "ty": tip[:, 1] if tip is not None else None,
        }
        tip_err_series = {
            "tx": (tip_e[:, 0]
                   if (tip_e is not None and tip_e.shape[1] > 0) else None),
            "ty": (tip_e[:, 1]
                   if (tip_e is not None and tip_e.shape[1] > 1) else None),
        }
        tip_mod_series = {
            "tx": tip_mod[:, 0] if tip_mod is not None else None,
            "ty": tip_mod[:, 1] if tip_mod is not None else None,
        }
        tip_colors = {"tx": "#1f77b4", "ty": "#d62728"}  # blue, red

        for ti, (lbl, tc_key, part) in enumerate(_TIP_PANELS):
            ax_t = ax_tip[ti]
            tv = tip_series.get(tc_key)
            if tv is None:
                ax_t.text(0.5, 0.5, "—", ha="center", va="center",
                          transform=ax_t.transAxes, color="0.55")
                continue
            col = tip_colors[tc_key]
            vals = np.real(tv) if part == "real" else np.imag(tv)
            te_v = tip_err_series.get(tc_key)
            if te_v is not None and show_error_bars:
                err_v = np.abs(te_v)
                err_v = err_v if np.any(err_v > 0) else None
            else:
                err_v = None
            ax_t.errorbar(per, vals, yerr=err_v, fmt="o", ms=3.5, color=col,
                          elinewidth=0.8, capsize=2, alpha=0.9, zorder=3)
            # model overlay
            tm_v = tip_mod_series.get(tc_key)
            if tm_v is not None and per_mod is not None:
                mod_vals = np.real(tm_v) if part == "real" else np.imag(tm_v)
                ax_t.plot(per_mod, mod_vals, "--", color=col,
                          lw=1.5, alpha=0.85, zorder=2)
            ax_t.axhline(0.0, color="0.65", lw=0.7, ls=":", zorder=1)
            ax_t.set_xscale("log")
            ax_t.set_ylim(*tipper_lim)
            ax_t.grid(True, which="both", alpha=0.2, lw=0.5)
            ax_t.grid(True, which="major", alpha=0.4, lw=0.7)
            ax_t.tick_params(labelsize=7)
            ax_t.set_xlabel("Period (s)", fontsize=8)
            ax_t.set_title(lbl, fontsize=8, pad=4)
            if ti == 0:
                ax_t.set_ylabel("Tipper", fontsize=8)
        # hide unused tipper panels (when n_comp < 4)
        for ti in range(4, len(ax_tip)):
            ax_tip[ti].set_visible(False)

    # ── 8. figure title ───────────────────────────────────────────────────
    fig.suptitle(
        title or st_name,
        fontsize=11, fontweight="bold", y=1.01,
    )
    return fig
