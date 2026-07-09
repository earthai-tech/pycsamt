# pycsamt/emtools/strike.py
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import re
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch as _Patch


import numpy as np
import pandas as pd

from ._core import (
    ensure_sites,
    _axes_list,
    _iter_items,
    _apply_each,
    _get_z_block,
    _get_t_block,
    _name,
    hide_polar_radius_labels,
)
from ..site import edit as _edit
from .tensor import build_phase_tensor_table
from ..api._rose_style import RoseStyle, resolve_rose_style, _UNSET
from ..api.labels import LOG10_PERIOD_LABEL
from ..api.station import PYCSAMT_STATION_RENDERING

# -------------------------- small helpers ------------------------------- #

def _rotmat(deg: float) -> np.ndarray:
    th = np.radians(float(deg))
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, s], [-s, c]], dtype=float)


def _rotate_tensor(z: np.ndarray, deg: float) -> np.ndarray:
    R = _rotmat(deg)
    Rt = R.T
    out = np.empty_like(z)
    for i in range(z.shape[0]):
        out[i] = R @ z[i] @ Rt
    return out


def _diag_offdiag_norm(z: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    d = np.sqrt(np.abs(z[:, 0, 0]) ** 2 + np.abs(z[:, 1, 1]) ** 2)
    o = np.sqrt(np.abs(z[:, 0, 1]) ** 2 + np.abs(z[:, 1, 0]) ** 2)
    return d, o


def _score(z: np.ndarray, kind: str) -> np.ndarray:
    d, o = _diag_offdiag_norm(z)
    if kind == "diag_ratio":
        return d / (o + 1e-12)
    if kind == "offdiag_neg":
        return -o
    if kind == "det_diag":
        det = np.abs(z[:, 0, 0] * z[:, 1, 1])
        return det / (o + 1e-12)
    return d / (o + 1e-12)


def _wrap90(a: np.ndarray) -> np.ndarray:
    x = ((a + 90.0) % 180.0) - 90.0
    return x


def _unwrap_deg_180(a: np.ndarray) -> np.ndarray:
    x = a.copy().astype(float)
    for i in range(1, x.size):
        d = x[i] - x[i - 1]
        if d > 90.0:
            x[i:] -= 180.0
        elif d < -90.0:
            x[i:] += 180.0
    return x


def _band_edges(p: np.ndarray, band: Optional[Tuple[float, float]]):
    if band is None:
        lo = max(1e-6, float(np.nanmin(p)))
        hi = float(np.nanmax(p))
        return lo, hi
    return float(band[0]), float(band[1])


def _site_lonlat(ed: Any) -> Tuple[Optional[float], Optional[float]]:
    """Return ``(lon, lat)`` for a Site/EDI object, or ``(None, None)``.

    Real ``Site`` objects expose coordinates only via ``.coords``
    (returning ``(lat, lon, elev)``), not flat ``.lon``/``.lat``
    attributes — checking the latter alone silently treats every real
    station as having no coordinates.
    """
    x = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
    y = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
    if x is None or y is None:
        coords = getattr(ed, "coords", None)
        if coords is not None and len(coords) >= 2:
            y = y if y is not None else coords[0]
            x = x if x is not None else coords[1]
    if x is None or y is None:
        _head = getattr(getattr(ed, "edi", None), "sections", {}).get("head")
        if _head is not None:
            y = y if y is not None else getattr(_head, "lat", None)
            x = x if x is not None else (
                getattr(_head, "long", None) or getattr(_head, "lon", None)
            )
    x = float(x) if x is not None else None
    y = float(y) if y is not None else None
    return x, y


# ----------------------- strike by sweep (Z) ----------------------------- #

def estimate_strike_sweep(
    sites: Any,
    *,
    angles: np.ndarray = np.arange(-90.0, 91.0, 1.0),
    metric: str = "diag_ratio",
    band: Optional[Tuple[float, float]] = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    rows: List[Dict[str, float]] = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        p = 1.0 / fr
        lo, hi = _band_edges(p, band)
        m = (p >= lo) & (p <= hi)
        if not np.any(m):
            continue
        zB = z[m]
        best = np.zeros(zB.shape[0], dtype=float)
        val = np.full(zB.shape[0], np.inf, dtype=float)
        for a in angles:
            zr = _rotate_tensor(zB, a)
            sc = _score(zr, metric)
            upd = sc < val
            val[upd] = sc[upd]
            best[upd] = a
        best_u = _unwrap_deg_180(best)
        ang = _wrap90(np.nanmedian(best_u))
        iqr = np.nanpercentile(best_u, 75) - np.nanpercentile(
            best_u, 25
        )
        rows.append(
            dict(
                station=st, ang=ang, iqr=iqr,
                lo=lo, hi=hi, n=int(zB.shape[0]),
            )
        )
    cols = ["station", "ang", "iqr", "lo", "hi", "n"]
    return pd.DataFrame.from_records(rows, columns=cols)


# ------------------ strike from phase tensor (theta) --------------------- #

def estimate_strike_phase_tensor(
    sites: Any,
    *,
    band: Optional[Tuple[float, float]] = None,
    robust: bool = True,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    df = build_phase_tensor_table(
        S, recursive=False, on_dup=on_dup,
        strict=False, verbose=verbose,
    )
    if df.empty:
        return pd.DataFrame(
            columns=["station", "ang", "iqr", "lo", "hi", "n"]
        )
    rows: List[Dict[str, float]] = []
    for st, sdf in df.groupby("station"):
        p = sdf["period"].to_numpy()
        lo, hi = _band_edges(p, band)
        m = (p >= lo) & (p <= hi)
        if not np.any(m):
            continue
        th = sdf.loc[m, "theta"].to_numpy(dtype=float)
        th = _unwrap_deg_180(th)
        if robust:
            ang = _wrap90(np.nanmedian(th))
            iqr = np.nanpercentile(th, 75) - np.nanpercentile(
                th, 25
            )
        else:
            ang = _wrap90(float(np.nanmean(th)))
            iqr = float(np.nanstd(th))
        rows.append(
            dict(station=st, ang=ang, iqr=iqr, lo=lo, hi=hi,
                 n=int(np.sum(m)))
        )
    return pd.DataFrame.from_records(
        rows, columns=["station", "ang", "iqr", "lo", "hi", "n"]
    )


# ---------------------- consensus strike (blend) ------------------------- #

def estimate_strike_consensus(
    sites: Any,
    *,
    band: Optional[Tuple[float, float]] = None,
    w_sweep: float = 0.5,
    w_pt: float = 0.5,
    angles: np.ndarray = np.arange(-90.0, 91.0, 1.0),
    metric: str = "diag_ratio",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    t1 = estimate_strike_sweep(
        sites,
        angles=angles, metric=metric,
        band=band, recursive=recursive,
        on_dup=on_dup, strict=strict, verbose=verbose,
    )
    t2 = estimate_strike_phase_tensor(
        sites,
        band=band, robust=True,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if t1.empty and t2.empty:
        return pd.DataFrame(
            columns=["station", "ang", "iqr", "lo", "hi", "n"]
        )
    df = pd.merge(
        t1, t2, on="station", how="outer",
        suffixes=("_sweep", "_pt"),
    )
    def pick(a, b, ws, wp):
        if np.isnan(a) and np.isnan(b):
            return np.nan
        a = 0.0 if np.isnan(a) else a
        b = 0.0 if np.isnan(b) else b
        u = _unwrap_deg_180(np.array([a, b]))
        u = _wrap90(u)
        return float(ws * u[0] + wp * u[1])
    ang = []
    iqr = []
    lo = []
    hi = []
    n = []
    for _, r in df.iterrows():
        ang.append(
            pick(r.get("ang_sweep", np.nan),
                 r.get("ang_pt", np.nan),
                 w_sweep, w_pt)
        )
        i1 = r.get("iqr_sweep", np.nan)
        i2 = r.get("iqr_pt", np.nan)
        iqr.append(
            float(np.nanmedian([i1, i2]))
        )
        lo.append(
            float(np.nanmin([
                r.get("lo_sweep", np.nan), r.get("lo_pt", np.nan)
            ]))
        )
        hi.append(
            float(np.nanmax([
                r.get("hi_sweep", np.nan), r.get("hi_pt", np.nan)
            ]))
        )
        n.append(int(np.nansum([r.get("n_sweep", 0), r.get("n_pt", 0)])))
    out = pd.DataFrame(
        dict(station=df["station"], ang=ang, iqr=iqr,
             lo=lo, hi=hi, n=n)
    )
    return out


# ----------------------- rotation applicators ---------------------------- #

def rotate_to_strike(
    sites: Any,
    *,
    method: str = "consensus",  # consensus|sweep|pt
    band: Optional[Tuple[float, float]] = None,
    angles: np.ndarray = np.arange(-90.0, 91.0, 1.0),
    metric: str = "diag_ratio",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if method == "sweep":
        TB = estimate_strike_sweep(
            S, angles=angles, metric=metric,
            band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    elif method == "pt":
        TB = estimate_strike_phase_tensor(
            S, band=band, robust=True,
            recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    else:
        TB = estimate_strike_consensus(
            S, band=band, angles=angles, metric=metric,
            recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    amap = {r.station: float(r.ang) for _, r in TB.iterrows()}

    def _one(Si):
        ed = next(_iter_items(Si))
        st = getattr(ed, "station", None) or getattr(ed, "name", None)
        ang = float(amap.get(st, 0.0))
        # _edit.rotate looks for a ``.Z`` section (the raw EDI layout);
        # a Site wrapper only exposes ``.z`` directly, so calling it on
        # `ed` (or on the `Si` collection) is a silent no-op for every
        # real station. Site.z reads through to edi.Z.z, so rotating
        # the underlying EDI in place also updates the Site wrapper.
        # _apply_each already handles inplace-vs-copy at the collection
        # level (it deep-copies each item first when inplace=False and
        # only keeps mutations made in place on that copy); forwarding
        # the outer *inplace* flag here would make _edit.rotate return a
        # *new* object whose mutation is then discarded, silently
        # leaving every station unrotated.
        edi = getattr(ed, "edi", None)
        _edit.rotate(edi if edi is not None else ed, ang, inplace=True)
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --------------------- per-frequency strike curve ------------------------ #

def strike_curve_sweep(
    sites: Any,
    *,
    angles: np.ndarray = np.arange(-90.0, 91.0, 1.0),
    metric: str = "diag_ratio",
    smooth: int = 5,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    rows: List[Dict[str, float]] = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        best = np.zeros(z.shape[0], dtype=float)
        val = np.full(z.shape[0], np.inf, dtype=float)
        for a in angles:
            zr = _rotate_tensor(z, a)
            sc = _score(zr, metric)
            upd = sc < val
            val[upd] = sc[upd]
            best[upd] = a
        best = _unwrap_deg_180(best)
        if smooth > 1 and best.size >= smooth:
            k = int(smooth)
            w = np.ones(k) / k
            u = np.convolve(best, w, mode="same")
            best = u
        for f, ang in zip(fr, _wrap90(best)):
            period = np.nan if not np.isfinite(f) or f == 0 else 1.0 / float(f)
            rows.append(
                dict(station=st, freq=float(f), period=float(period), ang=float(ang))
            )
    return pd.DataFrame.from_records(
        rows,
        columns=["station", "freq", "period", "ang"],
    )


def _auto_line(st: str) -> str:
    m = re.match(r"^([A-Za-z]+[0-9]+)", str(st))
    return m.group(1) if m else str(st)


def _axial_mean_deg(a: np.ndarray, w: np.ndarray) -> float:
    # axial mean: double angles, mean vector, halve back
    th = np.radians(2.0 * (a % 180.0))
    cw = np.cos(th) * w
    sw = np.sin(th) * w
    ang = 0.5 * np.degrees(np.arctan2(sw.sum(), cw.sum()))
    ang = (ang + 180.0) % 180.0
    return float(ang)

def plot_strike_rose_by_line(
    sites: Any,
    *,
    groups: dict[str, list[str]] | None = None,
    group_key: str | None = None,
    band: tuple[float, float] | None = None,
    method: str = "consensus",  # consensus|sweep|pt
    bins: int = 36,
    weight: str = "inv_iqr",  # inv_iqr|uniform
    axes=None,
    figsize: tuple[float, float] = (8.6, 4.6),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    # 1) per-station strike + iqr (weight proxy)
    if method == "sweep":
        TB = estimate_strike_sweep(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    elif method == "pt":
        TB = estimate_strike_phase_tensor(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    else:
        TB = estimate_strike_consensus(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    axes_given = _axes_list(axes, 1) if axes is not None else None
    if TB.empty:
        if axes_given is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111, polar=True)
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no strikes", ha="center", va="center")
        return fig
    TB = TB.copy()
    TB["ang"] = (TB["ang"] % 180.0).astype(float)
    TB["w"] = 1.0
    if weight == "inv_iqr":
        TB["w"] = 1.0 / (TB["iqr"].abs().to_numpy() + 1e-6)
    # 2) build group membership
    if groups is None:
        # from attribute on EDI, else from station prefix like E1
        lab = {}
        for i, ed in enumerate(_iter_items(S)):
            st = _name(ed, i)
            if group_key and hasattr(ed, group_key):
                lab[st] = str(getattr(ed, group_key))
            else:
                lab[st] = _auto_line(st)
        groups = {}
        for st, g in lab.items():
            groups.setdefault(g, []).append(st)
    # keep groups with at least 2 stations
    groups = {g: v for g, v in groups.items() if len(v) >= 2}
    if not groups:
        if axes_given is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111, polar=True)
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no groups", ha="center", va="center")
        return fig
    # 3) figure grid
    G = list(groups.keys())
    n = len(G)
    axes_given = _axes_list(axes, n) if axes is not None else None
    fig = plt.figure(figsize=figsize) if axes_given is None else axes_given[0].figure
    for i, g in enumerate(G, 1):
        ax = axes_given[i - 1] if axes_given is not None else fig.add_subplot(1, n, i, polar=True)
        subset = TB[TB["station"].isin(groups[g])]
        if subset.empty:
            ax.text(0.5, 0.5, "empty", ha="center", va="center")
            continue
        ang = subset["ang"].to_numpy(dtype=float)
        w = subset["w"].to_numpy(dtype=float)
        # histogram on 0..180 (axial), then mirror to 0..360
        bins = int(max(12, bins))
        edges = np.linspace(0.0, 180.0, bins + 1)
        h, _ = np.histogram(ang, bins=edges, weights=w)
        cen = 0.5 * (edges[1:] + edges[:-1])
        th = np.radians(np.concatenate([cen, cen + 180.0]))
        rr = np.concatenate([h, h])
        # color by height (nice gradient)
        vmax = rr.max() if rr.size else 1.0
        cols = plt.cm.YlOrRd(rr / (vmax + 1e-12))
        ax.bar(
            th, rr,
            width=np.radians(edges[1] - edges[0]),
            bottom=0.0,
            color=cols,
            edgecolor="none",
        )
        # axial mean + label box
        mu = _axial_mean_deg(ang, w)
        rmax = float(vmax) * 1.05
        for add in (0.0, 180.0):
            ax.plot(
                [np.radians(mu + add), np.radians(mu + add)],
                [0.0, rmax],
                color="crimson",
                lw=2.5,
            )
        ax.text(
            0.08, 0.90, f"{mu:.1f}°",
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.2",
                      fc="white", ec="0.2", lw=0.6),
        )
        # polar cosmetics
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        hide_polar_radius_labels(ax)
        ax.set_title(str(g), pad=12.0)
    fig.tight_layout()
    return fig


# ---- enhanced rose diagram ---------------------------------------------- #

def plot_strike_rose(
    sites: Any,
    *,
    # ── visual style ─────────────────────────────────────────────────────
    style: "str | RoseStyle | None" = "pycsamt",
    # ── data / algorithm ─────────────────────────────────────────────────
    groups: Optional[Dict[str, List[str]]] = None,
    group_key: Optional[str] = None,
    band: Optional[Tuple[float, float]] = None,
    freq_bands: Optional[List[Tuple[float, float]]] = None,
    band_labels: Optional[List[str]] = None,
    band_colors: Optional[List] = None,
    method: str = "consensus",
    bins: int = 36,
    weight: str = "inv_iqr",
    # ── visual overrides (None-like sentinel → taken from *style*) ───────
    bar_style=_UNSET,
    bar_color=_UNSET,
    bar_alpha=_UNSET,
    bar_edgecolor=_UNSET,
    bar_edgelw=_UNSET,
    cmap=_UNSET,
    outer_ring_lw=_UNSET,
    outer_ring_color=_UNSET,
    n_rings=_UNSET,
    ring_color=_UNSET,
    ring_ls=_UNSET,
    ring_lw=_UNSET,
    ring_labels=_UNSET,
    ring_label_angle=_UNSET,
    ring_label_fontsize=_UNSET,
    ring_label_color=_UNSET,
    ring_label_fmt=_UNSET,
    spoke_every=_UNSET,
    spoke_color=_UNSET,
    spoke_ls=_UNSET,
    spoke_lw=_UNSET,
    compass_labels=_UNSET,
    compass_fontsize=_UNSET,
    compass_color=_UNSET,
    compass_fontweight=_UNSET,
    show_mean=_UNSET,
    mean_color=_UNSET,
    mean_lw=_UNSET,
    mean_ls=_UNSET,
    show_secondary=_UNSET,
    secondary_color=_UNSET,
    secondary_ls=_UNSET,
    secondary_lw=_UNSET,
    show_annotation=_UNSET,
    annotation_pos=_UNSET,
    annotation_fontsize=_UNSET,
    annotation_bg=_UNSET,
    annotation_ec=_UNSET,
    show_n_stations=_UNSET,
    # ── layout ───────────────────────────────────────────────────────────
    subplot_size: float = 3.2,
    n_cols: Optional[int] = None,
    axes=None,
    figsize: Optional[Tuple[float, float]] = None,
    suptitle: str = "",
    suptitle_fontsize: float = 10.0,
    tight_layout: bool = True,
    # ── core ─────────────────────────────────────────────────────────────
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Publication-quality rose diagram of geoelectric strike direction.

    Each subplot shows the angular distribution of the estimated MT
    geoelectric strike for one station group (profile line).  Bars are
    drawn with axial symmetry — 0–180° mirrored to 180–360° — to reflect
    the inherent 180° ambiguity of geoelectric strike.

    Parameters
    ----------
    sites : any
        EDI paths, EDI objects, or SitesCollection accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    groups : dict[str, list[str]], optional
        Explicit map ``{group_label: [station_name, ...], ...}``.
        If *None*, stations are auto-grouped by profile prefix
        (e.g. ``"E1S01"`` → group ``"E1"``).
    group_key : str, optional
        EDI attribute to read as group label when *groups* is *None*.
    band : (float, float), optional
        Period band ``(lo_s, hi_s)`` in seconds for strike estimation.
        *None* uses all available frequencies.
    freq_bands : list of (float, float), optional
        Period sub-bands used for ``bar_style="bands"``.  Each tuple is
        ``(lo_s, hi_s)``; one histogram per band is stacked.
    band_labels : list[str], optional
        Legend labels matching *freq_bands* (one per band).
    band_colors : list, optional
        Bar colours for each *freq_bands* entry (any matplotlib colour
        spec).  Defaults to ``tab10`` palette.
    method : {"consensus", "sweep", "pt"}
        Strike estimation method — see :func:`estimate_strike_consensus`.
    bins : int
        Number of histogram bins over 0–180°, mirrored to 360°.
    weight : {"inv_iqr", "uniform"}
        Weighting scheme.  ``"inv_iqr"`` down-weights unstable sites.
    bar_style : {"gradient", "bands", "solid"}
        ``"gradient"`` — bars coloured by height via *cmap* (paper
        style); ``"bands"`` — stacked per-band bars with distinct
        colours; ``"solid"`` — uniform *bar_color*.
    bar_color : str
        Bar fill colour for ``bar_style="solid"``.
    bar_edgecolor : str
        Bar edge colour (``"none"`` → no edge).
    bar_edgelw : float
        Bar edge line-width.
    cmap : str
        Colormap name for ``bar_style="gradient"``.
    outer_ring_lw : float
        Line-width of the bold outer circle.
    outer_ring_color : str
        Colour of the outer circle.
    n_rings : int
        Number of concentric reference rings inside the plot.
    ring_color : str
        Colour of grid rings and radial spokes.
    ring_ls : str
        Line-style of grid rings and radial spokes.
    spoke_every : float
        Angular spacing (degrees) of radial spokes / tick marks.
    compass_labels : {"NESW", "degrees", "none"}
        Labels around the polar perimeter.  ``"NESW"`` shows cardinal
        directions; ``"degrees"`` shows degree values;
        ``"none"`` suppresses all labels.
    compass_fontsize : float
        Font size for compass / degree labels.
    compass_color : str
        Colour of compass labels.
    mean_color : str
        Colour of the mean-direction line.
    mean_lw : float
        Line-width of the mean-direction line.
    mean_ls : str
        Line-style of the mean-direction line.
    show_secondary : bool
        Draw the 180°-conjugate mean line (axial symmetry axis).
    secondary_color : str, optional
        Colour for the conjugate line; defaults to *mean_color*.
    secondary_ls : str
        Line-style for the conjugate line.
    secondary_lw : float, optional
        Line-width for the conjugate line; defaults to *mean_lw*.
    annotation_pos : (float, float)
        Axes-fraction ``(x, y)`` of the strike angle annotation box.
    annotation_fontsize : float
        Font size of the annotation text.
    annotation_bg : str
        Background colour of the annotation box.
    annotation_ec : str
        Edge colour of the annotation box.
    show_n_stations : bool
        Append the station count ``n = N`` to the annotation text.
    subplot_size : float
        Side length (inches) of each polar subplot.
    n_cols : int, optional
        Number of subplot columns.  Defaults to ``len(groups)``.
    figsize : (float, float), optional
        Override the auto-computed figure size.
    suptitle : str
        Figure-level super-title.
    suptitle_fontsize : float
        Font size of the super-title.
    tight_layout : bool
        Call ``fig.tight_layout()`` before returning.
    recursive : bool
        Passed to :func:`ensure_sites`.
    on_dup : str
        Duplicate-handling strategy for :func:`ensure_sites`.
    strict : bool
        Strict mode for :func:`ensure_sites`.
    verbose : int
        Verbosity level.

    Returns
    -------
    matplotlib.figure.Figure
        Figure with one polar axes per station group.

    Examples
    --------
    Basic usage — one rose per profile line, gradient style:

    >>> from pycsamt.emtools import plot_strike_rose
    >>> fig = plot_strike_rose("path/to/edis/")

    Frequency-band decomposition (short / long period stacked):

    >>> fig = plot_strike_rose(
    ...     sites,
    ...     bar_style="bands",
    ...     freq_bands=[(0.001, 0.1), (0.1, 100.0)],
    ...     band_labels=["Short period", "Long period"],
    ... )
    """
    # ── resolve style → fill _UNSET visual params ─────────────────────────
    rs = resolve_rose_style(style)
    def _v(val, attr): return getattr(rs, attr) if val is _UNSET else val
    bar_style        = _v(bar_style,        "bar_style")
    bar_color        = _v(bar_color,        "bar_color")
    bar_alpha        = _v(bar_alpha,        "bar_alpha")
    bar_edgecolor    = _v(bar_edgecolor,    "bar_edgecolor")
    bar_edgelw       = _v(bar_edgelw,       "bar_edgelw")
    cmap             = _v(cmap,             "cmap")
    outer_ring_lw    = _v(outer_ring_lw,    "outer_ring_lw")
    outer_ring_color = _v(outer_ring_color, "outer_ring_color")
    n_rings          = _v(n_rings,          "n_rings")
    ring_color       = _v(ring_color,       "ring_color")
    ring_ls          = _v(ring_ls,          "ring_ls")
    ring_lw          = _v(ring_lw,          "ring_lw")
    ring_labels      = _v(ring_labels,      "ring_labels")
    ring_label_angle = _v(ring_label_angle, "ring_label_angle")
    ring_label_fontsize = _v(ring_label_fontsize, "ring_label_fontsize")
    ring_label_color = _v(ring_label_color, "ring_label_color")
    ring_label_fmt   = _v(ring_label_fmt,   "ring_label_fmt")
    spoke_every      = _v(spoke_every,      "spoke_every")
    spoke_color      = _v(spoke_color,      "spoke_color")
    spoke_ls         = _v(spoke_ls,         "spoke_ls")
    spoke_lw         = _v(spoke_lw,         "spoke_lw")
    compass_labels   = _v(compass_labels,   "compass_labels")
    compass_fontsize = _v(compass_fontsize, "compass_fontsize")
    compass_color    = _v(compass_color,    "compass_color")
    compass_fontweight = _v(compass_fontweight, "compass_fontweight")
    show_mean        = _v(show_mean,        "show_mean")
    mean_color       = _v(mean_color,       "mean_color")
    mean_lw          = _v(mean_lw,          "mean_lw")
    mean_ls          = _v(mean_ls,          "mean_ls")
    show_secondary   = _v(show_secondary,   "show_secondary")
    secondary_color  = _v(secondary_color,  "secondary_color")
    secondary_ls     = _v(secondary_ls,     "secondary_ls")
    secondary_lw     = _v(secondary_lw,     "secondary_lw")
    show_annotation  = _v(show_annotation,  "show_annotation")
    annotation_pos   = _v(annotation_pos,   "annotation_pos")
    annotation_fontsize = _v(annotation_fontsize, "annotation_fontsize")
    annotation_bg    = _v(annotation_bg,    "annotation_bg")
    annotation_ec    = _v(annotation_ec,    "annotation_ec")
    show_n_stations  = _v(show_n_stations,  "show_n")

    S = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                     strict=strict, verbose=verbose)

    # ---- strike estimation -------------------------------------------------
    def _est(b: Optional[Tuple[float, float]]):
        if method == "sweep":
            return estimate_strike_sweep(
                S, band=b, recursive=False, on_dup=on_dup,
                strict=False, verbose=verbose)
        if method == "pt":
            return estimate_strike_phase_tensor(
                S, band=b, recursive=False, on_dup=on_dup,
                strict=False, verbose=verbose)
        return estimate_strike_consensus(
            S, band=b, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose)

    TB = _est(band)
    if TB.empty:
        axes_given = _axes_list(axes, 1) if axes is not None else None
        if axes_given is None:
            fig = plt.figure(figsize=figsize or (4.0, 4.0))
            ax = fig.add_subplot(111, polar=True)
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no strikes", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    TB = TB.copy()
    TB["ang"] = TB["ang"] % 180.0
    TB["w"] = (1.0 / (TB["iqr"].abs() + 1e-6)
               if weight == "inv_iqr" else np.ones(len(TB)))

    # ---- optional per-band tables (bar_style="bands") ----------------------
    use_bands = bar_style == "bands" and bool(freq_bands)
    TB_list: List[pd.DataFrame] = []
    _bc: List = []
    _bl: List[str] = []
    if use_bands:
        n_fb = len(freq_bands)  # type: ignore[arg-type]
        _bc = (list(band_colors) if band_colors is not None
               else list(plt.get_cmap("tab10")(np.linspace(0, 0.8, n_fb))))
        _bl = (list(band_labels) if band_labels is not None
               else [f"{lo:.4g}–{hi:.4g} s"
                     for lo, hi in freq_bands])  # type: ignore[union-attr]
        for fb in freq_bands:  # type: ignore[union-attr]
            tb = _est(fb)
            if not tb.empty:
                tb = tb.copy()
                tb["ang"] = tb["ang"] % 180.0
                tb["w"] = (1.0 / (tb["iqr"].abs() + 1e-6)
                           if weight == "inv_iqr" else np.ones(len(tb)))
            TB_list.append(tb)

    # ---- build groups ------------------------------------------------------
    if groups is None:
        lab: Dict[str, str] = {}
        for ii, ed in enumerate(_iter_items(S)):
            st = _name(ed, ii)
            if group_key and hasattr(ed, group_key):
                lab[st] = str(getattr(ed, group_key))
            else:
                lab[st] = _auto_line(st)
        groups = {}
        for st, g in lab.items():
            groups.setdefault(g, []).append(st)

    groups = {g: v for g, v in groups.items() if len(v) >= 2}
    if not groups:
        all_st = TB["station"].tolist()
        if all_st:
            groups = {"All": all_st}
        else:
            axes_given = _axes_list(axes, 1) if axes is not None else None
            if axes_given is None:
                fig = plt.figure(figsize=figsize or (4.0, 4.0))
                ax = fig.add_subplot(111, polar=True)
            else:
                ax = axes_given[0]
                fig = ax.figure
            ax.text(0.5, 0.5, "no groups", ha="center", va="center",
                    transform=ax.transAxes)
            return fig

    # ---- figure layout -----------------------------------------------------
    G = list(groups.keys())
    n_g = len(G)
    ncols = int(n_cols) if n_cols else n_g
    nrows = int(np.ceil(n_g / ncols))
    if figsize is None:
        figsize = (subplot_size * ncols + 0.4,
                   subplot_size * nrows + (0.6 if suptitle else 0.2))

    axes_given = _axes_list(axes, n_g) if axes is not None else None
    fig = plt.figure(figsize=figsize) if axes_given is None else axes_given[0].figure

    bins_ = int(max(12, bins))
    edges = np.linspace(0.0, 180.0, bins_ + 1)
    dw = np.radians(180.0 / bins_)

    for idx, g in enumerate(G):
        ax = axes_given[idx] if axes_given is not None else fig.add_subplot(nrows, ncols, idx + 1, polar=True)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)

        subset = TB[TB["station"].isin(groups[g])]
        n_sta = len(groups[g])

        if subset.empty:
            ax.text(0.5, 0.5, "empty", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(str(g), pad=10.0)
            continue

        ang = subset["ang"].to_numpy(dtype=float)
        w_arr = subset["w"].to_numpy(dtype=float)

        # ---- compute histograms ------------------------------------------
        h_main, _ = np.histogram(ang, bins=edges, weights=w_arr)
        cen = 0.5 * (edges[1:] + edges[:-1])
        th = np.radians(np.concatenate([cen, cen + 180.0]))

        if use_bands:
            h_stacks: List[np.ndarray] = []
            for tb_b in TB_list:
                if tb_b.empty:
                    h_stacks.append(np.zeros(bins_))
                else:
                    sub_b = tb_b[tb_b["station"].isin(groups[g])]
                    if sub_b.empty:
                        h_stacks.append(np.zeros(bins_))
                    else:
                        ang_b = sub_b["ang"].to_numpy(dtype=float)
                        w_b = sub_b["w"].to_numpy(dtype=float)
                        hb, _ = np.histogram(ang_b, bins=edges, weights=w_b)
                        h_stacks.append(hb)
            h_total = np.sum(h_stacks, axis=0)
        else:
            h_total = h_main

        rr = np.concatenate([h_total, h_total])
        rmax = max(float(rr.max()) if rr.size else 0.0, 1e-6)
        rline = rmax * 1.08  # mean-line tip just beyond tallest bar

        # ---- draw bars --------------------------------------------------
        if use_bands:
            bot = np.zeros(2 * bins_)
            for hb, bc in zip(h_stacks, _bc):
                rr_b = np.concatenate([hb, hb])
                ax.bar(th, rr_b, width=dw, bottom=bot,
                       color=bc, edgecolor=bar_edgecolor,
                       linewidth=bar_edgelw, align="center")
                bot = bot + rr_b
        elif bar_style == "gradient":
            col_vals = plt.get_cmap(cmap)(rr / (rmax + 1e-12))
            ax.bar(th, rr, width=dw, bottom=0.0, color=col_vals,
                   edgecolor=bar_edgecolor, linewidth=bar_edgelw,
                   align="center")
        else:
            ax.bar(th, rr, width=dw, bottom=0.0, color=bar_color,
                   edgecolor=bar_edgecolor, linewidth=bar_edgelw,
                   align="center")

        # ---- mean direction line ----------------------------------------
        mu = _axial_mean_deg(ang, w_arr)
        mu_rad = np.radians(mu)
        if show_mean:
            ax.plot([mu_rad, mu_rad], [0.0, rline],
                    color=mean_color, lw=mean_lw, ls=mean_ls,
                    solid_capstyle="round", zorder=5)
            if show_secondary:
                sc = secondary_color or mean_color
                sl = secondary_lw if secondary_lw is not None else mean_lw
                ax.plot([mu_rad + np.pi, mu_rad + np.pi], [0.0, rline],
                        color=sc, lw=sl, ls=secondary_ls,
                        solid_capstyle="round", zorder=5)

        # ---- annotation box (strike angle + optional station count) -----
        if show_annotation:
            txt = f"{mu:.1f}°"
            if show_n_stations:
                txt += f"\nn={n_sta}"
            ax.text(
                annotation_pos[0], annotation_pos[1], txt,
                transform=ax.transAxes,
                fontsize=annotation_fontsize, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.25",
                          fc=annotation_bg, ec=annotation_ec, lw=0.7),
                zorder=6,
            )

        # ---- radial scale (rings + optional count labels) ---------------
        ax.set_rmax(rline * 1.18)
        if ring_labels is not None:
            r_levels = [float(v) for v in ring_labels]
        else:
            step = rmax / max(1, n_rings)
            r_levels = [step * k for k in range(1, n_rings + 1)] if n_rings > 0 else []
        ax.set_yticks(r_levels)
        hide_polar_radius_labels(ax)

        # ---- angular ticks / compass labels -----------------------------
        spoke_angles = np.arange(0.0, 360.0, float(spoke_every))
        if compass_labels == "NESW":
            cpts = {0: "N", 90: "E", 180: "S", 270: "W"}
            lbls = [cpts.get(int(round(s)) % 360, "")
                    for s in spoke_angles]
        elif compass_labels == "degrees":
            lbls = [f"{int(s)}°" for s in spoke_angles]
        else:
            lbls = [""] * len(spoke_angles)
        ax.set_thetagrids(spoke_angles, labels=lbls)
        ax.tick_params(axis="x", labelsize=compass_fontsize,
                       labelcolor=compass_color, pad=4)
        for lbl in ax.get_xticklabels():
            lbl.set_fontweight(compass_fontweight)

        # ---- grid styling -----------------------------------------------
        ax.yaxis.grid(True, color=ring_color, linestyle=ring_ls,
                      linewidth=ring_lw, alpha=0.8)
        ax.xaxis.grid(True, color=spoke_color, linestyle=spoke_ls,
                      linewidth=spoke_lw, alpha=0.7)

        # ---- bold outer ring --------------------------------------------
        ax.spines["polar"].set_linewidth(outer_ring_lw)
        ax.spines["polar"].set_color(outer_ring_color)

        # ---- subplot title ----------------------------------------------
        ax.set_title(str(g), pad=14.0,
                     fontsize=annotation_fontsize + 1.5, fontweight="bold")

    # ---- figure-level band legend (bands style only) -----------------------
    if use_bands and _bl:
        handles = [_Patch(color=c, label=l) for c, l in zip(_bc, _bl)]
        fig.legend(handles=handles, loc="lower center",
                   ncol=min(len(_bl), 5), fontsize=annotation_fontsize - 0.5,
                   frameon=True, framealpha=0.9,
                   bbox_to_anchor=(0.5, -0.04))

    if suptitle:
        fig.suptitle(suptitle, fontsize=suptitle_fontsize,
                     y=1.02 if nrows == 1 else 1.01)
    if tight_layout:
        fig.tight_layout()
    return fig


# ---- STRIKE VIEWS: ribbon, profile, and map-sticks --------------------- #

def _hsv_rgb(h, s, v):
    hsv = np.stack([h, s, v], axis=-1)
    return mcolors.hsv_to_rgb(hsv)

def plot_strike_ribbon(
    sites: Any,
    *,
    method: str = "sweep",  # sweep|pt|consensus
    win: int = 5,
    show_colorbar: bool = True,
    cbar_ticks: list | None = None,   # None → [0, 45, 90, 135, 180]
    figsize: tuple[float, float] = (9.0, 4.2),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    if method == "sweep":
        df = strike_curve_sweep(
            sites,
            recursive=recursive, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )
    else:
        # consensus/pt → per-station single angle; expand flat
        tb = estimate_strike_consensus(
            sites,
            recursive=recursive, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )
        rows = []
        for _, r in tb.iterrows():
            # fake a thin band so it still renders
            for f in (1e-3, 1e3):
                rows.append(
                    dict(station=r.station, freq=f, ang=r.ang)
                )
        df = pd.DataFrame.from_records(rows)
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no strike", ha="center", va="center")
        return ax
    df = df.copy()
    df["lp"] = np.log10(1.0 / df["freq"].to_numpy())
    sts = list(df["station"].unique())
    X = []
    H = []
    for st in sts:
        s = df[df["station"] == st].sort_values("lp")
        th = s["ang"].to_numpy(dtype=float) % 180.0
        lp = s["lp"].to_numpy(dtype=float)
        h = th / 180.0
        k = max(3, int(win))
        if th.size >= k:
            vv = np.convolve(
                ((th - np.nanmean(th)) ** 2),
                np.ones(k) / k,
                mode="same",
            )
        else:
            vv = np.full_like(th, np.nan)
        v0 = (np.nanpercentile(vv, 5)
              if np.isfinite(vv).any() else 0.0)
        v1 = (np.nanpercentile(vv, 95)
              if np.isfinite(vv).any() else 1.0)
        s_sat = 1.0 - np.clip((vv - v0) / (v1 - v0 + 1e-12),
                              0.0, 1.0)
        H.append(np.vstack([h, s_sat, np.ones_like(h)]))
        X.append(lp)
    ygrid = np.unique(np.concatenate(X))
    img = np.zeros((ygrid.size, len(sts), 3))
    for j, (lp, hs) in enumerate(zip(X, H)):
        i = np.searchsorted(ygrid, lp)
        i = np.clip(i, 0, ygrid.size - 1)
        h, s, v = hs
        rgb = _hsv_rgb(h, s, v)
        for r, k in enumerate(i):
            img[k, j, :] = rgb[r]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.imshow(
        img,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_ylabel(LOG10_PERIOD_LABEL)
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(len(sts), dtype=float),
        sts,
        preset="pseudosection",
        xlim=(-0.5, len(sts) - 0.5),
    )
    yt = np.linspace(0, len(ygrid) - 1, num=min(8, len(ygrid)))
    yv = np.linspace(ygrid.min(), ygrid.max(), num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    if not ax.yaxis_inverted():
        ax.invert_yaxis()

    # ── colorbar: hue → strike angle (0–180°) ────────────────────────────
    if show_colorbar:
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        sm = ScalarMappable(
            cmap=plt.get_cmap("hsv"),
            norm=Normalize(vmin=0.0, vmax=180.0),
        )
        sm.set_array([])
        fig = ax.get_figure()
        cb = fig.colorbar(sm, ax=ax, shrink=0.82, pad=0.015, aspect=22)
        cb.set_label("Strike angle (°)", fontsize=8)
        tks = cbar_ticks if cbar_ticks is not None else [0, 45, 90, 135, 180]
        cb.set_ticks(tks)
        cb.ax.tick_params(labelsize=7)
        # saturation key: white = high variance, saturated = stable
        ax.text(
            1.18, 0.02,
            "Saturation → stability",
            transform=ax.transAxes,
            fontsize=6.5, color="0.40",
            ha="center", va="bottom",
            rotation=90,
        )

    return ax


def plot_strike_profile(
    sites: Any,
    *,
    method: str = "consensus",  # consensus|sweep|pt
    band: tuple[float, float] | None = None,
    sort_by: str = "auto",  # auto|lon|lat|name
    figsize: tuple[float, float] = (8.6, 3.8),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if method == "sweep":
        tb = estimate_strike_sweep(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    elif method == "pt":
        tb = estimate_strike_phase_tensor(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    else:
        tb = estimate_strike_consensus(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    if tb.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no strikes", ha="center", va="center")
        return ax
    def _key(st, ed):
        x, y = _site_lonlat(ed)
        if sort_by == "lon":
            return (1, st) if x is None else (0, float(x))
        if sort_by == "lat":
            return (1, st) if y is None else (0, float(y))
        if sort_by == "name":
            return (0, st)
        # auto: lon then name
        return (0, float(x)) if x is not None else (1, st)
    order = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        order.append((st, _key(st, ed)))
    order = [st for st, _ in sorted(order, key=lambda t: t[1])]
    tb = tb.set_index("station").reindex(order).reset_index()
    ang = tb["ang"].to_numpy(dtype=float) % 180.0
    iq = tb["iqr"].to_numpy(dtype=float)
    x = np.arange(tb.shape[0])
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.plot(x, ang, "-", lw=1.5)
    # IQR ribbon
    lo = (ang - 0.5 * iq)
    hi = (ang + 0.5 * iq)
    ax.fill_between(x, lo, hi, alpha=0.25)
    ax.set_ylim(-5.0, 185.0)
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_ylabel("Strike (deg)")
    ax.set_xlabel("Station")
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=90)
    ax.grid(True, alpha=0.2, which="both")
    return ax


def plot_strike_mapsticks(
    sites: Any,
    *,
    method: str = "consensus",  # consensus|sweep|pt
    band: tuple[float, float] | None = None,
    len_deg: float = 0.02,
    alpha_scale: float = 0.9,  # from confidence
    figsize: tuple[float, float] = (7.8, 6.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if method == "sweep":
        tb = estimate_strike_sweep(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    elif method == "pt":
        tb = estimate_strike_phase_tensor(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    else:
        tb = estimate_strike_consensus(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    if tb.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no strikes", ha="center", va="center")
        return ax
    segs = []
    alphas = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        row = tb[tb["station"] == st]
        if row.empty:
            continue
        lon, lat = _site_lonlat(ed)
        if lat is None or lon is None:
            continue
        a = float(row["ang"].iloc[0]) % 180.0
        c = 1.0 / (float(row["iqr"].iloc[0]) + 1e-6)
        # line segment centered at (lon,lat), axial symmetry
        th = np.radians(a)
        dx = 0.5 * len_deg * np.sin(th)
        dy = 0.5 * len_deg * np.cos(th)
        segs.append([(lon - dx, lat - dy), (lon + dx, lat + dy)])
        alphas.append(
            alpha_scale * np.clip(c / np.nanmax([c, 1.0]), 0.1, 1.0)
        )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    lc = LineCollection(
        segs, colors=[(0.1, 0.1, 0.1, a) for a in alphas],
        linewidths=2.0,
    )
    ax.add_collection(lc)
    xs = [s[0][0] for s in segs] + [s[1][0] for s in segs]
    ys = [s[0][1] for s in segs] + [s[1][1] for s in segs]
    if xs and ys:
        ax.set_xlim(min(xs) - len_deg, max(xs) + len_deg)
        ax.set_ylim(min(ys) - len_deg, max(ys) + len_deg)
    ax.set_xlabel("Lon")
    ax.set_ylabel("Lat")
    ax.set_aspect("equal", adjustable="box")
    return ax


# ---- shared rose-panel renderer ----------------------------------------- #

def _draw_rose_on_ax(
    ax: plt.Axes,
    angles_deg: np.ndarray,
    rs: "RoseStyle",
    *,
    bins: int = 36,
    cmap_override: Optional[str] = None,
    title: str = "",
    title_fc: str = "white",
    title_ec: str = "0.35",
) -> None:
    """Render one rose panel onto an existing polar Axes.

    Parameters
    ----------
    ax : polar Axes
    angles_deg : 1-D array of strike/azimuth angles (degrees, axial 0–180°).
    rs : :class:`~pycsamt.api._rose_style.RoseStyle`
    bins : histogram bins over 0–180°.
    cmap_override : optional colormap name; falls back to ``rs.cmap``.
    title : panel title text.
    title_fc / title_ec : facecolor / edgecolor of the title bbox.
    """
    ang = angles_deg[np.isfinite(angles_deg)] % 180.0
    cmap_name = cmap_override or rs.cmap

    bins_ = int(max(12, bins))
    edges = np.linspace(0.0, 180.0, bins_ + 1)
    dw = np.radians(180.0 / bins_)

    h, _ = np.histogram(ang, bins=edges)
    cen = 0.5 * (edges[1:] + edges[:-1])
    th = np.radians(np.concatenate([cen, cen + 180.0]))
    rr = np.concatenate([h, h])
    rmax = max(float(rr.max()), 1.0)
    rline = rmax * 1.08

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    # ── bars ─────────────────────────────────────────────────────────────────
    if rs.bar_style == "gradient":
        cm_ = plt.get_cmap(cmap_name)
        cols = cm_(rr / (rmax + 1e-12))
        ax.bar(th, rr, width=dw, color=cols,
               edgecolor=rs.bar_edgecolor, linewidth=rs.bar_edgelw,
               alpha=rs.bar_alpha, align="center")
    else:
        ax.bar(th, rr, width=dw, color=rs.bar_color,
               edgecolor=rs.bar_edgecolor, linewidth=rs.bar_edgelw,
               alpha=rs.bar_alpha, align="center")

    # ── concentric rings ─────────────────────────────────────────────────────
    if rs.ring_labels is not None:
        r_levels = [float(v) for v in rs.ring_labels]
    else:
        step = rmax / max(1, rs.n_rings)
        r_levels = [step * k for k in range(1, rs.n_rings + 1)]

    ax.set_rmax(rline * 1.18)
    ax.set_yticks(r_levels)
    hide_polar_radius_labels(ax)

    # ── spokes / compass labels ───────────────────────────────────────────────
    spoke_angles = np.arange(0.0, 360.0, float(rs.spoke_every))
    if rs.compass_labels == "NESW":
        _cpts = {0: "N", 90: "E", 180: "S", 270: "W"}
        lbls = [_cpts.get(int(round(s)) % 360, "") for s in spoke_angles]
    elif rs.compass_labels == "degrees":
        lbls = [f"{int(s)}°" for s in spoke_angles]
    else:
        lbls = [""] * len(spoke_angles)
    ax.set_thetagrids(spoke_angles, labels=lbls)
    ax.tick_params(axis="x", labelsize=rs.compass_fontsize,
                   labelcolor=rs.compass_color, pad=4)
    for lbl in ax.get_xticklabels():
        lbl.set_fontweight(rs.compass_fontweight)

    # ── grid styling ──────────────────────────────────────────────────────────
    ax.yaxis.grid(True, color=rs.ring_color, linestyle=rs.ring_ls,
                  linewidth=rs.ring_lw, alpha=0.8)
    ax.xaxis.grid(True, color=rs.spoke_color, linestyle=rs.spoke_ls,
                  linewidth=rs.spoke_lw, alpha=0.7)
    ax.spines["polar"].set_linewidth(rs.outer_ring_lw)
    ax.spines["polar"].set_color(rs.outer_ring_color)

    # ── mean direction + annotation ───────────────────────────────────────────
    if len(ang) == 0:
        ax.text(0.5, 0.5, "no data",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=rs.annotation_fontsize, color="0.55")
    if len(ang) > 0:
        mu = _axial_mean_deg(ang, np.ones(len(ang)))
        mu_rad = np.radians(mu)
        if rs.show_mean:
            ax.plot([mu_rad, mu_rad], [0.0, rline],
                    color=rs.mean_color, lw=rs.mean_lw, ls=rs.mean_ls,
                    solid_capstyle="round", zorder=5)
            if rs.show_secondary:
                sc = rs.secondary_color or rs.mean_color
                sl = rs.secondary_lw if rs.secondary_lw is not None else rs.mean_lw
                ax.plot([mu_rad + np.pi, mu_rad + np.pi], [0.0, rline],
                        color=sc, lw=sl, ls=rs.secondary_ls,
                        solid_capstyle="round", zorder=5)
        if rs.show_annotation:
            txt = f"{mu:.1f}°"
            if rs.show_n:
                txt += f"\nn={len(ang)}"
            ax.text(
                rs.annotation_pos[0], rs.annotation_pos[1], txt,
                transform=ax.transAxes,
                fontsize=rs.annotation_fontsize,
                va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.25",
                          fc=rs.annotation_bg, ec=rs.annotation_ec, lw=0.7),
                zorder=6,
            )

    # ── title with coloured box ───────────────────────────────────────────────
    ax.set_title(
        title,
        fontsize=rs.annotation_fontsize + 1.5,
        fontweight="bold",
        pad=14,
        bbox=dict(boxstyle="round,pad=0.3",
                  facecolor=title_fc, edgecolor=title_ec, lw=0.8),
    )


# ---- combined Strike-Analysis figure ------------------------------------ #

def plot_strike_analysis(
    sites: Any,
    *,
    # ── visual style ─────────────────────────────────────────────────────
    style: "str | RoseStyle | None" = "pycsamt",
    # ── data / algorithm ─────────────────────────────────────────────────
    band: Optional[Tuple[float, float]] = None,
    bins: int = 36,
    method: str = "sweep",
    # ── per-panel colormap overrides (None → rs.cmap from style) ─────────
    cmap_z: Optional[str] = None,
    cmap_pt: Optional[str] = None,
    cmap_tipper: Optional[str] = None,
    # ── title box colours ─────────────────────────────────────────────────
    title_fc_z: str = "#ffe0e0",
    title_fc_pt: str = "#ffffd0",
    title_fc_tipper: str = "#d5e8ff",
    title_ec: str = "0.35",
    # ── layout ───────────────────────────────────────────────────────────
    axes=None,
    figsize: Optional[Tuple[float, float]] = None,
    subplot_size: float = 3.8,
    suptitle: str = "",
    tight_layout: bool = True,
    # ── core ─────────────────────────────────────────────────────────────
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Three-panel rose diagram: Strike (Z), PT Azimuth, and Tipper Strike.

    Produces a publication-quality figure with one polar rose per analysis
    type, analogous to the MTPy ``StrikeAnalysis`` plot.  All three panels
    share the same :class:`~pycsamt.api._rose_style.RoseStyle` so colours
    remain visually consistent.  Each panel carries a coloured title box to
    distinguish the three quantities at a glance.

    Parameters
    ----------
    sites : any
        EDI paths, EDI objects, or
        :class:`~pycsamt.core.base.SitesCollection` accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    style : str, RoseStyle, or None
        Named style preset or :class:`~pycsamt.api._rose_style.RoseStyle`
        instance.  Default ``"pycsamt"`` uses the YlOrRd-gradient,
        crimson-mean-line paper style.
    band : (lo_s, hi_s) or None
        Period window in seconds applied to **all three** panels.
        ``None`` uses all available periods / frequencies.
    bins : int
        Number of histogram bins over 0–180°, mirrored to 0–360°.
        Default 36 → 5° bins.
    method : {"sweep", "pt", "consensus"}
        Strike estimation algorithm for the **Strike (Z)** panel.

        ``"sweep"`` — impedance-tensor rotation sweep
        (calls :func:`estimate_strike_sweep`);
        ``"pt"`` — phase-tensor θ median per station
        (calls :func:`estimate_strike_phase_tensor`);
        ``"consensus"`` — weighted blend of sweep and PT
        (calls :func:`estimate_strike_consensus`).
    cmap_z, cmap_pt, cmap_tipper : str or None
        Colormap name for each panel when ``bar_style="gradient"``.
        ``None`` falls back to the colormap in *style*.
    title_fc_z, title_fc_pt, title_fc_tipper : str
        Facecolour of the title annotation box for each panel.
    title_ec : str
        Edge colour shared by all title boxes.
    figsize : (float, float) or None
        Figure size.  Auto-derived from *subplot_size* when ``None``.
    subplot_size : float
        Side length (inches) of each polar panel when *figsize* is auto.
    suptitle : str
        Figure-level super-title.
    tight_layout : bool
        Call ``fig.tight_layout()`` before returning.
    recursive, on_dup, strict, verbose
        Passed to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    matplotlib.figure.Figure
        Figure with three polar axes: Strike (Z), PT Azimuth, Tipper Strike.

    Examples
    --------
    Default style, all periods:

    >>> from pycsamt.emtools import plot_strike_analysis
    >>> fig = plot_strike_analysis("path/to/edis/")
    >>> fig.savefig("strike_analysis.png", dpi=150, bbox_inches="tight")

    Short-period band, publication style:

    >>> fig = plot_strike_analysis(
    ...     sites,
    ...     band=(0.01, 1.0),
    ...     style="publication",
    ...     suptitle="WILLY AMT — short-period band",
    ... )
    """
    rs = resolve_rose_style(style)

    S = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                     strict=strict, verbose=verbose)

    # ── 1. Z-strike angles (one per station) ────────────────────────────────
    if method == "pt":
        df_z = estimate_strike_phase_tensor(
            S, band=band, recursive=False,
            on_dup=on_dup, strict=False, verbose=verbose,
        )
    elif method == "consensus":
        df_z = estimate_strike_consensus(
            S, band=band, recursive=False,
            on_dup=on_dup, strict=False, verbose=verbose,
        )
    else:
        df_z = estimate_strike_sweep(
            S, band=band, recursive=False,
            on_dup=on_dup, strict=False, verbose=verbose,
        )
    ang_z = (df_z["ang"].to_numpy(float) % 180.0
             if not df_z.empty else np.empty(0))

    # ── 2. PT azimuth angles (per frequency × station) ──────────────────────
    df_pt = build_phase_tensor_table(
        S, recursive=False, on_dup=on_dup, strict=False, verbose=verbose,
    )
    if not df_pt.empty:
        if band is not None:
            lo_, hi_ = float(band[0]), float(band[1])
            m_pt = (df_pt["period"] >= lo_) & (df_pt["period"] <= hi_)
            ang_pt = df_pt.loc[m_pt, "theta"].to_numpy(float) % 180.0
        else:
            ang_pt = df_pt["theta"].to_numpy(float) % 180.0
        ang_pt = ang_pt[np.isfinite(ang_pt)]
    else:
        ang_pt = np.empty(0)

    # ── 3. Tipper-strike angles (per frequency × station) ───────────────────
    _tip_list: List[float] = []
    for i, ed in enumerate(_iter_items(S)):
        _T, t, fr = _get_t_block(ed)
        if t is None or fr is None:
            continue
        per_t = 1.0 / np.where(fr == 0, np.nan, fr)
        mask_t = np.isfinite(per_t)
        if band is not None:
            lo_, hi_ = float(band[0]), float(band[1])
            mask_t &= (per_t >= lo_) & (per_t <= hi_)
        if not mask_t.any():
            continue
        tx = np.real(t[mask_t, 0])   # Re(Tzx)
        ty = np.real(t[mask_t, 1])   # Re(Tzy)
        az = np.degrees(np.arctan2(ty, tx)) % 180.0
        _tip_list.extend(az[np.isfinite(az)].tolist())
    ang_tipper = np.array(_tip_list, float)

    # ── figure ───────────────────────────────────────────────────────────────
    if figsize is None:
        figsize = (subplot_size * 3 + 0.6, subplot_size + 0.5)
    axes_given = _axes_list(axes, 3)
    if axes_given is None:
        fig, axes_arr = plt.subplots(
            1, 3, figsize=figsize, subplot_kw=dict(polar=True),
        )
    else:
        axes_arr = np.asarray(axes_given, dtype=object)
        fig = axes_arr[0].figure

    _draw_rose_on_ax(
        axes_arr[0], ang_z, rs,
        bins=bins, cmap_override=cmap_z,
        title="Strike (Z)",
        title_fc=title_fc_z, title_ec=title_ec,
    )
    _draw_rose_on_ax(
        axes_arr[1], ang_pt, rs,
        bins=bins, cmap_override=cmap_pt,
        title="PT Azimuth",
        title_fc=title_fc_pt, title_ec=title_ec,
    )
    _draw_rose_on_ax(
        axes_arr[2], ang_tipper, rs,
        bins=bins, cmap_override=cmap_tipper,
        title="Tipper Strike",
        title_fc=title_fc_tipper, title_ec=title_ec,
    )

    if suptitle:
        fig.suptitle(suptitle, fontsize=10.0, y=1.02)
    if tight_layout:
        fig.tight_layout()
    return fig
