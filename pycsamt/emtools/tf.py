from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import FancyArrowPatch

from ._core import (
    ensure_sites,
    _axes_list,
    _iter_items,
    _get_t_block,
    _get_z_block,
    _name,
    hide_polar_radius_labels,
)

from .tensor import build_phase_tensor_table
from ..api._rose_style import _UNSET, RoseStyle, resolve_rose_style
from ..api.style import PYCSAMT_STYLE
from ..api.control import PYCSAMT_CONTROL
from ..api.section import PYCSAMT_SECTION, SectionStyle
from ..api.labels import LOG10_PERIOD_LABEL
from ..api.station import PYCSAMT_STATION_RENDERING
from ..api.plot import add_colorbar, add_polar_colorbar


# ------------------------------- helpers -------------------------------- #

def _pick_station(S, station: Optional[str]) -> Tuple[str, Any]:
    pool = {}
    for i, ed in enumerate(_iter_items(S)):
        pool[_name(ed, i)] = ed
    if not pool:
        raise RuntimeError("no sites")
    if station is None:
        station = sorted(pool.keys())[0]
    ed = pool.get(station, None)
    if ed is None:
        raise RuntimeError("station not found")
    return station, ed


def _bands_from_periods(
    per: np.ndarray,
    *,
    bands: Optional[Sequence[Tuple[float, float]]] = None,
    n_bands: int = 3,
) -> List[np.ndarray]:
    if bands:
        ms = []
        for (lo, hi) in bands:
            ms.append((per >= float(lo)) & (per <= float(hi)))
        return ms
    # quantile bands
    q = np.linspace(0, 1, num=int(n_bands) + 1)
    bb = np.quantile(per, q)
    ms = []
    for i in range(len(bb) - 1):
        lo, hi = bb[i], bb[i + 1]
        # include low edge, exclude high (except last)
        if i == len(bb) - 2:
            ms.append((per >= lo) & (per <= hi))
        else:
            ms.append((per >= lo) & (per < hi))
    return ms


def _station_xy(ed: Any, i: int) -> Tuple[float, float]:
    for kx, ky in [
        ("east", "north"), ("easting", "northing"),
        ("x", "y"), ("lon", "lat"),
    ]:
        x = getattr(ed, kx, None); y = getattr(ed, ky, None)
        if isinstance(x, (int, float)) and isinstance(y, (int, float)):
            return float(x), float(y)
    return float(i), 0.0  # fallback: index on a line


def _nearest_idx(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    j = np.searchsorted(x, y)
    j = np.clip(j, 1, x.size - 1)
    use_left = np.abs(y - x[j - 1]) <= np.abs(y - x[j])
    j[use_left] -= 1
    return j


def _pt_angle(
    S, station: str, per: np.ndarray
) -> Optional[np.ndarray]:
    if build_phase_tensor_table is None:
        return None
    tb = build_phase_tensor_table(
        S, recursive=False, on_dup="replace",
        strict=False, verbose=0,
    )
    if getattr(tb, "empty", False):
        return None
    sdf = tb[tb["station"] == station]
    if sdf.empty:
        return None
    p_ref = sdf["period"].to_numpy(dtype=float)
    for col in ("azimuth", "strike", "phi", "theta"):
        if col in sdf.columns:
            phi = sdf[col].to_numpy(dtype=float)
            j = _nearest_idx(p_ref, per)
            return phi[j]
    return None


# ------------------------- 11) Tipper hodograms ------------------------- #

def plot_tipper_hodograms(
    sites: Any,
    *,
    station: Optional[str] = None,
    bands: Optional[Sequence[Tuple[float, float]]] = None,
    n_bands: int = 3,
    normalize: bool = False,
    colors: Optional[Sequence[str]] = None,
    marker: str = "o",
    ms: float = 3.0,
    lw: float = 1.0,
    ls: str = "-",
    unit_circle: bool = True,
    axes=None,
    figsize: Tuple[float, float] = (6.4, 3.2),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    st, ed = _pick_station(S, station)
    T, t, fr = _get_t_block(ed)
    axes_given = _axes_list(axes, 1) if axes is not None else None
    if T is None or t is None:
        if axes_given is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no tipper", ha="center", va="center")
        return fig
    per = 1.0 / fr
    Ms = _bands_from_periods(per, bands=bands, n_bands=n_bands)
    if colors is None:
        cols = [plt.cm.viridis(c) for c in np.linspace(0, 1, len(Ms))]
    else:
        cols = list(colors)[: len(Ms)]
    axes_given = _axes_list(axes, 2) if axes is not None else None
    if axes_given is None:
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(1, 2, wspace=0.25)
        axX = fig.add_subplot(gs[0]); axY = fig.add_subplot(gs[1])
    else:
        axX, axY = axes_given
        fig = axX.figure
    for k, m in enumerate(Ms):
        if not np.any(m):
            continue
        tx = t[m, 0]; ty = t[m, 1]
        X = np.real(tx); Y = np.imag(tx)
        U = np.real(ty); V = np.imag(ty)
        if normalize:
            s = np.nanpercentile(
                np.hypot(np.r_[X, U], np.r_[Y, V]), 95
            ) + 1e-24
            X, Y, U, V = X / s, Y / s, U / s, V / s
        col = cols[k]
        axX.plot(X, Y, ls=ls, lw=lw, color=col)
        axX.scatter(X, Y, s=12 * ms, color=col)
        axY.plot(U, V, ls=ls, lw=lw, color=col)
        axY.scatter(U, V, s=12 * ms, color=col)
    for ax, lab in [(axX, "Tx"), (axY, "Ty")]:
        ax.axhline(0, color="0.85", lw=0.8)
        ax.axvline(0, color="0.85", lw=0.8)
        if unit_circle:
            th = np.linspace(0, 2 * np.pi, 256)
            ax.plot(np.cos(th), np.sin(th), ":", color="0.7", lw=0.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Real")
        ax.set_ylabel("Imag")
        ax.set_title(f"{st} • {lab}")
        ax.grid(True, alpha=0.25)
    return fig


# --------- 12) Induction arrows + phase-tensor strike overlay ----------- #

def _arrow_from_tipper(
    t: np.ndarray,
    *,
    convention: str = "park",  # park|wiese|real|imag
) -> np.ndarray:
    tx = t[:, 0]; ty = t[:, 1]
    if convention == "real":
        vx = np.real(tx); vy = np.real(ty)
    else:
        # default: use imaginary (Parkinson-style)
        vx = -np.imag(tx); vy = -np.imag(ty)
        if convention == "wiese":
            # 90° rotation of Parkinson vector
            vx, vy = -vy, vx
        if convention == "imag":
            pass  # already imag-based
    return np.vstack([vx, vy]).T


def plot_induction_arrows(
    sites: Any,
    *,
    periods: Sequence[float] = (1.0,),
    convention: str = "park",   # park|wiese|real|imag
    scale: float = 1.0,         # quiver scale factor
    normalize: bool = True,
    strike_ticks: bool = True,
    tick_len: float = 0.25,
    figsize: Tuple[float, float] = (7.2, 3.4),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    # collect site positions and arrows per requested period
    sts, XY, AR = [], [], []
    per_layers: List[Tuple[float, np.ndarray]] = []
    for p in periods:
        xs, ys, u, v = [], [], [], []
        for i, ed in enumerate(_iter_items(S)):
            T, t, fr = _get_t_block(ed)
            if T is None or t is None:
                continue
            per = 1.0 / fr
            j = _nearest_idx(per, np.array([p], float))
            tx = t[j[0]]
            vec = _arrow_from_tipper(
                np.asarray([tx]), convention=convention
            )[0]
            x, y = _station_xy(ed, i)
            xs.append(x); ys.append(y)
            u.append(vec[0]); v.append(vec[1])
            if p == periods[0]:
                sts.append(_name(ed, i))
        per_layers.append(
            (p, np.vstack([np.array(xs), np.array(ys),
                           np.array(u), np.array(v)]))
        )
    if not per_layers:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no tipper", ha="center", va="center")
        return ax

    # normalization
    if normalize:
        all_mag = []
        for _, L in per_layers:
            mag = np.hypot(L[2], L[3])
            all_mag.append(mag)
        m95 = np.nanpercentile(np.hstack(all_mag), 95) + 1e-24
        for k in range(len(per_layers)):
            p, L = per_layers[k]
            L[2] /= m95; L[3] /= m95
            per_layers[k] = (p, L)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    # axis as map or profile depending on coords
    XY0 = per_layers[0][1]
    if np.unique(XY0[1]).size == 1:  # one line
        ax.set_xlabel("Station index / x")
        ax.set_ylabel("Arrow (arb.)")
    else:
        ax.set_xlabel("East / lon"); ax.set_ylabel("North / lat")

    # draw layers, one color per period
    cols = plt.cm.viridis(
        np.linspace(0.15, 0.85, len(per_layers))
    )
    for (c, (p, L)) in zip(cols, per_layers):
        ax.quiver(
            L[0], L[1], L[2] * scale, L[3] * scale,
            angles="xy", scale_units="xy", scale=1.0,
            color=c, width=0.003, headlength=4, headaxislength=3,
            minlength=0.0,
        )
    # optional strike ticks from PT
    if strike_ticks and build_phase_tensor_table is not None:
        # one average angle per site (first period layer’s sites)
        per = np.array(periods, dtype=float)
        # use the first period to sample strike
        p0 = float(per[0])
        labs = sts
        th_map = {}
        for i, ed in enumerate(_iter_items(S)):
            st = _name(ed, i)
            if st not in labs:
                continue
            Z, z, fr = _get_z_block(ed)[:3]
            if Z is None:
                continue
            th = _pt_angle(S, st, np.array([p0]))
            if th is None:
                continue
            th_map[st] = float(np.radians(th[0]))
        for i, ed in enumerate(_iter_items(S)):
            st = _name(ed, i)
            if st not in th_map:
                continue
            x, y = _station_xy(ed, i)
            th = th_map[st]
            dx = tick_len * np.cos(th)
            dy = tick_len * np.sin(th)
            ax.plot([x - dx, x + dx], [y - dy, y + dy],
                    "-", color="0.2", lw=1.2, alpha=0.9)

    # cosmetics
    if np.unique(XY0[1]).size == 1:
        # profile: center around zero vertically
        ax.axhline(0.0, color="0.8", lw=0.8)
    ax.grid(True, alpha=0.25)
    leg = [f"P={p:g}s" for p, _ in per_layers]
    fig = ax.figure
    fig.legend(
        [plt.Line2D([], [], color=c, lw=2.0) for c in cols],
        leg, ncol=min(len(leg), 4), frameon=False,
        loc="upper center", bbox_to_anchor=(0.5, 1.02), fontsize=9,
    )
    return ax


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers shared by the new functions
# ─────────────────────────────────────────────────────────────────────────────

def _spine_style(ax) -> None:
    ax.grid(True, which="both", ls=":", lw=0.4, color="0.75", zorder=0)
    ax.set_axisbelow(True)


def _set_map_aspect(ax, ys: np.ndarray, scale: float) -> None:
    """Set aspect ratio without wasteful whitespace.

    For a true 2-D map (stations spread in both x and y) use equal
    scaling.  For a profile (all stations on a near-horizontal line)
    switch to *auto* and enforce a tight y-window around the arrow
    extent so the arrow panel fills the axes without blank borders.
    """
    ax.autoscale(True)
    xlim = ax.get_xlim(); ylim = ax.get_ylim()
    xspan = xlim[1] - xlim[0] + 1e-12
    yspan = ylim[1] - ylim[0] + 1e-12

    if yspan < 0.25 * xspan:
        # Profile case: force a visible y-window ≈ 2× max arrow length
        margin = max(scale * 1.4, yspan * 0.5, xspan * 0.08)
        y_ctr  = 0.5 * (ylim[0] + ylim[1])
        ax.set_ylim(y_ctr - margin, y_ctr + margin)
        ax.set_aspect("auto")
    else:
        ax.set_aspect("equal", adjustable="box")


def _collect_tipper_spectrum(
    S: Any,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, Any]]:
    """Return dicts keyed by station: tipper (nf,2), freq, xy."""
    tip_dict: Dict[str, np.ndarray] = {}
    freq_dict: Dict[str, np.ndarray] = {}
    xy_dict:   Dict[str, Tuple[float, float]] = {}
    for i, ed in enumerate(_iter_items(S)):
        T, t, fr = _get_t_block(ed)
        if T is None or t is None:
            continue
        name = _name(ed, i)
        tip_dict[name]  = np.asarray(t, complex)
        freq_dict[name] = np.asarray(fr, float)
        xy_dict[name]   = _station_xy(ed, i)
    return tip_dict, freq_dict, xy_dict


def _tipper_from_spectra(
    sp_input: Any,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Extract tipper from a Spectra object or collection."""
    from ..seg.spectra import Spectra as _Spectra  # noqa: PLC0415

    if isinstance(sp_input, _Spectra):
        items = {"site": sp_input}
    elif isinstance(sp_input, dict):
        items = sp_input
    elif isinstance(sp_input, (list, tuple)):
        items = {getattr(s, "name", None) or f"site{k}": s
                 for k, s in enumerate(sp_input)}
    else:
        raise TypeError(type(sp_input))

    tip_dict: Dict[str, np.ndarray] = {}
    freq_dict: Dict[str, np.ndarray] = {}
    for name, sp in items.items():
        _, tip = sp.to_Z(estimate_error=False)
        if tip is None or tip.tipper is None:
            continue
        tip_dict[str(name)]  = tip.tipper[:, 0, :]   # (nf, 2)
        freq_dict[str(name)] = tip.freq
    return tip_dict, freq_dict


# ─────────────────────────────────────────────────────────────────────────────
# 13) Map view — real + imaginary induction arrows on station positions
# ─────────────────────────────────────────────────────────────────────────────

def plot_induction_map(
    sites: Any,
    *,
    period: float = 1.0,
    convention: str = "park",
    show_real: bool = True,
    show_imag: bool = True,
    scale: float = _UNSET,
    cmap: str = "plasma",
    clim: Optional[Tuple[float, float]] = None,
    show_colorbar: bool = True,
    reference_arrow: float = 0.1,
    station_labels: bool = True,
    title: str = "",
    figsize: Tuple[float, float] = (8, 7),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    r"""Map-view induction arrows at one period.

    Real (solid) and imaginary (dashed) Parkinson arrows at every
    station, coloured by |T| magnitude.

    Parameters
    ----------
    sites : Sites-like
    period : float
        Target period in seconds.
    convention : {'park', 'wiese', 'real', 'imag'}
    show_real, show_imag : bool
    scale : float or _UNSET
        Arrow scale factor (auto from station spacing).
    cmap : str
    clim : (vmin, vmax) or None
    show_colorbar : bool
    reference_arrow : float
        Length of the scale-bar reference arrow.
    station_labels : bool
    title : str
    figsize, ax : standard

    Returns
    -------
    ax : Axes
    """
    S = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                     strict=strict, verbose=verbose)

    names, xs, ys, re_vecs, im_vecs = [], [], [], [], []
    for i, ed in enumerate(_iter_items(S)):
        T, t, fr = _get_t_block(ed)
        if T is None or t is None:
            continue
        per = 1.0 / fr
        j   = _nearest_idx(per, np.array([float(period)]))[0]
        tx  = t[j, 0]; ty = t[j, 1]
        re_vecs.append([np.real(tx), np.real(ty)])
        im_vecs.append([np.imag(tx), np.imag(ty)])
        x, y = _station_xy(ed, i)
        xs.append(x); ys.append(y); names.append(_name(ed, i))

    if not names:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no tipper data", ha="center", va="center",
                transform=ax.transAxes)
        return ax

    xs = np.array(xs, float); ys = np.array(ys, float)
    re = np.array(re_vecs, float)
    im = np.array(im_vecs, float)

    if scale is _UNSET:
        dist  = (np.diff(xs) ** 2 + np.diff(ys) ** 2).mean() ** 0.5 if len(xs) > 1 else 1.0
        scale = float(dist * 0.4)

    mag  = np.hypot(re[:, 0], re[:, 1])
    norm = mcolors.Normalize(
        vmin=clim[0] if clim else mag.min(),
        vmax=clim[1] if clim else mag.max() + 1e-12,
    )
    cm   = plt.get_cmap(cmap)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.get_figure()

    for k in range(len(names)):
        c = cm(norm(mag[k]))
        if show_real:
            ax.annotate("",
                xy=(xs[k] + re[k, 0] * scale, ys[k] + re[k, 1] * scale),
                xytext=(xs[k], ys[k]),
                arrowprops=dict(arrowstyle="-|>", color=c,
                                lw=1.8, mutation_scale=10))
        if show_imag:
            ax.annotate("",
                xy=(xs[k] + im[k, 0] * scale, ys[k] + im[k, 1] * scale),
                xytext=(xs[k], ys[k]),
                arrowprops=dict(arrowstyle="-|>", color=c, lw=1.2,
                                linestyle="dashed", mutation_scale=8))
        ax.plot(xs[k], ys[k], "v", ms=5, color="0.3", zorder=5)
        if station_labels:
            ax.annotate(names[k], (xs[k], ys[k]),
                        xytext=(3, 4), textcoords="offset points",
                        fontsize=6.5, color="0.4")

    # Reference scale arrow
    x0 = xs.min()
    y0 = ys.min() - 0.15 * (np.ptp(ys) + 1e-6)
    ax.annotate("",
        xy=(x0 + reference_arrow * scale, y0), xytext=(x0, y0),
        arrowprops=dict(arrowstyle="-|>", color="0.2",
                        lw=1.8, mutation_scale=10))
    ax.text(x0 + 0.5 * reference_arrow * scale, y0,
            f"|T|={reference_arrow}", ha="center", va="top",
            fontsize=7, color="0.3")

    if show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        add_colorbar(sm, ax, label="|T|", size="3%", pad=0.04, max_ticks=5)

    from matplotlib.lines import Line2D
    handles = []
    if show_real:
        handles.append(Line2D([], [], color="0.4", lw=1.8,
                              label="Real (Parkinson)"))
    if show_imag:
        handles.append(Line2D([], [], color="0.4", lw=1.2, ls="--",
                              label="Imaginary"))
    if handles:
        ax.legend(handles=handles, fontsize=8, framealpha=0.8, loc="upper right")

    _set_map_aspect(ax, ys, scale)
    ax.set_xlabel("x / Easting  (m)", fontsize=9)
    ax.set_ylabel("y / Northing  (m)", fontsize=9)
    ax.set_title(
        title or f"Induction arrows  —  T = {period:g} s  [{convention}]",
        fontsize=10, pad=6)
    _spine_style(ax)
    return ax


# ─────────────────────────────────────────────────────────────────────────────
# 14) Period section — |T| pseudo-section (station × log-period)
# ─────────────────────────────────────────────────────────────────────────────

def plot_induction_section(
    sites: Any,
    *,
    component: str = "abs",
    n_periods: int = 20,
    cmap: str = "RdBu_r",
    clim: Optional[Tuple[float, float]] = None,
    section: Union[str, SectionStyle] = "pseudosection",
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Period × station pseudo-section coloured by |T| magnitude.

    Parameters
    ----------
    sites : Sites-like
    component : {'real', 'imag', 'abs'}
    n_periods : int
    cmap : str
    clim : (vmin, vmax) or None
    section : str or SectionStyle
    title, figsize, ax : standard

    Returns
    -------
    ax : Axes
    """
    S = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                     strict=strict, verbose=verbose)

    tip_dict, freq_dict, xy_dict = _collect_tipper_spectrum(S)
    if not tip_dict:
        if ax is None: _, ax = plt.subplots()
        ax.text(0.5, 0.5, "no tipper", ha="center", va="center",
                transform=ax.transAxes); return ax

    all_fr = list(freq_dict.values())
    f_min  = max(float(f.min()) for f in all_fr)
    f_max  = min(float(f.max()) for f in all_fr)
    f_grid = np.logspace(np.log10(f_min), np.log10(f_max), n_periods)
    per    = 1.0 / f_grid
    names  = list(tip_dict.keys())
    n_st   = len(names)

    mat = np.full((n_st, n_periods), np.nan)
    for si, name in enumerate(names):
        t, fr = tip_dict[name], freq_dict[name]
        for pi, f0 in enumerate(f_grid):
            j = np.argmin(np.abs(fr - f0))
            tx, ty = t[j, 0], t[j, 1]
            if component == "real":
                mat[si, pi] = np.hypot(np.real(tx), np.real(ty))
            elif component == "imag":
                mat[si, pi] = np.hypot(np.imag(tx), np.imag(ty))
            else:
                mat[si, pi] = np.hypot(np.abs(tx), np.abs(ty))

    y_log = np.log10(per)
    sty   = (section if isinstance(section, SectionStyle)
             else PYCSAMT_SECTION.style_for(str(section)).copy())
    if figsize is None:
        figsize = sty.figsize_for(n_stations=n_st, n_y=n_periods)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        fig = ax.get_figure()
    else:
        fig = ax.get_figure()

    st_x    = np.arange(n_st, dtype=float)
    x_edges = np.r_[st_x[0] - 0.5, st_x + 0.5]
    if len(y_log) > 1:
        dy = np.abs(np.diff(y_log)) / 2.0
        sgn = np.sign(np.diff(y_log))
        y_edges = np.r_[y_log[0] - dy[0],
                        y_log[:-1] + sgn * dy,
                        y_log[-1] + sgn[-1] * dy[-1]]
    else:
        y_edges = np.r_[y_log[0] - 0.2, y_log[0] + 0.2]

    vmin, vmax = clim if clim else (0.0, float(np.nanpercentile(mat, 95)) + 1e-12)
    pc = ax.pcolormesh(x_edges, y_edges, mat.T,
                       cmap=cmap, shading="flat", vmin=vmin, vmax=vmax)
    sty.add_colorbar(pc, ax, label=f"|T| {component}")

    if y_log[0] > y_log[-1]:
        ax.invert_yaxis()

    sty.apply_axis(ax, xlabel="Station", ylabel=LOG10_PERIOD_LABEL)
    sty.apply_stations(ax, st_x, names)
    ax.set_title(title or f"Tipper section  [{component}]",
                 fontsize=10, pad=6)
    _spine_style(ax)
    return ax


# ─────────────────────────────────────────────────────────────────────────────
# 15) Convention comparison — 2×2 panel
# ─────────────────────────────────────────────────────────────────────────────

def plot_induction_convention(
    sites: Any,
    *,
    period: float = 1.0,
    scale: float = _UNSET,
    station_labels: bool = True,
    title: str = "",
    axes=None,
    figsize: Tuple[float, float] = (11, 10),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> np.ndarray:
    r"""2×2 panel: Parkinson/Wiese × Real/Imaginary conventions.

    +-----------------------+-----------------------+
    | Parkinson — Real      | Parkinson — Imaginary |
    +-----------------------+-----------------------+
    | Wiese — Real          | Wiese — Imaginary     |
    +-----------------------+-----------------------+

    Parameters
    ----------
    sites : Sites-like
    period : float
    scale, station_labels, title, figsize : standard

    Returns
    -------
    axes : ndarray of Axes, shape (2, 2)
    """
    S = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                     strict=strict, verbose=verbose)

    names, xs, ys, txs, tys = [], [], [], [], []
    for i, ed in enumerate(_iter_items(S)):
        T, t, fr = _get_t_block(ed)
        if T is None or t is None: continue
        per = 1.0 / fr
        j   = _nearest_idx(per, np.array([float(period)]))[0]
        txs.append(t[j, 0]); tys.append(t[j, 1])
        x, y = _station_xy(ed, i)
        xs.append(x); ys.append(y); names.append(_name(ed, i))

    axes_given = _axes_list(axes, 4) if axes is not None else None
    if not names:
        if axes_given is None:
            fig, axs = plt.subplots(2, 2, figsize=figsize)
            axs_flat = axs.ravel()
        else:
            axs_flat = np.asarray(axes_given, dtype=object)
            fig = axs_flat[0].figure
        for ax in axs_flat:
            ax.text(0.5, 0.5, "no tipper", ha="center", va="center",
                    transform=ax.transAxes)
        return np.asarray(axs_flat, dtype=object).reshape(2, 2)

    xs  = np.array(xs, float); ys  = np.array(ys, float)
    txs = np.array(txs); tys = np.array(tys)

    if scale is _UNSET:
        dist  = (np.diff(xs)**2 + np.diff(ys)**2).mean()**0.5 if len(xs)>1 else 1.0
        scale = float(dist * 0.4)

    park_re = np.vstack([ np.real(txs),  np.real(tys)]).T
    park_im = np.vstack([ np.imag(txs),  np.imag(tys)]).T
    wies_re = np.vstack([-park_re[:,1],   park_re[:,0]]).T
    wies_im = np.vstack([-park_im[:,1],   park_im[:,0]]).T

    def _panel(ax, vecs, label):
        u, v = vecs[:,0]*scale, vecs[:,1]*scale
        mag  = np.hypot(u, v) / (scale+1e-24)
        norm = mcolors.Normalize(vmin=0, vmax=max(mag.max(), 1e-6))
        cm   = plt.get_cmap("plasma")
        for k in range(len(names)):
            c = cm(norm(mag[k]))
            ax.annotate("",
                xy=(xs[k]+u[k], ys[k]+v[k]), xytext=(xs[k], ys[k]),
                arrowprops=dict(arrowstyle="-|>", color=c,
                                lw=1.6, mutation_scale=10))
            ax.plot(xs[k], ys[k], "v", ms=5, color="0.3", zorder=5)
            if station_labels:
                ax.annotate(names[k], (xs[k], ys[k]),
                            xytext=(3,4), textcoords="offset points",
                            fontsize=6, color="0.4")
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_title(label, fontsize=9, pad=5)
        ax.set_xlabel("x  (m)", fontsize=8); ax.set_ylabel("y  (m)", fontsize=8)
        _spine_style(ax)

    axes_given = _axes_list(axes, 4) if axes is not None else None
    if axes_given is None:
        fig, axs = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)
    else:
        axs = np.asarray(axes_given, dtype=object).reshape(2, 2)
        fig = axs.ravel()[0].figure
    _panel(axs[0,0], park_re, "Parkinson — Real")
    _panel(axs[0,1], park_im, "Parkinson — Imaginary")
    _panel(axs[1,0], wies_re, "Wiese — Real")
    _panel(axs[1,1], wies_im, "Wiese — Imaginary")
    fig.suptitle(title or f"Induction arrow conventions  —  T = {period:g} s",
                 fontsize=11, y=1.01)
    return axs


# ─────────────────────────────────────────────────────────────────────────────
# 16) Polar tipper plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_tipper_polar(
    sites: Any,
    *,
    station: Optional[str] = None,
    component: str = "real",
    cmap: str = _UNSET,
    lw: float = _UNSET,
    alpha: float = _UNSET,
    title: str = "",
    figsize: Tuple[float, float] = (5.5, 5.5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Polar view: tipper azimuth (angle) and magnitude (radius) vs period.

    Each frequency is one scatter point; colour encodes log₁₀(period).
    North (0°) = up, clockwise positive, following geomagnetic convention.

    Parameters
    ----------
    sites : Sites-like
    station : str or None
    component : {'real', 'imag', 'abs'}
    cmap : str or _UNSET
    lw, alpha : float or _UNSET
    title, figsize, ax : standard

    Returns
    -------
    ax : polar Axes
    """
    _ml = PYCSAMT_STYLE.multiline
    if lw    is _UNSET: lw    = _ml.lw
    if alpha is _UNSET: alpha = _ml.alpha
    if cmap  is _UNSET: cmap  = "viridis"

    S  = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                      strict=strict, verbose=verbose)
    st, ed = _pick_station(S, station)
    T, t, fr = _get_t_block(ed)
    if T is None or t is None:
        if ax is None:
            fig = plt.figure(figsize=figsize); ax = fig.add_subplot(111, polar=True)
        hide_polar_radius_labels(ax)
        ax.set_title("no tipper"); return ax

    per = 1.0 / fr
    tx  = t[:, 0]; ty = t[:, 1]

    if component == "real":
        u, v = np.real(tx), np.real(ty)
    elif component == "imag":
        u, v = np.imag(tx), np.imag(ty)
    else:
        u, v = np.abs(tx), np.abs(ty)

    azimuth   = np.arctan2(v, u)
    magnitude = np.hypot(u, v)
    log_per   = np.log10(np.maximum(per, 1e-9))
    norm      = mcolors.Normalize(vmin=log_per.min(), vmax=log_per.max())

    if ax is None:
        fig = plt.figure(figsize=figsize); ax = fig.add_subplot(111, polar=True)

    sc = ax.scatter(azimuth, magnitude, c=log_per, cmap=cmap, norm=norm,
                    s=30, alpha=alpha, edgecolors="none", zorder=4)
    order = np.argsort(per)
    ax.plot(azimuth[order], magnitude[order],
            lw=lw*0.7, alpha=alpha*0.5, color="0.6", zorder=3)

    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.grid(True, alpha=0.3)
    hide_polar_radius_labels(ax)
    ax.set_title(title or f"{st} — tipper polar [{component}]",
                 fontsize=10, pad=10)
    cbar = add_polar_colorbar(
        sc,
        ax,
        label=r"$\log_{10}T$ (s)",
        pad=0.10,
        shrink=0.72,
        aspect=18,
        max_ticks=5,
    )
    cbar.ax.tick_params(labelsize=7)
    return ax


# ─────────────────────────────────────────────────────────────────────────────
# 17) Induction-arrow rose diagram
# ─────────────────────────────────────────────────────────────────────────────

def plot_induction_rose(
    sites: Any,
    *,
    component: str = "real",
    pband: Optional[Tuple[float, float]] = None,
    nbins: int = 36,
    style=_UNSET,
    title: str = "",
    figsize: Tuple[float, float] = (5, 5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Rose diagram of induction arrow azimuths (all stations & periods).

    Parameters
    ----------
    sites : Sites-like
    component : {'real', 'imag', 'abs'}
    pband : (T_min, T_max) or None
    nbins : int   (default 36 → 10° bins)
    style : RoseStyle or str or _UNSET
    title, figsize, ax : standard

    Returns
    -------
    ax : polar Axes
    """
    rs = (PYCSAMT_STYLE.rose.copy() if style is _UNSET
          else resolve_rose_style(style))
    S  = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                      strict=strict, verbose=verbose)

    azimuths: List[float] = []
    for i, ed in enumerate(_iter_items(S)):
        T, t, fr = _get_t_block(ed)
        if T is None or t is None: continue
        per  = 1.0 / fr
        mask = np.ones(fr.size, bool)
        if pband:
            lo, hi = float(pband[0]), float(pband[1])
            mask &= (per >= lo) & (per <= hi)
        if not mask.any(): continue
        tx, ty = t[mask, 0], t[mask, 1]
        if component == "real":   u, v = np.real(tx), np.real(ty)
        elif component == "imag": u, v = np.imag(tx), np.imag(ty)
        else:                     u, v = np.abs(tx), np.abs(ty)
        azimuths.extend((np.degrees(np.arctan2(v, u)) % 360.0).tolist())

    if not azimuths:
        if ax is None:
            fig = plt.figure(figsize=figsize); ax = fig.add_subplot(111, polar=True)
        hide_polar_radius_labels(ax)
        ax.set_title("no arrows"); return ax

    az = np.array(azimuths, float)
    edges  = np.linspace(0, 360, nbins + 1)
    counts, _ = np.histogram(az, bins=edges)
    counts = counts / counts.sum() * 100.0
    theta  = np.radians(0.5 * (edges[:-1] + edges[1:]))
    width  = np.radians(360.0 / nbins) * 0.92
    r_max  = counts.max() or 1.0

    if ax is None:
        fig = plt.figure(figsize=figsize); ax = fig.add_subplot(111, polar=True)

    if rs.bar_style == "gradient":
        _cm = plt.get_cmap(rs.cmap)
        cols = [_cm(c / (r_max + 1e-12)) for c in counts]
    else:
        cols = rs.bar_color

    ax.bar(theta, counts, width=width, color=cols,
           alpha=rs.bar_alpha, edgecolor=rs.bar_edgecolor,
           linewidth=rs.bar_edgelw, zorder=3)
    ax.plot(np.linspace(0, 2*np.pi, 256), np.full(256, r_max),
            lw=rs.outer_ring_lw, color=rs.outer_ring_color, zorder=4)
    if rs.show_mean and az.size > 0:
        mean_az = np.degrees(np.arctan2(
            np.sin(np.radians(az)).sum(), np.cos(np.radians(az)).sum())) % 360.0
        r_m = np.radians(mean_az)
        ax.plot([r_m, r_m+np.pi], [r_max*0.9]*2,
                lw=rs.mean_lw, color=rs.mean_color, ls=rs.mean_ls, zorder=5)

    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.grid(True, alpha=0.25)
    hide_polar_radius_labels(ax)
    ax.set_title(title or f"Induction arrow rose [{component}]",
                 fontsize=10, pad=12)
    return ax


# ─────────────────────────────────────────────────────────────────────────────
# 18) Spectra-direct wrappers (no Sites needed)
# ─────────────────────────────────────────────────────────────────────────────

def plot_induction_map_from_spectra(
    sp_input: Any,
    *,
    period: float = 1.0,
    coords: Optional[Dict[str, Tuple[float, float]]] = None,
    show_real: bool = True,
    show_imag: bool = True,
    scale: float = _UNSET,
    cmap: str = "plasma",
    station_labels: bool = True,
    title: str = "",
    figsize: Tuple[float, float] = (8, 5),
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Map-view induction arrows from :class:`~pycsamt.seg.spectra.Spectra`.

    Parameters
    ----------
    sp_input : Spectra, list, or dict[str, Spectra]
    period : float
    coords : dict[name → (x, y)] or None
        Station positions.  Equidistant line when ``None``.
    show_real, show_imag : bool
    scale, cmap, station_labels, title, figsize, ax : standard

    Returns
    -------
    ax : Axes
    """
    tip_dict, freq_dict = _tipper_from_spectra(sp_input)
    if not tip_dict:
        if ax is None: _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no tipper", ha="center", va="center",
                transform=ax.transAxes); return ax

    names = list(tip_dict.keys())
    if coords is None:
        coords = {n: (float(i), 0.0) for i, n in enumerate(names)}

    xs = np.array([coords.get(n, (i, 0.0))[0] for i, n in enumerate(names)], float)
    ys = np.array([coords.get(n, (i, 0.0))[1] for i, n in enumerate(names)], float)

    re_list, im_list = [], []
    for name in names:
        t  = tip_dict[name]
        fr = freq_dict[name]
        j  = np.argmin(np.abs(1.0/np.maximum(fr,1e-24) - period))
        tx, ty = t[j,0], t[j,1]
        re_list.append([np.real(tx), np.real(ty)])
        im_list.append([np.imag(tx), np.imag(ty)])
    re = np.array(re_list, float)
    im = np.array(im_list, float)

    if scale is _UNSET:
        dist  = abs(np.diff(xs)).mean() if len(xs) > 1 else 1.0
        scale = float(dist * 0.4) if dist > 1e-9 else 0.3

    mag  = np.hypot(re[:,0], re[:,1])
    norm = mcolors.Normalize(vmin=mag.min(), vmax=mag.max()+1e-12)
    cm   = plt.get_cmap(cmap)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.get_figure()

    for k, name in enumerate(names):
        c = cm(norm(mag[k]))
        if show_real:
            ax.annotate("",
                xy=(xs[k]+re[k,0]*scale, ys[k]+re[k,1]*scale),
                xytext=(xs[k], ys[k]),
                arrowprops=dict(arrowstyle="-|>", color=c,
                                lw=1.8, mutation_scale=10))
        if show_imag:
            ax.annotate("",
                xy=(xs[k]+im[k,0]*scale, ys[k]+im[k,1]*scale),
                xytext=(xs[k], ys[k]),
                arrowprops=dict(arrowstyle="-|>", color=c, lw=1.2,
                                linestyle="dashed", mutation_scale=8))
        ax.plot(xs[k], ys[k], "v", ms=5, color="0.3", zorder=5)
        if station_labels:
            ax.annotate(name, (xs[k], ys[k]),
                        xytext=(3,4), textcoords="offset points",
                        fontsize=6.5, color="0.4")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    add_colorbar(sm, ax, label="|T|", size="3%", pad=0.04, max_ticks=5)
    from matplotlib.lines import Line2D
    handles = []
    if show_real: handles.append(Line2D([],[],color="0.4",lw=1.8,label="Real"))
    if show_imag: handles.append(Line2D([],[],color="0.4",lw=1.2,ls="--",label="Imaginary"))
    if handles: ax.legend(handles=handles, fontsize=8, framealpha=0.8)
    _set_map_aspect(ax, ys, scale)
    ax.set_xlabel("Station", fontsize=9)
    ax.set_title(title or f"Induction arrows from spectra  —  T={period:g} s",
                 fontsize=10, pad=6)
    _spine_style(ax)
    return ax


def plot_tipper_polar_from_spectra(
    sp: Any,
    *,
    component: str = "real",
    cmap: str = "viridis",
    title: str = "",
    ax=None,
    figsize: Tuple[float, float] = (5.5, 5.5),
) -> plt.Axes:
    """Polar tipper from a :class:`~pycsamt.seg.spectra.Spectra` object.

    Parameters
    ----------
    sp : Spectra
    component : {'real', 'imag', 'abs'}
    cmap, title, figsize : standard

    Returns
    -------
    ax : polar Axes
    """
    tip_dict, freq_dict = _tipper_from_spectra(sp)
    if not tip_dict:
        if ax is None:
            fig = plt.figure(figsize=figsize); ax = fig.add_subplot(111,polar=True)
        else:
            fig = ax.figure
        hide_polar_radius_labels(ax)
        ax.set_title("no tipper"); return ax

    name  = list(tip_dict.keys())[0]
    t     = tip_dict[name]
    fr    = freq_dict[name]
    per   = 1.0 / np.maximum(fr, 1e-24)
    tx, ty = t[:,0], t[:,1]

    if component == "real":   u, v = np.real(tx), np.real(ty)
    elif component == "imag": u, v = np.imag(tx), np.imag(ty)
    else:                     u, v = np.abs(tx), np.abs(ty)

    azimuth   = np.arctan2(v, u)
    magnitude = np.hypot(u, v)
    log_per   = np.log10(np.maximum(per, 1e-9))
    norm      = mcolors.Normalize(vmin=log_per.min(), vmax=log_per.max())

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax  = fig.add_subplot(111, polar=True)
    else:
        fig = ax.figure
    sc  = ax.scatter(azimuth, magnitude, c=log_per, cmap=cmap, norm=norm,
                     s=30, alpha=0.85, edgecolors="none", zorder=4)
    order = np.argsort(per)
    ax.plot(azimuth[order], magnitude[order], lw=0.8, alpha=0.4, color="0.6")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.grid(True, alpha=0.3)
    hide_polar_radius_labels(ax)
    ax.set_title(title or f"{name} — polar tipper [{component}]",
                 fontsize=10, pad=10)
    add_polar_colorbar(
        sc,
        ax,
        label=r"$\log_{10}T$ (s)",
        pad=0.10,
        shrink=0.72,
        aspect=18,
        max_ticks=5,
    )
    return ax


def plot_induction_rose_from_spectra(
    sp_input: Any,
    *,
    component: str = "real",
    pband: Optional[Tuple[float, float]] = None,
    nbins: int = 36,
    style=_UNSET,
    title: str = "",
    figsize: Tuple[float, float] = (5, 5),
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Rose diagram of induction arrow directions from Spectra objects.

    Parameters
    ----------
    sp_input : Spectra, list, or dict[str, Spectra]
    component : {'real', 'imag', 'abs'}
    pband : (T_min, T_max) or None
    nbins : int
    style : RoseStyle or str or _UNSET
    title, figsize, ax : standard

    Returns
    -------
    ax : polar Axes
    """
    rs = (PYCSAMT_STYLE.rose.copy() if style is _UNSET
          else resolve_rose_style(style))
    tip_dict, freq_dict = _tipper_from_spectra(sp_input)

    azimuths: List[float] = []
    for name, t in tip_dict.items():
        fr  = freq_dict[name]
        per = 1.0 / np.maximum(fr, 1e-24)
        mask = np.ones(fr.size, bool)
        if pband:
            lo, hi = float(pband[0]), float(pband[1])
            mask &= (per >= lo) & (per <= hi)
        if not mask.any(): continue
        tx, ty = t[mask,0], t[mask,1]
        if component == "real":   u, v = np.real(tx), np.real(ty)
        elif component == "imag": u, v = np.imag(tx), np.imag(ty)
        else:                     u, v = np.abs(tx), np.abs(ty)
        azimuths.extend((np.degrees(np.arctan2(v, u)) % 360.0).tolist())

    if not azimuths:
        if ax is None:
            fig = plt.figure(figsize=figsize); ax = fig.add_subplot(111,polar=True)
        hide_polar_radius_labels(ax)
        ax.set_title("no arrows"); return ax

    az    = np.array(azimuths, float)
    edges = np.linspace(0, 360, nbins+1)
    counts, _ = np.histogram(az, bins=edges)
    counts    = counts / counts.sum() * 100.0
    theta     = np.radians(0.5*(edges[:-1]+edges[1:]))
    width     = np.radians(360.0/nbins) * 0.92
    r_max     = counts.max() or 1.0

    if ax is None:
        fig = plt.figure(figsize=figsize); ax = fig.add_subplot(111,polar=True)

    if rs.bar_style == "gradient":
        _cm = plt.get_cmap(rs.cmap)
        cols = [_cm(c/(r_max+1e-12)) for c in counts]
    else:
        cols = rs.bar_color

    ax.bar(theta, counts, width=width, color=cols,
           alpha=rs.bar_alpha, edgecolor=rs.bar_edgecolor,
           linewidth=rs.bar_edgelw, zorder=3)
    ax.plot(np.linspace(0,2*np.pi,256), np.full(256, r_max),
            lw=rs.outer_ring_lw, color=rs.outer_ring_color, zorder=4)
    if rs.show_mean and az.size > 0:
        mean_az = np.degrees(np.arctan2(
            np.sin(np.radians(az)).sum(),
            np.cos(np.radians(az)).sum())) % 360.0
        r_m = np.radians(mean_az)
        ax.plot([r_m, r_m+np.pi], [r_max*0.9]*2,
                lw=rs.mean_lw, color=rs.mean_color, ls=rs.mean_ls, zorder=5)

    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.grid(True, alpha=0.25)
    hide_polar_radius_labels(ax)
    ax.set_title(title or f"Induction rose (spectra) [{component}]",
                 fontsize=10, pad=12)
    return ax



# ─────────────────────────────────────────────────────────────────────────────
# 19) Multi-period induction vector map  (Boukhalfa et al. 2020 style)
# ─────────────────────────────────────────────────────────────────────────────

def _synthetic_background(
    x0: float, x1: float, y0: float, y1: float,
    *,
    nx: int = 120,
    ny: int = 100,
    seed: int = 42,
    elev_min: float = 1000.0,
    elev_max: float = 2200.0,
) -> Tuple[np.ndarray, Tuple[float, float, float, float]]:
    """Return a smooth synthetic terrain image (ny, nx) and its extent."""
    from scipy.ndimage import gaussian_filter as _gf
    rng  = np.random.default_rng(seed)
    raw  = rng.normal(size=(ny, nx))
    raw  = _gf(raw, sigma=12)
    raw  = (raw - raw.min()) / (np.ptp(raw) + 1e-12)
    bg   = raw * (elev_max - elev_min) + elev_min
    extent = (x0, x1, y0, y1)
    return bg, extent


def plot_induction_multiperiod_map(
    sites: Any,
    *,
    periods: Sequence[float] = (1.0, 10.0, 100.0, 1000.0),
    tipper_data: Optional[Dict[str, np.ndarray]] = None,
    convention: str = "park",
    panel_labels: Optional[Sequence[str]] = None,
    background: Optional[np.ndarray] = None,
    background_extent: Optional[Tuple[float, float, float, float]] = None,
    background_cmap: str = "terrain",
    background_alpha: float = 0.75,
    background_clim: Optional[Tuple[float, float]] = None,
    bg_colorbar_label: str = "Elevation  (m)",
    bg_colorbar_side: str = "right",
    show_background_cbar: bool = True,
    arrow_color: str = "black",
    arrow_scale: float = _UNSET,
    arrow_lw: float = 1.6,
    reference_arrow: float = 0.1,
    reference_panel: int = 0,
    reference_label: str = _UNSET,
    show_stations: bool = True,
    station_labels: bool = False,
    annotations: Optional[Dict[str, Any]] = None,
    annotation_color: str = "navy",
    annotation_fontsize: float = 8.0,
    title: str = "",
    axes=None,
    figsize: Any = _UNSET,
    panel_height: float = 3.0,
    panel_width:  float = 8.5,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Tuple[plt.Figure, np.ndarray]:
    r"""Stacked multi-period induction vector map.

    Produces one panel per period (stacked vertically), each showing the
    **real** Parkinson induction vectors on a background colour map
    (elevation, resistivity, or any 2-D raster).  The layout mimics
    the style of Boukhalfa et al. (2020, *GJI*) Fig. 7 and the paper
    figure described in the session above.

    Visual conventions
    ------------------
    * **Arrows**: black solid vectors using the Parkinson (1962) convention
      (real component pointing toward the conductor).
    * **Background**: smooth terrain coloured with *background_cmap*; a
      synthetic gradient is auto-generated when *background* is ``None``.
    * **Reference vector**: drawn in the first (or *reference_panel*)
      sub-plot; labelled with its normalised length.
    * **Shared colorbar**: placed on the **right** edge of the figure,
      spanning all panels.
    * **Panel labels**: ``"(A) 1 s"``, ``"(B) 10 s"``, … placed in the
      lower-left corner of each panel.
    * **Geological / site annotations**: optional blue text labels
      supplied via *annotations*.

    Parameters
    ----------
    sites : Sites-like
        EDI collection accepted by :func:`ensure_sites`.  When the EDIs
        carry no tipper (``Tipper.tipper`` all-zero), pass synthetic
        tipper via *tipper_data*.
    periods : sequence of float
        Target periods in seconds, one panel each.
    tipper_data : dict {period → (n_sites, 2) complex array}, optional
        Explicit tipper override.  Keys must match *periods* (nearest
        match is used).  Column 0 = T_x, column 1 = T_y.  When absent,
        the tipper is read from the EDIs.
    convention : str
        Arrow convention.  ``"park"`` (Parkinson real) is the only one
        used in published induction-vector maps; ``"wiese"`` is supported.
    panel_labels : sequence of str, optional
        Panel corner labels (e.g. ``["(A) 1 s", "(B) 10 s", ...]``).
        Auto-generated when ``None``.
    background : ndarray (ny, nx), optional
        Pre-computed background raster.  A smooth synthetic terrain is
        generated when ``None``.
    background_extent : (x_left, x_right, y_bottom, y_top), optional
        Geographic extent for the background raster.  Auto-inferred from
        station positions when ``None``.
    background_cmap : str
        Colormap for the background.  Default ``"terrain"`` (green→brown
        elevation appearance).
    background_alpha : float
        Background opacity (0–1).
    background_clim : (vmin, vmax) or None
    bg_colorbar_label : str
        Label for the shared colorbar.
    bg_colorbar_side : {'right', 'left', 'top', 'bottom'}, default 'right'
        Side of the figure where the shared background colorbar is
        placed.  Right-side placement is the package default because it
        keeps section axes and station labels visually grouped.
    show_background_cbar : bool
    arrow_color : str
        Single colour for all induction vectors.  Default ``"black"``.
    arrow_scale : float or _UNSET
        Multiplier for arrow length in data units.  Auto-computed from
        the typical station spacing when ``_UNSET``.
    arrow_lw : float
        Arrow shaft line width in points.
    reference_arrow : float
        Normalised length of the scale-reference arrow drawn in the
        *reference_panel* sub-plot.  Default 0.1.
    reference_panel : int
        Sub-plot index (0-based) where the reference arrow appears.
    reference_label : str or _UNSET
        Text for the reference arrow.  Defaults to
        ``"Vector length {reference_arrow}"``.
    show_stations : bool
        Draw a small marker (``▼``) at each station position.
    station_labels : bool
        Annotate station names next to markers.
    annotations : dict {label: (x, y)} or {label: (x, y, {fontsize, color, …})}, optional
        Geological or site annotations drawn in *annotation_color* on
        **every** panel.
    annotation_color : str
        Default colour for *annotations*.
    annotation_fontsize : float
    title : str
        Figure suptitle.
    figsize : (float, float) or _UNSET
        Auto-computed from *panel_height*, *panel_width*, and number
        of panels when omitted.
    panel_height, panel_width : float
        Per-panel size in inches used for auto *figsize*.
    recursive, on_dup, strict, verbose : standard ensure_sites kwargs.

    Returns
    -------
    fig : Figure
    axes : ndarray of Axes, shape (n_periods,)

    Examples
    --------
    Real data with tipper::

        fig, axs = plot_induction_multiperiod_map(
            "site.edi",
            periods=[1, 10, 100, 1000],
        )

    Synthetic tipper supplied explicitly::

        tipper = {1.0:   st_tips_1s,
                  10.0:  st_tips_10s,
                  100.0: st_tips_100s}
        fig, axs = plot_induction_multiperiod_map(
            "profile/*.edi",
            periods=[1, 10, 100],
            tipper_data=tipper,
        )
    """
    import string
    import matplotlib.gridspec as gridspec
    from scipy.ndimage import gaussian_filter as _gf

    S = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                     strict=strict, verbose=verbose)

    # ── Collect station positions ──────────────────────────────────────────
    all_items = list(_iter_items(S))
    st_names: List[str] = []
    st_xy:    List[Tuple[float, float]] = []
    for i, ed in enumerate(all_items):
        st_names.append(_name(ed, i))
        st_xy.append(_station_xy(ed, i))
    st_xy = np.array(st_xy, float)  # (n, 2)
    n_st  = len(st_names)

    axes_given = _axes_list(axes, 1) if axes is not None else None
    if n_st == 0:
        if axes_given is None:
            fig = plt.figure()
            ax  = fig.add_subplot(111)
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no sites", ha="center", va="center",
                transform=ax.transAxes)
        return fig, np.array([ax])

    xs, ys = st_xy[:, 0], st_xy[:, 1]

    # ── Auto scale from station spacing ───────────────────────────────────
    if arrow_scale is _UNSET:
        if n_st > 1:
            dists = np.sqrt(np.diff(xs)**2 + np.diff(ys)**2)
            med   = np.median(dists[dists > 1e-9]) if dists[dists > 1e-9].size else 1.0
        else:
            med = 1.0
        arrow_scale = float(med * 0.35)

    # ── Build background ──────────────────────────────────────────────────
    pad_x = max((np.ptp(xs) + 1e-9) * 0.12, arrow_scale * 2)
    pad_y = max((np.ptp(ys) + 1e-9) * 0.12, arrow_scale * 2)
    ext = (xs.min() - pad_x, xs.max() + pad_x,
           ys.min() - pad_y, ys.max() + pad_y)

    if background is None:
        background, _ = _synthetic_background(
            ext[0], ext[1], ext[2], ext[3],
            elev_min=1000.0, elev_max=2200.0,
        )
    if background_extent is None:
        background_extent = ext

    bg_vmin = (background_clim[0] if background_clim
               else float(np.nanpercentile(background, 2)))
    bg_vmax = (background_clim[1] if background_clim
               else float(np.nanpercentile(background, 98)))

    # ── Panel labels ──────────────────────────────────────────────────────
    n_per = len(periods)
    letters = list(string.ascii_uppercase)
    if panel_labels is None:
        panel_labels = [f"({letters[k]}) {p:g} s"
                        for k, p in enumerate(periods)]

    # ── Reference label ───────────────────────────────────────────────────
    if reference_label is _UNSET:
        reference_label = f"Vector length  {reference_arrow}"

    # ── Tipper per period ─────────────────────────────────────────────────
    def _get_tip(period: float) -> Optional[np.ndarray]:
        """Return (n_st, 2) complex tipper at *period*, or None."""
        if tipper_data is not None:
            # find nearest key
            keys = np.array(list(tipper_data.keys()), float)
            k    = int(np.argmin(np.abs(keys - period)))
            arr  = np.asarray(tipper_data[keys[k]], complex)
            if arr.shape[0] == n_st:
                return arr
        # fall back to EDI tipper
        vecs: List[complex] = []
        for i, ed in enumerate(all_items):
            T, t, fr = _get_t_block(ed)
            if T is None or t is None:
                vecs.append(0.0 + 0.0j)
                continue
            per = 1.0 / fr
            j   = _nearest_idx(per, np.array([float(period)]))[0]
            vecs.append(t[j, 0] + 0.0j)   # Tx only (first component)
        tx_arr  = np.array([v.real for v in vecs], float)
        ty_arr  = np.zeros(n_st, float)    # no Ty if single component
        return (tx_arr + 1j * ty_arr)[:, np.newaxis]   # fall-through

    # ── Figure layout ─────────────────────────────────────────────────────
    if figsize is _UNSET:
        figsize = (panel_width, panel_height * n_per + 0.6)

    axes_given = _axes_list(axes, n_per) if axes is not None else None
    if axes_given is None:
        fig = plt.figure(figsize=figsize)
        gs  = gridspec.GridSpec(
            n_per, 1,
            figure=fig,
            hspace=0.04,                   # minimal gap between panels
            left=0.12, right=0.97,
            top=0.95 if title else 0.97,
            bottom=0.06,
        )
        axs = [fig.add_subplot(gs[k]) for k in range(n_per)]
    else:
        axs = list(axes_given)
        fig = axs[0].figure
    # share axes
    for ax in axs[1:]:
        ax.sharex(axs[0])
        ax.sharey(axs[0])

    # ── Draw each panel ───────────────────────────────────────────────────
    im_handle = None
    for k, (period, ax) in enumerate(zip(periods, axs)):
        # background
        im = ax.imshow(
            background,
            extent=background_extent,
            origin="lower",
            cmap=background_cmap,
            alpha=background_alpha,
            vmin=bg_vmin, vmax=bg_vmax,
            aspect="auto",
            interpolation="bilinear",
            zorder=1,
        )
        if im_handle is None:
            im_handle = im

        # station markers
        if show_stations:
            ax.scatter(xs, ys, marker="v", s=18, color="0.15",
                       zorder=4, linewidths=0)
            if station_labels:
                for ki, (xi, yi) in enumerate(zip(xs, ys)):
                    ax.annotate(st_names[ki], (xi, yi),
                                xytext=(2, 3), textcoords="offset points",
                                fontsize=5, color="0.25")

        # induction vectors
        tip_arr = _get_tip(period)        # (n_st, 2) or (n_st,)
        if tip_arr is not None:
            tip_arr = np.asarray(tip_arr, complex)
            if tip_arr.ndim == 1:
                tip_arr = tip_arr[:, np.newaxis]
            for si in range(n_st):
                tx = complex(tip_arr[si, 0])
                ty = complex(tip_arr[si, 1]) if tip_arr.shape[1] > 1 else 0.0
                if convention == "wiese":
                    ux = -float(np.real(ty))
                    uy =  float(np.real(tx))
                else:   # Parkinson: real component
                    ux = float(np.real(tx))
                    uy = float(np.real(ty))
                mag = np.hypot(ux, uy)
                if mag < 1e-9:
                    continue
                ax.annotate(
                    "",
                    xy=(xs[si] + ux * arrow_scale,
                        ys[si] + uy * arrow_scale),
                    xytext=(xs[si], ys[si]),
                    arrowprops=dict(
                        arrowstyle="-|>",
                        color=arrow_color,
                        lw=arrow_lw,
                        mutation_scale=9,
                    ),
                    zorder=5,
                )

        # reference arrow (first panel only)
        if k == reference_panel:
            rx  = background_extent[0] + (background_extent[1] -
                                           background_extent[0]) * 0.05
            ry  = background_extent[2] + (background_extent[3] -
                                           background_extent[2]) * 0.80
            ax.annotate(
                "",
                xy=(rx + reference_arrow * arrow_scale, ry),
                xytext=(rx, ry),
                arrowprops=dict(arrowstyle="-|>", color="black",
                                lw=arrow_lw + 0.2, mutation_scale=11),
                zorder=6,
            )
            ax.text(rx + reference_arrow * arrow_scale * 0.5,
                    ry + (background_extent[3] - background_extent[2]) * 0.035,
                    reference_label,
                    ha="center", va="bottom", fontsize=7,
                    color="black", fontweight="bold", zorder=6)

        # geological / site annotations
        if annotations:
            for label, info in annotations.items():
                if isinstance(info, (list, tuple)) and len(info) >= 2:
                    ax_x, ax_y = float(info[0]), float(info[1])
                    kw = info[2] if len(info) > 2 and isinstance(info[2], dict) else {}
                    fc = kw.pop("color",     annotation_color)
                    fs = kw.pop("fontsize",  annotation_fontsize)
                    fw = kw.pop("fontweight", "bold")
                    ax.text(ax_x, ax_y, label,
                            color=fc, fontsize=fs, fontweight=fw,
                            ha="center", va="center", zorder=7, **kw)

        # panel label in lower-left corner
        ax.text(
            0.015, 0.06, panel_labels[k],
            transform=ax.transAxes,
            fontsize=9, color="white", fontweight="bold",
            va="bottom", ha="left", zorder=8,
            bbox=dict(boxstyle="round,pad=0.2", fc="black",
                      alpha=0.55, ec="none"),
        )

        # clean ticks: only bottom panel keeps x labels
        if k < n_per - 1:
            plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(labelsize=7)
        ax.set_xlim(background_extent[0], background_extent[1])
        ax.set_ylim(background_extent[2], background_extent[3])

    # x-axis label on bottom panel only
    axs[-1].set_xlabel("x  (m)", fontsize=8)
    for ax in axs:
        ax.set_ylabel("y  (m)", fontsize=8)

    # ── Shared background colorbar ────────────────────────────────────────
    if show_background_cbar and im_handle is not None:
        bg_colorbar_side = str(bg_colorbar_side).lower()
        if bg_colorbar_side not in {"right", "left", "top", "bottom"}:
            msg = (
                "bg_colorbar_side must be 'right', 'left', 'top', "
                "or 'bottom'."
            )
            raise ValueError(msg)
        vertical = bg_colorbar_side in {"right", "left"}
        cbar = fig.colorbar(
            im_handle,
            ax=axs,
            location=bg_colorbar_side,
            shrink=0.75,
            pad=0.04 if vertical else 0.08,
            fraction=0.025,
            aspect=30,
        )
        cbar.set_label(
            bg_colorbar_label,
            rotation=90 if vertical else 0,
            labelpad=10,
            fontsize=8,
        )
        cbar.ax.tick_params(labelsize=7)

    if title:
        fig.suptitle(title, fontsize=11, y=0.99)

    return fig, np.array(axs)
