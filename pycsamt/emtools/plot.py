# pycsamt/emtools/plot.py
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt

from ._core import (
    ensure_sites,
    _iter_items,
    _get_z_block,
    _name,
)

from ..utils.plot import plot_errorbar


# ------------------------------- helpers -------------------------------- #

_COL = {
    "xy": "#1f77b4",  # blue
    "yx": "#d62728",  # red
    "xx": "#2ca02c",  # green
    "yy": "#9467bd",  # purple
}

_IDX = {"x": 0, "y": 1}

def _comp_slice(z: np.ndarray, comp: str) -> np.ndarray:
    a, b = _IDX[comp[0]], _IDX[comp[1]]
    return z[:, a, b]

def _phase_deg(z: np.ndarray) -> np.ndarray:
    return np.degrees(np.angle(z))

def _wrap_phase(ph: np.ndarray, rng: Tuple[float, float]) -> np.ndarray:
    lo, hi = float(rng[0]), float(rng[1])
    w = hi - lo
    return (ph - lo) % w + lo

def _rhoa_from(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    return (0.2 * (np.abs(z) ** 2)) / (fr + 1e-24)

def _err_log10_rhoa(z: np.ndarray, ze: Optional[np.ndarray]) -> np.ndarray:
    if ze is None:
        return None
    # δ log10 ρ ≈ 2/(ln10) · σ/|Z|
    s = np.abs(ze)
    m = np.abs(z) + 1e-24
    return (2.0 / np.log(10.0)) * (s / m)

def _err_log10_mag(z: np.ndarray, ze: Optional[np.ndarray]) -> np.ndarray:
    if ze is None:
        return None
    # δ log10|Z| ≈ σ/(|Z| ln 10)
    s = np.abs(ze)
    m = np.abs(z) + 1e-24
    return s / (m * np.log(10.0))

def _err_phase_deg(z: np.ndarray, ze: Optional[np.ndarray]) -> np.ndarray:
    if ze is None:
        return None
    # small-angle approx: δφ ≈ σ/|Z| (rad)
    s = np.abs(ze)
    m = np.abs(z) + 1e-24
    return np.degrees(s / m)

def _pick(axs, i, j):
    return axs[i * 2 + j]


def _zblk_flex(ed: Any):
    try:
        return _get_z_block(ed, with_errors=True)
    except TypeError:
        return _get_z_block(ed)


# ------------------------------ main plot -------------------------------- #

def plot_sites_panels(
    sites: Any,
    *,
    components: Tuple[str, ...] = ("xy", "yx"),
    quantity: str = "rhoa",          # rhoa|impedance
    x_axis: str = "period",          # period|frequency
    phase_range: Tuple[float, float] | None = (-90.0, 90.0),
    stations: Optional[List[str]] = None,
    ncols: int = 6,
    wspace: float = 0.20,
    hspace: float = 0.08,
    height_ratio: Tuple[int, int] = (2, 1),
    figsize_scale: Tuple[float, float] = (2.6, 2.6),
    colors: Optional[Dict[str, str]] = None,
    marker: str = "o",
    ms: float = 2.5,
    lw: float = 1.2,
    ls: str = "-",
    show_error_bars: bool = True,
    show_legend: bool = False,
    title_fmt: str = "{station}",
    ylim_rhoa: Tuple[float, float] | None = None,
    ylim_phase: Tuple[float, float] | None = None,
    grid: bool = True,
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
    comps = tuple(c.lower() for c in components)
    cmap = {**_COL, **(colors or {})}
    # gather stations
    items = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        if stations and st not in stations:
            continue
        items.append((st, ed))
    if not items:
        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "no sites", ha="center", va="center")
        return fig
    n = len(items)
    ncols = max(1, int(ncols))
    nrows = (n + ncols - 1) // ncols
    W, H = figsize_scale
    fig = plt.figure(
        figsize=(W * ncols, H * nrows),
        constrained_layout=False,
    )
    outer = fig.add_gridspec(
        nrows, ncols, wspace=wspace, hspace=hspace
    )
    axs: List[plt.Axes] = []
    for k in range(nrows * ncols):
        r, c = divmod(k, ncols)
        if k >= n:
            # filler invisible axes
            gs = outer[r, c].subgridspec(
                2, 1, hspace=0.02, height_ratios=height_ratio
            )
            ax1 = fig.add_subplot(gs[0]); ax1.axis("off")
            ax2 = fig.add_subplot(gs[1]); ax2.axis("off")
            axs += [ax1, ax2]
            continue
        gs = outer[r, c].subgridspec(
            2, 1, hspace=0.02, height_ratios=height_ratio
        )
        axR = fig.add_subplot(gs[0])
        axP = fig.add_subplot(gs[1], sharex=axR)
        axs += [axR, axP]

    # draw each station
    for idx, (st, ed) in enumerate(items):
        Z, z, fr = _get_z_block(ed, return_errors=True)
        if isinstance(Z, tuple):
            # legacy tuple from helper fallback
            _, z, fr, ze = Z
        else:
            ze = getattr(getattr(ed, "Z", None), "z_err", None)
        if z is None or fr is None:
            continue
        x = (1.0 / fr) if x_axis == "period" else fr
        axR = _pick(axs, idx, 0)
        axP = _pick(axs, idx, 1)
        axR.set_xscale("log"); axP.set_xscale("log")
        if grid:
            axR.grid(True, alpha=0.25, which="both")
            axP.grid(True, alpha=0.25, which="both")
        for comp in comps:
            col = cmap.get(comp, "k")
            zz = _comp_slice(z, comp)
            ee = None
            if isinstance(ze, np.ndarray) and ze.shape == z.shape:
                ee = _comp_slice(ze, comp)
            # top: rhoa or |Z|
            if quantity == "impedance":
                ymag = np.abs(zz)
                y = np.log10(ymag)
                yerr = _err_log10_mag(zz, ee)
                ylab = "Log10|Z|"
            else:
                rho = _rhoa_from(zz, fr)
                y = np.log10(rho)
                yerr = _err_log10_rhoa(zz, ee)
                ylab = "Log10Rhoa (Ω·m)"
            plot_errorbar(
                axR, x, y, y_err=yerr, color=col,
                marker=marker, ms=ms, ls=ls, lw=lw,
                show_error_bars=show_error_bars,
            )
            # bottom: phase
            ph = _phase_deg(zz)
            if phase_range is not None:
                ph = _wrap_phase(ph, phase_range)
            plot_errorbar(
                axP, x, ph, y_err=_err_phase_deg(zz, ee),
                color=col, marker=marker, ms=ms, ls=ls, lw=lw,
                show_error_bars=False,  # phase bars optional; off
            )
        # cosmetics per panel
        if idx // ncols == nrows - 1:
            axP.set_xlabel(
                "Period (s)" if x_axis == "period" else "Freq (Hz)"
            )
        else:
            axP.set_xlabel("")
        axR.set_ylabel(ylab)
        axP.set_ylabel("Phase (°)")
        if ylim_rhoa:
            axR.set_ylim(*ylim_rhoa)
        if ylim_phase:
            axP.set_ylim(*ylim_phase)
        elif phase_range is not None:
            axP.set_ylim(*phase_range)
        if show_legend and idx == 0:
            labs = [c.upper() for c in comps]
            handles = []
            for c in comps:
                h = axR.plot([], [], color=cmap.get(c, "k"),
                             ls=ls, lw=lw)[0]
                handles.append(h)
            axR.legend(handles, labs, ncol=len(labs), fontsize=8)
        axR.set_title(title_fmt.format(station=st), pad=4)

    # tidy shared ticks
    for i in range(n):
        axR = _pick(axs, i, 0)
        axP = _pick(axs, i, 1)
        if (i % ncols) != 0:
            axR.set_yticklabels([])
            axP.set_yticklabels([])
        if (i // ncols) != (nrows - 1):
            axP.set_xticklabels([])
    return fig


def _pair_by_station(
    sites: Any, new_sites: Any | None
) -> List[Tuple[str, Any, Any | None]]:
    S0 = ensure_sites(sites, recursive=False, strict=False)
    M0: Dict[str, Any] = {}
    for i, ed in enumerate(_iter_items(S0)):
        M0[_name(ed, i)] = ed
    if new_sites is None:
        return [(k, v, None) for k, v in M0.items()]
    S1 = ensure_sites(new_sites, recursive=False, strict=False)
    M1: Dict[str, Any] = {}
    for i, ed in enumerate(_iter_items(S1)):
        M1[_name(ed, i)] = ed
    out: List[Tuple[str, Any, Any | None]] = []
    for st, ed in M0.items():
        out.append((st, ed, M1.get(st, None)))
    return out


def plot_sites_compare(
    sites: Any,
    new_sites: Any | None = None,
    *,
    components: Tuple[str, ...] = ("xy", "yx"),
    quantity: str = "rhoa",      # rhoa|impedance
    x_axis: str = "period",      # period|frequency
    phase_range: Tuple[float, float] | None = (-90.0, 90.0),
    stations: Optional[List[str]] = None,
    ncols_groups: int = 3,       # station groups per row
    group_gap: float = 0.35,     # space between station groups
    pair_wspace: float = 0.06,   # space between raw/after
    hspace: float = 0.06,
    height_ratio: Tuple[int, int] = (2, 1),
    figsize_scale: Tuple[float, float] = (3.0, 3.0),
    colors: Optional[Dict[str, str]] = None,
    marker: str = "o",
    ms: float = 2.5,
    lw: float = 1.2,
    ls: str = "-",
    show_error_bars: bool = True,
    labels: Tuple[str, str] = ("raw", "after"),
    title_group_fmt: str = "{station}",
    title_col_fmt: str = "{tag}",
    show_legend: bool = False,
    ylim_rhoa: Tuple[float, float] | None = None,
    ylim_phase: Tuple[float, float] | None = None,
    grid: bool = True,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    comps = tuple(c.lower() for c in components)
    cmap = {**_COL, **(colors or {})}
    pairs = _pair_by_station(sites, new_sites)
    if stations:
        keep = set(stations)
        pairs = [p for p in pairs if p[0] in keep]
    if not pairs:
        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "no stations",
                ha="center", va="center")
        return fig
    dual = any(p[2] is not None for p in pairs)
    cols_per_grp = 2 if dual else 1
    n = len(pairs)
    ncols_groups = max(1, int(ncols_groups))
    nrows = (n + ncols_groups - 1) // ncols_groups
    W, H = figsize_scale
    fig = plt.figure(
        figsize=(W * ncols_groups * cols_per_grp, H * nrows),
        constrained_layout=False,
    )
    outer = fig.add_gridspec(
        nrows, ncols_groups, wspace=group_gap, hspace=hspace
    )

    # prebuild inner axes
    AxR: List[List[plt.Axes]] = []
    AxP: List[List[plt.Axes]] = []
    for g in range(n):
        r, c = divmod(g, ncols_groups)
        gs = outer[r, c].subgridspec(
            2, cols_per_grp,
            hspace=0.02, wspace=pair_wspace,
            height_ratios=height_ratio,
        )
        axR = [fig.add_subplot(gs[0, j]) for j in range(cols_per_grp)]
        axP = [fig.add_subplot(gs[1, j], sharex=axR[0])
               for j in range(cols_per_grp)]
        AxR.append(axR); AxP.append(axP)

    # draw each group
    for g, (st, ed0, ed1) in enumerate(pairs):
        # group title
        AxR[g][0].set_title(
            title_group_fmt.format(station=st), pad=4
        )
        cols = [(labels[0], ed0)]
        if dual:
            cols.append((labels[1], ed1))
        for j, (tag, ed) in enumerate(cols):
            axR = AxR[g][j]
            axP = AxP[g][j]
            if ed is None:
                axR.axis("off"); axP.axis("off"); continue
            out = _zblk_flex(ed)
            if len(out) == 4:
                _, z, fr, ze = out
            else:
                _, z, fr = out[:3]; ze = None
            if z is None or fr is None:
                axR.axis("off"); axP.axis("off"); continue
            x = (1.0 / fr) if x_axis == "period" else fr
            axR.set_xscale("log"); axP.set_xscale("log")
            if grid:
                axR.grid(True, alpha=0.25, which="both")
                axP.grid(True, alpha=0.25, which="both")
            for comp in comps:
                col = cmap.get(comp, "k")
                zz = _comp_slice(z, comp)
                ee = None
                if isinstance(ze, np.ndarray) and ze.shape == z.shape:
                    ee = _comp_slice(ze, comp)
                if quantity == "impedance":
                    y = np.log10(np.abs(zz))
                    yerr = _err_log10_mag(zz, ee)
                    ylab = "Log10|Z|"
                else:
                    rho = _rhoa_from(zz, fr)
                    y = np.log10(rho)
                    yerr = _err_log10_rhoa(zz, ee)
                    ylab = "Log10Rhoa (Ω·m)"
                plot_errorbar(
                    axR, x, y, y_err=yerr, color=col,
                    marker=marker, ms=ms, ls=ls, lw=lw,
                    show_error_bars=show_error_bars,
                )
                ph = _phase_deg(zz)
                if phase_range is not None:
                    ph = _wrap_phase(ph, phase_range)
                plot_errorbar(
                    axP, x, ph,
                    y_err=_err_phase_deg(zz, ee),
                    color=col, marker=marker, ms=ms,
                    ls=ls, lw=lw, show_error_bars=False,
                )
            # per-column cosmetics
            axR.text(
                0.02, 0.96, title_col_fmt.format(tag=tag),
                ha="left", va="top", transform=axR.transAxes,
                fontsize=9,
            )
            if ylim_rhoa:
                axR.set_ylim(*ylim_rhoa)
            if ylim_phase:
                axP.set_ylim(*ylim_phase)
            elif phase_range is not None:
                axP.set_ylim(*phase_range)
            if (g // ncols_groups) == (nrows - 1):
                axP.set_xlabel(
                    "Period (s)" if x_axis == "period" else "Freq (Hz)"
                )
            else:
                axP.set_xlabel("")
            axR.set_ylabel(ylab); axP.set_ylabel("Phase (°)")

    # shared ticks/legend
    for g in range(n):
        for j in range(cols_per_grp):
            axR = AxR[g][j]; axP = AxP[g][j]
            # hide y labels on inner columns except leftmost
            if not (j == 0):
                axR.set_yticklabels([]); axP.set_yticklabels([])
            # hide x labels except bottom row (done above)
            if (g // ncols_groups) != (nrows - 1):
                axP.set_xticklabels([])
    if show_legend:
        labs = [c.upper() for c in comps]
        hs = [plt.Line2D([], [], color=cmap.get(c, "k"),
              lw=lw, ls=ls, marker=marker, ms=ms) for c in comps]
        fig.legend(
            hs, labs, loc="upper center", ncol=len(labs),
            frameon=False, bbox_to_anchor=(0.5, 1.02), fontsize=9,
        )
    return fig

# -------------------- measured vs predicted panels --------------------- #

def _pairs_meas_pred(
    sites: Any, pred_sites: Any
) -> List[Tuple[str, Any, Any]]:
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = ensure_sites(pred_sites, recursive=False, strict=False)
    m0 = { _name(ed, i): ed for i, ed in enumerate(_iter_items(S0)) }
    m1 = { _name(ed, i): ed for i, ed in enumerate(_iter_items(S1)) }
    out = []
    for st, ed in m0.items():
        if st in m1:
            out.append((st, ed, m1[st]))
    return sorted(out, key=lambda x: x[0])


def _nearest_idx(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    ix = np.searchsorted(x, y)
    ix = np.clip(ix, 1, x.size - 1)
    left = np.abs(y - x[ix - 1])
    right = np.abs(y - x[ix])
    pick_left = left <= right
    ix[pick_left] -= 1
    return ix


def _align_pred(fr_m: np.ndarray, fr_p: np.ndarray,
                z_p: np.ndarray) -> np.ndarray:
    j = _nearest_idx(fr_p, fr_m)
    return z_p[j]


def _rms_from(
    zm: np.ndarray,
    zp: np.ndarray,
    ze: Optional[np.ndarray],
    fr: np.ndarray,
    *,
    quantity: str,
) -> float:
    if quantity == "impedance":
        ym = np.log10(np.abs(zm))
        yp = np.log10(np.abs(zp))
        dy = ym - yp
        if ze is not None:
            ye = _err_log10_mag(zm, ze)
        else:
            ye = None
    else:
        rm = _rhoa_from(zm, fr)
        rp = _rhoa_from(zp, fr)
        dy = np.log10(rm) - np.log10(rp)
        if ze is not None:
            ye = _err_log10_rhoa(zm, ze)
        else:
            ye = None
    if ye is None:
        return float(np.sqrt(np.nanmean(dy * dy)))
    w = 1.0 / (np.square(ye) + 1e-12)
    return float(np.sqrt(np.nanmean(w * dy * dy)))

def plot_sites_fit_grid(
    sites: Any,
    pred_sites: Any,
    *,
    components: Tuple[str, ...] = ("xx", "xy", "yx", "yy"),
    quantity: str = "rhoa",      # rhoa|impedance
    x_axis: str = "period",      # period|frequency
    phase_range: Tuple[float, float] | None = (-180.0, 180.0),
    stations: Optional[List[str]] = None,
    ncols_groups: int = 2,       # station groups per row
    comp_wspace: float = 0.10,   # space between components
    group_hspace: float = 0.18,
    height_ratio: Tuple[int, int] = (2, 1),
    figsize_scale: Tuple[float, float] = (4.0, 3.0),
    colors_meas: Optional[Dict[str, str]] = None,
    color_fit_te: str = "#2ca02c",   # green
    color_fit_tm: str = "#e377c2",   # magenta
    marker: str = "o",
    ms: float = 2.5,
    lw: float = 1.4,
    ls_meas: str = "-",
    lw_fit: float = 2.0,
    ls_fit: str = "-",
    show_error_bars: bool = True,
    show_mode_legend: bool = True,
    title_group_fmt: str = "{station}",
    ylim_rhoa: Tuple[float, float] | None = None,
    ylim_phase: Tuple[float, float] | None = None,
    grid: bool = True,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    comps = tuple(c.lower() for c in components)
    cmap = {**_COL, **(colors_meas or {})}
    pairs = _pairs_meas_pred(sites, pred_sites)
    if stations:
        keep = set(stations)
        pairs = [p for p in pairs if p[0] in keep]
    if not pairs:
        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "no matching stations",
                ha="center", va="center")
        return fig
    n = len(pairs)
    ncols_groups = max(1, int(ncols_groups))
    nrows = (n + ncols_groups - 1) // ncols_groups
    W, H = figsize_scale
    fig = plt.figure(
        figsize=(W * ncols_groups, H * nrows),
        constrained_layout=False,
    )
    outer = fig.add_gridspec(
        nrows, ncols_groups, wspace=0.35, hspace=group_hspace
    )

    # build axes per group: 2 x ncomp
    AX = []
    for g in range(n):
        r, c = divmod(g, ncols_groups)
        gs = outer[r, c].subgridspec(
            2, len(comps),
            hspace=0.02, wspace=comp_wspace,
            height_ratios=height_ratio,
        )
        axR = [fig.add_subplot(gs[0, j]) for j in range(len(comps))]
        axP = [fig.add_subplot(gs[1, j], sharex=axR[0])
               for j in range(len(comps))]
        AX.append((axR, axP))

    # draw each station group
    for g, (st, edm, edp) in enumerate(pairs):
        axR, axP = AX[g]
        # station title centered over the row of components
        axR[0].set_title(title_group_fmt.format(station=st), pad=6)
        Zm = _zblk_flex(edm)
        if len(Zm) == 4:
            _, zm, frm, zem = Zm
        else:
            _, zm, frm = Zm[:3]; zem = None
        Zp = _zblk_flex(edp)
        if len(Zp) == 4:
            _, zp, frp, _ = Zp
        else:
            _, zp, frp = Zp[:3]; _ = None
        if zm is None or zp is None:
            for j in range(len(comps)):
                axR[j].axis("off"); axP[j].axis("off")
            continue
        x = (1.0 / frm) if x_axis == "period" else frm
        for j, comp in enumerate(comps):
            aR = axR[j]; aP = axP[j]
            aR.set_xscale("log"); aP.set_xscale("log")
            if grid:
                aR.grid(True, alpha=0.25, which="both")
                aP.grid(True, alpha=0.25, which="both")
            # measured
            zm_c = _comp_slice(zm, comp)
            ze_c = None
            if isinstance(zem, np.ndarray) and zem.shape == zm.shape:
                ze_c = _comp_slice(zem, comp)
            colm = cmap.get(comp, "k")
            if quantity == "impedance":
                ym = np.log10(np.abs(zm_c))
                yerr = _err_log10_mag(zm_c, ze_c)
                ylab = "Log10|Z|"
            else:
                rm = _rhoa_from(zm_c, frm)
                ym = np.log10(rm)
                yerr = _err_log10_rhoa(zm_c, ze_c)
                ylab = "Log10Rhoa (Ω·m)"
            plot_errorbar(
                aR, x, ym, y_err=yerr, color=colm,
                marker=marker, ms=ms, ls=ls_meas, lw=lw,
                show_error_bars=show_error_bars,
            )
            phm = _phase_deg(zm_c)
            if phase_range is not None:
                phm = _wrap_phase(phm, phase_range)
            plot_errorbar(
                aP, x, phm, y_err=_err_phase_deg(zm_c, ze_c),
                color=colm, marker=marker, ms=ms,
                ls=ls_meas, lw=lw, show_error_bars=False,
            )
            # predicted (aligned to measured freq)
            zp_c = _comp_slice(zp, comp)
            zp_c = _align_pred(frm, frp, zp_c)
            colf = color_fit_te if comp in ("xx", "xy") \
                   else color_fit_tm
            if quantity == "impedance":
                yp = np.log10(np.abs(zp_c))
            else:
                rp = _rhoa_from(zp_c, frm)
                yp = np.log10(rp)
            aR.plot(x, yp, ls_fit, color=colf, lw=lw_fit)
            php = _phase_deg(zp_c)
            if phase_range is not None:
                php = _wrap_phase(php, phase_range)
            aP.plot(x, php, ls_fit, color=colf, lw=lw_fit)
            # per-component header + RMS
            rmsc = _rms_from(zm_c, zp_c, ze_c, frm, quantity=quantity)
            aR.text(
                0.50, 1.02, f"Z{comp.upper()}",
                ha="center", va="bottom", transform=aR.transAxes,
                fontsize=9,
            )
            aR.text(
                0.02, 1.02, f"rms = {rmsc:.2f}",
                ha="left", va="bottom", transform=aR.transAxes,
                fontsize=8,
            )
            # axes cosmetics
            if ylim_rhoa:
                aR.set_ylim(*ylim_rhoa)
            if ylim_phase:
                aP.set_ylim(*ylim_phase)
            elif phase_range is not None:
                aP.set_ylim(*phase_range)
            if j == 0:
                aR.set_ylabel(ylab)
                aP.set_ylabel("Phase (°)")
            else:
                aR.set_yticklabels([]); aP.set_yticklabels([])
            if (g // ncols_groups) == (nrows - 1):
                aP.set_xlabel(
                    "Log10Period (s)" if x_axis == "period"
                    else "Log10Freq (Hz)"
                )
            else:
                aP.set_xlabel("")
    # global legend (TE/TM fits)
    if show_mode_legend:
        h_te = plt.Line2D([], [], color=color_fit_te, lw=lw_fit)
        h_tm = plt.Line2D([], [], color=color_fit_tm, lw=lw_fit)
        fig.legend(
            [h_te, h_tm],
            ["TE fit: xx/xy", "TM fit: yx/yy"],
            loc="upper center", ncol=2, frameon=False,
            bbox_to_anchor=(0.5, 1.02), fontsize=9,
        )
    return fig

