# pycsamt/emtools/plot.py
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple
import numpy as np
import matplotlib.pyplot as plt

from ._core import (
    ensure_sites,
    _iter_items,
    _get_t_block,
    _get_z_block,
    _name,
)

from ..utils.plot import plot_errorbar
from ..api._rose_style import _UNSET
from ..api.control import PYCSAMT_CONTROL, wrap_phase
from ..api.style import PYCSAMT_STYLE


# ------------------------------- helpers -------------------------------- #

def _style_col() -> dict:
    """Return component → color dict from the current PYCSAMT_STYLE.mt."""
    mt = PYCSAMT_STYLE.mt
    return {
        "xy": mt.xy.color,
        "yx": mt.yx.color,
        "xx": mt.xx.color,
        "yy": mt.yy.color,
    }


# legacy alias kept so that any direct caller outside this module still works
_COL = {
    "xy": "#1f77b4",
    "yx": "#d62728",
    "xx": "#2ca02c",
    "yy": "#9467bd",
}

_IDX = {"x": 0, "y": 1}

def _comp_slice(z: np.ndarray, comp: str) -> np.ndarray:
    a, b = _IDX[comp[0]], _IDX[comp[1]]
    return z[:, a, b]

def _phase_deg(z: np.ndarray) -> np.ndarray:
    return np.degrees(np.angle(z))

def _wrap_phase(ph: np.ndarray, rng: Tuple[float, float]) -> np.ndarray:
    return wrap_phase(ph, rng)

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


def _err_rhoa(z: np.ndarray, ze: Optional[np.ndarray],
              fr: np.ndarray) -> np.ndarray | None:
    """Return approximate apparent-resistivity uncertainty."""
    if ze is None:
        return None
    rho = _rhoa_from(z, fr)
    s = np.abs(ze)
    m = np.abs(z) + 1e-24
    return rho * 2.0 * s / m

def _pick(axs, i, j):
    return axs[i * 2 + j]


def _zblk_flex(ed: Any):
    try:
        return _get_z_block(ed, with_errors=True)
    except TypeError:
        return _get_z_block(ed)


def _control_or_default(control):
    """Return the active plotting-view control."""
    return control or PYCSAMT_CONTROL


def _site_items(
    sites: Any,
    *,
    recursive: bool,
    on_dup: str,
    strict: bool,
    verbose: int,
) -> list[Any]:
    """Return site items, falling back to raw iterables when needed."""
    try:
        normalized = ensure_sites(
            sites,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
        items = list(_iter_items(normalized))
        if items:
            return items
    except Exception:
        if strict:
            raise
    return list(_iter_items(sites))


def _x_from_freq(fr: np.ndarray, control) -> tuple[np.ndarray, str, bool]:
    """Return x values, label, and log-scale flag from frequencies."""
    xctrl = _control_or_default(control).x
    return xctrl.transform(fr), xctrl.label(), xctrl.use_log_scale()


def _rho_display(
    z: np.ndarray,
    ze: Optional[np.ndarray],
    fr: np.ndarray,
    control,
) -> tuple[np.ndarray, np.ndarray | None, str]:
    """Return apparent resistivity, errors, and axis label."""
    c = _control_or_default(control)
    rho = _rhoa_from(z, fr)
    rho_err = _err_rhoa(z, ze, fr)
    return c.rho.transform(rho), c.rho.error(rho, rho_err), c.rho.label()


def _phase_display(
    z: np.ndarray,
    ze: Optional[np.ndarray],
    control,
    phase_range: tuple[float, float] | None,
) -> tuple[np.ndarray, np.ndarray | None, str, tuple[float, float] | None]:
    """Return phase, errors, label, and display range."""
    c = _control_or_default(control)
    ph = _phase_deg(z)
    pc = c.phase
    if phase_range is not None:
        ph = wrap_phase(ph, phase_range)
        display_range = phase_range
    else:
        ph = pc.transform(ph)
        display_range = pc.range if pc.wrap else None
    err = _err_phase_deg(z, ze)
    if pc.unit.lower() == "radian":
        if err is not None:
            err = np.deg2rad(err)
        if display_range is not None:
            display_range = tuple(np.deg2rad(display_range))
    return ph, err, pc.label(), display_range


def _component_style(comp: str, raw: bool, force_style: bool):
    """Return style object for a component or raw data."""
    if raw and not force_style:
        return PYCSAMT_STYLE.raw
    return PYCSAMT_STYLE.mt.component(comp)


def _errorbar_from_style(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    yerr: np.ndarray | None,
    style,
    *,
    color: str | None = None,
    show_error_bars: bool = True,
):
    """Draw a styled errorbar curve."""
    kwargs = style.errorbar_kwargs() if hasattr(style, "errorbar_kwargs") else {}
    if color is not None:
        kwargs["color"] = color
        kwargs["ecolor"] = color
        kwargs["mec"] = color
    marker = kwargs.pop("marker", "o")
    ms = kwargs.pop("ms", kwargs.pop("markersize", 4.0))
    mfc = kwargs.pop("mfc", "white")
    mew = kwargs.pop("mew", 1.0)
    ls = kwargs.pop("ls", "-")
    lw = kwargs.pop("lw", 1.0)
    ecolor = kwargs.pop("ecolor", kwargs.get("color", "black"))
    capsize = kwargs.pop("capsize", 2.0)
    elinewidth = kwargs.pop("elinewidth", lw)
    label = kwargs.pop("label", None)
    if not show_error_bars:
        yerr = None
    return ax.errorbar(
        x,
        y,
        yerr=yerr,
        color=kwargs.pop("color", "black"),
        marker=marker,
        ms=ms,
        mfc=mfc,
        mec=ecolor,
        mew=mew,
        ls=ls,
        lw=lw,
        ecolor=ecolor,
        capsize=capsize,
        elinewidth=elinewidth,
        label=label,
        **kwargs,
    )


def _axes_group_bounds(axes: list[plt.Axes]) -> tuple[float, float, float, float]:
    """Return figure-coordinate bounds for a group of axes."""
    boxes = [ax.get_position() for ax in axes]
    left = min(box.x0 for box in boxes)
    right = max(box.x1 for box in boxes)
    bottom = min(box.y0 for box in boxes)
    top = max(box.y1 for box in boxes)
    return left, right, bottom, top


def _add_raw_group_labels(
    fig: plt.Figure,
    ax_rho: list[plt.Axes],
    ax_phase: list[plt.Axes],
    *,
    station: str,
    rho_label: str,
    phase_label: str,
    x_label: str,
    show_x_label: bool,
    x_label_pad: float = 0.078,
) -> None:
    """Draw shared labels for one raw 1-D station group."""
    left, right, bottom, top = _axes_group_bounds(ax_rho + ax_phase)
    xmid = 0.5 * (left + right)
    rho_y = 0.5 * (ax_rho[0].get_position().y0 + ax_rho[0].get_position().y1)
    phase_y = 0.5 * (
        ax_phase[0].get_position().y0 + ax_phase[0].get_position().y1
    )
    fig.text(
        xmid,
        min(top + 0.024, 0.985),
        station,
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )
    fig.text(
        max(left - 0.026, 0.006),
        rho_y,
        rho_label,
        ha="right",
        va="center",
        rotation=90,
        fontsize=8,
    )
    fig.text(
        max(left - 0.026, 0.006),
        phase_y,
        phase_label,
        ha="right",
        va="center",
        rotation=90,
        fontsize=8,
    )
    if show_x_label:
        fig.text(
            xmid,
            max(bottom - x_label_pad, 0.006),
            x_label,
            ha="center",
            va="top",
            fontsize=8,
        )


def _apply_raw_tick_style(
    ax: plt.Axes,
    *,
    show_x: bool,
    rotation: float,
    fontsize: int,
) -> None:
    """Apply compact tick-label styling to raw 1-D panels."""
    ax.tick_params(axis="both", labelsize=fontsize)
    if show_x:
        ax.tick_params(axis="x", labelbottom=True)
        for label in ax.get_xticklabels():
            label.set_rotation(rotation)
            label.set_ha("center" if abs(rotation) >= 80 else "right")
            label.set_fontsize(fontsize)
    else:
        ax.tick_params(axis="x", labelbottom=False)


def _raw_label_mode(
    label_mode: str | None,
    shared_group_labels: bool,
) -> str:
    """Resolve the raw-panel label mode."""
    if label_mode is None:
        return "shared" if shared_group_labels else "axis"
    mode = str(label_mode).lower()
    if mode not in {"shared", "axis"}:
        msg = "label_mode must be 'shared' or 'axis'."
        raise ValueError(msg)
    return mode


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
    colors: Optional[Dict[str, str]] = None,  # None → from PYCSAMT_STYLE.mt
    marker         = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.marker
    ms             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.ms
    lw             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.lw
    ls             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.ls
    show_error_bars: bool = True,
    show_legend: bool = False,
    title_fmt: str = "{station}",
    ylim_rhoa: Tuple[float, float] | None = None,
    ylim_phase: Tuple[float, float] | None = None,
    grid: bool = True,
    preserve_duplicates: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    # ── resolve visual style from PYCSAMT_STYLE.mt ───────────────────────
    _mt = PYCSAMT_STYLE.mt
    if marker is _UNSET: marker = _mt.xy.marker
    if ms     is _UNSET: ms     = _mt.xy.ms
    if lw     is _UNSET: lw     = _mt.xy.lw
    if ls     is _UNSET: ls     = _mt.xy.ls

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    comps = tuple(c.lower() for c in components)
    cmap = {**_style_col(), **(colors or {})}
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
        Z, z, fr, ze = _get_z_block(ed, with_errors=True)
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


def plot_raw_sites_1d(
    sites: Any,
    *,
    stations: Optional[List[str]] = None,
    components: Tuple[str, ...] = ("xx", "xy", "yx", "yy"),
    raw: bool = True,
    force_style: bool = False,
    control: Any | None = None,
    phase_range: Tuple[float, float] | None = None,
    ncols_groups: int = 3,
    comp_wspace: float = 0.12,
    group_hspace: float = 0.25,
    height_ratio: Tuple[int, int] = (2, 1),
    figsize_scale: Tuple[float, float] = (4.2, 3.1),
    colors: Optional[Dict[str, str]] = None,
    title_group_fmt: str = "{station}",
    title_comp_fmt: str = "Z{component}",
    shared_group_labels: bool = True,
    label_mode: str | None = None,
    shared_x_label_pad: float = 0.078,
    x_tick_rotation: float | None = None,
    tick_fontsize: int = 7,
    show_error_bars: bool = True,
    show_component_legend: bool = True,
    ylim_rhoa: Tuple[float, float] | None = None,
    ylim_phase: Tuple[float, float] | None = None,
    grid: bool = True,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    """Plot raw or processed 1-D rho/phase panels by station.

    The layout mirrors diagnostic raw-data figures used in AMT/MT
    workflows: every selected station is a group, every component is a
    column, and each component column contains apparent resistivity above
    phase.  When ``raw=True`` the raw-data style from
    :data:`pycsamt.api.style.PYCSAMT_STYLE.raw` is used automatically,
    producing black diagnostic traces unless ``force_style=True`` or
    explicit ``colors`` are provided.
    """
    ctl = _control_or_default(control)
    comps = tuple(c.lower() for c in components)
    resolved_label_mode = _raw_label_mode(label_mode, shared_group_labels)
    use_shared_labels = resolved_label_mode == "shared"
    if x_tick_rotation is None:
        x_tick_rotation = 45.0 if use_shared_labels else 90.0
    site_items = _site_items(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    items = []
    keep = set(stations or [])
    for i, ed in enumerate(site_items):
        st = _name(ed, i)
        if keep and st not in keep:
            continue
        items.append((st, ed))
    if not items:
        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "no stations", ha="center", va="center")
        return fig

    n = len(items)
    ncols_groups = max(1, int(ncols_groups))
    nrows = (n + ncols_groups - 1) // ncols_groups
    W, H = figsize_scale
    fig = plt.figure(
        figsize=(W * ncols_groups, H * nrows),
        constrained_layout=False,
    )
    outer = fig.add_gridspec(
        nrows,
        ncols_groups,
        wspace=0.35,
        hspace=group_hspace,
    )

    axes: list[tuple[list[plt.Axes], list[plt.Axes]]] = []
    for g in range(n):
        r, c = divmod(g, ncols_groups)
        gs = outer[r, c].subgridspec(
            2,
            len(comps),
            hspace=0.02,
            wspace=comp_wspace,
            height_ratios=height_ratio,
        )
        ax_rho = [fig.add_subplot(gs[0, j]) for j in range(len(comps))]
        ax_phase = [
            fig.add_subplot(gs[1, j], sharex=ax_rho[j])
            for j in range(len(comps))
        ]
        axes.append((ax_rho, ax_phase))

    for g, (station, ed) in enumerate(items):
        ax_rho, ax_phase = axes[g]
        station_label = title_group_fmt.format(station=station)
        if not use_shared_labels:
            ax_rho[0].set_title(station_label, pad=6)
        out = _zblk_flex(ed)
        if len(out) == 4:
            _, z, fr, ze = out
        else:
            _, z, fr = out[:3]
            ze = None
        if z is None or fr is None:
            for ax in ax_rho + ax_phase:
                ax.axis("off")
            continue
        fr = np.asarray(fr, dtype=float)
        x, x_label, log_x = _x_from_freq(fr, ctl)
        group_rho_label = ""
        group_phase_label = ""

        for j, comp in enumerate(comps):
            a_rho = ax_rho[j]
            a_phase = ax_phase[j]
            a_rho.text(
                0.5,
                0.96,
                title_comp_fmt.format(component=comp.upper()),
                ha="center",
                va="top",
                transform=a_rho.transAxes,
                fontsize=8,
                bbox={
                    "boxstyle": "round,pad=0.12",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.75,
                },
            )
            if log_x:
                a_rho.set_xscale("log")
                a_phase.set_xscale("log")
            if grid:
                a_rho.grid(True, alpha=0.25, which="both")
                a_phase.grid(True, alpha=0.25, which="both")

            zz = _comp_slice(z, comp)
            ee = None
            if isinstance(ze, np.ndarray) and ze.shape == z.shape:
                ee = _comp_slice(ze, comp)
            rho_y, rho_err, rho_label = _rho_display(zz, ee, fr, ctl)
            phase_y, phase_err, phase_label, display_range = _phase_display(
                zz,
                ee,
                ctl,
                phase_range,
            )
            group_rho_label = rho_label
            group_phase_label = phase_label
            style = _component_style(comp, raw, force_style)
            color = None
            if colors is not None:
                color = colors.get(comp)
            _errorbar_from_style(
                a_rho,
                x,
                rho_y,
                rho_err,
                style,
                color=color,
                show_error_bars=show_error_bars,
            )
            _errorbar_from_style(
                a_phase,
                x,
                phase_y,
                phase_err,
                style,
                color=color,
                show_error_bars=False,
            )
            if ylim_rhoa is not None:
                a_rho.set_ylim(*ylim_rhoa)
            if ylim_phase is not None:
                a_phase.set_ylim(*ylim_phase)
            elif display_range is not None:
                a_phase.set_ylim(*display_range)
            is_bottom_group = (g // ncols_groups) == (nrows - 1)
            if j == 0 and not use_shared_labels:
                a_rho.set_ylabel(rho_label)
                a_phase.set_ylabel(phase_label)
            else:
                a_rho.set_ylabel("")
                a_phase.set_ylabel("")
                if j != 0:
                    a_rho.set_yticklabels([])
                    a_phase.set_yticklabels([])
            if is_bottom_group and not use_shared_labels:
                a_phase.set_xlabel(x_label)
            else:
                a_phase.set_xlabel("")
            _apply_raw_tick_style(
                a_phase,
                show_x=is_bottom_group,
                rotation=float(x_tick_rotation),
                fontsize=int(tick_fontsize),
            )
            _apply_raw_tick_style(
                a_rho,
                show_x=False,
                rotation=float(x_tick_rotation),
                fontsize=int(tick_fontsize),
            )
        if use_shared_labels:
            _add_raw_group_labels(
                fig,
                ax_rho,
                ax_phase,
                station=station_label,
                rho_label=group_rho_label,
                phase_label=group_phase_label,
                x_label=x_label,
                show_x_label=(g // ncols_groups) == (nrows - 1),
                x_label_pad=shared_x_label_pad,
            )

    if show_component_legend:
        handles = []
        labels = []
        for comp in comps:
            style = _component_style(comp, raw, force_style)
            color = colors.get(comp) if colors else None
            kwargs = style.plot_kwargs() if hasattr(style, "plot_kwargs") else {}
            if color is not None:
                kwargs["color"] = color
            handles.append(plt.Line2D([], [], **kwargs))
            labels.append(comp.upper())
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=len(labels),
            frameon=False,
            bbox_to_anchor=(0.5, -0.015),
            fontsize=8,
        )
    return fig


def _add_response_tipper_group_labels(
    fig: plt.Figure,
    axes: Sequence[plt.Axes],
    row_axes: dict[str, plt.Axes],
    *,
    station: str,
    rho_label: str,
    phase_label: str,
    x_label: str,
    show_x_label: bool,
    x_label_pad: float = 0.070,
) -> None:
    """Draw shared labels for one response/tipper station group."""
    left, right, bottom, top = _axes_group_bounds(list(axes))
    xmid = 0.5 * (left + right)

    fig.text(
        xmid,
        min(top + 0.020, 0.985),
        station,
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )
    for key, label in (
        ("rho", rho_label),
        ("phase", phase_label),
        ("tx", r"$T_x$"),
        ("ty", r"$T_y$"),
    ):
        ax = row_axes.get(key)
        if ax is None or not label:
            continue
        box = ax.get_position()
        fig.text(
            max(left - 0.046, 0.006),
            0.5 * (box.y0 + box.y1),
            label,
            ha="right",
            va="center",
            rotation=90,
            fontsize=8,
        )
    if show_x_label:
        fig.text(
            xmid,
            max(bottom - x_label_pad, 0.006),
            x_label,
            ha="center",
            va="top",
            fontsize=8,
        )


def _plot_tipper_part(
    ax: plt.Axes,
    x: np.ndarray,
    values: np.ndarray,
    errors: np.ndarray | None,
    *,
    show_error_bars: bool,
    line_style: str | None = None,
    real_label: str = "real",
    imag_label: str = "imag",
) -> None:
    """Plot real and imaginary tipper parts on one compact axis."""
    real_style = PYCSAMT_STYLE.mt.xy
    imag_style = PYCSAMT_STYLE.mt.yx
    if line_style is not None:
        real_style = real_style.copy(ls=line_style)
        imag_style = imag_style.copy(ls=line_style)
    _errorbar_from_style(
        ax,
        x,
        np.real(values),
        errors,
        real_style,
        show_error_bars=show_error_bars,
    )
    imag_kw = imag_style.copy(marker="x", mfc=imag_style.color)
    _errorbar_from_style(
        ax,
        x,
        np.imag(values),
        errors,
        imag_kw,
        show_error_bars=show_error_bars,
    )
    ax.plot([], [], **real_style.plot_kwargs(label=real_label))
    ax.plot([], [], **imag_kw.plot_kwargs(label=imag_label))


def plot_response_tipper(
    sites: Any,
    *,
    stations: Optional[List[str]] = None,
    components: Tuple[str, ...] = ("xy", "yx"),
    tipper_components: Tuple[str, ...] = ("tx", "ty"),
    raw: bool = False,
    force_style: bool = False,
    control: Any | None = None,
    phase_range: Tuple[float, float] | None = None,
    ncols_groups: int = 3,
    comp_wspace: float = 0.12,
    group_hspace: float = 0.32,
    height_ratios: Tuple[float, ...] = (2.2, 1.1, 0.75, 0.75),
    figsize_scale: Tuple[float, float] = (4.8, 4.6),
    colors: Optional[Dict[str, str]] = None,
    tipper_span_group: bool = False,
    line_style: str | None = None,
    tipper_line_style: str | None = ":",
    label_component_x: bool = True,
    label_tipper_x: bool = True,
    title_group_fmt: str = "{station}",
    title_comp_fmt: str = "Z{component}",
    shared_group_labels: bool = True,
    shared_x_label_pad: float = 0.074,
    x_tick_rotation: float = 0.0,
    tick_fontsize: int = 7,
    show_error_bars: bool = True,
    show_tipper_error_bars: bool = False,
    show_component_legend: bool = True,
    show_tipper_legend: bool = True,
    ylim_rhoa: Tuple[float, float] | None = None,
    ylim_phase: Tuple[float, float] | None = None,
    ylim_tipper: Tuple[float, float] | None = (-0.6, 0.6),
    grid: bool = True,
    preserve_duplicates: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Plot impedance response panels with station-level tipper rows.

    The figure is designed for MT/AMT quality control where impedance and
    tipper behaviour must be inspected together.  Each station is a group:
    apparent resistivity and phase are shown per impedance component, while
    compact :math:`T_x` and :math:`T_y` panels span the full group width.

    Parameters
    ----------
    sites : Sites-like
        EDI path, collection of EDI files, ``Sites`` object, or iterable
        accepted by :func:`pycsamt.emtools._core.ensure_sites`.
    stations : list of str, optional
        Station names to display.  When omitted, all stations are used.
    components : tuple of str, default ("xy", "yx")
        Impedance tensor components plotted in the resistivity and phase
        rows.  Values are component keys such as ``"xy"``, ``"yx"``,
        ``"xx"``, or ``"yy"``.
    tipper_components : tuple of str, default ("tx", "ty")
        Tipper components to draw.  The usual complete diagnostic uses both
        ``"tx"`` and ``"ty"``.
    raw : bool, default False
        If ``True``, impedance response curves use the package raw-data
        style unless ``force_style`` is also true.  Tipper real/imaginary
        curves keep distinct package component colours for readability.
    force_style : bool, default False
        Use component colours even when ``raw=True``.
    control : object, optional
        Plot view control.  Defaults to
        :data:`pycsamt.api.control.PYCSAMT_CONTROL`.
    phase_range : tuple of float or None, optional
        Explicit phase display range.  If omitted, the active
        ``control.phase`` policy is used.
    ncols_groups : int, default 3
        Number of station groups per figure row.
    comp_wspace, group_hspace : float
        Spacing inside and between station groups.
    height_ratios : tuple of float, default (2.2, 1.1, 0.75, 0.75)
        Relative heights for rho, phase, Tx, and Ty rows.  If only one
        tipper component is requested, the unused extra ratio is ignored.
    figsize_scale : tuple of float, default (4.8, 4.6)
        Width and height multiplier for each station group row/column.
    colors : dict, optional
        Optional impedance component colour overrides.
    tipper_span_group : bool, default False
        If ``True``, each tipper row spans all impedance component columns
        in a station group.  If ``False``, Tx and Ty are repeated under
        each component column, giving a compact grid such as
        ``rho/phase/Tx/Ty`` for every component.
    line_style : str or None, optional
        Optional line style override for impedance curves.
    tipper_line_style : str or None, default ":"
        Optional line style override for real and imaginary tipper curves.
    label_component_x, label_tipper_x : bool, default True
        Put the active x-axis label under each impedance-component stack
        and under the bottom tipper row.
    title_group_fmt, title_comp_fmt : str
        Format strings for station and component labels.
    shared_group_labels : bool, default True
        Use group-level rho/phase/tipper labels instead of repeating labels
        on every small axis.
    shared_x_label_pad : float, default 0.074
        Figure-coordinate padding used for shared x labels.
    x_tick_rotation : float, default 0
        Rotation for bottom x tick labels.
    tick_fontsize : int, default 7
        Tick-label size for compact panels.
    show_error_bars, show_tipper_error_bars : bool
        Toggle impedance and tipper error bars independently.
    show_component_legend, show_tipper_legend : bool
        Toggle global legends.
    ylim_rhoa, ylim_phase, ylim_tipper : tuple or None
        Optional y-axis limits.
    grid : bool, default True
        Draw light panel grids.
    preserve_duplicates : bool, default False
        Preserve repeated in-memory EDI objects instead of normalizing them
        through :func:`ensure_sites`.  This is useful for diagnostic demos or
        before/after comparisons where two display stations intentionally
        share the same source station name.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    matplotlib.figure.Figure
        The assembled response/tipper figure.

    Examples
    --------
    >>> from pycsamt.emtools.plot import plot_response_tipper
    >>> fig = plot_response_tipper(
    ...     "data/AMT/TIPPER/HBH03_IMP.edi",
    ...     components=("xy", "yx"),
    ... )
    """
    ctl = _control_or_default(control)
    comps = tuple(c.lower() for c in components)
    tips = tuple(t.lower() for t in tipper_components)
    for comp in comps:
        if comp not in {"xx", "xy", "yx", "yy"}:
            msg = f"unknown impedance component {comp!r}."
            raise ValueError(msg)
    for tip in tips:
        if tip not in {"tx", "ty"}:
            msg = "tipper_components must contain only 'tx' and/or 'ty'."
            raise ValueError(msg)

    if (
        preserve_duplicates
        and isinstance(sites, (list, tuple))
        and not all(isinstance(item, (str, bytes)) for item in sites)
    ):
        site_items = list(sites)
    else:
        site_items = _site_items(
            sites,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
    keep = set(stations or [])
    items: list[tuple[str, Any]] = []
    seen_names: dict[str, int] = {}
    for i, ed in enumerate(site_items):
        base_st = _name(ed, i)
        count = seen_names.get(base_st, 0) + 1
        seen_names[base_st] = count
        st = base_st if count == 1 else f"{base_st}_{count}"
        if keep and st not in keep:
            continue
        items.append((st, ed))
    if not items:
        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "no stations", ha="center", va="center")
        return fig

    n = len(items)
    ncols_groups = max(1, int(ncols_groups))
    nrows = (n + ncols_groups - 1) // ncols_groups
    W, H = figsize_scale
    fig = plt.figure(
        figsize=(W * ncols_groups, H * nrows),
        constrained_layout=False,
    )
    outer = fig.add_gridspec(
        nrows,
        ncols_groups,
        wspace=0.35,
        hspace=group_hspace,
    )

    row_names = ("rho", "phase", *tips)
    ratios = tuple(height_ratios[: len(row_names)])
    if len(ratios) < len(row_names):
        ratios = ratios + tuple([0.75] * (len(row_names) - len(ratios)))

    groups: list[dict[str, Any]] = []
    for g in range(n):
        r, c = divmod(g, ncols_groups)
        gs = outer[r, c].subgridspec(
            len(row_names),
            len(comps),
            hspace=0.28,
            wspace=comp_wspace,
            height_ratios=ratios,
        )
        ax_rho = [fig.add_subplot(gs[0, j]) for j in range(len(comps))]
        ax_phase = [
            fig.add_subplot(gs[1, j], sharex=ax_rho[j])
            for j in range(len(comps))
        ]
        tip_axes: dict[str, plt.Axes | list[plt.Axes]] = {}
        for row, tip in enumerate(tips, start=2):
            if tipper_span_group:
                tip_axes[tip] = fig.add_subplot(gs[row, :], sharex=ax_rho[0])
            else:
                tip_axes[tip] = [
                    fig.add_subplot(gs[row, j], sharex=ax_rho[j])
                    for j in range(len(comps))
                ]
        groups.append(
            {"rho": ax_rho, "phase": ax_phase, "tip": tip_axes}
        )

    for g, (station, ed) in enumerate(items):
        group = groups[g]
        ax_rho = group["rho"]
        ax_phase = group["phase"]
        tip_axes = group["tip"]
        station_label = title_group_fmt.format(station=station)

        z_out = _zblk_flex(ed)
        if len(z_out) == 4:
            _, z, fr, ze = z_out
        else:
            _, z, fr = z_out[:3]
            ze = None
        if z is None or fr is None:
            for ax in ax_rho + ax_phase + list(tip_axes.values()):
                ax.axis("off")
            continue
        fr = np.asarray(fr, dtype=float)
        x, x_label, log_x = _x_from_freq(fr, ctl)
        group_rho_label = ""
        group_phase_label = ""
        is_bottom_group = (g // ncols_groups) == (nrows - 1)

        for j, comp in enumerate(comps):
            a_rho = ax_rho[j]
            a_phase = ax_phase[j]
            a_rho.text(
                0.5,
                0.96,
                title_comp_fmt.format(component=comp.upper()),
                ha="center",
                va="top",
                transform=a_rho.transAxes,
                fontsize=8,
                bbox={
                    "boxstyle": "round,pad=0.12",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.75,
                },
            )
            for ax in (a_rho, a_phase):
                if log_x:
                    ax.set_xscale("log")
                if grid:
                    ax.grid(True, alpha=0.25, which="both")
            zz = _comp_slice(z, comp)
            ee = None
            if isinstance(ze, np.ndarray) and ze.shape == z.shape:
                ee = _comp_slice(ze, comp)
            rho_y, rho_err, rho_label = _rho_display(zz, ee, fr, ctl)
            phase_y, phase_err, phase_label, display_range = _phase_display(
                zz,
                ee,
                ctl,
                phase_range,
            )
            group_rho_label = rho_label
            group_phase_label = phase_label
            style = _component_style(comp, raw, force_style)
            if line_style is not None and hasattr(style, "copy"):
                style = style.copy(ls=line_style)
            color = colors.get(comp) if colors else None
            _errorbar_from_style(
                a_rho,
                x,
                rho_y,
                rho_err,
                style,
                color=color,
                show_error_bars=show_error_bars,
            )
            _errorbar_from_style(
                a_phase,
                x,
                phase_y,
                phase_err,
                style,
                color=color,
                show_error_bars=False,
            )
            if ylim_rhoa is not None:
                a_rho.set_ylim(*ylim_rhoa)
            if ylim_phase is not None:
                a_phase.set_ylim(*ylim_phase)
            elif display_range is not None:
                a_phase.set_ylim(*display_range)
            if not shared_group_labels and j == 0:
                a_rho.set_ylabel(rho_label)
                a_phase.set_ylabel(phase_label)
            elif j != 0:
                a_rho.set_yticklabels([])
                a_phase.set_yticklabels([])
            a_rho.set_xlabel("")
            a_phase.set_xlabel(x_label if label_component_x else "")
            _apply_raw_tick_style(
                a_rho,
                show_x=False,
                rotation=x_tick_rotation,
                fontsize=tick_fontsize,
            )
            _apply_raw_tick_style(
                a_phase,
                show_x=bool(label_component_x),
                rotation=x_tick_rotation,
                fontsize=tick_fontsize,
            )

        tipper_x_label = None
        t_out = _get_t_block(ed, with_errors=True)
        if len(t_out) == 4:
            _, tipper, tfr, terr = t_out
        else:
            _, tipper, tfr = t_out[:3]
            terr = None
        if tipper is None or tfr is None:
            for tip, axes_obj in tip_axes.items():
                axes_list = axes_obj if isinstance(axes_obj, list) else [axes_obj]
                for ax in axes_list:
                    ax.text(
                        0.5,
                        0.5,
                        "no tipper",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                        fontsize=8,
                        color="0.4",
                    )
                    if grid:
                        ax.grid(True, alpha=0.25, which="both")
        else:
            tx, t_label, t_log_x = _x_from_freq(np.asarray(tfr, float), ctl)
            tipper = np.asarray(tipper, complex)
            terr_arr = terr if isinstance(terr, np.ndarray) else None
            for tip in tips:
                axes_obj = tip_axes[tip]
                axes_list = axes_obj if isinstance(axes_obj, list) else [axes_obj]
                col = 0 if tip == "tx" else 1
                err = None
                if terr_arr is not None and terr_arr.shape == tipper.shape:
                    err = np.abs(terr_arr[:, col])
                for ax in axes_list:
                    if t_log_x:
                        ax.set_xscale("log")
                    if grid:
                        ax.grid(True, alpha=0.25, which="both")
                    _plot_tipper_part(
                        ax,
                        tx,
                        tipper[:, col],
                        err,
                        show_error_bars=show_tipper_error_bars,
                        line_style=tipper_line_style,
                    )
                    if ylim_tipper is not None:
                        ax.set_ylim(*ylim_tipper)
                    ax.axhline(0.0, color="0.75", lw=0.7, zorder=0)
                tipper_x_label = t_label
            if tipper_x_label is not None:
                x_label = tipper_x_label

        flat_tip_axes: list[plt.Axes] = []
        for axes_obj in tip_axes.values():
            if isinstance(axes_obj, list):
                flat_tip_axes.extend(axes_obj)
            else:
                flat_tip_axes.append(axes_obj)
        all_axes = ax_rho + ax_phase + flat_tip_axes
        for ax in all_axes:
            _apply_raw_tick_style(
                ax,
                show_x=False,
                rotation=x_tick_rotation,
                fontsize=tick_fontsize,
            )
        if label_component_x:
            for ax in ax_phase:
                _apply_raw_tick_style(
                    ax,
                    show_x=True,
                    rotation=x_tick_rotation,
                    fontsize=tick_fontsize,
                )
                ax.set_xlabel(x_label, labelpad=2)
        bottom_axes_obj = list(tip_axes.values())[-1] if tip_axes else ax_phase
        bottom_axes = (
            bottom_axes_obj if isinstance(bottom_axes_obj, list)
            else [bottom_axes_obj]
        )
        for bottom_ax in bottom_axes:
            _apply_raw_tick_style(
                bottom_ax,
                show_x=bool(label_tipper_x),
                rotation=x_tick_rotation,
                fontsize=tick_fontsize,
            )
            bottom_ax.set_xlabel(x_label if label_tipper_x else "")
        if shared_group_labels:
            row_axes = {
                "rho": ax_rho[0],
                "phase": ax_phase[0],
            }
            if "tx" in tip_axes:
                tx_axes = tip_axes["tx"]
                row_axes["tx"] = tx_axes[0] if isinstance(tx_axes, list) else tx_axes
            if "ty" in tip_axes:
                ty_axes = tip_axes["ty"]
                row_axes["ty"] = ty_axes[0] if isinstance(ty_axes, list) else ty_axes
            _add_response_tipper_group_labels(
                fig,
                all_axes,
                row_axes,
                station=station_label,
                rho_label=group_rho_label,
                phase_label=group_phase_label,
                x_label=x_label,
                show_x_label=(
                    is_bottom_group
                    and not (label_component_x or label_tipper_x)
                ),
                x_label_pad=shared_x_label_pad,
            )

    if show_component_legend:
        handles = []
        labels = []
        for comp in comps:
            style = _component_style(comp, raw, force_style)
            color = colors.get(comp) if colors else None
            kwargs = style.plot_kwargs() if hasattr(style, "plot_kwargs") else {}
            if color is not None:
                kwargs["color"] = color
            handles.append(plt.Line2D([], [], **kwargs))
            labels.append(f"Z{comp.upper()}")
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=max(1, len(labels)),
            frameon=False,
            bbox_to_anchor=(0.38, -0.015),
            fontsize=8,
        )
    if show_tipper_legend:
        real_style = PYCSAMT_STYLE.mt.xy
        imag_style = PYCSAMT_STYLE.mt.yx.copy(
            ls="--",
            marker="x",
            mfc=PYCSAMT_STYLE.mt.yx.color,
        )
        handles = [
            plt.Line2D([], [], **real_style.plot_kwargs(label="real")),
            plt.Line2D([], [], **imag_style.plot_kwargs(label="imag")),
        ]
        fig.legend(
            handles,
            ["tipper real", "tipper imag"],
            loc="lower center",
            ncol=2,
            frameon=False,
            bbox_to_anchor=(0.68, -0.015),
            fontsize=8,
        )
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
    colors: Optional[Dict[str, str]] = None,  # None → from PYCSAMT_STYLE.mt
    marker         = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.marker
    ms             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.ms
    lw             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.lw
    ls             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.ls
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
    # ── resolve visual style from PYCSAMT_STYLE.mt ───────────────────────
    _mt = PYCSAMT_STYLE.mt
    if marker is _UNSET: marker = _mt.xy.marker
    if ms     is _UNSET: ms     = _mt.xy.ms
    if lw     is _UNSET: lw     = _mt.xy.lw
    if ls     is _UNSET: ls     = _mt.xy.ls

    comps = tuple(c.lower() for c in components)
    cmap = {**_style_col(), **(colors or {})}
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
    colors_meas: Optional[Dict[str, str]] = None,  # None → PYCSAMT_STYLE.mt
    color_fit_te   = _UNSET,   # default: PYCSAMT_STYLE.mt.te.color
    color_fit_tm   = _UNSET,   # default: PYCSAMT_STYLE.mt.tm.color
    marker         = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.marker
    ms             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.ms
    lw             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.lw
    ls_meas        = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.ls
    lw_fit: float  = 2.0,
    ls_fit: str    = "-",
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
    # ── resolve visual style from PYCSAMT_STYLE.mt ───────────────────────
    _mt = PYCSAMT_STYLE.mt
    if color_fit_te is _UNSET: color_fit_te = _mt.te.color
    if color_fit_tm is _UNSET: color_fit_tm = _mt.tm.color
    if marker       is _UNSET: marker       = _mt.xy.marker
    if ms           is _UNSET: ms           = _mt.xy.ms
    if lw           is _UNSET: lw           = _mt.xy.lw
    if ls_meas      is _UNSET: ls_meas      = _mt.xy.ls

    comps = tuple(c.lower() for c in components)
    cmap = {**_style_col(), **(colors_meas or {})}
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
