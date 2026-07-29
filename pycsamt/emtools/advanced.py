# pycsamt/emtools/advanced.py
"""
pycsamt.emtools.advanced
========================

Novel visualisations unique to pycsamt v2.  None of these plots exist in
MTPy, MARE2DEM, ModEM or other standard MT packages.

Functions
---------
plot_impedance_mohr_circles
    Double-Mohr-circle representation of the full impedance tensor: rotational
    trajectories in the complex (Re / Im) plane reveal 1-D / 2-D / 3-D character
    in a single glance.

plot_zt_argand
    Argand-space (complex-plane) trajectory of impedance components.  Period is
    the trajectory parameter; colour indicates depth penetration.

plot_survey_fingerprint
    Compact multi-metric "fingerprint" grid — station × log-period cells coloured
    by six simultaneous quantities.  The entire survey quality pattern fits on one
    page.

plot_sensitivity_depth_section
    Bostick-sensitivity-kernel pseudosection: each (station, period) cell is
    drawn as a coloured bar centred on its Bostick depth, width ∝ sensitivity
    window, colour ∝ ρa.  Shows *where* in depth the data are sensitive.
"""

from __future__ import annotations

from typing import Any

import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize

from ..api.plot import add_colorbar
from ..api.station import PYCSAMT_STATION_RENDERING
from ..api.style import PYCSAMT_STYLE
from ._core import (
    _get_t_block,
    _get_z_block,
    _iter_items,
    _name,
    ensure_sites,
    hide_polar_radius_labels,
)
from .tensor import build_phase_tensor_table

__all__ = [
    "plot_dimensionality_ternary",
    "plot_distortion_radar",
    "plot_tf_coherence_network",
    "plot_strike_stability_bands",
    "plot_impedance_mohr_circles",
    "plot_zt_argand",
    "plot_survey_fingerprint",
    "plot_sensitivity_depth_section",
    "plot_rho_phase_bode",
    "plot_pt_period_clock",
    "plot_apparent_resistivity_polar",
    "plot_apparent_anisotropy_section",
    "plot_dimensionality_depth_profile",
    "plot_mt_composite_section",
    "plot_snr_section",
    "plot_z_invariants_section",
]

_MU0 = 4.0 * np.pi * 1e-7  # magnetic permeability [H/m]
_COMP_IDX = {"xx": (0, 0), "xy": (0, 1), "yx": (1, 0), "yy": (1, 1)}


def _rot2(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s], [-s, c]])


def _z_at_theta(z0: np.ndarray, theta: float) -> np.ndarray:
    """Rotate 2×2 impedance tensor by *theta* radians."""
    R = _rot2(theta)
    return R @ z0 @ R.T


def _axes_list(axes: Any, n: int, *, label: str = "axes") -> list[Any] | None:
    """Return *n* flattened axes, or ``None`` when axes were not supplied."""
    if axes is None:
        return None
    if isinstance(axes, np.ndarray):
        out = list(axes.ravel())
    elif isinstance(axes, (list, tuple)):
        out = list(np.asarray(axes, dtype=object).ravel())
    else:
        out = [axes]
    if len(out) < n:
        raise ValueError(f"{label} must provide at least {n} axes; got {len(out)}.")
    return out[:n]


# ─────────────────────────────────────────────────────────────────────────────
# 1. Impedance Mohr Circles
# ─────────────────────────────────────────────────────────────────────────────


def plot_impedance_mohr_circles(
    sites: Any,
    *,
    station: str | None = None,
    periods: list[float] | None = None,
    n_periods: int = 8,
    n_theta: int = 360,
    components: tuple[str, str] = ("xx", "xy"),
    cmap: str = "plasma",
    alpha: float = 0.75,
    mark_zero: bool = True,
    axes=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Double-Mohr-circle diagram of the MT impedance tensor.

    For each period, the impedance tensor **Z** is rotated through all
    angles θ ∈ [0°, 360°) and the trajectory of a chosen pair of
    components is drawn in the complex plane — one circle per period.
    The two panels show the real and imaginary trajectories separately.

    **Physical interpretation** (Lilley 1998; Weaver et al. 2000):

    * **1-D** — every circle degenerates to a single point; all circles
      are centred at the same location (no off-diagonal trace).
    * **2-D** — circles are distinct but all pass through the origin.
    * **3-D** — circles do **not** pass through the origin; the distance
      of the centre from the origin indicates the degree of 3-D character.

    Parameters
    ----------
    sites : any
        EDI path(s) or Sites collection.
    station : str or None
        Station name.  ``None`` picks the first available.
    periods : list of float or None
        Target periods (s) to draw.  When ``None``, *n_periods* logarithmically
        spaced periods spanning the data range are chosen.
    n_periods : int, default 8
        Number of auto-selected periods.
    n_theta : int, default 360
        Angular resolution of each circle.
    components : (str, str), default ("xx", "xy")
        The two Z components traced on (x, y) axes of each panel.
    cmap : str, default "plasma"
        Colormap for period colour-coding.
    alpha : float, default 0.75
        Circle opacity.
    mark_zero : bool
        Mark the origin and θ=0° starting point on each circle.
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    References
    ----------
    Lilley F.E.M., 1998.  *Magnetotelluric tensor decomposition: Part I,
    Theory for a basic procedure.*  Geophysics, 63, 1885–1897.

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_impedance_mohr_circles
    >>> fig = plot_impedance_mohr_circles(sites, station="S12")
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    st_name = None
    ed_pick = None
    for ii, ed in enumerate(_iter_items(S)):
        nm = _name(ed, ii)
        if station is None or nm == station:
            st_name, ed_pick = nm, ed
            break
    if ed_pick is None:
        raise RuntimeError(f"Station {station!r} not found.")

    Z_obj, z, fr = _get_z_block(ed_pick)
    if z is None:
        raise RuntimeError("No impedance data for selected station.")
    per = 1.0 / np.where(fr == 0, np.nan, fr)

    # pick target period indices
    if periods is not None:
        pidx = [int(np.nanargmin(np.abs(per - p))) for p in periods]
    else:
        lo, hi = np.nanmin(per), np.nanmax(per)
        tgt = np.logspace(np.log10(lo), np.log10(hi), n_periods)
        pidx = [int(np.nanargmin(np.abs(per - p))) for p in tgt]
    pidx = sorted(set(pidx))

    ri0, ci0 = _COMP_IDX[components[0]]
    ri1, ci1 = _COMP_IDX[components[1]]
    theta_arr = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)

    # colour scaling by log10(period)
    per_vals = per[pidx]
    cmap_obj = plt.get_cmap(cmap)
    norm_p = mcolors.LogNorm(vmin=per_vals.min(), vmax=per_vals.max())

    axes_given = _axes_list(axes, 2)
    if axes_given is None:
        if figsize is None:
            figsize = (10, 5)
        fig, (ax_re, ax_im) = plt.subplots(1, 2, figsize=figsize)
    else:
        ax_re, ax_im = axes_given
        fig = ax_re.figure

    for pi in pidx:
        z0 = z[pi]
        col = cmap_obj(norm_p(per[pi]))
        traj_re_x, traj_re_y = [], []
        traj_im_x, traj_im_y = [], []
        for th in theta_arr:
            Zr = _z_at_theta(z0, th)
            traj_re_x.append(Zr[ri0, ci0].real)
            traj_re_y.append(Zr[ri1, ci1].real)
            traj_im_x.append(Zr[ri0, ci0].imag)
            traj_im_y.append(Zr[ri1, ci1].imag)
        # close the loop
        traj_re_x.append(traj_re_x[0])
        traj_re_y.append(traj_re_y[0])
        traj_im_x.append(traj_im_x[0])
        traj_im_y.append(traj_im_y[0])

        ax_re.plot(traj_re_x, traj_re_y, color=col, lw=1.2, alpha=alpha)
        ax_im.plot(traj_im_x, traj_im_y, color=col, lw=1.2, alpha=alpha)
        if mark_zero:
            ax_re.plot(traj_re_x[0], traj_re_y[0], "o", ms=4, color=col, zorder=5)
            ax_im.plot(traj_im_x[0], traj_im_y[0], "o", ms=4, color=col, zorder=5)

    for ax, half in ((ax_re, "Re"), (ax_im, "Im")):
        ax.axhline(0, color="0.55", lw=0.7, ls=":")
        ax.axvline(0, color="0.55", lw=0.7, ls=":")
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlabel(
            f"$\\{half}\\,(Z_{{\\rm {components[0].upper()}}})$  (Ω)",
            fontsize=9,
        )
        ax.set_ylabel(
            f"$\\{half}\\,(Z_{{\\rm {components[1].upper()}}})$  (Ω)",
            fontsize=9,
        )
        ax.set_title(f"{half}(Z) Mohr circles", fontsize=9, fontweight="bold")
        ax.grid(True, alpha=0.2, lw=0.5)
        ax.tick_params(labelsize=8)

    fig.subplots_adjust(left=0.08, right=0.87, top=0.90, bottom=0.10, wspace=0.35)
    cax = fig.add_axes([0.895, 0.15, 0.018, 0.68])
    sm = ScalarMappable(cmap=cmap_obj, norm=norm_p)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("Period (s)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    fig.suptitle(
        title or f"Impedance Mohr circles — {st_name}",
        fontsize=10,
        fontweight="bold",
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 2. Z-T Argand Trajectory
# ─────────────────────────────────────────────────────────────────────────────


def plot_zt_argand(
    sites: Any,
    *,
    station: str | None = None,
    components: tuple[str, ...] = ("xy", "yx"),
    period_range: tuple[float, float] | None = None,
    cmap: str = "viridis",
    lw: float = 1.6,
    ms: float = 4.5,
    arrow_every: int = 4,
    normalize: bool = False,
    axes=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Argand-space trajectory of MT impedance components.

    Each component **Z_ij** is plotted as a curve in the complex plane
    (Re vs Im), parametrised by period.  The curve is colour-coded from
    short periods (shallow) to long periods (deep) using *cmap*.  Arrows
    along the curve show the direction of increasing period.

    **Physical interpretation**:

    * A **1-D** subsurface produces a straight line at 45° through the
      origin for Z_xy (pure resistive + inductive).
    * **2-D** structures bend and rotate the trajectory.
    * **3-D** structures produce strongly curved or looping trajectories.
    * The **winding number** of the trajectory around the origin is related
      to the number of layer interfaces penetrated.

    Parameters
    ----------
    sites : any
        EDI path(s) or Sites collection.
    station : str or None
        Station name.  ``None`` → first station.
    components : tuple of {"xx","xy","yx","yy"}, default ("xy","yx")
        Which Z components to trace.  One sub-panel per component.
    period_range : (T_min, T_max) or None
    cmap : str, default "viridis"
        Colormap for period (shallow = dark, deep = bright).
    lw : float
        Trajectory line-width.
    ms : float
        Marker size at each period point.
    arrow_every : int
        Draw a direction arrow every *arrow_every* points.
    normalize : bool
        If ``True``, normalise each trajectory to unit magnitude for
        shape comparison.
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_zt_argand
    >>> fig = plot_zt_argand(
    ...     sites, station="S12", components=("xy", "yx"), normalize=False
    ... )
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    st_name = None
    ed_pick = None
    for ii, ed in enumerate(_iter_items(S)):
        nm = _name(ed, ii)
        if station is None or nm == station:
            st_name, ed_pick = nm, ed
            break
    if ed_pick is None:
        raise RuntimeError(f"Station {station!r} not found.")

    Z_obj, z, fr = _get_z_block(ed_pick)
    if z is None:
        raise RuntimeError("No impedance data.")
    per = 1.0 / np.where(fr == 0, np.nan, fr)

    # period mask
    mask = np.isfinite(per)
    if period_range is not None:
        lo, hi = float(period_range[0]), float(period_range[1])
        mask &= (per >= lo) & (per <= hi)
    per = per[mask]
    z = z[mask]
    order = np.argsort(per)
    per = per[order]
    z = z[order]

    n_comp = len(components)
    axes_given = _axes_list(axes, n_comp)
    if axes_given is None:
        if figsize is None:
            figsize = (4.5 * n_comp, 4.5)
        fig, axes_list = plt.subplots(1, n_comp, figsize=figsize)
        if n_comp == 1:
            axes_list = [axes_list]
        else:
            axes_list = list(axes_list)
    else:
        axes_list = axes_given
        fig = axes_list[0].figure

    cmap_obj = plt.get_cmap(cmap)
    lper = np.log10(per)
    lper_norm = (lper - lper.min()) / (lper.max() - lper.min() + 1e-12)

    for ci, comp in enumerate(components):
        ri, cj = _COMP_IDX[comp.lower()]
        ax = axes_list[ci]
        vals = z[:, ri, cj]
        if normalize:
            scale = np.abs(vals).max() + 1e-30
            vals = vals / scale
        re_v = vals.real
        im_v = vals.imag
        cmap_obj(lper_norm)

        # gradient coloured line segments
        for k in range(len(per) - 1):
            ax.plot(
                re_v[k : k + 2],
                im_v[k : k + 2],
                color=cmap_obj((lper_norm[k] + lper_norm[k + 1]) / 2),
                lw=lw,
                solid_capstyle="round",
            )
        ax.scatter(re_v, im_v, c=lper_norm, cmap=cmap, s=ms**2, zorder=4)

        # direction arrows
        for k in range(0, len(per) - 1, arrow_every):
            dx = re_v[k + 1] - re_v[k]
            dy = im_v[k + 1] - im_v[k]
            ax.annotate(
                "",
                xy=(re_v[k] + 0.6 * dx, im_v[k] + 0.6 * dy),
                xytext=(re_v[k], im_v[k]),
                arrowprops=dict(
                    arrowstyle="-|>", color="0.3", lw=0.6, mutation_scale=8
                ),
            )

        # mark start (short T) and end (long T)
        ax.plot(
            re_v[0],
            im_v[0],
            "*",
            ms=9,
            color="gold",
            zorder=6,
            markeredgecolor="0.3",
            markeredgewidth=0.5,
        )
        ax.plot(
            re_v[-1],
            im_v[-1],
            "^",
            ms=7,
            color="red",
            zorder=6,
            markeredgecolor="0.3",
            markeredgewidth=0.5,
        )

        ax.axhline(0, color="0.55", lw=0.7, ls=":")
        ax.axvline(0, color="0.55", lw=0.7, ls=":")
        ax.set_aspect("equal", adjustable="datalim")
        unit = "" if normalize else "  (Ω)"
        ax.set_xlabel(f"Re($Z_{{\\rm {comp.upper()}}}$){unit}", fontsize=9)
        ax.set_ylabel(f"Im($Z_{{\\rm {comp.upper()}}}$){unit}", fontsize=9)
        ax.set_title(
            f"$Z_{{\\rm {comp.upper()}}}$ trajectory",
            fontsize=9,
            fontweight="bold",
        )
        ax.grid(True, alpha=0.2, lw=0.5)
        ax.tick_params(labelsize=8)

    fig.subplots_adjust(left=0.08, right=0.87, top=0.90, bottom=0.10, wspace=0.35)
    cax = fig.add_axes([0.895, 0.15, 0.018, 0.68])
    sm = ScalarMappable(cmap=cmap_obj, norm=Normalize(vmin=lper.min(), vmax=lper.max()))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("$\\log_{10}$(Period / s)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle(
        title or f"Impedance Argand trajectory — {st_name}",
        fontsize=10,
        fontweight="bold",
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 3. Survey Fingerprint
# ─────────────────────────────────────────────────────────────────────────────

_FINGERPRINT_QUANTITIES = {
    "skew": dict(label="Skew β (°)", cmap="RdBu_r", sym=True, pct=(5, 95)),
    "ellipt": dict(label="Ellipticity λ", cmap="viridis", sym=False, pct=(2, 98)),
    "theta": dict(label="Strike θ (°)", cmap="hsv", sym=False, pct=(2, 98)),
    "s1": dict(label="φ_max", cmap="plasma", sym=False, pct=(5, 95)),
    "s2": dict(label="φ_min", cmap="plasma", sym=False, pct=(5, 95)),
    "beta": dict(label="|β| (°)", cmap="Reds", sym=False, pct=(5, 95)),
}


def plot_survey_fingerprint(
    sites: Any,
    *,
    quantities: list[str] | None = None,
    period_range: tuple[float, float] | None = None,
    station_order: list[str] | None = None,
    cell_aspect: float = 1.0,
    axes=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Compact multi-metric survey fingerprint grid.

    Plots all stations × periods as a colour-coded image for several
    simultaneous physical quantities derived from the phase tensor.  The
    result is a compact "fingerprint" that reveals spatial and frequency
    patterns across the entire survey on a single page.

    Each row of panels corresponds to one quantity; each column of pixels
    within a panel is one station; each row of pixels is one period
    (log-spaced, top = short period = shallow, bottom = long period = deep).

    **Quantities available** (all derived from
    :func:`~pycsamt.emtools.tensor.build_phase_tensor_table`):

    * ``"skew"``   — skewness β: the 3-D indicator
    * ``"ellipt"`` — ellipticity λ = φ_min/φ_max
    * ``"theta"``  — principal-axis strike angle
    * ``"s1"``     — φ_max (maximum phase)
    * ``"s2"``     — φ_min (minimum phase)
    * ``"beta"``   — |β| (absolute skew)

    Parameters
    ----------
    sites : any
    quantities : list of str or None
        Subset of the available quantities.
        Default: ``["skew", "ellipt", "theta", "s1"]``.
    period_range : (T_min, T_max) or None
    station_order : list of str or None
        Explicit station order along the x-axis.  Auto from data when None.
    cell_aspect : float
        Aspect ratio of individual cells (width / height per cell).
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_survey_fingerprint
    >>> fig = plot_survey_fingerprint(all_sites, period_range=(1e-4, 1.0))
    """

    if quantities is None:
        quantities = ["skew", "ellipt", "theta", "s1"]

    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    axes_given = _axes_list(axes, 1) if axes is not None else None
    if df.empty:
        if axes_given is None:
            fig, ax = plt.subplots()
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no phase tensor data", ha="center", va="center")
        return fig

    if period_range is not None:
        lo, hi = float(period_range[0]), float(period_range[1])
        df = df[(df["period"] >= lo) & (df["period"] <= hi)]

    # station order
    all_stations = (
        list(df["station"].unique()) if station_order is None else station_order
    )
    n_sta = len(all_stations)

    # common log-period grid
    all_per = df["period"].to_numpy(float)
    p_lo = float(np.nanmin(all_per))
    p_hi = float(np.nanmax(all_per))
    n_grid = 40
    per_grid = np.logspace(np.log10(max(p_lo, 1e-6)), np.log10(p_hi), n_grid)

    n_q = len(quantities)
    if figsize is None:
        figsize = (max(8, n_sta * 0.3 + 2.5), 2.8 * n_q + 0.8)

    axes_given = _axes_list(axes, n_q) if axes is not None else None
    if axes_given is None:
        fig = plt.figure(figsize=figsize, constrained_layout=False)
        gs = gridspec.GridSpec(
            n_q,
            1,
            figure=fig,
            hspace=0.35,
            left=0.12,
            right=0.88,
            top=0.92,
            bottom=0.06,
        )
    else:
        fig = axes_given[0].figure
        gs = None

    for qi, qty in enumerate(quantities):
        if qty not in _FINGERPRINT_QUANTITIES:
            continue
        qmeta = _FINGERPRINT_QUANTITIES[qty]
        col = qty if qty in df.columns else ("beta" if qty == "|beta|" else None)
        if col is None or col not in df.columns:
            continue

        ax = axes_given[qi] if axes_given is not None else fig.add_subplot(gs[qi])

        # build station × period image
        img = np.full((n_grid, n_sta), np.nan)
        for si, st in enumerate(all_stations):
            sdf = df[df["station"] == st].sort_values("period")
            if sdf.empty:
                continue
            p_s = sdf["period"].to_numpy(float)
            v_s = sdf[col].to_numpy(float)
            for gi, pg in enumerate(per_grid):
                j = int(np.argmin(np.abs(p_s - pg)))
                if np.abs(np.log10(p_s[j] + 1e-30) - np.log10(pg + 1e-30)) < 0.25:
                    img[gi, si] = v_s[j]

        # colour limits
        vv = img[np.isfinite(img)]
        if vv.size == 0:
            continue
        lo_p, hi_p = (
            float(np.percentile(vv, qmeta["pct"][0])),
            float(np.percentile(vv, qmeta["pct"][1])),
        )
        if qmeta["sym"] and lo_p != hi_p:
            vmax = max(abs(lo_p), abs(hi_p))
            vmin, vmax = -vmax, vmax
        else:
            vmin, vmax = lo_p, hi_p
        if vmin == vmax:
            vmax = vmin + 1.0

        cmap_obj = plt.get_cmap(qmeta["cmap"])
        im = ax.imshow(
            img,
            aspect="auto",
            origin="upper",
            extent=[0, n_sta, np.log10(per_grid[-1]), np.log10(per_grid[0])],
            cmap=cmap_obj,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )

        # axes decoration
        ax.set_ylabel(qmeta["label"], fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        if qi == 0:
            PYCSAMT_STATION_RENDERING.apply(
                ax,
                np.arange(n_sta, dtype=float) + 0.5,
                all_stations,
                preset="pseudosection",
                xlim=(0, n_sta),
            )

        # y axis: log period labels
        y_ticks = [
            t
            for t in np.arange(
                int(np.log10(per_grid[0])) - 1,
                int(np.log10(per_grid[-1])) + 2,
            )
            if np.log10(per_grid[0]) <= t <= np.log10(per_grid[-1])
        ]
        ax.set_yticks(y_ticks)
        ax.set_yticklabels([f"$10^{{{int(t)}}}$" for t in y_ticks], fontsize=7)

        # colorbar
        cax = ax.inset_axes([1.01, 0.05, 0.015, 0.9])
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=6)

    if title:
        fig.suptitle(title, fontsize=10, fontweight="bold")
    else:
        fig.suptitle(
            "Survey fingerprint  (phase-tensor metrics)",
            fontsize=10,
            fontweight="bold",
        )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 4. Sensitivity Depth Section
# ─────────────────────────────────────────────────────────────────────────────


def plot_sensitivity_depth_section(
    sites: Any,
    *,
    component: str = "xy",
    period_range: tuple[float, float] | None = None,
    depth_max: float | None = None,
    depth_unit: str = "km",
    cmap: str = "jet_r",
    alpha_bar: float = 0.55,
    bar_width_fraction: float = 0.70,
    rho_lim: tuple[float, float] | None = None,
    station_order: list[str] | None = None,
    show_bostick_depth: bool = True,
    ax=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Bostick-sensitivity-kernel pseudosection.

    For every (station, period) cell the Bostick penetration depth is
    computed and the cell is rendered as a **vertical bar** centred on
    that depth:

    .. math::

        d_B = \\sqrt{\\frac{\\rho_a}{\\mu_0 \\, 2 \\pi f}}

    The bar's colour encodes the apparent resistivity ρa; its vertical
    extent reflects the sensitivity window Δd ≈ d_B × (d log ρa /
    d log T + 1) / 2.  Overlapping bars from multiple periods build a
    natural depth-smoothed image.

    Unlike a simple Bostick pseudosection (which maps period → depth but
    loses the sensitivity context), this plot explicitly shows **where**
    in depth and **how broadly** each datum is sensitive.

    Parameters
    ----------
    sites : any
    component : {"xy", "yx"}, default "xy"
        Impedance component used for ρa and penetration depth.
    period_range : (T_min, T_max) or None
    depth_max : float or None
        Clip the depth axis at this value (km or m depending on *depth_unit*).
    depth_unit : {"km", "m"}, default "km"
    cmap : str, default "jet_r"
        Colormap for ρa colour-coding.
    alpha_bar : float
        Opacity of each kernel bar.
    bar_width_fraction : float
        Fraction of the inter-station spacing used as bar width.
    rho_lim : (vmin, vmax) or None
        Explicit ρa colour limits.  ``None`` → 5th–95th percentile.
    station_order : list of str or None
    show_bostick_depth : bool
        Overlay a thin line at the Bostick depth (without the bar width).
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_sensitivity_depth_section
    >>> fig = plot_sensitivity_depth_section(
    ...     sites, component="xy", depth_max=5.0
    ... )
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    ri, ci = _COMP_IDX[component.lower()]

    rows = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z_o, z, fr = _get_z_block(ed)
        if z is None:
            continue
        rho = getattr(ed, "rho", None)
        if rho is None:
            continue
        per = 1.0 / np.where(fr == 0, np.nan, fr)
        rho_c = rho[:, ri, ci]
        for k in range(len(fr)):
            if not np.isfinite(per[k]) or not (rho_c[k] > 0):
                continue
            if period_range is not None:
                lo_, hi_ = float(period_range[0]), float(period_range[1])
                if per[k] < lo_ or per[k] > hi_:
                    continue
            d_b = float(np.sqrt(rho_c[k] / (_MU0 * 2 * np.pi * fr[k])))
            # sensitivity half-width via finite-difference d(log rho)/d(log T)
            if k > 0 and k < len(fr) - 1:
                dlr = np.log10(max(rho_c[k + 1], 1e-6)) - np.log10(
                    max(rho_c[k - 1], 1e-6)
                )
                dlt = np.log10(max(per[k + 1], 1e-30)) - np.log10(
                    max(per[k - 1], 1e-30)
                )
                slope = (dlr / dlt) if dlt != 0 else 0.0
            else:
                slope = 0.0
            sens_half = d_b * abs(slope + 1.0) / 2.0 * 0.7
            rows.append(
                dict(
                    station=st,
                    period=per[k],
                    rho=rho_c[k],
                    depth=d_b,
                    sens_half=sens_half,
                )
            )

    if not rows:
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize or (10, 5))
        else:
            fig = ax.figure
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return fig

    import pandas as _pd

    df = _pd.DataFrame(rows)

    # station order
    all_st = (
        station_order if station_order is not None else list(df["station"].unique())
    )
    n_st = len(all_st)
    st_pos = {st: i for i, st in enumerate(all_st)}
    bar_w = bar_width_fraction

    d_scale = 1.0 if depth_unit == "m" else 1e-3
    df["depth_plot"] = df["depth"] * d_scale
    df["sens_half_plot"] = df["sens_half"] * d_scale

    rho_all = df["rho"].to_numpy(float)
    if rho_lim is not None:
        vmin, vmax = float(rho_lim[0]), float(rho_lim[1])
    else:
        vmin = float(np.percentile(rho_all, 5))
        vmax = float(np.percentile(rho_all, 95))
    vmin = max(vmin, 1e-6)

    cmap_obj = plt.get_cmap(cmap)
    log_norm = LogNorm(vmin=vmin, vmax=vmax)

    d_all_plot = df["depth_plot"].to_numpy(float)
    d_max_plot = float(np.nanmax(d_all_plot)) * 1.05
    if depth_max is not None:
        d_max_plot = float(depth_max)

    if ax is None:
        if figsize is None:
            figsize = (max(9, n_st * 0.35 + 2), 5.5)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    for _, row in df.iterrows():
        st = row["station"]
        if st not in st_pos:
            continue
        xi = float(st_pos[st])
        d = float(row["depth_plot"])
        dh = max(float(row["sens_half_plot"]), d * 0.05)
        rho = float(row["rho"])
        col = cmap_obj(log_norm(np.clip(rho, vmin, vmax)))
        d_lo = max(0.0, d - dh)
        d_hi = min(d_max_plot, d + dh)
        ax.add_patch(
            plt.Rectangle(
                (xi - bar_w / 2, d_lo),
                bar_w,
                d_hi - d_lo,
                facecolor=col,
                edgecolor="none",
                alpha=alpha_bar,
            )
        )
        if show_bostick_depth:
            ax.plot(
                [xi - bar_w / 2, xi + bar_w / 2],
                [d, d],
                color="k",
                lw=0.5,
                alpha=0.4,
                zorder=3,
            )

    ax.set_ylim(d_max_plot, 0)
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(n_st, dtype=float),
        all_st,
        preset="pseudosection",
        xlim=(-0.5, n_st - 0.5),
    )
    ax.set_ylabel(f"Depth  ({depth_unit})", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(True, which="both", axis="y", alpha=0.2, lw=0.5)

    sm = ScalarMappable(cmap=cmap_obj, norm=log_norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label(r"$\rho_a$  (Ω·m)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    ax.set_title(
        title
        or f"Bostick sensitivity-kernel section  (Z$_{{\\rm {component.upper()}}}$)",
        fontsize=10,
        fontweight="bold",
        pad=8,
    )
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 5. Dimensionality Ternary Diagram
# ─────────────────────────────────────────────────────────────────────────────


def _ternary_coords(u1, u2, u3):
    """Map (1D, 2D, 3D) soft memberships to 2-D Cartesian ternary coords."""
    # vertices: 1D=(0,0), 2D=(1,0), 3D=(0.5, sqrt(3)/2)
    x = u2 + 0.5 * u3
    y = (np.sqrt(3) / 2.0) * u3
    return x, y


def _draw_ternary_frame(ax, labels=("1-D", "2-D", "3-D"), gridlines=(0.25, 0.50, 0.75)):
    """Draw the ternary triangle, axis labels, and grid lines."""
    tri_x = [0.0, 1.0, 0.5, 0.0]
    tri_y = [0.0, 0.0, np.sqrt(3) / 2, 0.0]
    ax.plot(tri_x, tri_y, "k-", lw=1.5, zorder=5)

    # grid lines parallel to each axis
    for frac in gridlines:
        # parallel to 1D–2D edge (constant 3D fraction)
        x0, y0 = _ternary_coords(frac, 1 - frac, 0)
        x1, y1 = _ternary_coords(0, 1 - frac, frac)
        ax.plot([x0, x1], [y0, y1], ":", color="0.72", lw=0.6, zorder=2)
        # parallel to 2D–3D edge (constant 1D fraction)
        x0, y0 = _ternary_coords(0, frac, 1 - frac)
        x1, y1 = _ternary_coords(1 - frac, 0, frac)
        ax.plot([x0, x1], [y0, y1], ":", color="0.72", lw=0.6, zorder=2)
        # parallel to 1D–3D edge (constant 2D fraction)
        x0, y0 = _ternary_coords(frac, 0, 1 - frac)
        x1, y1 = _ternary_coords(0, frac, 1 - frac)
        ax.plot([x0, x1], [y0, y1], ":", color="0.72", lw=0.6, zorder=2)

    # corner labels
    ax.text(
        -0.06,
        -0.04,
        labels[0],
        ha="center",
        va="top",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        1.06,
        -0.04,
        labels[1],
        ha="center",
        va="top",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        0.50,
        np.sqrt(3) / 2 + 0.04,
        labels[2],
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )

    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.10, np.sqrt(3) / 2 + 0.12)
    ax.set_aspect("equal")
    ax.axis("off")


def plot_dimensionality_ternary(
    sites: Any,
    *,
    beta_thresh: float = 5.0,
    ellipt_thresh: float = 0.1,
    period_range: tuple[float, float] | None = None,
    color_by: str = "period",
    cmap: str = "plasma",
    ms: float = 4.0,
    alpha: float = 0.65,
    add_density: bool = True,
    ax=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """1-D / 2-D / 3-D ternary classification diagram.

    Every (station, period) observation is mapped to a **ternary triangle**
    whose vertices represent pure 1-D, 2-D, and 3-D character.  The position
    within the triangle is determined by two continuous membership functions
    derived from the phase tensor:

    .. math::

        u_{3D} = \\min(1, |\\beta| / \\beta_{thresh})

        u_{1D} = (1 - u_{3D}) \\cdot \\max(0,\\, 1 - \\lambda / \\lambda_{thresh})

        u_{2D} = 1 - u_{1D} - u_{3D}

    where β is the phase-tensor skewness and λ is the ellipticity.

    Unlike the standard traffic-light grid, the ternary diagram reveals the
    **continuous spread** of the dataset's dimensionality: where the cloud
    sits, how tightly clustered it is, and whether it bridges two regimes.

    Parameters
    ----------
    sites : any
    beta_thresh : float, default 5.0
        Skewness |β| (°) above which the 3-D membership saturates to 1.
    ellipt_thresh : float, default 0.1
        Ellipticity λ above which the 2-D membership saturates to 1.
    period_range : (T_min, T_max) or None
    color_by : {"period", "station", "skew", "ellipt"}
        Quantity mapped to point colour.
    cmap : str, default "plasma"
    ms : float
        Scatter marker size.
    alpha : float
        Scatter opacity.
    add_density : bool
        Overlay a hexbin density map behind the scatter.
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_dimensionality_ternary
    >>> fig = plot_dimensionality_ternary(
    ...     sites, beta_thresh=5.0, ellipt_thresh=0.1
    ... )
    """
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize or (6, 6))
        else:
            fig = ax.figure
        ax.text(0.5, 0.5, "no phase tensor data", ha="center", va="center")
        return fig

    if period_range is not None:
        lo, hi = float(period_range[0]), float(period_range[1])
        df = df[(df["period"] >= lo) & (df["period"] <= hi)]

    ellipt = df["ellipt"].to_numpy(float)
    beta = np.abs(df["beta"].to_numpy(float))
    per = df["period"].to_numpy(float)

    # soft memberships (sum = 1)
    u3d = np.clip(beta / float(beta_thresh), 0.0, 1.0)
    u1d = (1 - u3d) * np.clip(1.0 - ellipt / float(ellipt_thresh), 0.0, 1.0)
    u2d = 1.0 - u1d - u3d

    x, y = _ternary_coords(u1d, u2d, u3d)

    # colour values
    if color_by == "period":
        c_vals = np.log10(np.clip(per, 1e-10, None))
        cb_label = "$\\log_{10}$(Period / s)"
        cm = cmap
    elif color_by == "skew":
        c_vals = np.abs(df["skew"].to_numpy(float))
        cb_label = "|skew β| (°)"
        cm = "Reds"
    elif color_by == "ellipt":
        c_vals = ellipt
        cb_label = "Ellipticity λ"
        cm = "viridis"
    else:
        # station → integer index
        sts = list(df["station"].unique())
        c_vals = np.array([sts.index(s) for s in df["station"]], float)
        cb_label = "Station index"
        cm = cmap

    if ax is None:
        if figsize is None:
            figsize = (7.5, 7.0)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    _draw_ternary_frame(ax)

    if add_density:
        ax.hexbin(
            x,
            y,
            gridsize=20,
            cmap="Greys",
            alpha=0.25,
            zorder=1,
            extent=[0, 1, 0, np.sqrt(3) / 2],
        )

    sc = ax.scatter(
        x, y, c=c_vals, cmap=cm, s=ms**2, alpha=alpha, linewidths=0, zorder=4
    )
    cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02, shrink=0.7)
    cb.set_label(cb_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # threshold boundary lines
    # |beta| = beta_thresh ↔ u3d = 1 → y = sqrt(3)/2 line segment
    # line where u3d = 0.5 (|beta| = beta_thresh/2)
    u3_half = 0.5
    x_line = np.linspace(0, 1 - u3_half, 50)
    u1_line = np.clip(1 - u3_half - x_line, 0, 1)
    x_l, y_l = _ternary_coords(u1_line, x_line, np.full_like(x_line, u3_half))
    ax.plot(
        x_l,
        y_l,
        "--",
        color="tomato",
        lw=1.2,
        alpha=0.7,
        label=f"|β| = {beta_thresh / 2:.1f}°",
    )
    # line where ellipt = ellipt_thresh (given u3d = 0)
    u2_line = np.linspace(0, 1, 50)
    u1_line = 1 - u2_line
    x_l, y_l = _ternary_coords(u1_line, u2_line, np.zeros_like(u2_line))
    # mark the point where u2=ellipt_thresh
    xi, yi = _ternary_coords(
        np.array([1 - ellipt_thresh / 1.0]),
        np.array([min(ellipt_thresh / 1.0, 1.0)]),
        np.array([0.0]),
    )
    ax.plot(
        xi,
        yi,
        "^",
        ms=9,
        color="steelblue",
        zorder=6,
        label=f"λ = {ellipt_thresh:.2f}",
    )

    ax.legend(
        fontsize=8,
        loc="lower center",
        framealpha=0.85,
        bbox_to_anchor=(0.5, -0.05),
        ncol=2,
    )
    ax.set_title(
        title
        or f"Dimensionality ternary  (β$_{{thresh}}$={beta_thresh}°, "
        f"λ$_{{thresh}}$={ellipt_thresh})",
        fontsize=10,
        fontweight="bold",
        pad=12,
    )
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 6. Galvanic Distortion Decomposition Radar
# ─────────────────────────────────────────────────────────────────────────────

_RADAR_LABELS = [
    "Swift ν",
    "Bahr η",
    "Phase asym.",
    "|β| (skew)",
    "1 − λ",
    "Strike IQR",
]


def _distortion_params(
    z: np.ndarray, rho: np.ndarray, phase: np.ndarray, fr: np.ndarray
) -> np.ndarray:
    """Return 6 distortion proxies in [0,1] for a single station."""
    a, b = z[:, 0, 0], z[:, 0, 1]
    c, d = z[:, 1, 0], z[:, 1, 1]

    # 1. Swift skew ν = |(Zxx + Zyy)| / |(Zxy − Zyx)|  → normalise with ν/(1+ν)
    swift_num = np.abs(a + d)
    swift_den = np.abs(b - c) + 1e-30
    swift = swift_num / swift_den
    p_swift = float(np.nanmedian(swift / (1 + swift)))

    # 2. Bahr η = |(Zxy + Zyx)| / |(Zxy − Zyx)|  → normalise
    bahr_num = np.abs(b + c)
    bahr_den = np.abs(b - c) + 1e-30
    bahr_eta = bahr_num / bahr_den
    p_bahr = float(np.nanmedian(bahr_eta / (1 + bahr_eta)))

    # 3. Phase asymmetry = |φ_xy + φ_yx − 180°| / 90°
    phi_xy = phase[:, 0, 1]
    phi_yx = phase[:, 1, 0]
    p_asym = float(np.nanmedian(np.abs(phi_xy + phi_yx - 180.0) / 90.0))
    p_asym = np.clip(p_asym, 0.0, 1.0)

    # 4. |β| (skewness) from PT — approximate as arctan(|Zxx+Zyy|/|Zxy-Zyx|)/90
    p_beta = (
        float(np.nanmedian(np.degrees(np.arctan(swift_num / (swift_den + 1e-30)))))
        / 45.0
    )
    p_beta = np.clip(p_beta, 0.0, 1.0)

    # 5. 1 − λ (ellipticity complement) from ρa components
    rho_xy = rho[:, 0, 1]
    rho_yx = rho[:, 1, 0]
    rho_xy + rho_yx + 1e-30
    ellipt = np.minimum(rho_xy, rho_yx) / np.maximum(rho_xy, rho_yx)
    p_ellipt = float(1.0 - np.nanmedian(ellipt))
    p_ellipt = np.clip(p_ellipt, 0.0, 1.0)

    # 6. Strike stability: IQR of optimal rotation angle normalised by 90°
    scores = []
    th_arr = np.linspace(-np.pi / 2, np.pi / 2, 91)
    for zi in z:
        best_val = np.inf
        best_ang = 0.0
        for th in th_arr:
            c_, s_ = np.cos(th), np.sin(th)
            R = np.array([[c_, s_], [-s_, c_]])
            Zr = R @ zi @ R.T
            score = float(np.abs(Zr[0, 0]) + np.abs(Zr[1, 1]))
            if score < best_val:
                best_val = score
                best_ang = float(np.degrees(th))
        scores.append(best_ang)
    iqr = float(np.percentile(scores, 75) - np.percentile(scores, 25))
    p_iqr = np.clip(iqr / 90.0, 0.0, 1.0)

    return np.array([p_swift, p_bahr, p_asym, p_beta, p_ellipt, p_iqr], dtype=float)


def plot_distortion_radar(
    sites: Any,
    *,
    stations: list[str] | None = None,
    max_stations: int = 8,
    period_range: tuple[float, float] | None = None,
    fill_alpha: float = 0.18,
    line_alpha: float = 0.85,
    lw: float = 1.5,
    cmap: str = "tab10",
    ax=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Galvanic-distortion decomposition radar chart.

    Each station is rendered as a filled polygon on a **six-axis radar
    (spider) chart**.  The axes encode six independent distortion proxies
    derived from the impedance tensor and apparent resistivity:

    1. **Swift ν** — |(Zxx + Zyy)| / |(Zxy − Zyx)|, normalised.
       Measures departure from the 2-D "anti-diagonal" structure.
    2. **Bahr η** — |(Zxy + Zyx)| / |(Zxy − Zyx)|, normalised.
       Quantifies the symmetric off-diagonal contamination.
    3. **Phase asymmetry** — |φ_xy + φ_yx − 180°| / 90°.
       Zero for 1-D/2-D; increases with galvanic mixing.
    4. **|β| skewness** — scaled phase-tensor skewness.
    5. **1 − λ** — ellipticity complement; high = near-isotropic (1-D).
    6. **Strike IQR** — interquartile range of the sweep-optimal strike
       over all frequencies; large = unstable strike = 3-D/distorted.

    Pure 2-D galvanic distortion (twist + shear only) produces a
    characteristic **narrow polygon aligned with axes 1–3**; true 3-D
    structures also inflate axes 4–5.

    Parameters
    ----------
    sites : any
    stations : list of str or None
        Station names to display.  Auto-selects *max_stations* when None.
    max_stations : int, default 8
    period_range : (T_min, T_max) or None
    fill_alpha, line_alpha, lw : visual parameters.
    cmap : str, default "tab10"
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_distortion_radar
    >>> fig = plot_distortion_radar(sites, stations=["S05", "S12", "S20"])
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    eds_all = [(i, _name(ed, i), ed) for i, ed in enumerate(_iter_items(S))]
    if stations is not None:
        eds_all = [(i, nm, ed) for i, nm, ed in eds_all if nm in stations]
    if len(eds_all) > max_stations:
        step = max(1, len(eds_all) // max_stations)
        eds_all = eds_all[::step][:max_stations]

    if not eds_all:
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize or (6, 6))
        else:
            fig = ax.figure
        ax.text(0.5, 0.5, "no stations", ha="center", va="center")
        return fig

    n_axes = len(_RADAR_LABELS)
    angles = np.linspace(0, 2 * np.pi, n_axes, endpoint=False).tolist()
    angles_plot = angles + [angles[0]]

    cmap_obj = plt.get_cmap(cmap)
    if ax is None:
        if figsize is None:
            figsize = (7.5, 7.0)
        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))
    else:
        fig = ax.figure

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.degrees(angles), labels=_RADAR_LABELS, fontsize=9)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    hide_polar_radius_labels(ax)
    ax.grid(True, alpha=0.3, lw=0.7)

    for k, (_, nm, ed) in enumerate(eds_all):
        Z_o, z, fr = _get_z_block(ed)
        if z is None:
            continue
        rho = getattr(ed, "rho", None)
        phase = getattr(ed, "phase", None)
        if rho is None or phase is None:
            continue
        per = 1.0 / np.where(fr == 0, np.nan, fr)
        if period_range is not None:
            lo_, hi_ = float(period_range[0]), float(period_range[1])
            mask = np.isfinite(per) & (per >= lo_) & (per <= hi_)
        else:
            mask = np.isfinite(per)
        if not mask.any():
            continue

        params = _distortion_params(z[mask], rho[mask], phase[mask], fr[mask])
        vals = params.tolist() + [params[0]]
        col = cmap_obj(k / max(len(eds_all) - 1, 1))

        ax.plot(angles_plot, vals, lw=lw, color=col, alpha=line_alpha, label=nm)
        ax.fill(angles_plot, vals, color=col, alpha=fill_alpha)
        # mark each vertex
        ax.scatter(angles, params, s=25, color=col, zorder=5, alpha=0.9)

    ax.legend(
        loc="upper right",
        bbox_to_anchor=(1.35, 1.15),
        fontsize=8,
        framealpha=0.85,
    )
    ax.set_title(
        title or "Galvanic distortion radar",
        fontsize=10,
        fontweight="bold",
        pad=20,
    )
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 7. Transfer-Function Coherence Network
# ─────────────────────────────────────────────────────────────────────────────


def plot_tf_coherence_network(
    sites: Any,
    *,
    component: str = "xy",
    period_range: tuple[float, float] | None = None,
    threshold: float = 0.85,
    max_edges: int = 120,
    node_c_by: str = "skew",
    node_cmap: str = "RdBu_r",
    edge_cmap: str = "YlOrRd",
    node_ms: float = 8.0,
    lw_max: float = 2.5,
    alpha_edge: float = 0.65,
    ax=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Inter-station transfer-function coherence network map.

    Stations are placed at their geographic positions.  Pairs of stations
    whose log₁₀(ρa) curves are **correlated above** *threshold* (Pearson r)
    are connected by an edge: edge width ∝ r, edge colour ∝ r.  Isolated
    nodes (no connections above *threshold*) represent data outliers or
    strongly localised 3-D anomalies.

    The node colour encodes a per-station summary quantity (*node_c_by*):

    * ``"skew"``   — median |β| (high = 3-D)
    * ``"ellipt"`` — median ellipticity
    * ``"rho"``    — median log₁₀(ρa)

    Parameters
    ----------
    sites : any
    component : {"xy", "yx"}, default "xy"
        ρa component used to compute pairwise correlation.
    period_range : (T_min, T_max) or None
    threshold : float, default 0.85
        Minimum Pearson r to draw an edge.
    max_edges : int, default 120
        Maximum edges drawn (highest-r first).
    node_c_by : {"skew", "ellipt", "rho"}
    node_cmap : str, default "RdBu_r"
    edge_cmap : str, default "YlOrRd"
    node_ms : float
    lw_max, alpha_edge : edge visual parameters.
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_tf_coherence_network
    >>> fig = plot_tf_coherence_network(all_sites, threshold=0.90)
    """
    from scipy.stats import pearsonr as _pearsonr

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    ri, ci = _COMP_IDX[component.lower()]

    # 1. collect per-station rho curves + coords
    sta_data = {}
    for i, ed in enumerate(_iter_items(S)):
        nm = _name(ed, i)
        c = getattr(ed, "coords", None)
        # A missing EDI >HEAD LAT/LONG (e.g. KAP03) surfaces as NaN, not
        # None; "is None" alone lets NaN through and eventually crashes
        # ax.set_aspect() below.
        if (
            c is None
            or c[0] is None
            or c[1] is None
            or not np.isfinite(c[0])
            or not np.isfinite(c[1])
        ):
            continue
        lat, lon = float(c[0]), float(c[1])
        Z_o, z, fr = _get_z_block(ed)
        rho = getattr(ed, "rho", None)
        if rho is None or z is None:
            continue
        per = 1.0 / np.where(fr == 0, np.nan, fr)
        mask = np.isfinite(per)
        if period_range is not None:
            lo_, hi_ = float(period_range[0]), float(period_range[1])
            mask &= (per >= lo_) & (per <= hi_)
        if not mask.any():
            continue
        rho_c = np.log10(np.clip(rho[mask, ri, ci], 1e-6, None))
        per_c = per[mask]
        sta_data[nm] = dict(lat=lat, lon=lon, rho=rho_c, per=per_c)

    if len(sta_data) < 2:
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize or (8, 6))
        else:
            fig = ax.figure
        ax.text(0.5, 0.5, "insufficient coord data", ha="center", va="center")
        return fig

    # 2. common period grid for interpolation
    all_per = np.concatenate([v["per"] for v in sta_data.values()])
    lo_g = float(np.nanmin(all_per))
    hi_g = float(np.nanmax(all_per))
    per_grid = np.logspace(np.log10(max(lo_g, 1e-6)), np.log10(hi_g), 40)

    names = list(sta_data.keys())
    n_st = len(names)
    rho_mat = np.full((n_st, len(per_grid)), np.nan)
    for k, nm in enumerate(names):
        v = sta_data[nm]
        idx = np.argsort(v["per"])
        p_s = v["per"][idx]
        r_s = v["rho"][idx]
        finite = np.isfinite(p_s) & np.isfinite(r_s)
        if finite.sum() >= 2:
            rho_mat[k] = np.interp(per_grid, p_s[finite], r_s[finite])

    # 3. pairwise Pearson r
    edges = []
    for i in range(n_st):
        for j in range(i + 1, n_st):
            row_i = rho_mat[i]
            row_j = rho_mat[j]
            ok = np.isfinite(row_i) & np.isfinite(row_j)
            if ok.sum() < 4:
                continue
            try:
                r, _ = _pearsonr(row_i[ok], row_j[ok])
            except Exception:
                continue
            if r >= float(threshold):
                edges.append((i, j, float(r)))

    # keep strongest edges
    edges.sort(key=lambda e: -e[2])
    edges = edges[:max_edges]

    # 4. node colour quantity
    df_pt = build_phase_tensor_table(
        S, recursive=False, on_dup=on_dup, strict=False, verbose=0
    )
    node_vals = {}
    for nm in names:
        sdf = df_pt[df_pt["station"] == nm] if not df_pt.empty else None
        if sdf is not None and not sdf.empty:
            if node_c_by == "skew":
                node_vals[nm] = float(np.nanmedian(np.abs(sdf["skew"])))
            elif node_c_by == "rho":
                node_vals[nm] = float(
                    np.nanmedian([sta_data[nm]["rho"]] if nm in sta_data else [np.nan])
                )
            else:
                node_vals[nm] = float(np.nanmedian(sdf["ellipt"]))
        else:
            node_vals[nm] = np.nan

    nv_arr = np.array([node_vals.get(nm, np.nan) for nm in names], float)
    nv_fin = nv_arr[np.isfinite(nv_arr)]
    nv0 = float(np.nanpercentile(nv_fin, 5)) if nv_fin.size else 0.0
    nv1 = float(np.nanpercentile(nv_fin, 95)) if nv_fin.size else 1.0

    lats = np.array([sta_data[nm]["lat"] for nm in names], float)
    lons = np.array([sta_data[nm]["lon"] for nm in names], float)
    lat_r = float(lats.max() - lats.min())
    lon_r = float(lons.max() - lons.min())
    # Ratio of the two spans. A near-linear survey (one span tiny
    # relative to the other -- a north-south AMT line, or a regional
    # profile spanning tens of degrees at depressive angle) is treated
    # as a "profile": true geographic aspect is skipped (it would
    # collapse the plot into an unreadable sliver) and the figsize
    # ratio is capped at 4x so the canvas stays a sane, allocatable
    # size regardless of whether stations are metres or hundreds of
    # kilometres apart.
    raw_ratio = lat_r / lon_r if lon_r > 1e-12 else np.inf
    is_profile = not (0.25 <= raw_ratio <= 4.0)
    geo_ratio = min(max(raw_ratio, 0.25), 4.0) if np.isfinite(raw_ratio) else 4.0

    if figsize is None:
        base = 8.0
        if geo_ratio >= 1.0:
            figsize = (base, base * geo_ratio)
        else:
            figsize = (base / geo_ratio, base)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # draw edges
    r_vals = [e[2] for e in edges]
    r_min = min(r_vals) if r_vals else float(threshold)
    r_max = max(r_vals) if r_vals else 1.0
    ec_obj = plt.get_cmap(edge_cmap)
    er_norm = Normalize(vmin=r_min, vmax=r_max)

    for i_idx, j_idx, r in edges:
        xi, yi = lons[i_idx], lats[i_idx]
        xj, yj = lons[j_idx], lats[j_idx]
        lw = lw_max * (r - r_min) / (r_max - r_min + 1e-12)
        col = ec_obj(er_norm(r))
        ax.plot(
            [xi, xj],
            [yi, yj],
            color=col,
            lw=lw,
            alpha=alpha_edge,
            zorder=2,
            solid_capstyle="round",
        )

    # draw nodes
    plt.get_cmap(node_cmap)
    Normalize(vmin=nv0, vmax=nv1)
    sc = ax.scatter(
        lons,
        lats,
        c=nv_arr,
        cmap=node_cmap,
        vmin=nv0,
        vmax=nv1,
        s=node_ms**2,
        zorder=4,
        linewidths=0.5,
        edgecolors="0.3",
    )

    add_colorbar(sc, ax, label=node_c_by, side="right", size="3%", pad=0.06)

    # edge colorbar
    sm_e = ScalarMappable(cmap=ec_obj, norm=er_norm)
    sm_e.set_array([])
    cb_e = fig.colorbar(sm_e, ax=ax, fraction=0.02, pad=0.12, orientation="vertical")
    cb_e.set_label("Pearson r", fontsize=8)
    cb_e.ax.tick_params(labelsize=7)

    ax.set_xlabel("Longitude", fontsize=9)
    ax.set_ylabel("Latitude", fontsize=9)
    lat_mid = float(lats.mean())
    # A true geographic aspect ratio is right for a 2-D-spread survey,
    # but forcing it on a near-linear profile (one of lat_r/lon_r tiny
    # relative to the other, e.g. a north-south line where longitude
    # barely varies) shrinks the axes box down to an unreadable sliver
    # regardless of the requested figsize -- fall back to "auto" there
    # so the plot actually fills the figure.
    if is_profile:
        ax.set_aspect("auto")
    else:
        ax.set_aspect(1.0 / max(np.cos(np.radians(lat_mid)), 1e-6))
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.15, lw=0.5)
    ax.set_title(
        title
        or f"ρa coherence network  (Z$_{{\\rm {component.upper()}}}$, "
        f"r ≥ {threshold:.2f})",
        fontsize=10,
        fontweight="bold",
    )
    from matplotlib.ticker import ScalarFormatter

    for _axis in (ax.xaxis, ax.yaxis):
        _fmt = ScalarFormatter(useOffset=False)
        _fmt.set_scientific(False)
        _axis.set_major_formatter(_fmt)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 8. Strike Stability Bands
# ─────────────────────────────────────────────────────────────────────────────


def plot_strike_stability_bands(
    sites: Any,
    *,
    methods: tuple[str, ...] = ("sweep", "pt", "tipper"),
    period_range: tuple[float, float] | None = None,
    n_period_bins: int = 30,
    agreement_tol: float = 10.0,
    smooth_window: int = 3,
    method_colors: dict | None = None,
    fill_alpha: float = 0.25,
    line_alpha: float = 0.90,
    consensus_alpha: float = 0.18,
    ax=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Multi-method strike stability band diagram.

    Three strike-estimation methods are evaluated at every period and
    displayed as **coloured ribbons** (median ± 0.5 × IQR across all
    stations).  Where all active methods agree within *agreement_tol*
    degrees, a grey **consensus zone** is shaded.

    The plot answers: *"At which periods is the estimated strike reliable,
    and do different methods agree?"*

    Methods:

    * ``"sweep"``  — impedance-rotation sweep (Z-based), per frequency.
    * ``"pt"``     — phase-tensor principal axis θ.
    * ``"tipper"`` — induction-arrow azimuth, real component
      (skipped silently when no tipper data are found).

    Parameters
    ----------
    sites : any
    methods : tuple, default ("sweep", "pt", "tipper")
    period_range : (T_min, T_max) or None
    n_period_bins : int, default 30
        Number of log-spaced period bins used for ribbon statistics.
    agreement_tol : float, default 10.0
        Strike agreement window in degrees (consensus shading threshold).
    smooth_window : int, default 3
        Moving-average smoothing kernel applied to the median ribbon.
    method_colors : dict or None
        ``{"sweep": color, "pt": color, "tipper": color}``.
    fill_alpha, line_alpha, consensus_alpha : float
        Opacity parameters.
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_strike_stability_bands
    >>> fig = plot_strike_stability_bands(all_sites, agreement_tol=10.0)
    """
    from .strike import strike_curve_sweep as _sweep

    _colors = {"sweep": "#1f77b4", "pt": "#d62728", "tipper": "#2ca02c"}
    if method_colors:
        _colors.update(method_colors)

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    # ── collect raw per-(station, period) strike angles ───────────────────
    data = {}  # method → list of (period, angle_deg) pairs

    if "sweep" in methods:
        try:
            df_sw = _sweep(
                S,
                recursive=False,
                on_dup=on_dup,
                strict=False,
                verbose=verbose,
            )
            if not df_sw.empty:
                per_sw = 1.0 / np.where(
                    df_sw["freq"] == 0, np.nan, df_sw["freq"].to_numpy(float)
                )
                ang_sw = df_sw["ang"].to_numpy(float) % 180.0
                data["sweep"] = np.column_stack([per_sw, ang_sw])
        except Exception:
            pass

    if "pt" in methods:
        df_pt = build_phase_tensor_table(
            S,
            recursive=False,
            on_dup=on_dup,
            strict=False,
            verbose=0,
        )
        if not df_pt.empty:
            per_pt = df_pt["period"].to_numpy(float)
            ang_pt = df_pt["theta"].to_numpy(float) % 180.0
            data["pt"] = np.column_stack([per_pt, ang_pt])

    if "tipper" in methods:
        tip_rows = []
        for _i, ed in enumerate(_iter_items(S)):
            T_o, t, fr = _get_t_block(ed)
            if t is None or fr is None:
                continue
            per_t = 1.0 / np.where(fr == 0, np.nan, fr)
            tx, ty = np.real(t[:, 0]), np.real(t[:, 1])
            az = np.degrees(np.arctan2(ty, tx)) % 180.0
            for p, a in zip(per_t, az):
                if np.isfinite(p) and np.isfinite(a):
                    tip_rows.append([p, a])
        if tip_rows:
            data["tipper"] = np.array(tip_rows, float)

    if not data:
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize or (10, 4))
        else:
            fig = ax.figure
        ax.text(0.5, 0.5, "no strike data", ha="center", va="center")
        return fig

    # ── period grid ───────────────────────────────────────────────────────
    all_per = np.concatenate([v[:, 0] for v in data.values()])
    all_per = all_per[np.isfinite(all_per)]
    p_lo = float(np.nanmin(all_per))
    p_hi = float(np.nanmax(all_per))
    if period_range is not None:
        p_lo = max(p_lo, float(period_range[0]))
        p_hi = min(p_hi, float(period_range[1]))
    per_grid = np.logspace(np.log10(max(p_lo, 1e-8)), np.log10(p_hi), n_period_bins)

    def _smooth(arr, k):
        if k <= 1 or arr.size < k:
            return arr
        w = np.ones(k) / k
        return np.convolve(arr, w, mode="same")

    # ── compute median + IQR ribbons ──────────────────────────────────────
    stats = {}
    for meth, arr in data.items():
        per_m = arr[:, 0]
        ang_m = arr[:, 1]
        med_arr = np.full(len(per_grid), np.nan)
        lo_arr = np.full(len(per_grid), np.nan)
        hi_arr = np.full(len(per_grid), np.nan)
        for gi, pg in enumerate(per_grid):
            # pick points within 0.3 decade of grid node
            dist = np.abs(np.log10(per_m + 1e-30) - np.log10(pg + 1e-30))
            sel = dist < 0.3
            a = ang_m[sel]
            if a.size < 2:
                continue
            # axial circular median (double & halve)
            rad2 = np.radians(2.0 * a)
            med2 = 0.5 * np.degrees(
                np.arctan2(np.sin(rad2).mean(), np.cos(rad2).mean())
            )
            med_arr[gi] = med2 % 180.0
            lo_arr[gi] = float(np.percentile(a, 25))
            hi_arr[gi] = float(np.percentile(a, 75))
        stats[meth] = dict(
            med=_smooth(med_arr, smooth_window),
            lo=_smooth(lo_arr, smooth_window),
            hi=_smooth(hi_arr, smooth_window),
        )

    # ── figure ────────────────────────────────────────────────────────────
    if ax is None:
        if figsize is None:
            figsize = (11, 4.5)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    for meth, st in stats.items():
        col = _colors.get(meth, "C0")
        valid = np.isfinite(st["med"])
        ax.fill_between(
            per_grid[valid],
            st["lo"][valid],
            st["hi"][valid],
            color=col,
            alpha=fill_alpha,
            zorder=2,
        )
        ax.plot(
            per_grid[valid],
            st["med"][valid],
            color=col,
            lw=2.0,
            alpha=line_alpha,
            label=meth,
            zorder=3,
        )

    # ── consensus zone ────────────────────────────────────────────────────
    if len(stats) >= 2:
        meds = np.vstack([st["med"] for st in stats.values()])
        max_spread = np.nanmax(meds, axis=0) - np.nanmin(meds, axis=0)
        agree = max_spread <= float(agreement_tol)
        # shade periods of agreement
        for gi in range(len(per_grid) - 1):
            if agree[gi]:
                ax.axvspan(
                    per_grid[gi],
                    per_grid[gi + 1],
                    color="0.55",
                    alpha=consensus_alpha,
                    zorder=1,
                )

    ax.set_xscale("log")
    ax.set_ylim(0, 180)
    ax.set_xlabel("Period (s)", fontsize=9)
    ax.set_ylabel("Strike angle (°)", fontsize=9)
    ax.set_yticks([0, 45, 90, 135, 180])
    ax.grid(True, which="both", alpha=0.2, lw=0.5)
    ax.grid(True, which="major", alpha=0.45, lw=0.7)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=8, framealpha=0.85, loc="best")

    from matplotlib.patches import Patch

    ax.legend(
        handles=(
            [
                plt.Line2D([0], [0], color=_colors.get(m, "C0"), lw=2, label=m)
                for m in stats
            ]
            + [
                Patch(
                    color="0.55",
                    alpha=0.5,
                    label=f"consensus (±{agreement_tol:.0f}°)",
                )
            ]
        ),
        fontsize=8,
        framealpha=0.85,
    )
    ax.set_title(
        title or f"Strike stability bands  ({', '.join(stats.keys())})",
        fontsize=10,
        fontweight="bold",
    )
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 9. Rho-Phase Bode Consistency Check
# ─────────────────────────────────────────────────────────────────────────────


def plot_rho_phase_bode(
    sites: Any,
    *,
    station: str | None = None,
    component: str = "xy",
    period_range: tuple[float, float] | None = None,
    smooth_window: int = 0,
    axes=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Bode consistency diagram: observed ρa/φ vs Bostick-predicted φ.

    For a minimum-phase medium the phase can be predicted from ρa via

    .. math::

        \\phi_{Bode}(T) \\approx \\frac{\\pi}{4}
                        \\left(1 + \\frac{d\\,\\ln\\rho_a}{d\\,\\ln T}\\right)

    Significant departure of observed φ from φ_Bode indicates galvanic
    distortion or near-field source effects.

    Parameters
    ----------
    sites : any
    station : str or None
    component : {"xy", "yx"}, default "xy"
    period_range : (T_min, T_max) or None
    smooth_window : int
        Half-width (points) of a centred moving-average smoothing kernel.
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_rho_phase_bode
    >>> fig = plot_rho_phase_bode(sites, component="xy")
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    st_name = None
    ed_pick = None
    for ii, ed in enumerate(_iter_items(S)):
        nm = _name(ed, ii)
        if station is None or nm == station:
            st_name, ed_pick = nm, ed
            break
    if ed_pick is None:
        raise RuntimeError(f"Station {station!r} not found.")

    Z_obj, z, fr = _get_z_block(ed_pick)
    if z is None:
        raise RuntimeError("No impedance data.")
    ri, ci = _COMP_IDX[component.lower()]
    per = 1.0 / np.where(fr == 0, np.nan, fr)

    rho_raw = getattr(ed_pick, "rho", None)
    pha_raw = getattr(ed_pick, "phase", None)
    rho = (
        rho_raw[:, ri, ci]
        if rho_raw is not None
        else (0.2 / np.where(fr == 0, np.nan, fr)) * np.abs(z[:, ri, ci]) ** 2
    )
    phi = (
        pha_raw[:, ri, ci]
        if pha_raw is not None
        else np.degrees(np.angle(z[:, ri, ci]))
    )

    mask = np.isfinite(per) & (rho > 0) & np.isfinite(phi)
    if period_range is not None:
        lo, hi = float(period_range[0]), float(period_range[1])
        mask &= (per >= lo) & (per <= hi)
    per, rho, phi = per[mask], rho[mask], phi[mask]
    order = np.argsort(per)
    per, rho, phi = per[order], rho[order], phi[order]

    def _smooth(arr, w):
        if w <= 0 or arr.size < 2 * w + 1:
            return arr
        return np.convolve(arr, np.ones(2 * w + 1) / (2 * w + 1), mode="same")

    rho = _smooth(rho, smooth_window)
    phi = _smooth(phi, smooth_window)

    # Bostick predicted phase
    log_rho = np.log10(np.clip(rho, 1e-10, None))
    log_per = np.log10(np.clip(per, 1e-30, None))
    d_log_rho = np.gradient(log_rho, log_per)
    phi_bode = 45.0 * (1.0 + d_log_rho)

    mt = PYCSAMT_STYLE.mt
    comp_upper = component.upper()
    axes_given = _axes_list(axes, 2)
    if axes_given is None:
        if figsize is None:
            figsize = (9, 6)
        fig, (ax_rho, ax_phi) = plt.subplots(2, 1, figsize=figsize, sharex=True)
    else:
        ax_rho, ax_phi = axes_given
        fig = ax_rho.figure

    ax_rho.loglog(
        per,
        rho,
        color=mt.xy.color,
        lw=1.8,
        marker="o",
        ms=3.5,
        label=f"$\\rho_{{a,{comp_upper}}}$",
    )
    ax_rho.set_ylabel(r"$\rho_a$  (Ω·m)", fontsize=9)
    ax_rho.grid(True, which="both", alpha=0.25, lw=0.5)
    ax_rho.legend(fontsize=8, framealpha=0.85, loc="best")
    ax_rho.tick_params(labelsize=8)
    ax_rho.set_title(
        title or f"Bode consistency — {st_name}  (Z$_{{\\rm {comp_upper}}}$)",
        fontsize=10,
        fontweight="bold",
    )

    ax_phi.semilogx(
        per,
        phi,
        color=mt.xy.color,
        lw=1.8,
        marker="o",
        ms=3.5,
        label="$\\phi_{obs}$",
    )
    ax_phi.semilogx(
        per,
        phi_bode,
        color="tomato",
        lw=1.5,
        ls="--",
        label="$\\phi_{Bode}$ (predicted)",
    )
    ax_phi.fill_between(
        per, phi, phi_bode, color="tomato", alpha=0.15, label="discrepancy"
    )
    ax_phi.axhline(45, color="0.55", lw=0.8, ls=":")
    ax_phi.set_ylim(0, 90)
    ax_phi.set_ylabel("Phase  (°)", fontsize=9)
    ax_phi.set_xlabel("Period (s)", fontsize=9)
    ax_phi.grid(True, which="both", alpha=0.25, lw=0.5)
    ax_phi.legend(fontsize=8, framealpha=0.85, loc="best")
    ax_phi.tick_params(labelsize=8)

    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 10. Phase-Tensor Period Clock
# ─────────────────────────────────────────────────────────────────────────────


def plot_pt_period_clock(
    sites: Any,
    *,
    station: str | None = None,
    n_rings: int = 6,
    period_range: tuple[float, float] | None = None,
    cmap: str = "plasma",
    ax=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Phase-tensor period clock: PT ellipses on concentric rings.

    Each concentric ring = one period (inner = short = shallow; outer =
    long = deep).  The PT ellipse is placed at the top of each ring,
    oriented by θ (strike) and elongated by λ (ellipticity).  When
    *station* is None, the survey-wide median θ and λ are used.

    Parameters
    ----------
    sites : any
    station : str or None
    n_rings : int, default 6
    period_range : (T_min, T_max) or None
    cmap : str
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_pt_period_clock
    >>> fig = plot_pt_period_clock(all_sites, n_rings=6)
    """
    from matplotlib.patches import Ellipse as _Ell

    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize or (6, 6))
        else:
            fig = ax.figure
        ax.text(0.5, 0.5, "no PT data", ha="center", va="center")
        return fig

    if period_range is not None:
        lo, hi = float(period_range[0]), float(period_range[1])
        df = df[(df["period"] >= lo) & (df["period"] <= hi)]
    if station is not None:
        df = df[df["station"] == station]

    all_per = df["period"].to_numpy(float)
    p_lo = float(np.nanmin(all_per))
    p_hi = float(np.nanmax(all_per))
    per_rings = np.logspace(np.log10(max(p_lo, 1e-10)), np.log10(p_hi), n_rings)

    thetas, lambdas = [], []
    for pg in per_rings:
        dist = np.abs(
            np.log10(df["period"].to_numpy(float) + 1e-30) - np.log10(pg + 1e-30)
        )
        sub = df[dist < 0.3]
        thetas.append(float(np.nanmedian(sub["theta"])) if not sub.empty else np.nan)
        lambdas.append(float(np.nanmedian(sub["ellipt"])) if not sub.empty else np.nan)

    cmap_obj = plt.get_cmap(cmap)
    lp_min = np.log10(per_rings[0] + 1e-30)
    lp_max = np.log10(per_rings[-1] + 1e-30)
    per_norm = Normalize(vmin=lp_min, vmax=lp_max)

    if ax is None:
        if figsize is None:
            figsize = (7.5, 7.5)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_aspect("equal")
    ax.axis("off")

    R_min, R_max = 0.15, 1.0
    r_arr = np.linspace(R_min, R_max, n_rings)
    dr = (R_max - R_min) / max(n_rings - 1, 1)
    a_size = dr * 0.42  # ellipse semi-major half-width in data units

    for r, pg, th, lam in zip(r_arr, per_rings, thetas, lambdas):
        ax.add_patch(plt.Circle((0, 0), r, fill=False, color="0.78", lw=0.7, zorder=1))
        if not (np.isfinite(th) and np.isfinite(lam)):
            continue
        col = cmap_obj(per_norm(np.log10(pg + 1e-30)))
        a = a_size
        b = a * max(1.0 - float(lam), 0.06)
        mpl_ang = 90.0 - float(th)  # geoelectric N-CW → mpl CCW-from-E
        ell = _Ell(
            (0.0, r),
            width=2 * b,
            height=2 * a,
            angle=mpl_ang,
            facecolor=col,
            edgecolor="k",
            linewidth=0.6,
            alpha=0.80,
            zorder=4,
        )
        ax.add_patch(ell)
        ax.text(
            r + 0.04,
            0.0,
            f"$10^{{{np.log10(pg + 1e-30):.1f}}}$ s",
            ha="left",
            va="center",
            fontsize=6.5,
            color=col,
        )

    for ang_d, lbl in ((90, "N"), (0, "E"), (270, "S"), (180, "W")):
        ang_r = np.radians(ang_d)
        rx = (R_max + 0.14) * np.cos(ang_r)
        ry = (R_max + 0.14) * np.sin(ang_r)
        ax.text(
            rx,
            ry,
            lbl,
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
        )

    sm = ScalarMappable(cmap=cmap_obj, norm=per_norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02, shrink=0.65)
    cb.set_label("$\\log_{10}$(Period / s)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    lbl_sta = station or "survey median"
    ax.set_title(
        title or f"Phase-tensor period clock — {lbl_sta}",
        fontsize=10,
        fontweight="bold",
    )
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 11. Apparent Resistivity Polar Diagram
# ─────────────────────────────────────────────────────────────────────────────


def plot_apparent_resistivity_polar(
    sites: Any,
    *,
    station: str | None = None,
    n_periods: int = 8,
    period_range: tuple[float, float] | None = None,
    normalize: bool = True,
    cmap: str = "plasma",
    ax=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Apparent-resistivity polar diagram: ρa(θ) petals per period.

    The impedance tensor is rotated through θ ∈ [0°, 360°) and
    ρa_xy(θ) is computed at each angle.  One petal per period is drawn
    on a polar axes, colour-coded by log(period).

    Parameters
    ----------
    sites : any
    station : str or None
    n_periods : int, default 8
    period_range : (T_min, T_max) or None
    normalize : bool
        Normalise each petal to its maximum so shapes are comparable.
    cmap : str
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_apparent_resistivity_polar
    >>> fig = plot_apparent_resistivity_polar(
    ...     sites, n_periods=8, normalize=True
    ... )
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    st_name = None
    ed_pick = None
    for ii, ed in enumerate(_iter_items(S)):
        nm = _name(ed, ii)
        if station is None or nm == station:
            st_name, ed_pick = nm, ed
            break
    if ed_pick is None:
        raise RuntimeError(f"Station {station!r} not found.")

    Z_obj, z, fr = _get_z_block(ed_pick)
    if z is None:
        raise RuntimeError("No impedance data.")
    per = 1.0 / np.where(fr == 0, np.nan, fr)

    mask = np.isfinite(per)
    if period_range is not None:
        lo, hi = float(period_range[0]), float(period_range[1])
        mask &= (per >= lo) & (per <= hi)
    per, z, fr = per[mask], z[mask], fr[mask]

    # select n_periods log-spaced
    lo_p, hi_p = float(np.nanmin(per)), float(np.nanmax(per))
    tgt = np.logspace(np.log10(lo_p), np.log10(hi_p), n_periods)
    pidx = sorted({int(np.nanargmin(np.abs(per - t))) for t in tgt})

    cmap_obj = plt.get_cmap(cmap)
    per_vals = per[pidx]
    pnorm = mcolors.LogNorm(vmin=per_vals.min(), vmax=per_vals.max())

    theta_arr = np.linspace(0, 2 * np.pi, 361)
    n_th = len(theta_arr)

    if ax is None:
        if figsize is None:
            figsize = (7, 7)
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})
    else:
        fig = ax.figure
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    hide_polar_radius_labels(ax)

    for pi in pidx:
        z0 = z[pi]
        fi = float(fr[pi])
        if fi <= 0:
            continue
        rho_th = np.zeros(n_th)
        for k, th in enumerate(theta_arr):
            Zr = _z_at_theta(z0, th)
            rho_th[k] = (0.2 / fi) * float(np.abs(Zr[0, 1]) ** 2)
        if normalize:
            mx = rho_th.max()
            if mx > 0:
                rho_th /= mx
        col = cmap_obj(pnorm(per[pi]))
        ax.plot(theta_arr, rho_th, color=col, lw=1.5, alpha=0.80)
        ax.fill(theta_arr, rho_th, color=col, alpha=0.12)

    sm = ScalarMappable(cmap=cmap_obj, norm=pnorm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.08, shrink=0.65)
    cb.set_label("Period (s)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    ax.set_title(
        title or f"ρa polar diagram — {st_name}",
        fontsize=10,
        fontweight="bold",
        pad=16,
    )
    ax.tick_params(labelsize=8)
    hide_polar_radius_labels(ax)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 12. Apparent Anisotropy Pseudosection
# ─────────────────────────────────────────────────────────────────────────────


def _build_psection_image(
    S: Any,
    extract_fn,
    *,
    n_grid: int = 50,
    period_range: tuple[float, float] | None = None,
    station_order: list[str] | None = None,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Build a (n_grid × n_sta) pseudosection image using *extract_fn*.

    extract_fn(ed, fr, per) → 1-D value array aligned with fr/per.
    """
    rows: dict = {}
    for i, ed in enumerate(_iter_items(S)):
        nm = _name(ed, i)
        Z_obj, z, fr = _get_z_block(ed)
        if z is None:
            continue
        per = 1.0 / np.where(fr == 0, np.nan, fr)
        try:
            vals = extract_fn(ed, z, fr, per)
        except Exception:
            continue
        if vals is None:
            continue
        mask = np.isfinite(per) & np.isfinite(vals)
        if period_range is not None:
            lo_, hi_ = float(period_range[0]), float(period_range[1])
            mask &= (per >= lo_) & (per <= hi_)
        if not mask.any():
            continue
        rows[nm] = (per[mask], vals[mask])

    all_st = station_order if station_order is not None else list(rows.keys())
    n_st = len(all_st)

    if not rows or n_st == 0:
        return np.full((n_grid, 1), np.nan), np.zeros(n_grid), ["—"]

    all_per = np.concatenate([rows[s][0] for s in all_st if s in rows])
    p_lo = float(np.nanmin(all_per))
    p_hi = float(np.nanmax(all_per))
    per_grid = np.logspace(np.log10(max(p_lo, 1e-10)), np.log10(p_hi), n_grid)

    img = np.full((n_grid, n_st), np.nan)
    for si, st in enumerate(all_st):
        if st not in rows:
            continue
        p_s, v_s = rows[st]
        idx_order = np.argsort(p_s)
        p_s, v_s = p_s[idx_order], v_s[idx_order]
        for gi, pg in enumerate(per_grid):
            j = int(np.argmin(np.abs(p_s - pg)))
            if np.abs(np.log10(p_s[j] + 1e-30) - np.log10(pg + 1e-30)) < 0.3:
                img[gi, si] = v_s[j]
    return img, per_grid, all_st


def plot_apparent_anisotropy_section(
    sites: Any,
    *,
    period_range: tuple[float, float] | None = None,
    show_pt_arrows: bool = False,
    arrow_every: int = 4,
    cmap: str = "RdBu_r",
    vmax: float = 1.0,
    station_order: list[str] | None = None,
    ax=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Apparent-anisotropy pseudosection: log₁₀(ρa_XY / ρa_YX).

    Positive (warm) cells indicate ρa_XY > ρa_YX; negative (cool) cells
    the reverse.  Zero corresponds to isotropic response at that station
    and period.

    Parameters
    ----------
    sites : any
    period_range : (T_min, T_max) or None
    show_pt_arrows : bool
        Overlay phase-tensor principal-axis arrows.
    arrow_every : int
        Draw PT arrows every this many stations.
    cmap : str, default "RdBu_r"
    vmax : float
        Symmetric colour limit in log₁₀ units.
    station_order : list of str or None
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_apparent_anisotropy_section
    >>> fig = plot_apparent_anisotropy_section(sites, period_range=(1e-4, 1.0))
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _extract(ed, z, fr, per):
        rho = getattr(ed, "rho", None)
        if rho is None:
            return None
        r_xy = rho[:, 0, 1]
        r_yx = rho[:, 1, 0]
        return np.log10(np.clip(r_xy, 1e-6, None) / np.clip(r_yx, 1e-6, None))

    img, per_grid, all_st = _build_psection_image(
        S,
        _extract,
        period_range=period_range,
        station_order=station_order,
    )
    n_st = len(all_st)
    n_grid = img.shape[0]

    if ax is None:
        if figsize is None:
            figsize = (max(9, n_st * 0.35 + 2), 5.5)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    lp = np.log10(per_grid + 1e-30)
    ax.imshow(
        img,
        aspect="auto",
        origin="upper",
        extent=(-0.5, n_st - 0.5, lp[-1], lp[0]),
        cmap=cmap,
        vmin=-float(vmax),
        vmax=float(vmax),
        interpolation="nearest",
    )

    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(n_st, dtype=float),
        all_st,
        preset="pseudosection",
        xlim=(-0.5, n_st - 0.5),
    )
    ax.set_ylabel("$\\log_{10}T$ (s)", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)

    sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=-float(vmax), vmax=float(vmax)))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("$\\log_{10}(\\rho_{XY}/\\rho_{YX})$", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    if show_pt_arrows:
        df_pt = build_phase_tensor_table(
            S, recursive=False, on_dup=on_dup, strict=False, verbose=0
        )
        if not df_pt.empty:
            arrow_len = (
                0.55 * (lp[0] - lp[-1]) / max(n_grid - 1, 1) if n_grid > 1 else 0.3
            )
            for si in range(0, n_st, max(1, int(arrow_every))):
                sub = df_pt[df_pt["station"] == all_st[si]]
                if sub.empty:
                    continue
                p_pt = sub["period"].to_numpy(float)
                th_pt = sub["theta"].to_numpy(float)
                for gi, pg in enumerate(per_grid):
                    j = int(
                        np.nanargmin(
                            np.abs(np.log10(p_pt + 1e-30) - np.log10(pg + 1e-30))
                        )
                    )
                    if abs(np.log10(p_pt[j] + 1e-30) - np.log10(pg + 1e-30)) >= 0.3:
                        continue
                    th = th_pt[j]
                    if not np.isfinite(th):
                        continue
                    rad = np.radians(90.0 - th)  # geoelectric N-CW -> plot CCW-from-E
                    dx = arrow_len * np.cos(rad) * 0.6
                    dy = arrow_len * np.sin(rad)
                    ax.plot(
                        [si - dx, si + dx],
                        [lp[gi] - dy, lp[gi] + dy],
                        color="k",
                        lw=0.9,
                        alpha=0.75,
                        solid_capstyle="round",
                        zorder=5,
                    )

    ax.set_title(
        title or "Apparent anisotropy pseudosection",
        fontsize=10,
        fontweight="bold",
        pad=8,
    )
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 13. Dimensionality Depth Profile
# ─────────────────────────────────────────────────────────────────────────────


def plot_dimensionality_depth_profile(
    sites: Any,
    *,
    component: str = "xy",
    beta_thresh: float = 5.0,
    ellipt_thresh: float = 0.1,
    period_range: tuple[float, float] | None = None,
    depth_max: float | None = None,
    depth_unit: str = "km",
    cmap: str = "RdYlGn_r",
    station_order: list[str] | None = None,
    ax=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Dimensionality classification mapped to Bostick depth space.

    Each (station, period) datum is placed at its Bostick penetration
    depth and coloured by the 3-D membership u₃D ∈ [0, 1] (red = 3-D,
    green = 1-D/2-D).

    Parameters
    ----------
    sites : any
    component : {"xy", "yx"}, default "xy"
    beta_thresh : float
    ellipt_thresh : float
    period_range : (T_min, T_max) or None
    depth_max : float or None
    depth_unit : {"km", "m"}
    cmap : str
    station_order : list of str or None
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_dimensionality_depth_profile
    >>> fig = plot_dimensionality_depth_profile(sites, depth_max=5.0)
    """
    df_pt = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    S = ensure_sites(
        sites, recursive=False, on_dup=on_dup, strict=strict, verbose=verbose
    )

    ri, ci = _COMP_IDX[component.lower()]
    d_scale = 1e-3 if depth_unit == "km" else 1.0

    rows = []
    for i, ed in enumerate(_iter_items(S)):
        nm = _name(ed, i)
        Z_obj, z, fr = _get_z_block(ed)
        if z is None:
            continue
        rho_raw = getattr(ed, "rho", None)
        if rho_raw is None:
            continue
        per = 1.0 / np.where(fr == 0, np.nan, fr)
        rho_c = rho_raw[:, ri, ci]
        for k in range(len(fr)):
            if not np.isfinite(per[k]) or not (rho_c[k] > 0):
                continue
            if period_range is not None:
                lo_, hi_ = float(period_range[0]), float(period_range[1])
                if per[k] < lo_ or per[k] > hi_:
                    continue
            d_b = float(np.sqrt(rho_c[k] / (_MU0 * 2 * np.pi * fr[k]))) * d_scale
            # PT membership for this station + period
            sub = df_pt[(df_pt["station"] == nm)] if not df_pt.empty else None
            u3d = np.nan
            if sub is not None and not sub.empty:
                dist = np.abs(
                    np.log10(sub["period"].to_numpy(float) + 1e-30)
                    - np.log10(per[k] + 1e-30)
                )
                j = int(np.argmin(dist))
                if dist[j] < 0.3:
                    beta_v = float(abs(sub["beta"].iloc[j]))
                    float(sub["ellipt"].iloc[j])
                    u3d = min(1.0, beta_v / float(beta_thresh))
            rows.append(dict(station=nm, depth=d_b, u3d=u3d))

    if not rows:
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize or (10, 5))
        else:
            fig = ax.figure
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return fig

    import pandas as _pd

    df = _pd.DataFrame(rows)
    all_st = (
        station_order if station_order is not None else list(df["station"].unique())
    )
    n_st = len(all_st)
    st_pos = {s: i for i, s in enumerate(all_st)}

    d_max_plot = float(df["depth"].quantile(0.97)) * 1.05
    if depth_max is not None:
        d_max_plot = float(depth_max)

    if ax is None:
        if figsize is None:
            figsize = (max(9, n_st * 0.35 + 2), 5.5)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    cmap_obj = plt.get_cmap(cmap)
    u3d_vals = df["u3d"].to_numpy(float)
    u3d_vals[np.isfinite(u3d_vals)]

    for _, row in df.iterrows():
        st = row["station"]
        if st not in st_pos:
            continue
        xi = float(st_pos[st])
        di = float(row["depth"])
        if di > d_max_plot:
            continue
        uv = float(row["u3d"])
        col = cmap_obj(uv) if np.isfinite(uv) else "0.65"
        ax.plot(xi, di, "o", ms=3.5, color=col, alpha=0.75, zorder=3)

    ax.set_ylim(d_max_plot, 0)
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(n_st, dtype=float),
        all_st,
        preset="pseudosection",
        xlim=(-0.5, n_st - 0.5),
    )
    ax.set_ylabel(f"Bostick depth  ({depth_unit})", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(True, axis="y", alpha=0.2, lw=0.5)

    sm = ScalarMappable(cmap=cmap_obj, norm=Normalize(vmin=0.0, vmax=1.0))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("3-D membership $u_{3D}$", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.set_title(
        title or f"Dimensionality depth profile  (β$_{{thresh}}$={beta_thresh}°)",
        fontsize=10,
        fontweight="bold",
        pad=8,
    )
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 14. MT Composite Section
# ─────────────────────────────────────────────────────────────────────────────

_COMPOSITE_META: dict = {
    "rho": dict(
        label=r"$\log_{10}\rho_a$  (Ω·m)",
        cmap="jet_r",
        sym=False,
        log=True,
        pct=(5, 95),
    ),
    "phase": dict(label="Phase  (°)", cmap="plasma", sym=False, log=False, pct=(5, 95)),
    "skew": dict(label="|β| skew  (°)", cmap="Reds", sym=False, log=False, pct=(5, 95)),
    "theta": dict(label="Strike θ  (°)", cmap="hsv", sym=False, log=False, pct=(2, 98)),
    "snr": dict(label="SNR", cmap="RdYlGn", sym=False, log=False, pct=(5, 95)),
}


def plot_mt_composite_section(
    sites: Any,
    *,
    component: str = "xy",
    quantities: list[str] | None = None,
    period_range: tuple[float, float] | None = None,
    station_order: list[str] | None = None,
    axes=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Multi-row aligned pseudosection: ρa, φ, |β|, θ, SNR.

    Up to five aligned rows share the same station axis.  Station labels
    and site-triangle markers appear at the top of the first row via the
    package API.

    Parameters
    ----------
    sites : any
    component : {"xy", "yx"}, default "xy"
    quantities : list of {"rho","phase","skew","theta","snr"} or None
        Default: all five.
    period_range : (T_min, T_max) or None
    station_order : list of str or None
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_mt_composite_section
    >>> fig = plot_mt_composite_section(
    ...     sites, component="xy", quantities=["rho", "phase", "skew"]
    ... )
    """
    if quantities is None:
        quantities = list(_COMPOSITE_META.keys())

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    ri, ci = _COMP_IDX[component.lower()]

    # collect raw data per station
    rho_d: dict = {}
    phi_d: dict = {}
    snr_d: dict = {}
    for i, ed in enumerate(_iter_items(S)):
        nm = _name(ed, i)
        Z_obj, z, fr, ze = _get_z_block(ed, with_errors=True)
        if z is None:
            continue
        per = 1.0 / np.where(fr == 0, np.nan, fr)
        rho_raw = getattr(ed, "rho", None)
        pha_raw = getattr(ed, "phase", None)
        rho_v = (
            rho_raw[:, ri, ci]
            if rho_raw is not None
            else (0.2 / np.where(fr == 0, np.nan, fr)) * np.abs(z[:, ri, ci]) ** 2
        )
        phi_v = (
            pha_raw[:, ri, ci]
            if pha_raw is not None
            else np.degrees(np.angle(z[:, ri, ci]))
        )
        if ze is not None and ze.shape == z.shape:
            snr_v = np.abs(z[:, ri, ci]) / (np.abs(ze[:, ri, ci]) + 1e-30)
        else:
            snr_v = np.full(len(fr), np.nan)
        mask = np.isfinite(per)
        if period_range is not None:
            lo_, hi_ = float(period_range[0]), float(period_range[1])
            mask &= (per >= lo_) & (per <= hi_)
        rho_d[nm] = (per[mask], rho_v[mask])
        phi_d[nm] = (per[mask], phi_v[mask])
        snr_d[nm] = (per[mask], snr_v[mask])

    df_pt = build_phase_tensor_table(
        S, recursive=False, on_dup=on_dup, strict=False, verbose=0
    )

    all_st = station_order if station_order is not None else list(rho_d.keys())
    n_st = len(all_st)
    axes_given = _axes_list(axes, 1) if axes is not None else None
    if n_st == 0:
        if axes_given is None:
            fig, ax = plt.subplots()
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return fig

    all_per = np.concatenate([rho_d[s][0] for s in all_st if s in rho_d])
    p_lo = float(np.nanmin(all_per))
    p_hi = float(np.nanmax(all_per))
    n_grid = 50
    per_grid = np.logspace(np.log10(max(p_lo, 1e-10)), np.log10(p_hi), n_grid)

    def _to_img(data_dict):
        img = np.full((n_grid, n_st), np.nan)
        for si, st in enumerate(all_st):
            if st not in data_dict:
                continue
            p_s, v_s = data_dict[st]
            idx_o = np.argsort(p_s)
            p_s, v_s = p_s[idx_o], v_s[idx_o]
            for gi, pg in enumerate(per_grid):
                j = int(np.argmin(np.abs(p_s - pg)))
                if np.abs(np.log10(p_s[j] + 1e-30) - np.log10(pg + 1e-30)) < 0.3:
                    img[gi, si] = v_s[j]
        return img

    def _pt_img(col):
        if df_pt.empty or col not in df_pt.columns:
            return np.full((n_grid, n_st), np.nan)
        d: dict = {}
        for st in all_st:
            sub = df_pt[df_pt["station"] == st]
            if sub.empty:
                continue
            d[st] = (sub["period"].to_numpy(float), sub[col].to_numpy(float))
        return _to_img(d)

    imgs: dict = {}
    for qty in quantities:
        if qty == "rho":
            img = _to_img(rho_d)
            imgs[qty] = np.log10(np.clip(img, 1e-6, None))
        elif qty == "phase":
            imgs[qty] = _to_img(phi_d)
        elif qty == "skew":
            imgs[qty] = np.abs(_pt_img("beta"))
        elif qty == "theta":
            imgs[qty] = _pt_img("theta")
        elif qty == "snr":
            imgs[qty] = _to_img(snr_d)

    n_q = len([q for q in quantities if q in imgs])
    if n_q == 0:
        if axes_given is None:
            fig, ax = plt.subplots()
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return fig

    axes_given = _axes_list(axes, n_q) if axes is not None else None
    if axes_given is None:
        if figsize is None:
            figsize = (max(9, n_st * 0.35 + 2), 2.5 * n_q + 0.6)
        fig = plt.figure(figsize=figsize, constrained_layout=False)
        gs = gridspec.GridSpec(
            n_q,
            1,
            figure=fig,
            hspace=0.30,
            left=0.10,
            right=0.88,
            top=0.92,
            bottom=0.04,
        )
    else:
        fig = axes_given[0].figure
        gs = None

    lp = np.log10(per_grid + 1e-30)
    valid_qs = [q for q in quantities if q in imgs]
    for qi, qty in enumerate(valid_qs):
        meta = _COMPOSITE_META[qty]
        img = imgs[qty]
        vv = img[np.isfinite(img)]
        if vv.size == 0:
            continue
        vmin = float(np.percentile(vv, meta["pct"][0]))
        vmax_ = float(np.percentile(vv, meta["pct"][1]))
        if vmin == vmax_:
            vmax_ = vmin + 1.0

        ax = axes_given[qi] if axes_given is not None else fig.add_subplot(gs[qi])
        im = ax.imshow(
            img,
            aspect="auto",
            origin="upper",
            extent=(-0.5, n_st - 0.5, lp[-1], lp[0]),
            cmap=meta["cmap"],
            vmin=vmin,
            vmax=vmax_,
            interpolation="nearest",
        )
        ax.set_ylabel(meta["label"], fontsize=7.5)
        ax.tick_params(axis="y", labelsize=7)
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        if qi == 0:
            PYCSAMT_STATION_RENDERING.apply(
                ax,
                np.arange(n_st, dtype=float),
                all_st,
                preset="pseudosection",
                xlim=(-0.5, n_st - 0.5),
            )
        cax = ax.inset_axes([1.01, 0.05, 0.015, 0.9])
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=6)

    comp_upper = component.upper()
    fig.suptitle(
        title or f"MT composite section  (Z$_{{\\rm {comp_upper}}}$)",
        fontsize=10,
        fontweight="bold",
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 15. SNR Pseudosection
# ─────────────────────────────────────────────────────────────────────────────


def plot_snr_section(
    sites: Any,
    *,
    components: tuple[str, ...] = ("xy", "yx"),
    period_range: tuple[float, float] | None = None,
    snr_thresh: float = 3.0,
    cmap: str = "RdYlGn",
    vmax: float = 10.0,
    station_order: list[str] | None = None,
    axes=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Signal-to-noise ratio pseudosection: SNR = |Z| / |Z_err|.

    Each panel shows one impedance component.  A contour at *snr_thresh*
    separates acceptable (green) from poor (red) quality cells.

    Parameters
    ----------
    sites : any
    components : tuple, default ("xy", "yx")
    period_range : (T_min, T_max) or None
    snr_thresh : float, default 3.0
    cmap : str
    vmax : float
        Upper SNR colour limit.
    station_order : list of str or None
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_snr_section
    >>> fig = plot_snr_section(sites, components=("xy", "yx"), snr_thresh=3.0)
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    n_comp = len(components)
    # collect SNR per component
    snr_dicts: list = [{} for _ in range(n_comp)]
    all_st_order: list = []
    for i, ed in enumerate(_iter_items(S)):
        nm = _name(ed, i)
        if nm not in all_st_order:
            all_st_order.append(nm)
        Z_obj, z, fr, ze = _get_z_block(ed, with_errors=True)
        if z is None:
            continue
        per = 1.0 / np.where(fr == 0, np.nan, fr)
        mask = np.isfinite(per)
        if period_range is not None:
            lo_, hi_ = float(period_range[0]), float(period_range[1])
            mask &= (per >= lo_) & (per <= hi_)
        for ci_idx, comp in enumerate(components):
            ri, cj = _COMP_IDX[comp.lower()]
            if ze is not None and ze.shape == z.shape:
                snr_v = np.abs(z[:, ri, cj]) / (np.abs(ze[:, ri, cj]) + 1e-30)
            else:
                snr_v = np.full(len(fr), np.nan)
            snr_dicts[ci_idx][nm] = (per[mask], snr_v[mask])

    all_st = station_order if station_order is not None else all_st_order
    n_st = len(all_st)

    axes_given = _axes_list(axes, 1) if axes is not None else None
    if n_st == 0:
        if axes_given is None:
            fig, ax = plt.subplots()
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return fig

    all_per = np.concatenate([snr_dicts[0][s][0] for s in all_st if s in snr_dicts[0]])
    p_lo = float(np.nanmin(all_per))
    p_hi = float(np.nanmax(all_per))
    n_grid = 50
    per_grid = np.logspace(np.log10(max(p_lo, 1e-10)), np.log10(p_hi), n_grid)

    def _img(dct):
        img = np.full((n_grid, n_st), np.nan)
        for si, st in enumerate(all_st):
            if st not in dct:
                continue
            p_s, v_s = dct[st]
            idx_o = np.argsort(p_s)
            p_s, v_s = p_s[idx_o], v_s[idx_o]
            for gi, pg in enumerate(per_grid):
                j = int(np.argmin(np.abs(p_s - pg)))
                if np.abs(np.log10(p_s[j] + 1e-30) - np.log10(pg + 1e-30)) < 0.3:
                    img[gi, si] = v_s[j]
        return img

    axes_given = _axes_list(axes, n_comp) if axes is not None else None
    if axes_given is None:
        if figsize is None:
            figsize = (max(9, n_st * 0.35 + 2), 4.5 * n_comp)
        fig = plt.figure(figsize=figsize, constrained_layout=False)
        gs = gridspec.GridSpec(
            n_comp,
            1,
            figure=fig,
            hspace=0.30,
            left=0.10,
            right=0.88,
            top=0.92,
            bottom=0.04,
        )
    else:
        fig = axes_given[0].figure
        gs = None

    lp = np.log10(per_grid + 1e-30)
    for pi, comp in enumerate(components):
        img = _img(snr_dicts[pi])
        ax = axes_given[pi] if axes_given is not None else fig.add_subplot(gs[pi])
        im = ax.imshow(
            np.clip(img, 0, float(vmax)),
            aspect="auto",
            origin="upper",
            extent=(-0.5, n_st - 0.5, lp[-1], lp[0]),
            cmap=cmap,
            vmin=0.0,
            vmax=float(vmax),
            interpolation="nearest",
        )
        # threshold contour
        finite_mask = np.isfinite(img)
        if finite_mask.any():
            try:
                ax.contour(
                    np.arange(n_st),
                    lp[::-1],
                    np.where(finite_mask, img, 0.0)[::-1, :],
                    levels=[float(snr_thresh)],
                    colors="k",
                    linewidths=0.8,
                    linestyles="--",
                )
            except Exception:
                pass

        ax.set_ylabel(f"SNR  Z$_{{\\rm {comp.upper()}}}$", fontsize=8.5)
        ax.tick_params(axis="y", labelsize=7)
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        if pi == 0:
            PYCSAMT_STATION_RENDERING.apply(
                ax,
                np.arange(n_st, dtype=float),
                all_st,
                preset="pseudosection",
                xlim=(-0.5, n_st - 0.5),
            )
        cax = ax.inset_axes([1.01, 0.05, 0.015, 0.9])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("SNR", fontsize=7)
        cb.ax.tick_params(labelsize=6)

    fig.suptitle(
        title or "SNR pseudosection  (|Z| / |Z$_{err}$|)",
        fontsize=10,
        fontweight="bold",
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 16. Z Rotation-Invariants Section
# ─────────────────────────────────────────────────────────────────────────────


def plot_z_invariants_section(
    sites: Any,
    *,
    period_range: tuple[float, float] | None = None,
    station_order: list[str] | None = None,
    axes=None,
    figsize: tuple[float, float] | None = None,
    title: str = "",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Four-panel impedance rotation-invariants pseudosection.

    Each panel is a station × log-period image of one invariant:

    1. **Swift ν** = |Zxx + Zyy| / |Zxy − Zyx|  (0 = ideal 2-D)
    2. **Bahr μ** = |Zxy + Zyx| / |Zxy − Zyx|  (0 = no galvanic mixing)
    3. **|det Z|^½** = √|ZxxZyy − ZxyZyx|  (distortion-invariant ρa proxy)
    4. **|tr Z| / ||Zxy| − |Zyx||**  (anisotropy proxy: small when the two
       off-diagonal magnitudes are close, large as they diverge)

    Parameters
    ----------
    sites : any
    period_range : (T_min, T_max) or None
    station_order : list of str or None
    figsize : (float, float) or None
    title : str

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools.advanced import plot_z_invariants_section
    >>> fig = plot_z_invariants_section(sites, period_range=(1e-4, 1.0))
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    _INV_META = [
        dict(label="Swift ν", cmap="Reds", pct=(5, 95), sym=False),
        dict(label="Bahr μ", cmap="Oranges", pct=(5, 95), sym=False),
        dict(label=r"$|\det Z|^{1/2}$", cmap="viridis", pct=(5, 95), sym=False),
        dict(
            label=r"$|\mathrm{tr}\,Z| / |{|Z_{xy}|-|Z_{yx}|}|$",
            cmap="plasma",
            pct=(5, 95),
            sym=False,
        ),
    ]

    def _extract_all(ed, fr, per):
        Z_obj, z, _fr = _get_z_block(ed)
        if z is None:
            return None
        a, b = z[:, 0, 0], z[:, 0, 1]
        c, d = z[:, 1, 0], z[:, 1, 1]
        dz = np.abs(b - c) + 1e-30
        swift = np.abs(a + d) / dz
        bahr = np.abs(b + c) / dz
        det_z = np.sqrt(np.abs(a * d - b * c))
        tr_dz = np.abs(a + d) / dz
        return np.column_stack([swift, bahr, det_z, tr_dz])

    # collect per station
    rows: dict = {}
    all_st_order: list = []
    for i, ed in enumerate(_iter_items(S)):
        nm = _name(ed, i)
        if nm not in all_st_order:
            all_st_order.append(nm)
        Z_obj, z, fr = _get_z_block(ed)
        if z is None:
            continue
        per = 1.0 / np.where(fr == 0, np.nan, fr)
        a, b = z[:, 0, 0], z[:, 0, 1]
        c, d_t = z[:, 1, 0], z[:, 1, 1]
        dz = np.abs(b - c) + 1e-30
        # |tr Z| / dmag: dmag = ||Zxy| - |Zyx|| (difference of off-diagonal
        # *magnitudes*), not |Zxy - Zyx| (dz, the complex difference used
        # by Swift/Bahr above) -- using dz here made this column an exact
        # duplicate of the Swift-nu column. dmag makes this a genuine
        # anisotropy proxy: it is small when |Zxy| and |Zyx| are close
        # (near-isotropic response) and grows as they diverge.
        dmag = np.abs(np.abs(b) - np.abs(c)) + 1e-30
        inv4 = np.column_stack(
            [
                np.abs(a + d_t) / dz,
                np.abs(b + c) / dz,
                np.sqrt(np.abs(a * d_t - b * c)),
                np.abs(a + d_t) / dmag,
            ]
        )
        mask = np.isfinite(per)
        if period_range is not None:
            lo_, hi_ = float(period_range[0]), float(period_range[1])
            mask &= (per >= lo_) & (per <= hi_)
        rows[nm] = (per[mask], inv4[mask])

    all_st = station_order if station_order is not None else all_st_order
    n_st = len(all_st)

    axes_given = _axes_list(axes, 1) if axes is not None else None
    if n_st == 0:
        if axes_given is None:
            fig, ax = plt.subplots()
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return fig

    all_per = np.concatenate([rows[s][0] for s in all_st if s in rows])
    p_lo = float(np.nanmin(all_per))
    p_hi = float(np.nanmax(all_per))
    n_grid = 50
    per_grid = np.logspace(np.log10(max(p_lo, 1e-10)), np.log10(p_hi), n_grid)

    imgs = []
    for inv_idx in range(4):
        img = np.full((n_grid, n_st), np.nan)
        for si, st in enumerate(all_st):
            if st not in rows:
                continue
            p_s, v4 = rows[st]
            v_s = v4[:, inv_idx]
            idx_o = np.argsort(p_s)
            p_s, v_s = p_s[idx_o], v_s[idx_o]
            for gi, pg in enumerate(per_grid):
                j = int(np.argmin(np.abs(p_s - pg)))
                if np.abs(np.log10(p_s[j] + 1e-30) - np.log10(pg + 1e-30)) < 0.3:
                    img[gi, si] = v_s[j]
        imgs.append(img)

    axes_given = _axes_list(axes, 4) if axes is not None else None
    if axes_given is None:
        if figsize is None:
            figsize = (max(9, n_st * 0.35 + 2), 2.8 * 4 + 0.6)
        fig = plt.figure(figsize=figsize, constrained_layout=False)
        gs = gridspec.GridSpec(
            4,
            1,
            figure=fig,
            hspace=0.30,
            left=0.10,
            right=0.88,
            top=0.92,
            bottom=0.04,
        )
    else:
        fig = axes_given[0].figure
        gs = None
    lp = np.log10(per_grid + 1e-30)

    for inv_idx, (img, meta) in enumerate(zip(imgs, _INV_META)):
        vv = img[np.isfinite(img)]
        if vv.size == 0:
            continue
        vmin = float(np.percentile(vv, meta["pct"][0]))
        vmax_ = float(np.percentile(vv, meta["pct"][1]))
        if vmin == vmax_:
            vmax_ = vmin + 1.0

        ax = (
            axes_given[inv_idx]
            if axes_given is not None
            else fig.add_subplot(gs[inv_idx])
        )
        im = ax.imshow(
            img,
            aspect="auto",
            origin="upper",
            extent=(-0.5, n_st - 0.5, lp[-1], lp[0]),
            cmap=meta["cmap"],
            vmin=vmin,
            vmax=vmax_,
            interpolation="nearest",
        )
        ax.set_ylabel(meta["label"], fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        if inv_idx == 0:
            PYCSAMT_STATION_RENDERING.apply(
                ax,
                np.arange(n_st, dtype=float),
                all_st,
                preset="pseudosection",
                xlim=(-0.5, n_st - 0.5),
            )
        cax = ax.inset_axes([1.01, 0.05, 0.015, 0.9])
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=6)

    fig.suptitle(
        title or "Z rotation-invariants section",
        fontsize=10,
        fontweight="bold",
    )
    return fig
