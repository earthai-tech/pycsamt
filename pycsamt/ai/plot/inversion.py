# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Publication-quality visualisation for 2-D EM inversion results.

:func:`plot_inversion_result_2d`
    Multi-panel summary figure for :class:`~pycsamt.ai.inversion.EMInverter2D`
    (U-Net) output.  Up to five panels are arranged in two rows:

    * **(a) True model** — shown when ``log_true`` is supplied.
    * **(b) Predicted model** — always shown.
    * **(c) Signed misfit** ``Δlog₁₀ρ`` — shown when ``show_misfit=True`` and
      ``log_true`` is available.
    * **(d) Training convergence** — shown when ``train_loss`` is supplied.
    * **(e) Per-station RMSE** — shown when ``show_rmse=True``.

    Every aspect of the figure is controllable through keyword arguments.
    Pre-built :class:`~matplotlib.axes.Axes` can be injected via the ``axes``
    parameter for embedding in a larger composite figure.

Usage
-----
>>> from pycsamt.ai.plot import plot_inversion_result_2d
>>> fig = plot_inversion_result_2d(
...     log_pred,
...     log_true=log_true,
...     depths=depths,
...     stations=stations,
...     train_loss=train_loss,
...     val_loss=val_loss,
...     fault_positions=[9.75],
...     suptitle="EMInverter2D — profile L22PLT",
... )
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from ._style import (
    EMStyle, EM_COLORS, EM_CMAPS, EM_FIGSIZE, StationTickConfig,
)

__all__ = ["plot_inversion_result_2d"]

# ── default figsize table (n_top columns → width, fixed height per row) ─────
_TOP_WIDTHS  = {1: 6.0,  2: 10.0, 3: 14.0}
_TOP_HEIGHT  = 4.5
_BOT_HEIGHT  = 2.5


# ─────────────────────────────────────────────────────────────────────────────
# Private helpers
# ─────────────────────────────────────────────────────────────────────────────

def _cell_edges(centres: np.ndarray) -> np.ndarray:
    """Return (n+1,) cell-boundary array from (n,) cell-centre array."""
    if len(centres) < 2:
        half = 0.5
        return np.array([centres[0] - half, centres[0] + half])
    d = np.diff(centres)
    mids = centres[:-1] + d / 2.0
    lo   = centres[0]  - d[0]  / 2.0
    hi   = centres[-1] + d[-1] / 2.0
    return np.concatenate([[lo], mids, [hi]])


def _safe_vlim(
    arr: np.ndarray,
    pct_lo: float = 2.0,
    pct_hi: float = 98.0,
) -> Tuple[float, float]:
    lo = float(np.nanpercentile(arr, pct_lo))
    hi = float(np.nanpercentile(arr, pct_hi))
    if np.isclose(lo, hi):
        lo, hi = lo - 0.5, hi + 0.5
    return lo, hi


def _draw_section(
    ax: Axes,
    data: np.ndarray,               # (n_depth, n_stations)
    s_edges: np.ndarray,
    d_edges: np.ndarray,
    vmin: float,
    vmax: float,
    cmap: str,
    title: str,
    xlabel: str,
    ylabel: str,
    show_sites: bool,
    site_marker: str,
    site_color: str,
    site_ms: float,
    stations: np.ndarray,
    station_labels: Optional[List[str]],
    tick_cfg: StationTickConfig,
) -> Any:
    """Render one 2-D section panel; return the pcolormesh mappable."""
    im = ax.pcolormesh(
        s_edges, d_edges, data,
        cmap=cmap, vmin=vmin, vmax=vmax, shading="flat",
    )
    ax.invert_yaxis()
    if title:
        ax.set_title(title, fontsize=9, pad=4)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7)

    if show_sites:
        ax.plot(
            stations,
            np.full(len(stations), float(d_edges[0])),
            site_marker, ms=site_ms, color=site_color, zorder=5, clip_on=False,
        )

    if station_labels is not None:
        tick_cfg.apply(ax, stations, station_labels, xlabel=xlabel)
    return im


def _add_faults(
    ax: Axes,
    fault_positions: List[float],
    color: str,
    lw: float,
    ls: str,
    alpha: float,
    label_fontsize: float,
    y_label_frac: float = 0.92,
) -> None:
    ylim = ax.get_ylim()
    y_label = ylim[0] + (ylim[1] - ylim[0]) * (1.0 - y_label_frac)
    for x in fault_positions:
        ax.axvline(x, color=color, lw=lw, ls=ls, alpha=alpha, zorder=3)
        ax.text(
            x + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.008,
            y_label, "Fault",
            fontsize=label_fontsize, color=color,
            rotation=90, va="bottom", alpha=alpha,
        )


def _add_annotations(ax: Axes, annotations: List[dict]) -> None:
    """
    Apply a list of annotation dicts to *ax*.

    Each dict is forwarded to ``ax.annotate`` with sensible defaults.
    Required keys: ``text``, ``xy``.  Optional: ``xytext``, ``fontsize``,
    ``color``, ``arrowprops``, ``bbox``.
    """
    _default_arrow = dict(arrowstyle="->", color="white", lw=0.8)
    _default_bbox  = dict(fc="#1a1a1a", ec="none", alpha=0.55, pad=2)
    for ann in annotations:
        kw = dict(
            fontsize  = ann.get("fontsize", 7),
            color     = ann.get("color", "white"),
            arrowprops= ann.get("arrowprops", _default_arrow),
            bbox      = ann.get("bbox", _default_bbox),
        )
        if "xytext" in ann:
            kw["xytext"] = ann["xytext"]
        ax.annotate(ann["text"], xy=ann["xy"], **kw)


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_inversion_result_2d(
    log_pred: np.ndarray,
    *,
    # ── optional true model ──────────────────────────────────────────────────
    log_true: Optional[np.ndarray] = None,
    # ── grid geometry ────────────────────────────────────────────────────────
    depths: Optional[np.ndarray] = None,
    stations: Optional[np.ndarray] = None,
    station_labels: Optional[List[str]] = None,
    depth_max: float = 1500.0,
    station_spacing: float = 0.5,
    # ── resistivity colour scale ─────────────────────────────────────────────
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    cmap: Optional[str] = None,
    vmin_percentile: float = 2.0,
    vmax_percentile: float = 98.0,
    # ── misfit panel ─────────────────────────────────────────────────────────
    show_misfit: bool = True,
    misfit_cmap: Optional[str] = None,
    misfit_vlim: Optional[float] = None,
    misfit_percentile: float = 98.0,
    # ── training convergence panel ───────────────────────────────────────────
    show_convergence: bool = True,
    train_loss: Optional[np.ndarray] = None,
    val_loss: Optional[np.ndarray] = None,
    convergence_epochs: Optional[np.ndarray] = None,
    convergence_log_scale: bool = True,
    convergence_target: Optional[float] = None,
    convergence_target_color: str = "#aaaaaa",
    convergence_target_ls: str = ":",
    convergence_target_lw: float = 0.8,
    convergence_target_label: str = "target",
    # ── per-station RMSE panel ───────────────────────────────────────────────
    show_rmse: bool = True,
    rmse: Optional[np.ndarray] = None,
    rmse_threshold: float = 0.10,
    rmse_bar_color: Optional[str] = None,
    rmse_error_color: Optional[str] = None,
    rmse_threshold_color: Optional[str] = None,
    rmse_threshold_ls: str = "--",
    rmse_threshold_lw: float = 1.1,
    rmse_threshold_label: str = "target RMSE",
    # ── site markers ─────────────────────────────────────────────────────────
    show_sites: bool = True,
    site_marker: str = "v",
    site_color: str = "#333333",
    site_ms: float = 3.5,
    # ── fault overlays ───────────────────────────────────────────────────────
    fault_positions: Optional[List[float]] = None,
    fault_color: str = "#eeeeee",
    fault_lw: float = 0.9,
    fault_ls: str = "--",
    fault_alpha: float = 0.75,
    fault_label_fontsize: float = 6.5,
    # ── geology annotations (on leftmost section panel) ──────────────────────
    annotations: Optional[List[dict]] = None,
    # ── colourbar labels ─────────────────────────────────────────────────────
    colorbar_label: str = r"$\log_{10}(\rho)$  (Ω·m)",
    misfit_colorbar_label: str = r"$\Delta\log_{10}(\rho)$",
    # ── panel titles ─────────────────────────────────────────────────────────
    title_true: str = "(a) True model",
    title_pred: Optional[str] = None,
    title_misfit: str = r"(c) Misfit  $\Delta\log_{10}\rho$",
    title_convergence: str = "(d) Training convergence",
    title_rmse: str = "(e) Per-station RMSE",
    xlabel: str = "Distance (km)",
    ylabel: str = "Depth (m)",
    suptitle: str = "",
    # ── layout ───────────────────────────────────────────────────────────────
    figsize: Optional[Tuple[float, float]] = None,
    top_height_ratio: float = 2.5,
    bottom_height_ratio: float = 1.2,
    hspace: float = 0.40,
    wspace: float = 0.36,
    axes: Optional[Dict[str, Axes]] = None,
    # ── station tick spacing ─────────────────────────────────────────────────
    tick_every: Union[int, str] = "auto",
    tick_label_rotation: float = 0.0,
    tick_fontsize: int = 7,
    station_tick_config: Optional[StationTickConfig] = None,
) -> Figure:
    """
    Multi-panel 2-D inversion result figure.

    Parameters
    ----------
    log_pred : ndarray, shape (n_depth, n_stations)
        Predicted :math:`\\log_{10}(\\rho)` section from the network.
    log_true : ndarray (n_depth, n_stations) or None
        Ground-truth :math:`\\log_{10}(\\rho)` section.  When supplied,
        panel **(a)** is drawn and the misfit panel **(c)** is enabled.
    depths : ndarray (n_depth,) or None
        Depth node positions in metres.  When ``None``, equally spaced
        between 0 and ``depth_max``.
    stations : ndarray (n_stations,) or None
        Station positions.  When ``None``, spaced by ``station_spacing``.
    station_labels : list of str or None
        Custom labels for the x-axis tick marks (e.g. station names).
        Uses :class:`~pycsamt.ai.plot._style.StationTickConfig` to prevent
        overlap automatically.
    depth_max : float
        Maximum depth when ``depths`` is ``None``.
    station_spacing : float
        Station interval when ``stations`` is ``None``.
    vmin, vmax : float or None
        Colour-scale limits for the resistivity panels.  Auto-computed
        from ``vmin_percentile`` / ``vmax_percentile`` of *log_true* when
        ``None``.
    cmap : str or None
        Colourmap for resistivity panels.  Defaults to
        ``EM_CMAPS['resistivity']`` (``'RdYlBu_r'``).
    vmin_percentile, vmax_percentile : float
        Percentile clipping for automatic colour-scale calculation.
    show_misfit : bool
        Draw the signed-misfit panel **(c)**.  Silently ignored when
        ``log_true`` is ``None``.
    misfit_cmap : str or None
        Colourmap for the misfit panel.  Defaults to
        ``EM_CMAPS['misfit']`` (``'RdBu_r'``).
    misfit_vlim : float or None
        Symmetric limit for the misfit colorbar.  Auto-computed from
        ``misfit_percentile`` when ``None``.
    misfit_percentile : float
        Percentile used to clip the symmetric misfit colour scale.
    show_convergence : bool
        Draw the training-convergence panel **(d)**.  Silently ignored
        when ``train_loss`` is ``None``.
    train_loss : ndarray (n_epochs,) or None
        Per-epoch training loss values.
    val_loss : ndarray (n_epochs,) or None
        Per-epoch validation loss values.  Omitted when ``None``.
    convergence_epochs : ndarray (n_epochs,) or None
        Epoch indices.  Defaults to ``1 … len(train_loss)``.
    convergence_log_scale : bool
        Use a log y-axis on the convergence panel.
    convergence_target : float or None
        Draw a horizontal reference line at this loss value.
    convergence_target_color, convergence_target_ls, convergence_target_lw : str / float
        Style of the target reference line.
    convergence_target_label : str
        Legend label for the target line.
    show_rmse : bool
        Draw the per-station RMSE bar chart **(e)**.
    rmse : ndarray (n_stations,) or None
        Pre-computed per-station RMSE.  Computed from ``log_pred - log_true``
        when ``None`` and ``log_true`` is available; otherwise zeros.
    rmse_threshold : float
        Dashed reference line on the RMSE panel.
    rmse_bar_color, rmse_error_color : str or None
        Bar colours for stations below / above ``rmse_threshold``.
    rmse_threshold_color : str or None
        Colour of the RMSE threshold line.
    rmse_threshold_ls, rmse_threshold_lw : str / float
        Line style / width of the RMSE threshold line.
    rmse_threshold_label : str
        Legend label for the threshold line.
    show_sites : bool
        Draw site-marker triangles at the surface of each section panel.
    site_marker, site_color, site_ms : str / float
        Marker style, colour, and size.
    fault_positions : list of float or None
        x-positions (same units as ``stations``) where fault lines are
        drawn on every section panel.
    fault_color, fault_lw, fault_ls, fault_alpha : str / float
        Style of the fault overlay.
    fault_label_fontsize : float
        Font size of the ``"Fault"`` text label.
    annotations : list of dict or None
        Geology annotation dicts applied to the leftmost section panel.
        Each dict must contain ``text`` and ``xy``; optional keys are
        ``xytext``, ``fontsize``, ``color``, ``arrowprops``, ``bbox``.
    colorbar_label : str
        Label for the shared resistivity colorbar.
    misfit_colorbar_label : str
        Label for the misfit colorbar.
    title_true, title_pred, title_misfit : str
        Panel titles.  ``title_pred`` defaults to ``"(a) Predicted model"``
        when ``log_true`` is ``None``, else ``"(b) Predicted model"``.
    title_convergence, title_rmse : str
        Panel titles for the bottom-row panels.
    xlabel, ylabel : str
        Axis labels for every section panel.
    suptitle : str
        Figure-level super-title.
    figsize : (width, height) or None
        Figure size.  Auto-derived from the number of top panels when
        ``None``.
    top_height_ratio, bottom_height_ratio : float
        Relative heights of the top and bottom grid rows.
    hspace, wspace : float
        GridSpec spacing parameters.
    axes : dict {str → Axes} or None
        Pre-built axes keyed by ``"true"``, ``"pred"``, ``"misfit"``,
        ``"convergence"``, ``"rmse"``.  When supplied, figure and
        GridSpec creation is skipped entirely and the function draws
        into the provided axes.
    tick_every : int or ``"auto"``
        Station-tick spacing for ``station_labels`` x-axis.  Forwarded
        to :class:`~pycsamt.ai.plot._style.StationTickConfig`.
    tick_label_rotation : float
        Tick-label rotation for station labels.
    tick_fontsize : int
        Tick-label font size for station labels.
    station_tick_config : StationTickConfig or None
        Pre-built tick config; overrides ``tick_every`` / rotation /
        fontsize when provided.

    Returns
    -------
    fig : Figure

    Examples
    --------
    Minimal (prediction only, no training curves):

    >>> fig = plot_inversion_result_2d(log_pred, depths=depths, stations=stations)

    Full five-panel result:

    >>> fig = plot_inversion_result_2d(
    ...     log_pred,
    ...     log_true=log_true,
    ...     depths=depths,
    ...     stations=stations,
    ...     train_loss=train_loss,
    ...     val_loss=val_loss,
    ...     fault_positions=[9.75],
    ...     annotations=[{
    ...         "text": "Conductive\\nbasin",
    ...         "xy": (6.0, 450.0),
    ...         "xytext": (3.5, 270.0),
    ...     }],
    ...     suptitle="EMInverter2D — profile L22PLT",
    ... )
    """
    log_pred = np.asarray(log_pred, dtype=float)
    n_depth, n_stations = log_pred.shape

    if log_true is not None:
        log_true = np.asarray(log_true, dtype=float)

    # ── grid geometry ────────────────────────────────────────────────────────
    if depths is None:
        depths = np.linspace(0.0, depth_max, n_depth)
    if stations is None:
        stations = np.arange(n_stations) * float(station_spacing)

    d_edges = _cell_edges(depths)
    s_edges = _cell_edges(stations)

    # ── colour limits ────────────────────────────────────────────────────────
    ref = log_true if log_true is not None else log_pred
    if vmin is None:
        vmin, _ = _safe_vlim(ref, vmin_percentile, 50.0)
        vmin = float(np.nanpercentile(ref, vmin_percentile))
    if vmax is None:
        vmax = float(np.nanpercentile(ref, vmax_percentile))
    if cmap is None:
        cmap = EM_CMAPS["resistivity"]

    # ── misfit ───────────────────────────────────────────────────────────────
    _show_misfit = show_misfit and (log_true is not None)
    misfit = (log_pred - log_true) if log_true is not None else None
    if misfit_cmap is None:
        misfit_cmap = EM_CMAPS["misfit"]
    if _show_misfit:
        if misfit_vlim is None:
            misfit_vlim = float(np.nanpercentile(np.abs(misfit), misfit_percentile))

    # ── RMSE ─────────────────────────────────────────────────────────────────
    _show_rmse = show_rmse and (train_loss is not None or rmse is not None
                                or log_true is not None)
    if rmse is None:
        if log_true is not None:
            rmse = np.sqrt(np.nanmean((log_pred - log_true) ** 2, axis=0))
        else:
            _show_rmse = False
    if rmse is not None:
        rmse = np.asarray(rmse, dtype=float)

    if rmse_bar_color is None:
        rmse_bar_color = EM_COLORS["primary"]
    if rmse_error_color is None:
        rmse_error_color = EM_COLORS["error"]
    if rmse_threshold_color is None:
        rmse_threshold_color = EM_COLORS["secondary"]

    # ── convergence ──────────────────────────────────────────────────────────
    _show_conv = show_convergence and (train_loss is not None)
    if train_loss is not None:
        train_loss = np.asarray(train_loss, dtype=float)
    if val_loss is not None:
        val_loss = np.asarray(val_loss, dtype=float)
    if convergence_epochs is None and train_loss is not None:
        convergence_epochs = np.arange(1, len(train_loss) + 1)

    # ── panel list & auto title for pred ────────────────────────────────────
    has_true = log_true is not None
    panels_top: List[str] = []
    if has_true:
        panels_top.append("true")
    panels_top.append("pred")
    if _show_misfit:
        panels_top.append("misfit")
    n_top = len(panels_top)

    has_bottom = _show_conv or _show_rmse
    panels_bot: List[str] = []
    if _show_conv:
        panels_bot.append("convergence")
    if _show_rmse:
        panels_bot.append("rmse")

    if title_pred is None:
        letter = chr(ord("a") + panels_top.index("pred"))
        title_pred = f"({letter}) Predicted model"
    # override title letters to match position
    _titles = {
        "true":        title_true,
        "pred":        title_pred,
        "misfit":      title_misfit,
        "convergence": title_convergence,
        "rmse":        title_rmse,
    }

    # ── station tick config ──────────────────────────────────────────────────
    if station_tick_config is None:
        station_tick_config = StationTickConfig(
            every=tick_every,
            rotation=tick_label_rotation,
            fontsize=tick_fontsize,
        )

    # ── figure / axes creation ───────────────────────────────────────────────
    using_external_axes = axes is not None
    ax: Dict[str, Axes] = {}

    if using_external_axes:
        ax = dict(axes)
        fig = next(iter(ax.values())).get_figure()
    else:
        if figsize is None:
            w = _TOP_WIDTHS.get(n_top, n_top * 5.0)
            h = _TOP_HEIGHT + ((_BOT_HEIGHT + 0.5) if has_bottom else 0.0)
            figsize = (w, h)

        fig = plt.figure(figsize=figsize)

        if has_bottom:
            gs_outer = fig.add_gridspec(
                2, 1,
                height_ratios=[top_height_ratio, bottom_height_ratio],
                hspace=hspace,
            )
            gs_top = gs_outer[0].subgridspec(1, n_top, wspace=wspace)
            gs_bot = gs_outer[1].subgridspec(1, n_top, wspace=wspace)
        else:
            gs_top = fig.add_gridspec(1, n_top, wspace=wspace)
            gs_bot = None

        for i, key in enumerate(panels_top):
            ax[key] = fig.add_subplot(gs_top[0, i])

        if gs_bot is not None:
            n_bot = len(panels_bot)
            if n_bot == 2:
                ax["convergence"] = fig.add_subplot(
                    gs_bot[0, : max(1, n_top - 1)]
                )
                ax["rmse"] = fig.add_subplot(gs_bot[0, n_top - 1])
            elif n_bot == 1:
                ax[panels_bot[0]] = fig.add_subplot(gs_bot[0, :])

    # ── draw section panels ──────────────────────────────────────────────────
    leftmost_sec_key = panels_top[0]   # where annotations + faults go
    im_res = None

    for key in ("true", "pred"):
        if key not in ax:
            continue
        data = log_true if key == "true" else log_pred
        im = _draw_section(
            ax[key], data, s_edges, d_edges,
            vmin, vmax, cmap,
            _titles[key], xlabel, ylabel,
            show_sites, site_marker, site_color, site_ms,
            stations, station_labels, station_tick_config,
        )
        im_res = im

    # misfit panel
    im_mis = None
    if "misfit" in ax:
        im_mis = _draw_section(
            ax["misfit"], misfit, s_edges, d_edges,
            -misfit_vlim, misfit_vlim, misfit_cmap,
            _titles["misfit"], xlabel, ylabel,
            show_sites, site_marker, site_color, site_ms,
            stations, station_labels, station_tick_config,
        )

    # ── shared resistivity colorbar ──────────────────────────────────────────
    if im_res is not None and not using_external_axes:
        res_axes = [ax[k] for k in ("true", "pred") if k in ax]
        cb_r = fig.colorbar(
            im_res, ax=res_axes,
            shrink=0.88, pad=0.02, fraction=0.022,
        )
        cb_r.set_label(colorbar_label, fontsize=8)
        cb_r.ax.tick_params(labelsize=7)

    if im_mis is not None and not using_external_axes:
        cb_m = fig.colorbar(
            im_mis, ax=ax["misfit"],
            shrink=0.88, pad=0.02, fraction=0.044,
        )
        cb_m.set_label(misfit_colorbar_label, fontsize=8)
        cb_m.ax.tick_params(labelsize=7)

    # ── fault overlays ───────────────────────────────────────────────────────
    if fault_positions:
        for key in ("true", "pred", "misfit"):
            if key in ax:
                _add_faults(
                    ax[key], fault_positions,
                    fault_color, fault_lw, fault_ls,
                    fault_alpha, fault_label_fontsize,
                )

    # ── geology annotations ──────────────────────────────────────────────────
    if annotations and leftmost_sec_key in ax:
        _add_annotations(ax[leftmost_sec_key], annotations)

    # ── convergence panel ────────────────────────────────────────────────────
    if "convergence" in ax:
        ax_c = ax["convergence"]
        plot_fn = ax_c.semilogy if convergence_log_scale else ax_c.plot
        plot_fn(
            convergence_epochs, train_loss,
            color=EM_COLORS["primary"], lw=1.6, label="Train loss",
        )
        if val_loss is not None:
            plot_fn(
                convergence_epochs, val_loss,
                color=EM_COLORS["secondary"], lw=1.6, ls="--", label="Val loss",
            )
        if convergence_target is not None:
            ax_c.axhline(
                convergence_target,
                color=convergence_target_color,
                lw=convergence_target_lw,
                ls=convergence_target_ls,
                label=convergence_target_label,
            )
            xlim = ax_c.get_xlim()
            ax_c.text(
                xlim[1] * 0.96, convergence_target * 1.08,
                convergence_target_label,
                fontsize=6.5, color=convergence_target_color,
                ha="right", va="bottom",
            )
        ax_c.set_xlabel("Epoch", fontsize=8)
        ax_c.set_ylabel("MSE (log scale)" if convergence_log_scale else "MSE", fontsize=8)
        ax_c.set_title(_titles["convergence"], fontsize=9)
        ax_c.legend(fontsize=7, framealpha=0.7)
        ax_c.tick_params(labelsize=7)
        if convergence_epochs is not None:
            ax_c.set_xlim(convergence_epochs[0], convergence_epochs[-1])

    # ── RMSE bar chart ────────────────────────────────────────────────────────
    if "rmse" in ax and rmse is not None:
        ax_r = ax["rmse"]
        bar_c = [
            rmse_error_color if v > rmse_threshold else rmse_bar_color
            for v in rmse
        ]
        dx = float(np.diff(s_edges).mean())
        ax_r.bar(
            stations, rmse, width=dx * 0.65,
            color=bar_c, edgecolor="none", alpha=0.85,
        )
        ax_r.axhline(
            rmse_threshold,
            color=rmse_threshold_color,
            lw=rmse_threshold_lw,
            ls=rmse_threshold_ls,
            label=rmse_threshold_label,
        )
        ax_r.set_xlabel(xlabel, fontsize=8)
        ax_r.set_ylabel(r"RMSE  ($\Delta\log_{10}\rho$)", fontsize=8)
        ax_r.set_title(_titles["rmse"], fontsize=9)
        ax_r.legend(fontsize=7, framealpha=0.7)
        ax_r.tick_params(labelsize=7)
        ax_r.set_xlim(s_edges[0], s_edges[-1])

    # ── suptitle ─────────────────────────────────────────────────────────────
    if suptitle:
        fig.suptitle(suptitle, fontsize=10, fontweight="bold", y=1.002)

    return fig
