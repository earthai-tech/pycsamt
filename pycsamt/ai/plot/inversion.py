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

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ._style import (
    EM_CMAPS,
    EM_COLORS,
    EMStyle,
    StationTickConfig,
)

__all__ = ["plot_inversion_result_2d", "plot_inversion_result_3d"]

# ── default figsize table (n_top columns → width, fixed height per row) ─────
_TOP_WIDTHS = {1: 6.0, 2: 10.0, 3: 14.0}
_TOP_HEIGHT = 4.5
_BOT_HEIGHT = 2.5


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
    lo = centres[0] - d[0] / 2.0
    hi = centres[-1] + d[-1] / 2.0
    return np.concatenate([[lo], mids, [hi]])


def _safe_vlim(
    arr: np.ndarray,
    pct_lo: float = 2.0,
    pct_hi: float = 98.0,
) -> tuple[float, float]:
    lo = float(np.nanpercentile(arr, pct_lo))
    hi = float(np.nanpercentile(arr, pct_hi))
    if np.isclose(lo, hi):
        lo, hi = lo - 0.5, hi + 0.5
    return lo, hi


def _draw_section(
    ax: Axes,
    data: np.ndarray,  # (n_depth, n_stations)
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
    station_labels: list[str] | None,
    tick_cfg: StationTickConfig,
) -> Any:
    """Render one 2-D section panel; return the pcolormesh mappable."""
    im = ax.pcolormesh(
        s_edges,
        d_edges,
        data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="flat",
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
            site_marker,
            ms=site_ms,
            color=site_color,
            zorder=5,
            clip_on=False,
        )

    if station_labels is not None:
        tick_cfg.apply(ax, stations, station_labels, xlabel=xlabel)
    return im


def _add_faults(
    ax: Axes,
    fault_positions: list[float],
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
            y_label,
            "Fault",
            fontsize=label_fontsize,
            color=color,
            rotation=90,
            va="bottom",
            alpha=alpha,
        )


def _add_annotations(ax: Axes, annotations: list[dict]) -> None:
    """
    Apply a list of annotation dicts to *ax*.

    Each dict is forwarded to ``ax.annotate`` with sensible defaults.
    Required keys: ``text``, ``xy``.  Optional: ``xytext``, ``fontsize``,
    ``color``, ``arrowprops``, ``bbox``.
    """
    _default_arrow = dict(arrowstyle="->", color="white", lw=0.8)
    _default_bbox = dict(fc="#1a1a1a", ec="none", alpha=0.55, pad=2)
    for ann in annotations:
        kw = dict(
            fontsize=ann.get("fontsize", 7),
            color=ann.get("color", "white"),
            arrowprops=ann.get("arrowprops", _default_arrow),
            bbox=ann.get("bbox", _default_bbox),
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
    log_true: np.ndarray | None = None,
    # ── grid geometry ────────────────────────────────────────────────────────
    depths: np.ndarray | None = None,
    stations: np.ndarray | None = None,
    station_labels: list[str] | None = None,
    depth_max: float = 1500.0,
    station_spacing: float = 0.5,
    # ── resistivity colour scale ─────────────────────────────────────────────
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | None = None,
    vmin_percentile: float = 2.0,
    vmax_percentile: float = 98.0,
    # ── misfit panel ─────────────────────────────────────────────────────────
    show_misfit: bool = True,
    misfit_cmap: str | None = None,
    misfit_vlim: float | None = None,
    misfit_percentile: float = 98.0,
    # ── training convergence panel ───────────────────────────────────────────
    show_convergence: bool = True,
    train_loss: np.ndarray | None = None,
    val_loss: np.ndarray | None = None,
    convergence_epochs: np.ndarray | None = None,
    convergence_log_scale: bool = True,
    convergence_target: float | None = None,
    convergence_target_color: str = "#aaaaaa",
    convergence_target_ls: str = ":",
    convergence_target_lw: float = 0.8,
    convergence_target_label: str = "target",
    # ── per-station RMSE panel ───────────────────────────────────────────────
    show_rmse: bool = True,
    rmse: np.ndarray | None = None,
    rmse_threshold: float = 0.10,
    rmse_bar_color: str | None = None,
    rmse_error_color: str | None = None,
    rmse_threshold_color: str | None = None,
    rmse_threshold_ls: str = "--",
    rmse_threshold_lw: float = 1.1,
    rmse_threshold_label: str = "target RMSE",
    # ── site markers ─────────────────────────────────────────────────────────
    show_sites: bool = True,
    site_marker: str = "v",
    site_color: str = "#333333",
    site_ms: float = 3.5,
    # ── fault overlays ───────────────────────────────────────────────────────
    fault_positions: list[float] | None = None,
    fault_color: str = "#eeeeee",
    fault_lw: float = 0.9,
    fault_ls: str = "--",
    fault_alpha: float = 0.75,
    fault_label_fontsize: float = 6.5,
    # ── geology annotations (on leftmost section panel) ──────────────────────
    annotations: list[dict] | None = None,
    # ── colourbar labels ─────────────────────────────────────────────────────
    colorbar_label: str = r"$\log_{10}(\rho)$  (Ω·m)",
    misfit_colorbar_label: str = r"$\Delta\log_{10}(\rho)$",
    # ── panel titles ─────────────────────────────────────────────────────────
    title_true: str = "(a) True model",
    title_pred: str | None = None,
    title_misfit: str = r"(c) Misfit  $\Delta\log_{10}\rho$",
    title_convergence: str = "(d) Training convergence",
    title_rmse: str = "(e) Per-station RMSE",
    xlabel: str = "Distance (km)",
    ylabel: str = "Depth (m)",
    suptitle: str = "",
    # ── layout ───────────────────────────────────────────────────────────────
    figsize: tuple[float, float] | None = None,
    top_height_ratio: float = 2.5,
    bottom_height_ratio: float = 1.2,
    hspace: float = 0.40,
    wspace: float = 0.36,
    axes: dict[str, Axes] | None = None,
    # ── station tick spacing ─────────────────────────────────────────────────
    tick_every: int | str = "auto",
    tick_label_rotation: float = 0.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
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

    >>> fig = plot_inversion_result_2d(
    ...     log_pred, depths=depths, stations=stations
    ... )

    Full five-panel result:

    >>> fig = plot_inversion_result_2d(
    ...     log_pred,
    ...     log_true=log_true,
    ...     depths=depths,
    ...     stations=stations,
    ...     train_loss=train_loss,
    ...     val_loss=val_loss,
    ...     fault_positions=[9.75],
    ...     annotations=[
    ...         {
    ...             "text": "Conductive\\nbasin",
    ...             "xy": (6.0, 450.0),
    ...             "xytext": (3.5, 270.0),
    ...         }
    ...     ],
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
    _show_rmse = show_rmse and (
        train_loss is not None or rmse is not None or log_true is not None
    )
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
    panels_top: list[str] = []
    if has_true:
        panels_top.append("true")
    panels_top.append("pred")
    if _show_misfit:
        panels_top.append("misfit")
    n_top = len(panels_top)

    has_bottom = _show_conv or _show_rmse
    panels_bot: list[str] = []
    if _show_conv:
        panels_bot.append("convergence")
    if _show_rmse:
        panels_bot.append("rmse")

    if title_pred is None:
        letter = chr(ord("a") + panels_top.index("pred"))
        title_pred = f"({letter}) Predicted model"
    # override title letters to match position
    _titles = {
        "true": title_true,
        "pred": title_pred,
        "misfit": title_misfit,
        "convergence": title_convergence,
        "rmse": title_rmse,
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
    ax: dict[str, Axes] = {}

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
                2,
                1,
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
                ax["convergence"] = fig.add_subplot(gs_bot[0, : max(1, n_top - 1)])
                ax["rmse"] = fig.add_subplot(gs_bot[0, n_top - 1])
            elif n_bot == 1:
                ax[panels_bot[0]] = fig.add_subplot(gs_bot[0, :])

    # ── draw section panels ──────────────────────────────────────────────────
    leftmost_sec_key = panels_top[0]  # where annotations + faults go
    im_res = None

    for key in ("true", "pred"):
        if key not in ax:
            continue
        data = log_true if key == "true" else log_pred
        im = _draw_section(
            ax[key],
            data,
            s_edges,
            d_edges,
            vmin,
            vmax,
            cmap,
            _titles[key],
            xlabel,
            ylabel,
            show_sites,
            site_marker,
            site_color,
            site_ms,
            stations,
            station_labels,
            station_tick_config,
        )
        im_res = im

    # misfit panel
    im_mis = None
    if "misfit" in ax:
        im_mis = _draw_section(
            ax["misfit"],
            misfit,
            s_edges,
            d_edges,
            -misfit_vlim,
            misfit_vlim,
            misfit_cmap,
            _titles["misfit"],
            xlabel,
            ylabel,
            show_sites,
            site_marker,
            site_color,
            site_ms,
            stations,
            station_labels,
            station_tick_config,
        )

    # ── shared resistivity colorbar ──────────────────────────────────────────
    if im_res is not None and not using_external_axes:
        res_axes = [ax[k] for k in ("true", "pred") if k in ax]
        cb_r = fig.colorbar(
            im_res,
            ax=res_axes,
            shrink=0.88,
            pad=0.02,
            fraction=0.022,
        )
        cb_r.set_label(colorbar_label, fontsize=8)
        cb_r.ax.tick_params(labelsize=7)

    if im_mis is not None and not using_external_axes:
        cb_m = fig.colorbar(
            im_mis,
            ax=ax["misfit"],
            shrink=0.88,
            pad=0.02,
            fraction=0.044,
        )
        cb_m.set_label(misfit_colorbar_label, fontsize=8)
        cb_m.ax.tick_params(labelsize=7)

    # ── fault overlays ───────────────────────────────────────────────────────
    if fault_positions:
        for key in ("true", "pred", "misfit"):
            if key in ax:
                _add_faults(
                    ax[key],
                    fault_positions,
                    fault_color,
                    fault_lw,
                    fault_ls,
                    fault_alpha,
                    fault_label_fontsize,
                )

    # ── geology annotations ──────────────────────────────────────────────────
    if annotations and leftmost_sec_key in ax:
        _add_annotations(ax[leftmost_sec_key], annotations)

    # ── convergence panel ────────────────────────────────────────────────────
    if "convergence" in ax:
        ax_c = ax["convergence"]
        plot_fn = ax_c.semilogy if convergence_log_scale else ax_c.plot
        plot_fn(
            convergence_epochs,
            train_loss,
            color=EM_COLORS["primary"],
            lw=1.6,
            label="Train loss",
        )
        if val_loss is not None:
            plot_fn(
                convergence_epochs,
                val_loss,
                color=EM_COLORS["secondary"],
                lw=1.6,
                ls="--",
                label="Val loss",
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
                xlim[1] * 0.96,
                convergence_target * 1.08,
                convergence_target_label,
                fontsize=6.5,
                color=convergence_target_color,
                ha="right",
                va="bottom",
            )
        ax_c.set_xlabel("Epoch", fontsize=8)
        ax_c.set_ylabel(
            "MSE (log scale)" if convergence_log_scale else "MSE", fontsize=8
        )
        ax_c.set_title(_titles["convergence"], fontsize=9)
        ax_c.legend(fontsize=7, framealpha=0.7)
        ax_c.tick_params(labelsize=7)
        if convergence_epochs is not None:
            ax_c.set_xlim(convergence_epochs[0], convergence_epochs[-1])

    # ── RMSE bar chart ────────────────────────────────────────────────────────
    if "rmse" in ax and rmse is not None:
        ax_r = ax["rmse"]
        bar_c = [
            rmse_error_color if v > rmse_threshold else rmse_bar_color for v in rmse
        ]
        dx = float(np.diff(s_edges).mean())
        ax_r.bar(
            stations,
            rmse,
            width=dx * 0.65,
            color=bar_c,
            edgecolor="none",
            alpha=0.85,
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


# ═════════════════════════════════════════════════════════════════════════════
# Helpers – 3-D panels
# ═════════════════════════════════════════════════════════════════════════════


def _per_station_to_volume(
    station_models: np.ndarray,  # (n_stations, n_depth)
    station_xy: np.ndarray,  # (n_stations, 2)
    x_coords: np.ndarray,  # (n_x,)
    y_coords: np.ndarray,  # (n_y,)
    method: str = "linear",
) -> np.ndarray:  # (n_depth, n_y, n_x)
    """
    Interpolate scattered per-station 1-D models onto a regular (y, x) grid.

    Requires :mod:`scipy`.  Falls back to nearest-neighbour if the
    requested *method* raises an error (e.g. too few points for cubic).
    """
    from scipy.interpolate import griddata

    _, n_depth = station_models.shape
    n_y, n_x = len(y_coords), len(x_coords)
    XX, YY = np.meshgrid(x_coords, y_coords)  # (n_y, n_x)
    pts = station_xy.astype(float)
    volume = np.empty((n_depth, n_y, n_x), dtype=float)
    for d in range(n_depth):
        vals = station_models[:, d]
        try:
            layer = griddata(pts, vals, (XX, YY), method=method, fill_value=np.nan)
        except Exception:
            layer = griddata(pts, vals, (XX, YY), method="nearest")
        volume[d] = layer
    return volume


def _draw_map_panel(
    ax: Axes,
    data: np.ndarray,  # (n_y, n_x)
    x_edges: np.ndarray,  # (n_x+1,)
    y_edges: np.ndarray,  # (n_y+1,)
    vmin: float,
    vmax: float,
    cmap: str,
    title: str,
    xlabel: str,
    ylabel: str,
    aspect: str = "auto",
) -> Any:
    """Render a horizontal depth-slice map; return the pcolormesh mappable."""
    im = ax.pcolormesh(
        x_edges,
        y_edges,
        data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="flat",
    )
    ax.set_aspect(aspect)
    if title:
        ax.set_title(title, fontsize=9, pad=4)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7)
    return im


def _draw_rmse_scatter(
    ax: Axes,
    station_xy: np.ndarray,  # (n_stations, 2)
    rmse: np.ndarray,  # (n_stations,)
    cmap: str,
    vmax: float,
    marker: str,
    ms: float,
    threshold: float | None,
    title: str,
    xlabel: str,
    ylabel: str,
    aspect: str = "auto",
) -> Any:
    """Scatter map of per-station RMSE; return scatter mappable."""
    sc = ax.scatter(
        station_xy[:, 0],
        station_xy[:, 1],
        c=rmse,
        cmap=cmap,
        vmin=0.0,
        vmax=vmax,
        s=ms,
        marker=marker,
        edgecolors="none",
        zorder=4,
    )
    if threshold is not None:
        hi = rmse > threshold
        if hi.any():
            ax.scatter(
                station_xy[hi, 0],
                station_xy[hi, 1],
                facecolors="none",
                edgecolors="#cc0000",
                s=ms * 1.6,
                lw=0.9,
                zorder=5,
                label=f"RMSE > {threshold:.2f}",
            )
    ax.set_aspect(aspect)
    if title:
        ax.set_title(title, fontsize=9, pad=4)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7)
    return sc


def _overlay_cut_lines(
    ax: Axes,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
    inline_idx: int | None,
    crossline_idx: int | None,
    color: str = "#ffffff",
    lw: float = 1.0,
    ls: str = "-",
    alpha: float = 0.75,
    fontsize: float = 6.5,
) -> None:
    """Draw inline / crossline position indicators on a map-view axes."""
    x0, x1 = float(x_coords[0]), float(x_coords[-1])
    y0, y1 = float(y_coords[0]), float(y_coords[-1])
    dx, dy = x1 - x0, y1 - y0

    if inline_idx is not None and 0 <= inline_idx < len(y_coords):
        yc = float(y_coords[inline_idx])
        ax.axhline(yc, color=color, lw=lw, ls=ls, alpha=alpha, zorder=3)
        ax.text(
            x0 + dx * 0.01,
            yc + dy * 0.015,
            "Inline",
            fontsize=fontsize,
            color=color,
            alpha=alpha,
            ha="left",
            va="bottom",
        )

    if crossline_idx is not None and 0 <= crossline_idx < len(x_coords):
        xc = float(x_coords[crossline_idx])
        ax.axvline(xc, color=color, lw=lw, ls=ls, alpha=alpha, zorder=3)
        ax.text(
            xc + dx * 0.012,
            y0 + dy * 0.03,
            "X-line",
            fontsize=fontsize,
            color=color,
            alpha=alpha,
            ha="left",
            va="bottom",
            rotation=90,
        )


def _draw_graph_edges(
    ax: Axes,
    station_xy: np.ndarray,  # (n_stations, 2)
    adjacency: np.ndarray,  # (n_stations, n_stations)
    color: str,
    lw: float,
    alpha: float,
    threshold: float = 1e-2,
) -> None:
    """Draw GCN graph edges between connected stations on a map axes."""
    n = len(station_xy)
    for i in range(n):
        for j in range(i + 1, n):
            if float(adjacency[i, j]) > threshold:
                ax.plot(
                    [station_xy[i, 0], station_xy[j, 0]],
                    [station_xy[i, 1], station_xy[j, 1]],
                    color=color,
                    lw=lw,
                    alpha=alpha,
                    zorder=2,
                )


def _render_convergence(
    ax: Axes,
    epochs: np.ndarray,
    train_loss: np.ndarray,
    val_loss: np.ndarray | None,
    log_scale: bool,
    target: float | None,
    target_color: str,
    target_ls: str,
    target_lw: float,
    target_label: str,
    title: str,
) -> None:
    """Draw training convergence curves on *ax*."""
    plot_fn = ax.semilogy if log_scale else ax.plot
    plot_fn(
        epochs,
        train_loss,
        color=EM_COLORS["primary"],
        lw=1.6,
        label="Train loss",
    )
    if val_loss is not None:
        plot_fn(
            epochs,
            val_loss,
            color=EM_COLORS["secondary"],
            lw=1.6,
            ls="--",
            label="Val loss",
        )
    if target is not None:
        ax.axhline(
            target,
            color=target_color,
            lw=target_lw,
            ls=target_ls,
            label=target_label,
        )
        xl = ax.get_xlim()
        ax.text(
            xl[1] * 0.96,
            target * 1.10,
            target_label,
            fontsize=6.5,
            color=target_color,
            ha="right",
            va="bottom",
        )
    ax.set_xlabel("Epoch", fontsize=8)
    ax.set_ylabel("MSE (log scale)" if log_scale else "MSE", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7, framealpha=0.7)
    ax.tick_params(labelsize=7)
    ax.set_xlim(epochs[0], epochs[-1])


# ═════════════════════════════════════════════════════════════════════════════
# plot_inversion_result_3d
# ═════════════════════════════════════════════════════════════════════════════


@EMStyle()
def plot_inversion_result_3d(
    log_pred: np.ndarray,
    *,
    # ── optional true model ──────────────────────────────────────────────────
    log_true: np.ndarray | None = None,
    # ── grid geometry ────────────────────────────────────────────────────────
    depths: np.ndarray | None = None,
    x_coords: np.ndarray | None = None,
    y_coords: np.ndarray | None = None,
    depth_max: float = 1500.0,
    x_spacing: float = 0.5,
    y_spacing: float = 0.5,
    # ── per-station → volume interpolation (when log_pred is 2-D) ────────────
    station_xy: np.ndarray | None = None,
    interp_method: str = "linear",
    interp_nx: int = 30,
    interp_ny: int = 30,
    # ── slice planes ─────────────────────────────────────────────────────────
    depth_slice_idx: int | None = None,
    depth_slice_m: float | None = None,
    inline_idx: int | None = None,
    crossline_idx: int | None = None,
    # ── panel visibility ─────────────────────────────────────────────────────
    show_depth_slice: bool = True,
    show_inline: bool = True,
    show_crossline: bool = True,
    show_misfit: bool = True,
    show_convergence: bool = True,
    show_rmse_map: bool = True,
    show_section_lines: bool = True,
    # ── resistivity colour scale ─────────────────────────────────────────────
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | None = None,
    vmin_percentile: float = 2.0,
    vmax_percentile: float = 98.0,
    map_aspect: str = "auto",
    # ── misfit ───────────────────────────────────────────────────────────────
    misfit_cmap: str | None = None,
    misfit_vlim: float | None = None,
    misfit_percentile: float = 98.0,
    # ── convergence ──────────────────────────────────────────────────────────
    train_loss: np.ndarray | None = None,
    val_loss: np.ndarray | None = None,
    convergence_epochs: np.ndarray | None = None,
    convergence_log_scale: bool = True,
    convergence_target: float | None = None,
    convergence_target_color: str = "#aaaaaa",
    convergence_target_ls: str = ":",
    convergence_target_lw: float = 0.8,
    convergence_target_label: str = "target",
    # ── per-station RMSE scatter map ─────────────────────────────────────────
    station_rmse: np.ndarray | None = None,
    rmse_cmap: str = "YlOrRd",
    rmse_vmax: float | None = None,
    rmse_marker: str = "o",
    rmse_marker_size: float = 60.0,
    rmse_threshold: float | None = None,
    rmse_colorbar_label: str = r"RMSE  ($\Delta\log_{10}\rho$)",
    # ── station overlays on sections and map ─────────────────────────────────
    show_stations: bool = True,
    station_marker: str = "^",
    station_color: str = "#333333",
    station_ms: float = 4.0,
    # ── GCN graph edges on map ────────────────────────────────────────────────
    adjacency: np.ndarray | None = None,
    show_graph_edges: bool = True,
    edge_color: str = "#bbbbbb",
    edge_lw: float = 0.4,
    edge_alpha: float = 0.5,
    # ── fault planes ─────────────────────────────────────────────────────────
    fault_x_positions: list[float] | None = None,
    fault_y_positions: list[float] | None = None,
    fault_color: str = "#eeeeee",
    fault_lw: float = 0.9,
    fault_ls: str = "--",
    fault_alpha: float = 0.75,
    fault_label_fontsize: float = 6.5,
    # ── annotations (per panel type) ─────────────────────────────────────────
    annotations_map: list[dict] | None = None,
    annotations_inline: list[dict] | None = None,
    annotations_crossline: list[dict] | None = None,
    # ── colourbar labels ─────────────────────────────────────────────────────
    colorbar_label: str = r"$\log_{10}(\rho)$  (Ω·m)",
    misfit_colorbar_label: str = r"$\Delta\log_{10}(\rho)$",
    # ── panel titles ─────────────────────────────────────────────────────────
    title_map_true: str = "(a) True — depth slice",
    title_map_pred: str | None = None,
    title_map_misfit: str = r"(c) Misfit — depth slice",
    title_inline_true: str = "(d) True — inline section",
    title_inline_pred: str | None = None,
    title_crossline_pred: str | None = None,
    title_rmse_map: str = "Per-station RMSE",
    title_convergence: str = "Training convergence",
    # ── axis labels ──────────────────────────────────────────────────────────
    xlabel_ew: str = "E–W distance (km)",
    ylabel_ns: str = "N–S distance (km)",
    xlabel_inline: str = "E–W distance (km)",
    xlabel_crossline: str = "N–S distance (km)",
    ylabel_section: str = "Depth (m)",
    suptitle: str = "",
    # ── layout ───────────────────────────────────────────────────────────────
    figsize: tuple[float, float] | None = None,
    hspace: float = 0.40,
    wspace: float = 0.32,
    map_height_ratio: float = 1.0,
    section_height_ratio: float = 1.4,
    convergence_height_ratio: float = 0.55,
    axes: dict[str, Axes] | None = None,
) -> Figure:
    """
    Multi-panel 3-D inversion result figure for
    :class:`~pycsamt.ai.inversion.GCNInverter3D` output.

    Up to seven panels are arranged in two or three rows:

    **Comparison mode** (when *log_true* is supplied):

    * **(a) True depth slice** — plan view at *depth_slice_idx*.
    * **(b) Predicted depth slice** — same scale.
    * **(c) Depth-slice misfit** ``Δlog₁₀ρ`` — diverging colourmap.
    * **(d) True inline section** — vertical x-z section at *inline_idx*.
    * **(e) Predicted inline section** — same scale.
    * **(f) Per-station RMSE scatter map** — 2-D spatial error distribution.
    * **(g) Training convergence** — optional full-width bottom row.

    **Prediction mode** (no *log_true*):

    * **(a) Depth slice** — plan view.
    * **(b) Inline section** — x-z plane.
    * **(c) Crossline section** — y-z plane.
    * **(d) Training convergence** — spans 2 columns.
    * **(e) Per-station RMSE map** — last column.

    Parameters
    ----------
    log_pred : ndarray, shape (n_depth, n_y, n_x)  **or**  (n_stations, n_depth)
        Predicted :math:`\\log_{10}(\\rho)` volume.  When the array is
        2-D, it is interpreted as per-station 1-D models and interpolated
        onto a regular grid using *station_xy*, *interp_nx*, *interp_ny*,
        and *interp_method* (requires :mod:`scipy`).
    log_true : ndarray (n_depth, n_y, n_x) or (n_stations, n_depth) or None
        Ground-truth volume.  Enables comparison mode.
    depths : ndarray (n_depth,) or None
        Depth node positions in metres.
    x_coords : ndarray (n_x,) or None
        E–W coordinates.  Defaults to equally spaced with *x_spacing*.
    y_coords : ndarray (n_y,) or None
        N–S coordinates.  Defaults to equally spaced with *y_spacing*.
    depth_max : float
        Maximum depth used when *depths* is ``None``.
    x_spacing, y_spacing : float
        Node spacing when *x_coords* / *y_coords* are ``None``.
    station_xy : ndarray (n_stations, 2) or None
        Station (x, y) positions — required for RMSE scatter map and
        graph-edge overlay, and for 2-D ``log_pred`` interpolation.
    interp_method : ``"linear"``, ``"nearest"``, or ``"cubic"``
        Interpolation method for per-station → volume conversion.
    interp_nx, interp_ny : int
        Target grid size when auto-building from per-station data.
    depth_slice_idx : int or None
        Index along the depth axis for the horizontal map panels.
        Auto-set to ``n_depth // 4`` when ``None``.
    depth_slice_m : float or None
        Alternative: specify depth in metres (converted to index internally).
    inline_idx : int or None
        y-index for the inline section (x-z plane).  Defaults to midpoint.
    crossline_idx : int or None
        x-index for the crossline section (y-z plane).  Defaults to midpoint.
    show_depth_slice, show_inline, show_crossline : bool
        Toggle individual section panels.
    show_misfit : bool
        Show misfit panels (requires *log_true*).
    show_convergence : bool
        Show convergence panel (requires *train_loss*).
    show_rmse_map : bool
        Show the per-station RMSE scatter map (requires *station_xy* and
        *station_rmse* or pre-computable misfit).
    show_section_lines : bool
        Draw inline / crossline position indicators on depth-slice panels.
    vmin, vmax : float or None
        Shared colour-scale limits for resistivity panels.
    cmap : str or None
        Colourmap for resistivity panels.  Defaults to ``'RdYlBu_r'``.
    vmin_percentile, vmax_percentile : float
        Auto colour-scale percentile clipping.
    map_aspect : str
        Aspect ratio for map panels (``"equal"`` or ``"auto"``).
    misfit_cmap : str or None
        Colourmap for misfit panels.  Defaults to ``'RdBu_r'``.
    misfit_vlim : float or None
        Symmetric misfit colour-scale limit.
    misfit_percentile : float
        Percentile for auto misfit limit.
    train_loss, val_loss : ndarray (n_epochs,) or None
        Training and validation loss curves.
    convergence_epochs : ndarray (n_epochs,) or None
        Epoch indices (defaults to ``1 … n``).
    convergence_log_scale : bool
        Log y-axis on convergence panel.
    convergence_target : float or None
        Optional horizontal reference line on the convergence panel.
    station_rmse : ndarray (n_stations,) or None
        Pre-computed per-station RMSE values for the scatter map.
    rmse_cmap : str
        Colourmap for RMSE scatter.  Defaults to ``'YlOrRd'``.
    rmse_vmax : float or None
        Upper colour limit for RMSE.  Auto-computed when ``None``.
    rmse_marker, rmse_marker_size : str / float
        Scatter marker style and size.
    rmse_threshold : float or None
        Stations with RMSE above this value get a red outline ring.
    show_stations : bool
        Overlay station markers on section and map panels.
    station_marker, station_color, station_ms : str / float
        Station overlay marker style.
    adjacency : ndarray (n_stations, n_stations) or None
        Normalised adjacency matrix for GCN graph-edge drawing on map.
    show_graph_edges : bool
        Draw graph edges on the map panel.
    edge_color, edge_lw, edge_alpha : str / float
        Graph-edge visual style.
    fault_x_positions : list of float or None
        x-coordinates for vertical fault indicators drawn on map and inline
        section panels.
    fault_y_positions : list of float or None
        y-coordinates for fault indicators on map and crossline panels.
    fault_color, fault_lw, fault_ls, fault_alpha : str / float
        Fault-line style.
    annotations_map : list of dict or None
        Annotation dicts applied to the leftmost depth-slice map panel.
    annotations_inline : list of dict or None
        Annotation dicts applied to the inline section panel.
    annotations_crossline : list of dict or None
        Annotation dicts applied to the crossline section panel.
    colorbar_label : str
        Shared resistivity colourbar label.
    misfit_colorbar_label : str
        Misfit colourbar label.
    title_map_true, title_map_pred, title_map_misfit : str
        Map-row panel titles.  *title_map_pred* defaults to ``"(a)"`` or
        ``"(b)"`` automatically.
    title_inline_true, title_inline_pred : str
        Inline section panel titles.  Auto-set when ``None``.
    title_crossline_pred : str
        Crossline section title.  Auto-set when ``None``.
    title_rmse_map, title_convergence : str
        Diagnostic panel titles.
    xlabel_ew, ylabel_ns : str
        Axis labels for map panels.
    xlabel_inline, xlabel_crossline, ylabel_section : str
        Axis labels for vertical section panels.
    suptitle : str
        Figure-level super-title.
    figsize : (width, height) or None
        Auto-derived from layout mode when ``None``.
    hspace, wspace : float
        GridSpec spacing.
    map_height_ratio : float
        Relative height of map rows vs. section rows.
    section_height_ratio : float
        Relative height of section rows.
    convergence_height_ratio : float
        Relative height of the optional convergence row.
    axes : dict {str → Axes} or None
        Pre-built axes keyed by ``"map_true"``, ``"map_pred"``,
        ``"map_misfit"``, ``"inline_true"``, ``"inline_pred"``,
        ``"crossline_pred"``, ``"rmse_map"``, ``"convergence"``.
        Skips figure/GridSpec creation when provided.

    Returns
    -------
    fig : Figure

    Examples
    --------
    Comparison mode with per-station GCN output:

    >>> fig = plot_inversion_result_3d(
    ...     gcn_pred,  # (n_stations, n_layers)
    ...     log_true=true_volume,  # (n_depth, n_y, n_x)
    ...     station_xy=coords,
    ...     depths=depths,
    ...     train_loss=train_loss,
    ...     val_loss=val_loss,
    ...     station_rmse=rmse_per_stn,
    ...     adjacency=A,
    ...     depth_slice_m=350.0,
    ...     fault_x_positions=[7.5],
    ...     suptitle="GCNInverter3D — survey block A",
    ... )

    Prediction-only mode (3 sections + convergence + RMSE map):

    >>> fig = plot_inversion_result_3d(
    ...     log_pred_volume,
    ...     depths=depths,
    ...     x_coords=xv,
    ...     y_coords=yv,
    ...     station_xy=sxy,
    ...     station_rmse=rmse,
    ...     train_loss=tl,
    ...     val_loss=vl,
    ... )
    """
    # ── 1. input coercion ─────────────────────────────────────────────────────
    log_pred = np.asarray(log_pred, dtype=float)

    # Per-station 2-D input: (n_stations, n_depth) → volume
    if log_pred.ndim == 2:
        if station_xy is None:
            raise ValueError(
                "station_xy= is required to interpolate a 2-D per-station "
                "log_pred (n_stations, n_depth) onto a 3-D volume."
            )
        sxy = np.asarray(station_xy, dtype=float)
        n_st, n_d = log_pred.shape
        xi = np.linspace(sxy[:, 0].min(), sxy[:, 0].max(), interp_nx)
        yi = np.linspace(sxy[:, 1].min(), sxy[:, 1].max(), interp_ny)
        if x_coords is None:
            x_coords = xi
        if y_coords is None:
            y_coords = yi
        log_pred = _per_station_to_volume(
            log_pred, sxy, x_coords, y_coords, interp_method
        )
        if log_true is not None:
            log_true = np.asarray(log_true, dtype=float)
            if log_true.ndim == 2:
                log_true = _per_station_to_volume(
                    log_true, sxy, x_coords, y_coords, interp_method
                )
    elif log_pred.ndim != 3:
        raise ValueError(
            f"log_pred must be 3-D (n_depth, n_y, n_x) or 2-D "
            f"(n_stations, n_depth); got shape {log_pred.shape}."
        )
    else:
        if log_true is not None:
            log_true = np.asarray(log_true, dtype=float)

    n_depth, n_y, n_x = log_pred.shape

    # ── 2. grid geometry ──────────────────────────────────────────────────────
    if depths is None:
        depths = np.linspace(0.0, depth_max, n_depth)
    if x_coords is None:
        x_coords = np.arange(n_x) * float(x_spacing)
    if y_coords is None:
        y_coords = np.arange(n_y) * float(y_spacing)

    d_edges = _cell_edges(depths)
    x_edges = _cell_edges(x_coords)
    y_edges = _cell_edges(y_coords)

    # ── 3. slice indices ─────────────────────────────────────────────────────
    if depth_slice_idx is None:
        if depth_slice_m is not None:
            depth_slice_idx = int(np.argmin(np.abs(depths - depth_slice_m)))
        else:
            depth_slice_idx = max(0, n_depth // 4)
    depth_slice_idx = int(np.clip(depth_slice_idx, 0, n_depth - 1))

    if inline_idx is None:
        inline_idx = n_y // 2
    inline_idx = int(np.clip(inline_idx, 0, n_y - 1))

    if crossline_idx is None:
        crossline_idx = n_x // 2
    crossline_idx = int(np.clip(crossline_idx, 0, n_x - 1))

    # ── 4. colour limits ──────────────────────────────────────────────────────
    ref = log_true if log_true is not None else log_pred
    if vmin is None:
        vmin = float(np.nanpercentile(ref, vmin_percentile))
    if vmax is None:
        vmax = float(np.nanpercentile(ref, vmax_percentile))
    if cmap is None:
        cmap = EM_CMAPS["resistivity"]

    # ── 5. misfit ─────────────────────────────────────────────────────────────
    has_true = log_true is not None
    _show_mis = show_misfit and has_true
    misfit_3d: np.ndarray | None = None
    if _show_mis:
        misfit_3d = log_pred - log_true
        if misfit_cmap is None:
            misfit_cmap = EM_CMAPS["misfit"]
        if misfit_vlim is None:
            misfit_vlim = float(np.nanpercentile(np.abs(misfit_3d), misfit_percentile))
    else:
        misfit_cmap = misfit_cmap or EM_CMAPS["misfit"]

    # ── 6. RMSE scatter ───────────────────────────────────────────────────────
    if station_xy is not None:
        station_xy = np.asarray(station_xy, dtype=float)
    _show_rmse = show_rmse_map and station_xy is not None and station_rmse is not None
    if station_rmse is not None:
        station_rmse = np.asarray(station_rmse, dtype=float)
        if rmse_vmax is None:
            rmse_vmax = float(np.nanpercentile(station_rmse, 97))

    # ── 7. convergence ────────────────────────────────────────────────────────
    _show_conv = show_convergence and (train_loss is not None)
    if train_loss is not None:
        train_loss = np.asarray(train_loss, dtype=float)
    if val_loss is not None:
        val_loss = np.asarray(val_loss, dtype=float)
    if convergence_epochs is None and train_loss is not None:
        convergence_epochs = np.arange(1, len(train_loss) + 1)

    # ── 8. auto panel titles ──────────────────────────────────────────────────
    # comparison mode: true panels on col 0, pred on col 1, misfit on col 2
    if comparison_mode := (has_true and _show_mis and show_depth_slice):
        if title_map_pred is None:
            title_map_pred = "(b) Predicted — depth slice"
        if title_inline_pred is None:
            title_inline_pred = "(e) Predicted — inline section"
        if title_crossline_pred is None:
            title_crossline_pred = "(e) Crossline section"
    else:
        letter = "a"
        if title_map_pred is None:
            title_map_pred = f"({letter}) Depth slice — predicted"
        if title_inline_pred is None:
            title_inline_pred = (
                "(b) Inline section — predicted"
                if show_depth_slice
                else "(a) Inline section"
            )
        if title_crossline_pred is None:
            title_crossline_pred = "(c) Crossline section — predicted"

    # ── 9. figure / axes creation ─────────────────────────────────────────────
    using_external = axes is not None
    ax: dict[str, Axes] = {}

    if using_external:
        ax = dict(axes)
        fig = next(iter(ax.values())).get_figure()
    else:
        if comparison_mode:
            # Row 0: map_true | map_pred | map_misfit
            # Row 1: inline_true | inline_pred | rmse_map
            # Row 2 (opt): convergence
            n_rows = 2 + (1 if _show_conv else 0)
            hr = [map_height_ratio, section_height_ratio] + (
                [convergence_height_ratio] if _show_conv else []
            )
            if figsize is None:
                h = 4.5 * map_height_ratio + 5.0 * section_height_ratio
                h += 2.5 * convergence_height_ratio if _show_conv else 0.0
                figsize = (14.0, h)
        else:
            # Row 0: up to 3 top sections
            top_keys: list[str] = []
            if show_depth_slice:
                top_keys.append("map_pred")
            if show_inline:
                top_keys.append("inline_pred")
            if show_crossline:
                top_keys.append("crossline_pred")
            if not top_keys:
                top_keys = ["map_pred"]
            n_top = len(top_keys)

            has_bot = _show_conv or _show_rmse
            n_rows = 1 + (1 if has_bot else 0)
            hr = [section_height_ratio] + (
                [convergence_height_ratio] if has_bot else []
            )
            if figsize is None:
                w = _TOP_WIDTHS.get(n_top, n_top * 5.0)
                h = 5.0 + (2.5 if has_bot else 0.0)
                figsize = (w, h)

        fig = plt.figure(figsize=figsize)

        if comparison_mode:
            gs_outer = fig.add_gridspec(n_rows, 1, height_ratios=hr, hspace=hspace)
            gs_top = gs_outer[0].subgridspec(1, 3, wspace=wspace)
            gs_mid = gs_outer[1].subgridspec(1, 3, wspace=wspace)
            ax["map_true"] = fig.add_subplot(gs_top[0, 0])
            ax["map_pred"] = fig.add_subplot(gs_top[0, 1])
            ax["map_misfit"] = fig.add_subplot(gs_top[0, 2])
            ax["inline_true"] = fig.add_subplot(gs_mid[0, 0])
            ax["inline_pred"] = fig.add_subplot(gs_mid[0, 1])
            ax["rmse_map"] = fig.add_subplot(gs_mid[0, 2])
            if _show_conv:
                ax["convergence"] = fig.add_subplot(gs_outer[2])
        else:
            gs_outer = fig.add_gridspec(n_rows, 1, height_ratios=hr, hspace=hspace)
            gs_top = gs_outer[0].subgridspec(1, n_top, wspace=wspace)
            for i, key in enumerate(top_keys):
                ax[key] = fig.add_subplot(gs_top[0, i])
            if has_bot:
                gs_bot = gs_outer[1].subgridspec(1, n_top, wspace=wspace)
                if _show_conv and _show_rmse:
                    ax["convergence"] = fig.add_subplot(gs_bot[0, : max(1, n_top - 1)])
                    ax["rmse_map"] = fig.add_subplot(gs_bot[0, n_top - 1])
                elif _show_conv:
                    ax["convergence"] = fig.add_subplot(gs_bot[0, :])
                elif _show_rmse:
                    ax["rmse_map"] = fig.add_subplot(gs_bot[0, :])

    # ── 10. render section panels ────────────────────────────────────────────
    depth_label_str = f"{depths[depth_slice_idx]:.0f} m"

    def _station_surface(sec_ax, along_coords):
        """Draw site markers at the surface of a section panel."""
        if show_stations and station_xy is not None:
            sec_ax.plot(
                along_coords,
                np.zeros(len(along_coords)),
                station_marker,
                ms=station_ms,
                color=station_color,
                zorder=5,
                clip_on=False,
            )

    im_res: Any | None = None  # last resistivity mappable (for shared cbar)

    # ── depth-slice maps ──────────────────────────────────────────────────────
    for key, data, ttl in [
        (
            "map_true",
            log_true,
            title_map_true if log_true is not None else None,
        ),
        ("map_pred", log_pred, title_map_pred),
        ("map_misfit", misfit_3d, title_map_misfit),
    ]:
        if key not in ax or data is None:
            continue
        slc = data[depth_slice_idx, :, :]
        _v0 = -misfit_vlim if key == "map_misfit" else vmin
        _v1 = misfit_vlim if key == "map_misfit" else vmax
        _cm = misfit_cmap if key == "map_misfit" else cmap
        _ttl = (ttl or "").replace("depth slice", f"depth slice  ({depth_label_str})")
        im = _draw_map_panel(
            ax[key],
            slc,
            x_edges,
            y_edges,
            _v0,
            _v1,
            _cm,
            _ttl,
            xlabel_ew,
            ylabel_ns,
            aspect=map_aspect,
        )
        if key != "map_misfit":
            im_res = im

        if show_stations and station_xy is not None:
            ax[key].scatter(
                station_xy[:, 0],
                station_xy[:, 1],
                s=station_ms**2 * 0.8,
                marker=station_marker,
                c=station_color,
                zorder=5,
                edgecolors="none",
            )

        if (
            show_graph_edges
            and adjacency is not None
            and station_xy is not None
            and key != "map_misfit"
        ):
            _draw_graph_edges(
                ax[key],
                station_xy,
                adjacency,
                edge_color,
                edge_lw,
                edge_alpha,
            )

        if show_section_lines and key != "map_misfit":
            _overlay_cut_lines(
                ax[key],
                x_coords,
                y_coords,
                inline_idx,
                crossline_idx,
            )

        if fault_x_positions and key != "map_misfit":
            for xf in fault_x_positions:
                ax[key].axvline(
                    xf,
                    color=fault_color,
                    lw=fault_lw,
                    ls=fault_ls,
                    alpha=fault_alpha,
                    zorder=3,
                )
        if fault_y_positions and key != "map_misfit":
            for yf in fault_y_positions:
                ax[key].axhline(
                    yf,
                    color=fault_color,
                    lw=fault_lw,
                    ls=fault_ls,
                    alpha=fault_alpha,
                    zorder=3,
                )

    if annotations_map:
        leftmost_map = next((k for k in ("map_true", "map_pred") if k in ax), None)
        if leftmost_map:
            _add_annotations(ax[leftmost_map], annotations_map)

    # ── inline sections (x-z plane at y = inline_idx) ────────────────────────
    for key, data, ttl in [
        ("inline_true", log_true, title_inline_true),
        ("inline_pred", log_pred, title_inline_pred),
    ]:
        if key not in ax or data is None:
            continue
        slc = data[:, inline_idx, :]  # (n_depth, n_x)
        im = _draw_section(
            ax[key],
            slc,
            x_edges,
            d_edges,
            vmin,
            vmax,
            cmap,
            ttl + f"  (y = {y_coords[inline_idx]:.1f} km)",
            xlabel_inline,
            ylabel_section,
            show_sites=False,
            site_marker=station_marker,
            site_color=station_color,
            site_ms=station_ms,
            stations=x_coords,
            station_labels=None,
            tick_cfg=StationTickConfig(),
        )
        im_res = im

        if show_stations and station_xy is not None:
            # project stations onto this section line
            ax[key].plot(
                station_xy[:, 0],
                np.zeros(len(station_xy)),
                station_marker,
                ms=station_ms,
                color=station_color,
                zorder=5,
                clip_on=False,
            )

        if fault_x_positions:
            _add_faults(
                ax[key],
                fault_x_positions,
                fault_color,
                fault_lw,
                fault_ls,
                fault_alpha,
                fault_label_fontsize,
            )

        if annotations_inline:
            _add_annotations(ax[key], annotations_inline)

    # ── crossline section (y-z plane at x = crossline_idx) ───────────────────
    for key, data, ttl in [
        ("crossline_pred", log_pred, title_crossline_pred),
    ]:
        if key not in ax or data is None:
            continue
        slc = data[:, :, crossline_idx]  # (n_depth, n_y)
        im = _draw_section(
            ax[key],
            slc,
            y_edges,
            d_edges,
            vmin,
            vmax,
            cmap,
            ttl + f"  (x = {x_coords[crossline_idx]:.1f} km)",
            xlabel_crossline,
            ylabel_section,
            show_sites=False,
            site_marker=station_marker,
            site_color=station_color,
            site_ms=station_ms,
            stations=y_coords,
            station_labels=None,
            tick_cfg=StationTickConfig(),
        )
        im_res = im

        if fault_y_positions:
            _add_faults(
                ax[key],
                fault_y_positions,
                fault_color,
                fault_lw,
                fault_ls,
                fault_alpha,
                fault_label_fontsize,
            )

        if annotations_crossline:
            _add_annotations(ax[key], annotations_crossline)

    # ── RMSE scatter map ──────────────────────────────────────────────────────
    im_rmse = None
    if "rmse_map" in ax and _show_rmse and station_rmse is not None:
        im_rmse = _draw_rmse_scatter(
            ax["rmse_map"],
            station_xy,
            station_rmse,
            rmse_cmap,
            rmse_vmax,
            rmse_marker,
            rmse_marker_size,
            rmse_threshold,
            title_rmse_map,
            xlabel_ew,
            ylabel_ns,
            aspect=map_aspect,
        )
        if show_section_lines:
            _overlay_cut_lines(
                ax["rmse_map"],
                x_coords,
                y_coords,
                inline_idx,
                crossline_idx,
                color="#444444",
            )

    # ── convergence panel ─────────────────────────────────────────────────────
    if "convergence" in ax and _show_conv:
        _render_convergence(
            ax["convergence"],
            convergence_epochs,
            train_loss,
            val_loss,
            convergence_log_scale,
            convergence_target,
            convergence_target_color,
            convergence_target_ls,
            convergence_target_lw,
            convergence_target_label,
            title_convergence,
        )

    # ── 11. shared colorbars ──────────────────────────────────────────────────
    if not using_external:
        # resistivity: shared across all non-misfit panels
        if im_res is not None:
            res_ax_list = [
                ax[k]
                for k in (
                    "map_true",
                    "map_pred",
                    "inline_true",
                    "inline_pred",
                    "crossline_pred",
                )
                if k in ax
            ]
            if res_ax_list:
                cb_r = fig.colorbar(
                    im_res,
                    ax=res_ax_list,
                    shrink=0.80,
                    pad=0.02,
                    fraction=0.015,
                )
                cb_r.set_label(colorbar_label, fontsize=8)
                cb_r.ax.tick_params(labelsize=7)

        # misfit colorbar
        if "map_misfit" in ax and misfit_3d is not None:
            im_mf = (
                ax["map_misfit"].collections[0]
                if ax["map_misfit"].collections
                else None
            )
            if im_mf is not None:
                cb_m = fig.colorbar(
                    im_mf,
                    ax=ax["map_misfit"],
                    shrink=0.88,
                    pad=0.02,
                    fraction=0.044,
                )
                cb_m.set_label(misfit_colorbar_label, fontsize=8)
                cb_m.ax.tick_params(labelsize=7)

        # RMSE colourbar
        if im_rmse is not None and "rmse_map" in ax:
            cb_rr = fig.colorbar(
                im_rmse,
                ax=ax["rmse_map"],
                shrink=0.88,
                pad=0.02,
                fraction=0.044,
            )
            cb_rr.set_label(rmse_colorbar_label, fontsize=8)
            cb_rr.ax.tick_params(labelsize=7)

    # ── 12. suptitle ─────────────────────────────────────────────────────────
    if suptitle:
        fig.suptitle(suptitle, fontsize=10, fontweight="bold", y=1.002)

    return fig
