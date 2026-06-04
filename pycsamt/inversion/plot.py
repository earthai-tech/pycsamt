# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Plotting helpers for :mod:`pycsamt.inversion` results.

The functions in this module follow the shared pyCSAMT plotting API:
section-like views use :data:`pycsamt.api.section.PYCSAMT_SECTION`,
station labels use :data:`pycsamt.api.station.PYCSAMT_STATION_RENDERING`,
multi-line diagnostics use :data:`pycsamt.api.style.PYCSAMT_STYLE`, and
optional saving goes through :func:`pycsamt.api.plot.save_fig`.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..api.plot import save_fig
from ..api.section import PYCSAMT_SECTION, SectionStyle
from ..api.style import PYCSAMT_STYLE
from .results import InversionResult

__all__ = ["plot_model", "plot_rms"]


def plot_model(
    result: InversionResult,
    ax: Any = None,
    *,
    log_rho: bool = True,
    cmap: str = "jet_r",
    colorbar: bool = True,
    show_stations: bool = True,
    section: str | SectionStyle = "inversion",
    title: str | None = None,
    savepath: str | None = None,
    savefig_kw: dict[str, Any] | None = None,
):
    """Plot a recovered 1-D or 2-D resistivity model.

    Parameters
    ----------
    result : InversionResult
        Result produced by :mod:`pycsamt.inversion`.
    ax : matplotlib Axes, optional
        Existing axes. If omitted, a new figure and axes are created.
    log_rho : bool, default True
        Plot log10 resistivity. If ``False``, values are converted to ohm-m.
    cmap : str, default "jet_r"
        Matplotlib colormap for 2-D sections.
    colorbar : bool, default True
        Add a colorbar for 2-D plots.
    show_stations : bool, default True
        Draw station markers and labels using the section station preset.
    section : str or SectionStyle, default "inversion"
        Shared section style preset or explicit style object.
    title : str, optional
        Axes title. Defaults to a backend/method summary.
    savepath : str, optional
        If given, save the figure using :func:`pycsamt.api.plot.save_fig`.
    savefig_kw : dict, optional
        Extra keyword arguments forwarded to ``save_fig``.
    """
    import matplotlib.pyplot as plt

    section_style = _resolve_section_style(section)
    if ax is None:
        # Figure size is refined after the model is loaded below.
        fig = None
    else:
        fig = ax.get_figure()
    model = result.to_resistivity_model()
    rho = np.asarray(model.rho_2d, dtype=float)
    values = rho if log_rho else 10.0 ** rho
    x = np.asarray(model.x_centers, dtype=float)
    z = np.asarray(model.z_centers, dtype=float)
    labels = list(model.station_names)

    if ax is None:
        figsize = section_style.figsize_for(
            n_stations=max(values.shape[1], len(labels)),
            n_y=values.shape[0],
            labels=labels,
            colorbar=colorbar and values.shape[1] > 1,
        )
        fig, ax = plt.subplots(
            figsize=figsize,
            constrained_layout=section_style.figure.constrained,
        )

    if values.shape[1] == 1:
        y = values[:, 0]
        ax.step(y, z, where="mid", **_model_line_kwargs())
        xlabel = r"$\log_{10}\rho$ ($\Omega\cdot$m)" if log_rho else r"$\rho$ ($\Omega\cdot$m)"
        section_style.apply_axis(
            ax,
            xlabel=xlabel,
            ylabel="Depth (m)",
            title=_title(result, title),
        )
    else:
        im = ax.pcolormesh(
            _edges(x),
            _edges(z),
            values,
            shading="auto",
            cmap=cmap,
        )
        section_style.apply_axis(
            ax,
            xlabel="Station",
            ylabel="Depth (m)",
            title=_title(result, title),
        )
        if colorbar:
            label = (
                r"$\log_{10}\rho$ ($\Omega\cdot$m)"
                if log_rho
                else r"$\rho$ ($\Omega\cdot$m)"
            )
            section_style.add_colorbar(im, ax, label=label)
        if show_stations:
            station_x = np.asarray(model.station_x, dtype=float)
            if station_x.size:
                station_labels = labels or [f"S{i:03d}" for i in range(station_x.size)]
                section_style.apply_stations(
                    ax,
                    station_x,
                    station_labels,
                    xlim=(float(_edges(x)[0]), float(_edges(x)[-1])),
                )

    _finish_figure(ax, section_style)
    if savepath:
        save_fig(ax, savepath, **(savefig_kw or {}))
    return ax


def plot_rms(
    result: InversionResult,
    ax: Any = None,
    *,
    title: str = "Inversion misfit",
    savepath: str | None = None,
    savefig_kw: dict[str, Any] | None = None,
    **kwargs: Any,
):
    """Plot station RMS values when available, otherwise one global value."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(6.5, 3.2))
    station_rms = result.metadata.get("station_rms")
    if station_rms is not None:
        y = np.asarray(station_rms, dtype=float)
        x = np.arange(y.size)
        line_kw = PYCSAMT_STYLE.multiline.line_kwargs(0, 1, marker="o")
        line_kw.update(kwargs)
        ax.plot(x, y, **line_kw)
        ax.set_xlabel("Station index")
    else:
        y = np.asarray([result.rms], dtype=float)
        bar_kw = {"color": PYCSAMT_STYLE.mt.xy.color, "alpha": 0.85}
        bar_kw.update(kwargs)
        ax.bar([0], y, **bar_kw)
        ax.set_xticks([0], ["global"])
    ax.set_ylabel("Weighted RMS")
    ax.set_title(title)
    ax.grid(False)
    if savepath:
        save_fig(ax, savepath, **(savefig_kw or {}))
    return ax


def _resolve_section_style(section: str | SectionStyle) -> SectionStyle:
    if isinstance(section, SectionStyle):
        return section.copy()
    return PYCSAMT_SECTION.style_for(str(section)).copy()


def _title(result: InversionResult, title: str | None) -> str | None:
    if title is not None:
        return title
    return f"{result.backend} {result.method.upper()} {result.dimension} model"


def _model_line_kwargs() -> dict[str, Any]:
    st = PYCSAMT_STYLE.mt.xy
    return {
        "color": st.color,
        "lw": st.lw,
        "alpha": st.alpha,
    }


def _finish_figure(ax: Any, section_style: SectionStyle) -> None:
    fig = ax.get_figure()
    if section_style.figure.tight:
        try:
            fig.tight_layout()
        except Exception:
            pass


def _edges(centers: np.ndarray) -> np.ndarray:
    centers = np.asarray(centers, dtype=float)
    if centers.size == 1:
        dx = max(abs(float(centers[0])) * 0.1, 1.0)
        return np.array([centers[0] - dx, centers[0] + dx], dtype=float)
    mids = 0.5 * (centers[:-1] + centers[1:])
    first = centers[0] - (mids[0] - centers[0])
    last = centers[-1] + (centers[-1] - mids[-1])
    return np.r_[first, mids, last]
