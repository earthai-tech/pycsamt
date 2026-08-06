# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Plotting helpers for :mod:`pycsamt.inversion` results.

This module provides compact quick-look plots for
:class:`pycsamt.inversion.results.InversionResult` objects. The public
functions are intentionally small and backend-neutral: each plot consumes the
common result API rather than SimPEG, pyGIMLi, Occam2D, ModEM, or built-in
solver internals.

The plotting helpers follow the shared pyCSAMT plotting API. Section-like views
use :data:`pycsamt.api.section.PYCSAMT_SECTION`, diagnostic line styling uses
:data:`pycsamt.api.style.PYCSAMT_STYLE`, and optional saving goes through
:func:`pycsamt.api.plot.save_fig`.

Available plots
---------------
``plot_model``
    Plot a recovered 1-D layered model as a depth profile, or a 2-D inversion
    section as a profile/depth color mesh.
``plot_rms``
    Plot station RMS values when available, otherwise plot the global weighted
    RMS value.

Coordinate and value conventions
--------------------------------
Depth is positive downward. Model values are stored as
``log10(rho / ohm m)`` by the interpretation model container, and plotting uses
that log scale by default. Pass ``log_rho=False`` to display linear
resistivity in ohm metres.

Examples
--------
Plot a recovered inversion section::

    >>> from pycsamt.inversion.plot import plot_model
    >>> ax = plot_model(result, section="compact")  # doctest: +SKIP
    >>> ax.get_ylabel()  # doctest: +SKIP
    'Depth (m)'

Plot RMS diagnostics::

    >>> from pycsamt.inversion.plot import plot_rms
    >>> ax = plot_rms(result)  # doctest: +SKIP
    >>> ax.get_ylabel()  # doctest: +SKIP
    'Weighted RMS'

Save figures using global pyCSAMT plot settings::

    >>> from pycsamt.inversion.plot import plot_model
    >>> plot_model(result, savepath="figures/inversion_model")  # doctest: +SKIP

See Also
--------
pycsamt.inversion.results.InversionResult
    Backend-neutral result object consumed by the plotting helpers.
pycsamt.inversion.export
    File export helpers for CSV, NPZ, GeoJSON, VTK, GeoTIFF, and ZIP products.
pycsamt.api.plot.save_fig
    Shared figure-saving helper used by this module.

References
----------
.. [module-1] Hunter, J. D. (2007). Matplotlib: A 2D graphics environment.
   *Computing in Science & Engineering*, 9(3), 90-95.
.. [module-2] Tufte, E. R. (2001). *The Visual Display of Quantitative Information*,
   2nd edition. Graphics Press.
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

    ``plot_model`` is the quick-look model visualizer for
    :class:`pycsamt.inversion.results.InversionResult`. It first converts the
    result through ``result.to_resistivity_model()`` and then chooses the plot
    type from the recovered grid shape:

    * one model column -> a depth profile using ``Axes.step``;
    * multiple columns -> a profile/depth section using ``Axes.pcolormesh``.

    The depth axis is positive downward, matching
    :class:`pycsamt.interp.ResistivityModel` and the rest of the interpretation
    API. Section sizing, station labels, colorbar style, and optional saving use
    the shared pyCSAMT plotting configuration.

    Parameters
    ----------
    result : InversionResult
        Result produced by :mod:`pycsamt.inversion`. The result must be
        convertible to a 2-D :class:`pycsamt.interp.ResistivityModel` through
        ``result.to_resistivity_model()``.
    ax : matplotlib Axes, optional
        Existing axes to draw into. If omitted, a new figure and axes are
        created using the selected section style.
    log_rho : bool, default True
        Plot ``log10(rho / ohm m)`` values. If ``False``, values are converted
        to linear resistivity in ohm metres before plotting.
    cmap : str, default "jet_r"
        Matplotlib colormap for 2-D sections. Ignored for single-column 1-D
        depth profiles.
    colorbar : bool, default True
        Add a colorbar for 2-D section plots. Ignored for single-column 1-D
        depth profiles.
    show_stations : bool, default True
        Draw station markers and labels using the section station preset when
        station positions are available.
    section : str or SectionStyle, default "inversion"
        Shared section style preset name or explicit
        :class:`pycsamt.api.section.SectionStyle` object. Names are resolved
        through :data:`pycsamt.api.section.PYCSAMT_SECTION`.
    title : str, optional
        Axes title. If omitted, a backend/method/dimension summary is used.
    savepath : str, optional
        If given, save the figure using :func:`pycsamt.api.plot.save_fig`.
        The path may omit the extension when global pyCSAMT plot formats are
        configured.
    savefig_kw : dict, optional
        Extra keyword arguments forwarded to ``save_fig``.

    Returns
    -------
    matplotlib.axes.Axes
        Axes containing the model plot.

    Notes
    -----
    The function does not call ``matplotlib.pyplot.show``. This keeps it safe
    for scripts, notebooks, test suites, and batch figure generation. Use the
    returned axes to further customize labels, limits, annotations, or overlays.

    Examples
    --------
    Plot an inversion result returned by a workflow::

        >>> from pycsamt.inversion.plot import plot_model
        >>> ax = plot_model(result, section="compact", colorbar=False)  # doctest: +SKIP
        >>> ax.get_ylabel()  # doctest: +SKIP
        'Depth (m)'

    Save a publication copy using the shared pyCSAMT plot settings::

        >>> from pycsamt.inversion.plot import plot_model
        >>> plot_model(result, savepath="figures/inversion_model")  # doctest: +SKIP

    Draw linear resistivity instead of log10 resistivity::

        >>> from pycsamt.inversion.plot import plot_model
        >>> plot_model(result, log_rho=False, cmap="viridis")  # doctest: +SKIP

    References
    ----------
    .. [plot-model-1] Tufte, E. R. (2001). *The Visual Display of Quantitative
       Information*, 2nd edition. Graphics Press.
    .. [plot-model-2] Hunter, J. D. (2007). Matplotlib: A 2D graphics environment.
       *Computing in Science & Engineering*, 9(3), 90-95.
    .. [plot-model-3] Chave, A. D. and Jones, A. G. (2012). *The Magnetotelluric Method:
       Theory and Practice*. Cambridge University Press.
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
    values = rho if log_rho else 10.0**rho
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
        xlabel = (
            r"$\log_{10}\rho$ ($\Omega\cdot$m)"
            if log_rho
            else r"$\rho$ ($\Omega\cdot$m)"
        )
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
                station_labels = labels or [
                    f"S{i:03d}" for i in range(station_x.size)
                ]
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
    """Plot inversion RMS misfit diagnostics.

    ``plot_rms`` visualizes the weighted root-mean-square misfit stored on an
    :class:`pycsamt.inversion.results.InversionResult`. When
    ``result.metadata["station_rms"]`` is available, the function draws one
    marker per station. Otherwise it draws a single global RMS bar from
    ``result.rms``.

    Parameters
    ----------
    result : InversionResult
        Inversion result containing a global ``rms`` value and optionally
        ``metadata["station_rms"]``.
    ax : matplotlib Axes, optional
        Existing axes to draw into. If omitted, a new compact figure and axes
        are created.
    title : str, default "Inversion misfit"
        Axes title.
    savepath : str, optional
        If given, save the figure using :func:`pycsamt.api.plot.save_fig`.
    savefig_kw : dict, optional
        Extra keyword arguments forwarded to ``save_fig``.
    **kwargs
        Additional Matplotlib keyword arguments forwarded to ``Axes.plot`` for
        station RMS curves or ``Axes.bar`` for a global RMS bar.

    Returns
    -------
    matplotlib.axes.Axes
        Axes containing the RMS plot.

    Notes
    -----
    RMS near 1 is often interpreted as consistency with the supplied data-error
    model, but the correct target depends on data quality, error floors,
    regularization strength, and backend conventions. This plot is therefore a
    diagnostic companion to model plots rather than a standalone quality
    guarantee.

    Examples
    --------
    Plot the global RMS or station RMS diagnostics::

        >>> from pycsamt.inversion.plot import plot_rms
        >>> ax = plot_rms(result)  # doctest: +SKIP
        >>> ax.get_ylabel()  # doctest: +SKIP
        'Weighted RMS'

    Customize line/bar appearance and save the plot::

        >>> from pycsamt.inversion.plot import plot_rms
        >>> plot_rms(result, color="black", savepath="figures/rms")  # doctest: +SKIP

    References
    ----------
    .. [plot-rms-1] Aster, R. C., Borchers, B. and Thurber, C. H. (2018). *Parameter
       Estimation and Inverse Problems*, 3rd edition. Elsevier.
    .. [plot-rms-2] Constable, S. C., Parker, R. L. and Constable, C. G. (1987).
       Occam's inversion: A practical algorithm for generating smooth models
       from electromagnetic sounding data. *Geophysics*, 52(3), 289-300.
    """
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
