# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Plotting helpers for TEM data and transformed apparent resistivity."""

from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np

__all__ = ["PlotDecayCurve", "PlotTransformedRho"]

# Brand colours (consistent with pyCSAMT theme)
_ORANGE = "#fbb040"
_BLUE   = "#3e65b0"
_RED    = "#f15a29"


class PlotDecayCurve:
    r"""
    Log–log plot of the TEM transient decay curve.

    Plots :math:`|\partial B_z/\partial t|` vs time for one or
    more :class:`~pycsamt.tdem.TEMSounding` objects.  Error bars are
    drawn if ``sounding.error`` is available.

    Parameters
    ----------
    soundings : TEMSounding or list of TEMSounding
        Input data.
    title : str, optional
        Axes title.  Defaults to the station name.
    figsize : tuple, optional
        Figure size in inches.  Default ``(7, 5)``.

    Examples
    --------
    >>> from pycsamt.tdem import TEMSounding
    >>> from pycsamt.tdem.plot import PlotDecayCurve
    >>> import numpy as np
    >>> t = np.logspace(-5, -2, 30)
    >>> snd = TEMSounding(t, 1e-3 * t**(-5./2.), current=8.0, tx_area=1e4,
    ...                   station_name="S01")
    >>> PlotDecayCurve(snd).plot()
    """

    def __init__(
        self,
        soundings,
        *,
        title: Optional[str] = None,
        figsize: tuple = (7, 5),
    ) -> None:
        from ._base import TEMSounding as _TEM
        if isinstance(soundings, _TEM):
            soundings = [soundings]
        self.soundings = list(soundings)
        self.title = title
        self.figsize = figsize

    def plot(self, ax=None):
        """
        Draw the decay curves.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Existing axes to draw into.  Creates a new figure if not
            supplied.

        Returns
        -------
        matplotlib.axes.Axes
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _fig, ax = plt.subplots(figsize=self.figsize)

        colors = [_BLUE, _ORANGE, _RED, "green", "purple"]
        for i, snd in enumerate(self.soundings):
            c = colors[i % len(colors)]
            dBdt = np.abs(snd.dBdt())
            t = snd.time_gates
            label = snd.station_name or f"Sounding {i + 1}"

            if snd.error is not None:
                ax.errorbar(
                    t * 1e3, dBdt,
                    yerr=np.abs(snd.error),
                    fmt="o", color=c, markersize=4,
                    ecolor=c, elinewidth=0.8, capsize=2,
                    label=label,
                )
            else:
                ax.plot(t * 1e3, dBdt, "o-", color=c,
                        markersize=4, linewidth=1.2, label=label)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel(r"$|\partial B_z / \partial t|$ (T s$^{-1}$)")
        ax.set_title(
            self.title
            or (self.soundings[0].station_name if len(self.soundings) == 1 else "TEM decay")
        )
        ax.legend(fontsize=8)
        ax.grid(True, which="both", linestyle="--", alpha=0.4)

        plt.tight_layout()
        return ax


class PlotTransformedRho:
    r"""
    Plot TEM-derived apparent resistivity and phase vs pseudo-frequency.

    Accepts either raw :class:`~pycsamt.tdem.TEMSounding` objects (and
    runs :class:`~pycsamt.tdem.transform.LateTimeTransform` internally)
    or pre-computed result dicts from ``LateTimeTransform.transform()``.

    Parameters
    ----------
    data : TEMSounding, dict, or list
        Input data.  A ``TEMSounding`` is automatically transformed
        using the :class:`~pycsamt.tdem.transform.LateTimeTransform`
        with default settings.
    show_phase : bool, optional
        Draw a second panel with the phase curve.  Default ``True``.
    freq_convention : str, optional
        Passed to :class:`~pycsamt.tdem.transform.LateTimeTransform`
        when the input is a :class:`~pycsamt.tdem.TEMSounding`.
        Default ``"skin_depth"``.
    phase_mode : str, optional
        Passed to :class:`~pycsamt.tdem.transform.LateTimeTransform`.
        Default ``"weidelt"``.
    figsize : tuple, optional
        Figure size in inches.  Default ``(7, 6)`` (single panel) or
        ``(7, 9)`` (with phase).

    Examples
    --------
    >>> from pycsamt.tdem import TEMSounding
    >>> from pycsamt.tdem.plot import PlotTransformedRho
    >>> import numpy as np
    >>> t = np.logspace(-5, -2, 30)
    >>> M = 8.0 * 1e4
    >>> dBdt = (4e-7/np.pi) * ((4e-7)**2 * M / (20 * 10.**1.5 * t**2.5))**(2/3)
    >>> snd = TEMSounding(t, dBdt, current=8.0, tx_area=1e4, station_name="S01")
    >>> PlotTransformedRho(snd, phase_mode="weidelt").plot()
    """

    def __init__(
        self,
        data,
        *,
        show_phase: bool = True,
        freq_convention: str = "skin_depth",
        phase_mode: str = "weidelt",
        figsize: Optional[tuple] = None,
    ) -> None:
        self.show_phase = show_phase
        self.freq_convention = freq_convention
        self.phase_mode = phase_mode
        self.figsize = figsize

        from ._base import TEMSounding as _TEM
        if isinstance(data, _TEM):
            data = [data]
        if isinstance(data, dict):
            data = [data]

        results = []
        for item in data:
            if isinstance(item, _TEM):
                from .transform import LateTimeTransform
                tr = LateTimeTransform(
                    freq_convention=freq_convention,
                    phase_mode=phase_mode,
                )
                results.append(tr.transform(item))
            elif isinstance(item, dict) and "freq" in item:
                results.append(item)
            else:
                raise TypeError(
                    f"Expected TEMSounding or result dict, got {type(item)}"
                )
        self.results = results

    def plot(self, axes=None):
        """
        Draw the apparent-resistivity (and optionally phase) curves.

        Parameters
        ----------
        axes : matplotlib.axes.Axes or array of Axes, optional
            Pre-existing axes.  For ``show_phase=True`` supply a
            two-element array ``[ax_rho, ax_phase]``.

        Returns
        -------
        tuple of matplotlib.axes.Axes
        """
        import matplotlib.pyplot as plt

        n_panels = 2 if self.show_phase else 1
        if axes is None:
            if self.figsize is None:
                fsize = (7, 9) if self.show_phase else (7, 5)
            else:
                fsize = self.figsize
            if self.show_phase:
                fig, axes = plt.subplots(2, 1, figsize=fsize, sharex=True)
            else:
                fig, ax = plt.subplots(figsize=fsize)
                axes = [ax]
        else:
            if not hasattr(axes, "__len__"):
                axes = [axes]

        ax_rho = axes[0]
        ax_ph = axes[1] if self.show_phase else None

        colors = [_BLUE, _ORANGE, _RED, "green", "purple"]
        for i, res in enumerate(self.results):
            c = colors[i % len(colors)]
            label = res.get("station_name") or f"Site {i + 1}"
            freq = res["freq"]
            rho = res["rho_a"]
            phase = res.get("phase_xy")

            ax_rho.loglog(freq, rho, "o-", color=c,
                          markersize=4, linewidth=1.2, label=label)

            if ax_ph is not None and phase is not None:
                ax_ph.semilogx(freq, phase, "s--", color=c,
                               markersize=4, linewidth=1.0)

        ax_rho.set_ylabel(r"Apparent resistivity $\rho_a$ (Ω·m)")
        ax_rho.legend(fontsize=8)
        ax_rho.grid(True, which="both", linestyle="--", alpha=0.4)
        ax_rho.set_title("TEM-derived apparent resistivity")

        if ax_ph is not None:
            ax_ph.set_xlabel("Pseudo-frequency (Hz)")
            ax_ph.set_ylabel("Phase Z_xy (°)")
            ax_ph.set_ylim(0, 90)
            ax_ph.axhline(45, color="gray", linestyle=":", linewidth=0.8)
            ax_ph.grid(True, which="both", linestyle="--", alpha=0.4)
        else:
            ax_rho.set_xlabel("Pseudo-frequency (Hz)")

        plt.tight_layout()
        return tuple(axes[:n_panels])
