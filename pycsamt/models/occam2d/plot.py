# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Occam2D plotting — Python replacement for the MATLAB Occam2DMT scripts.

This module reimplements the functionality of the five original MATLAB
scripts in pure Python using matplotlib, and adds three diagnostic views:

    plotOccam2DMT.m             → PlotModel
    plotOccam2DMTResponse.m     → PlotResponse
    plotOccam2DMTPseudo.m       → PlotPseudo
    plotOccamIterMisfit.m       → PlotMisfit
    ExtractOccam2DMTProfile.m   → PlotModel.extract_profile()

Additional diagnostic views
---------------------------
    PlotSounding1D   — 1-D ρ–depth profiles extracted from the 2-D model
                       at each station column position.
    PlotSiteMisfit   — Per-site RMS bar chart + normalised-residual
                       pseudosection map.
    PlotResponseGrid — Compact grid of observed vs predicted curves for
                       all stations with per-site RMS annotations.

All plot classes accept an ``InversionResult`` object and expose a
``plot()`` method that returns a ``matplotlib.figure.Figure``.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np

from .base import OccamBase

__all__ = [
    "PlotModel",
    "PlotResponse",
    "PlotPseudo",
    "PlotMisfit",
    "PlotSounding1D",
    "PlotSiteMisfit",
    "PlotResponseGrid",
]

# Data-type code metadata: code → (short_name, colorbar_label, is_rho)
_TYPE_INFO = {
    1: ("RhoTE", r"$\rho_a$ TE (Ω·m)",  True),
    2: ("PhsTE", "Phase TE (°)",          False),
    5: ("RhoTM", r"$\rho_a$ TM (Ω·m)",  True),
    6: ("PhsTM", "Phase TM (°)",          False),
}
_RHO_CODES = frozenset({1, 5})
_TE_RHO, _TE_PHS = 1, 2
_TM_RHO, _TM_PHS = 5, 6


def _edges(centers: np.ndarray) -> np.ndarray:
    """Return cell-edge positions from 1-D cell-centre array."""
    c = np.asarray(centers, dtype=float)
    e = np.empty(len(c) + 1)
    if len(c) > 1:
        e[1:-1] = (c[:-1] + c[1:]) / 2.0
        e[0]    = c[0]  - (c[1]  - c[0])  / 2.0
        e[-1]   = c[-1] + (c[-1] - c[-2]) / 2.0
    elif len(c) == 1:
        e[0], e[1] = c[0] - 0.5, c[0] + 0.5
    return e


class _OccamPlotBase(OccamBase):
    """Shared initialisation for all Occam2D plot classes."""

    def __init__(
        self,
        result=None,
        figsize: Optional[Tuple[float, float]] = None,
        cmap: str = "jet_r",
        dpi: int = 100,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.result  = result
        self.figsize = figsize
        self.cmap    = cmap
        self.dpi     = dpi

    def plot(self):
        raise NotImplementedError(
            f"{self.__class__.__name__}.plot — not yet implemented"
        )


# ---------------------------------------------------------------------------
# PlotMisfit
# ---------------------------------------------------------------------------

class PlotMisfit(_OccamPlotBase):
    r"""Plot Occam2D convergence metrics by iteration.

    ``PlotMisfit`` visualizes the convergence history stored
    in an :class:`~pycsamt.models.occam2d.OccamLog`
    attached to an
    :class:`~pycsamt.models.occam2d.InversionResult`.
    It replaces the MATLAB ``plotOccamIterMisfit.m`` view.

    The main curve is the normalized root-mean-square data
    misfit:

    .. math::

        \mathrm{RMS} =
        \sqrt{\frac{1}{N}\sum_{i=1}^{N} r_i^2},

    Here :math:`r_i` is the weighted residual for datum ``i``.
    A run is commonly acceptable when the RMS approaches the
    target value of 1.0 [1]_.

    Parameters
    ----------
    result : InversionResult
        Loaded inversion result. It must expose ``log`` with
        ``iterations``, ``rms``, ``roughness``, ``lagrange``,
        and ``n_iter`` attributes.
    show_roughness : bool, default True
        If ``True``, add a secondary y-axis for roughness.
        Roughness is plotted on a log scale because it can
        vary by several orders of magnitude.
    show_lagrange : bool, default False
        If ``True``, add a lower panel for the accepted
        Lagrange multiplier at each iteration.
    target_line : bool, default True
        If ``True``, draw a dashed line at RMS equal
        to 1.0.
    figsize : tuple of float, optional
        Matplotlib figure size passed through the shared base.
        If omitted, the number of panels controls the size.
    cmap : str, default "jet_r"
        Stored for consistency with other plot classes.
        It is not used by this line plot.
    dpi : int, default 100
        Figure resolution in dots per inch.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the convergence plot.

    Raises
    ------
    RuntimeError
        If ``result.log`` is missing or has no iterations.

    See Also
    --------
    OccamLog
        Parses convergence values from the Occam log file.
    InversionResult.plot_misfit
        Convenience wrapper that instantiates this class.

    Examples
    --------
    >>> from pycsamt.models.occam2d import InversionResult
    >>> from pycsamt.models.occam2d import PlotMisfit
    >>> result = InversionResult("occam_run")
    >>> fig = PlotMisfit(result, show_lagrange=True).plot()

    References
    ----------
    .. [1] Constable, S. C., Parker, R. L., and Constable,
       C. G., "Occam's inversion: A practical algorithm for
       generating smooth models from electromagnetic sounding
       data", Geophysics, 52(3), 289-300, 1987.
    """

    def __init__(
        self,
        result=None,
        show_roughness: bool = True,
        show_lagrange:  bool = False,
        target_line:    bool = True,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.show_roughness = show_roughness
        self.show_lagrange  = show_lagrange
        self.target_line    = target_line

    def plot(self):
        """Return a convergence Figure (RMS ± roughness vs iteration).

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker

        log = getattr(self.result, "log", None)
        if log is None or log.n_iter == 0:
            raise RuntimeError("PlotMisfit: no log data available")

        iters = log.iterations
        rms   = log.rms
        rough = log.roughness
        lagr  = log.lagrange

        n_panels = 2 if self.show_lagrange else 1
        fig, axes = plt.subplots(
            n_panels, 1,
            figsize=self.figsize or (8, 4 * n_panels),
            dpi=self.dpi,
            squeeze=False,
        )
        ax_rms = axes[0, 0]

        ax_rms.plot(iters, rms, "o-", color="C0", label="RMS misfit", zorder=3)
        if self.target_line:
            ax_rms.axhline(1.0, color="k", ls="--", lw=0.8, label="Target  1.0")
        ax_rms.set_xlabel("Iteration")
        ax_rms.set_ylabel("RMS misfit", color="C0")
        ax_rms.tick_params(axis="y", labelcolor="C0")
        ax_rms.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax_rms.set_title("Occam2D convergence")
        ax_rms.legend(loc="upper right", fontsize="small")

        if self.show_roughness and rough.size > 0:
            valid = np.isfinite(rough)
            if valid.any():
                ax2 = ax_rms.twinx()
                ax2.semilogy(
                    iters[valid], rough[valid], "s--",
                    color="C1", label="Roughness", zorder=2,
                )
                ax2.set_ylabel("Roughness", color="C1")
                ax2.tick_params(axis="y", labelcolor="C1")

        if self.show_lagrange and lagr.size > 0:
            ax_lag = axes[1, 0]
            valid  = np.isfinite(lagr)
            if valid.any():
                ax_lag.semilogy(iters[valid], lagr[valid], "^-", color="C2")
            ax_lag.set_xlabel("Iteration")
            ax_lag.set_ylabel("Lagrange multiplier")
            ax_lag.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# PlotModel
# ---------------------------------------------------------------------------

class PlotModel(_OccamPlotBase):
    r"""Plot a two-dimensional Occam resistivity model.

    ``PlotModel`` displays the selected iteration model from
    :class:`~pycsamt.models.occam2d.InversionResult` as a
    depth section. It replaces ``plotOccam2DMT.m``, and the
    companion profile extractor replaces
    ``ExtractOccam2DMTProfile.m``.

    Occam iteration files store model parameters as
    :math:`\log_{10}` resistivity. The plot converts them back
    to ohm metres before drawing:

    .. math::

        \rho(x, z) = 10^{m(x, z)}.

    The mesh is centered around the profile midpoint.
    Depth is positive downward.

    Parameters
    ----------
    result : InversionResult
        Loaded result containing ``rho_2d`` and ``mesh``.
        Station markers are drawn when ``result.data.offsets``
        is available.
    rho_min, rho_max : float, default 1.0, 1000.0
        Color-scale limits in ohm metres. They are passed to
        :class:`matplotlib.colors.LogNorm`; both values must
        be positive.
    depth_max : float, optional
        Maximum display depth in metres. If omitted, the full
        mesh depth is shown.
    show_stations : bool, default True
        If ``True``, overlay station triangles at the surface
        using the offsets stored in the data file.
    profile_distance_unit : {"m", "km"}, default "km"
        Unit used on the horizontal axis and by
        :meth:`extract_profile`.
    figsize : tuple of float, optional
        Matplotlib figure size. The default suits a profile
        section.
    cmap : str, default "jet_r"
        Colormap used for the resistivity image.
    dpi : int, default 100
        Figure resolution in dots per inch.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the resistivity model section.

    Raises
    ------
    RuntimeError
        If the result does not contain ``rho_2d``.

    Methods
    -------
    plot()
        Return the model section as a Matplotlib figure.
    extract_profile(x0, x1)
        Extract centered coordinates, depth centers, and the
        log10-resistivity grid between ``x0`` and ``x1``.

    See Also
    --------
    InversionResult
        Reconstructs ``rho_2d`` from Occam output files.
    PlotSounding1D
        Extracts model columns at station positions.

    Examples
    --------
    >>> from pycsamt.models.occam2d import InversionResult
    >>> from pycsamt.models.occam2d import PlotModel
    >>> result = InversionResult("occam_run")
    >>> fig = PlotModel(result, depth_max=2000).plot()
    >>> plotter = PlotModel(result)
    >>> x, z, log_rho = plotter.extract_profile(-1.0, 1.0)

    References
    ----------
    .. [1] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    """

    def __init__(
        self,
        result=None,
        rho_min: float = 1.0,
        rho_max: float = 1000.0,
        depth_max: Optional[float] = None,
        show_stations: bool = True,
        profile_distance_unit: str = "km",
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.rho_min               = rho_min
        self.rho_max               = rho_max
        self.depth_max             = depth_max
        self.show_stations         = show_stations
        self.profile_distance_unit = profile_distance_unit

    def plot(self):
        """Return a pcolormesh Figure of the 2-D resistivity model.

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        result = self.result
        if result is None or result.rho_2d is None:
            raise RuntimeError("PlotModel: no rho_2d available")

        mesh = result.mesh
        rho  = result.rho_2d   # (n_zcells, n_xcells), log10(Ω·m)

        with np.errstate(over="ignore"):
            rho_lin = np.where(np.isfinite(rho), 10.0 ** rho, np.nan)

        xn = mesh.x_nodes.copy()
        zn = mesh.z_nodes.copy()
        xc = (xn[:-1] + xn[1:]) / 2.0
        x_shift = float(xc.mean())
        xn_c    = xn - x_shift

        if self.profile_distance_unit == "km":
            xn_plot     = xn_c / 1000.0
            x_off_scale = 1000.0
            xlabel      = "Distance (km)"
        else:
            xn_plot     = xn_c
            x_off_scale = 1.0
            xlabel      = "Distance (m)"

        depth_max = self.depth_max if self.depth_max is not None else float(zn[-1])
        zi_max    = int(np.searchsorted(zn, depth_max))
        zi_max    = min(zi_max + 1, len(zn))

        figsize = self.figsize or (12, 5)
        fig, ax  = plt.subplots(figsize=figsize, dpi=self.dpi)

        rho_sub = np.ma.masked_invalid(rho_lin[:zi_max - 1, :])
        im = ax.pcolormesh(
            xn_plot, zn[:zi_max], rho_sub,
            cmap=self.cmap,
            norm=LogNorm(vmin=self.rho_min, vmax=self.rho_max),
            shading="flat",
        )

        ax.set_ylim(float(zn[zi_max - 1]), 0.0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Depth (m)")
        iter_n = getattr(getattr(result, "best_iter", None), "iteration", "?")
        ax.set_title(f"Occam2D resistivity model  [iter {iter_n}]")

        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label("Resistivity (Ω·m)")

        if self.show_stations:
            data_obj = getattr(result, "data", None)
            if data_obj is not None and data_obj.offsets.size > 0:
                sta_x = (data_obj.offsets - x_shift) / x_off_scale
                ax.plot(sta_x, np.zeros_like(sta_x), "v", color="k",
                        ms=4, clip_on=False, zorder=5)

        fig.tight_layout()
        return fig

    def extract_profile(
        self, x0: float, x1: float
    ) -> tuple:
        """Extract (x_centers, z_centers, rho_subset) between x0 and x1.

        Coordinates use ``profile_distance_unit``.
        ``x0`` and ``x1`` are in the centered profile frame.

        Returns
        -------
        tuple
            ``(x_centers, z_centers, rho_2d_subset)`` where
            ``rho_2d_subset`` has shape
            ``(n_zcells, n_cols_in_range)``.
        """
        result = self.result
        if result is None or result.rho_2d is None:
            raise RuntimeError("extract_profile: no rho_2d available")

        mesh = result.mesh
        xn   = mesh.x_nodes
        zn   = mesh.z_nodes
        xc   = (xn[:-1] + xn[1:]) / 2.0
        zc   = (zn[:-1] + zn[1:]) / 2.0

        xc_c  = xc - float(xc.mean())
        scale = 1000.0 if self.profile_distance_unit == "km" else 1.0
        xc_s  = xc_c / scale

        mask = (xc_s >= x0) & (xc_s <= x1)
        return xc_s[mask], zc, result.rho_2d[:, mask]


# ---------------------------------------------------------------------------
# PlotResponse
# ---------------------------------------------------------------------------

class PlotResponse(_OccamPlotBase):
    r"""Plot observed and modeled Occam response curves.

    ``PlotResponse`` compares observed data from the Occam
    file with modeled values from an Occam ``.resp`` file. It
    replaces ``plotOccam2DMTResponse.m``.

    Apparent-resistivity rows are stored as log10 values.
    They are converted before plotting:

    .. math::

        \rho_a = 10^{d_\rho}.

    Phase rows are plotted in degrees. The frequency index in
    the table is mapped back to physical frequency when
    corresponding :class:`OccamData` object is available.

    Parameters
    ----------
    result : InversionResult
        Loaded result containing ``response`` and ideally
        ``data``. The response must expose the seven-column
        table.
    stations : list of str, list of int, or None, default None
        Stations to plot. Strings are matched against
        ``result.data.sites``. Integers are one-based site
        indices when names are unavailable. If ``None``,
        stations are sampled from the response table.
    modes : list of str, optional
        Electromagnetic modes to draw. Supported values are
        ``"TE"`` and ``"TM"``.
    period_axis : bool, default True
        Reserved for interface clarity. The implementation
        uses period when frequencies are available and indices
        otherwise.
    max_stations : int, default 9
        Maximum number of station columns when ``stations`` is
        ``None``.
    figsize : tuple of float, optional
        Figure size. If omitted, width scales with station
        count.
    cmap : str, default "jet_r"
        Stored for a consistent interface. It is not used
        by this curve plot.
    dpi : int, default 100
        Figure resolution in dots per inch.

    Returns
    -------
    matplotlib.figure.Figure
        Figure with apparent-resistivity and phase panels.

    Raises
    ------
    RuntimeError
        If response data are missing, modes are absent,
        or no stations can be selected.

    See Also
    --------
    PlotResponseGrid
        Compact version designed for many stations.
    OccamResponse
        Reader for the ``.resp`` file used by this plot.

    Examples
    --------
    >>> from pycsamt.models.occam2d import InversionResult
    >>> from pycsamt.models.occam2d import PlotResponse
    >>> result = InversionResult("occam_run")
    >>> fig = PlotResponse(result, modes=["TM"]).plot()

    References
    ----------
    .. [1] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    """

    def __init__(
        self,
        result=None,
        stations=None,
        modes: Optional[List[str]] = None,
        period_axis: bool = True,
        max_stations: int = 9,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.stations     = stations
        self.modes        = modes or ["TE", "TM"]
        self.period_axis  = period_axis
        self.max_stations = max_stations

    def plot(self):
        """Return a Figure with rho_a / phase subplots per station.

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt

        result   = self.result
        response = getattr(result, "response", None)
        data_obj = getattr(result, "data", None)

        if response is None or response.n_data == 0:
            raise RuntimeError("PlotResponse: no response data available")

        # Build mode → {rho_code, phs_code}
        _mode_map = {
            "TE": {"rho": _TE_RHO, "phs": _TE_PHS},
            "TM": {"rho": _TM_RHO, "phs": _TM_PHS},
        }
        mode_codes = {m: _mode_map[m] for m in self.modes if m in _mode_map}

        present = set(response.data[:, 2].astype(int))
        mode_codes = {
            m: c for m, c in mode_codes.items()
            if c["rho"] in present or c["phs"] in present
        }
        if not mode_codes:
            raise RuntimeError(
                "PlotResponse: none of the requested modes found in response"
            )

        # Station index selection
        all_site_idx = np.unique(response.data[:, 0].astype(int))
        if self.stations is not None:
            if data_obj is not None and data_obj.sites:
                idx_map  = {name: i + 1 for i, name in enumerate(data_obj.sites)}
                sel_idx  = np.array(
                    [idx_map[s] for s in self.stations if s in idx_map], dtype=int
                )
            else:
                sel_idx = np.asarray(self.stations, dtype=int)
        else:
            sel_idx = all_site_idx

        if len(sel_idx) > self.max_stations:
            step    = max(1, len(sel_idx) // self.max_stations)
            sel_idx = sel_idx[::step][: self.max_stations]

        n_sel = len(sel_idx)
        if n_sel == 0:
            raise RuntimeError("PlotResponse: no stations to plot")

        fig, axes = plt.subplots(
            2, n_sel,
            figsize=self.figsize or (4 * n_sel, 6),
            dpi=self.dpi,
            squeeze=False,
        )

        colors = {"TE": "C0", "TM": "C1"}

        if data_obj is not None and data_obj.frequencies.size > 0:
            freq_arr = data_obj.frequencies
        else:
            freq_arr = None

        for si, site_idx in enumerate(sel_idx):
            ax_rho = axes[0, si]
            ax_phs = axes[1, si]

            site_mask = response.data[:, 0].astype(int) == site_idx

            for mode_name, codes in mode_codes.items():
                color = colors.get(mode_name, "C0")

                for ax, key, code in [
                    (ax_rho, "rho", codes["rho"]),
                    (ax_phs, "phs", codes["phs"]),
                ]:
                    mask = site_mask & (response.data[:, 2].astype(int) == code)
                    if not mask.any():
                        continue
                    sub = response.data[mask]
                    fi  = sub[:, 1].astype(int)   # 1-based frequency indices
                    obs = sub[:, 4]
                    mod = sub[:, 5]

                    if freq_arr is not None:
                        idx0   = np.clip(fi - 1, 0, len(freq_arr) - 1)
                        x_vals = 1.0 / freq_arr[idx0]
                    else:
                        x_vals = fi.astype(float)

                    order   = np.argsort(x_vals)
                    x_vals  = x_vals[order]
                    obs     = obs[order]
                    mod     = mod[order]

                    if key == "rho":
                        obs_v = np.where(np.isfinite(obs), 10.0 ** obs, np.nan)
                        mod_v = np.where(np.isfinite(mod), 10.0 ** mod, np.nan)
                        ax.loglog(x_vals, obs_v, "o", ms=4, color=color,
                                  label=f"obs {mode_name}")
                        ax.loglog(x_vals, mod_v, "-", lw=1, color=color,
                                  label=f"mod {mode_name}")
                    else:
                        ax.semilogx(x_vals, obs, "o", ms=4, color=color,
                                    label=f"obs {mode_name}")
                        ax.semilogx(x_vals, mod, "-", lw=1, color=color,
                                    label=f"mod {mode_name}")

            site_name = (
                data_obj.sites[site_idx - 1]
                if (data_obj is not None and data_obj.sites
                    and site_idx - 1 < len(data_obj.sites))
                else f"S{site_idx:03d}"
            )
            ax_rho.set_title(site_name, fontsize="small")
            ax_rho.set_ylabel(r"$\rho_a$ (Ω·m)")
            ax_phs.set_ylabel("Phase (°)")
            x_label = "Period (s)" if freq_arr is not None else "Freq. index"
            ax_phs.set_xlabel(x_label)
            if si == 0:
                ax_rho.legend(fontsize="x-small", loc="best")

        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# PlotPseudo
# ---------------------------------------------------------------------------

class PlotPseudo(_OccamPlotBase):
    r"""Plot an Occam observed-data pseudosection.

    ``PlotPseudo`` displays one data component from the Occam
    data file as a station-period view. It is the Python
    replacement for ``plotOccam2DMTPseudo.m``.

    The horizontal axis is station offset in kilometres when
    offsets are available. The vertical axis is
    :math:`\log_{10}(T)`, where :math:`T = 1/f` is period in
    seconds. Apparent resistivity is converted from log10
    storage to ohm metres; phase data remain in degrees.

    Parameters
    ----------
    result : InversionResult
        Loaded result containing an :class:`OccamData` object
        with data blocks, offsets, and frequencies.
    mode : {"TE", "TM"}, default "TM"
        Electromagnetic mode to display. ``"TE"`` maps to
        codes 1 and 2; ``"TM"`` maps to codes 5 and 6.
    data_type : {"rho", "phase"}, default "rho"
        Quantity to display. ``"rho"`` selects apparent
        resistivity; ``"phase"`` selects phase.
    figsize : tuple of float, optional
        Matplotlib figure size.
    cmap : str, default "jet_r"
        Colormap used for the pseudosection.
    dpi : int, default 100
        Figure resolution in dots per inch.

    Returns
    -------
    matplotlib.figure.Figure
        Pseudosection figure.

    Raises
    ------
    RuntimeError
        If no data blocks are available or the selected type
        is not present.
    ValueError
        If ``mode`` and ``data_type`` are unsupported.

    See Also
    --------
    OccamData
        Provides the data block for the pseudosection.
    PlotSiteMisfit
        Builds a residual pseudosection from response values.

    Examples
    --------
    >>> from pycsamt.models.occam2d import InversionResult
    >>> from pycsamt.models.occam2d import PlotPseudo
    >>> result = InversionResult("occam_run")
    >>> fig = PlotPseudo(result, mode="TE").plot()
    """

    def __init__(
        self,
        result=None,
        mode: str = "TM",
        data_type: str = "rho",
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.mode      = mode
        self.data_type = data_type

    def plot(self):
        """Return a pseudosection pcolormesh Figure.

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        result   = self.result
        data_obj = getattr(result, "data", None)
        if data_obj is None or data_obj.data_blocks.size == 0:
            raise RuntimeError("PlotPseudo: no data blocks available")

        _code_map = {
            ("TE", "rho"):   1,
            ("TE", "phase"): 2,
            ("TM", "rho"):   5,
            ("TM", "phase"): 6,
        }
        code = _code_map.get((self.mode.upper(), self.data_type.lower()))
        if code is None:
            raise ValueError(
                f"Unknown mode/data_type: ({self.mode!r}, {self.data_type!r})"
            )

        db   = data_obj.data_blocks
        mask = db[:, 2].astype(int) == code
        if not mask.any():
            raise RuntimeError(
                f"PlotPseudo: no data with type code {code} "
                f"(mode={self.mode}, data_type={self.data_type})"
            )
        sub      = db[mask]
        site_idx = sub[:, 0].astype(int) - 1   # 0-based
        freq_idx = sub[:, 1].astype(int) - 1   # 0-based
        values   = sub[:, 3]

        n_sites = data_obj.n_sites
        n_freqs = data_obj.n_frequencies
        grid    = np.full((n_freqs, n_sites), np.nan)
        ok      = (
            (site_idx >= 0) & (site_idx < n_sites) &
            (freq_idx >= 0) & (freq_idx < n_freqs)
        )
        grid[freq_idx[ok], site_idx[ok]] = values[ok]

        is_rho = code in _RHO_CODES
        if is_rho:
            with np.errstate(over="ignore"):
                grid_plot = np.where(np.isfinite(grid), 10.0 ** grid, np.nan)
        else:
            grid_plot = grid.copy()

        # x: station offsets converted to km
        if data_obj.offsets.size == n_sites:
            x_km = data_obj.offsets / 1000.0
        else:
            x_km = np.arange(n_sites, dtype=float)

        # y: log10(period)
        if data_obj.frequencies.size == n_freqs:
            periods     = 1.0 / data_obj.frequencies
            log_periods = np.log10(periods)
        else:
            log_periods = np.arange(n_freqs, dtype=float)

        x_edges = _edges(x_km)
        y_edges = _edges(log_periods)

        figsize = self.figsize or (12, 5)
        fig, ax  = plt.subplots(figsize=figsize, dpi=self.dpi)

        masked = np.ma.masked_invalid(grid_plot)
        if is_rho:
            finite = grid_plot[np.isfinite(grid_plot)]
            vmin   = float(finite.min()) if finite.size else 1.0
            vmax   = float(finite.max()) if finite.size else 1000.0
            norm   = LogNorm(vmin=max(vmin, 1e-6), vmax=max(vmax, vmin * 10.0))
            im = ax.pcolormesh(
                x_edges, y_edges, masked,
                cmap=self.cmap, norm=norm, shading="flat",
            )
        else:
            im = ax.pcolormesh(
                x_edges, y_edges, masked,
                cmap=self.cmap, shading="flat",
            )

        _, cb_label, _ = _TYPE_INFO.get(code, ("", str(code), False))
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label(cb_label)

        ax.set_xlabel("Distance (km)")
        ax.set_ylabel(r"$\log_{10}$(Period [s])")
        ax.set_title(
            f"Pseudosection — {self.mode} {self.data_type}  "
            f"({n_sites} sites, {n_freqs} freqs)"
        )

        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# PlotSounding1D
# ---------------------------------------------------------------------------

class PlotSounding1D(_OccamPlotBase):
    r"""Plot station-centered 1-D profiles from a 2-D Occam model.

    ``PlotSounding1D`` samples the reconstructed 2-D
    resistivity grid at the mesh column nearest each station.
    The result is a set of resistivity-depth curves to compare
    vertical structure below stations.

    The plotted resistivity is converted from the log10 grid:

    .. math::

        \rho(z) = 10^{m(z)}.

    Air layers are omitted using ``result.mesh.n_airlayers``.

    Parameters
    ----------
    result : InversionResult
        Loaded result containing ``rho_2d``, ``mesh``, and
        ``data.offsets``.
    stations : list of str or None, default None
        Station names to plot. If ``None``, stations are
        sampled from all available offsets.
    max_stations : int, default 16
        Maximum station profiles when ``stations`` is
        ``None``.
    depth_max : float, optional
        Maximum plotted depth. If omitted, mesh depth
        controls the lower limit.
    rho_min, rho_max : float, default 1.0, 1000.0
        Horizontal resistivity-axis limits in ohm metres.
    overlay : bool, default False
        If ``True``, draw selected profiles on one axis. If
        ``False``, draw one panel per station.
    figsize : tuple of float, optional
        Figure size. Defaults depend on station count.
    cmap : str, default "jet_r"
        Stored for interface consistency. Overlay plots use a
        tabular Matplotlib colormap internally.
    dpi : int, default 100
        Figure resolution in dots per inch.

    Returns
    -------
    matplotlib.figure.Figure
        Resistivity-depth profile figure.

    Raises
    ------
    RuntimeError
        If ``rho_2d`` is missing, station offsets are missing,
        or station selection is empty.

    See Also
    --------
    PlotModel.extract_profile
        Extracts a horizontal interval from the same grid.
    InversionResult
        Provides the reconstructed model grid.

    Examples
    --------
    >>> from pycsamt.models.occam2d import InversionResult
    >>> from pycsamt.models.occam2d import PlotSounding1D
    >>> result = InversionResult("occam_run")
    >>> fig = PlotSounding1D(result, overlay=True).plot()
    """

    def __init__(
        self,
        result=None,
        stations=None,
        max_stations: int = 16,
        depth_max: Optional[float] = None,
        rho_min: float = 1.0,
        rho_max: float = 1000.0,
        overlay: bool = False,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.stations     = stations
        self.max_stations = max_stations
        self.depth_max    = depth_max
        self.rho_min      = rho_min
        self.rho_max      = rho_max
        self.overlay      = overlay

    def plot(self):
        """Return a Figure of 1-D ρ–depth soundings.

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm

        result = self.result
        if result is None or result.rho_2d is None:
            raise RuntimeError("PlotSounding1D: no rho_2d available")

        mesh     = result.mesh
        data_obj = getattr(result, "data", None)
        rho_2d   = result.rho_2d   # (n_zcells, n_xcells), log10(Ω·m)

        xn = mesh.x_nodes
        zn = mesh.z_nodes
        xc = (xn[:-1] + xn[1:]) / 2.0
        zc = (zn[:-1] + zn[1:]) / 2.0

        if data_obj is None or data_obj.offsets.size == 0:
            raise RuntimeError("PlotSounding1D: no station offsets in result.data")

        offsets    = data_obj.offsets
        site_names = (
            list(data_obj.sites)
            if data_obj.sites else
            [f"S{i+1:03d}" for i in range(len(offsets))]
        )

        if self.stations is not None:
            sel = [(n, o) for n, o in zip(site_names, offsets) if n in self.stations]
        else:
            sel = list(zip(site_names, offsets))

        if len(sel) > self.max_stations:
            step = max(1, len(sel) // self.max_stations)
            sel  = sel[::step][: self.max_stations]
        if not sel:
            raise RuntimeError("PlotSounding1D: no stations to plot")

        n_air     = mesh.n_airlayers
        depth_max = self.depth_max if self.depth_max is not None else float(zn[-1])
        zi_max    = min(int(np.searchsorted(zn, depth_max)) + 1, len(zn) - 1)
        zc_sub    = zc[n_air:zi_max]

        iter_n = getattr(getattr(result, "best_iter", None), "iteration", "?")

        if self.overlay:
            fig, ax = plt.subplots(
                figsize=self.figsize or (6, 8), dpi=self.dpi
            )
            colors = cm.tab20(np.linspace(0, 1, len(sel)))
            for (name, off), color in zip(sel, colors):
                col_idx = int(np.argmin(np.abs(xc - off)))
                prof    = rho_2d[n_air:zi_max, col_idx]
                rho_lin = np.where(np.isfinite(prof), 10.0 ** prof, np.nan)
                ax.semilogx(rho_lin, zc_sub, "-", lw=1.5, color=color, label=name)

            ax.set_xlim(self.rho_min, self.rho_max)
            ax.set_ylim(float(zc_sub[-1]) if zc_sub.size else 1.0, 0.0)
            ax.set_xlabel("Resistivity (Ω·m)")
            ax.set_ylabel("Depth (m)")
            ax.legend(fontsize="x-small", bbox_to_anchor=(1.02, 1), loc="upper left")
            ax.set_title(f"1-D soundings from 2-D model  [iter {iter_n}]")
            ax.grid(True, which="both", lw=0.4, alpha=0.4)
            fig.tight_layout()
            return fig

        n_sta  = len(sel)
        n_cols = min(4, n_sta)
        n_rows = (n_sta + n_cols - 1) // n_cols

        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=self.figsize or (3.5 * n_cols, 4.0 * n_rows),
            dpi=self.dpi,
            squeeze=False,
        )

        for si, (name, off) in enumerate(sel):
            row, col = divmod(si, n_cols)
            ax       = axes[row, col]
            col_idx  = int(np.argmin(np.abs(xc - off)))
            prof     = rho_2d[n_air:zi_max, col_idx]
            rho_lin  = np.where(np.isfinite(prof), 10.0 ** prof, np.nan)

            ax.semilogx(rho_lin, zc_sub, "-", lw=1.5, color="C0")
            ax.set_xlim(self.rho_min, self.rho_max)
            ax.set_ylim(float(zc_sub[-1]) if zc_sub.size else 1.0, 0.0)
            ax.set_title(name, fontsize="small")
            if col == 0:
                ax.set_ylabel("Depth (m)", fontsize="small")
            if row == n_rows - 1:
                ax.set_xlabel("ρ (Ω·m)", fontsize="small")
            ax.grid(True, which="both", lw=0.4, alpha=0.4)

        for si in range(n_sta, n_rows * n_cols):
            row, col = divmod(si, n_cols)
            axes[row, col].set_visible(False)

        fig.suptitle(
            f"1-D soundings from 2-D model  [iter {iter_n}]", fontsize="medium"
        )
        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# PlotSiteMisfit
# ---------------------------------------------------------------------------

class PlotSiteMisfit(_OccamPlotBase):
    r"""Plot per-site Occam response misfit diagnostics.

    ``PlotSiteMisfit`` summarizes the fit between observed and
    modeled values at each station. The top panel is a
    bar chart of RMS residual by station and data type. The
    optional lower panel is a residual pseudosection.

    Residuals are normalized by the error column in the Occam
    response table:

    .. math::

        r_i = \frac{d_i^{obs} - d_i^{pred}}{\sigma_i}.

    Per-site RMS values use these normalized residuals. If a
    response error is non-positive, that residual is ignored.

    Parameters
    ----------
    result : InversionResult
        Loaded result containing ``response``. Station labels
        and frequencies come from ``result.data``.
    modes : list of str, optional
        Modes included in the summary. Supported values are
        ``"TE"`` and ``"TM"``. If omitted, both are used.
    show_residual_map : bool, default True
        If ``True``, draw the normalized residual map
        under the bar chart.
    rms_target : float, default 1.0
        Target RMS value drawn in the bar panel.
        Use ``None`` to omit the target line.
    figsize : tuple of float, optional
        Figure size. Defaults scale with station count.
    cmap : str, default "jet_r"
        Stored for the shared interface. The residual map uses
        a diverging colormap.
    dpi : int, default 100
        Figure resolution in dots per inch.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing per-site RMS diagnostics.

    Raises
    ------
    RuntimeError
        If response data are missing or requested type codes
        are absent.

    See Also
    --------
    OccamResponse.misfit_per_site
        Returns a simpler per-site RMS dictionary.
    PlotResponseGrid
        Shows observed and modeled curves for many stations.

    Examples
    --------
    >>> from pycsamt.models.occam2d import InversionResult
    >>> from pycsamt.models.occam2d import PlotSiteMisfit
    >>> result = InversionResult("occam_run")
    >>> fig = PlotSiteMisfit(result).plot()
    """

    def __init__(
        self,
        result=None,
        modes: Optional[List[str]] = None,
        show_residual_map: bool = True,
        rms_target: float = 1.0,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.modes             = modes or ["TE", "TM"]
        self.show_residual_map = show_residual_map
        self.rms_target        = rms_target

    def plot(self):
        """Return a per-site misfit Figure.

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt

        result   = self.result
        response = getattr(result, "response", None)
        data_obj = getattr(result, "data", None)

        if response is None or response.n_data == 0:
            raise RuntimeError("PlotSiteMisfit: no response data available")

        # response columns: [site(1-b), freq(1-b), type_code, err, obs, pred]
        rd      = response.data
        site_1b = rd[:, 0].astype(int)
        freq_1b = rd[:, 1].astype(int)
        tcodes  = rd[:, 2].astype(int)
        errors  = rd[:, 3]
        obs     = rd[:, 4]
        pred    = rd[:, 5]

        safe_e = np.where(errors > 0, errors, 1.0)
        resid  = np.where(errors > 0, (obs - pred) / safe_e, np.nan)

        _mode_codes = {
            "TE": [_TE_RHO, _TE_PHS],
            "TM": [_TM_RHO, _TM_PHS],
        }
        present      = set(np.unique(tcodes))
        wanted_codes = [
            c for m in self.modes
            for c in _mode_codes.get(m, [])
            if c in present
        ]
        if not wanted_codes:
            raise RuntimeError("PlotSiteMisfit: no matching data type codes")

        n_sites     = int(site_1b.max())
        site_labels = (
            list(data_obj.sites)[:n_sites]
            if (data_obj and data_obj.sites) else
            [f"S{i+1:03d}" for i in range(n_sites)]
        )

        rms_table = np.full((n_sites, len(wanted_codes)), np.nan)
        for ci, code in enumerate(wanted_codes):
            for si in range(n_sites):
                m   = (site_1b == si + 1) & (tcodes == code)
                r   = resid[m]
                fin = r[np.isfinite(r)]
                if fin.size:
                    rms_table[si, ci] = float(np.sqrt(np.mean(fin ** 2)))

        code_labels = [_TYPE_INFO.get(c, (str(c), "", False))[0] for c in wanted_codes]
        _bar_colors = ["C0", "C0", "C1", "C1"]
        _bar_alpha  = [1.0,  0.55, 1.0,  0.55]

        n_panels = 2 if self.show_residual_map else 1
        fig, axes = plt.subplots(
            n_panels, 1,
            figsize=self.figsize or (max(10, n_sites * 0.7), 5 * n_panels),
            dpi=self.dpi,
            squeeze=False,
        )

        # ---- top: bar chart ----
        ax_bar  = axes[0, 0]
        n_codes = len(wanted_codes)
        x_pos   = np.arange(n_sites)
        bar_w   = 0.8 / n_codes

        for ci, (label, clr, alp) in enumerate(
            zip(code_labels, _bar_colors[:n_codes], _bar_alpha[:n_codes])
        ):
            offset = (ci - n_codes / 2.0 + 0.5) * bar_w
            ax_bar.bar(
                x_pos + offset, rms_table[:, ci],
                width=bar_w, label=label,
                color=clr, alpha=alp, edgecolor="none",
            )

        overall = np.sqrt(np.nanmean(rms_table ** 2, axis=1))
        ax_bar.plot(x_pos, overall, "kD", ms=5, zorder=5, label="Overall")

        if self.rms_target is not None:
            ax_bar.axhline(
                self.rms_target, color="k", ls="--", lw=0.8,
                label=f"Target {self.rms_target:.1f}",
            )

        ax_bar.set_xticks(x_pos)
        ax_bar.set_xticklabels(
            site_labels, rotation=45, ha="right", fontsize="x-small"
        )
        ax_bar.set_ylabel("RMS misfit")
        ax_bar.set_title("Per-site RMS misfit by data type")
        ax_bar.legend(fontsize="small", loc="upper right")

        if not self.show_residual_map:
            fig.tight_layout()
            return fig

        # ---- bottom: normalised-residual pseudosection ----
        ax_map = axes[1, 0]

        if data_obj is not None and data_obj.frequencies.size > 0:
            n_freqs     = data_obj.n_frequencies
            freq_arr    = data_obj.frequencies
            log_periods = np.log10(1.0 / freq_arr)
        else:
            n_freqs     = int(freq_1b.max())
            log_periods = np.arange(n_freqs, dtype=float)

        x_km = (
            data_obj.offsets / 1000.0
            if (data_obj is not None and data_obj.offsets.size > 0)
            else np.arange(n_sites, dtype=float)
        )

        resid_sum = np.zeros((n_freqs, n_sites))
        count_map = np.zeros((n_freqs, n_sites), dtype=int)

        for code in wanted_codes:
            cmask = tcodes == code
            s_idx = site_1b[cmask] - 1
            f_idx = freq_1b[cmask] - 1
            rv    = resid[cmask]
            ok    = (
                np.isfinite(rv) &
                (s_idx >= 0) & (s_idx < n_sites) &
                (f_idx >= 0) & (f_idx < n_freqs)
            )
            np.add.at(resid_sum, (f_idx[ok], s_idx[ok]), rv[ok])
            np.add.at(count_map, (f_idx[ok], s_idx[ok]), 1)

        with np.errstate(invalid="ignore"):
            resid_map = np.where(count_map > 0, resid_sum / count_map, np.nan)

        x_edges = _edges(x_km)
        y_edges = _edges(log_periods)
        vmax    = max(3.0, float(np.nanpercentile(np.abs(resid_map), 95)))

        im = ax_map.pcolormesh(
            x_edges, y_edges,
            np.ma.masked_invalid(resid_map),
            cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="flat",
        )
        cb = fig.colorbar(im, ax=ax_map, pad=0.02)
        cb.set_label(r"(obs $-$ pred) / error  [$\sigma$]")

        ax_map.set_xlabel("Distance (km)")
        ax_map.set_ylabel(r"$\log_{10}$(Period [s])")
        ax_map.set_title("Normalised residuals (pseudosection)")

        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# PlotResponseGrid
# ---------------------------------------------------------------------------

class PlotResponseGrid(_OccamPlotBase):
    r"""Plot a compact grid of observed and modeled responses.

    ``PlotResponseGrid`` is designed to scan many stations.
    Each station uses two axes: apparent resistivity above and
    phase below. Observed values are drawn as points and
    modeled values as lines. Titles include per-site RMS when
    response errors are available.

    Apparent-resistivity values are converted from log10 Occam
    storage before plotting:

    .. math::

        \rho_a = 10^{d_\rho}.

    Parameters
    ----------
    result : InversionResult
        Loaded result with ``response`` and ideally ``data``.
    stations : list of str or None, default None
        Station names to include. If omitted, station indices
        sampled from the response table.
    n_cols : int, default 5
        Maximum number of station columns in each grid row.
    modes : list of str, optional
        Modes to draw. Use ``"TE"``, ``"TM"``, or both.
    max_stations : int, default 25
        Maximum station count included when ``stations`` is
        ``None``.
    figsize : tuple of float, optional
        Figure size. Defaults scale with columns and rows.
    cmap : str, default "jet_r"
        Stored for interface consistency. It is unused by this
        curve plot.
    dpi : int, default 100
        Figure resolution in dots per inch.

    Returns
    -------
    matplotlib.figure.Figure
        Compact response-grid figure.

    Raises
    ------
    RuntimeError
        If response data or station choices are missing.

    See Also
    --------
    PlotResponse
        Larger response panels for a smaller station subset.
    PlotSiteMisfit
        Per-site residual summary from the response table.

    Examples
    --------
    >>> from pycsamt.models.occam2d import InversionResult
    >>> from pycsamt.models.occam2d import PlotResponseGrid
    >>> result = InversionResult("occam_run")
    >>> fig = PlotResponseGrid(result, n_cols=4).plot()
    """

    def __init__(
        self,
        result=None,
        stations=None,
        n_cols: int = 5,
        modes: Optional[List[str]] = None,
        max_stations: int = 25,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.stations     = stations
        self.n_cols       = n_cols
        self.modes        = modes or ["TE", "TM"]
        self.max_stations = max_stations

    def plot(self):
        """Return a compact response-grid Figure.

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt

        result   = self.result
        response = getattr(result, "response", None)
        data_obj = getattr(result, "data", None)

        if response is None or response.n_data == 0:
            raise RuntimeError("PlotResponseGrid: no response data available")

        _mode_map = {
            "TE": {"rho": _TE_RHO, "phs": _TE_PHS},
            "TM": {"rho": _TM_RHO, "phs": _TM_PHS},
        }
        mode_codes = {m: _mode_map[m] for m in self.modes if m in _mode_map}

        rd      = response.data
        all_idx = np.unique(rd[:, 0].astype(int))

        if self.stations is not None and data_obj is not None and data_obj.sites:
            idx_map = {name: i + 1 for i, name in enumerate(data_obj.sites)}
            sel_idx = np.array(
                [idx_map[s] for s in self.stations if s in idx_map], dtype=int
            )
        else:
            sel_idx = all_idx

        if len(sel_idx) > self.max_stations:
            step    = max(1, len(sel_idx) // self.max_stations)
            sel_idx = sel_idx[::step][: self.max_stations]

        if len(sel_idx) == 0:
            raise RuntimeError("PlotResponseGrid: no stations to plot")

        freq_arr = (
            data_obj.frequencies
            if (data_obj is not None and data_obj.frequencies.size > 0)
            else None
        )

        errors = rd[:, 3]
        obs    = rd[:, 4]
        pred   = rd[:, 5]
        safe_e = np.where(errors > 0, errors, 1.0)
        resid  = np.where(errors > 0, (obs - pred) / safe_e, np.nan)

        def _site_rms(s1b: int) -> float:
            m = rd[:, 0].astype(int) == s1b
            r = resid[m]
            r = r[np.isfinite(r)]
            return float(np.sqrt(np.mean(r ** 2))) if r.size else np.nan

        n_sta  = len(sel_idx)
        n_cols = min(self.n_cols, n_sta)
        n_rows = (n_sta + n_cols - 1) // n_cols

        fig, axes = plt.subplots(
            2 * n_rows, n_cols,
            figsize=self.figsize or (3.0 * n_cols, 4.0 * n_rows),
            dpi=self.dpi,
            squeeze=False,
        )

        colors = {"TE": "C0", "TM": "C1"}

        for si, site_1b in enumerate(sel_idx):
            pg_row, col = divmod(si, n_cols)
            ax_rho = axes[2 * pg_row,     col]
            ax_phs = axes[2 * pg_row + 1, col]
            smask  = rd[:, 0].astype(int) == site_1b

            for mode_name, codes in mode_codes.items():
                color = colors.get(mode_name, "C0")
                for ax, key, code in [
                    (ax_rho, "rho", codes["rho"]),
                    (ax_phs, "phs", codes["phs"]),
                ]:
                    mask = smask & (rd[:, 2].astype(int) == code)
                    if not mask.any():
                        continue
                    sub = rd[mask]
                    fi  = sub[:, 1].astype(int)
                    ov  = sub[:, 4]
                    mv  = sub[:, 5]

                    if freq_arr is not None:
                        idx0   = np.clip(fi - 1, 0, len(freq_arr) - 1)
                        x_vals = 1.0 / freq_arr[idx0]
                    else:
                        x_vals = fi.astype(float)

                    order  = np.argsort(x_vals)
                    x_vals = x_vals[order]
                    ov     = ov[order]
                    mv     = mv[order]

                    if key == "rho":
                        ov_lin = np.where(np.isfinite(ov), 10.0 ** ov, np.nan)
                        mv_lin = np.where(np.isfinite(mv), 10.0 ** mv, np.nan)
                        ax.loglog(x_vals, ov_lin, ".", ms=3, color=color)
                        ax.loglog(x_vals, mv_lin, "-", lw=1, color=color)
                    else:
                        ax.semilogx(x_vals, ov, ".", ms=3, color=color)
                        ax.semilogx(x_vals, mv, "-", lw=1, color=color)

            rms_val = _site_rms(site_1b)
            name    = (
                data_obj.sites[site_1b - 1]
                if (data_obj and data_obj.sites
                    and site_1b - 1 < len(data_obj.sites))
                else f"S{site_1b:03d}"
            )
            rms_str = f"{rms_val:.2f}" if np.isfinite(rms_val) else "n/a"
            ax_rho.set_title(f"{name}  [{rms_str}]", fontsize="xx-small")
            ax_rho.tick_params(labelsize=6, pad=1)
            ax_phs.tick_params(labelsize=6, pad=1)
            if col == 0:
                ax_rho.set_ylabel(r"$\rho_a$", fontsize="x-small")
                ax_phs.set_ylabel("φ (°)",     fontsize="x-small")

        for si in range(n_sta, n_rows * n_cols):
            pg_row, col = divmod(si, n_cols)
            axes[2 * pg_row,     col].set_visible(False)
            axes[2 * pg_row + 1, col].set_visible(False)

        iter_n = getattr(getattr(result, "best_iter", None), "iteration", "?")
        fig.suptitle(f"Response grid  [iter {iter_n}]", fontsize="small")
        fig.tight_layout()
        return fig
