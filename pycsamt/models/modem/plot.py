# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Diagnostic plotting helpers for ModEM inversion results.

The classes in this module turn an
:class:`~pycsamt.models.modem.results.InversionResult` into
matplotlib figures for checking convergence, inspecting 2-D
and 3-D resistivity models, and viewing apparent-resistivity
or phase responses.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from .base import ModEmBase
from .results import InversionResult

__all__ = [
    "PlotMisfit",
    "PlotModel2D",
    "PlotModel3D",
    "PlotResponse",
    "PlotPseudo",
]


# ----------------------------------------------------------------------
# Shared base
# ----------------------------------------------------------------------

class _ModEmPlotBase(ModEmBase):
    """Shared result handling for ModEM plotting helpers."""

    def __init__(self, result: InversionResult | None = None, **kwargs):
        """Initialize a plotter with an optional result object."""
        super().__init__(**kwargs)
        self.result = result

    def _check_result(self) -> InversionResult:
        """Return the attached result or raise a clear error."""
        if self.result is None:
            msg = (
                "No InversionResult attached. Pass result= "
                "to the constructor."
            )
            raise ValueError(msg)
        return self.result


# ----------------------------------------------------------------------
# PlotMisfit
# ----------------------------------------------------------------------

class PlotMisfit(_ModEmPlotBase):
    """Plot RMS misfit as a function of inversion iteration.

    ``PlotMisfit`` visualizes the convergence history parsed
    from :class:`~pycsamt.models.modem.log.ModEmLog`. The plot
    shows normalized RMS misfit against iteration number and,
    optionally, marks the lowest RMS value. A reference line at
    :math:`RMS=1` is included because values near one generally
    indicate a data fit comparable to the assigned errors.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded ModEM inversion result. The result must contain
        a parsed ``log`` attribute.
    show_best : bool, default True
        Whether to mark the iteration with the lowest parsed
        RMS value.

    Examples
    --------
    >>> from pycsamt.models.modem.results import InversionResult
    >>> from pycsamt.models.modem.plot import PlotMisfit
    >>> result = InversionResult("modem_run")
    >>> fig = PlotMisfit(result=result).plot()
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        show_best: bool = True,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.show_best = show_best

    def plot(self):
        """Return a matplotlib figure containing the RMS curve.

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing one axes with RMS history.

        Raises
        ------
        ValueError
            If no result is attached or the result has no
            parsed log.
        """
        import matplotlib.pyplot as plt

        r = self._check_result()
        if r.log is None:
            raise ValueError("InversionResult has no log loaded.")

        iters = r.log.iterations
        rms = r.log.rms

        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(iters, rms, "k.-", lw=1.5, ms=5, label="RMS")

        if self.show_best and rms.size:
            bi = int(np.argmin(rms))
            ax.scatter(
                [iters[bi]],
                [rms[bi]],
                color="tab:red",
                zorder=5,
                label=f"best={rms[bi]:.3f}",
            )

        ax.axhline(1.0, color="gray", ls="--", lw=0.8, label="RMS=1")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("RMS misfit")
        ax.set_title(f"ModEM inversion - mode={r.mode}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        return fig


# ----------------------------------------------------------------------
# PlotModel2D
# ----------------------------------------------------------------------

class PlotModel2D(_ModEmPlotBase):
    """Plot a 2-D ModEM resistivity cross-section.

    The plot displays either the initial or final 2-D model as
    a depth section. Resistivity is shown on a logarithmic
    colour scale because magnetotelluric models commonly span
    several orders of magnitude.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded inversion result containing ``model_initial`` or
        ``model_final``.
    which : {"final", "initial"}, default "final"
        Which model to display.
    depth_max : float, optional
        Maximum depth in metres to display. If omitted, all
        layers are plotted.
    rho_min, rho_max : float, default 1.0, 1000.0
        Resistivity colour-scale limits in ohm metres.
    cmap : str, default "jet_r"
        Matplotlib colourmap name.

    Examples
    --------
    >>> from pycsamt.models.modem.plot import PlotModel2D
    >>> fig = PlotModel2D(result=result, depth_max=5000).plot()
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        which: str = "final",
        depth_max: float | None = None,
        rho_min: float = 1.0,
        rho_max: float = 1000.0,
        cmap: str = "jet_r",
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.which = which
        self.depth_max = depth_max
        self.rho_min = rho_min
        self.rho_max = rho_max
        self.cmap = cmap

    def plot(self):
        """Return a matplotlib figure containing the model section."""
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt

        r = self._check_result()
        model = (
            r.model_final if self.which == "final" else r.model_initial
        )
        if model is None:
            msg = f"No {self.which} model in InversionResult."
            raise ValueError(msg)

        x_nodes = model.x_nodes
        z_nodes = model.z_nodes
        rho = model.rho_linear

        x_km = x_nodes / 1e3
        z_km = z_nodes / 1e3

        # depth crop
        nz = model.nz
        if self.depth_max is not None:
            iz_max = int(np.searchsorted(z_nodes, self.depth_max))
            iz_max = min(iz_max, nz)
        else:
            iz_max = nz

        fig, ax = plt.subplots(figsize=(10, 4))
        norm = mcolors.LogNorm(vmin=self.rho_min, vmax=self.rho_max)
        pm = ax.pcolormesh(
            x_km,
            z_km[:iz_max + 1],
            rho[:iz_max, :],
            norm=norm,
            cmap=self.cmap,
            shading="flat",
        )
        fig.colorbar(pm, ax=ax, label="Resistivity (ohm m)", pad=0.02)
        ax.set_xlabel("Offset (km)")
        ax.set_ylabel("Depth (km)")
        ax.set_title(f"ModEM 2D model - {self.which}")
        ax.invert_yaxis()
        fig.tight_layout()
        return fig


# ----------------------------------------------------------------------
# PlotModel3D
# ----------------------------------------------------------------------

class PlotModel3D(_ModEmPlotBase):
    """Plot horizontal slices through a 3-D ModEM model.

    ``PlotModel3D`` extracts model layers nearest to requested
    depths and renders each as an x-y resistivity map. The
    selected depths are interpreted in metres below the model
    top, including the depth coordinate convention stored in
    the model object.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded inversion result containing a 3-D model.
    depths : sequence of float, optional
        Depths in metres at which to extract slices. If
        omitted, the first four active earth-layer centres are
        used.
    which : {"final", "initial"}, default "final"
        Which model to display.
    rho_min, rho_max : float, default 1.0, 1000.0
        Resistivity colour-scale limits in ohm metres.
    cmap : str, default "jet_r"
        Matplotlib colourmap name.
    n_cols : int, default 2
        Number of columns in the subplot grid.
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        depths: Sequence[float] | None = None,
        which: str = "final",
        rho_min: float = 1.0,
        rho_max: float = 1000.0,
        cmap: str = "jet_r",
        n_cols: int = 2,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.depths = depths
        self.which = which
        self.rho_min = rho_min
        self.rho_max = rho_max
        self.cmap = cmap
        self.n_cols = n_cols

    def plot(self):
        """Return a matplotlib figure containing model slices."""
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt

        r = self._check_result()
        model = (
            r.model_final if self.which == "final" else r.model_initial
        )
        if model is None:
            msg = f"No {self.which} model in InversionResult."
            raise ValueError(msg)

        z_nodes = model.z_nodes
        z_centres = 0.5 * (z_nodes[:-1] + z_nodes[1:])
        n_air = model.n_air

        if self.depths is None:
            active = z_centres[n_air:n_air + 4]
        else:
            active = np.asarray(self.depths, dtype=float)

        # find layer index for each requested depth
        iz_list = [
            int(np.argmin(np.abs(z_centres - d))) for d in active
        ]
        n_slices = len(iz_list)
        n_cols = min(self.n_cols, n_slices)
        n_rows = int(np.ceil(n_slices / n_cols))

        x_km = model.x_nodes / 1e3
        y_km = model.y_nodes / 1e3
        norm = mcolors.LogNorm(vmin=self.rho_min, vmax=self.rho_max)

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(5 * n_cols, 4 * n_rows),
            squeeze=False,
        )

        for k, iz in enumerate(iz_list):
            row, col = divmod(k, n_cols)
            ax = axes[row][col]
            rho_slice = model.rho_linear[iz, :, :]
            pm = ax.pcolormesh(
                x_km,
                y_km,
                rho_slice,
                norm=norm,
                cmap=self.cmap,
                shading="flat",
            )
            depth_m = float(z_centres[iz])
            ax.set_title(f"z = {depth_m/1e3:.2f} km", fontsize=9)
            ax.set_xlabel("X (km)", fontsize=8)
            ax.set_ylabel("Y (km)", fontsize=8)
            ax.set_aspect("equal")
            fig.colorbar(pm, ax=ax, label="ohm m")

        # hide unused axes
        for k in range(n_slices, n_rows * n_cols):
            row, col = divmod(k, n_cols)
            axes[row][col].set_visible(False)

        fig.suptitle(f"ModEM 3D model slices - {self.which}", y=1.01)
        fig.tight_layout()
        return fig


# ----------------------------------------------------------------------
# PlotResponse
# ----------------------------------------------------------------------

class PlotResponse(_ModEmPlotBase):
    """Plot apparent-resistivity and phase curves by station.

    The current implementation plots observed data from
    ``result.data_obs``. It groups rows by station and draws
    apparent resistivity and phase for each component present
    in the data rows.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded result containing observed data.
    stations : list of str, optional
        Site names to plot. Defaults to all stations up to
        ``max_stations``.
    max_stations : int, default 16
        Maximum number of stations to show.
    n_cols : int, default 4
        Number of station columns in the subplot grid.
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        stations: Sequence[str] | None = None,
        max_stations: int = 16,
        n_cols: int = 4,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.stations = stations
        self.max_stations = max_stations
        self.n_cols = n_cols

    def plot(self):
        """Return a matplotlib figure containing response curves."""
        import matplotlib.pyplot as plt

        r = self._check_result()
        if r.data_obs is None:
            raise ValueError("InversionResult has no data_obs loaded.")

        data = r.data_obs
        names = self.stations or data.site_names
        names = list(names)[: self.max_stations]
        n_sites = len(names)
        n_cols = min(self.n_cols, n_sites)
        n_rows = int(np.ceil(n_sites / n_cols)) * 2

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(4 * n_cols, 3 * n_rows),
            squeeze=False,
        )

        mu0 = 4 * np.pi * 1e-7

        for idx, name in enumerate(names):
            pg_row, col = divmod(idx, n_cols)
            ax_rho = axes[2 * pg_row][col]
            ax_phs = axes[2 * pg_row + 1][col]

            si = (
                data.site_names.index(name)
                if name in data.site_names
                else None
            )
            if si is None:
                ax_rho.set_visible(False)
                ax_phs.set_visible(False)
                continue

            for blk in data.blocks:
                rows = [row for row in blk["rows"] if row[1] == si]
                if not rows:
                    continue
                p_vals = np.array([row[0] for row in rows])
                real = np.array([row[6] for row in rows])
                imag = np.array([row[7] for row in rows])
                z_abs2 = real ** 2 + imag ** 2
                rho_a = z_abs2 / (mu0 * (2 * np.pi / p_vals))
                phase = np.degrees(np.arctan2(imag, real))
                comp = rows[0][5]
                ax_rho.loglog(p_vals, rho_a, ".", ms=4, label=comp)
                ax_phs.semilogx(p_vals, phase, ".", ms=4)

            ax_rho.set_title(name, fontsize=8)
            ax_rho.set_ylabel("rho_a (ohm m)", fontsize=7)
            ax_phs.set_xlabel("Period (s)", fontsize=7)
            ax_phs.set_ylabel("Phase (deg)", fontsize=7)
            ax_phs.set_ylim(0, 90)
            ax_rho.legend(fontsize=6, loc="best")

        for idx in range(n_sites, (n_rows // 2) * n_cols):
            pg_row, col = divmod(idx, n_cols)
            axes[2 * pg_row][col].set_visible(False)
            axes[2 * pg_row + 1][col].set_visible(False)

        fig.suptitle("ModEM observed data", fontsize=10)
        fig.tight_layout()
        return fig


# ----------------------------------------------------------------------
# PlotPseudo
# ----------------------------------------------------------------------

class PlotPseudo(_ModEmPlotBase):
    """Plot apparent-resistivity and phase pseudo-sections.

    ``PlotPseudo`` selects one component from observed data and
    arranges apparent resistivity and phase on station offset
    versus period grids. This is a quick survey-scale view of
    lateral and period-dependent response variation.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded result containing observed data.
    component : str
        Data component to display (e.g. ``'TE'``, ``'ZXY'``).
    rho_min, rho_max : float, default 1.0, 1000.0
        Apparent-resistivity colour-scale limits in ohm
        metres.
    cmap : str, default "jet_r"
        Matplotlib colourmap used for apparent resistivity.
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        component: str = "TE",
        rho_min: float = 1.0,
        rho_max: float = 1000.0,
        cmap: str = "jet_r",
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.component = component
        self.rho_min = rho_min
        self.rho_max = rho_max
        self.cmap = cmap

    def plot(self):
        """Return a matplotlib figure containing pseudo-sections."""
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt

        r = self._check_result()
        if r.data_obs is None:
            raise ValueError("InversionResult has no data_obs loaded.")

        data = r.data_obs
        periods = data.periods
        site_names = data.site_names
        offsets = data.offsets / 1e3

        n_per = len(periods)
        n_site = len(site_names)
        rho_mat = np.full((n_per, n_site), np.nan)
        phs_mat = np.full((n_per, n_site), np.nan)
        mu0 = 4 * np.pi * 1e-7

        for blk in data.blocks:
            for row in blk["rows"]:
                comp = row[5]
                if comp != self.component:
                    continue
                period = row[0]
                si = row[1]
                real, imag = row[6], row[7]
                pi_idx = np.argmin(np.abs(periods - period))
                z2 = real ** 2 + imag ** 2
                rho_mat[pi_idx, si] = z2 / (mu0 * 2 * np.pi / period)
                phs_mat[pi_idx, si] = np.degrees(np.arctan2(imag, real))

        log_p = np.log10(periods)

        fig, (ax_rho, ax_phs) = plt.subplots(
            2,
            1,
            figsize=(10, 7),
            sharex=True,
        )

        norm_rho = mcolors.LogNorm(vmin=self.rho_min, vmax=self.rho_max)
        pm_r = ax_rho.pcolormesh(
            offsets,
            log_p,
            rho_mat,
            norm=norm_rho,
            cmap=self.cmap,
            shading="nearest",
        )
        fig.colorbar(pm_r, ax=ax_rho, label="rho_a (ohm m)")
        ax_rho.set_ylabel("log10 Period (s)")
        ax_rho.set_title(
            f"ModEM pseudo-section - {self.component} rho_a"
        )

        pm_p = ax_phs.pcolormesh(
            offsets,
            log_p,
            phs_mat,
            vmin=0,
            vmax=90,
            cmap="RdBu_r",
            shading="nearest",
        )
        fig.colorbar(pm_p, ax=ax_phs, label="Phase (deg)")
        ax_phs.set_ylabel("log10 Period (s)")
        ax_phs.set_xlabel("Offset (km)")
        ax_phs.set_title(
            f"ModEM pseudo-section - {self.component} phase"
        )

        fig.tight_layout()
        return fig
