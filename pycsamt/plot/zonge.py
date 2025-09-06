# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

"""
This module provides the Plot class for visualizing Zonge AVG
data. It offers high-level methods for creating common geophysical
plots like sounding curves and pseudosections, mirroring the
capabilities of Zonge's plotting utilities.
"""
from __future__ import annotations

from typing import ( 
    TYPE_CHECKING, 
    Optional, 
    Union, 
    Tuple, 
    List, 
    Literal
)
import warnings 

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (
    MaxNLocator, NullFormatter, NullLocator
)

from scipy.interpolate import UnivariateSpline
import numpy as np
import pandas as pd 

from ..constants import PI  
from ..exceptions import ProcessingError
from ..utils.validation import has_read
from ..utils.generic import get_valid_kwargs
from ..utils.plot import set_axis_grid 
from ..zonge.config import Zonge
from ..zonge.proc_utils import ( 
    get_skew, 
    get_strike, 
    prepare_strike_frame
)

from .base import BasePlot

if TYPE_CHECKING:
    from ..zonge.avg import BaseAVG, AMTAVG
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes

__all__ = ["AVGPlot"]


class AVGPlot(Zonge, BasePlot):
    r"""A class for visualizing Zonge AVG data.

    This class provides a high-level, object-oriented interface for
    creating standard geophysical plots from a loaded
    :class:`~.avg.AVG` or :class:`~.avg.AMTAVG` object. It serves
    as the primary entry point for all visualization tasks.
    
    Parameters
    ----------
    avg_data : :class:`~.avg.BaseAVG`, optional
        An initialized and loaded AVG object. If provided, the
        plotter is immediately ready for use, otherwise the
        :meth:`read` method must be called.
    **kws
        Additional keyword arguments are passed to the
        :class:`~pycsamt.plot.base.BasePlot` parent class to control
        the styling of the figures and axes (e.g., `figsize`,
        `cmap`, `font_size`).
    
    Attributes
    ----------
    avg : :class:`~.avg.BaseAVG` or None
        The AVG data object that the plotter is operating on.
    
    Methods
    -------
    read(avg_data)
        Loads an AVG data object into the plotter.
    plot_sounding_curves(...)
        Plots apparent resistivity and phase for all components.
    plot_tensor_soundings(...)
        Creates a detailed, multi-panel plot for specific stations.
    plot_pseudosection(...)
        Generates a single, detailed pseudosection plot.
    plot_pseudosections(...)
        Creates a figure with multiple pseudosection subplots.
    plot_location_map(...)
        Visualizes station coordinates and elevation.
    plot_strike(...)
        Creates a rose diagram of geoelectric strike angles.
    
    Notes
    -----
    The `Plot` class uses a composition approach, holding an
    instance of a loaded `AVG` object. It inherits from two base
    classes:
    - :class:`~.config.Zonge`: For common attributes like `verbose`
      and a standardized logger.
    - :class:`~.plot.base.BasePlot`: For a rich set of styling
      attributes that control the appearance of the plots.
    
    This design separates the data representation (`AVG`) from the
    visualization logic (`Plot`), creating a clean and maintainable
    architecture.
    
    Examples
    --------
    The typical workflow is to first load data using `AMTAVG` and
    then pass the loaded object to the `Plot` class.
    
    >>> from pycsamt.zonge import AMTAVG
    >>> from pycsamt.plot import AVGPlot
    >>> # 1. Load the data
    >>> avg = AMTAVG.from_file('data/avg/K2.avg')
    >>>
    >>> # 2. Create a plotter and generate a plot
    >>> plotter = AVGPlot(avg, figsize=(8, 6), cmap='viridis')
    >>> fig, ax = plotter.plot_pseudosection(var='rho')
    >>> fig.savefig('rho_pseudosection.png')
    
    See Also
    --------
    pycsamt.plot.base.BasePlot : The base class for plot styling.
    pycsamt.zonge.avg.AMTAVG : The primary data container class.
    """
    def __init__(
        self,
        avg_data: Optional["BaseAVG"] = None,
        *,
        verbose: bool = False,
        **kws
    ):
        Zonge.__init__(self, verbose = verbose,  **kws)
        BasePlot.__init__(self, **kws)
        
        self.avg: Union["BaseAVG", "AMTAVG", None] = avg_data
        if avg_data is not None:
            self.read(avg_data)

    def read(self, avg_data: "BaseAVG") -> "AVGPlot":
        r"""
        Load an AVG data object into the plotter.

        This is the primary method for associating a dataset with
        the Plot class.

        Parameters
        ----------
        avg_data : :class:`~.avg.BaseAVG`
            An `AVG` or `AMTAVG` object that has already been
            loaded with data.

        Returns
        -------
        self : Plot
            The method returns the instance of the class, allowing
            for convenient method chaining.
        """
        has_read(avg_data)
        self.avg = avg_data
        if self.verbose:
            src_name = (
                self.avg._source_path.name
                if self.avg._source_path
                else 'DataFrame'
            )
            self._logger.info(
                f"Plotter initialized with data from '{src_name}'."
            )
        return self
    
    def plot_sounding_curves(
        self,
        station_id: Optional[
            Union[int, float, List[Union[int, float]]]] = None,
        *,
        x_axis: str = 'frequency',
        phase_in_degrees: bool = False,
        corrected_df: Optional[pd.DataFrame] = None,
        **kwargs
    ) -> Tuple["Figure", np.ndarray]:
        """
        Plot apparent resistivity and phase sounding curves.

        This method generates a two-panel figure with a log-log plot
        of apparent resistivity vs. frequency and a semi-log plot
        of phase vs. frequency for each station and component.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)

        df = self.avg.df
        if station_id is not None:
            s_ids = [station_id] if isinstance(
                station_id, (int, float)) else station_id
            df = df[df['station'].isin(s_ids)]
            if df.empty:
                raise ValueError(f"No data found for station(s): {s_ids}")

        components = sorted(df['comp'].unique())
        if not components:
            warnings.warn("No components found to plot.")
            return plt.subplots(2, 1, figsize=self.fig_size)

        n_comps = len(components)
        fig, axes = plt.subplots(
            2, n_comps, figsize=(n_comps * 4, 8),
            sharex=True, squeeze=False,
            gridspec_kw={'height_ratios': [2, 1]}
        )

        off_diag_color = 'blue'
        diag_color = 'red'

        for i, comp in enumerate(components):
            ax_rho, ax_phase = axes[0, i], axes[1, i]
            data = df[df['comp'] == comp].sort_values('freq')
            color = off_diag_color if comp in ['ExHy', 'EyHx'] else diag_color

            x_values = 1. / data['freq'] if x_axis == 'period' else data['freq']
            phase_values = (
                data['phase'] * (180 / (PI * 1000)) if phase_in_degrees
                else data['phase']
            )

            ax_rho.loglog(
                x_values, data['rho'], 'o', ms=self.ms,
                color=color, label=f'Original {comp}'
            )
            ax_phase.semilogx(
                x_values, phase_values, 'o', ms=self.ms, color=color
            )

            if corrected_df is not None:
                corr = corrected_df[
                    (corrected_df['station'].isin(df.station.unique())) &
                    (corrected_df['comp'] == comp)
                ].sort_values('freq')
                if not corr.empty:
                    corr_x = 1. / corr['freq'] if x_axis == 'period' else corr['freq']
                    corr_phase = (
                        corr['phase'] * (180 / (PI * 1000))
                        if phase_in_degrees else corr['phase']
                    )
                    ax_rho.loglog(
                        corr_x, corr['rho'], '-', color=color, lw=self.lw,
                        label=f'Corrected {comp}'
                    )
                    ax_phase.semilogx(
                        corr_x, corr_phase, '-', color=color, lw=self.lw
                    )

            ax_rho.set_title(comp)
            ax_rho.grid(True, which='both', ls=':', alpha=0.6)
            ax_phase.grid(True, which='both', ls=':', alpha=0.6)

        # Common styling
        x_label = "Period (s)" if x_axis == 'period' else "Frequency (Hz)"
        phase_label = "Phase (deg)" if phase_in_degrees else "Phase (mrad)"
        for ax in axes[1, :]:
            ax.set_xlabel(x_label)
        if x_axis == 'frequency':
            axes[1, 0].invert_xaxis()
        axes[0, 0].set_ylabel("App. Res. ($\Omega \cdot m$)")
        axes[1, 0].set_ylabel(phase_label)
        axes[0, 0].legend(fontsize='small')

        title = self.fig_title or "Tensor Sounding Curves"
        fig.suptitle(
            title, fontsize=self.font_size + 6, fontweight='bold'
        )
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)
        return fig, axes

    def plot_pseudosection(
        self,
        var: str = 'rho',
        comp: str = 'ExHy',
        *,
        station_names: Optional[List[str]] = None,
        tick_step: int = 1,
        cbar_label: Optional[str] = None,
        ax: Optional["Axes"] = None,
        **kwargs
    ) -> Tuple["Figure", "Axes"]:
        """
        Plot a single, detailed pseudosection with dual x-axes.
        """
        has_read(self, attributes="avg")
        kwargs = get_valid_kwargs(self.update, kwargs)
        self.update(**kwargs)

        df_comp = self.avg.df[self.avg.df['comp'] == comp]
        if df_comp.empty:
            warnings.warn(f"No data for component '{comp}'. Plot is empty.")
            return self._fig_ax()

        df_pivot = df_comp.pivot_table(
            index='freq', columns='station', values=var
        )
        x_coords = df_pivot.columns.values
        y_coords = df_pivot.index.values

        fig, ax = (plt.gcf(), ax) if ax else self._fig_ax()

        norm = (
            mcolors.LogNorm() if 'rho' in var.lower()
            else mcolors.Normalize()
        )

        if self.plt_style == 'imshow':
            c = ax.imshow(
                df_pivot.values, aspect='auto', cmap=self.cmap,
                norm=norm, origin='lower',
                extent=[x_coords.min(), x_coords.max(),
                        y_coords.min(), y_coords.max()]
            )
        else: # Default to pcolormesh
            c = ax.pcolormesh(
                x_coords, y_coords, df_pivot.values,
                cmap=self.cmap, shading='auto', norm=norm
            )
        
        if cbar_label is None:
            if 'rho' in var:
                cbar_label = r"$\rho_a$ ($\Omega \cdot m$)"
            elif 'phase' in var:
                cbar_label = "Phase (mrad)"
            else:
                cbar_label = var.replace('_', ' ').capitalize()

        cbar = fig.colorbar(c, ax=ax, **self.cb_props)
        cbar.set_label(cbar_label)

        ax.set_yscale('log')
        ax.invert_yaxis()
        ax.set_xlabel(self.xlabel or "Station Distance (m)")
        ax.set_ylabel(self.ylabel or "Frequency (Hz)")

        ax.set_xticks(x_coords[::tick_step])
        ax.tick_params(axis='x')

        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        ax_top.set_xticks(x_coords[::tick_step])
        s_names = station_names or [
            f"S{i+1}" for i in range(len(x_coords))
        ]
        ax_top.set_xticklabels(s_names[::tick_step])
        ax_top.tick_params(axis='x', labelsize='small', rotation=45)
        ax_top.set_xlabel("Station")

        set_axis_grid(
            ax, show_grid=self.show_grid, grid_props=self.grid_props
        )

        ax.set_title(
            self.fig_title or f"{comp} {var.capitalize()}",
            fontsize=self.font_size + 4, fontweight='bold'
        )
        return fig, ax

    def plot_pseudosections(
        self,
        vars: Union[str, List[str]] = 'rho',
        comp: str = 'ExHy',
        station_names: Optional[List[str]] = None,
        tick_step: int = 1,
        cbar_label: Optional[str] = None,
        max_cols: int = 2,
        **kwargs
    ) -> Tuple["Figure", np.ndarray]:
        """
        Plot one or more pseudosections in a flexible grid.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)

        if isinstance(vars, str):
            vars = [vars]

        for var in vars:
            if var not in self.avg.df.columns:
                raise ProcessingError(
                    f"Variable '{var}' not found in the DataFrame."
                )

        n_vars = len(vars)
        ncols = min(n_vars, max_cols)
        nrows = int(np.ceil(n_vars / ncols))

        fig, axes = plt.subplots(
            nrows, ncols, figsize=(ncols * 7, nrows * 5),
            sharex=True, sharey=True, squeeze=False
        )
        axes = axes.flatten()

        for i, var in enumerate(vars):
            self.plot_pseudosection(
                var=var, comp=comp, ax=axes[i],
                station_names=station_names, 
                cbar_label=cbar_label, 
                tick_step = tick_step, 
            )

        # Hide unused subplots
        for i in range(n_vars, len(axes)):
            axes[i].set_visible(False)

        fig.suptitle(
            self.fig_title or f"Pseudosections for {comp}",
            fontsize=self.font_size + 6, fontweight='bold'
        )
        fig.tight_layout(rect=[0, 0.03, 1, 0.96])

        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)

        return fig, axes
    

    def plot_tensor_soundings(
        self,
        station_id: Optional[
            Union[int, float, List[Union[int, float]]]
        ] = None,
        *,
        tensor: str = "rho",
        x_axis: str = "frequency",
        todeg: bool = False,
        corrected_df: Optional[pd.DataFrame] = None,
        show_fit: bool = False,
        **kwargs,
    ) -> Tuple["Figure", np.ndarray]:

        # ensure data are loaded and sync config with kwargs
        has_read(self, attributes="avg")
        self.update(**kwargs)
    
        # filter stations if a subset is requested
        df = self.avg.df
        s_ids = df["station"].unique()
        if station_id is not None:
            s_ids = (
                [station_id]
                if isinstance(station_id, (int, float))
                else station_id
            )
            df = df[df["station"].isin(s_ids)]
            if df.empty:
                raise ValueError(
                    f"No data found for station(s): {s_ids}"
                )
    
        # collect tensor components present in the data
        comps = sorted(df["comp"].unique())
        if not comps:
            warnings.warn("No components found to plot.")
            return plt.subplots(2, 1, figsize=self.fig_size)
    
        # make a 2 x N grid of subplots
        n = len(comps)
        fig, axes = plt.subplots(
            2,
            n,
            figsize=(n * 4.5, 8.0),
            sharex=True,
            squeeze=False,
            gridspec_kw={"height_ratios": [2, 1]},
        )
    
        # color scheme for off-diagonal vs diagonal
        off_diag = "blue"
        diag = "red"
    
        # iterate over components and plot each station
        for i, comp in enumerate(comps):
            ax1, ax2 = axes[0, i], axes[1, i]
    
            for stn in s_ids:
                sel = (
                    (df["comp"] == comp)
                    & (df["station"] == stn)
                )
                data = df[sel].sort_values("freq")
                if data.empty:
                    continue
    
                color = (
                    off_diag
                    if comp in ["ExHy", "EyHx"]
                    else diag
                )
                label = f"Stn {stn}"
    
                # choose x-domain
                x_vals = (
                    1.0 / data["freq"]
                    if x_axis == "period"
                    else data["freq"]
                )
    
                # choose amplitude source (rho or |Z|)
                if tensor == "z":
                    y1_vals = np.abs(self.avg.z.z[data.index])
                else:
                    y1_vals = data["rho"]
    
                # convert phase if requested
                if todeg:
                    # phase is in mrad → convert to degrees
                    y2_vals = (
                        data["phase"] * 180.0 / (np.pi * 1000.0)
                    )
                else:
                    y2_vals = data["phase"]
    
                # scatter original points
                ax1.loglog(
                    x_vals,
                    y1_vals,
                    "o",
                    ms=self.ms,
                    color=color,
                    label=label,
                )
                ax2.semilogx(
                    x_vals,
                    y2_vals,
                    "o",
                    ms=self.ms,
                    color=color,
                )
    
                # optional smooth fit for trend visualization
                if show_fit and len(data) >= 5:
                    # fit on log-freq for smoother trends
                    lf = np.log(data["freq"])
    
                    s1 = UnivariateSpline(
                        lf,
                        np.log(y1_vals),
                        s=0.1,
                    )
                    s2 = UnivariateSpline(
                        lf,
                        y2_vals,
                        s=0.1,
                    )
    
                    fx = np.logspace(
                        np.log10(x_vals.min()),
                        np.log10(x_vals.max()),
                        100,
                    )
                    # map fx to frequency depending on x_axis
                    if x_axis == "period":
                        f_fit = 1.0 / fx
                    else:
                        f_fit = fx
    
                    ax1.loglog(
                        fx,
                        np.exp(s1(np.log(f_fit))),
                        "--",
                        color=color,
                        lw=self.lw,
                    )
                    ax2.semilogx(
                        fx,
                        s2(np.log(f_fit)),
                        "--",
                        color=color,
                        lw=self.lw,
                    )
    
            # panel styling for this component
            ax1.set_title(comp)
            ax1.grid(True, which="both", ls=":", alpha=0.6)
            ax2.grid(True, which="both", ls=":", alpha=0.6)
    
        # common labels and legend
        x_lab = (
            "Period (s)"
            if x_axis == "period"
            else "Frequency (Hz)"
        )
        y1_lab = (
            r"$|Z|$ ($\Omega$)"
            if tensor == "z"
            else r"App. Res. ($\Omega \cdot m$)"
        )
        y2_lab = "Phase (deg)" if todeg else "Phase (mrad)"
    
        for ax in axes[1, :]:
            ax.set_xlabel(x_lab)
    
        if x_axis == "frequency":
            axes[1, 0].invert_xaxis()
    
        axes[0, 0].set_ylabel(y1_lab)
        axes[1, 0].set_ylabel(y2_lab)
        axes[0, 0].legend(fontsize="small")
    
        # figure title and layout
        title = self.fig_title or "Tensor Sounding Curves"
        fig.suptitle(
            title,
            fontsize=self.font_size + 6,
            fontweight="bold",
        )
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    
        # optional saving
        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)
    
        return fig, axes
 
    def plot_phase_tensor_ellipses(
        self,
        station_id: Optional[Union[int, float]] = None,
        *,
        frequencies: Optional[List[float]] = None,
        **kwargs
    ) -> Tuple["Figure", np.ndarray]:
        """
        Visualize the phase tensor as ellipses for a station.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)

        if station_id is None:
            station_id = self.avg.df['station'].unique()[0]

        df_station = self.avg.df[
            self.avg.df['station'] == station_id
        ].copy() # Use a copy to add the 'z' column
        if df_station.empty:
            raise ValueError(
                f"Station ID '{station_id}' not found in data."
            )
        
        # --- FIX: Calculate and add the complex 'z' column ---
        df_station['z'] = self.avg.z.z.loc[df_station.index]
        # -----------------------------------------------------

        if frequencies is None:
            n_freqs = min(9, df_station['freq'].nunique())
            freqs = np.logspace(
                np.log10(df_station['freq'].min()),
                np.log10(df_station['freq'].max()),
                n_freqs
            ).round(2)
        else:
            freqs = frequencies

        n_freqs = len(freqs)
        ncols = min(n_freqs, 4)
        nrows = int(np.ceil(n_freqs / ncols))
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(ncols * 3, nrows * 3),
            squeeze=False
        )
        axes = axes.flatten()

        for i, freq in enumerate(freqs):
            ax = axes[i]
            freq_data = df_station.iloc[
                (df_station['freq'] - freq).abs().argsort()[:4]
            ]
            if freq_data.empty: continue

            Z = np.zeros((2, 2), dtype=complex)
            for _, row in freq_data.iterrows():
                if row['comp'] in ('ExHx', 'Zxx'): Z[0, 0] = row['z']
                if row['comp'] in ('ExHy', 'Zxy'): Z[0, 1] = row['z']
                if row['comp'] in ('EyHx', 'Zyx'): Z[1, 0] = row['z']
                if row['comp'] in ('EyHy', 'Zyy'): Z[1, 1] = row['z']

            if np.linalg.det(Z.real) == 0:
                if self.verbose:
                    self._logger.warning(f"Skipping freq {freq} due to "
                                         "singular real impedance tensor.")
                continue

            M = np.linalg.inv(Z.real)
            Phi = np.dot(M, Z.imag)

            eigvals, eigvecs = np.linalg.eigh(Phi)
            phi_maj, phi_min = eigvals
            angle = np.rad2deg(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
            skew = np.abs(Z[0, 0] + Z[1, 1]) / np.abs(Z[0, 1] - Z[1, 0])

            ellipse = Ellipse(
                (0, 0), width=phi_maj*2, height=phi_min*2,
                angle=angle,
                facecolor=plt.cm.viridis(skew / 0.5),
                alpha=0.7
            )
            ax.add_patch(ellipse)

            ax.set_title(f"{freq:.1f} Hz", fontsize='small')
            ax.set_aspect('equal', 'box')
            lim = np.max(np.abs(eigvals)) * 1.2
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.grid(True, ls=':')

        for i in range(n_freqs, len(axes)):
            axes[i].set_visible(False)

        fig.suptitle(
            f"Phase Tensor Ellipses for Station {station_id}",
            fontsize=self.font_size + 6, fontweight='bold'
        )
        fig.tight_layout(rect=[0, 0, 1, 0.95])

        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)

        return fig, axes

    def plot_remediation(
        self,
        station_id: Union[int, float],
        corrected_df: pd.DataFrame,
        **kwargs
    ) -> Tuple["Figure", "Axes"]:
        """
        Visualize the effect of a remedial action on a sounding.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)

        original_df = self.avg.df[
            self.avg.df['station'] == station_id
        ]
        if original_df.empty:
            raise ValueError(
                f"Station ID '{station_id}' not found in original data."
            )

        fig, axes = plt.subplots(
            2, 1, figsize=self.fig_size, sharex=True)
        ax_rho, ax_phase = axes

        for comp, group in original_df.groupby('comp'):
            # Plot original data with error bars
            ax_rho.errorbar(
                group['freq'], group['rho'],
                yerr=(group['rho'] * group.get('pc_rho', 0) / 100),
                fmt='o', ms=self.ms, label=f'Original {comp}',
                capsize=3
            )
            ax_phase.errorbar(
                group['freq'], group['phase'],
                yerr=group.get('s_phz', 0),
                fmt='o', ms=self.ms, capsize=3
            )

            # Plot corrected data as a line
            corr_group = corrected_df[
                (corrected_df['station'] == station_id) &
                (corrected_df['comp'] == comp)
            ]
            if not corr_group.empty:
                ax_rho.plot(
                    corr_group['freq'], corr_group['rho'],
                    ls='-', lw=self.lw,
                    label=f'Corrected {comp}'
                )
                ax_phase.plot(
                    corr_group['freq'], corr_group['phase'],
                    ls='-', lw=self.lw
                )

        # Styling
        ax_rho.set_yscale('log')
        ax_rho.set_ylabel("App. Res. ($\Omega \cdot m$)")
        ax_rho.legend(fontsize='small')
        ax_rho.grid(True, which='both', ls=':', alpha=0.6)

        ax_phase.set_xscale('log')
        ax_phase.set_xlabel("Frequency (Hz)")
        ax_phase.set_ylabel("Phase (mrad)")
        ax_phase.invert_xaxis()
        ax_phase.grid(True, which='both', ls=':', alpha=0.6)

        fig.suptitle(
            f"Remediation Plot for Station {station_id}",
            fontsize=self.font_size + 6, fontweight='bold'
        )
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)

        return fig, axes

    def plot_station(
        self,
        station_id: Union[int, float, List[Union[int, float]]],
        *,
        x_axis: str = "frequency",
        todeg: bool = False,
        corrected_df: Optional[pd.DataFrame] = None,
        **kwargs,
    ) -> "Figure":
        """
        Create a comprehensive multi-panel plot for stations.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)
    
        station_ids = (
            [station_id]
            if isinstance(station_id, (int, float))
            else station_id
        )
        df_station = self.avg.df[
            self.avg.df["station"].isin(station_ids)
        ]
        if df_station.empty:
            raise ValueError(
                f"No data found for station(s): {station_ids}"
            )
    
        fig = plt.figure(figsize=self.fig_size)
        gs = gridspec.GridSpec(3, 3, figure=fig)
    
        ax_rho = fig.add_subplot(gs[0, :2])
        ax_emag = fig.add_subplot(gs[0, 2], sharey=ax_rho)
        ax_phase = fig.add_subplot(gs[1, :2], sharex=ax_rho)
        ax_hmag = fig.add_subplot(gs[1, 2], sharey=ax_phase)
        ax_skew = fig.add_subplot(gs[2, :2], sharex=ax_rho)
    
        for stn_id in station_ids:
            stn_data = df_station[df_station.station == stn_id]
            if todeg: 
                stn_data['phase'] = stn_data['phase'] * (
                    180 / (np.pi * 1000)) 
                
            self._plot_single_station(
                stn_data,
                ax_rho,
                ax_phase,
                ax_emag,
                ax_hmag,
                corrected_df,
                x_axis
            )
    
        # Add the 'z' column before calling get_skew 
        df_for_skew = df_station[
            df_station.station == station_ids[0]
        ].copy()
        df_for_skew['z'] = self.avg.z.z.loc[df_for_skew.index]
        
        skew_data = get_skew(df_for_skew)
        
        if not skew_data.empty:
            ax_skew.semilogx(
                skew_data["freq"],
                skew_data["skew"],
                marker="d",
                color="purple",
            )
    
        self._style_station_plot(
            fig,
            ax_rho,
            ax_phase,
            ax_emag,
            ax_hmag,
            ax_skew,
            x_axis,
            todeg,
        )
        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)
   
        return fig
    
    
    def _plot_single_station(
        self,
        df: pd.DataFrame,
        ax_rho,
        ax_phase,
        ax_emag,
        ax_hmag,
        corrected_df: Optional[pd.DataFrame],
        x_axis
    ) -> None:
        """Helper to plot data for one station on the axes."""
        station_id = df["station"].iloc[0]
        components = sorted(df["comp"].unique())
    
        for comp in components:
            data = df[df["comp"] == comp].sort_values("freq")
            if data.empty:
                continue
    
            x_values = (
                1.0 / data["freq"]
                if x_axis == "period"
                else data["freq"]
            )
            color = (
                "blue"
                if comp in ["ExHy", "EyHx"]
                else "red"
            )
            label = f"Stn {station_id} ({comp})"
    
            ax_rho.loglog(
                x_values,
                data["rho"],
                "o",
                color=color,
                label=label,
            )
            ax_phase.semilogx(
                x_values,
                data["phase"],
                "o",
                color=color,
            )
            ax_emag.loglog(
                x_values,
                data["emag"],
                "s",
                color=color,
            )
            ax_hmag.loglog(
                x_values,
                data["hmag"],
                "^",
                color=color,
            )
    
            if corrected_df is not None:
                corr = corrected_df[
                    (corrected_df["station"] == station_id)
                    & (corrected_df["comp"] == comp)
                ].sort_values("freq")
                if not corr.empty:
                    corr_x = (
                        1.0 / corr["freq"]
                        if x_axis == "period"
                        else corr["freq"]
                    )
                    ax_rho.loglog(
                        corr_x,
                        corr["rho"],
                        "-",
                        color=color,
                    )
                    ax_phase.semilogx(
                        corr_x,
                        corr["phase"],
                        "-",
                        color=color,
                    )
    
    
    def _style_station_plot(
        self,
        fig,
        ax_rho,
        ax_phase,
        ax_emag,
        ax_hmag,
        ax_skew,
        x_axis: str,
        todeg: bool,
    ) -> None:
        """Apply final styling to the station plot."""
        x_label = (
            "Period (s)" if x_axis == "period" else "Frequency (Hz)"
        )
        phase_label = (
            "Phase (deg)" if todeg else "Phase (mrad)"
        )
    
        ax_rho.set_ylabel(r"App. Res. ($\Omega \cdot m$)")
        ax_rho.grid(True, which="both", ls=":", alpha=0.6)
        ax_rho.legend(fontsize="small")
        plt.setp(ax_rho.get_xticklabels(), visible=False)
    
        ax_phase.set_ylabel(phase_label)
        ax_phase.grid(True, which="both", ls=":", alpha=0.6)
    
        ax_emag.set_ylabel("E-Magnitude (nV/Am)")
        ax_emag.yaxis.tick_right()
        ax_emag.yaxis.set_label_position("right")
    
        ax_hmag.set_ylabel("H-Magnitude (pT/A)")
        ax_hmag.yaxis.tick_right()
        ax_hmag.yaxis.set_label_position("right")
    
        ax_skew.set_ylabel("Skew")
        ax_skew.set_xlabel(x_label)
    
        if x_axis == "frequency":
            ax_skew.invert_xaxis()
    
        fig.tight_layout(
            rect=(0, 0.03, 1, 0.95)
        )
        
    def plot_profile(
        self,
        *,
        station_names: Optional[Union[str, List[str]]] = None,
        tick_step: int = 1,
        right_axis: str = "mask",   # {"mask","show","off"}
        top_max_ticks: int = 40,
        top_rotation: float = 0.0,
        top_label: Optional[str] = None,  # None → no label
        **kwargs,
    ) -> Tuple["Figure", "Axes"]:
        r"""
        Plot elevation vs distance with a masked top axis.
    
        - Top axis spines are invisible by default.
        - If ``station_names`` is provided (or inferred), only the
          **labels** are drawn (no line, no ticks).
        - ``station_names='mask'`` or ``None`` hides all labels.
        - ``right_axis='mask'`` fully hides the right y-axis.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)
    
        if not getattr(self.avg, "topo", None):
            raise ProcessingError(
                "Topography not found. Call "
                "`avg.add_topography(stn_file=...)` first."
            )
    
        topo = self.avg.topo
        dist = np.asarray(topo.stations, dtype=float)
        elev = np.asarray(topo.elevations, dtype=float)
    
        fig, ax = plt.subplots(1, 1, figsize=self.fig_size)
    
        ax.plot(
            dist, elev, marker="o", ls="-",
            **getattr(self, "plt_kws", {}),
        )
    
        ax.set_xlabel(self.xlabel or "Distance Along Line (m)")
        ax.set_ylabel(self.ylabel or "Elevation (m)")
        ax.set_title(self.fig_title or "Station Elevation Profile")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
        ax.grid(True, which="both", ls=":", alpha=0.7)
    
        # ---------------- right axis control ---------------- #
        if right_axis == "show":
            ax_r = ax.twinx()
            ax_r.set_ylim(ax.get_ylim())
            ax_r.set_ylabel(self.ylabel or "Elevation (m)")
        elif right_axis == "mask":
            ax_r = ax.twinx()
            # hide *everything* on the right axis
            ax_r.set_visible(False)
        # if "off": do not create a twin axis
    
        # ---------------- top axis with names --------------- #
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
    
        # hide ALL spines on the top axis
        for sp in ax_top.spines.values():
            sp.set_visible(False)
    
        # remove tick marks by default
        ax_top.tick_params(
            axis="x", which="both", length=0, width=0
        )
    
        # Decide tick stride to avoid MAXTICKS warnings
        n = len(dist)
        auto = max(1, int(np.ceil(n / float(top_max_ticks))))
        step = max(1, int(tick_step), auto)
        ticks = dist[::step]
    
        def _mask_top_axis() -> None:
            ax_top.set_xticks([])
            ax_top.xaxis.set_major_locator(NullLocator())
            ax_top.set_xlabel("" if top_label is None else top_label)
            # make sure no ticklabels are drawn
            ax_top.set_xticklabels([])
    
        # station_names logic
        if isinstance(station_names, str):
            if station_names.lower() == "mask":
                _mask_top_axis()
            else:
                names = [station_names] * len(dist)
                labs = names[::step]
                ax_top.set_xticks(ticks)
                ax_top.set_xticklabels(
                    labs, rotation=top_rotation, ha="center",
                )
                if top_label is not None:
                    ax_top.set_xlabel(top_label)
        else:
            # infer from topo or use provided list
            if station_names is None:
                # default: show labels (short form), not mask
                names = np.asarray(topo.stations).astype(str)
            else:
                names = np.asarray(station_names).astype(str)
    
            if names.shape[0] != dist.shape[0]:
                # mismatch → mask the top axis entirely
                _mask_top_axis()
            else:
                labs = names[::step]
                ax_top.set_xticks(ticks)
                ax_top.set_xticklabels(
                    labs, rotation=top_rotation, ha="center",
                )
                if top_label is not None:
                    ax_top.set_xlabel(top_label)
    
        # ensure top axis draws only labels; no ticks nor frames
        ax_top.xaxis.set_tick_params(length=0)
        ax_top.yaxis.set_major_locator(NullLocator())
        ax_top.yaxis.set_major_formatter(NullFormatter())
    
        fig.tight_layout()
    
        if getattr(self, "savefig", None):
            fig.savefig(self.savefig, dpi=self.fig_dpi)
    
        return fig, ax

    def plot_location_map(
        self,
        *,
        kind: Literal["scatter", "contour", "contourf"] = "scatter",
        n_levels: int = 10,
        cbar_label: Optional[str] = None,
        tick_step: int = 1,
        **kwargs
    ) -> Tuple["Figure", "Axes"]:
        r"""Plot a map or profile of station locations.

        This method creates a visualization of the survey geometry.
        It can generate a 2D plan-view map or a 2D elevation
        profile along the survey line.

        Parameters
        ----------
        kind : {'scatter', 'contour', 'contourf', 'profile'}, default 'scatter'
            The type of plot to generate.
            - 'scatter': A 2D plan-view map with points colored by
              elevation.
            - 'contour', 'contourf': A 2D plan-view contour map of
              the topography. Note: This requires non-collinear
              station locations.
            - 'profile': A 2D plot of elevation versus the
              cumulative distance along the survey line.
        n_levels : int, default 10
            The number of contour levels for 'contour' plots.
        cbar_label : str, optional
            A custom label for the color bar.
        tick_step : int, default 1
            The step for displaying station name labels on the top
            axis of the profile plot to prevent overlap.
            
        **kwargs
            Additional keyword arguments passed to the underlying
            matplotlib plot function.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)

        if not hasattr(self.avg, 'topo') or self.avg.topo is None:
            raise ProcessingError(
                "Topography data not found. Use the "
                "`avg.add_topography()` method to load a .stn file first."
            )

        topo = self.avg.topo
        fig, ax = self._fig_ax()
    
        # Handle map-view plots
        if kind == 'scatter':
            sc = ax.scatter(
                topo.eastings, topo.northings, c=topo.elevations,
                cmap=self.cmap, s=self.s, **self.plt_kws
            )
            cbar = fig.colorbar(sc, ax=ax, **self.cb_props)
            cbar.set_label(cbar_label or "Elevation (m)")
        else: # contour or contourf
            try:
                grid_x, grid_y, grid_z = topo.to_grid()
                levels = np.linspace(
                    np.nanmin(grid_z), np.nanmax(grid_z), n_levels
                )
                if kind == 'contourf':
                    cf = ax.contourf(
                        grid_x, grid_y, grid_z, levels=levels,
                        cmap=self.cmap, **self.plt_kws
                    )
                    cbar = fig.colorbar(cf, ax=ax, **self.cb_props)
                    cbar.set_label(cbar_label or "Elevation (m)")
                else: # contour
                    cs = ax.contour(
                        grid_x, grid_y, grid_z, levels=levels,
                        cmap=self.cmap, **self.plt_kws
                    )
                    ax.clabel(cs, inline=True, fontsize='small')
            except Exception as e:
                warnings.warn(f"Could not create contour plot: {e}")

        ax.set_aspect('equal', 'box')
        ax.set_xlabel(self.xlabel or "Easting (m)")
        ax.set_ylabel(self.ylabel or "Northing (m)")
        ax.set_title(
            self.fig_title or "Station Location Map"
        )
        # Set a reasonable number of ticks to avoid warnings
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=10))
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=10))
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")


        set_axis_grid(
            ax, show_grid=self.show_grid, grid_props=self.grid_props
        )
        fig.tight_layout()

        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)

        return fig, ax   
    

    def plot_strike(
        self,
        *,
        num_bins: int = 36,
        corrected_df: Optional[pd.DataFrame] = None,
        **kwargs,
    ) -> Tuple["Figure", "Axes"]:
        r"""
        Plot a rose diagram of geoelectric strike angles.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)
    
        # --------- build strike input (prefers complex z) -------- #
        zf = getattr(getattr(self.avg, "z", None), "frame", None)
        try:
            strike_in = prepare_strike_frame(
                z_frame=zf,
                df=self.avg.df,
                prefer="z",
                phase_unit="auto",
                drop_na=True,
                na_policy="any",
            )
        except ProcessingError as exc:
            raise ProcessingError(
                f"Unable to prepare strike input: {exc}"
            ) from exc
    
        strike_df = get_strike(strike_in)
        if strike_df.empty:
            warnings.warn(
                "No strike angles could be calculated.",
                stacklevel=2,
            )
            return plt.subplots(
                figsize=self.fig_size,
                subplot_kw={"projection": "polar"},
            )
    
        # ------------- histogram on polar axes (rose) ------------ #
        ang = strike_df["strike_angle"].dropna().to_numpy(float)
        ang = np.concatenate([ang, ang + 180.0])
        ang_rad = np.deg2rad(ang)
    
        bins = np.linspace(0.0, 2.0 * np.pi, num_bins + 1)
        counts, bin_edges = np.histogram(ang_rad, bins=bins)
        widths = np.diff(bin_edges)
        centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    
        fig, ax = plt.subplots(
            figsize=self.fig_size,
            subplot_kw={"projection": "polar"},
        )
    
        bars = ax.bar(
            bin_edges[:-1],
            counts,
            width=widths,
            edgecolor="k",
            lw=0.5,
            alpha=0.7,
            label="Original",
        )
        norm = mcolors.Normalize(vmin=0, vmax=max(1, counts.max()))
        cmap = getattr(self, "cmap", "viridis")
        cm = plt.cm.get_cmap(cmap)
        for b, c in zip(bars, counts):
            b.set_facecolor(cm(norm(c)))
    
        # ------------------- optional corrected ------------------ #
        if corrected_df is not None:
            try:
                corr_in = prepare_strike_frame(
                    z_frame=None,
                    df=corrected_df,
                    prefer="df",
                    phase_unit="auto",
                    drop_na=True,
                    na_policy="any",
                )
                strike_corr = get_strike(corr_in)
            except Exception as exc:
                warnings.warn(
                    f"Corrected strike skipped: {exc}",
                    stacklevel=2,
                )
                strike_corr = pd.DataFrame()
    
            if not strike_corr.empty:
                a2 = strike_corr["strike_angle"].dropna()
                a2 = np.concatenate([a2, a2 + 180.0])
                a2r = np.deg2rad(a2)
                cnt2, _ = np.histogram(a2r, bins=bins)
                ax.plot(
                    centers,
                    cnt2,
                    color="red",
                    lw=1.5,
                    label="Corrected",
                    drawstyle="steps-mid",
                )
                ax.legend()
    
        # -------------------------- style ------------------------ #
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_title(
            self.fig_title or "Geoelectric Strike Distribution",
            fontsize=self.font_size + 4,
            fontweight="bold",
            pad=20,
        )
    
        if getattr(self, "savefig", None):
            fig.savefig(self.savefig, dpi=self.fig_dpi)
    
        return fig, ax


