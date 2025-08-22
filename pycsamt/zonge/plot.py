# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
pycsamt.zonge.plot
------------------

This module provides the Plot class for visualizing Zonge AVG
data. It offers high-level methods for creating common geophysical
plots like sounding curves and pseudosections, mirroring the
capabilities of Zonge's plotting utilities.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union, Tuple, List
import warnings 

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

import numpy as np
import pandas as pd 

from ..exceptions import ProcessingError
from ..utils.validation import has_read
from ..plot.base import BasePlot
 
from .config import Zonge
from ..constants import PI 
from .proc_utils import get_skew, get_strike

if TYPE_CHECKING:
    from .avg import BaseAVG, AMTAVG
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes

__all__ = ["Plot"]


class Plot(Zonge, BasePlot):
    r"""A class for visualizing Zonge AVG data.

    This class provides a suite of methods for creating standard
    geophysical plots from a loaded :class:`~.avg.AVG` or
    :class:`~.avg.AMTAVG` object.

    Parameters
    ----------
    avg_data : :class:`~.avg.BaseAVG`, optional
        An initialized and loaded AVG object. If provided, the
        plotter is immediately ready for use.
    verbose : bool, default False
        If ``True``, log messages will be printed to the console
        during plotting operations.

    Examples
    --------
    >>> from pycsamt.zonge import AMTAVG, Plot
    >>> # 1. Load the data
    >>> avg = AMTAVG.from_file('data/avg/K2.avg')
    >>>
    >>> # 2. Create a plotter and generate a plot
    >>> plotter = Plot(avg)
    >>> fig, ax = plotter.plot_sounding_curves()
    >>> fig.show()
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

    def read(self, avg_data: "BaseAVG") -> "Plot":
        """
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
        station_id: Optional[Union[int, float, List[Union[int, float]]]] = None,
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
        **kwargs
    ) -> Tuple["Figure", "Axes"]:
        """
        Plot a pseudosection of a given data variable.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)

        if var not in self.avg.df.columns:
            raise ProcessingError(
                f"Variable '{var}' not found in the DataFrame."
            )

        # Filter data for the specified component
        df_comp = self.avg.df[self.avg.df['comp'] == comp]
        if df_comp.empty:
            warnings.warn(
                f"No data found for component '{comp}'. "
                "Plot will be empty."
            )
            df_pivot = pd.DataFrame()
        else:
            df_pivot = df_comp.pivot_table(
                index='freq', columns='station', values=var
            )

        fig, ax = self._fig_ax()
        ax.grid (False) 
        
        if not df_pivot.empty:
            # Create the color mesh plot
            c = ax.pcolormesh(
                df_pivot.columns,
                df_pivot.index,
                df_pivot.values,
                cmap=self.cmap,
                shading='auto'
            )
            # Add a colorbar
            cbar = fig.colorbar(c, ax=ax)
            cbar.set_label(
                f"{var.capitalize()} ({self.avg.info.meta.get(f'Unit.{var}', '')})"
            )

        # --- Styling ---
        ax.set_yscale('log')
        ax.invert_yaxis()  # Depth increases downwards

        ax.set_xlabel(self.xlabel or "Station")
        ax.set_ylabel(self.ylabel or "Frequency (Hz)")
        ax.set_title(
            self.fig_title or f"{comp} {var.capitalize()} Pseudosection",
            fontsize=self.font_size + 4,
            fontweight='bold'
        )

        if self.show_grid:
            ax.grid(True, which="both", ls=":", alpha=0.5)

        fig.tight_layout()

        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)
            self._logger.info(f"Figure saved to {self.savefig}")

        return fig, ax
    
    def plot_pseudosections(
        self,
        vars: Union[str, List[str]] = 'rho',
        comp: str = 'ExHy',
        **kwargs
    ) -> Tuple["Figure", np.ndarray]:
        """
        Plot one or more pseudosections.
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

        df_comp = self.avg.df[self.avg.df['comp'] == comp]
        if df_comp.empty:
            warnings.warn(f"No data for component '{comp}'.")
            return plt.subplots(len(vars), 1, squeeze=False)

        fig, axes = plt.subplots(
            len(vars), 1, figsize=(10, len(vars) * 4),
            sharex=True, squeeze=False
        )
        axes = axes.flatten()

        for i, var in enumerate(vars):
            ax = axes[i]
            df_pivot = df_comp.pivot_table(
                index='freq', columns='station', values=var
            )

            # Use logarithmic color scale for resistivity
            norm = (
                mcolors.LogNorm() if 'rho' in var.lower()
                else mcolors.Normalize()
            )

            c = ax.pcolormesh(
                df_pivot.columns, df_pivot.index, df_pivot.values,
                cmap=self.cmap, shading='auto', norm=norm
            )
            cbar = fig.colorbar(c, ax=ax)
            cbar.set_label(var)

            ax.set_yscale('log')
            ax.set_ylabel("Frequency (Hz)")

        ax.invert_yaxis()
        ax.set_xlabel("Station")
        fig.suptitle(
            self.fig_title or f"Pseudosections for {comp}",
            fontsize=self.font_size + 6, fontweight='bold'
        )
        fig.tight_layout(rect=[0, 0, 1, 0.96])

        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)

        return fig, axes

    def plot_tensor_soundings(
        self,
        station_id: Optional[Union[int, float]] = None,
        *,
        corrected_df: Optional[pd.DataFrame] = None,
        **kwargs
    ) -> Tuple["Figure", np.ndarray]:
        """
        Plot resistivity and phase soundings for all tensor
        components of a single station.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)

        if station_id is None:
            station_id = self.avg.df['station'].unique()[0]

        df_station = self.avg.df[
            self.avg.df['station'] == station_id
        ]
        if df_station.empty:
            raise ValueError(
                f"Station ID '{station_id}' not found in data."
            )

        components = sorted(df_station['comp'].unique())
        if not components:
            warnings.warn("No components found for the station.")
            return plt.subplots(2, 1, figsize=self.fig_size)

        n_comps = len(components)
        fig, axes = plt.subplots(
            2, n_comps, figsize=(n_comps * 4, 8),
            sharex=True, squeeze=False
        )

        off_diagonal_color = 'blue'
        diagonal_color = 'red'

        for i, comp in enumerate(components):
            ax_rho = axes[0, i]
            ax_phase = axes[1, i]

            data = df_station[df_station['comp'] == comp].sort_values(
                'freq'
            )
            color = (
                off_diagonal_color if comp in ['ExHy', 'EyHx']
                else diagonal_color
            )

            # Plot original data
            ax_rho.loglog(
                data['freq'], data['rho'], marker='o', ls='none',
                color=color, label='Original'
            )
            ax_phase.semilogx(
                data['freq'], data['phase'], marker='o', ls='none',
                color=color
            )

            # Plot corrected data if provided
            if corrected_df is not None:
                corr_data = corrected_df[
                    (corrected_df['station'] == station_id) &
                    (corrected_df['comp'] == comp)
                ].sort_values('freq')
                if not corr_data.empty:
                    ax_rho.loglog(
                        corr_data['freq'], corr_data['rho'],
                        color=color, ls='-', lw=1.5, label='Corrected'
                    )
                    ax_phase.semilogx(
                        corr_data['freq'], corr_data['phase'],
                        color=color, ls='-', lw=1.5
                    )

            ax_rho.set_title(comp)
            ax_rho.grid(True, which='both', ls=':', alpha=0.6)
            ax_phase.grid(True, which='both', ls=':', alpha=0.6)

        # Common styling
        for ax in axes[1, :]:
            ax.set_xlabel("Frequency (Hz)")
            ax.invert_xaxis()
        for ax in axes[0, :]:
            ax.set_ylabel("App. Res. ($\Omega \cdot m$)")
        for ax in axes[1, :]:
            ax.set_ylabel("Phase (mrad)")

        axes[0, 0].legend()
        fig.suptitle(
            self.fig_title or f"Tensor Soundings for Station {station_id}",
            fontsize=self.font_size + 6, fontweight='bold'
        )
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)
            self._logger.info(f"Figure saved to {self.savefig}")

        return fig, axes
    
    def plot_phase_tensor(
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
        ]
        if df_station.empty:
            raise ValueError(
                f"Station ID '{station_id}' not found in data."
            )

        if frequencies is None:
            # Auto-select a few representative frequencies
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
            # Find closest frequency in data
            freq_data = df_station.iloc[
                (df_station['freq'] - freq).abs().argsort()[:4]
            ]
            if freq_data.empty: continue

            # Reconstruct the 2x2 impedance tensor for this freq
            Z = np.zeros((2, 2), dtype=complex)
            for _, row in freq_data.iterrows():
                if row['comp'] in ('ExHx', 'Zxx'): Z[0, 0] = row['z']
                if row['comp'] in ('ExHy', 'Zxy'): Z[0, 1] = row['z']
                if row['comp'] in ('EyHx', 'Zyx'): Z[1, 0] = row['z']
                if row['comp'] in ('EyHy', 'Zyy'): Z[1, 1] = row['z']

            # Calculate the real induction tensor and phase tensor
            M = np.linalg.inv(Z.real)
            N = Z.imag
            Phi = np.dot(M, N)

            # Eigenvalue decomposition to get ellipse parameters
            eigvals, eigvecs = np.linalg.eigh(Phi)
            phi_maj, phi_min = eigvals
            angle = np.rad2deg(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))

            # Calculate skew for coloring
            skew = np.abs(Z[0, 0] + Z[1, 1]) / np.abs(Z[0, 1] - Z[1, 0])
            
            ellipse = Ellipse(
                (0, 0), width=phi_maj*2, height=phi_min*2,
                angle=angle,
                facecolor=plt.cm.viridis(skew / 0.5), # Normalize skew
                alpha=0.7
            )
            ax.add_patch(ellipse)

            ax.set_title(f"{freq:.1f} Hz", fontsize='small')
            ax.set_aspect('equal', 'box')
            lim = np.max(np.abs(eigvals)) * 1.2
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.grid(True, ls=':')

        # Hide unused subplots
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

    
    
    def plot_stations(
        self,
        station_id: Union[int, float, List[Union[int, float]]],
        *,
        x_axis: str = "frequency",
        phase_in_degrees: bool = False,
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
            self._plot_single_station(
                stn_data,
                ax_rho,
                ax_phase,
                ax_emag,
                ax_hmag,
                corrected_df,
            )
    
        skew_data = get_skew(
            df_station[df_station.station == station_ids[0]]
        )
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
            phase_in_degrees,
        )
        return fig
    
    
    def _plot_single_station(
        self,
        df: pd.DataFrame,
        ax_rho,
        ax_phase,
        ax_emag,
        ax_hmag,
        corrected_df: Optional[pd.DataFrame],
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
                if self.x_axis == "period"
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
                        if self.x_axis == "period"
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
        phase_in_degrees: bool,
    ) -> None:
        """Apply final styling to the station plot."""
        x_label = (
            "Period (s)" if x_axis == "period" else "Frequency (Hz)"
        )
        phase_label = (
            "Phase (deg)" if phase_in_degrees else "Phase (mrad)"
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
    
    def plot_location_map(
        self,
        **kwargs
    ) -> Tuple["Figure", "Axes"]:
        """
        Plot a 2D map of station locations with elevation.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)

        required_cols = ['easting', 'northing', 'elevation']
        if not all(c in self.avg.df.columns for c in required_cols):
            raise ProcessingError(
                "DataFrame must contain 'easting', 'northing', and "
                "'elevation' columns. Load a .stn file first."
            )

        # Get unique station locations
        locations = self.avg.df[
            ['station', 'easting', 'northing', 'elevation']
        ].drop_duplicates('station').sort_values('station')

        fig, ax = self._fig_ax()

        sc = ax.scatter(
            locations['easting'],
            locations['northing'],
            c=locations['elevation'],
            cmap=self.cmap,
            s=self.s
        )

        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label("Elevation (m)")

        ax.set_aspect('equal', 'box')
        ax.set_xlabel(self.xlabel or "Easting (m)")
        ax.set_ylabel(self.ylabel or "Northing (m)")
        ax.set_title(
            self.fig_title or "Station Location Map",
            fontsize=self.font_size + 4,
            fontweight='bold'
        )

        if self.show_grid:
            ax.grid(True, ls=':', alpha=0.6)

        fig.tight_layout()

        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)

        return fig, ax

    def plot_strike_rose(
        self,
        *,
        num_bins: int = 36,
        corrected_df: Optional[pd.DataFrame] = None,
        **kwargs
    ) -> Tuple["Figure", "Axes"]:
        """
        Plot a rose diagram of geoelectric strike angles.
        """
        has_read(self, attributes="avg")
        self.update(**kwargs)

        if not hasattr(self.avg, 'z'):
            raise ProcessingError(
                "Strike plot requires an AMTAVG object with a Z component."
            )

        strike_df = get_strike(self.avg.z.frame)
        if strike_df.empty:
            warnings.warn("No strike angles could be calculated.")
            return plt.subplots(subplot_kw={'projection': 'polar'})

        # Handle 180-degree ambiguity for original data
        angles = strike_df['strike_angle'].dropna()
        angles_rad = np.deg2rad(
            np.concatenate([angles, angles + 180])
        )

        bins = np.linspace(0, 2 * np.pi, num_bins + 1)
        counts, bin_edges = np.histogram(angles_rad, bins=bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        fig, ax = plt.subplots(
            figsize=self.fig_size,
            subplot_kw={'projection': 'polar'}
        )

        # Plot the original data as bars
        widths = np.diff(bin_edges)
        bars = ax.bar(
            bin_edges[:-1], counts, width=widths,
            edgecolor='k', lw=0.5, alpha=0.7, label="Original"
        )
        norm = mcolors.Normalize(vmin=0, vmax=counts.max())
        for bar, count in zip(bars, counts):
            bar.set_facecolor(plt.cm.get_cmap(self.cmap)(norm(count)))

        # Plot corrected data if provided
        if corrected_df is not None:
            strike_corr_df = get_strike(corrected_df)
            if not strike_corr_df.empty:
                angles_corr = strike_corr_df['strike_angle'].dropna()
                angles_rad_corr = np.deg2rad(
                    np.concatenate([angles_corr, angles_corr + 180])
                )
                counts_corr, _ = np.histogram(angles_rad_corr, bins=bins)
                ax.plot(
                    bin_centers, counts_corr, color='red',
                    lw=1.5, label='Corrected', drawstyle='steps-mid'
                )
            ax.legend()

        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_title(
            self.fig_title or "Geoelectric Strike Distribution",
            fontsize=self.font_size + 4, fontweight='bold', pad=20
        )

        if self.savefig:
            fig.savefig(self.savefig, dpi=self.fig_dpi)

        return fig, ax