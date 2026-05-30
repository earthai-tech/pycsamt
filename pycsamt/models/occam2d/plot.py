# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Occam2D plotting — Python replacement for the MATLAB Occam2DMT scripts.

This module reimplements the functionality of the five original MATLAB
scripts in pure Python using matplotlib:

    plotOccam2DMT.m       → PlotModel
    plotOccam2DMTResponse.m → PlotResponse
    plotOccam2DMTPseudo.m → PlotPseudo
    plotOccamIterMisfit.m → PlotMisfit
    ExtractOccam2DMTProfile.m → PlotModel.extract_profile()

All plot classes accept an ``InversionResult`` object and expose a
``plot()`` method that returns a ``matplotlib.figure.Figure``.
"""

from __future__ import annotations

from typing import Optional, Tuple, Union

from .base import OccamBase

__all__ = [
    "PlotModel",
    "PlotResponse",
    "PlotPseudo",
    "PlotMisfit",
]


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
        raise NotImplementedError(f"{self.__class__.__name__}.plot — not yet implemented")


class PlotModel(_OccamPlotBase):
    """2-D resistivity model plot.

    Equivalent to ``plotOccam2DMT.m``.

    Parameters
    ----------
    result : InversionResult
    rho_min, rho_max : float
        Colour-scale limits (Ω·m).
    depth_max : float
        Maximum depth to display (m).
    show_stations : bool
        Overlay station triangles on the surface.
    profile_distance_unit : str
        ``'m'`` or ``'km'``.

    Methods
    -------
    plot() → Figure
    extract_profile(x0, x1) → (x, z, rho)
        Extract a vertical profile between horizontal positions x0 and x1.
        Reimplements ``ExtractOccam2DMTProfile.m``.
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
        # TODO: implement
        # 1. Extract mesh x/z grids from result.mesh
        # 2. Map result.rho_2d onto pcolormesh
        # 3. Overlay station positions if show_stations
        # 4. Add colourbar (log scale), axes labels, title
        raise NotImplementedError("PlotModel.plot — not yet implemented")

    def extract_profile(
        self, x0: float, x1: float
    ) -> tuple:
        """Extract (x, z, rho) column between profile positions *x0* and *x1*."""
        # TODO: implement — replaces ExtractOccam2DMTProfile.m
        raise NotImplementedError("PlotModel.extract_profile — not yet implemented")


class PlotResponse(_OccamPlotBase):
    """Observed vs predicted response plot.

    Equivalent to ``plotOccam2DMTResponse.m``.

    Plots apparent resistivity and phase as a function of frequency for
    selected stations.  Both TE and TM modes shown side-by-side.

    Parameters
    ----------
    result : InversionResult
    stations : list[str] or None
        Station names to plot.  ``None`` → all stations.
    modes : list[str]
        Subset of ``["TE", "TM"]``.
    period_axis : bool
        Use period (s) instead of frequency (Hz) on the x-axis.
    """

    def __init__(
        self,
        result=None,
        stations=None,
        modes: Optional[list] = None,
        period_axis: bool = True,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.stations   = stations
        self.modes      = modes or ["TE", "TM"]
        self.period_axis = period_axis

    def plot(self):
        # TODO: implement
        # 1. Pull observed & predicted from result.response
        # 2. One subplot grid per station: rho_a (top) / phase (bottom)
        # 3. Symbols for observed, lines for predicted
        raise NotImplementedError("PlotResponse.plot — not yet implemented")


class PlotPseudo(_OccamPlotBase):
    """Pseudosection of observed apparent resistivity / phase.

    Equivalent to ``plotOccam2DMTPseudo.m``.

    Parameters
    ----------
    result : InversionResult
    mode : str
        ``'TE'`` or ``'TM'``.
    data_type : str
        ``'rho'`` (apparent resistivity) or ``'phase'``.
    """

    def __init__(
        self,
        result=None,
        mode: str = "TE",
        data_type: str = "rho",
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.mode      = mode
        self.data_type = data_type

    def plot(self):
        # TODO: implement
        # 1. Reshape data_blocks into (n_freq, n_sites) grid
        # 2. pcolormesh with log scale for rho, linear for phase
        # 3. x-axis = station offsets, y-axis = log period
        raise NotImplementedError("PlotPseudo.plot — not yet implemented")


class PlotMisfit(_OccamPlotBase):
    """RMS misfit vs iteration plot.

    Equivalent to ``plotOccamIterMisfit.m``.

    Parameters
    ----------
    result : InversionResult
    show_roughness : bool
        Add a secondary y-axis for roughness.
    show_lagrange : bool
        Add a third panel for the Lagrange multiplier.
    target_line : bool
        Draw a horizontal dashed line at misfit = 1.0.
    """

    def __init__(
        self,
        result=None,
        show_roughness: bool = True,
        show_lagrange: bool  = False,
        target_line: bool    = True,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.show_roughness = show_roughness
        self.show_lagrange  = show_lagrange
        self.target_line    = target_line

    def plot(self):
        # TODO: implement
        # 1. Pull log.iterations, log.rms, log.roughness, log.lagrange
        # 2. Plot RMS vs iteration (left y-axis)
        # 3. Optionally overlay roughness (right y-axis)
        # 4. Optionally add lagrange panel below
        # 5. Dashed line at y=1 if target_line
        raise NotImplementedError("PlotMisfit.plot — not yet implemented")
