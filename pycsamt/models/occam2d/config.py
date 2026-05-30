# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Default configuration for Occam2D inputs.

``OccamConfig`` holds every tuneable parameter that ``InputBuilder`` and
``OccamRunner`` need.  Override individual fields when constructing or pass
the config object explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

__all__ = ["OccamConfig"]


@dataclass
class OccamConfig:
    """Occam2D run configuration.

    Data-file options
    -----------------
    modes : list[str]
        TE/TM/Tipper components to include.  Default ``["TE", "TM"]``.
    error_floor_rho : float
        Minimum relative error on apparent resistivity (fraction).
    error_floor_phase : float
        Minimum absolute error on phase (degrees).
    freq_min, freq_max : float or None
        Frequency band limits (Hz).  ``None`` means use all EDI frequencies.

    Mesh options
    ------------
    n_layers : int
        Number of model parameter layers.
    n_airlayers : int
        Air layers above topography.
    cell_size_horizontal : float
        Horizontal cell width at station positions (m).
    cell_size_vertical_top : float
        Top-layer vertical thickness (m).
    depth_scale : float
        Geometric increase factor for layer thickness with depth.
    n_padding_x : int
        Number of horizontal padding cells on each side.

    Startup options
    ---------------
    max_iterations : int
        Maximum Occam iterations.
    target_misfit : float
        Desired normalised RMS misfit.
    roughness_type : int
        1 = gradient (default), 2 = curvature.
    diagonal_penalties : int
        0 = off (default).
    stepsize_cut_count : int
        Lagrange stepsize reduction limit.
    debug_level : int
        Occam binary debug output verbosity.
    initial_rho : float
        Starting halfspace resistivity (Ω·m).
    lagrange_start : float
        Initial Lagrange multiplier.
    """

    # --- data file ---
    modes: List[str] = field(default_factory=lambda: ["TE", "TM"])
    error_floor_rho: float = 0.05
    error_floor_phase: float = 0.5
    freq_min: Optional[float] = None
    freq_max: Optional[float] = None

    # --- mesh ---
    n_layers: int = 30
    n_airlayers: int = 5
    cell_size_horizontal: float = 100.0
    cell_size_vertical_top: float = 10.0
    depth_scale: float = 1.2
    n_padding_x: int = 7

    # --- startup ---
    max_iterations: int = 100
    target_misfit: float = 1.0
    roughness_type: int = 1
    diagonal_penalties: int = 0
    stepsize_cut_count: int = 8
    debug_level: int = 1
    initial_rho: float = 100.0
    lagrange_start: float = 5.0

    # --- file names (relative to workdir) ---
    data_file: str = "OccamDataFile.dat"
    mesh_file: str = "Occam2DMesh"
    model_file: str = "Occam2DModel"
    startup_file: str = "Startup"

    # --- binary ---
    binary_name: str = "Occam2D"
