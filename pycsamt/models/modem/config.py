# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Default configuration for ModEM inputs.

``ModEmConfig`` holds every tuneable parameter that ``InputBuilder`` and
``ModEmRunner`` need.  Override individual fields when constructing or pass
the config object explicitly.

ModEM supports two dimensionalities controlled by the ``mode`` field:

* ``'2d'`` — runs ``Mod2DMT``; data uses TE/TM impedance components.
* ``'3d'`` — runs ``Mod3DMT``; data uses full or off-diagonal impedance tensor.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

__all__ = ["ModEmConfig"]


@dataclass
class ModEmConfig:
    """ModEM run configuration.

    Dimensionality
    --------------
    mode : str
        ``'2d'`` or ``'3d'``.  Controls which binary is used, which file
        formats apply, and whether the covariance file is written.

    Data options
    ------------
    component_type : str
        Data component to write.

        * 2D: ``'TE_Impedance'``, ``'TM_Impedance'``
        * 3D: ``'Full_Impedance'``, ``'Off_Diagonal_Impedance'``,
          ``'Determinant_Impedance'``, ``'Full_Vertical_Components'``

    sign_convention : str
        Time-harmonic convention: ``'exp(+i\\\\omega t)'`` (default) or
        ``'exp(-i\\\\omega t)'``.
    units : str
        Impedance units written in the data-file header.
        ``'[mV/km]/[nT]'`` (default) or ``'[V/m]/[T]'``.
    error_floor_z : float
        Minimum relative error on |Z| as a fraction of |Z| (default 0.05).
    error_floor_z_floor : float
        Absolute minimum error value — applied after the relative floor.
    freq_min, freq_max : float or None
        Frequency band limits (Hz).

    2-D grid options
    ----------------
    nx_2d : int      Horizontal (profile) cells.
    nz_2d : int      Vertical cells.
    n_airlayers_2d : int
    cell_size_h_2d : float   Station-zone cell width (m).
    cell_size_v_top_2d : float   Top layer thickness (m).
    depth_scale_2d : float   Geometric depth growth factor.
    n_padding_x_2d : int   Padding cells each side of station zone.

    3-D grid options
    ----------------
    nx : int        N–S (x) cells (excluding padding).
    ny : int        E–W (y) cells (excluding padding).
    nz : int        Vertical cells.
    n_airlayers : int
    cell_size_h : float   Nominal station-zone cell size in x and y (m).
    cell_size_v_top : float   Top layer thickness (m).
    depth_scale : float   Geometric depth growth factor.
    n_padding_xy : int   Padding cells each side in x and y.

    Covariance options (3D only)
    ----------------------------
    smooth_x, smooth_y, smooth_z : float
        Smoothing coefficients (0–1).  Typical value: 0.1.
    n_smooth_iter : int
        Number of times the smoothing is applied (≥ 0).

    Inversion control
    -----------------
    max_iterations : int
    target_rms : float
        Exit when RMS drops below this value (default 1.05).
    initial_lambda : float
        Initial damping factor λ.
    lambda_divisor : float
        Divide λ by this at each success (default 100).
    initial_alpha : float
        Initial line-search step (in model units).
    rms_diff_tol : float
        Restart step when RMS change is less than this (default 5e-4).
    lambda_exit : float
        Exit when λ is less than this (default 1e-4).

    Initial model
    -------------
    initial_rho : float
        Starting half-space resistivity (Ω·m).

    File names
    ----------
    data_file, model_file, covariance_file, control_file, log_file : str
        Relative to the working directory.
    output_stem : str
        Stem used by ModEM for output files (matches *Model and data output
        file name* in the control file).

    Binary / MPI
    ------------
    binary_2d, binary_3d : str
        Executable names (located via PATH or absolute path).
    use_mpi : bool
        Wrap invocation with ``mpirun`` (3D only, requires MPI build).
    n_procs : int
        Number of MPI processes.
    mpi_command : str
        MPI launcher (default ``'mpirun'``).
    """

    # ---- dimensionality ----
    mode: str = "3d"

    # ---- data ----
    component_type: str = "Full_Impedance"
    sign_convention: str = "exp(+i\\omega t)"
    units: str = "[mV/km]/[nT]"
    error_floor_z: float = 0.05
    error_floor_z_floor: float = 0.0
    freq_min: Optional[float] = None
    freq_max: Optional[float] = None

    # ---- 2D grid ----
    nx_2d: int = 100
    nz_2d: int = 50
    n_airlayers_2d: int = 5
    cell_size_h_2d: float = 100.0
    cell_size_v_top_2d: float = 10.0
    depth_scale_2d: float = 1.2
    n_padding_x_2d: int = 7

    # ---- 3D grid ----
    nx: int = 20
    ny: int = 20
    nz: int = 30
    n_airlayers: int = 5
    cell_size_h: float = 500.0
    cell_size_v_top: float = 10.0
    depth_scale: float = 1.2
    n_padding_xy: int = 7

    # ---- covariance (3D) ----
    smooth_x: float = 0.1
    smooth_y: float = 0.1
    smooth_z: float = 0.1
    n_smooth_iter: int = 2

    # ---- inversion control ----
    max_iterations: int = 100
    target_rms: float = 1.05
    initial_lambda: float = 10.0
    lambda_divisor: float = 100.0
    initial_alpha: float = 10.0
    rms_diff_tol: float = 5.0e-4
    lambda_exit: float = 1.0e-4

    # ---- initial model ----
    initial_rho: float = 100.0

    # ---- file names ----
    data_file: str = "ModEMData.dat"
    model_file: str = "ModEM_Model.rho"
    covariance_file: str = "ModEM.cov"
    control_file: str = "ModEM.inv"
    log_file: str = "Modular_NLCG.log"
    output_stem: str = "ModEM_out"

    # ---- binary ----
    binary_2d: str = "Mod2DMT"
    binary_3d: str = "Mod3DMT"

    # ---- MPI ----
    use_mpi: bool = False
    n_procs: int = 4
    mpi_command: str = "mpirun"

    # ---- derived helpers ----
    @property
    def is_3d(self) -> bool:
        return self.mode.strip().lower() == "3d"

    @property
    def binary_name(self) -> str:
        return self.binary_3d if self.is_3d else self.binary_2d
