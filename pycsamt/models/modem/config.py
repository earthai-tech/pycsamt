# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Configuration defaults for ModEM input and run workflows.

The module exposes :class:`ModEmConfig`, a dataclass that
collects data-selection, model-grid, covariance, inversion,
file-name, binary, and MPI settings used throughout the ModEM
subpackage.
"""

from __future__ import annotations

from dataclasses import dataclass

from .doc import _modem_param_docs as _params

__all__ = ["ModEmConfig"]


@dataclass
class ModEmConfig:
    # ---- dimensionality ----
    mode: str = "3d"

    # ---- data ----
    component_type: str = "Full_Impedance"
    sign_convention: str = "exp(+i\\omega t)"
    units: str = "[mV/km]/[nT]"
    error_floor_z: float = 0.05
    error_floor_z_floor: float = 0.0
    freq_min: float | None = None
    freq_max: float | None = None

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
        """Return ``True`` when :attr:`mode` selects 3-D ModEM."""
        return self.mode.strip().lower() == "3d"

    @property
    def binary_name(self) -> str:
        """Return the executable name implied by :attr:`mode`."""
        return self.binary_3d if self.is_3d else self.binary_2d


ModEmConfig.__doc__ = rf"""
Collect settings that define a ModEM run.

``ModEmConfig`` is the central configuration object for the
ModEM v2 subpackage. It is a plain dataclass, so values may be
set at construction time, changed before building files, or
shared across data, model, covariance, control, runner, and
result objects. The class does not perform file I/O or launch
the ModEM executable; it records the choices used by those
objects.

The configuration covers six workflow areas:

* dimensionality and data-component selection;
* impedance units, sign convention, and error floors;
* 2-D and 3-D starting-model geometry;
* 3-D covariance smoothing and active-cell behavior;
* nonlinear inversion-control values;
* file names, executable names, and MPI launch settings.

The starting-model classes store resistivity in logarithmic
form. For a positive half-space resistivity :math:`\rho_0`,
the initial model value is commonly

.. math::

   m_0 = \ln(\rho_0).

ModEM then updates this model while balancing data fit and
regularization [1]_, [2]_.

Dimensionality
--------------
{_params.config.mode}

Data Options
------------
{_params.config.component_type}
{_params.config.sign_convention}
{_params.config.units}
{_params.config.error_floor_z}
{_params.config.error_floor_z_floor}
{_params.config.freq_min}
{_params.config.freq_max}

2-D Grid Options
----------------
{_params.model.nx_2d}
{_params.model.nz_2d}
{_params.model.n_airlayers_2d}
{_params.model.cell_size_h_2d}
{_params.model.cell_size_v_top_2d}
{_params.model.depth_scale_2d}
{_params.model.n_padding_x_2d}

3-D Grid Options
----------------
{_params.model.nx}
{_params.model.ny}
{_params.model.nz}
{_params.model.n_airlayers}
{_params.model.cell_size_h}
{_params.model.cell_size_v_top}
{_params.model.depth_scale}
{_params.model.n_padding_xy}

Covariance Options
------------------
{_params.covariance.smooth_x}
{_params.covariance.smooth_y}
{_params.covariance.smooth_z}
{_params.covariance.n_smooth_iter}

Inversion Control
-----------------
{_params.control.max_iterations}
{_params.control.target_rms}
{_params.control.initial_lambda}
{_params.control.lambda_divisor}
{_params.control.initial_alpha}
{_params.control.rms_diff_tol}
{_params.control.lambda_exit}

Initial Model
-------------
{_params.model.initial_rho}

File Names
----------
data_file : str, default "ModEMData.dat"
    Default observed-data filename used by workflows that
    consume configuration file names. The path is interpreted
    relative to the run working directory unless callers
    provide an absolute path.
model_file : str, default "ModEM_Model.rho"
    Default model filename used by configured runner or result
    workflows. Builders may override this with mode-specific
    names such as ``"m0.ws"`` or ``"m0.rho"``.
covariance_file : str, default "ModEM.cov"
    Default 3-D covariance filename. It identifies the file
    containing smoothing coefficients and active-cell masks.
control_file : str, default "ModEM.inv"
    Default inversion-control filename. It stores nonlinear
    solver settings and output naming controls.
log_file : str, default "Modular_NLCG.log"
    Default ModEM nonlinear conjugate-gradient log filename.
    Result loaders use this name when scanning run output.
{_params.control.output_stem}

Binary And MPI
--------------
{_params.runner.binary_2d}
{_params.runner.binary_3d}
{_params.runner.use_mpi}
{_params.runner.n_procs}
{_params.runner.mpi_command}

Attributes
----------
is_3d : bool
    Derived property that returns ``True`` when ``mode`` is
    ``"3d"`` after whitespace stripping and lower-casing.
binary_name : str
    Derived property returning ``binary_3d`` for 3-D mode and
    ``binary_2d`` otherwise.

Notes
-----
The dataclass is intentionally permissive. It records user
choices but does not validate every field at construction
time. Downstream builders, readers, writers, and runners check
the subset of fields they need. This keeps quick configuration
experiments lightweight while preserving clear failure points
near the operation that requires a value.

For 2-D workflows, choose ``component_type`` from TE/TM
impedance families and use the ``*_2d`` grid fields. For 3-D
workflows, use full, off-diagonal, determinant, or vertical
component families and the 3-D grid/covariance fields.

See Also
--------
InputBuilder
    Consumes this configuration while writing ModEM input
    files.
ModEmData.from_edi
    Uses data options to select components, units, error
    floors, and frequency limits.
ModEmModel2D.halfspace
    Uses 2-D grid and initial-resistivity settings.
ModEmModel3D.halfspace
    Uses 3-D grid and initial-resistivity settings.
ModEmCovariance.from_model
    Uses covariance smoothing settings for 3-D runs.
ModEmControl.from_config
    Converts inversion-control fields into a ModEM control
    object.
ModEmRunner
    Uses binary, MPI, mode, and file-name settings.

Examples
--------
Create a default 3-D configuration:

>>> from pycsamt.models.modem.config import ModEmConfig
>>> cfg = ModEmConfig()
>>> cfg.is_3d
True
>>> cfg.binary_name
'Mod3DMT'

Configure a 2-D TE workflow:

>>> cfg = ModEmConfig(
...     mode="2d",
...     component_type="TE_Impedance",
...     nx_2d=120,
...     nz_2d=60,
... )
>>> cfg.binary_name
'Mod2DMT'

Set data weighting and frequency limits:

>>> cfg = ModEmConfig(
...     error_floor_z=0.07,
...     error_floor_z_floor=1e-6,
...     freq_min=0.01,
...     freq_max=1000.0,
... )

Configure an MPI-enabled 3-D run:

>>> cfg = ModEmConfig(
...     mode="3d",
...     use_mpi=True,
...     n_procs=16,
...     binary_3d="Mod3DMT_MPI",
... )

References
----------
.. [1] Egbert, G. D., and Kelbert, A., "Computational
   recipes for electromagnetic inverse problems", Geophysical
   Journal International, 189(1), 251-267, 2012,
   doi:10.1111/j.1365-246X.2011.05347.x.
.. [2] Kelbert, A., Meqbel, N., Egbert, G. D., and Tandon,
   K., "ModEM: A modular system for inversion of
   electromagnetic geophysical data", Computers and
   Geosciences, 66, 40-53, 2014,
   doi:10.1016/j.cageo.2014.01.010.
"""
