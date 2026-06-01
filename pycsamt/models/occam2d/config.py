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
    r"""Collect settings that define an Occam2D run.

    ``OccamConfig`` groups the options shared by the Occam2D
    builder, file containers, runner, and result loaders.
    It is a plain dataclass, so users may set fields at
    construction time or mutate them before calling
    :class:`InputBuilder`.

    The configuration controls four parts of the workflow:

    * data selection and error floors;
    * mesh geometry and depth discretization;
    * startup and inversion-control values;
    * file names and executable discovery.

    Data Options
    ------------
    modes : list of str, default ["TE", "TM"]
        Electromagnetic modes written to the Occam data file.
        ``"TE"`` selects the :math:`Z_{xy}` component and
        ``"TM"`` selects :math:`Z_{yx}`. Each selected mode
        writes apparent-resistivity and phase rows.
    error_floor_rho : float
        Relative apparent-resistivity error floor. A value of
        ``0.05`` means five percent. Because Occam stores
        apparent resistivity as :math:`\log_{10}(\rho_a)`,
        the builder converts this floor before writing data.
    error_floor_phase : float
        Absolute phase error floor in degrees. This keeps
        phase rows with unrealistically small source errors
        from dominating the normalized data misfit.
    freq_min : float or None
        Lower frequency limit in hertz. Frequencies below this
        value are excluded when built from EDI sources.
        ``None`` leaves the lower bound open.
    freq_max : float or None
        Upper frequency limit in hertz. Frequencies above this
        value are excluded when built from EDI sources.
        ``None`` leaves the upper bound open.

    Mesh Options
    ------------
    n_layers : int
        Number of active earth layers below the air layers.
        Larger values represent more vertical structure but
        increase the parameter count.
    n_airlayers : int
        Number of air layers above the earth model. These
        layers stabilize finite-element boundaries near
        the surface.
    cell_size_horizontal : float
        Target horizontal cell width in metres near stations.
        Smaller values give finer lateral detail and bigger
        meshes.
    cell_size_vertical_top : float
        Thickness, in metres, of the top earth layer and air
        layers used by the current mesh builder.
    depth_scale : float
        Geometric multiplier applied to layer thickness with
        depth. Values greater than 1 make layers progressively
        thicker.
    n_padding_x : int
        Number of horizontal padding cells added on each side.
        Padding moves side boundaries away from the survey
        profile.

    Startup Options
    ---------------
    max_iterations : int
        Maximum number of Occam iterations requested in the
        startup file.
    target_misfit : float
        Target normalized RMS misfit. Values near 1.0 are
        typical when data errors are realistic.
    roughness_type : int
        Roughness penalty type passed to Occam. A value of
        ``1`` selects the standard gradient penalty; ``2``
        selects curvature when supported by the executable.
    diagonal_penalties : int
        Flag controlling diagonal roughness penalties in the
        startup file. ``0`` disables them.
    stepsize_cut_count : int
        Maximum number of Lagrange step-size reductions
        allowed during a line-search stage.
    debug_level : int
        Debug verbosity passed to the Occam executable.
    initial_rho : float
        Starting half-space resistivity in ohm metres. The
        startup vector is initialized as :math:`\log_{10}` of
        this value.
    lagrange_start : float
        Initial Lagrange multiplier written to startup.

    File Options
    ------------
    data_file : str
        Data filename written inside the run directory.
    mesh_file : str
        Mesh filename written inside the run directory.
    model_file : str
        Model filename written inside the run directory.
    startup_file : str
        Startup filename passed to the Occam executable.
    binary_name : str
        Executable name searched by :class:`OccamRunner`.

    Notes
    -----
    ``InputBuilder.build`` accepts one-shot overrides for
    common data and mesh fields. Overrides update the same
    ``OccamConfig`` instance stored on the builder.

    See Also
    --------
    InputBuilder
        Consumes this configuration while writing input files.
    OccamData.from_edi
        Uses data options to select modes, frequencies, and
        errors.
    OccamMesh.from_data
        Uses mesh options to build finite-element geometry.
    OccamStartup.from_model
        Uses startup options to initialize inversion controls.
    OccamRunner
        Uses file and binary settings during execution.

    Examples
    --------
    Create a standard configuration for the builder:

    >>> from pycsamt.models.occam2d import OccamConfig
    >>> from pycsamt.models.occam2d import InputBuilder
    >>> cfg = OccamConfig(n_layers=32, target_misfit=1.0)
    >>> cfg.error_floor_rho = 0.07
    >>> builder = InputBuilder([], workdir="run", config=cfg)

    Restrict the frequency band and use only TM mode:

    >>> cfg = OccamConfig(modes=["TM"])
    >>> cfg.freq_min = 0.1
    >>> cfg.freq_max = 1000.0

    Configure a finer near-station mesh:

    >>> cfg = OccamConfig()
    >>> cfg.cell_size_horizontal = 50.0
    >>> cfg.cell_size_vertical_top = 5.0
    >>> cfg.depth_scale = 1.15

    References
    ----------
    .. [1] Constable, S. C., Parker, R. L., and Constable,
       C. G., "Occam's inversion: A practical algorithm for
       generating smooth models from electromagnetic sounding
       data", Geophysics, 52(3), 289-300, 1987.
    .. [2] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
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
