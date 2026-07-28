# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Default configuration for Occam2D inputs.

``OccamConfig`` holds every tuneable parameter that
``InputBuilder`` and ``OccamRunner`` need. Override individual
fields when constructing or pass the config object explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..config_io import (
    ConfigParameter,
    read_config_file,
    write_config_template,
)

__all__ = ["OccamConfig"]


_OCCAM_CONFIG_SCHEMA = [
    ConfigParameter(
        "modes",
        "Electromagnetic modes written to the Occam data file. "
        "Use 'TE' for Zxy, 'TM' for Zyx, or both modes when "
        "both apparent resistivity and phase curves should be "
        "exported.",
        "Data Options",
    ),
    ConfigParameter(
        "error_floor_rho",
        "Relative apparent-resistivity error floor. A value of "
        "0.05 means five percent and is converted to the log10 "
        "data units used by Occam.",
        "Data Options",
    ),
    ConfigParameter(
        "error_floor_phase",
        "Minimum absolute phase uncertainty in degrees. This "
        "prevents phase rows with very small source errors from "
        "dominating the normalized misfit.",
        "Data Options",
    ),
    ConfigParameter(
        "freq_min",
        "Optional lower frequency limit in hertz. Frequencies "
        "below this value are excluded when building data from "
        "EDI-like sources.",
        "Data Options",
    ),
    ConfigParameter(
        "freq_max",
        "Optional upper frequency limit in hertz. Frequencies "
        "above this value are excluded when building data from "
        "EDI-like sources.",
        "Data Options",
    ),
    ConfigParameter(
        "n_layers",
        "Number of active earth layers in the Occam model. "
        "Larger values allow more vertical structure but "
        "increase the inversion parameter count.",
        "Mesh Options",
    ),
    ConfigParameter(
        "n_airlayers",
        "Number of air layers above the earth model. These "
        "layers help represent the air-earth boundary.",
        "Mesh Options",
    ),
    ConfigParameter(
        "cell_size_horizontal",
        "Target horizontal cell width in metres near stations. "
        "Smaller values provide finer lateral detail and larger "
        "meshes.",
        "Mesh Options",
    ),
    ConfigParameter(
        "cell_size_vertical_top",
        "Thickness in metres of the top earth layer and current " "air-layer cells.",
        "Mesh Options",
    ),
    ConfigParameter(
        "depth_scale",
        "Geometric multiplier applied to layer thickness with "
        "depth. Values greater than one make deeper layers "
        "progressively thicker.",
        "Mesh Options",
    ),
    ConfigParameter(
        "n_padding_x",
        "Number of horizontal padding cells added on each side "
        "of the survey profile.",
        "Mesh Options",
    ),
    ConfigParameter(
        "max_iterations",
        "Maximum number of Occam inversion iterations requested "
        "in the startup file.",
        "Startup Options",
    ),
    ConfigParameter(
        "target_misfit",
        "Target normalized RMS misfit. Values near one are "
        "typical when data errors are realistic.",
        "Startup Options",
    ),
    ConfigParameter(
        "roughness_type",
        "Roughness penalty type passed to Occam. The standard "
        "gradient penalty commonly uses value 1.",
        "Startup Options",
    ),
    ConfigParameter(
        "diagonal_penalties",
        "Flag controlling diagonal roughness penalties in the "
        "startup file. A value of 0 disables them.",
        "Startup Options",
    ),
    ConfigParameter(
        "stepsize_cut_count",
        "Maximum number of Lagrange step-size reductions "
        "allowed during a line-search stage.",
        "Startup Options",
    ),
    ConfigParameter(
        "debug_level",
        "Debug verbosity level passed to the Occam executable.",
        "Startup Options",
    ),
    ConfigParameter(
        "initial_rho",
        "Starting half-space resistivity in ohm metres. The "
        "startup vector is initialized as log10(initial_rho).",
        "Startup Options",
    ),
    ConfigParameter(
        "lagrange_start",
        "Initial Lagrange multiplier written to the startup file.",
        "Startup Options",
    ),
    ConfigParameter(
        "data_file",
        "Data filename written inside the Occam2D run directory.",
        "File Options",
    ),
    ConfigParameter(
        "mesh_file",
        "Mesh filename written inside the Occam2D run directory.",
        "File Options",
    ),
    ConfigParameter(
        "model_file",
        "Model-definition filename written inside the run directory.",
        "File Options",
    ),
    ConfigParameter(
        "startup_file",
        "Startup filename passed to the Occam executable.",
        "File Options",
    ),
    ConfigParameter(
        "binary_name",
        "Occam2D executable name or path used by the runner.",
        "Binary",
    ),
]


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

    Source-Of-Truth Files
    ---------------------
    Users can generate an editable configuration file before
    building an Occam2D run. Python is the default template
    format because it supports rich inline comments and can be
    read safely by :meth:`from_file` using literal parsing.
    YAML templates also keep comments. JSON templates store
    explanations in a ``"_schema"`` metadata block and editable
    values under ``"config"`` because standard JSON has no
    comment syntax.

    The recommended workflow is:

    1. Generate a template with :meth:`write_template`.
    2. Edit values in the generated file.
    3. Load the edited file with :meth:`from_file` or
       :meth:`read`.
    4. Pass the resulting configuration to builders and
       runners.

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

    Generate a documented source-of-truth template:

    >>> path = OccamConfig.write_template("occam2d_config.py")
    >>> cfg = OccamConfig.from_file(path)
    >>> cfg.binary_name
    'Occam2D'

    Use YAML when the configuration will be edited outside
    Python:

    >>> OccamConfig.write_template("occam2d_config.yml")
    PosixPath('occam2d_config.yml')

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
    modes: list[str] = field(default_factory=lambda: ["TE", "TM"])
    error_floor_rho: float = 0.05
    error_floor_phase: float = 0.5
    freq_min: float | None = None
    freq_max: float | None = None

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

    def to_template(
        self,
        path: str | Path = "occam2d_config.py",
        *,
        fmt: str | None = None,
    ) -> Path:
        """Write this configuration as an editable template.

        Parameters
        ----------
        path : path-like, default "occam2d_config.py"
            Destination file. If the path has no suffix, the
            suffix is inferred from ``fmt`` and defaults to
            ``.py``.
        fmt : {"py", "json", "yml", "yaml"}, optional
            Template format. Python and YAML templates include
            comments. JSON templates include a ``"_schema"``
            metadata block because standard JSON does not
            support comments.

        Returns
        -------
        pathlib.Path
            Path of the generated source-of-truth file.
        """
        return write_config_template(
            path,
            self,
            _OCCAM_CONFIG_SCHEMA,
            fmt=fmt,
            title="Occam2D source-of-truth configuration",
        )

    @classmethod
    def write_template(
        cls,
        path: str | Path = "occam2d_config.py",
        *,
        fmt: str | None = None,
    ) -> Path:
        """Write a default editable Occam2D configuration file.

        Parameters
        ----------
        path : path-like, default "occam2d_config.py"
            Destination file. Suffixes ``.py``, ``.json``,
            ``.yml``, and ``.yaml`` select the output format.
        fmt : {"py", "json", "yml", "yaml"}, optional
            Explicit output format. When omitted, the suffix of
            ``path`` is used; paths without a suffix produce a
            Python template.

        Returns
        -------
        pathlib.Path
            Path of the generated template.
        """
        return cls().to_template(path, fmt=fmt)

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        *,
        strict: bool = True,
    ) -> OccamConfig:
        """Create a configuration from a source-of-truth file.

        Parameters
        ----------
        path : path-like
            Python, JSON, YML, or YAML configuration file
            generated by :meth:`write_template` or following
            the same structure.
        strict : bool, default True
            If ``True``, unknown editable keys raise
            :class:`ValueError`. If ``False``, unknown keys are
            ignored. Metadata keys beginning with ``"_"`` are
            always ignored.

        Returns
        -------
        OccamConfig
            Configuration populated from edited file values.
        """
        values = read_config_file(path, cls, strict=strict)
        return cls(**values)

    read = from_file
