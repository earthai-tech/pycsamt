# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Configuration defaults for MARE2DEM input and run workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..config_io import (
    ConfigParameter,
    read_config_file,
    write_config_template,
)
from .doc import _mare2dem_param_docs as _params

__all__ = ["Mare2DEMConfig"]


_MARE2DEM_CONFIG_SCHEMA = [
    # ---- source management ----
    ConfigParameter(
        "source_dir",
        "Explicit path to the MARE2DEM Fortran source tree. "
        "Leave blank (None) to let SourceManager auto-resolve "
        "the location through the environment variable, bundled "
        "_source/ directory, or platform user-data directory.",
        "Source Management",
    ),
    ConfigParameter(
        "fc_compiler",
        "MPI-Fortran compiler for building MARE2DEM. Auto-detects "
        "mpiifort (Intel oneAPI) then mpifort when None.",
        "Source Management",
    ),
    ConfigParameter(
        "cc_compiler",
        "MPI-C compiler for building Triangle and ScaLAPACK. "
        "Auto-detects mpiicc then mpicc when None.",
        "Source Management",
    ),
    # ---- binary and MPI ----
    ConfigParameter(
        "binary",
        "Name or absolute path of the compiled MARE2DEM executable.",
        "Binary And MPI",
    ),
    ConfigParameter(
        "use_mpi",
        "Whether to launch MARE2DEM through an MPI wrapper. "
        "MARE2DEM requires MPI; set False only for special "
        "single-process debug builds.",
        "Binary And MPI",
    ),
    ConfigParameter(
        "n_procs",
        "Number of MPI processes requested when use_mpi is True.",
        "Binary And MPI",
    ),
    ConfigParameter(
        "mpi_command",
        "MPI launcher command, commonly mpirun or mpiexec.",
        "Binary And MPI",
    ),
    # ---- inversion ----
    ConfigParameter(
        "max_iterations",
        "Maximum number of Occam inversion iterations.",
        "Inversion Control",
    ),
    ConfigParameter(
        "target_rms",
        "Normalized RMS misfit target for stopping the inversion.",
        "Inversion Control",
    ),
    ConfigParameter(
        "initial_rho",
        "Starting half-space resistivity in ohm-metres.",
        "Initial Model",
    ),
    # ---- file names ----
    ConfigParameter(
        "data_file",
        "Default MARE2DEM data filename (.emdata extension).",
        "File Names",
    ),
    ConfigParameter(
        "resistivity_file",
        "Default MARE2DEM resistivity model filename (.resistivity).",
        "File Names",
    ),
    ConfigParameter(
        "settings_file",
        "Default MARE2DEM inversion-settings filename (.settings).",
        "File Names",
    ),
]


@dataclass
class Mare2DEMConfig:
    # ---- source management ----
    source_dir: str | Path | None = None
    fc_compiler: str | None = None
    cc_compiler: str | None = None

    # ---- binary and MPI ----
    binary: str = "MARE2DEM"
    use_mpi: bool = True
    n_procs: int = 4
    mpi_command: str = "mpirun"

    # ---- inversion control ----
    max_iterations: int = 150
    target_rms: float = 1.0

    # ---- initial model ----
    initial_rho: float = 1.0

    # ---- file names ----
    data_file: str = "mare2dem.emdata"
    resistivity_file: str = "mare2dem.resistivity"
    settings_file: str = "mare2dem.settings"

    # ---- derived helpers ----
    @property
    def resistivity_stem(self) -> str:
        """Return the stem of :attr:`resistivity_file` without extension."""
        return Path(self.resistivity_file).stem

    def to_template(
        self,
        path: str | Path = "mare2dem_config.py",
        *,
        fmt: str | None = None,
    ) -> Path:
        """Write this configuration as an editable template file.

        Parameters
        ----------
        path : path-like, default "mare2dem_config.py"
            Destination file.
        fmt : {"py", "json", "yml", "yaml"}, optional
            Output format. Defaults to ``.py`` when the path has no
            recognized suffix.

        Returns
        -------
        pathlib.Path
            Path of the generated template.
        """
        return write_config_template(
            path,
            self,
            _MARE2DEM_CONFIG_SCHEMA,
            fmt=fmt,
            title="MARE2DEM source-of-truth configuration",
        )

    @classmethod
    def write_template(
        cls,
        path: str | Path = "mare2dem_config.py",
        *,
        fmt: str | None = None,
    ) -> Path:
        """Write a default editable MARE2DEM configuration file.

        Parameters
        ----------
        path : path-like, default "mare2dem_config.py"
            Destination file. Suffixes ``.py``, ``.json``,
            ``.yml``, and ``.yaml`` select the output format.
        fmt : {"py", "json", "yml", "yaml"}, optional
            Explicit output format.

        Returns
        -------
        pathlib.Path
            Path of the generated source-of-truth file.

        Examples
        --------
        >>> from pycsamt.models.mare2dem.config import Mare2DEMConfig
        >>> path = Mare2DEMConfig.write_template("mare2dem_config.py")
        >>> path.name
        'mare2dem_config.py'
        """
        return cls().to_template(path, fmt=fmt)

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        *,
        strict: bool = True,
    ) -> Mare2DEMConfig:
        """Create a configuration from a source-of-truth file.

        Parameters
        ----------
        path : path-like
            Python, JSON, YML, or YAML configuration file
            generated by :meth:`write_template`.
        strict : bool, default True
            If ``True``, unknown keys raise :class:`ValueError`.

        Returns
        -------
        Mare2DEMConfig
            Configuration populated from the edited file.

        Examples
        --------
        >>> from pycsamt.models.mare2dem.config import Mare2DEMConfig
        >>> Mare2DEMConfig.write_template("mare2dem_config.json")
        PosixPath('mare2dem_config.json')
        >>> cfg = Mare2DEMConfig.from_file("mare2dem_config.json")
        >>> cfg.binary
        'MARE2DEM'
        """
        values = read_config_file(path, cls, strict=strict)
        return cls(**values)

    read = from_file


Mare2DEMConfig.__doc__ = rf"""
Collect settings that define a MARE2DEM run.

``Mare2DEMConfig`` is the central configuration object for the
MARE2DEM wrapper. It is a plain dataclass whose fields cover
source management, MPI execution, inversion control, starting
model, and default file names.

The configuration is shared by :class:`SourceManager`,
:class:`InputBuilder`, :class:`Mare2DEMRunner`, and
:class:`InversionResult` so that all components of a workflow
use consistent parameters.

Source Management
-----------------
{_params.config.source_dir}
{_params.config.fc_compiler}
{_params.config.cc_compiler}

Binary And MPI
--------------
{_params.config.binary}
{_params.config.use_mpi}
{_params.config.n_procs}
{_params.config.mpi_command}

Inversion Control
-----------------
{_params.config.max_iterations}
{_params.config.target_rms}

Initial Model
-------------
{_params.config.initial_rho}

File Names
----------
{_params.config.data_file}
{_params.config.resistivity_file}
{_params.config.settings_file}

Attributes
----------
resistivity_stem : str
    The stem of :attr:`resistivity_file` without its extension.
    MARE2DEM receives this stem as its only positional argument
    and derives the data and settings filenames from it.

Notes
-----
MARE2DEM implements 2.5-D finite-element MT and CSEM forward
modelling in the frequency domain with an Occam-style
regularized inversion [Mare2DEMConfig-1]_. The inversion minimizes an
objective of the form

.. math::

    \Phi(m) =
    \| W_d (F(m) - d) \|_2^2 +
    \\lambda \| \\nabla m \|_2^2,

where :math:`F(m)` is the 2.5-D FEM forward operator,
:math:`d` is the observed data vector, and
:math:`W_d` is the data-weighting operator. The configuration
fields ``target_rms`` and ``max_iterations`` bound the
inversion iteration.

Source-Of-Truth Files
---------------------
The recommended workflow is:

1. Generate a template with :meth:`write_template`.
2. Edit the values in the generated file.
3. Load the edited file with :meth:`from_file`.
4. Pass the configuration to :class:`SourceManager`,
   :class:`InputBuilder`, :class:`Mare2DEMRunner`.

See Also
--------
SourceManager
    Download and compile the MARE2DEM Fortran source.
InputBuilder
    Write MARE2DEM resistivity model, data, and settings files.
Mare2DEMRunner
    Launch the MARE2DEM binary subprocess.
InversionResult
    Load MARE2DEM inversion output files.

Examples
--------
Create a default configuration:

>>> from pycsamt.models.mare2dem.config import Mare2DEMConfig
>>> cfg = Mare2DEMConfig()
>>> cfg.resistivity_stem
'mare2dem'

Configure a parallel run with 8 MPI processes:

>>> cfg = Mare2DEMConfig(use_mpi=True, n_procs=8)

Point to a custom source directory:

>>> cfg = Mare2DEMConfig(source_dir="/opt/mare2dem_source")

References
----------
.. [Mare2DEMConfig-1] Key, K. (2016). MARE2DEM: A 2-D inversion code for
   controlled-source electromagnetic and magnetotelluric data.
   *Geophysical Journal International*, 207(1), 571-588.
   doi:10.1093/gji/ggw290.
"""
