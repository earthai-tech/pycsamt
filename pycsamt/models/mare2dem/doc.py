# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Docstring fragments for :mod:`pycsamt.models.mare2dem`."""

from __future__ import annotations

from ...api.docs import DocstringComponents

__all__ = ["_mare2dem_param_docs"]


_mare2dem_common_params = dict(
    path="""
path : path-like
    Path to a MARE2DEM input or output file. The value may be a
    string, :class:`pathlib.Path`, or any object accepted by
    :class:`pathlib.Path`. Relative paths are interpreted from
    the current working directory unless the caller resolves
    them against a MARE2DEM working directory.
""",
    workdir="""
workdir : path-like, default "."
    Directory that contains, or will receive, a MARE2DEM run.
    The builder writes the resistivity model, data, and settings
    files here. The runner executes the MARE2DEM binary from
    this directory so the input file stem is resolved relative
    to the run folder. The directory is created before output is
    written.
""",
    config="""
config : Mare2DEMConfig, optional
    Configuration object controlling the resistivity model,
    data component selection, inversion parameters, source
    management, executable name, and MPI settings. If omitted,
    a default :class:`Mare2DEMConfig` is created. Pass an
    explicit configuration when several MARE2DEM objects must
    share exactly the same run parameters.
""",
    verbose="""
verbose : int or bool, default 0
    Verbosity level. ``0`` or ``False`` keeps the object quiet.
    Positive values enable diagnostic messages through the
    instance logger. Larger values may be used to request more
    detailed run, parsing, or build information.
""",
    logger="""
logger : logging.Logger, optional
    Logger for progress and diagnostic messages. If omitted, a
    class-specific PyCSAMT logger is created automatically.
    Provide a logger when integrating MARE2DEM workflows into an
    application-wide logging configuration.
""",
)

_mare2dem_config_params = dict(
    source_dir="""
source_dir : path-like or None, default None
    Explicit path to the MARE2DEM Fortran source tree. When
    ``None``, :class:`SourceManager` applies a four-level
    fallback: the ``PYCSAMT_MARE2DEM_SOURCE`` environment
    variable, the bundled ``_source/`` directory inside the
    package (writable dev installs), and finally the platform
    user-data directory (PyPI installs). Set this field when the
    source tree lives at an unconventional location.
""",
    binary="""
binary : str, default "MARE2DEM"
    Name or absolute path of the compiled MARE2DEM executable.
    The runner first checks whether the name is on ``PATH``, then
    looks inside the resolved source directory, and finally
    checks the user-data directory. Use an absolute path to pin
    a specific build.
""",
    use_mpi="""
use_mpi : bool, default True
    Whether to launch MARE2DEM through an MPI wrapper. MARE2DEM
    is an MPI-parallel code; serial execution is not supported
    by the official distribution. Set to ``False`` only when
    testing with a special single-process build.
""",
    n_procs="""
n_procs : int, default 4
    Number of MPI processes requested when ``use_mpi`` is
    ``True``. Each process handles a subset of the spatial
    wavenumbers. Typical values range from 4 to the number of
    physical cores available on the node.
""",
    mpi_command="""
mpi_command : str, default "mpirun"
    MPI launcher command. Common alternatives are ``"mpiexec"``
    and ``"srun"`` (SLURM). The value is prepended to the
    MARE2DEM command line when ``use_mpi`` is ``True``.
""",
    fc_compiler="""
fc_compiler : str or None, default None
    MPI-Fortran compiler used to build MARE2DEM from source. When
    ``None``, :meth:`SourceManager.build` auto-detects in the
    order ``mpiifort`` (Intel oneAPI), ``mpifort`` (generic
    OpenMPI/MPICH). Intel compilers are preferred because the
    MARE2DEM Makefile uses Intel-specific pre-processor flags and
    the Intel MKL is required.
""",
    cc_compiler="""
cc_compiler : str or None, default None
    MPI-C compiler used for the Triangle mesh library and
    ScaLAPACK. Auto-detects ``mpiicc`` then ``mpicc`` when
    ``None``.
""",
    max_iterations="""
max_iterations : int, default 150
    Maximum number of Occam inversion iterations MARE2DEM is
    allowed to perform before stopping regardless of the misfit
    target.
""",
    target_rms="""
target_rms : float, default 1.0
    Normalized RMS misfit target. MARE2DEM stops iterating when
    the data misfit falls below this value. Lower values demand
    a tighter fit and may produce more structured models.
""",
    initial_rho="""
initial_rho : float, default 1.0
    Starting half-space resistivity in ohm-metres for the
    initial homogeneous resistivity model. MARE2DEM internally
    uses log10-resistivity, so the value must be positive.
""",
    data_file="""
data_file : str, default "mare2dem.emdata"
    Default MARE2DEM data filename. The file uses the ``.emdata``
    extension and records source–receiver geometry, observed data,
    and data uncertainties for one or more EM methods (MT, CSEM,
    or combined).
""",
    resistivity_file="""
resistivity_file : str, default "mare2dem.resistivity"
    Default MARE2DEM resistivity model filename. The file
    describes the 2-D finite-element mesh and layer resistivities
    and is the primary input file stem passed to the binary.
""",
    settings_file="""
settings_file : str, default "mare2dem.settings"
    Default MARE2DEM inversion-settings filename. It controls
    inversion type (Occam vs NLCG), iteration limits, target
    misfit, output verbosity, and optional regularization
    parameters.
""",
)

_mare2dem_runner_params = dict(
    resistivity_stem="""
resistivity_stem : str or path-like
    File stem of the MARE2DEM resistivity model passed to the
    binary as the sole positional argument. MARE2DEM derives the
    data filename by replacing the stem extension with
    ``.emdata`` and the settings filename with ``.settings``.
    The stem may be absolute or relative to ``workdir``.
""",
)

_mare2dem_source_params = dict(
    repo_url="""
repo_url : str, default "https://bitbucket.org/mare2dem/mare2dem_source"
    Bitbucket URL of the MARE2DEM source repository used by
    :meth:`SourceManager.download` when ``method="git"``.
""",
    archive_url="""
archive_url : str
    Direct URL of the MARE2DEM source archive used when
    ``method="archive"`` or when ``git`` is not available.
    The URL points to the current tip-of-master archive on
    Bitbucket.
""",
    method="""
method : {"git", "archive", "auto"}, default "auto"
    Download strategy. ``"git"`` clones the repository with the
    system ``git`` executable. ``"archive"`` downloads and
    extracts a ``.tar.gz`` archive using :mod:`requests`.
    ``"auto"`` tries ``"git"`` first and falls back to
    ``"archive"`` if ``git`` is not on ``PATH``.
""",
    force="""
force : bool, default False
    When ``True``, re-download even if the source tree is
    already present. When ``False``, :meth:`download` is a
    no-op if the source tree is already populated.
""",
    inc_file="""
inc_file : path-like or None, default None
    Path to a custom Make include file passed to
    ``make INCLUDE=<inc_file>``. When ``None``,
    :meth:`build` auto-generates a platform-specific include
    file by detecting the available MPI-Fortran/C compilers and
    the Intel MKL root. Provide an explicit file to override
    compiler flags for a specific cluster or toolchain.
""",
)

_mare2dem_param_docs = DocstringComponents.from_nested_components(
    common=DocstringComponents(_mare2dem_common_params),
    config=DocstringComponents(_mare2dem_config_params),
    runner=DocstringComponents(_mare2dem_runner_params),
    source=DocstringComponents(_mare2dem_source_params),
)
