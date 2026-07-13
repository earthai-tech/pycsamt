"""Docstring fragments for :mod:`pycsamt.models.occam2d`."""

from __future__ import annotations

from ...api.docs import DocstringComponents

__all__ = ["_occam_param_docs"]


_occam_common_params = dict(
    path="""
path : path-like
    Path to an Occam2D input or output file. The value may be
    a string, :class:`pathlib.Path`, or any object accepted by
    :class:`pathlib.Path`. Relative paths are interpreted from
    the current working directory of the Python process.
    Readers require the file to exist, while writers create
    parent directories when the owning method supports output.
""",
    workdir="""
workdir : path-like, default "."
    Directory that contains, or will receive, the Occam2D run
    files. Builders write ``OccamDataFile.dat``,
    ``Occam2DMesh``, ``Occam2DModel``, and ``Startup`` inside
    this directory. Runners execute the compiled binary from
    this location and capture standard output and error logs
    there. The directory is created when output is written.
""",
    config="""
config : OccamConfig, optional
    Configuration object controlling data selection, mesh
    geometry, startup controls, file names, and executable
    discovery. It centralizes choices such as modes, error
    floors, frequency limits, layer counts, cell sizes, target
    misfit, starting resistivity, and Occam file names. If
    omitted, a default :class:`OccamConfig` is created.
""",
    verbose="""
verbose : int or bool, default 0
    Verbosity level for progress reporting. ``0`` or ``False``
    keeps the object quiet. Positive values enable progress
    messages through the instance logger; larger values may be
    used by callers to request more diagnostic detail.
""",
    logger="""
logger : logging.Logger, optional
    Logger used for progress and diagnostic messages. If
    omitted, a class-specific PyCSAMT logger is created
    automatically. Pass an explicit logger when integrating
    Occam2D objects into an application-level logging setup.
""",
)

_occam_data_params = dict(
    source="""
source : Sites, EDICollection, or iterable
    EDI-derived survey source used to build the Occam data
    file. Accepted inputs include :class:`pycsamt.site.Sites`,
    an EDI collection, or any iterable of site-like objects.
    Each item must expose frequency, apparent resistivity, and
    phase arrays. Coordinates are strongly preferred because
    they allow station offsets to be ordered along profile.
    When coordinates are absent, fallback spacing is used.
""",
    modes="""
modes : list of str, optional
    Electromagnetic modes written to the data file. Supported
    values are ``"TE"`` for the :math:`Z_{xy}` component and
    ``"TM"`` for the :math:`Z_{yx}` component. Both apparent
    resistivity and phase rows are written for each selected
    mode. If omitted, ``config.modes`` is used.
""",
    title="""
title : str, default "pycsamt Occam2D data file"
    Free-text title written into the Occam data-file header.
    Use it to record survey name, processing version,
    inversion purpose, or other provenance attached to
    the generated ``OccamDataFile.dat``.
""",
    error_floor_rho="""
error_floor_rho : float, optional
    Relative apparent-resistivity error floor used when data
    built from EDI sources. The value is a fraction, for
    example ``0.05`` for five percent. Occam stores apparent
    resistivity in log10 units, so the floor is converted as
    :math:`\\sigma_d=\\sigma_{\\rho}/\\ln(10)`. This prevents
    very small resistivity errors from dominating inversion
    objective.
""",
    error_floor_phase="""
error_floor_phase : float, optional
    Minimum absolute phase uncertainty in degrees. This value
    is applied to phase rows when source errors are missing or
    smaller than the floor. It stabilizes phase weighting
    data in the Occam objective function.
""",
    freq_min="""
freq_min : float, optional
    Lower frequency limit in hertz. Frequencies below this
    value are excluded when building data from EDI sources.
    Use this to remove low-frequency samples that are noisy,
    poorly sampled, or outside the intended depth range.
""",
    freq_max="""
freq_max : float, optional
    Upper frequency limit in hertz. Frequencies above this
    value are excluded when building data from EDI sources.
    Use this to remove high-frequency samples affected by
    near-surface noise, static effects, or processing limits.
""",
)

_occam_mesh_params = dict(
    n_layers="""
n_layers : int, optional
    Number of active subsurface parameter layers in the Occam
    model. Larger values allow the inversion to represent more
    vertical structure, but they also increase the number of
    model parameters and can make the problem less stable
    without adequate data support.
""",
    cell_size="""
cell_size : float, optional
    Horizontal cell width in metres near station positions.
    Smaller values provide finer lateral resolution around the
    profile but increase mesh size and runtime. This value
    overrides ``config.cell_size_horizontal`` for one build
    only.
""",
    mesh="""
mesh : OccamMesh
    Populated mesh object defining the horizontal and vertical
    finite-element cells used by the Occam2D forward solver.
    The mesh must provide node coordinates, cell widths,
    air-layer count, and cell-type rows before mapping
    to model parameters.
""",
    model="""
model : OccamModel
    Populated model-definition object that maps mesh cells to
    Occam inversion parameters. The model must contain at
    one parameter layer and a positive parameter count before
    startup files or iteration vectors can be built from it.
""",
)

_occam_runner_params = dict(
    binary_path="""
binary_path : path-like, optional
    Explicit path to a compiled Occam2D executable. Use this
    when the binary is stored outside the run directory or is
    not on ``PATH``. If omitted, the runner searches
    ``workdir``, the system ``PATH``, and finally the bundled
    Fortran source directory when automatic compilation is
    enabled.
""",
    startup_file="""
startup_file : str, default "Startup"
    Name of the startup file passed to the executable. It is
    resolved relative to ``workdir``. The file should be an
    ``OCCAMITER_FLEX`` startup file with ``Iteration: 0`` and
    a parameter vector matching the model parameter
    count.
""",
    auto_compile="""
auto_compile : bool, default True
    If ``True``, attempt to compile the bundled Fortran source
    when no executable is found during discovery. Compilation
    requires a working Fortran compiler and ``make``. Set to
    ``False`` when the binary must be supplied explicitly.
""",
    max_iter="""
max_iter : int, optional
    Temporary override for the ``Iterations to run`` field in
    the startup file before launching the executable. This
    changes the run control file for the current run. It is
    useful for continuation or shortened test inversions.
""",
    target_misfit="""
target_misfit : float, optional
    Temporary override for startup ``Target Misfit``. Occam
    seeks a normalized RMS near this value, commonly close
    to 1.0 when data errors are well estimated. Lower values
    demand tighter data fit and may increase model roughness.
""",
)

_occam_result_params = dict(
    iteration="""
iteration : int or None, default None
    Iteration number to load from ``workdir``. If ``None``,
    the highest-numbered ``.iter`` file is selected. If the
    requested iteration is not found, result loaders may fall
    back to the last available iteration, depending on caller.
""",
    output_file="""
output_file : path-like
    Destination path for an exported ASCII result file. Parent
    directories are created when necessary. Existing files at
    the same path may be overwritten by the writer.
""",
)

_occam_param_docs = DocstringComponents.from_nested_components(
    common=DocstringComponents(_occam_common_params),
    data=DocstringComponents(_occam_data_params),
    mesh=DocstringComponents(_occam_mesh_params),
    runner=DocstringComponents(_occam_runner_params),
    result=DocstringComponents(_occam_result_params),
)
