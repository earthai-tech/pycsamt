"""Docstring fragments for :mod:`pycsamt.models.occam2d`."""

from __future__ import annotations

from ...api.docs import DocstringComponents

__all__ = ["_occam_param_docs"]


_occam_common_params = dict(
    path="""
path : path-like
    Path to an Occam2D input or output file. The value may be
    a string or :class:`pathlib.Path`; it is converted to
    :class:`pathlib.Path` before use.
""",
    workdir="""
workdir : path-like, default "."
    Directory that contains, or will receive, the Occam2D run
    files. Builders and runners create it when writing output.
""",
    config="""
config : OccamConfig, optional
    Configuration for data selection, mesh geometry, startup
    values, file names, and binary discovery. If omitted, a
    default :class:`OccamConfig` instance is created.
""",
    verbose="""
verbose : int or bool, default 0
    Verbosity level. Zero keeps the object quiet. Positive
    values emit progress messages through the instance logger.
""",
    logger="""
logger : logging.Logger, optional
    Logger used by the object. If omitted, a class-specific
    PyCSAMT logger is created automatically.
""",
)

_occam_data_params = dict(
    source="""
source : Sites, EDICollection, or iterable
    EDI-derived source used to build the data file. Iteration
    must yield site-like objects exposing frequency, apparent
    resistivity, phase, and preferably coordinates.
""",
    modes="""
modes : list of str, optional
    Modes to include. Supported values are ``"TE"`` for
    :math:`Z_{xy}` and ``"TM"`` for :math:`Z_{yx}`. When
    omitted, ``config.modes`` is used.
""",
    title="""
title : str, default "pycsamt Occam2D data file"
    Free-text title written into the Occam data-file header.
    Use it to record survey or processing provenance.
""",
    error_floor_rho="""
error_floor_rho : float, optional
    Relative apparent-resistivity error floor for data rows.
    It is converted to log10 units as
    :math:`\\sigma_d=\\sigma_{\\rho}/\\ln(10)`.
""",
    error_floor_phase="""
error_floor_phase : float, optional
    Minimum absolute phase uncertainty in degrees. This lower
    bound keeps small phase errors from dominating inversion.
""",
    freq_min="""
freq_min : float, optional
    Lower frequency limit in hertz. Frequencies below this
    value are excluded when building data from EDI sources.
""",
    freq_max="""
freq_max : float, optional
    Upper frequency limit in hertz. Frequencies above this
    value are excluded when building data from EDI sources.
""",
)

_occam_mesh_params = dict(
    n_layers="""
n_layers : int, optional
    Number of active subsurface layers. Larger values add
    vertical structure but increase inversion size.
""",
    cell_size="""
cell_size : float, optional
    Horizontal cell size in metres near stations. Overrides
    ``config.cell_size_horizontal`` for one build call.
""",
    mesh="""
mesh : OccamMesh
    Populated mesh defining horizontal and vertical finite
    elements for the Occam2D forward problem.
""",
    model="""
model : OccamModel
    Populated model-definition object that maps mesh cells to
    Occam inversion parameters.
""",
)

_occam_runner_params = dict(
    binary_path="""
binary_path : path-like, optional
    Path to a compiled Occam2D executable. If omitted, the
    runner searches ``workdir``, ``PATH``, and bundled source.
""",
    startup_file="""
startup_file : str, default "Startup"
    Name of the startup file passed to the executable. It is
    resolved relative to ``workdir``.
""",
    auto_compile="""
auto_compile : bool, default True
    If ``True``, compile bundled source when no binary is
    found during discovery.
""",
    max_iter="""
max_iter : int, optional
    Temporary override for startup iterations before launching
    the binary.
""",
    target_misfit="""
target_misfit : float, optional
    Temporary override for startup target misfit. Occam seeks
    a normalized RMS near this value.
""",
)

_occam_result_params = dict(
    iteration="""
iteration : int or None, default None
    Iteration number to load from ``workdir``. If ``None``,
    the highest-numbered ``.iter`` file is selected.
""",
    output_file="""
output_file : path-like
    Destination path for an exported ASCII result file. Parent
    directories are created when necessary.
""",
)

_occam_param_docs = (
    DocstringComponents.from_nested_components(
        common=DocstringComponents(_occam_common_params),
        data=DocstringComponents(_occam_data_params),
        mesh=DocstringComponents(_occam_mesh_params),
        runner=DocstringComponents(_occam_runner_params),
        result=DocstringComponents(_occam_result_params),
    )
)
