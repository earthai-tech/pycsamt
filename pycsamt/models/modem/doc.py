"""Docstring fragments for :mod:`pycsamt.models.modem`."""

from __future__ import annotations

from ...api.docs import DocstringComponents

__all__ = ["_modem_param_docs"]


_modem_common_params = dict(
    path="""
path : path-like
    Path to a ModEM input or output file. The value may be a
    string, :class:`pathlib.Path`, or any object accepted by
    :class:`pathlib.Path`. Relative paths are interpreted from
    the current working directory of the Python process unless
    the caller resolves them against a ModEM working directory.
    Reader methods expect the file to exist. Writer methods
    create parent directories when the owning object supports
    file output.
""",
    workdir="""
workdir : path-like, default "."
    Directory that contains, or will receive, a ModEM run. The
    builder writes the data file, starting model, covariance
    file, and inversion-control file here. The runner executes
    the ModEM binary from this directory so relative file names
    in command-line arguments and control files refer to the
    same run folder. The directory is created before output is
    written.
""",
    config="""
config : ModEmConfig, optional
    Configuration object controlling dimensionality, data
    components, error floors, grid geometry, covariance
    smoothing, inversion controls, file names, executable
    names, and MPI settings. If omitted, a default
    :class:`ModEmConfig` is created. Pass an explicit
    configuration when several ModEM objects must use exactly
    the same run parameters.
""",
    verbose="""
verbose : int or bool, default 0
    Verbosity level used for progress reporting. ``0`` or
    ``False`` keeps the object quiet. Positive values enable
    diagnostic messages through the instance logger. Larger
    values may be used by callers to request more detailed
    run, parsing, or export information.
""",
    logger="""
logger : logging.Logger, optional
    Logger used for progress and diagnostic messages. If
    omitted, a class-specific PyCSAMT logger is created
    automatically. Provide a logger when integrating ModEM
    workflows into an application-wide logging configuration.
""",
)

_modem_config_params = dict(
    mode="""
mode : {"2d", "3d"}, default "3d"
    Dimensionality of the ModEM workflow. ``"2d"`` selects
    two-dimensional input formats and the ``Mod2DMT`` binary.
    ``"3d"`` selects three-dimensional formats, writes a
    covariance file, and uses the ``Mod3DMT`` binary by
    default. The value also controls which model class is
    built from survey data.
""",
    component_type="""
component_type : str, default "Full_Impedance"
    ModEM data component family written to the data file.
    Typical 2-D choices are ``"TE_Impedance"`` and
    ``"TM_Impedance"``. Common 3-D choices include
    ``"Full_Impedance"``, ``"Off_Diagonal_Impedance"``,
    ``"Determinant_Impedance"``, and
    ``"Full_Vertical_Components"``. The selected component
    controls which impedance tensor entries are exported.
""",
    sign_convention="""
sign_convention : str, default "exp(+i\\\\omega t)"
    Time-harmonic sign convention recorded in the ModEM data
    header. Use ``"exp(+i\\\\omega t)"`` or
    ``"exp(-i\\\\omega t)"`` to match the convention used by
    the impedance estimates. A mismatch changes the sign of
    imaginary components and can lead to inconsistent phase
    responses.
""",
    units="""
units : str, default "[mV/km]/[nT]"
    Impedance units written to the data-file header. The
    default corresponds to common magnetotelluric field units.
    Use ``"[V/m]/[T]"`` when the impedance tensors are stored
    in SI units. Unit consistency is important because ModEM
    interprets data values directly from the file.
""",
    error_floor_z="""
error_floor_z : float, default 0.05
    Relative impedance-error floor expressed as a fraction of
    :math:`|Z|`. A value of ``0.05`` enforces a five percent
    minimum uncertainty on each impedance component. The floor
    prevents very small formal errors from dominating the
    objective function and stabilizes inversion weighting.
""",
    error_floor_z_floor="""
error_floor_z_floor : float, default 0.0
    Absolute lower bound applied to impedance errors after the
    relative floor. Use this when some components have very
    small amplitudes and :math:`|Z|`-scaled errors alone would
    still be too small for a stable inversion.
""",
    freq_min="""
freq_min : float, optional
    Lower frequency limit in hertz. Frequencies below this
    value are excluded when building ModEM data from EDI-like
    sources. Use it to remove low-frequency samples that are
    noisy, sparsely sampled, or outside the desired depth
    range.
""",
    freq_max="""
freq_max : float, optional
    Upper frequency limit in hertz. Frequencies above this
    value are excluded when building ModEM data from EDI-like
    sources. Use it to remove high-frequency samples affected
    by near-surface noise, instrument limits, or processing
    artefacts.
""",
)

_modem_data_params = dict(
    source="""
source : iterable of site-like objects
    Survey source used to build a ModEM data file. Each item
    must expose station name, coordinates, frequency samples,
    impedance tensor values, and impedance errors in the form
    accepted by :meth:`ModEmData.from_edi`. EDI collections,
    site containers, and custom objects with compatible
    attributes can be used. Coordinates are required for
    building station locations in the ModEM coordinate system.
""",
    data="""
data : ModEmData
    Populated ModEM data object containing station names,
    station coordinates, periods, component blocks, observed
    complex values, and data errors. Builders use this object
    to derive model extents and to write the observed-data
    file. Plotting and result objects may also use it as the
    observed or predicted response container.
""",
    comment="""
comment : str, optional
    Free-text comment written near the top of a ModEM data
    file. Use it to record survey name, processing version,
    data-selection rules, or other provenance that should
    travel with the exported inversion input.
""",
    periods="""
periods : array-like of float
    Periods in seconds represented by the data object. Periods
    are the inverse of frequency, :math:`T=1/f`, and are used
    by ModEM to order response blocks. Values should be
    positive and should correspond to the impedance samples
    stored in the component rows.
""",
    site_names="""
site_names : sequence of str
    Ordered station names used in the data file and in plots.
    The order must match rows in ``site_coords`` and station
    indices in the component blocks. Stable names make it
    easier to compare observed data, predicted data, and model
    responses across inversion iterations.
""",
    site_coords="""
site_coords : array-like, shape (n_sites, 3)
    Station coordinates in metres in the ModEM coordinate
    system. Columns represent northing, easting, and elevation
    or compatible local coordinates used by the current data
    writer. The coordinate order must match ``site_names``.
""",
    blocks="""
blocks : list of dict
    Parsed or assembled ModEM component rows. Each block
    stores period, station, coordinates, component name, real
    value, imaginary value, and error. The block structure is
    used internally by readers, writers, result loaders, and
    response plotting helpers.
""",
)

_modem_model_params = dict(
    initial_rho="""
initial_rho : float, default 100.0
    Starting half-space resistivity in ohm metres. Builders
    fill active earth cells with this value before inversion.
    A positive value is required because model writers store
    resistivity in logarithmic form for ModEM-compatible files.
""",
    model="""
model : ModEmModel2D or ModEmModel3D
    Model object containing grid cell widths, origin
    coordinates, air-layer count, and resistivity values.
    Two-dimensional models represent profile and depth cells.
    Three-dimensional models represent northing, easting, and
    depth cells. The object must be populated before writing,
    plotting, or deriving a covariance file.
""",
    nx_2d="""
nx_2d : int, default 100
    Number of core horizontal cells in a 2-D ModEM model. The
    cells describe the inversion region along profile before
    lateral padding is added. Larger values can represent more
    lateral structure but increase model size and run time.
""",
    nz_2d="""
nz_2d : int, default 50
    Number of active earth layers in a 2-D ModEM model. This
    count excludes air layers. More layers allow finer depth
    variation but require enough period coverage to constrain
    the additional parameters.
""",
    n_airlayers_2d="""
n_airlayers_2d : int, default 5
    Number of air layers placed above the earth in 2-D
    models. Air layers allow the forward solver to represent
    the air-earth boundary. Their resistivity is normally fixed
    to a very high value and is not interpreted geologically.
""",
    cell_size_h_2d="""
cell_size_h_2d : float, default 100.0
    Nominal horizontal cell width in metres near the 2-D
    station zone. Smaller values increase near-station
    resolution and file size. The value should be chosen with
    station spacing and shortest useful period in mind.
""",
    cell_size_v_top_2d="""
cell_size_v_top_2d : float, default 10.0
    Thickness in metres of the shallowest earth layer in a
    2-D model. Subsequent layers grow geometrically according
    to ``depth_scale_2d``. Choose a value fine enough to
    represent shallow sensitivity without over-refining the
    mesh.
""",
    depth_scale_2d="""
depth_scale_2d : float, default 1.2
    Geometric growth factor for 2-D earth-layer thicknesses.
    Values slightly greater than one create gradually thicker
    cells with depth. Larger values reach great depths with
    fewer layers but reduce vertical resolution.
""",
    n_padding_x_2d="""
n_padding_x_2d : int, default 7
    Number of lateral padding cells added to each side of the
    2-D station zone. Padding moves artificial boundaries away
    from the survey line and helps reduce edge effects in
    forward responses.
""",
    nx="""
nx : int, default 20
    Number of core cells along the ModEM 3-D x direction,
    commonly local northing. Padding cells are added outside
    this core region. Increase the value when station coverage
    or expected structure requires finer north-south detail.
""",
    ny="""
ny : int, default 20
    Number of core cells along the ModEM 3-D y direction,
    commonly local easting. Padding cells are added outside
    this core region. Increase the value when the survey has
    dense east-west coverage or strong lateral gradients.
""",
    nz="""
nz : int, default 30
    Number of active earth layers in a 3-D ModEM model. The
    count excludes air layers. The depth range and vertical
    resolution are controlled jointly by ``cell_size_v_top``,
    ``depth_scale``, and this layer count.
""",
    n_airlayers="""
n_airlayers : int, default 5
    Number of air layers placed above the earth in 3-D
    models. These layers help represent the air-earth
    boundary in the forward solver and are usually assigned
    very high resistivity values.
""",
    cell_size_h="""
cell_size_h : float, default 500.0
    Nominal horizontal cell width in metres for 3-D core
    cells in both x and y directions. It should be consistent
    with station spacing, expected target size, and the
    shortest periods retained in the data file.
""",
    cell_size_v_top="""
cell_size_v_top : float, default 10.0
    Thickness in metres of the shallowest 3-D earth layer.
    Deeper layers grow according to ``depth_scale``. Smaller
    values improve shallow resolution but increase the number
    of cells needed to reach the same depth.
""",
    depth_scale="""
depth_scale : float, default 1.2
    Geometric growth factor for 3-D earth-layer thicknesses.
    The thickness of deeper layers approximately follows
    :math:`h_k=h_0 s^k`, where :math:`h_0` is the top-layer
    thickness and :math:`s` is this scale factor.
""",
    n_padding_xy="""
n_padding_xy : int, default 7
    Number of horizontal padding cells added on each side of
    the 3-D core grid in both x and y directions. Padding
    expands the computational domain so boundary conditions
    are farther from the survey area.
""",
)

_modem_covariance_params = dict(
    covariance="""
covariance : ModEmCovariance, optional
    Three-dimensional ModEM covariance object describing
    smoothing weights, smoothing iterations, and active-cell
    masks. It is required for standard 3-D inversions and is
    normally derived from the starting model by
    :meth:`ModEmCovariance.from_model`.
""",
    smooth_x="""
smooth_x : float, default 0.1
    Smoothing coefficient applied between neighbouring cells
    in the ModEM x direction. Larger values impose stronger
    model smoothness along this direction. The value should be
    balanced with data coverage and expected geological strike.
""",
    smooth_y="""
smooth_y : float, default 0.1
    Smoothing coefficient applied between neighbouring cells
    in the ModEM y direction. For 3-D inversions this controls
    lateral regularization perpendicular to ``smooth_x``.
""",
    smooth_z="""
smooth_z : float, default 0.1
    Smoothing coefficient applied between neighbouring cells
    in depth. Increasing this value discourages rapid vertical
    resistivity changes, while smaller values allow sharper
    layering when supported by the data.
""",
    n_smooth_iter="""
n_smooth_iter : int, default 2
    Number of times ModEM applies the covariance smoothing
    operator. ``0`` disables repeated smoothing. Higher values
    strengthen regularization and can produce smoother, more
    conservative model updates.
""",
)

_modem_control_params = dict(
    max_iterations="""
max_iterations : int, default 100
    Maximum number of nonlinear conjugate-gradient iterations
    requested in the ModEM control file. The run may stop
    earlier when the target RMS, lambda-exit threshold, or
    convergence criteria are reached.
""",
    target_rms="""
target_rms : float, default 1.05
    Target normalized root-mean-square misfit. When data
    errors are realistic, values near :math:`1` indicate a fit
    comparable to the assigned uncertainties. Smaller targets
    demand tighter data fit and may increase model roughness.
""",
    initial_lambda="""
initial_lambda : float, default 10.0
    Initial damping or regularization trade-off parameter used
    by ModEM. Larger values favour smoother, smaller model
    updates at the start of inversion. The value is reduced
    during successful iterations according to
    ``lambda_divisor``.
""",
    lambda_divisor="""
lambda_divisor : float, default 100.0
    Factor by which ModEM reduces lambda after successful
    steps. Larger divisors make lambda decrease faster, while
    smaller divisors keep stronger regularization for more
    iterations.
""",
    initial_alpha="""
initial_alpha : float, default 10.0
    Initial line-search step length used by the inversion
    control file. It controls the first trial update along the
    search direction and may be adjusted internally by ModEM
    during the nonlinear solve.
""",
    rms_diff_tol="""
rms_diff_tol : float, default 5.0e-4
    RMS-change tolerance used as a convergence or restart
    criterion. When successive RMS values differ by less than
    this threshold, the inversion is considered to have made
    little progress.
""",
    lambda_exit="""
lambda_exit : float, default 1.0e-4
    Lower lambda threshold used to stop the inversion. Once
    the regularization trade-off parameter is smaller than
    this value, additional reductions are unlikely to improve
    the solution in a stable way.
""",
    output_stem="""
output_stem : str, default "ModEM_out"
    Stem used by ModEM for generated model, response, and log
    outputs. It corresponds to the model and data output-name
    field in the control file. Use a distinctive stem when
    several inversions share one directory.
""",
)

_modem_runner_params = dict(
    model="""
model : path-like
    Starting or current ModEM model file passed to the runner.
    The path may be absolute or relative to ``workdir``. In
    2-D runs it usually points to a ``.rho`` file; in 3-D runs
    it usually points to a ``.ws`` model file. The grid,
    dimensionality, and resistivity encoding must be
    compatible with the selected executable.
""",
    model_file="""
model_file : path-like
    Starting or current ModEM model file passed to the runner.
    The path may be absolute or relative to ``workdir``. In
    2-D runs it usually points to a ``.rho`` file; in 3-D runs
    it usually points to a ``.ws`` model file.
""",
    data="""
data : path-like
    Observed-data file passed to ModEM. The file may be
    absolute or relative to ``workdir``. It should match the
    dimensionality, component selection, sign convention,
    period list, and units expected by the selected executable
    and control settings.
""",
    data_file="""
data_file : path-like
    Observed-data file passed to ModEM. The file should match
    the dimensionality, component selection, sign convention,
    and units expected by the selected executable and control
    settings.
""",
    control="""
control : path-like
    ModEM inversion-control file passed to inversion runs.
    It contains iteration limits, target misfit, lambda
    controls, output stem, line-search settings, and related
    nonlinear solver parameters. The path may be absolute or
    relative to ``workdir``.
""",
    control_file="""
control_file : path-like
    ModEM inversion-control file. It contains iteration
    limits, target misfit, lambda controls, output stem, and
    related nonlinear solver settings. The runner passes this
    file to the executable for inversion runs.
""",
    covariance="""
covariance : path-like, optional
    ModEM covariance file passed to 3-D inversion runs. The
    file describes smoothing strengths, smoothing iteration
    count, and active-cell masks. It is commonly omitted for
    2-D workflows and supplied for standard 3-D inversions.
""",
    covariance_file="""
covariance_file : path-like, optional
    ModEM covariance file used by 3-D inversion runs. The file
    describes smoothing and active-cell masks. It is commonly
    omitted for 2-D workflows and required for standard 3-D
    inversions.
""",
    binary_2d="""
binary_2d : str, default "Mod2DMT"
    Name or path of the ModEM 2-D executable. If only a name
    is supplied, the runner searches the working directory and
    the system ``PATH``. Use an absolute path when the binary
    is installed in a non-standard location.
""",
    binary_3d="""
binary_3d : str, default "Mod3DMT"
    Name or path of the ModEM 3-D executable. If only a name
    is supplied, the runner searches the working directory and
    the system ``PATH``. MPI-enabled builds can be launched
    with ``use_mpi`` and ``n_procs``.
""",
    use_mpi="""
use_mpi : bool, default False
    Whether to launch the 3-D executable through an MPI
    command. This requires a ModEM binary compiled with MPI
    support. When ``False``, the runner starts the executable
    directly as a serial process.
""",
    n_procs="""
n_procs : int, default 4
    Number of MPI processes requested when ``use_mpi`` is
    true. The value is passed to the MPI launcher and should
    be compatible with the machine, scheduler allocation, and
    ModEM build.
""",
    mpi_command="""
mpi_command : str, default "mpirun"
    MPI launcher used when ``use_mpi`` is true. Common values
    include ``"mpirun"`` and ``"mpiexec"``. Cluster
    environments may require a site-specific wrapper.
""",
    extra_args="""
extra_args : sequence of str, optional
    Additional command-line arguments appended to the ModEM
    executable invocation. Use this for advanced executable
    flags while keeping model, data, control, and covariance
    file handling under the runner.
""",
    timeout="""
timeout : float, optional
    Maximum run time in seconds for the ModEM process. If the
    process exceeds this duration, it is terminated by the
    caller. Leave as ``None`` for no Python-side timeout.
""",
    load_result="""
load_result : bool, default True
    Whether to scan ``workdir`` and return an
    :class:`InversionResult` after the executable finishes.
    Set to ``False`` when the caller only needs process
    completion and will inspect output files separately.
""",
)

_modem_builder_params = dict(
    data_filename="""
data_filename : str, default "data.dat"
    Name of the observed-data file written by
    :meth:`InputBuilder.build`. The name is resolved relative
    to ``workdir``. Use a descriptive value when several data
    selections or component sets are written into the same
    parent directory.
""",
    model_filename="""
model_filename : str, optional
    Name of the starting-model file written by the builder.
    If omitted, the builder selects ``"m0.ws"`` for 3-D runs
    and ``"m0.rho"`` for 2-D runs. The extension should
    remain compatible with the selected ModEM executable and
    model writer.
""",
    cov_filename="""
cov_filename : str, default "covariance.cov"
    Name of the covariance file written for 3-D inversions.
    The file stores smoothing weights, smoothing iteration
    count, and active-cell masks derived from the starting
    model. It is not written for 2-D builder workflows.
""",
    ctrl_filename="""
ctrl_filename : str, default "control.inv"
    Name of the ModEM inversion-control file written by the
    builder. The file contains nonlinear solver settings,
    target RMS, lambda controls, and output-stem information
    derived from :class:`ModEmConfig`.
""",
)

_modem_result_params = dict(
    result="""
result : InversionResult, optional
    Loaded ModEM inversion result containing parsed logs,
    initial and final models, iteration models, observed data,
    predicted data, and covariance where available. Plotting
    classes use this object as their shared data source.
""",
    model_initial="""
model_initial : ModEmModel2D or ModEmModel3D, optional
    Initial model loaded from the run directory. It represents
    the starting half-space or user-supplied model before
    inversion updates were applied.
""",
    model_final="""
model_final : ModEmModel2D or ModEmModel3D, optional
    Final model loaded from the run directory. When several
    iteration models are present, this is the highest-numbered
    or otherwise final model detected by the result scanner.
""",
    models="""
models : dict
    Mapping from iteration labels to parsed ModEM model
    objects. The dictionary allows callers to inspect model
    evolution across iterations rather than only the final
    inversion result.
""",
    data_obs="""
data_obs : ModEmData, optional
    Observed data loaded from the run directory. Response and
    pseudo-section plots compare this object with predicted
    data when both are available.
""",
    data_pred="""
data_pred : ModEmData, optional
    Predicted response data loaded from ModEM output. It
    should share stations, periods, and component choices with
    ``data_obs`` so residuals and response plots are
    meaningful.
""",
)

_modem_iotool_params = dict(
    impedance="""
impedance : ImpedanceFile
    Parsed impedance container storing periods, station
    locations, component names, complex impedance values, and
    errors. Conversion helpers use this object to translate
    between older block formats and newer list-style ModEM
    impedance files.
""",
    input_file="""
input_file : path-like
    Source file read by a conversion, interpolation, or export
    helper. The required format depends on the function, for
    example an old-style impedance file, Mackie model, or
    ModEM model file.
""",
    output_file="""
output_file : path-like
    Destination file written by a conversion, interpolation,
    or export helper. Parent directories are created by the
    caller when supported by that helper. Existing files may
    be overwritten.
""",
    bg_rho="""
bg_rho : float or array-like, default 100.0
    Background resistivity in ohm metres used when
    interpolating outside the source model domain. Scalars
    create a uniform background. Arrays allow depth-dependent
    or cell-dependent background values.
""",
    log_type="""
log_type : {"loge", "log10", "linear"}, default "loge"
    Resistivity encoding used by model conversion helpers.
    ``"loge"`` stores :math:`\\ln(\\rho)`, ``"log10"``
    stores :math:`\\log_{10}(\\rho)`, and ``"linear"``
    stores resistivity directly in ohm metres.
""",
)

_modem_param_docs = (
    DocstringComponents.from_nested_components(
        common=DocstringComponents(_modem_common_params),
        config=DocstringComponents(_modem_config_params),
        data=DocstringComponents(_modem_data_params),
        model=DocstringComponents(_modem_model_params),
        covariance=DocstringComponents(_modem_covariance_params),
        control=DocstringComponents(_modem_control_params),
        builder=DocstringComponents(_modem_builder_params),
        runner=DocstringComponents(_modem_runner_params),
        result=DocstringComponents(_modem_result_params),
        iotool=DocstringComponents(_modem_iotool_params),
    )
)
