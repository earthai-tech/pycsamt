"""Docstring fragments for :mod:`pycsamt.inversion`."""

from __future__ import annotations

from ..api.docs import DocstringComponents

__all__ = ["_inversion_examples", "_inversion_param_docs"]


_inversion_common_params = dict(
    method="""
method : {"mt", "amt", "csamt", "emap", "tdem"}, default "mt"
    Electromagnetic method to invert. Natural-source methods use frequency,
    apparent-resistivity, and phase observations. TDEM uses time-gate decay
    observations. The value is normalized to lower case and is used with
    ``dimension`` to select the backend execution path.
""",
    dimension="""
dimension : {"1d", "2d", "3d"}, default "1d"
    Inversion dimensionality. One-dimensional runs recover layered-earth
    models. Two-dimensional runs produce section models either by stitched
    station inversions or by the built-in finite-difference profile mode.
    Three-dimensional execution is delegated to optional or external engines.
""",
    backend="""
backend : {"builtin", "simpeg", "pygimli", "occam2d", "modem"}, default "builtin"
    Inversion engine. ``"builtin"`` is dependency-light and runnable for
    local smoke studies. ``"simpeg"`` and ``"pygimli"`` are optional physics
    engines imported only when selected. ``"occam2d"`` and ``"modem"`` prepare,
    validate, and optionally run external workflows.
""",
    data="""
data : mapping, object, sequence, or path-like, optional
    Observed EM data. Mappings can contain ``freqs``/``frequencies``,
    ``rho_a``/``phase``, or ``times``/``values`` for TDEM. Station collections,
    response-like objects, and survey objects exposing ``to_soundings()`` are
    coerced through :class:`pycsamt.inversion.EMData`.
""",
    workdir="""
workdir : path-like, default "pycsamt_inversion"
    Directory used for backend files, logs, and prepared external-run inputs.
    Local backends store provenance here when they create products. External
    adapters resolve runner files relative to this folder.
""",
    backend_options="""
backend_options : dict, optional
    Backend-specific options. Common entries include mesh controls, component
    selection, ``profile_mode``, error-model overrides, SimPEG/pyGIMLi knobs,
    and external runner file names or executable settings. Unknown keys are
    passed through to the selected backend.
""",
)

_inversion_workflow_params = dict(
    config="""
config : InversionConfig or dict, optional
    Complete inversion configuration. Dictionaries are converted to
    :class:`pycsamt.inversion.config.InversionConfig`. Keyword arguments
    override matching configuration fields. If omitted, keyword arguments are
    used to build a default configuration.
""",
    data_override="""
data : mapping, object, sequence, or path-like, optional
    Optional data override for one run. When omitted, ``config.data`` is used.
    Accepted values are the same inputs handled by
    :class:`pycsamt.inversion.data.EMData`.
""",
    keyword_overrides="""
**kw
    Configuration fields that override values supplied by ``config``. Common
    keys include ``method``, ``dimension``, ``backend``, ``data``,
    ``starting_model``, ``max_iter``, ``regularization``, and
    ``backend_options``.
""",
    examples="""
Examples
--------
Run from an explicit configuration object::

    >>> from pycsamt.inversion.config import InversionConfig
    >>> from pycsamt.inversion.workflow import InversionWorkflow
    >>> cfg = InversionConfig(
    ...     method="mt",
    ...     dimension="1d",
    ...     backend="builtin",
    ...     data={"freqs": [1.0, 10.0],
    ...           "rho_a": [100.0, 120.0],
    ...           "phase": [45.0, 47.0]},
    ...     max_iter=4,
    ... )
    >>> result = InversionWorkflow(cfg).run()  # doctest: +SKIP

Run from a dictionary in one line::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion({
    ...     "method": "tdem",
    ...     "dimension": "1d",
    ...     "backend": "builtin",
    ...     "data": {"times": [1e-5, 1e-4],
    ...              "values": [1e-8, 2e-9]},
    ...     "max_iter": 4,
    ... })  # doctest: +SKIP
""",
    references="""
References
----------
.. [1] Aster, R. C., Borchers, B. and Thurber, C. H. (2018). *Parameter
       Estimation and Inverse Problems*, 3rd edition. Elsevier.
.. [2] Chave, A. D. and Jones, A. G. (2012). *The Magnetotelluric Method:
       Theory and Practice*. Cambridge University Press.
""",
)

_inversion_model_params = dict(
    n_layers="""
n_layers : int, default 4
    Number of layered-earth cells for 1-D style parameterizations, including
    the final halfspace. Built-in 1-D and stitched 2-D workflows use this value
    unless a supplied starting model defines its own layer count.
""",
    starting_model="""
starting_model : StartingModel, mapping, or object, optional
    Initial layered model. Resistivities are linear ohm-m values and
    thicknesses are metres. Mappings may use singular or plural keys such as
    ``resistivity``/``resistivities`` and ``thickness``/``thicknesses``.
""",
    reference_model="""
reference_model : StartingModel, mapping, or object, optional
    Model used as the reference for damped or smooth regularization when the
    selected backend supports it. If omitted, solvers regularize against their
    starting model or backend-native default.
""",
    regularization="""
regularization : {"none", "smooth", "damped", "blocky"}, default "smooth"
    Regularization family used by runnable backends. Shared helpers map this
    vocabulary to built-in residual penalties, pyGIMLi ``lam`` values, and
    SimPEG weighted least-squares settings where possible.
""",
    resistivities="""
resistivities : array-like of float
    Layer resistivities in ohm metres. Values must be strictly positive. The
    last value represents the bottom halfspace.
""",
    thicknesses="""
thicknesses : array-like of float
    Layer thicknesses in metres. Values must be strictly positive and the array
    length must be ``len(resistivities) - 1`` because the final layer is a
    halfspace.
""",
    name="""
name : str, optional
    Human-readable model label recorded with the object and passed to forward
    model adapters where supported.
""",
    metadata="""
metadata : dict, optional
    Free-form provenance metadata attached to the model.
""",
    examples="""
Examples
--------
Build a three-layer model with a conductive target over a resistive basement::

    >>> from pycsamt.inversion.model import StartingModel
    >>> start = StartingModel([100.0, 30.0, 800.0], [200.0, 900.0])
    >>> start.n_layers
    3
    >>> start.depths.tolist()
    [0.0, 200.0, 1100.0]

Coerce a mapping supplied in a configuration file::

    >>> from pycsamt.inversion.model import StartingModel
    >>> start = StartingModel.coerce({
    ...     "resistivities": [80.0, 250.0],
    ...     "thicknesses": [500.0],
    ... })
    >>> start.resistivities.tolist()
    [80.0, 250.0]
""",
    references="""
References
----------
.. [1] Constable, S. C., Parker, R. L. and Constable, C. G. (1987). Occam's
       inversion: A practical algorithm for generating smooth models from
       electromagnetic sounding data. *Geophysics*, 52(3), 289-300.
.. [2] Aster, R. C., Borchers, B. and Thurber, C. H. (2018). *Parameter
       Estimation and Inverse Problems*, 3rd edition. Elsevier.
""",
)

_inversion_error_params = dict(
    error_floor="""
error_floor : float, default 0.05
    Relative error floor for amplitude-like components when explicit errors
    are unavailable. A value of ``0.05`` means five percent. More detailed
    rho, phase, impedance, TDEM, and station-mask rules can be supplied through
    ``backend_options["error_model"]``.
""",
    phase_error="""
phase_error : float, default 3.0
    Absolute phase uncertainty in degrees used by MT/AMT/CSAMT objectives when
    phase errors are unavailable.
""",
    include_phase="""
include_phase : bool, default True
    Whether built-in natural-source objectives include phase observations in
    addition to apparent resistivity.
""",
    error_model_examples="""
Examples
--------
Build an error model directly::

    >>> from pycsamt.inversion.objective import ErrorModel
    >>> model = ErrorModel(rho_relative=0.05, phase_absolute=2.0)
    >>> model.errors([100.0, 200.0], component="rho").tolist()
    [5.0, 10.0]

Use component masks to down-weight stations or samples::

    >>> from pycsamt.inversion.objective import ErrorModel
    >>> model = ErrorModel(masks={"station": [True, False]})
    >>> model.mask([[1.0, 2.0], [3.0, 4.0]], component="rho").tolist()
    [[True, True], [False, False]]
""",
    error_model_references="""
References
----------
.. [1] Tarantola, A. (2005). *Inverse Problem Theory and Methods for Model
       Parameter Estimation*. SIAM.
.. [2] Aster, R. C., Borchers, B. and Thurber, C. H. (2018). *Parameter
       Estimation and Inverse Problems*, 3rd edition. Elsevier.
.. [3] Chave, A. D. and Jones, A. G. (2012). *The Magnetotelluric Method:
       Theory and Practice*. Cambridge University Press.
""",
)

_inversion_solver_params = dict(
    max_iter="""
max_iter : int, default 80
    Maximum local optimizer iterations or backend iteration request. External
    adapters may write this value into native control files or record it as
    lifecycle metadata depending on the backend.
""",
    tol="""
tol : float, default 1e-5
    Local optimizer termination tolerance used by runnable backends.
""",
    bounds="""
bounds : dict, optional
    Optional parameter bounds. Built-in layered inversions understand
    ``log10_rho`` and ``log10_thickness`` bounds. Other engines may interpret
    backend-specific bound keys through ``backend_options``.
""",
    run_external="""
run_external : bool, default False
    Whether external adapters are allowed to launch compiled executables.
    When ``False``, Occam2D and ModEM prepare and validate inputs, assemble
    commands, and return a ready/prepared result without running the binary.
""",
)

_inversion_regularization_params = dict(
    kind="""
kind : {"none", "smooth", "damped", "blocky"}, default "smooth"
    Regularization family. ``"none"`` disables the penalty, ``"smooth"``
    penalizes model roughness, ``"damped"`` combines smallness/reference and
    roughness terms, and ``"blocky"`` applies an edge-preserving normalized
    gradient penalty.
""",
    alpha_s="""
alpha_s : float, default 1.0
    Smallness or reference-model weight. This controls residual terms of the
    form :math:`m - m_{ref}` when damping/reference regularization is active.
""",
    alpha_x="""
alpha_x : float, default 1.0
    Lateral roughness weight applied along the profile or X axis.
""",
    alpha_z="""
alpha_z : float, default 1.0
    Vertical roughness weight applied along the depth or Z axis.
""",
    reference_weight="""
reference_weight : float, default 0.0
    Extra multiplier for the reference-model term. When a reference model is
    provided, values below one are promoted to one for damped/smallness terms.
""",
    metadata="""
metadata : dict, optional
    Free-form provenance metadata attached to the regularization settings.
""",
    residual_examples="""
Examples
--------
Build a smoothness penalty for a 2-D log-resistivity model::

    >>> import numpy as np
    >>> from pycsamt.inversion.regularization import Regularization
    >>> from pycsamt.inversion.regularization import regularization_residual
    >>> reg = Regularization(kind="smooth", alpha_x=2.0, alpha_z=1.0)
    >>> residual = regularization_residual(np.ones((3, 4)), regularization=reg)
    >>> residual.size
    17

Read shared settings from an inversion config::

    >>> from pycsamt.inversion.config import InversionConfig
    >>> from pycsamt.inversion.regularization import regularization_from_config
    >>> cfg = InversionConfig(regularization="damped",
    ...     backend_options={"alpha_s": 0.5, "regularization_weight": 2.0})
    >>> regularization_from_config(cfg).alpha_s
    0.5
""",
    references="""
References
----------
.. [1] Tikhonov, A. N. and Arsenin, V. Y. (1977). *Solutions of Ill-Posed
       Problems*. Winston.
.. [2] Constable, S. C., Parker, R. L. and Constable, C. G. (1987). Occam's
       inversion: A practical algorithm for generating smooth models from
       electromagnetic sounding data. *Geophysics*, 52(3), 289-300.
.. [3] Farquharson, C. G. and Oldenburg, D. W. (1998). Non-linear inversion
       using general measures of data misfit and model structure. *Geophysical
       Journal International*, 134(1), 213-227.
""",
)

_inversion_result_params = dict(
    result="""
result : InversionResult
    Backend-neutral inversion output. It stores the recovered model, mesh,
    predicted response, RMS/objective values, files, warnings, native backend
    objects, uncertainty diagnostics, convergence history, and metadata.
""",
    exports="""
exports : module
    :mod:`pycsamt.inversion.export` writes common products from any
    :class:`InversionResult`, including CSV, NPZ, GeoJSON, VTK, optional
    GeoTIFF, and portable archive snapshots.
""",
    history_examples="""
Examples
--------
Build convergence history from backend records::

    >>> from pycsamt.inversion.results import InversionHistory
    >>> history = InversionHistory(records=[
    ...     {"iteration": 0, "objective": 5.0, "rms": 2.0},
    ...     {"iteration": 1, "objective": 2.5, "rms": 1.2},
    ... ])
    >>> history.arrays()["objective"].tolist()
    [5.0, 2.5]
""",
    uncertainty_examples="""
Examples
--------
Store model confidence and sensitivity maps::

    >>> from pycsamt.inversion.results import InversionUncertainty
    >>> uncertainty = InversionUncertainty(
    ...     confidence=[[0.9, 0.8]],
    ...     sensitivity=[[1.0, 0.5]],
    ... )
    >>> uncertainty.confidence.shape
    (1, 2)
""",
    result_examples="""
Examples
--------
Create a minimal 2-D result and convert it to an interpretation model::

    >>> import numpy as np
    >>> from pycsamt.inversion.results import InversionResult
    >>> result = InversionResult(
    ...     method="mt",
    ...     dimension="2d",
    ...     backend="builtin",
    ...     model={
    ...         "rho_2d": np.array([[2.0, 2.2]]),
    ...         "x_centers": np.array([0.0, 500.0]),
    ...         "z_centers": np.array([100.0]),
    ...     },
    ... )
    >>> result.to_resistivity_model().rho_2d.shape
    (1, 2)

Summarize a result in logs or notebooks::

    >>> from pycsamt.inversion.results import InversionResult
    >>> InversionResult("mt", "1d", "builtin", rms=1.2).summary()
    "InversionResult(method='mt', dimension='1d', backend='builtin', status='success', rms=1.2)"
""",
    references="""
References
----------
.. [1] Constable, S. C., Parker, R. L. and Constable, C. G. (1987). Occam's
       inversion: A practical algorithm for generating smooth models from
       electromagnetic sounding data. *Geophysics*, 52(3), 289-300.
.. [2] Aster, R. C., Borchers, B. and Thurber, C. H. (2018). *Parameter
       Estimation and Inverse Problems*, 3rd edition. Elsevier.
""",
)

_inversion_data_params = dict(
    method="""
method : {"mt", "amt", "csamt", "emap", "tdem"}, default "mt"
    Survey method represented by the observations. The value is normalized to
    lower case. Natural-source methods require frequency-domain responses;
    ``"tdem"`` requires time-gate decay values.
""",
    frequencies="""
frequencies : array-like of float, optional
    Frequency samples in hertz for MT, AMT, CSAMT, and EMAP responses. Values
    must be strictly positive. The aliases ``freqs``, ``freq``, ``frequency``,
    and period arrays via ``periods`` are accepted by coercion helpers.
""",
    times="""
times : array-like of float, optional
    Time-gate samples in seconds for TDEM data. Values must be strictly
    positive. Sounding objects may expose these as ``time_gates``, ``times``,
    or ``time``.
""",
    rho_a="""
rho_a : array-like of float, optional
    Apparent resistivity in ohm metres. One-dimensional arrays represent one
    station. Two-dimensional arrays represent ``(n_stations, n_samples)``; if a
    response object provides ``(n_samples, n_stations)``, coercion transposes it
    when the frequency count makes the orientation unambiguous. Common aliases
    include ``rhoa``, ``app_res``, and ``apparent_resistivity``.
""",
    phase="""
phase : array-like of float, optional
    Impedance phase in degrees, using the same shape convention as ``rho_a``.
    The alias ``phi`` is accepted by object coercion.
""",
    values="""
values : array-like of float, optional
    Generic observed values. This is primarily used for TDEM decay data such as
    ``dBz_dt``/``dbdt`` arrays. One-dimensional arrays represent one sounding;
    two-dimensional arrays represent ``(n_soundings, n_times)``.
""",
    errors="""
errors : array-like of float, optional
    Absolute observation uncertainties with the same sample axis as the
    corresponding data. If omitted, backends derive errors from
    :class:`InversionConfig` error floors and component-specific
    ``ErrorModel`` settings.
""",
    station_names="""
station_names : sequence of str, optional
    Station or sounding labels. The length must match the number of station
    rows in the data. Mapping aliases include ``stations``; response-object
    aliases include ``station_names``, ``stations``, ``site_names``, and
    ``names``.
""",
    station_x="""
station_x : array-like of float, optional
    Profile-coordinate positions in metres. The length must match the number
    of station rows. Object coercion also checks aliases such as ``x``,
    ``x_stations``, ``easting``, ``longitude``, and ``lon``.
""",
    source="""
source : object, optional
    Original data object, collection, survey, or path retained for provenance
    and backend adapters. It is not modified by :class:`EMData`.
""",
    metadata="""
metadata : dict, optional
    Free-form provenance metadata. Object readers copy ``metadata`` or ``meta``
    dictionaries when available and add a ``reader`` tag describing the coercion
    path.
""",
    mapping="""
data : mapping
    Mapping accepted by :meth:`EMData.from_dict`. Keys may include
    ``method``, ``freqs``/``frequencies``, ``periods``, ``rho_a``, ``phase``,
    ``times``, ``values``, ``errors``, ``stations``/``station_names``,
    ``station_x``, and ``metadata``.
""",
    coerce="""
data : EMData, mapping, response object, sounding object, survey object, or sequence
    Input converted to :class:`EMData`. Supported objects include natural-source
    response objects with frequency and rho/phase attributes; station
    collections of such objects; TDEM sounding objects exposing time gates and
    decay values; and survey objects exposing ``to_soundings()``.
""",
    examples="""
Examples
--------
Build one MT/AMT/CSAMT station from a mapping::

    >>> from pycsamt.inversion.data import EMData
    >>> data = EMData.from_dict({
    ...     "method": "mt",
    ...     "freqs": [1.0, 10.0, 100.0],
    ...     "rho_a": [100.0, 120.0, 140.0],
    ...     "phase": [45.0, 47.0, 49.0],
    ...     "stations": ["S01"],
    ... })
    >>> data.has_mt_response
    True

Build a stitched profile from station rows::

    >>> from pycsamt.inversion.data import EMData
    >>> profile = EMData.coerce({
    ...     "method": "amt",
    ...     "freqs": [1.0, 10.0],
    ...     "rho_a": [[100.0, 120.0], [90.0, 115.0]],
    ...     "phase": [[45.0, 47.0], [43.0, 46.0]],
    ...     "station_x": [0.0, 500.0],
    ...     "station_names": ["A01", "A02"],
    ... })
    >>> profile.n_stations
    2

Build TDEM data from time gates and decay values::

    >>> from pycsamt.inversion.data import EMData
    >>> tdem = EMData.coerce({
    ...     "method": "tdem",
    ...     "times": [1e-5, 1e-4, 1e-3],
    ...     "values": [1e-8, 2e-9, 5e-10],
    ... })
    >>> tdem.has_tdem_response
    True
""",
    references="""
References
----------
.. [1] Simpson, F. and Bahr, K. (2005). *Practical Magnetotellurics*.
       Cambridge University Press.
.. [2] Chave, A. D. and Jones, A. G. (2012). *The Magnetotelluric Method:
       Theory and Practice*. Cambridge University Press.
.. [3] Nabighian, M. N. and Macnae, J. C. (1991). Time domain electromagnetic
       prospecting methods. In *Electromagnetic Methods in Applied Geophysics*,
       volume 2, SEG.
""",
)

_inversion_param_docs = DocstringComponents.from_nested_components(
    common=DocstringComponents(_inversion_common_params),
    data=DocstringComponents(_inversion_data_params),
    model=DocstringComponents(_inversion_model_params),
    errors=DocstringComponents(_inversion_error_params),
    regularization=DocstringComponents(_inversion_regularization_params),
    solver=DocstringComponents(_inversion_solver_params),
    workflow=DocstringComponents(_inversion_workflow_params),
    result=DocstringComponents(_inversion_result_params),
)

_inversion_examples = DocstringComponents(
    dict(
        mt1d="""
Examples
--------
>>> import numpy as np
>>> from pycsamt.inversion import InversionConfig, InversionWorkflow, StartingModel
>>> cfg = InversionConfig(
...     method="mt",
...     dimension="1d",
...     backend="builtin",
...     data={"freqs": np.logspace(-1, 2, 8),
...           "rho_a": np.full(8, 100.0),
...           "phase": np.full(8, 45.0)},
...     starting_model=StartingModel([80.0, 200.0], [500.0]),
...     max_iter=4,
... )
>>> result = InversionWorkflow(cfg).run()  # doctest: +SKIP
""",
        tdem1d="""
Examples
--------
>>> from pycsamt.inversion.tdem import TDEM1DInversion
>>> inv = TDEM1DInversion(
...     {"times": [1e-5, 1e-4, 1e-3], "values": [1e-8, 3e-9, 8e-10]},
...     n_layers=3,
...     max_iter=4,
... )
>>> result = inv.run()  # doctest: +SKIP
""",
        exports="""
Examples
--------
>>> from pycsamt.inversion import export
>>> export.to_csv(result, "model.csv")       # doctest: +SKIP
>>> export.to_geojson(result, "model.geojson")  # doctest: +SKIP
>>> export.to_archive(result, "run_snapshot.zip")  # doctest: +SKIP
""",
    )
)
