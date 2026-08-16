# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Dual-Uncertainty Hybrid Inversion preparation for Occam2D.

This module connects two uncertainty descriptions to solver-native
Occam2D inputs:

* observation reliability modifies each datum error in data space;
* AI predictive uncertainty modifies model prejudice in parameter
  space.

For datum error :math:`\sigma_i` and reliability :math:`c_i`, the
effective error is

.. math::

    \sigma_{i,\mathrm{eff}}
    =
    \frac{\sigma_i}
    {\sqrt{\max(c_i, c_{\min})}}.

For AI parameter mean :math:`\mu_j`, standard deviation
:math:`\sigma_j`, uncertainty floor :math:`\sigma_0`, and AI weight
:math:`\lambda_{\mathrm{AI}}`, the Occam prejudice amplitude is

.. math::

    w_j
    =
    \frac{\lambda_{\mathrm{AI}}}
    {\sqrt{\sigma_j^2 + \sigma_0^2}}.

The :class:`DUHIInverter2D` object prepares these inputs but does not
execute the external solver. Execution remains the responsibility of
:class:`pycsamt.models.occam2d.OccamRunner` or the common Occam2D
inversion backend.

Entry points
------------
``apply_observation_reliability(errors, reliability)``
    Convert nominal errors to reliability-weighted effective errors.
``map_ai_grid_to_occam(grid, model, mesh)``
    Map a coordinated 2-D AI grid to Occam parameter order.
``DUHIInverter2D.prepare(builder, ...)``
    Write a completed Occam2D project with both DUHI branches.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np

from ...api.property import PyCSAMTObject
from ...compat.sklearn import validate_params
from ...models.occam2d import (
    InputBuilder,
    OccamMesh,
    OccamModel,
    OccamPrejudice,
)
from .mapping2d import map_ai_grid_to_occam
from .schema import (
    DUHI_INIT_SCHEMA,
    DUHI_PREPARE_SCHEMA,
    OBSERVATION_RELIABILITY_SCHEMA,
)

GridMapper = Callable[
    [
        np.ndarray,
        OccamModel,
        OccamMesh,
        Optional[np.ndarray],
        Optional[np.ndarray],
    ],
    np.ndarray,
]

__all__ = [
    "DUHIPreparation",
    "DUHIInverter2D",
    "apply_observation_reliability",
    "map_ai_grid_to_occam",
]


@validate_params(OBSERVATION_RELIABILITY_SCHEMA)
def apply_observation_reliability(
    errors: np.ndarray,
    reliability: np.ndarray,
    *,
    reliability_floor: float = 1.0e-6,
) -> np.ndarray:
    r"""Return reliability-weighted effective datum errors.

    Nominal errors are divided by the square root of observation
    reliability. Low-reliability data therefore receive larger
    effective errors and exert less influence on the normalized data
    misfit. Inputs are not modified.

    Parameters
    ----------
    errors : array-like of float
        Positive nominal standard errors. The returned array has the
        same shape.
    reliability : array-like of float
        Reliability values in the closed interval ``[0, 1]``. The
        values must be broadcastable to the shape of ``errors``.
    reliability_floor : float, default 1e-6
        Smallest reliability used in the denominator. It must lie in
        ``(0, 1]`` and prevents division by zero for rejected data.

    Returns
    -------
    numpy.ndarray of float
        Effective errors with the same shape as ``errors``.

    Raises
    ------
    ValueError
        Raised when inputs cannot be broadcast together, errors are
        non-finite or non-positive, reliability values are outside
        ``[0, 1]``, or the floor is invalid.

    Notes
    -----
    A reliability of one leaves the nominal error unchanged. A
    reliability of zero uses ``reliability_floor`` and therefore has
    a finite but potentially very small influence.

    See Also
    --------
    DUHIInverter2D.prepare
        Applies the transformation to an Occam data table.

    Examples
    --------
    >>> from pycsamt.ai.inversion.duhi2d import (
    ...     apply_observation_reliability,
    ... )
    >>> apply_observation_reliability([2.0, 2.0], [1.0, 0.25]).tolist()
    [2.0, 4.0]
    """
    sigma = np.asarray(errors, dtype=float)
    confidence = np.asarray(reliability, dtype=float)
    if sigma.size == 0:
        raise ValueError("errors must contain at least one value")
    if not np.isfinite(reliability_floor):
        raise ValueError("reliability_floor must be finite")
    if not 0 < reliability_floor <= 1:
        raise ValueError("reliability_floor must lie in (0, 1]")
    try:
        confidence = np.broadcast_to(confidence, sigma.shape)
    except ValueError as exc:
        raise ValueError("reliability is not broadcastable to errors") from exc
    if not np.all(np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError("errors must be finite and strictly positive")
    invalid_confidence = not np.all(np.isfinite(confidence)) or np.any(
        (confidence < 0) | (confidence > 1)
    )
    if invalid_confidence:
        raise ValueError("reliability must be finite and lie in [0, 1]")
    return sigma / np.sqrt(np.maximum(confidence, reliability_floor))


@dataclass(frozen=True)
class DUHIPreparation(PyCSAMTObject):
    """Summarize one completed DUHI Occam2D preparation.

    Parameters
    ----------
    workdir : pathlib.Path
        Occam2D project directory modified by the preparation.
    prejudice_file : pathlib.Path
        Written sparse prejudice file.
    data_file : pathlib.Path
        Rewritten Occam data file containing effective errors.
    model_file : pathlib.Path
        Rewritten model file referencing ``prejudice_file``.
    startup_file : pathlib.Path
        Rewritten startup file, optionally containing the AI mean.
    n_data : int
        Number of Occam data rows weighted by reliability.
    n_params : int
        Number of Occam model parameters receiving mapped AI values.
    effective_error_min : float
        Minimum effective error after reliability weighting.
    effective_error_max : float
        Maximum effective error after reliability weighting.
    prejudice_weight_min : float
        Minimum mapped AI prejudice weight.
    prejudice_weight_max : float
        Maximum mapped AI prejudice weight.
    ai_initialized : bool
        Whether the AI mean replaced the startup parameter vector.

    Attributes
    ----------
    files : dict of str to pathlib.Path
        Mapping of solver input roles to generated paths.

    Examples
    --------
    A preparation result is returned by
    :meth:`DUHIInverter2D.prepare`::

        result = inverter.prepare(
            builder,
            ai_mean=mean,
            ai_std=std,
            observation_reliability=reliability,
        )
        result.files["prejudice"]
    """

    workdir: Path
    prejudice_file: Path
    data_file: Path
    model_file: Path
    startup_file: Path
    n_data: int
    n_params: int
    effective_error_min: float
    effective_error_max: float
    prejudice_weight_min: float
    prejudice_weight_max: float
    ai_initialized: bool

    @property
    def files(self) -> dict[str, Path]:
        """Return generated solver-input paths by role.

        Returns
        -------
        dict of str to pathlib.Path
            Keys are ``"data"``, ``"model"``, ``"startup"``, and
            ``"prejudice"``.

        Examples
        --------
        >>> sorted(result.files)  # doctest: +SKIP
        ['data', 'model', 'prejudice', 'startup']
        """
        return {
            "data": self.data_file,
            "model": self.model_file,
            "startup": self.startup_file,
            "prejudice": self.prejudice_file,
        }


class DUHIInverter2D(PyCSAMTObject):
    r"""Prepare an Occam2D project for DUHI physics refinement.

    ``DUHIInverter2D`` converts observation reliability and ensemble
    AI uncertainty into solver-native Occam2D inputs. The object is a
    preparation-stage inversion component: it does not train a network
    and does not execute Occam2D.

    The preparation modifies a completed :class:`InputBuilder` project
    in place. It performs four operations:

    1. replace nominal datum errors by reliability-weighted errors;
    2. map AI mean and standard-deviation grids to Occam parameters;
    3. optionally replace the startup vector with the AI mean; and
    4. write uncertainty-dependent prejudice records and reference
       them from the Occam model file.

    Parameters
    ----------
    lambda_ai : float, default 1.0
        Non-negative global multiplier for AI prejudice amplitudes.
        A value of zero writes no active prejudice records.
    sigma_ai_floor : float, default 0.05
        Positive model-uncertainty floor in log10 resistivity. It
        prevents unbounded prejudice weights where ensemble spread is
        very small.
    reliability_floor : float, default 1e-6
        Smallest observation reliability used in effective errors. It
        must lie in ``(0, 1]``.
    prejudice_filename : str, default "DUHIPrejudice"
        Solver-local prejudice filename. It must be a non-empty file
        name without directory components.
    grid_mapper : callable, optional
        Function accepting ``(grid, model, mesh, x_coordinates,
        z_coordinates)`` and returning a one-dimensional vector of
        length ``model.n_params``. If omitted, the geometry-aware
        :func:`map_ai_grid_to_occam` mapper is used.
    verbose : int or bool, default 0
        Verbosity level. Positive values print a compact completion
        message after preparation.

    Attributes
    ----------
    lambda_ai : float
        Global AI prejudice multiplier.
    sigma_ai_floor : float
        Model-space uncertainty floor.
    reliability_floor : float
        Data-space reliability floor.
    prejudice_filename : str
        Native Occam prejudice filename.
    grid_mapper : callable
        Active AI-grid to Occam-parameter mapping function.
    verbose : int
        Integer verbosity level.
    is_prepared : bool
        Whether :meth:`prepare` completed successfully.
    preparation : DUHIPreparation
        Most recent completed preparation result.

    Notes
    -----
    Instances are intentionally one-shot because the builder is
    modified in place. Reapplying reliability weights to an already
    prepared data table would compound error inflation. Create a new
    ``DUHIInverter2D`` and a fresh ``InputBuilder`` for each run.

    When AI coordinates are omitted, the default mapper infers uniform
    AI cell centres spanning the complete Occam horizontal and earth
    domains. Explicit coordinates are preferred for archived runs.

    See Also
    --------
    EMInverter2D
        Produces learned 2-D model proposals.
    EnsembleInverter
        Produces ensemble mean and predictive uncertainty.
    OccamPrejudice
        Encodes model targets and weights for the native solver.
    pycsamt.models.occam2d.OccamRunner
        Executes a prepared Occam2D project.

    Examples
    --------
    Prepare a completed Occam2D project:

    >>> from pycsamt.ai.inversion import DUHIInverter2D
    >>> inverter = DUHIInverter2D(
    ...     lambda_ai=1.0,
    ...     sigma_ai_floor=0.05,
    ... )
    >>> result = inverter.prepare(  # doctest: +SKIP
    ...     builder,
    ...     ai_mean=ensemble_mean,
    ...     ai_std=ensemble_std,
    ...     observation_reliability=reliability,
    ... )
    >>> result.prejudice_file  # doctest: +SKIP
    PosixPath('occam_run/DUHIPrejudice')

    Supply a geometry-aware mapper:

    >>> inverter = DUHIInverter2D(  # doctest: +SKIP
    ...     grid_mapper=physical_mesh_mapper,
    ... )

    References
    ----------
    .. [DUHIInverter2D-1] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional models
       from magnetotelluric data", Geophysics, 55(12), 1613-1624,
       1990.
    """

    __repr_fields__ = (
        "lambda_ai",
        "sigma_ai_floor",
        "reliability_floor",
        "prejudice_filename",
        "is_prepared",
    )

    @validate_params(DUHI_INIT_SCHEMA)
    def __init__(
        self,
        *,
        lambda_ai: float = 1.0,
        sigma_ai_floor: float = 0.05,
        reliability_floor: float = 1.0e-6,
        prejudice_filename: str = "DUHIPrejudice",
        grid_mapper: GridMapper | None = None,
        verbose: int | bool = 0,
    ) -> None:
        self.lambda_ai = float(lambda_ai)
        self.sigma_ai_floor = float(sigma_ai_floor)
        self.reliability_floor = float(reliability_floor)
        self.prejudice_filename = str(prejudice_filename)
        self.grid_mapper: GridMapper = (
            map_ai_grid_to_occam if grid_mapper is None else grid_mapper
        )
        self.verbose = int(verbose)
        self._is_prepared = False
        self._preparation: DUHIPreparation | None = None
        self._prejudice: OccamPrejudice | None = None
        self._ai_mean_parameters: np.ndarray | None = None
        self._ai_std_parameters: np.ndarray | None = None
        self._prejudice_weights: np.ndarray | None = None
        self.validate()

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    def validate(self) -> None:
        """Validate DUHI preparation hyperparameters.

        Raises
        ------
        TypeError
            Raised when ``grid_mapper`` is not callable.
        ValueError
            Raised when numerical controls are non-finite or outside
            their accepted ranges, or when ``prejudice_filename`` is
            empty or contains directory components.

        See Also
        --------
        DUHIInverter2D.prepare
            Validates run-specific inputs after these controls.

        Examples
        --------
        >>> DUHIInverter2D(lambda_ai=0.5).validate()
        """
        if not np.isfinite(self.lambda_ai) or self.lambda_ai < 0:
            raise ValueError("lambda_ai must be finite and non-negative")
        if not np.isfinite(self.sigma_ai_floor) or self.sigma_ai_floor <= 0:
            raise ValueError("sigma_ai_floor must be finite and positive")
        if (
            not np.isfinite(self.reliability_floor)
            or not 0 < self.reliability_floor <= 1
        ):
            raise ValueError(
                "reliability_floor must be finite and lie in (0, 1]"
            )
        filename = Path(self.prejudice_filename)
        if (
            not self.prejudice_filename.strip()
            or filename.name != self.prejudice_filename
            or "/" in self.prejudice_filename
            or "\\" in self.prejudice_filename
        ):
            raise ValueError(
                "prejudice_filename must be a non-empty local file name"
            )
        if not callable(self.grid_mapper):
            raise TypeError("grid_mapper must be callable")

    # ------------------------------------------------------------------
    # Preparation
    # ------------------------------------------------------------------
    @validate_params(DUHI_PREPARE_SCHEMA)
    def prepare(
        self,
        builder: InputBuilder,
        *,
        ai_mean: np.ndarray,
        ai_std: np.ndarray,
        observation_reliability: np.ndarray,
        ai_initialize: bool = True,
        ai_x: np.ndarray | None = None,
        ai_z: np.ndarray | None = None,
    ) -> DUHIPreparation:
        """Apply both DUHI branches to completed Occam2D inputs.

        Parameters
        ----------
        builder : InputBuilder
            Completed Occam2D input builder. ``builder.is_ready`` must
            be ``True`` and its data, model, startup, configuration,
            and work directory must be populated. Compatible builder
            objects exposing the same interface are also accepted.
        ai_mean : array-like of float, shape (n_depth, n_horizontal)
            Finite ensemble-mean log10-resistivity grid.
        ai_std : array-like of float, shape (n_depth, n_horizontal)
            Finite non-negative predictive standard-deviation grid.
            It must have the same shape as ``ai_mean``.
        observation_reliability : array-like of float, shape (n_data,)
            One reliability value in ``[0, 1]`` for each Occam data
            row, in exactly the same order as ``data_blocks``.
        ai_initialize : bool, default True
            If ``True``, replace ``Startup.param_values`` with the
            mapped AI mean. If ``False``, retain the existing startup
            vector while still writing the DUHI prejudice.
        ai_x : array-like of float, optional
            Horizontal coordinates of AI grid columns in the Occam
            mesh coordinate system. The length must equal
            ``ai_mean.shape[1]``. If omitted, uniform centres spanning
            the complete mesh width are inferred.
        ai_z : array-like of float, optional
            AI grid depth coordinates, positive downward from the
            earth surface. The length must equal ``ai_mean.shape[0]``.
            If omitted, uniform centres spanning the earth mesh depth
            are inferred.

        Returns
        -------
        DUHIPreparation
            Immutable summary containing generated paths, dimensions,
            effective-error bounds, prejudice-weight bounds, and the
            initialization choice.

        Raises
        ------
        RuntimeError
            Raised when this one-shot object has already prepared a
            project.
        TypeError
            Raised when ``builder`` does not expose the expected
            Occam2D input interface.
        ValueError
            Raised when the builder is incomplete, AI grids are
            incompatible, reliability length is wrong, mapped vectors
            are invalid, or model/startup dimensions disagree.

        Notes
        -----
        The builder and its data, model, and startup objects are
        modified in place. The files are then rewritten sequentially
        using their standard Occam writers.

        See Also
        --------
        apply_observation_reliability
            Implements the data-space transformation.
        OccamPrejudice.from_dense
            Builds the model-space sparse constraint.
        DUHIInverter2D.preparation
            Returns the same result after successful preparation.

        Examples
        --------
        >>> result = inverter.prepare(  # doctest: +SKIP
        ...     builder,
        ...     ai_mean=mean,
        ...     ai_std=std,
        ...     observation_reliability=reliability,
        ... )
        >>> result.n_params  # doctest: +SKIP
        336
        """
        if self._is_prepared:
            raise RuntimeError(
                "DUHIInverter2D instances prepare one project only"
            )
        self._validate_builder(builder)

        mean_grid = np.asarray(ai_mean, dtype=float)
        std_grid = np.asarray(ai_std, dtype=float)
        if mean_grid.shape != std_grid.shape:
            raise ValueError("ai_mean and ai_std must have identical shape")
        if std_grid.ndim != 2 or std_grid.size == 0:
            raise ValueError("ai_mean and ai_std must be non-empty 2-D arrays")
        if not np.all(np.isfinite(mean_grid)):
            raise ValueError("ai_mean must contain only finite values")
        if not np.all(np.isfinite(std_grid)) or np.any(std_grid < 0):
            raise ValueError("ai_std must be finite and non-negative")

        data = builder.data
        mesh = builder.mesh
        model = builder.model
        startup = builder.startup
        reliability = np.asarray(
            observation_reliability,
            dtype=float,
        ).reshape(-1)
        n_data = int(data.data_blocks.shape[0])
        if reliability.size != n_data:
            raise ValueError(
                "one reliability value is required per Occam datum"
            )

        effective_errors = apply_observation_reliability(
            data.data_blocks[:, 4],
            reliability,
            reliability_floor=self.reliability_floor,
        )
        x_coordinates = None if ai_x is None else np.asarray(ai_x, dtype=float)
        z_coordinates = None if ai_z is None else np.asarray(ai_z, dtype=float)
        mean_parameters = self._map_grid(
            mean_grid,
            model,
            mesh,
            x_coordinates,
            z_coordinates,
            "ai_mean",
        )
        std_parameters = self._map_grid(
            std_grid,
            model,
            mesh,
            x_coordinates,
            z_coordinates,
            "ai_std",
        )
        if np.any(std_parameters < 0):
            raise ValueError("mapped ai_std parameters must be non-negative")
        if startup.param_values.size != int(model.n_params):
            raise ValueError(
                "startup parameter count does not match Occam model"
            )

        weights = self.lambda_ai / np.sqrt(
            std_parameters**2 + self.sigma_ai_floor**2
        )
        prejudice = OccamPrejudice.from_dense(
            mean_parameters,
            weights,
            config=builder.config,
        )
        prejudice.validate_parameter_count(model.n_params)

        workdir = Path(builder.workdir)
        data_file = workdir / builder.config.data_file
        model_file = workdir / builder.config.model_file
        startup_file = workdir / builder.config.startup_file
        prejudice_file = workdir / self.prejudice_filename

        data.data_blocks[:, 4] = effective_errors
        model.prejudice_file = self.prejudice_filename
        if ai_initialize:
            startup.param_values = mean_parameters.copy()

        prejudice.write(prejudice_file)
        data.write(data_file)
        model.write(model_file)
        startup.write(startup_file)

        result = DUHIPreparation(
            workdir=workdir,
            prejudice_file=prejudice_file,
            data_file=data_file,
            model_file=model_file,
            startup_file=startup_file,
            n_data=n_data,
            n_params=int(model.n_params),
            effective_error_min=float(effective_errors.min()),
            effective_error_max=float(effective_errors.max()),
            prejudice_weight_min=float(weights.min()),
            prejudice_weight_max=float(weights.max()),
            ai_initialized=bool(ai_initialize),
        )
        self._preparation = result
        self._prejudice = prejudice
        self._ai_mean_parameters = mean_parameters.copy()
        self._ai_std_parameters = std_parameters.copy()
        self._prejudice_weights = weights.copy()
        self._is_prepared = True

        if self.verbose:
            print(
                "DUHIInverter2D prepared "
                f"{result.n_data} data and {result.n_params} parameters "
                f"in {result.workdir}"
            )
        return result

    # ------------------------------------------------------------------
    # Prepared state
    # ------------------------------------------------------------------
    @property
    def is_prepared(self) -> bool:
        """Return whether Occam2D inputs were prepared successfully.

        Returns
        -------
        bool
            ``True`` after :meth:`prepare` completes without error.

        Examples
        --------
        >>> DUHIInverter2D().is_prepared
        False
        """
        return self._is_prepared

    @property
    def preparation(self) -> DUHIPreparation:
        """Return the completed preparation summary.

        Returns
        -------
        DUHIPreparation
            Immutable result from the successful preparation.

        Raises
        ------
        RuntimeError
            Raised when :meth:`prepare` has not completed.

        Examples
        --------
        >>> inverter.preparation  # doctest: +SKIP
        DUHIPreparation(...)
        """
        self._check_prepared()
        return self._preparation

    @property
    def prejudice(self) -> OccamPrejudice:
        """Return the generated sparse Occam prejudice object.

        Returns
        -------
        OccamPrejudice
            Generated model-space target and weight definition.

        Raises
        ------
        RuntimeError
            Raised when :meth:`prepare` has not completed.
        """
        self._check_prepared()
        return self._prejudice

    @property
    def ai_mean_parameters(self) -> np.ndarray:
        """Return a copy of the mapped AI mean parameter vector.

        Returns
        -------
        numpy.ndarray of float, shape (n_params,)
            AI mean in Occam layer-major parameter order.

        Raises
        ------
        RuntimeError
            Raised when :meth:`prepare` has not completed.
        """
        self._check_prepared()
        return self._ai_mean_parameters.copy()

    @property
    def ai_std_parameters(self) -> np.ndarray:
        """Return a copy of mapped AI standard deviations.

        Returns
        -------
        numpy.ndarray of float, shape (n_params,)
            Predictive standard deviations in Occam parameter order.

        Raises
        ------
        RuntimeError
            Raised when :meth:`prepare` has not completed.
        """
        self._check_prepared()
        return self._ai_std_parameters.copy()

    @property
    def prejudice_weights(self) -> np.ndarray:
        """Return a copy of mapped uncertainty-dependent weights.

        Returns
        -------
        numpy.ndarray of float, shape (n_params,)
            Dense prejudice weights in Occam parameter order.

        Raises
        ------
        RuntimeError
            Raised when :meth:`prepare` has not completed.
        """
        self._check_prepared()
        return self._prejudice_weights.copy()

    def summary(self, *, max_fields: int | None = None) -> str:
        """Return a compact DUHI configuration and state summary.

        Parameters
        ----------
        max_fields : int, optional
            Accepted for compatibility with
            :class:`PyCSAMTObject`. DUHI uses a fixed scientific
            summary and therefore ignores this value.

        Returns
        -------
        str
            One-line summary containing uncertainty controls and the
            preparation state.

        Examples
        --------
        >>> "unprepared" in DUHIInverter2D().summary()
        True
        """
        state = "prepared" if self.is_prepared else "unprepared"
        return (
            "DUHIInverter2D("
            f"lambda_ai={self.lambda_ai:g}, "
            f"sigma_ai_floor={self.sigma_ai_floor:g}, "
            f"reliability_floor={self.reliability_floor:g}, "
            f"{state})"
        )

    # ------------------------------------------------------------------
    # Internal validation
    # ------------------------------------------------------------------
    @staticmethod
    def _validate_builder(builder: Any) -> None:
        """Validate the minimal completed InputBuilder interface."""
        required = (
            "is_ready",
            "workdir",
            "config",
            "data",
            "mesh",
            "model",
            "startup",
        )
        missing = [name for name in required if not hasattr(builder, name)]
        if missing:
            raise TypeError(
                "builder does not provide required attributes: "
                + ", ".join(missing)
            )
        if not builder.is_ready:
            raise ValueError("builder must contain completed Occam2D inputs")
        if not hasattr(builder.data, "data_blocks"):
            raise TypeError("builder.data must provide data_blocks")
        blocks = np.asarray(builder.data.data_blocks)
        if blocks.ndim != 2 or blocks.shape[1] < 5 or blocks.shape[0] < 1:
            raise ValueError(
                "builder.data_blocks must have shape (n_data, >=5)"
            )
        if not hasattr(builder.model, "n_params"):
            raise TypeError("builder.model must provide n_params")
        if int(builder.model.n_params) < 1:
            raise ValueError("builder.model.n_params must be positive")
        if not hasattr(builder.startup, "param_values"):
            raise TypeError("builder.startup must provide param_values")

    def _map_grid(
        self,
        grid: np.ndarray,
        model: OccamModel,
        mesh: OccamMesh,
        x_coordinates: np.ndarray | None,
        z_coordinates: np.ndarray | None,
        name: str,
    ) -> np.ndarray:
        """Map and validate one AI grid against an Occam model."""
        mapped = np.asarray(
            self.grid_mapper(
                grid,
                model,
                mesh,
                x_coordinates,
                z_coordinates,
            ),
            dtype=float,
        )
        if mapped.ndim != 1:
            raise ValueError(f"mapped {name} must be one-dimensional")
        if mapped.size != int(model.n_params):
            raise ValueError(
                f"mapped {name} length does not match model.n_params"
            )
        if not np.all(np.isfinite(mapped)):
            raise ValueError(f"mapped {name} must contain finite values")
        return mapped

    def _check_prepared(self) -> None:
        """Raise when prepared-state outputs are unavailable."""
        if not self._is_prepared:
            raise RuntimeError(
                "Call prepare() before accessing DUHI preparation outputs"
            )
