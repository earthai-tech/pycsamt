# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Validated source-of-truth configuration for Occam1D inversions.

The configuration contains only user-controlled scientific and operational
settings. Derived values, such as the geometric layer-growth factor, are
computed from the authoritative fields rather than stored independently.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, fields, replace
from numbers import Integral, Real
from pathlib import Path
from types import MappingProxyType
from typing import Any

from ..config_io import (
    ConfigParameter,
    read_config_file,
    write_config_template,
)

__all__ = ["Occam1DConfig"]

_MODES = frozenset({"xy", "yx", "te", "tm", "determinant"})

_CONFIG_SCHEMA = [
    ConfigParameter(
        "mode",
        "Sounding response used by the 1-D inversion. 'xy' and 'TE' "
        "select Zxy; 'yx' and 'TM' select Zyx. 'determinant' uses "
        "sqrt(-Zxy*Zyx) and therefore requires the complex tensor.",
        "Data Selection",
    ),
    ConfigParameter(
        "error_floor_rho",
        "Minimum relative apparent-resistivity uncertainty. The value is "
        "dimensionless: 0.05 means five percent. It is converted to the "
        "base-10 logarithmic units used by the solver data file.",
        "Data Selection",
    ),
    ConfigParameter(
        "error_floor_phase",
        "Minimum absolute phase uncertainty in degrees.",
        "Data Selection",
    ),
    ConfigParameter(
        "freq_min",
        "Optional inclusive lower frequency bound in hertz. Null leaves "
        "the lower edge open.",
        "Data Selection",
    ),
    ConfigParameter(
        "freq_max",
        "Optional inclusive upper frequency bound in hertz. Null leaves "
        "the upper edge open.",
        "Data Selection",
    ),
    ConfigParameter(
        "n_layers",
        "Number of earth layers including the bottom half-space. There are "
        "n_layers - 1 finite layer thicknesses.",
        "Layer Model",
    ),
    ConfigParameter(
        "first_thickness",
        "Thickness of the shallowest finite earth layer in metres.",
        "Layer Model",
    ),
    ConfigParameter(
        "depth_max",
        "Top depth of the bottom half-space in metres. Together with "
        "n_layers and first_thickness this uniquely determines geometric "
        "layer growth.",
        "Layer Model",
    ),
    ConfigParameter(
        "starting_resistivity",
        "Uniform starting half-space resistivity in ohm metres. Startup "
        "parameters are written as log10 of this physical value.",
        "Layer Model",
    ),
    ConfigParameter(
        "max_iterations",
        "Maximum number of inversion iterations requested from Occam1D.",
        "Inversion Control",
    ),
    ConfigParameter(
        "target_misfit",
        "Target normalized RMS data misfit. A value near one is appropriate "
        "when observational uncertainties are realistic.",
        "Inversion Control",
    ),
    ConfigParameter(
        "roughness_type",
        "Native Occam roughness penalty code. Supported values are 1 for "
        "first-difference roughness and 2 for second-difference roughness.",
        "Inversion Control",
    ),
    ConfigParameter(
        "stepsize_cut_count",
        "Maximum number of line-search step-size reductions. Zero disables "
        "additional reductions.",
        "Inversion Control",
    ),
    ConfigParameter(
        "debug_level",
        "Non-negative diagnostic level passed to the native solver.",
        "Inversion Control",
    ),
    ConfigParameter(
        "lagrange_start",
        "Strictly positive initial Lagrange multiplier.",
        "Inversion Control",
    ),
    ConfigParameter(
        "lagrange_min",
        "Absolute non-negative lower bound for every Lagrange multiplier "
        "evaluated by the native pyCSAMT solver.",
        "Inversion Control",
    ),
    ConfigParameter(
        "lagrange_max",
        "Absolute positive upper bound for every Lagrange multiplier "
        "evaluated by the native pyCSAMT solver.",
        "Inversion Control",
    ),
    ConfigParameter(
        "data_file",
        "Local filename for the native EMData_1.x input file.",
        "Native Files",
    ),
    ConfigParameter(
        "model_file",
        "Local filename for the Resistivity1DMod_1.x layer model.",
        "Native Files",
    ),
    ConfigParameter(
        "startup_file",
        "Local filename for the OCCAMITER_FLEX startup controls.",
        "Native Files",
    ),
    ConfigParameter(
        "binary_name",
        "Executable name searched in the run directory and system PATH.",
        "Native Files",
    ),
]


@dataclass
class Occam1DConfig:
    r"""Define one reproducible Occam1D inversion setup.

    All physical values use explicit SI-derived field units: frequency in
    hertz, depth and thickness in metres, phase uncertainty in degrees, and
    resistivity in ohm metres. Apparent-resistivity uncertainty is relative
    and dimensionless.

    Parameters
    ----------
    mode : {"xy", "yx", "te", "tm", "determinant"}, default="determinant"
        Electromagnetic response selected from each station. Values are
        canonicalized to lowercase during construction.
    error_floor_rho : float, default=0.05
        Relative apparent-resistivity uncertainty floor in ``(0, 1]``.
    error_floor_phase : float, default=2.0
        Absolute phase uncertainty floor in degrees, in ``(0, 90]``.
    freq_min, freq_max : float or None
        Inclusive frequency bounds in hertz. Each finite bound must be
        strictly positive and ``freq_min <= freq_max``.
    n_layers : int, default=40
        Number of earth layers including the bottom half-space. Must be at
        least two.
    first_thickness : float, default=5.0
        Thickness of the shallowest finite layer in metres.
    depth_max : float, default=5000.0
        Top depth of the bottom half-space in metres. It must be at least
        ``(n_layers - 1) * first_thickness`` so finite layers can remain
        non-decreasing in thickness.
    max_iterations : int, default=30
        Positive maximum inversion iteration count.
    target_misfit : float, default=1.0
        Strictly positive target normalized RMS.
    roughness_type : {1, 2}, default=1
        First- or second-difference model roughness penalty.
    stepsize_cut_count : int, default=8
        Non-negative maximum line-search cut count.
    debug_level : int, default=1
        Non-negative solver diagnostic level.
    starting_resistivity : float, default=100.0
        Strictly positive uniform starting resistivity in ohm metres.
    lagrange_start : float, default=5.0
        Strictly positive initial Lagrange multiplier.
    lagrange_min : float, default=1e-12
        Non-negative absolute lower multiplier bound.
    lagrange_max : float, default=1e12
        Positive absolute upper multiplier bound. It must not be below
        ``lagrange_min``, and ``lagrange_start`` must lie inside the interval.
    data_file, model_file, startup_file : str
        Distinct local native filenames. Absolute paths, parent traversal,
        empty names, and directory components are rejected.
    binary_name : str, default="Occam1D"
        Non-empty executable name. An explicit binary path belongs to
        :class:`~pycsamt.models.occam1d.Occam1DRunner`, not this portable
        configuration.

    Raises
    ------
    TypeError
        If a value has the wrong scalar type. Boolean values are not accepted
        as integers or real-valued scientific parameters.
    ValueError
        If a scalar is outside its allowed domain or fields are inconsistent.

    Notes
    -----
    Validation runs at construction time and again at public workflow
    boundaries. Direct attribute assignment remains possible for ergonomic
    interactive use, but callers must invoke :meth:`validate` after a group of
    mutations. Prefer :meth:`updated` when atomic validation is desirable.

    ``layer_growth_factor`` is derived rather than stored. This prevents the
    contradictory four-parameter layer definitions that arise when layer
    count, first thickness, maximum depth, and growth factor are all editable.

    Examples
    --------
    Create a determinant inversion configuration:

    >>> cfg = Occam1DConfig(
    ...     mode="determinant",
    ...     n_layers=40,
    ...     first_thickness=5.0,
    ...     depth_max=5000.0,
    ... )
    >>> cfg.layer_growth_factor > 1.0
    True

    Apply validated changes without mutating the original:

    >>> shallow = cfg.updated(depth_max=2000.0, n_layers=30)
    >>> shallow.depth_max
    2000.0
    >>> cfg.depth_max
    5000.0

    Generate and reload a documented YAML source of truth:

    >>> path = cfg.to_template("occam1d.yml")  # doctest: +SKIP
    >>> restored = Occam1DConfig.from_file(path)  # doctest: +SKIP
    """

    mode: str = "determinant"
    error_floor_rho: float = 0.05
    error_floor_phase: float = 2.0
    freq_min: float | None = None
    freq_max: float | None = None

    n_layers: int = 40
    first_thickness: float = 5.0
    depth_max: float = 5000.0
    starting_resistivity: float = 100.0

    max_iterations: int = 30
    target_misfit: float = 1.0
    roughness_type: int = 1
    stepsize_cut_count: int = 8
    debug_level: int = 1
    lagrange_start: float = 5.0
    lagrange_min: float = 1e-12
    lagrange_max: float = 1e12

    data_file: str = "Occam1DData"
    model_file: str = "Occam1DModel"
    startup_file: str = "Startup"
    binary_name: str = "Occam1D"

    def __post_init__(self) -> None:
        """Canonicalize text fields and validate all invariants."""
        if isinstance(self.mode, str):
            self.mode = self.mode.strip().lower()
        self.validate()

    @staticmethod
    def _real(name: str, value: Any, *, positive: bool = False) -> float:
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError(f"{name} must be a real number.")
        result = float(value)
        if not math.isfinite(result):
            raise ValueError(f"{name} must be finite.")
        if positive and result <= 0:
            raise ValueError(f"{name} must be strictly positive.")
        return result

    @staticmethod
    def _integer(name: str, value: Any, *, minimum: int) -> int:
        if isinstance(value, bool) or not isinstance(value, Integral):
            raise TypeError(f"{name} must be an integer.")
        result = int(value)
        if result < minimum:
            raise ValueError(f"{name} must be at least {minimum}.")
        return result

    @staticmethod
    def _local_filename(name: str, value: Any) -> str:
        if not isinstance(value, str):
            raise TypeError(f"{name} must be a string.")
        text = value.strip()
        if not text:
            raise ValueError(f"{name} cannot be empty.")
        path = Path(text)
        if path.is_absolute() or path.name != text or text in {".", ".."}:
            raise ValueError(f"{name} must be a local filename: {value!r}.")
        return text

    def validate(self) -> None:
        """Validate scalar domains and relationships between fields.

        Returns
        -------
        None
            Validation succeeds silently.

        Raises
        ------
        TypeError
            If a field has an unsupported scalar type.
        ValueError
            If a field is outside its domain or conflicts with another field.
        """
        if not isinstance(self.mode, str):
            raise TypeError("mode must be a string.")
        mode = self.mode.strip().lower()
        if mode not in _MODES:
            allowed = ", ".join(sorted(_MODES))
            raise ValueError(f"mode must be one of: {allowed}.")
        self.mode = mode

        rho_floor = self._real(
            "error_floor_rho", self.error_floor_rho, positive=True
        )
        if rho_floor > 1:
            raise ValueError("error_floor_rho cannot exceed 1.0.")
        phase_floor = self._real(
            "error_floor_phase", self.error_floor_phase, positive=True
        )
        if phase_floor > 90:
            raise ValueError("error_floor_phase cannot exceed 90 degrees.")

        bounds = {}
        for name in ("freq_min", "freq_max"):
            value = getattr(self, name)
            bounds[name] = (
                None
                if value is None
                else self._real(name, value, positive=True)
            )
        if (
            bounds["freq_min"] is not None
            and bounds["freq_max"] is not None
            and bounds["freq_min"] > bounds["freq_max"]
        ):
            raise ValueError("freq_min cannot exceed freq_max.")

        n_layers = self._integer("n_layers", self.n_layers, minimum=2)
        first = self._real(
            "first_thickness", self.first_thickness, positive=True
        )
        maximum = self._real("depth_max", self.depth_max, positive=True)
        minimum_depth = (n_layers - 1) * first
        tolerance = 1e-12 * max(1.0, minimum_depth)
        if maximum + tolerance < minimum_depth:
            raise ValueError(
                "depth_max must be at least (n_layers - 1) * "
                "first_thickness for non-decreasing finite layers."
            )
        self._real(
            "starting_resistivity",
            self.starting_resistivity,
            positive=True,
        )
        self._integer("max_iterations", self.max_iterations, minimum=1)
        self._real("target_misfit", self.target_misfit, positive=True)
        roughness = self._integer(
            "roughness_type", self.roughness_type, minimum=1
        )
        if roughness not in {1, 2}:
            raise ValueError("roughness_type must be 1 or 2.")
        self._integer(
            "stepsize_cut_count", self.stepsize_cut_count, minimum=0
        )
        self._integer("debug_level", self.debug_level, minimum=0)
        lagrange_start = self._real(
            "lagrange_start",
            self.lagrange_start,
            positive=True,
        )
        lagrange_min = self._real("lagrange_min", self.lagrange_min)
        lagrange_max = self._real(
            "lagrange_max",
            self.lagrange_max,
            positive=True,
        )
        if lagrange_min < 0:
            raise ValueError("lagrange_min must be non-negative.")
        if lagrange_min > lagrange_max:
            raise ValueError("lagrange_min cannot exceed lagrange_max.")
        if not lagrange_min <= lagrange_start <= lagrange_max:
            raise ValueError(
                "lagrange_start must lie between lagrange_min and "
                "lagrange_max."
            )

        native = [
            self._local_filename("data_file", self.data_file),
            self._local_filename("model_file", self.model_file),
            self._local_filename("startup_file", self.startup_file),
        ]
        if len(set(name.casefold() for name in native)) != len(native):
            raise ValueError("data, model, and startup filenames must differ.")
        self._local_filename("binary_name", self.binary_name)

        self.error_floor_rho = rho_floor
        self.error_floor_phase = phase_floor
        self.freq_min = bounds["freq_min"]
        self.freq_max = bounds["freq_max"]
        self.n_layers = n_layers
        self.first_thickness = first
        self.depth_max = maximum
        self.starting_resistivity = float(self.starting_resistivity)
        self.max_iterations = int(self.max_iterations)
        self.target_misfit = float(self.target_misfit)
        self.roughness_type = roughness
        self.stepsize_cut_count = int(self.stepsize_cut_count)
        self.debug_level = int(self.debug_level)
        self.lagrange_start = lagrange_start
        self.lagrange_min = lagrange_min
        self.lagrange_max = lagrange_max
        self.data_file, self.model_file, self.startup_file = native
        self.binary_name = self.binary_name.strip()

    @property
    def canonical_component(self) -> str:
        """Return ``xy``, ``yx``, or ``determinant`` for data extraction."""
        return {"te": "xy", "tm": "yx"}.get(self.mode, self.mode)

    @property
    def frequency_bounds(self) -> tuple[float | None, float | None]:
        """Inclusive ``(minimum, maximum)`` frequency bounds in hertz."""
        return self.freq_min, self.freq_max

    @property
    def minimum_model_depth(self) -> float:
        """Minimum depth compatible with non-decreasing thicknesses."""
        return float((self.n_layers - 1) * self.first_thickness)

    @property
    def layer_growth_factor(self) -> float:
        r"""Geometric thickness factor implied by the layer settings.

        The returned scale ``r >= 1`` satisfies

        .. math::

           d_{max} = h_0 \sum_{k=0}^{N-2} r^k,

        where ``h_0`` is ``first_thickness`` and ``N`` is ``n_layers``.
        """
        count = self.n_layers - 1
        if math.isclose(self.depth_max, self.minimum_model_depth):
            return 1.0

        def total(scale: float) -> float:
            return self.first_thickness * sum(
                scale**power for power in range(count)
            )

        lower, upper = 1.0, 2.0
        while total(upper) < self.depth_max:
            upper *= 2.0
        for _ in range(80):
            middle = 0.5 * (lower + upper)
            if total(middle) < self.depth_max:
                lower = middle
            else:
                upper = middle
        return 0.5 * (lower + upper)

    @property
    def native_filenames(self) -> tuple[str, str, str]:
        """Data, model, and startup filenames in dependency order."""
        return self.data_file, self.model_file, self.startup_file

    def run_paths(self, workdir: str | Path) -> MappingProxyType:
        """Return read-only absolute paths for native run inputs.

        Parameters
        ----------
        workdir : path-like
            Run-directory root. It is not created by this method.
        """
        root = Path(workdir).expanduser().resolve(strict=False)
        paths = {
            "workdir": root,
            "data": root / self.data_file,
            "model": root / self.model_file,
            "startup": root / self.startup_file,
        }
        return MappingProxyType(paths)

    def to_dict(self) -> dict[str, Any]:
        """Return a plain editable mapping of authoritative fields."""
        return asdict(self)

    def copy(self):
        """Return an independently validated shallow configuration copy."""
        return replace(self)

    def updated(self, **changes: Any):
        """Return a validated copy containing ``changes``.

        Unknown names raise :class:`TypeError` through
        :func:`dataclasses.replace`.
        """
        return replace(self, **changes)

    def update(self, **changes: Any):
        """Atomically validate and apply changes in place.

        The current instance is unchanged when the proposed configuration is
        invalid.
        """
        candidate = self.updated(**changes)
        for field in fields(self):
            setattr(self, field.name, getattr(candidate, field.name))
        return self

    def to_template(
        self,
        path: str | Path = "occam1d_config.py",
        *,
        fmt: str | None = None,
    ) -> Path:
        """Write this instance as a documented editable configuration."""
        self.validate()
        return write_config_template(
            path,
            self,
            _CONFIG_SCHEMA,
            fmt=fmt,
            title="Occam1D source-of-truth configuration",
        )

    @classmethod
    def write_template(
        cls,
        path: str | Path = "occam1d_config.py",
        *,
        fmt: str | None = None,
    ) -> Path:
        """Write a documented template containing validated defaults."""
        return cls().to_template(path, fmt=fmt)

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        *,
        strict: bool = True,
    ):
        """Load and validate Python, JSON, or YAML configuration values."""
        values = read_config_file(path, cls, strict=strict)
        return cls(**values)

    read = from_file
