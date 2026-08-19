# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Vertical-earth discretization and native Occam1D model I/O.

The scientific model stores layer-top depths and physical resistivity. Native
regularization columns are retained explicitly so reading and writing a model
does not silently discard solver controls.
"""

from __future__ import annotations

import math
import re
from numbers import Real
from pathlib import Path
from typing import Any

import numpy as np

from ...compat.sklearn import validate_params
from .base import Occam1DBase
from .schema import MODEL_BUILD_SCHEMA, MODEL_SCHEMA, PATH_SCHEMA

__all__ = ["Occam1DModel"]

_FORMAT = "Resistivity1DMod_1.0"
_FORMAT_TOKEN = "resistivity1dmod_1."
_LAYERS = re.compile(r"^\s*#\s*layers\s*:\s*(.*?)\s*$", re.I)
_FORTRAN_EXPONENT = re.compile(r"(?<=\d)[Dd](?=[-+]?\d)")


def _native_float(token: str, *, line_number: int, field: str) -> float:
    """Parse one native float, including Fortran ``D`` exponents."""
    try:
        return float(_FORTRAN_EXPONENT.sub("E", token))
    except ValueError as error:
        raise ValueError(
            f"Invalid {field} value {token!r} at line {line_number}."
        ) from error


class Occam1DModel(Occam1DBase):
    r"""Represent an Occam1D layered-earth parameterization.

    Parameters
    ----------
    depth : array-like of shape (n_layers,)
        Layer-top depths in metres. The first value must be zero and all later
        values must increase strictly. The final top begins the half-space.
    resistivity : array-like of shape (n_layers,)
        Layer resistivities in ohm metres. Finite values must be positive.
        ``nan`` may represent an unset starting value in memory, but
        :meth:`write` requires a fully parameterized model.
    penalty : array-like of shape (n_layers,), optional
        Non-negative native regularization penalty values. Defaults to zero.
    preference : array-like of shape (n_layers,), optional
        Finite native preferred model values. Defaults to zero.
    weight : array-like of shape (n_layers,), optional
        Non-negative native preference weights. Defaults to zero.
    verbose, logger, path, metadata, stream
        Shared options inherited from :class:`Occam1DBase`.

    Attributes
    ----------
    depth, resistivity, penalty, preference, weight : ndarray
        Independent one-dimensional arrays aligned by layer.

    Raises
    ------
    TypeError
        If a supplied array is not real-valued.
    ValueError
        If geometry, shape, resistivity, or regularization controls violate
        the model invariants.

    Notes
    -----
    There are ``n_layers - 1`` finite-thickness layers. The final layer is a
    half-space and therefore has infinite thickness. ``depth_max`` denotes
    its top, not a finite bottom of the earth model.

    Examples
    --------
    >>> model = Occam1DModel.build(4, 5.0, 500.0, resistivity=100.0)
    >>> model.n_layers
    4
    >>> model.thickness[-1]
    inf
    >>> model.is_parameterized
    True
    """

    @validate_params(MODEL_SCHEMA)
    def __init__(
        self,
        depth,
        resistivity,
        *,
        penalty=None,
        preference=None,
        weight=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.depth = self._vector(depth, "depth")
        self.resistivity = self._vector(resistivity, "resistivity")
        size = self.depth.size
        self.penalty = self._control(penalty, "penalty", size)
        self.preference = self._control(
            preference, "preference", size
        )
        self.weight = self._control(weight, "weight", size)
        self._validate_state()

    @staticmethod
    def _vector(values, name: str) -> np.ndarray:
        """Return a defensive one-dimensional floating-point array."""
        try:
            result = np.array(values, dtype=float, copy=True)
        except (TypeError, ValueError) as error:
            raise TypeError(f"{name} must contain real numbers.") from error
        if result.ndim != 1:
            raise ValueError(f"{name} must be one-dimensional.")
        return result

    @classmethod
    def _control(cls, values, name: str, size: int) -> np.ndarray:
        """Normalize an optional native control vector."""
        if values is None:
            return np.zeros(size, dtype=float)
        result = cls._vector(values, name)
        if result.size != size:
            raise ValueError(
                f"{name} must contain exactly {size} layer values."
            )
        return result

    def _validate_state(self) -> None:
        """Validate vertical geometry and all aligned layer properties."""
        arrays = {
            "depth": self.depth,
            "resistivity": self.resistivity,
            "penalty": self.penalty,
            "preference": self.preference,
            "weight": self.weight,
        }
        if any(values.ndim != 1 for values in arrays.values()):
            raise ValueError("Occam1D model arrays must be one-dimensional.")
        if len({values.shape for values in arrays.values()}) != 1:
            raise ValueError("All Occam1D model arrays must have equal shape.")
        if self.depth.size < 2:
            raise ValueError("An Occam1D model requires at least two layers.")
        if np.any(~np.isfinite(self.depth)):
            raise ValueError("depth must contain only finite values.")
        tolerance = np.finfo(float).eps * max(1.0, abs(self.depth[-1]))
        if not math.isclose(self.depth[0], 0.0, abs_tol=tolerance):
            raise ValueError("depth must start at zero metres.")
        self.depth[0] = 0.0
        if np.any(np.diff(self.depth) <= 0):
            raise ValueError("depth must increase strictly.")

        finite_rho = np.isfinite(self.resistivity)
        if np.any(np.isinf(self.resistivity)):
            raise ValueError("resistivity cannot contain infinite values.")
        if np.any(self.resistivity[finite_rho] <= 0):
            raise ValueError("Finite resistivity values must be positive.")
        for name in ("penalty", "preference", "weight"):
            values = getattr(self, name)
            if np.any(~np.isfinite(values)):
                raise ValueError(f"{name} must contain finite values.")
        if np.any(self.penalty < 0):
            raise ValueError("penalty values must be non-negative.")
        if np.any(self.weight < 0):
            raise ValueError("weight values must be non-negative.")

    @property
    def n_layers(self) -> int:
        """Number of layers, including the bottom half-space."""
        return int(self.depth.size)

    @property
    def n_finite_layers(self) -> int:
        """Number of finite-thickness layers above the half-space."""
        return self.n_layers - 1

    @property
    def thickness(self) -> np.ndarray:
        """Layer thicknesses in metres; the half-space value is infinite."""
        return np.r_[np.diff(self.depth), np.inf]

    @property
    def finite_thickness(self) -> np.ndarray:
        """Thicknesses in metres excluding the bottom half-space."""
        return np.diff(self.depth)

    @property
    def depth_max(self) -> float:
        """Top depth of the half-space in metres."""
        return float(self.depth[-1])

    @property
    def is_parameterized(self) -> bool:
        """Whether every layer has a finite physical resistivity."""
        return bool(np.all(np.isfinite(self.resistivity)))

    @property
    def n_parameters(self) -> int:
        """Number of finite layer-resistivity parameters."""
        return int(np.count_nonzero(np.isfinite(self.resistivity)))

    @property
    def resistivity_bounds(self) -> tuple[float, float] | None:
        """Minimum and maximum finite resistivity in ohm metres."""
        values = self.resistivity[np.isfinite(self.resistivity)]
        if not values.size:
            return None
        return float(np.min(values)), float(np.max(values))

    @property
    def conductance(self) -> np.ndarray:
        r"""Finite-layer longitudinal conductance in siemens.

        Each value is :math:`h_i / \rho_i`. Unset resistivity produces
        ``nan``. The half-space is excluded because its thickness is infinite.
        """
        return np.divide(
            self.finite_thickness,
            self.resistivity[:-1],
            out=np.full(self.n_finite_layers, np.nan),
            where=np.isfinite(self.resistivity[:-1]),
        )

    def with_resistivity(self, resistivity) -> Occam1DModel:
        """Return an independent model with replacement resistivity values."""
        return type(self)(
            self.depth,
            resistivity,
            penalty=self.penalty,
            preference=self.preference,
            weight=self.weight,
            metadata=self.metadata,
            verbose=self.verbose,
        )

    def filled(self, value: float) -> Occam1DModel:
        """Return a copy with unset resistivities replaced by ``value``.

        Parameters
        ----------
        value : float
            Finite, strictly positive replacement in ohm metres.
        """
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError("value must be a real number in ohm metres.")
        if not math.isfinite(value) or value <= 0:
            raise ValueError("value must be finite and strictly positive.")
        resistivity = np.where(
            np.isfinite(self.resistivity), self.resistivity, float(value)
        )
        return self.with_resistivity(resistivity)

    def to_records(self) -> list[dict[str, float]]:
        """Return one JSON-friendly physical/native record per layer."""
        names = (
            "depth",
            "resistivity",
            "penalty",
            "preference",
            "weight",
        )
        return [
            {name: float(getattr(self, name)[index]) for name in names}
            for index in range(self.n_layers)
        ]

    @classmethod
    @validate_params(MODEL_BUILD_SCHEMA)
    def build(
        cls,
        n_layers,
        first_thickness,
        depth_max,
        resistivity=100.0,
        **kwargs,
    ):
        r"""Construct a geometric layered mesh with uniform resistivity.

        ``n_layers - 1`` positive thicknesses form a geometric progression.
        Its first value is exactly ``first_thickness`` and its sum is exactly
        ``depth_max`` within floating-point precision. The unique growth
        factor is solved by bounded bisection; a factor of one gives uniform
        thickness.

        Returns
        -------
        Occam1DModel
            Fully parameterized model whose final depth is the half-space top.
        """
        finite_count = n_layers - 1
        minimum_depth = finite_count * first_thickness
        tolerance = np.finfo(float).eps * max(1.0, depth_max) * 16
        if finite_count == 1 and not math.isclose(
            depth_max, first_thickness, abs_tol=tolerance
        ):
            raise ValueError(
                "For a two-layer model, depth_max must equal "
                "first_thickness."
            )
        if depth_max + tolerance < minimum_depth:
            raise ValueError(
                "depth_max is too small for n_layers and first_thickness."
            )
        if math.isclose(depth_max, minimum_depth, abs_tol=tolerance):
            increments = np.full(finite_count, first_thickness)
        else:
            growth = cls._solve_growth(
                finite_count, first_thickness, depth_max
            )
            powers = np.arange(finite_count, dtype=float)
            increments = first_thickness * np.power(growth, powers)
            increments[-1] += depth_max - float(np.sum(increments))
        if np.any(~np.isfinite(increments)) or np.any(increments <= 0):
            raise ValueError(
                "The requested geometry exceeds floating-point limits."
            )
        depth = np.r_[0.0, np.cumsum(increments)]
        depth[-1] = depth_max
        return cls(
            depth,
            np.full(n_layers, resistivity, dtype=float),
            **kwargs,
        )

    @staticmethod
    def _solve_growth(count: int, first: float, total: float) -> float:
        """Solve the positive geometric-series growth factor."""

        def thickness_sum(growth: float) -> float:
            with np.errstate(over="ignore", invalid="ignore"):
                return float(
                    first
                    * np.sum(np.power(growth, np.arange(count, dtype=float)))
                )

        lower, upper = 1.0, 2.0
        for _ in range(64):
            if thickness_sum(upper) >= total:
                break
            upper *= 2.0
        else:
            raise ValueError("Could not bracket a finite layer-growth factor.")
        for _ in range(100):
            middle = (lower + upper) / 2.0
            if thickness_sum(middle) < total:
                lower = middle
            else:
                upper = middle
        return (lower + upper) / 2.0

    @validate_params(PATH_SCHEMA)
    def write(self, path) -> Path:
        """Write a complete ``Resistivity1DMod_1.0`` native model.

        Raises
        ------
        ValueError
            If any resistivity remains unset. Use :meth:`filled` to explicitly
            supply a physical starting value before serialization.
        """
        self.validate()
        if not self.is_parameterized:
            raise ValueError(
                "Cannot write a model with unset resistivity; call filled() "
                "with a physical starting value first."
            )
        target = self._prepare_output_file(path)
        with target.open("w", encoding="utf8", newline="\n") as stream:
            stream.write(f"Format: {_FORMAT}\n")
            stream.write(f"# Layers: {self.n_layers}\n")
            stream.write(
                "! Top_Depth Resistivity Penalty Preference Weight\n"
            )
            for record in zip(
                self.depth,
                self.resistivity,
                self.penalty,
                self.preference,
                self.weight,
            ):
                stream.write(" ".join(f"{value:.12g}" for value in record))
                stream.write("\n")
        self._bind_path(target)
        self.logger.debug(
            "Wrote %d-layer Occam1D model to %s.", self.n_layers, target
        )
        return target

    @classmethod
    @validate_params(PATH_SCHEMA)
    def read(cls, path, **kwargs):
        """Read a ``Resistivity1DMod_1.x`` native model.

        The reader accepts legacy two-column rows and fills their absent
        control columns with zero. Rows containing three or four columns are
        likewise completed. Extra columns are rejected to expose ambiguity.
        """
        probe = cls.__new__(cls)
        Occam1DBase.__init__(probe, **kwargs)
        source = probe._require_input_file(path, purpose="Occam1D model")
        lines = source.read_text(
            encoding="utf8", errors="replace"
        ).splitlines()
        if not any(
            _FORMAT_TOKEN in line.lower() for line in lines[:5]
        ):
            raise ValueError(f"Not an Occam1D model file: {source}")
        marker = None
        count = None
        for index, line in enumerate(lines):
            match = _LAYERS.match(line)
            if match:
                marker = index
                try:
                    count = int(match.group(1))
                except ValueError as error:
                    raise ValueError(
                        f"Invalid # Layers count in {source}."
                    ) from error
                break
        if marker is None or count is None:
            raise ValueError(f"Missing # Layers header in {source}.")
        if count < 2:
            raise ValueError("# Layers must declare at least two layers.")

        rows = []
        for index, line in enumerate(lines[marker + 1 :], marker + 2):
            text = line.strip()
            if not text or text.startswith(("!", "#")):
                continue
            fields = text.split()
            if not 2 <= len(fields) <= 5:
                raise ValueError(
                    f"Layer row at {source}:{index} must have two to five "
                    "numeric columns."
                )
            values = [
                _native_float(
                    field,
                    line_number=index,
                    field=f"layer column {column}",
                )
                for column, field in enumerate(fields, 1)
            ]
            values.extend([0.0] * (5 - len(values)))
            rows.append(values)
        if len(rows) != count:
            raise ValueError(
                f"Declared # Layers count is {count}, but {len(rows)} "
                f"layer rows were found in {source}."
            )
        columns = np.asarray(rows, dtype=float).T
        obj = cls(
            columns[0],
            columns[1],
            penalty=columns[2],
            preference=columns[3],
            weight=columns[4],
            **kwargs,
        )
        obj._bind_path(source)
        obj.logger.debug(
            "Read %d-layer Occam1D model from %s.", obj.n_layers, source
        )
        return obj

    def summary(self) -> str:
        """Return a concise human-readable model summary."""
        bounds = self.resistivity_bounds
        resistivity = (
            "unset"
            if bounds is None
            else f"{bounds[0]:.6g}..{bounds[1]:.6g} ohm m"
        )
        return (
            f"Occam1DModel(layers={self.n_layers}, "
            f"halfspace_top={self.depth_max:.6g} m, "
            f"resistivity={resistivity}, "
            f"parameterized={self.is_parameterized})"
        )

    def diagnostics(self) -> dict[str, Any]:
        """Extend lifecycle diagnostics with geometry and parameter ranges."""
        values = super().diagnostics()
        values.update(
            {
                "n_layers": self.n_layers,
                "n_finite_layers": self.n_finite_layers,
                "n_parameters": self.n_parameters,
                "depth_max_m": self.depth_max,
                "first_thickness_m": float(self.finite_thickness[0]),
                "resistivity_bounds_ohm_m": self.resistivity_bounds,
                "is_parameterized": self.is_parameterized,
            }
        )
        return values
