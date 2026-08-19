# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Strict ``OCCAMITER_FLEX`` startup and iteration-file support.

Startup files and solver-produced iteration files share the same native
format. This module normalizes their control fields, preserves additional
solver diagnostics, and keeps model parameters in their native
:math:`\log_{10}(\rho)` representation.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from numbers import Integral, Real
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from ...compat.sklearn import validate_params
from .base import Occam1DBase
from .config import Occam1DConfig
from .model import Occam1DModel
from .schema import PATH_SCHEMA, STARTUP_SCHEMA

__all__ = ["Occam1DStartup"]

_FORMAT = "OCCAMITER_FLEX"
_FORTRAN_EXPONENT = re.compile(r"(?<=\d)[Dd](?=[-+]?\d)")
_ALIASES = {
    "format": "format",
    "model file": "model_file",
    "data file": "data_file",
    "iteration": "iteration",
    "target misfit": "target_misfit",
    "max iterations": "max_iterations",
    "max iter": "max_iterations",
    "iterations to run": "max_iterations",
    "roughness type": "roughness_type",
    "stepsize cut count": "stepsize_cut_count",
    "step size cut count": "stepsize_cut_count",
    "debug level": "debug_level",
    "lagrange value": "lagrange_start",
    "lagrange start": "lagrange_start",
    "param count": "parameter_count",
}
_REQUIRED = frozenset(
    {
        "format",
        "model_file",
        "data_file",
        "iteration",
        "target_misfit",
        "max_iterations",
        "roughness_type",
        "stepsize_cut_count",
        "debug_level",
        "lagrange_start",
        "parameter_count",
    }
)
_INTEGER_FIELDS = frozenset(
    {
        "iteration",
        "max_iterations",
        "roughness_type",
        "stepsize_cut_count",
        "debug_level",
        "parameter_count",
    }
)
_FLOAT_FIELDS = frozenset({"target_misfit", "lagrange_start"})


def _native_float(token: str, *, field: str, line_number: int) -> float:
    """Parse a native float, including Fortran ``D`` exponents."""
    try:
        return float(_FORTRAN_EXPONENT.sub("E", token))
    except ValueError as error:
        raise ValueError(
            f"Invalid {field} value {token!r} at line {line_number}."
        ) from error


class Occam1DStartup(Occam1DBase):
    r"""Represent native Occam controls and log-resistivity parameters.

    Parameters
    ----------
    model_file, data_file : path-like
        Non-empty native file references passed to the executable. Newlines
        are forbidden; relative and absolute references are otherwise kept.
    parameters : array-like of shape (n_parameters,)
        Finite :math:`\log_{10}` resistivity values. At least two parameters
        are required because an Occam1D model includes a half-space.
    iteration : int, default=0
        Non-negative iteration number. Zero denotes initial startup controls.
    target_misfit : float, default=1.0
        Finite, strictly positive normalized RMS target.
    max_iterations : int, default=30
        Positive number of solver iterations to run.
    roughness_type : {1, 2}, default=1
        Native first- or second-difference roughness operator.
    stepsize_cut_count : int, default=8
        Non-negative maximum line-search step cuts.
    debug_level : int, default=1
        Non-negative native diagnostic level.
    lagrange_start : float, default=5.0
        Finite, strictly positive starting/current Lagrange multiplier.
    extra_fields : mapping, optional
        Additional native headers, such as ``Roughness Value``, ``Misfit
        Value``, ``Description``, or ``Date/Time``. They are preserved during
        a read/write round trip.
    verbose, logger, path, metadata, stream
        Shared options inherited from :class:`Occam1DBase`.

    Raises
    ------
    TypeError
        If controls or parameters use unsupported types.
    ValueError
        If any control is outside its scientific/native domain.

    Notes
    -----
    ``parameters`` are not physical ohm metres. Use
    :attr:`physical_resistivity` or :meth:`apply_to_model` for physical model
    values. Exponentiation overflow raises rather than returning infinity.

    Examples
    --------
    >>> startup = Occam1DStartup(
    ...     "Model", "Data", [2.0, 2.0, 2.5], target_misfit=1.0,
    ... )
    >>> startup.is_initial, startup.n_parameters
    (True, 3)
    >>> startup.physical_resistivity.tolist()
    [100.0, 100.0, 316.22776601683796]
    """

    @validate_params(STARTUP_SCHEMA)
    def __init__(
        self,
        model_file,
        data_file,
        parameters,
        iteration=0,
        target_misfit=1.0,
        max_iterations=30,
        roughness_type=1,
        stepsize_cut_count=8,
        debug_level=1,
        lagrange_start=5.0,
        extra_fields=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.model_file = self._file_reference("model_file", model_file)
        self.data_file = self._file_reference("data_file", data_file)
        self.parameters = self._parameter_vector(parameters)
        self.iteration = self._integer("iteration", iteration, minimum=0)
        self.target_misfit = self._positive_real(
            "target_misfit", target_misfit
        )
        self.max_iterations = self._integer(
            "max_iterations", max_iterations, minimum=1
        )
        self.roughness_type = self._integer(
            "roughness_type", roughness_type, minimum=1
        )
        self.stepsize_cut_count = self._integer(
            "stepsize_cut_count", stepsize_cut_count, minimum=0
        )
        self.debug_level = self._integer(
            "debug_level", debug_level, minimum=0
        )
        self.lagrange_start = self._positive_real(
            "lagrange_start", lagrange_start
        )
        self.extra_fields = self._extras(extra_fields)
        self._validate_state()

    @staticmethod
    def _file_reference(name: str, value) -> str:
        text = str(value).strip()
        if not text:
            raise ValueError(f"{name} cannot be empty.")
        if "\n" in text or "\r" in text:
            raise ValueError(f"{name} cannot contain a newline.")
        return text

    @staticmethod
    def _parameter_vector(values) -> np.ndarray:
        try:
            result = np.array(values, dtype=float, copy=True)
        except (TypeError, ValueError) as error:
            raise TypeError("parameters must contain real numbers.") from error
        if result.ndim != 1:
            raise ValueError("parameters must be one-dimensional.")
        return result

    @staticmethod
    def _integer(name: str, value, *, minimum: int) -> int:
        if isinstance(value, bool) or not isinstance(value, Integral):
            raise TypeError(f"{name} must be an integer.")
        result = int(value)
        if result < minimum:
            raise ValueError(f"{name} must be at least {minimum}.")
        return result

    @staticmethod
    def _positive_real(name: str, value) -> float:
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError(f"{name} must be a real number.")
        result = float(value)
        if not math.isfinite(result) or result <= 0:
            raise ValueError(f"{name} must be finite and strictly positive.")
        return result

    @staticmethod
    def _extras(values) -> dict[str, str]:
        if values is None:
            return {}
        if not isinstance(values, Mapping):
            raise TypeError("extra_fields must be a mapping or None.")
        result = {}
        reserved = set(_ALIASES)
        for key, value in values.items():
            name = str(key).strip()
            text = str(value).strip()
            if not name or ":" in name or "\n" in name or "\r" in name:
                raise ValueError(f"Invalid extra field name: {key!r}.")
            if name.casefold() in reserved:
                raise ValueError(
                    f"extra_fields cannot override native field {name!r}."
                )
            if "\n" in text or "\r" in text:
                raise ValueError(
                    f"Extra field {name!r} cannot contain a newline."
                )
            result[name] = text
        return result

    def _validate_state(self) -> None:
        """Validate parameter and native-control invariants after mutation."""
        if self.parameters.ndim != 1 or self.parameters.size < 2:
            raise ValueError(
                "parameters must contain at least two log10 model values."
            )
        if np.any(~np.isfinite(self.parameters)):
            raise ValueError("parameters must contain finite log10 values.")
        self._file_reference("model_file", self.model_file)
        self._file_reference("data_file", self.data_file)
        self._integer("iteration", self.iteration, minimum=0)
        self._positive_real("target_misfit", self.target_misfit)
        self._integer("max_iterations", self.max_iterations, minimum=1)
        roughness = self._integer(
            "roughness_type", self.roughness_type, minimum=1
        )
        if roughness not in {1, 2}:
            raise ValueError("roughness_type must be either 1 or 2.")
        self._integer(
            "stepsize_cut_count", self.stepsize_cut_count, minimum=0
        )
        self._integer("debug_level", self.debug_level, minimum=0)
        self._positive_real("lagrange_start", self.lagrange_start)
        self._extras(self.extra_fields)

    @property
    def n_parameters(self) -> int:
        """Number of log-resistivity model parameters."""
        return int(self.parameters.size)

    @property
    def is_initial(self) -> bool:
        """Whether this object represents iteration-zero startup controls."""
        return self.iteration == 0

    @property
    def is_iteration(self) -> bool:
        """Whether this object represents a solver-produced later iteration."""
        return self.iteration > 0

    @property
    def extra_fields_view(self):
        """Read-only live view of additional native header fields."""
        return MappingProxyType(self.extra_fields)

    @property
    def physical_resistivity(self) -> np.ndarray:
        """Model resistivities in ohm metres."""
        with np.errstate(over="ignore", invalid="ignore"):
            values = np.power(10.0, self.parameters)
        if np.any(~np.isfinite(values)) or np.any(values <= 0):
            raise ValueError(
                "Log10 model parameters exceed physical resistivity range."
            )
        return values

    @property
    def resistivity_bounds(self) -> tuple[float, float]:
        """Minimum and maximum physical resistivity in ohm metres."""
        values = self.physical_resistivity
        return float(np.min(values)), float(np.max(values))

    def apply_to_model(self, model) -> Occam1DModel:
        """Return a layered model populated by this parameter vector.

        Raises
        ------
        TypeError
            If ``model`` is not an :class:`Occam1DModel`.
        ValueError
            If the layer and parameter counts differ.
        """
        if not isinstance(model, Occam1DModel):
            raise TypeError("model must be an Occam1DModel instance.")
        if model.n_layers != self.n_parameters:
            raise ValueError(
                "Startup parameter count does not match model layer count."
            )
        return model.with_resistivity(self.physical_resistivity)

    @classmethod
    def from_model(cls, model, config=None, **kwargs):
        """Create iteration-zero controls from a layered model.

        Unset model resistivities are explicitly filled using
        ``config.starting_resistivity`` before conversion to log10 space.
        """
        if not isinstance(model, Occam1DModel):
            raise TypeError("model must be an Occam1DModel instance.")
        if config is not None and not isinstance(config, Occam1DConfig):
            raise TypeError("config must be an Occam1DConfig or None.")
        cfg = config or Occam1DConfig()
        cfg.validate()
        rho = np.where(
            np.isfinite(model.resistivity),
            model.resistivity,
            cfg.starting_resistivity,
        )
        metadata = dict(kwargs.pop("metadata", {}) or {})
        metadata.update(
            {
                "source_model_layers": model.n_layers,
                "filled_resistivity_count": int(
                    np.count_nonzero(~np.isfinite(model.resistivity))
                ),
            }
        )
        return cls(
            cfg.model_file,
            cfg.data_file,
            np.log10(rho),
            target_misfit=cfg.target_misfit,
            max_iterations=cfg.max_iterations,
            roughness_type=cfg.roughness_type,
            stepsize_cut_count=cfg.stepsize_cut_count,
            debug_level=cfg.debug_level,
            lagrange_start=cfg.lagrange_start,
            metadata=metadata,
            **kwargs,
        )

    @validate_params(PATH_SCHEMA)
    def write(self, path) -> Path:
        """Write a solver-canonical ``OCCAMITER_FLEX`` file.

        The parameter block begins immediately after ``Param Count``, matching
        the native executable. Four log10 parameters are written per line.
        """
        self.validate()
        target = self._prepare_output_file(path)
        fields = [
            ("Format", _FORMAT),
            ("Model File", self.model_file),
            ("Data File", self.data_file),
            ("Iterations to run", self.max_iterations),
            ("Target Misfit", self.target_misfit),
            ("Roughness Type", self.roughness_type),
            ("Stepsize Cut Count", self.stepsize_cut_count),
            ("Debug Level", self.debug_level),
            ("Iteration", self.iteration),
            ("Lagrange Value", self.lagrange_start),
        ]
        fields.extend(self.extra_fields.items())
        fields.append(("Param Count", self.n_parameters))
        with target.open("w", encoding="utf8", newline="\n") as stream:
            for key, value in fields:
                stream.write(f"{key}: {value}\n")
            for start in range(0, self.n_parameters, 4):
                values = self.parameters[start : start + 4]
                stream.write(" ".join(f"{value:.12g}" for value in values))
                stream.write("\n")
        self._bind_path(target)
        self.logger.debug(
            "Wrote Occam1D iteration %d with %d parameters to %s.",
            self.iteration,
            self.n_parameters,
            target,
        )
        return target

    @classmethod
    @validate_params(PATH_SCHEMA)
    def read(cls, path, **kwargs):
        """Read a startup or solver iteration ``OCCAMITER_FLEX`` file.

        Field names are case-insensitive. Both ``Iterations to run`` and the
        legacy ``Max Iterations`` spelling are accepted. Parameter blocks may
        start immediately after ``Param Count`` or after a legacy
        ``Model Values:`` marker. Unknown headers are retained in
        :attr:`extra_fields`.
        """
        probe = cls.__new__(cls)
        Occam1DBase.__init__(probe, **kwargs)
        source = probe._require_input_file(
            path, purpose="OCCAMITER_FLEX startup/iteration"
        )
        lines = source.read_text(
            encoding="utf8", errors="replace"
        ).splitlines()
        parsed = {}
        extras = {}
        parameter_tokens = []
        reading_parameters = False
        for line_number, raw in enumerate(lines, start=1):
            text = raw.strip()
            if not text or text.startswith(("!", "#")):
                continue
            if reading_parameters:
                if text.casefold() == "model values:":
                    continue
                parameter_tokens.extend(
                    (token, line_number) for token in text.split()
                )
                continue
            if ":" not in text:
                raise ValueError(
                    f"Malformed OCCAMITER_FLEX header at "
                    f"{source}:{line_number}."
                )
            raw_key, raw_value = text.split(":", 1)
            key = raw_key.strip()
            value = raw_value.strip()
            canonical = _ALIASES.get(key.casefold())
            if canonical is None:
                if key in extras:
                    raise ValueError(
                        f"Duplicate extra field {key!r} in {source}."
                    )
                extras[key] = value
                continue
            if canonical in parsed:
                raise ValueError(
                    f"Duplicate native field {key!r} in {source}."
                )
            parsed[canonical] = cls._parse_header_value(
                canonical, value, line_number
            )
            if canonical == "parameter_count":
                reading_parameters = True

        missing = sorted(_REQUIRED.difference(parsed))
        if missing:
            names = ", ".join(missing)
            raise ValueError(
                f"Missing required field(s) in {source}: {names}."
            )
        if str(parsed["format"]).upper() != _FORMAT:
            raise ValueError(
                f"Expected format {_FORMAT!r}, got {parsed['format']!r}."
            )
        count = parsed.pop("parameter_count")
        if count < 2:
            raise ValueError("Param Count must declare at least two values.")
        parameters = np.asarray(
            [
                _native_float(
                    token,
                    field="model parameter",
                    line_number=line_number,
                )
                for token, line_number in parameter_tokens
            ],
            dtype=float,
        )
        if parameters.size != count:
            raise ValueError(
                f"Declared Param Count is {count}, but {parameters.size} "
                f"model values were found in {source}."
            )
        parsed.pop("format")
        obj = cls(
            parameters=parameters,
            extra_fields=extras,
            **parsed,
            **kwargs,
        )
        obj._bind_path(source)
        obj.logger.debug(
            "Read Occam1D iteration %d with %d parameters from %s.",
            obj.iteration,
            obj.n_parameters,
            source,
        )
        return obj

    @staticmethod
    def _parse_header_value(field, value, line_number):
        """Coerce one recognized native header value."""
        if field in {"format", "model_file", "data_file"}:
            return value
        if field in _INTEGER_FIELDS:
            number = _native_float(
                value, field=field, line_number=line_number
            )
            if not math.isfinite(number) or number != math.floor(number):
                raise ValueError(
                    f"{field} must be an integer at line {line_number}."
                )
            return int(number)
        if field in _FLOAT_FIELDS:
            return _native_float(
                value, field=field, line_number=line_number
            )
        return value

    def summary(self) -> str:
        """Return a concise human-readable startup/iteration summary."""
        low, high = self.resistivity_bounds
        kind = "startup" if self.is_initial else "iteration"
        return (
            f"Occam1DStartup(kind={kind!r}, iteration={self.iteration}, "
            f"parameters={self.n_parameters}, target={self.target_misfit:g}, "
            f"resistivity={low:.6g}..{high:.6g} ohm m)"
        )

    def diagnostics(self) -> dict[str, Any]:
        """Extend lifecycle diagnostics with inversion-control statistics."""
        values = super().diagnostics()
        low, high = self.resistivity_bounds
        values.update(
            {
                "model_file": self.model_file,
                "data_file": self.data_file,
                "iteration": self.iteration,
                "is_initial": self.is_initial,
                "n_parameters": self.n_parameters,
                "target_misfit": self.target_misfit,
                "max_iterations": self.max_iterations,
                "roughness_type": self.roughness_type,
                "stepsize_cut_count": self.stepsize_cut_count,
                "debug_level": self.debug_level,
                "lagrange_start": self.lagrange_start,
                "resistivity_min_ohm_m": low,
                "resistivity_max_ohm_m": high,
                "extra_fields": dict(self.extra_fields),
            }
        )
        return values
