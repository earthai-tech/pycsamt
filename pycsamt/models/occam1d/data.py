# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Scientific container and native I/O for one Occam1D sounding.

Occam1D represents apparent resistivity in :math:`\log_{10}` space while
pyCSAMT exposes it in physical ohm metres. This module owns that conversion,
preserves unavailable observations as ``nan``, and validates observation and
uncertainty alignment before native files reach an inversion.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import numpy as np

from ...compat.sklearn import validate_params
from .base import Occam1DBase
from .schema import DATA_SCHEMA, MODE_OPTIONS, PATH_SCHEMA

__all__ = ["Occam1DData"]

_FORMAT = "EMData_1.1"
_RHO_CODE = 103
_PHASE_CODE = 104
_LOG10_SCALE = math.log(10.0)
_FORTRAN_EXPONENT = re.compile(r"(?<=\d)[Dd](?=[-+]?\d)")
_HEADER = re.compile(r"^\s*#\s*([^:]+)\s*:\s*(.*?)\s*$")


def _native_float(token: str, *, line_number: int, field: str) -> float:
    """Parse a native decimal, including a Fortran ``D`` exponent."""
    try:
        return float(_FORTRAN_EXPONENT.sub("E", token))
    except ValueError as error:
        raise ValueError(
            f"Invalid {field} value {token!r} at line {line_number}."
        ) from error


class Occam1DData(Occam1DBase):
    r"""Represent one validated MT resistivity/phase sounding.

    Parameters
    ----------
    frequency : array-like of shape (n_frequencies,)
        Unique, strictly positive, finite frequencies in hertz. Ascending and
        descending order are both accepted.
    resistivity : array-like of shape (n_frequencies,)
        Apparent resistivity in ohm metres. Finite values must be positive;
        unavailable observations must be ``nan``.
    phase : array-like of shape (n_frequencies,)
        Impedance phase in degrees. Unavailable observations must be ``nan``.
    resistivity_error : array-like of shape (n_frequencies,)
        Positive relative resistivity uncertainty as a fraction. Thus,
        ``0.05`` means five percent.
    phase_error : array-like of shape (n_frequencies,)
        Positive absolute phase uncertainty in degrees.
    mode : {"xy", "yx", "te", "tm", "determinant"}, default="determinant"
        Tensor component or invariant represented by the sounding.
    station : str, default="site"
        Non-empty station identifier stored as native metadata.
    verbose, logger, path, metadata, stream
        Shared options inherited from :class:`Occam1DBase`.

    Raises
    ------
    TypeError
        If an input cannot be represented as a numerical vector.
    ValueError
        If arrays are empty, misaligned, or scientifically invalid.

    Notes
    -----
    Missing data and physical zero are distinct. Missing observations remain
    ``nan`` and are omitted from the native file. An uncertainty supplied
    without its observation is normalized to ``nan`` and recorded as a
    warning. Every finite observation requires a finite positive uncertainty.

    Examples
    --------
    >>> data = Occam1DData(
    ...     [100.0, 10.0], [50.0, 80.0], [42.0, 45.0],
    ...     [0.05, 0.05], [1.5, 1.5], station="S01",
    ... )
    >>> data.n_frequencies, data.n_data
    (2, 4)
    >>> data.period.tolist()
    [0.01, 0.1]
    """

    @validate_params(DATA_SCHEMA)
    def __init__(
        self,
        frequency,
        resistivity,
        phase,
        resistivity_error,
        phase_error,
        mode="determinant",
        station="site",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.frequency = self._vector(frequency, "frequency")
        self.resistivity = self._vector(resistivity, "resistivity")
        self.phase = self._vector(phase, "phase")
        self.resistivity_error = self._vector(
            resistivity_error, "resistivity_error"
        )
        self.phase_error = self._vector(phase_error, "phase_error")
        self.mode = mode.lower()
        self.station = station.strip()
        if not self.station:
            raise ValueError("station must be a non-empty string.")
        self._normalize_orphan_errors()
        self._validate_state()

    @staticmethod
    def _vector(values, name: str) -> np.ndarray:
        """Return an independent one-dimensional floating-point array."""
        try:
            array = np.array(values, dtype=float, copy=True)
        except (TypeError, ValueError) as error:
            raise TypeError(
                f"{name} must contain real numeric values."
            ) from error
        if array.ndim != 1:
            raise ValueError(f"{name} must be one-dimensional.")
        return array

    def _normalize_orphan_errors(self) -> None:
        """Discard uncertainties that have no corresponding observation."""
        for observation_name, error_name in (
            ("resistivity", "resistivity_error"),
            ("phase", "phase_error"),
        ):
            observation = getattr(self, observation_name)
            errors = getattr(self, error_name)
            if observation.shape != errors.shape:
                continue
            orphan = ~np.isfinite(observation) & np.isfinite(errors)
            if np.any(orphan):
                errors[orphan] = np.nan
                self.add_warning(
                    f"Discarded {int(np.count_nonzero(orphan))} {error_name} "
                    f"value(s) without a finite {observation_name}."
                )

    def _validate_state(self) -> None:
        """Validate array alignment and scientific sounding invariants."""
        arrays = {
            "frequency": self.frequency,
            "resistivity": self.resistivity,
            "phase": self.phase,
            "resistivity_error": self.resistivity_error,
            "phase_error": self.phase_error,
        }
        if any(values.ndim != 1 for values in arrays.values()):
            raise ValueError("Occam1D data arrays must be one-dimensional.")
        if len({values.shape for values in arrays.values()}) != 1:
            raise ValueError("All Occam1D data arrays must have equal shape.")
        if not self.frequency.size:
            raise ValueError("At least one frequency is required.")
        if np.any(~np.isfinite(self.frequency)) or np.any(
            self.frequency <= 0
        ):
            raise ValueError("frequency must contain finite positive values.")
        if np.unique(self.frequency).size != self.frequency.size:
            raise ValueError("frequency values must be unique.")

        rho = np.isfinite(self.resistivity)
        phi = np.isfinite(self.phase)
        if np.any(self.resistivity[rho] <= 0):
            raise ValueError("Finite resistivity values must be positive.")
        if not np.any(rho | phi):
            raise ValueError("At least one finite observation is required.")
        for observed, name, errors in (
            (rho, "resistivity_error", self.resistivity_error),
            (phi, "phase_error", self.phase_error),
        ):
            if np.any(~np.isfinite(errors[observed])):
                raise ValueError(
                    f"Every finite observation requires a finite {name}."
                )
            if np.any(errors[observed] <= 0):
                raise ValueError(f"Finite {name} values must be positive.")
        if not isinstance(self.station, str) or not self.station.strip():
            raise ValueError("station must be a non-empty string.")
        if self.mode not in MODE_OPTIONS:
            choices = ", ".join(sorted(MODE_OPTIONS))
            raise ValueError(f"mode must be one of: {choices}.")

    @property
    def n_frequencies(self) -> int:
        """Number of unique sampled frequencies."""
        return int(self.frequency.size)

    @property
    def resistivity_mask(self) -> np.ndarray:
        """Boolean mask of usable apparent-resistivity observations."""
        return np.isfinite(self.resistivity)

    @property
    def phase_mask(self) -> np.ndarray:
        """Boolean mask of usable phase observations."""
        return np.isfinite(self.phase)

    @property
    def n_resistivity(self) -> int:
        """Number of usable apparent-resistivity observations."""
        return int(np.count_nonzero(self.resistivity_mask))

    @property
    def n_phase(self) -> int:
        """Number of usable phase observations."""
        return int(np.count_nonzero(self.phase_mask))

    @property
    def n_data(self) -> int:
        """Total number of scalar observations written for inversion."""
        return self.n_resistivity + self.n_phase

    @property
    def period(self) -> np.ndarray:
        """Periods in seconds, aligned with :attr:`frequency`."""
        return np.reciprocal(self.frequency)

    @property
    def frequency_bounds(self) -> tuple[float, float]:
        """Minimum and maximum frequency in hertz."""
        return float(np.min(self.frequency)), float(np.max(self.frequency))

    def to_records(self) -> list[dict[str, float]]:
        """Return one JSON-friendly physical-unit record per frequency."""
        names = (
            "frequency",
            "resistivity",
            "phase",
            "resistivity_error",
            "phase_error",
        )
        return [
            {name: float(getattr(self, name)[index]) for name in names}
            for index in range(self.n_frequencies)
        ]

    @validate_params(PATH_SCHEMA)
    def write(self, path) -> Path:
        r"""Write the sounding using the native ``EMData_1.1`` dialect.

        Resistivity is serialized as :math:`\log_{10}(\rho_a)` and relative
        uncertainty as ``error / ln(10)``. The object is bound only after the
        write succeeds.
        """
        self.validate()
        target = self._prepare_output_file(path)
        rows = []
        for offset in range(self.n_frequencies):
            index = offset + 1
            if self.resistivity_mask[offset]:
                rows.append(
                    (
                        _RHO_CODE,
                        index,
                        math.log10(self.resistivity[offset]),
                        self.resistivity_error[offset] / _LOG10_SCALE,
                    )
                )
            if self.phase_mask[offset]:
                rows.append(
                    (
                        _PHASE_CODE,
                        index,
                        self.phase[offset],
                        self.phase_error[offset],
                    )
                )
        with target.open("w", encoding="utf8", newline="\n") as stream:
            stream.write(f"Format: {_FORMAT}\n")
            stream.write(f"! Station: {self.station}; Mode: {self.mode}\n")
            stream.write(f"# Frequencies: {self.n_frequencies}\n")
            for value in self.frequency:
                stream.write(f"{value:.12g}\n")
            stream.write("# Receivers: 1\n")
            stream.write("0 0 0 0 0 0\n")
            stream.write(f"# Data: {len(rows)}\n")
            stream.write("! Type Freq# Tx# Rx# Data Std_Error\n")
            for code, index, value, error in rows:
                stream.write(
                    f"{code:d} {index:d} 0 1 "
                    f"{value:.12g} {error:.12g}\n"
                )
        self._bind_path(target)
        self.logger.debug(
            "Wrote %d observations for station %s to %s.",
            self.n_data,
            self.station,
            target,
        )
        return target

    @classmethod
    @validate_params(PATH_SCHEMA)
    def read(cls, path, **kwargs):
        r"""Read and validate an ``EMData_1.x`` native data file.

        Malformed rows, invalid indices, duplicate physical rows, and
        declared-count mismatches raise descriptive errors. Unknown type
        codes are ignored with an object warning for extension compatibility.
        """
        probe = cls.__new__(cls)
        Occam1DBase.__init__(probe, **kwargs)
        source = probe._require_input_file(path, purpose="Occam1D data")
        lines = source.read_text(
            encoding="utf8", errors="replace"
        ).splitlines()
        if not any("emdata_1." in line.lower() for line in lines[:5]):
            raise ValueError(f"Not an Occam1D data file: {source}")

        headers = cls._headers(lines)
        frequency_marker, count_text = cls._required_header(
            headers, "frequencies", source
        )
        count = cls._positive_count(count_text, "Frequencies", source)
        if frequency_marker + count >= len(lines):
            raise ValueError(
                f"Frequency block is truncated in Occam1D data file: {source}"
            )
        frequency = np.asarray(
            [
                _native_float(
                    lines[index].strip(),
                    line_number=index + 1,
                    field="frequency",
                )
                for index in range(
                    frequency_marker + 1, frequency_marker + count + 1
                )
            ],
            dtype=float,
        )
        data_marker, declared_text = cls._required_header(
            headers, "data", source
        )
        declared = cls._nonnegative_count(declared_text, "Data", source)
        parsed, unknown = cls._parse_rows(
            lines[data_marker + 1 :], count, data_marker + 2, source
        )
        if len(parsed) + unknown != declared:
            raise ValueError(
                f"Declared # Data count is {declared}, but "
                f"{len(parsed) + unknown} rows were found in {source}."
            )

        station, mode = cls._metadata(lines)
        values = {
            "resistivity": np.full(count, np.nan),
            "phase": np.full(count, np.nan),
            "resistivity_error": np.full(count, np.nan),
            "phase_error": np.full(count, np.nan),
        }
        for code, index, value, error in parsed:
            if code == _RHO_CODE:
                try:
                    resistivity = 10.0**value
                except OverflowError as exception:
                    raise ValueError(
                        "Log-resistivity overflows physical units at "
                        f"frequency index {index + 1} in {source}."
                    ) from exception
                if not math.isfinite(resistivity) or resistivity <= 0:
                    raise ValueError(
                        "Log-resistivity is outside representable physical "
                        f"units at frequency index {index + 1} in {source}."
                    )
                values["resistivity"][index] = resistivity
                values["resistivity_error"][index] = error * _LOG10_SCALE
            else:
                values["phase"][index] = value
                values["phase_error"][index] = error
        obj = cls(
            frequency,
            values["resistivity"],
            values["phase"],
            values["resistivity_error"],
            values["phase_error"],
            mode=mode,
            station=station,
            **kwargs,
        )
        obj._bind_path(source)
        if unknown:
            obj.add_warning(
                f"Ignored {unknown} native row(s) with unsupported data "
                "type codes."
            )
        obj.logger.debug(
            "Read %d observations for station %s from %s.",
            obj.n_data,
            obj.station,
            source,
        )
        return obj

    @staticmethod
    def _headers(lines: list[str]) -> dict[str, tuple[int, str]]:
        """Index native count headers case-insensitively."""
        headers = {}
        for index, line in enumerate(lines):
            match = _HEADER.match(line)
            if match:
                headers[match.group(1).strip().lower()] = (
                    index,
                    match.group(2).strip(),
                )
        return headers

    @staticmethod
    def _required_header(headers, name: str, source: Path):
        try:
            return headers[name]
        except KeyError as error:
            raise ValueError(
                f"Missing # {name.title()} header in {source}."
            ) from error

    @staticmethod
    def _nonnegative_count(text: str, name: str, source: Path) -> int:
        try:
            value = int(text)
        except ValueError as error:
            raise ValueError(
                f"Invalid # {name} count {text!r} in {source}."
            ) from error
        if value < 0:
            raise ValueError(f"# {name} count cannot be negative in {source}.")
        return value

    @classmethod
    def _positive_count(cls, text: str, name: str, source: Path) -> int:
        value = cls._nonnegative_count(text, name, source)
        if value == 0:
            raise ValueError(f"# {name} count must be positive in {source}.")
        return value

    @staticmethod
    def _metadata(lines: list[str]) -> tuple[str, str]:
        station, mode = "site", "determinant"
        for line in lines:
            if not line.lstrip().lower().startswith("! station:"):
                continue
            fields = line.lstrip()[1:].split(";")
            station = fields[0].split(":", 1)[1].strip() or "site"
            for field in fields[1:]:
                key, separator, value = field.partition(":")
                if separator and key.strip().lower() == "mode":
                    mode = value.strip().lower()
        return station, mode

    @staticmethod
    def _parse_rows(lines, count, first_line, source):
        parsed = []
        unknown = 0
        occupied = set()
        for offset, line in enumerate(lines):
            line_number = first_line + offset
            text = line.strip()
            if not text or text.startswith("!") or text.startswith("#"):
                continue
            fields = text.split()
            if len(fields) < 6:
                raise ValueError(
                    f"Malformed data row at {source}:{line_number}; "
                    "expected at least six fields."
                )
            try:
                code = int(fields[0])
                index = int(fields[1]) - 1
            except ValueError as error:
                raise ValueError(
                    f"Invalid type or frequency index at "
                    f"{source}:{line_number}."
                ) from error
            if not 0 <= index < count:
                raise ValueError(
                    f"Frequency index {index + 1} at {source}:{line_number} "
                    f"is outside 1..{count}."
                )
            value = _native_float(
                fields[4], line_number=line_number, field="data"
            )
            error = _native_float(
                fields[5], line_number=line_number, field="standard error"
            )
            if not math.isfinite(value) or not math.isfinite(error):
                raise ValueError(
                    f"Non-finite data or standard error at "
                    f"{source}:{line_number}."
                )
            if code not in {_RHO_CODE, _PHASE_CODE}:
                unknown += 1
                continue
            key = code, index
            if key in occupied:
                raise ValueError(
                    f"Duplicate type {code} at frequency index {index + 1} "
                    f"in {source}."
                )
            occupied.add(key)
            parsed.append((code, index, value, error))
        return parsed, unknown

    def summary(self) -> str:
        """Return a concise human-readable sounding summary."""
        low, high = self.frequency_bounds
        return (
            f"Occam1DData(station={self.station!r}, mode={self.mode!r}, "
            f"frequencies={self.n_frequencies}, data={self.n_data}, "
            f"range={low:.6g}..{high:.6g} Hz)"
        )

    def diagnostics(self) -> dict[str, Any]:
        """Extend lifecycle diagnostics with sounding coverage statistics."""
        values = super().diagnostics()
        low, high = self.frequency_bounds
        values.update(
            {
                "station": self.station,
                "mode": self.mode,
                "n_frequencies": self.n_frequencies,
                "n_resistivity": self.n_resistivity,
                "n_phase": self.n_phase,
                "n_data": self.n_data,
                "frequency_min_hz": low,
                "frequency_max_hz": high,
            }
        )
        return values
