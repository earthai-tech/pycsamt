# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Parse and query native Occam1D convergence histories.

Occam executables exist in several closely related text dialects. Some report
accepted RMS using ``AND IS =`` while others use ``RMS:``, ``R.M.S.``, or
``MISFIT``. Fortran ``D`` exponents are also common. This module normalizes
those representations into immutable iteration records without silently
replacing missing diagnostics by physical zeros.
"""

from __future__ import annotations

import csv
import math
import re
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from numbers import Real
from pathlib import Path
from typing import Any

import numpy as np

from ...compat.sklearn import validate_params
from .base import Occam1DBase
from .schema import (
    LOG_CONVERGENCE_SCHEMA,
    LOG_ITERATION_SCHEMA,
    LOG_READ_SCHEMA,
    LOG_WRITE_SCHEMA,
)

__all__ = ["Occam1DLog", "Occam1DLogRecord"]

_NAN = float("nan")
_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
_ITERATION = re.compile(r"\bITERATION\s*(?:=|:|#)?\s*(\d+)\b", re.I)

_METRIC_PATTERNS = {
    "rms": (
        re.compile(r"\bAND\s+IS\s*(?:=|:)?\s*(" + _NUMBER + r")", re.I),
        re.compile(r"\bR\.?M\.?S\.?\s*(?:=|:)?\s*(" + _NUMBER + r")", re.I),
        re.compile(r"\bMISFIT\s*(?:=|:)?\s*(" + _NUMBER + r")", re.I),
    ),
    "roughness": (
        re.compile(
            r"\bROUGHNESS(?:\s+IS)?\s*(?:=|:)?\s*(" + _NUMBER + r")",
            re.I,
        ),
    ),
    "lagrange": (
        re.compile(
            r"\bLAGRANGE(?:\s+VALUE)?\s*(?:=|:)?\s*(" + _NUMBER + r")",
            re.I,
        ),
        re.compile(r"\bMU\s*(?:=|:)?\s*(" + _NUMBER + r")", re.I),
    ),
    "stepsize": (
        re.compile(
            r"\bSTEP\s*SIZE(?:\s+IS)?\s*(?:=|:)?\s*(" + _NUMBER + r")",
            re.I,
        ),
    ),
}


def _native_float(value: str) -> float:
    """Convert a native decimal or Fortran exponent token."""
    return float(value.replace("D", "E").replace("d", "e"))


def _metric(line: str, name: str) -> float | None:
    """Return the first metric value recognized in ``line``."""
    for pattern in _METRIC_PATTERNS[name]:
        match = pattern.search(line)
        if match is not None:
            return _native_float(match.group(1))
    return None


@dataclass(frozen=True)
class Occam1DLogRecord:
    r"""Immutable convergence diagnostics for one Occam iteration.

    Parameters
    ----------
    iteration : int
        Non-negative native iteration number.
    rms : float, default=nan
        Accepted normalized root-mean-square data misfit.
    roughness : float, default=nan
        Model roughness penalty reported by the solver.
    lagrange : float, default=nan
        Accepted Lagrange multiplier :math:`\mu`.
    stepsize : float, default=nan
        Accepted line-search step size.

    Notes
    -----
    Missing metrics remain ``nan``. They are never converted to zero because
    zero has a distinct numerical meaning for roughness and step size.
    """

    iteration: int
    rms: float = _NAN
    roughness: float = _NAN
    lagrange: float = _NAN
    stepsize: float = _NAN

    def __post_init__(self) -> None:
        if isinstance(self.iteration, bool) or not isinstance(
            self.iteration, (int, np.integer)
        ):
            raise TypeError("iteration must be a non-negative integer.")
        if self.iteration < 0:
            raise ValueError("iteration must be non-negative.")
        object.__setattr__(self, "iteration", int(self.iteration))
        for name in ("rms", "roughness", "lagrange", "stepsize"):
            value = float(getattr(self, name))
            if math.isinf(value):
                raise ValueError(f"{name} cannot be infinite.")
            if math.isfinite(value) and value < 0:
                raise ValueError(f"{name} cannot be negative.")
            object.__setattr__(self, name, value)

    @property
    def is_complete(self) -> bool:
        """Whether RMS and roughness are both available."""
        return math.isfinite(self.rms) and math.isfinite(self.roughness)

    def to_dict(self) -> dict[str, int | float]:
        """Return a plain record mapping."""
        return asdict(self)


class Occam1DLog(Occam1DBase):
    r"""Represent a parsed Occam1D convergence log.

    Parameters
    ----------
    records : iterable of Occam1DLogRecord or mappings, optional
        Iteration records in source order. Mapping values must be accepted by
        :class:`Occam1DLogRecord`.
    target_misfit : float, default=1.0
        Positive reference RMS used by :attr:`has_converged` and summaries.
    verbose, logger, path, metadata, stream
        Shared options inherited from :class:`Occam1DBase`.

    Attributes
    ----------
    iterations : ndarray of int, shape (n_iterations,)
        Iteration numbers in source order.
    rms, roughness, lagrange, stepsize : ndarray of float
        Aligned diagnostic vectors. Missing native metrics are ``nan``.
    target_misfit : float
        Default normalized RMS convergence target.

    Raises
    ------
    TypeError
        If records or target values have unsupported types.
    ValueError
        If records contain negative or infinite diagnostics, or if the target
        is not finite and strictly positive.

    Notes
    -----
    The best iteration minimizes finite RMS. It is not necessarily the final
    iteration because later Occam steps may trade data fit for smoothness or
    encounter line-search difficulties.

    Repeated iteration numbers are preserved because restarted solver logs may
    append a second attempt. :meth:`get_iteration` returns the last occurrence
    by default and can return all occurrences when requested.

    Examples
    --------
    Parse a native log and select its best fit:

    >>> log = Occam1DLog.read("run/occam1d.log")  # doctest: +SKIP
    >>> log.best_iteration  # doctest: +SKIP
    12
    >>> log.converged(target=1.0)  # doctest: +SKIP
    True

    Construct records directly for testing or external adapters:

    >>> log = Occam1DLog([
    ...     {"iteration": 1, "rms": 2.0, "roughness": 10.0},
    ...     {"iteration": 2, "rms": 0.95, "roughness": 12.0},
    ... ])
    >>> log.best_iteration
    2
    >>> log.has_converged
    True
    """

    def __init__(
        self,
        records: Iterable[Occam1DLogRecord | dict[str, Any]] | None = None,
        *,
        target_misfit: float = 1.0,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        self.target_misfit = self._target(target_misfit)
        self._records = self._coerce_records(records)
        self._synchronize_arrays()
        self._validate_state()

    @staticmethod
    def _target(value: Real) -> float:
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError("target_misfit must be a real number.")
        target = float(value)
        if not math.isfinite(target) or target <= 0:
            raise ValueError("target_misfit must be finite and positive.")
        return target

    @staticmethod
    def _coerce_records(records):
        if records is None:
            return tuple()
        if isinstance(records, (str, bytes)):
            raise TypeError(
                "records must be an iterable of iteration records."
            )
        converted = []
        for record in records:
            if isinstance(record, Occam1DLogRecord):
                converted.append(record)
            elif isinstance(record, Mapping):
                converted.append(Occam1DLogRecord(**record))
            else:
                raise TypeError(
                    "records must contain Occam1DLogRecord objects or "
                    "mappings."
                )
        return tuple(converted)

    def _synchronize_arrays(self) -> None:
        """Build aligned NumPy views from immutable records."""
        self.iterations = np.asarray(
            [record.iteration for record in self._records], dtype=int
        )
        for name in ("rms", "roughness", "lagrange", "stepsize"):
            setattr(
                self,
                name,
                np.asarray(
                    [getattr(record, name) for record in self._records],
                    dtype=float,
                ),
            )

    @classmethod
    @validate_params(LOG_READ_SCHEMA)
    def read(cls, path, *, strict=True, **kwargs):
        r"""Parse an Occam1D convergence log.

        Parameters
        ----------
        path : path-like
            Native text log.
        strict : bool, default=True
            If ``True``, a file with no recognizable iteration blocks raises
            :class:`ValueError`. If ``False``, an empty bound log object is
            returned with a warning record.
        **kwargs
            Forwarded to the constructor, including ``target_misfit`` and base
            logging options.

        Returns
        -------
        Occam1DLog
            Parsed convergence history bound to ``path``.

        Raises
        ------
        FileNotFoundError
            If ``path`` does not exist.
        IsADirectoryError
            If ``path`` is not a regular file.
        ValueError
            If strict parsing finds no iteration blocks or a recognized metric
            is scientifically invalid.
        """
        probe = cls(**kwargs)
        source = probe._require_input_file(path, purpose="Occam1D log")
        parsed: list[Occam1DLogRecord] = []
        current: dict[str, int | float] | None = None

        def save() -> None:
            if current is not None:
                parsed.append(Occam1DLogRecord(**current))

        for line_number, raw in enumerate(
            source.read_text(encoding="utf8", errors="replace").splitlines(),
            start=1,
        ):
            header = _ITERATION.search(raw)
            if header is not None:
                save()
                current = {
                    "iteration": int(header.group(1)),
                    "rms": _NAN,
                    "roughness": _NAN,
                    "lagrange": _NAN,
                    "stepsize": _NAN,
                }
            if current is None:
                continue
            for name in ("rms", "roughness", "lagrange", "stepsize"):
                value = _metric(raw, name)
                if value is not None:
                    if not math.isfinite(value) or value < 0:
                        raise ValueError(
                            f"Invalid {name} at {source}:{line_number}: "
                            f"{value!r}."
                        )
                    current[name] = value
        save()

        if not parsed and strict:
            raise ValueError(
                f"No Occam iteration blocks were recognized in {source}."
            )
        obj = cls(parsed, **kwargs)
        obj._bind_path(source)
        if not parsed:
            obj.add_warning("No iteration blocks were recognized.")
        obj._record_parse_warnings()
        return obj

    def _record_parse_warnings(self) -> None:
        """Record incomplete and non-monotonic history diagnostics."""
        incomplete = [
            record.iteration
            for record in self._records
            if not record.is_complete
        ]
        if incomplete:
            values = ", ".join(str(value) for value in incomplete)
            self.add_warning(
                "Incomplete RMS/roughness diagnostics at iteration(s): "
                f"{values}."
            )
        if self.iterations.size > 1 and np.any(np.diff(self.iterations) <= 0):
            self.add_warning(
                "Iteration numbers are repeated or non-monotonic; the log may "
                "contain a restarted run."
            )

    @property
    def records(self) -> tuple[Occam1DLogRecord, ...]:
        """Immutable parsed records in native source order."""
        return self._records

    @property
    def n_iterations(self) -> int:
        """Number of parsed iteration blocks, including restart duplicates."""
        return len(self._records)

    n_iter = n_iterations

    @property
    def best_index(self) -> int | None:
        """Zero-based record index with minimum finite RMS."""
        finite = np.isfinite(self.rms)
        if not np.any(finite):
            return None
        return int(np.argmin(np.where(finite, self.rms, np.inf)))

    @property
    def best_iteration(self) -> int:
        """Iteration number with minimum finite RMS, or zero if unavailable."""
        index = self.best_index
        return 0 if index is None else int(self.iterations[index])

    @property
    def best_rms(self) -> float:
        """Minimum finite RMS, or ``nan`` when no fit metric is available."""
        index = self.best_index
        return _NAN if index is None else float(self.rms[index])

    @property
    def initial_rms(self) -> float:
        """First finite RMS in source order."""
        finite = self.rms[np.isfinite(self.rms)]
        return float(finite[0]) if finite.size else _NAN

    @property
    def final_rms(self) -> float:
        """Last finite RMS in source order."""
        finite = self.rms[np.isfinite(self.rms)]
        return float(finite[-1]) if finite.size else _NAN

    @property
    def has_converged(self) -> bool:
        """Whether any iteration reaches ``target_misfit``."""
        return self.converged(self.target_misfit)

    @validate_params(LOG_CONVERGENCE_SCHEMA)
    def converged(self, target=1.0) -> bool:
        """Return whether a finite RMS is at or below ``target``."""
        finite = self.rms[np.isfinite(self.rms)]
        return bool(finite.size and np.min(finite) <= target)

    @validate_params(LOG_ITERATION_SCHEMA)
    def get_iteration(self, iteration, *, all_matches=False):
        """Return the last or every record matching an iteration number.

        Raises :class:`KeyError` when no matching record exists.
        """
        matches = tuple(
            record for record in self._records if record.iteration == iteration
        )
        if not matches:
            raise KeyError(f"Iteration {iteration} is not present in the log.")
        return matches if all_matches else matches[-1]

    def to_records(self) -> list[dict[str, int | float]]:
        """Return independent dictionaries suitable for JSON serialization."""
        return [record.to_dict() for record in self._records]

    def to_array(self) -> np.ndarray:
        """Return ``(n_iterations, 5)`` numeric convergence data."""
        if not self.n_iterations:
            return np.empty((0, 5), dtype=float)
        return np.column_stack(
            (
                self.iterations,
                self.rms,
                self.roughness,
                self.lagrange,
                self.stepsize,
            )
        )

    @validate_params(LOG_WRITE_SCHEMA)
    def to_csv(self, path, *, include_missing=True) -> Path:
        """Export normalized convergence records as UTF-8 CSV."""
        target = self._prepare_output_file(path)
        with target.open("w", newline="", encoding="utf8") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                ["iteration", "rms", "roughness", "lagrange", "stepsize"]
            )
            for record in self._records:
                values = (
                    record.iteration,
                    record.rms,
                    record.roughness,
                    record.lagrange,
                    record.stepsize,
                )
                if include_missing or all(
                    math.isfinite(value) for value in values[1:]
                ):
                    writer.writerow(values)
        return target

    def summary(self) -> str:
        """Return a concise human-readable convergence summary."""
        if not self.n_iterations:
            return "Occam1DLog(no iteration records)"
        return (
            f"Occam1DLog(iterations={self.n_iterations}, "
            f"initial_rms={self.initial_rms:.6g}, "
            f"final_rms={self.final_rms:.6g}, "
            f"best_iteration={self.best_iteration}, "
            f"best_rms={self.best_rms:.6g}, "
            f"converged={self.has_converged})"
        )

    def diagnostics(self) -> dict[str, Any]:
        """Extend base diagnostics with convergence statistics."""
        values = super().diagnostics()
        values.update(
            {
                "n_iterations": self.n_iterations,
                "target_misfit": self.target_misfit,
                "initial_rms": self.initial_rms,
                "final_rms": self.final_rms,
                "best_iteration": self.best_iteration,
                "best_rms": self.best_rms,
                "has_converged": self.has_converged,
            }
        )
        return values

    def _validate_state(self) -> None:
        """Ensure record and array views remain aligned."""
        sizes = {
            self.iterations.size,
            self.rms.size,
            self.roughness.size,
            self.lagrange.size,
            self.stepsize.size,
            len(self._records),
        }
        if len(sizes) != 1:
            raise ValueError("Occam1D log diagnostic arrays are not aligned.")
