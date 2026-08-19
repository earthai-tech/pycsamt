# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Normalized forward responses produced by an Occam1D executable.

Two solver dialects are supported: canonical seven-column response tables and
``EMData_1.x``-style rows beginning with data type ``103`` or ``104``. Both
are normalized to one seven-column representation before scientific queries
or plotting are performed.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Any

import numpy as np

from ...compat.sklearn import validate_params
from .base import Occam1DBase
from .schema import (
    RESPONSE_READ_SCHEMA,
    RESPONSE_SELECT_SCHEMA,
    RESPONSE_WRITE_SCHEMA,
)

__all__ = ["Occam1DResponse"]

_N_COLUMNS = 7
_RHO_CODE = 103
_PHASE_CODE = 104
_KNOWN_CODES = frozenset({_RHO_CODE, _PHASE_CODE})
_COLUMN_NAMES = (
    "site_index",
    "frequency_index",
    "type_code",
    "error",
    "observed",
    "modeled",
    "weighted_residual",
)
_FORTRAN_EXPONENT = re.compile(r"(?<=\d)[Dd](?=[-+]?\d)")


def _native_float(token: str, line_number: int) -> float:
    """Parse a response number, including a Fortran ``D`` exponent."""
    try:
        return float(_FORTRAN_EXPONENT.sub("E", token))
    except ValueError as error:
        raise ValueError(
            f"Invalid numeric token {token!r} at line {line_number}."
        ) from error


class Occam1DResponse(Occam1DBase):
    r"""Represent observed, modeled, and weighted Occam1D responses.

    Parameters
    ----------
    rows : array-like of shape (n_data, 7), optional
        Normalized rows containing ``site_index``, ``frequency_index``,
        ``type_code``, solver-unit ``error``, solver-unit ``observed``,
        solver-unit ``modeled``, and ``weighted_residual``. An omitted value
        creates a valid empty response container.
    verbose, logger, path, metadata, stream
        Shared options inherited from :class:`Occam1DBase`.

    Attributes
    ----------
    rows : ndarray of float, shape (n_data, 7)
        Defensive copy of the normalized response table.

    Raises
    ------
    TypeError
        If response rows cannot be converted to real numbers.
    ValueError
        If table shape, indices, codes, errors, or numeric values violate
        response invariants.

    Notes
    -----
    Type ``103`` stores apparent resistivity and its prediction as
    :math:`\log_{10}(\rho_a)`; type ``104`` stores phase in degrees. Weighted
    residuals are normally ``(observed - modeled) / error``. Canonical solver
    rows retain their reported residual because additional internal weighting
    may have been applied.

    Missing numerical results remain ``nan``. Infinite values are rejected.
    Site and frequency indices are one-based.

    Examples
    --------
    >>> response = Occam1DResponse([
    ...     [1, 1, 103, 0.1, 2.0, 2.1, -1.0],
    ...     [1, 1, 104, 2.0, 45.0, 44.0, 0.5],
    ... ])
    >>> response.n_data, response.rms
    (2, 0.7905694150420949)
    >>> response.physical_values()[0].tolist()
    [100.0, 45.0]
    """

    def __init__(self, rows=None, **kwargs):
        super().__init__(**kwargs)
        self.rows = self._table(rows)
        self._validate_state()
        self._record_warnings()

    @staticmethod
    def _table(rows) -> np.ndarray:
        """Return an independent normalized response table."""
        if rows is None:
            return np.empty((0, _N_COLUMNS), dtype=float)
        try:
            result = np.array(rows, dtype=float, copy=True)
        except (TypeError, ValueError) as error:
            raise TypeError(
                "response rows must contain real numbers."
            ) from error
        if result.size == 0:
            return np.empty((0, _N_COLUMNS), dtype=float)
        if result.ndim != 2 or result.shape[1] != _N_COLUMNS:
            raise ValueError("response rows must have shape (n_data, 7).")
        return result

    def _validate_state(self) -> None:
        """Validate normalized table shape and scientific row invariants."""
        if self.rows.ndim != 2 or self.rows.shape[1] != _N_COLUMNS:
            raise ValueError("response rows must have shape (n_data, 7).")
        if np.any(np.isinf(self.rows)):
            raise ValueError("response rows cannot contain infinite values.")
        if not self.n_data:
            return
        for column, name in ((0, "site"), (1, "frequency"), (2, "type")):
            values = self.rows[:, column]
            if np.any(~np.isfinite(values)) or np.any(values < 1):
                raise ValueError(
                    f"{name} indices/codes must be finite positive integers."
                )
            if np.any(values != np.floor(values)):
                raise ValueError(
                    f"{name} indices/codes must be integer-valued."
                )
        errors = self.rows[:, 3]
        finite_error = np.isfinite(errors)
        if np.any(errors[finite_error] < 0):
            raise ValueError("Finite response errors must be non-negative.")
        keys = self.rows[:, :3].astype(int)
        if np.unique(keys, axis=0).shape[0] != self.n_data:
            raise ValueError(
                "Response rows must have unique site/frequency/type keys."
            )

    def _record_warnings(self) -> None:
        """Record recoverable completeness and type-code diagnostics."""
        if not self.n_data:
            return
        unknown = sorted(set(self.type_codes).difference(_KNOWN_CODES))
        if unknown:
            values = ", ".join(str(value) for value in unknown)
            self.add_warning(f"Unsupported response type code(s): {values}.")
        incomplete = np.any(~np.isfinite(self.rows[:, 3:]), axis=1)
        if np.any(incomplete):
            self.add_warning(
                f"Response contains {int(np.count_nonzero(incomplete))} "
                "incomplete row(s)."
            )

    @classmethod
    @validate_params(RESPONSE_READ_SCHEMA)
    def read(cls, path, *, strict=True, **kwargs):
        """Read and normalize a solver response file.

        Parameters
        ----------
        path : path-like
            Response file produced for one solver iteration.
        strict : bool, default=True
            Raise on numeric-looking malformed rows. When ``False``, malformed
            rows are skipped and recorded as an object warning.
        **kwargs
            Shared constructor options.

        Returns
        -------
        Occam1DResponse
            Parsed response bound to the resolved input path.

        Raises
        ------
        FileNotFoundError, IsADirectoryError
            If the input is missing or not a regular file.
        ValueError
            If no response rows are found or a strict parse encounters a
            malformed numeric row.
        """
        probe = cls(**kwargs)
        source = probe._require_input_file(path, purpose="Occam1D response")
        rows = []
        skipped = 0
        for line_number, raw in enumerate(
            source.read_text(encoding="utf8", errors="replace").splitlines(),
            start=1,
        ):
            text = raw.strip()
            if not text or text.startswith(("!", "#")):
                continue
            fields = text.split()
            try:
                _native_float(fields[0], line_number)
            except ValueError:
                continue
            try:
                values = [
                    _native_float(field, line_number) for field in fields
                ]
                row = cls._normalize_row(values, line_number)
            except ValueError:
                if strict:
                    raise
                skipped += 1
                continue
            rows.append(row)
        if not rows:
            raise ValueError(f"No response rows found in {source}.")
        obj = cls(rows, **kwargs)
        obj._bind_path(source)
        if skipped:
            obj.add_warning(
                f"Skipped {skipped} malformed numeric response row(s)."
            )
        obj.logger.debug(
            "Read %d Occam1D response rows from %s.", obj.n_data, source
        )
        return obj

    @staticmethod
    def _normalize_row(values, line_number):
        """Normalize one canonical or EMData-style numeric row."""
        if not values:
            raise ValueError(
                f"Empty numeric response row at line {line_number}."
            )
        first = values[0]
        if (
            math.isfinite(first)
            and first == math.floor(first)
            and int(first) in _KNOWN_CODES
        ):
            if len(values) < 7:
                raise ValueError(
                    f"EMData response row at line {line_number} requires "
                    "at least seven columns."
                )
            code = int(first)
            frequency_index = values[1]
            observed, error, modeled = values[4:7]
            residual = (
                (observed - modeled) / error
                if math.isfinite(error) and error > 0
                else float("nan")
            )
            return [
                1.0,
                frequency_index,
                float(code),
                error,
                observed,
                modeled,
                residual,
            ]
        if len(values) != _N_COLUMNS:
            raise ValueError(
                f"Canonical response row at line {line_number} must have "
                "exactly seven columns."
            )
        return values

    @property
    def n_data(self) -> int:
        """Number of normalized scalar response rows."""
        return int(self.rows.shape[0])

    @property
    def site_index(self) -> np.ndarray:
        """One-based site indices for every response row."""
        return self.rows[:, 0].astype(int)

    site_indices = site_index

    @property
    def frequency_index(self) -> np.ndarray:
        """One-based frequency indices for every response row."""
        return self.rows[:, 1].astype(int)

    frequency_indices = frequency_index

    @property
    def type_code(self) -> np.ndarray:
        """Occam data type code for every response row."""
        return self.rows[:, 2].astype(int)

    @property
    def type_codes(self) -> np.ndarray:
        """Sorted unique Occam data type codes present."""
        return np.unique(self.type_code)

    @property
    def error(self) -> np.ndarray:
        """Observation errors in solver units."""
        return self.rows[:, 3]

    @property
    def observed(self) -> np.ndarray:
        """Observed values in solver units."""
        return self.rows[:, 4]

    @property
    def modeled(self) -> np.ndarray:
        """Modeled values in solver units."""
        return self.rows[:, 5]

    @property
    def residuals(self) -> np.ndarray:
        """Weighted residuals reported or derived from the solver row."""
        return self.rows[:, 6]

    @property
    def finite_residual_mask(self) -> np.ndarray:
        """Boolean mask identifying usable weighted residuals."""
        return np.isfinite(self.residuals)

    @property
    def n_residuals(self) -> int:
        """Number of finite weighted residuals contributing to RMS."""
        return int(np.count_nonzero(self.finite_residual_mask))

    @property
    def rms(self) -> float:
        r"""Root-mean-square finite weighted residual."""
        finite = self.residuals[self.finite_residual_mask]
        return float(np.sqrt(np.mean(finite**2))) if finite.size else np.nan

    @property
    def mean_residual(self) -> float:
        """Mean finite weighted residual, useful for bias detection."""
        finite = self.residuals[self.finite_residual_mask]
        return float(np.mean(finite)) if finite.size else np.nan

    @property
    def mean_absolute_residual(self) -> float:
        """Mean absolute finite weighted residual."""
        finite = self.residuals[self.finite_residual_mask]
        return float(np.mean(np.abs(finite))) if finite.size else np.nan

    @property
    def chi_square(self) -> float:
        """Sum of squared finite weighted residuals."""
        finite = self.residuals[self.finite_residual_mask]
        return float(np.sum(finite**2)) if finite.size else np.nan

    def physical_values(self) -> tuple[np.ndarray, np.ndarray]:
        """Return observed and modeled responses in physical units.

        Type ``103`` values are converted from log10 resistivity to ohm
        metres. Phase and unknown type codes pass through unchanged.
        """
        observed = self.observed.copy()
        modeled = self.modeled.copy()
        rho = self.type_code == _RHO_CODE
        with np.errstate(over="ignore", invalid="ignore"):
            observed[rho] = np.power(10.0, observed[rho])
            modeled[rho] = np.power(10.0, modeled[rho])
        if np.any(np.isinf(observed)) or np.any(np.isinf(modeled)):
            raise ValueError(
                "Log-resistivity response exceeds physical float range."
            )
        return observed, modeled

    def physical_errors(self) -> np.ndarray:
        """Return approximate one-sigma errors in physical response units.

        For resistivity, the exact symmetric log-space interval is not
        symmetric in physical units. This method returns the differential
        approximation ``rho * ln(10) * error`` used by response plots.
        """
        observed, _ = self.physical_values()
        errors = self.error.copy()
        rho = self.type_code == _RHO_CODE
        errors[rho] = observed[rho] * math.log(10.0) * errors[rho]
        return errors

    def computed_residuals(self) -> np.ndarray:
        """Compute ``(observed - modeled) / error`` in solver units."""
        return np.divide(
            self.observed - self.modeled,
            self.error,
            out=np.full(self.n_data, np.nan),
            where=np.isfinite(self.error) & (self.error > 0),
        )

    @validate_params(RESPONSE_SELECT_SCHEMA)
    def select(self, *, type_code=None, frequency_index=None):
        """Return an independent response containing matching rows."""
        mask = np.ones(self.n_data, dtype=bool)
        if type_code is not None:
            mask &= self.type_code == type_code
        if frequency_index is not None:
            mask &= self.frequency_index == frequency_index
        return type(self)(
            self.rows[mask],
            metadata=self.metadata,
            verbose=self.verbose,
        )

    def misfit_per_frequency(self) -> dict[int, float]:
        """Return finite-residual RMS grouped by one-based frequency index."""
        return self._grouped_rms(self.frequency_index)

    def misfit_per_type(self) -> dict[int, float]:
        """Return finite-residual RMS grouped by Occam data type code."""
        return self._grouped_rms(self.type_code)

    def _grouped_rms(self, groups) -> dict[int, float]:
        result = {}
        for value in np.unique(groups):
            residuals = self.residuals[
                (groups == value) & self.finite_residual_mask
            ]
            result[int(value)] = (
                float(np.sqrt(np.mean(residuals**2)))
                if residuals.size
                else np.nan
            )
        return result

    def validate_against(
        self,
        data,
        *,
        rtol=1e-6,
        atol=1e-9,
    ) -> None:
        """Validate response indices and observations against sounding data.

        Parameters
        ----------
        data : Occam1DData
            Sounding used by the inversion.
        rtol, atol : float
            Relative and absolute tolerances for solver-unit observations.

        Raises
        ------
        TypeError
            If ``data`` is not an :class:`Occam1DData`.
        ValueError
            If indices are out of range, a type is unsupported, or observed
            values do not match the sounding at their indexed frequencies.
        """
        from .data import Occam1DData

        if not isinstance(data, Occam1DData):
            raise TypeError("data must be an Occam1DData instance.")
        if isinstance(rtol, bool) or not isinstance(rtol, (int, float)):
            raise TypeError("rtol must be a non-negative real number.")
        if isinstance(atol, bool) or not isinstance(atol, (int, float)):
            raise TypeError("atol must be a non-negative real number.")
        if not math.isfinite(rtol) or rtol < 0:
            raise ValueError("rtol must be finite and non-negative.")
        if not math.isfinite(atol) or atol < 0:
            raise ValueError("atol must be finite and non-negative.")
        indices = self.frequency_index - 1
        if np.any(indices < 0) or np.any(indices >= data.n_frequencies):
            raise ValueError(
                "Response frequency index is outside the sounding range."
            )
        expected = np.full(self.n_data, np.nan)
        rho = self.type_code == _RHO_CODE
        phase = self.type_code == _PHASE_CODE
        if np.any(~(rho | phase)):
            raise ValueError(
                "Response contains types unsupported by Occam1D data."
            )
        with np.errstate(divide="ignore", invalid="ignore"):
            expected[rho] = np.log10(data.resistivity[indices[rho]])
        expected[phase] = data.phase[indices[phase]]
        if np.any(~np.isfinite(expected)):
            raise ValueError(
                "Response references a missing sounding observation."
            )
        if not np.allclose(
            self.observed, expected, rtol=rtol, atol=atol, equal_nan=False
        ):
            mismatch = np.flatnonzero(
                ~np.isclose(
                    self.observed,
                    expected,
                    rtol=rtol,
                    atol=atol,
                    equal_nan=False,
                )
            )[0]
            raise ValueError(
                "Response observation does not match sounding data at row "
                f"{int(mismatch) + 1}."
            )

    def to_records(self, *, physical=False) -> list[dict[str, int | float]]:
        """Return independent response dictionaries.

        Set ``physical=True`` to express observed, modeled, and error values
        in ohm metres or degrees.
        """
        table = self.rows.copy()
        if physical:
            table[:, 4], table[:, 5] = self.physical_values()
            table[:, 3] = self.physical_errors()
        records = []
        for row in table:
            record = dict(zip(_COLUMN_NAMES, row))
            for name in _COLUMN_NAMES[:3]:
                record[name] = int(record[name])
            records.append(record)
        return records

    @validate_params(RESPONSE_WRITE_SCHEMA)
    def to_csv(self, path, *, physical=True) -> Path:
        """Export a normalized response table as UTF-8 CSV."""
        target = self._prepare_output_file(path)
        with target.open("w", newline="", encoding="utf8") as stream:
            writer = csv.DictWriter(stream, fieldnames=_COLUMN_NAMES)
            writer.writeheader()
            writer.writerows(self.to_records(physical=physical))
        return target

    def summary(self) -> str:
        """Return a concise human-readable response summary."""
        return (
            f"Occam1DResponse(data={self.n_data}, "
            f"finite_residuals={self.n_residuals}, rms={self.rms:.6g}, "
            f"mean_residual={self.mean_residual:.6g}, "
            f"types={self.type_codes.tolist()})"
        )

    def diagnostics(self) -> dict[str, Any]:
        """Extend lifecycle diagnostics with response-fit statistics."""
        values = super().diagnostics()
        values.update(
            {
                "n_data": self.n_data,
                "n_residuals": self.n_residuals,
                "type_codes": self.type_codes.tolist(),
                "frequency_indices": np.unique(
                    self.frequency_index
                ).tolist(),
                "rms": self.rms,
                "mean_residual": self.mean_residual,
                "mean_absolute_residual": self.mean_absolute_residual,
                "chi_square": self.chi_square,
            }
        )
        return values
