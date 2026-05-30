# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamResponse — read Occam2D response (.resp) files.

A response file is produced by the Fortran binary after each iteration.
It contains seven whitespace-delimited columns per data point — no header:

::

    site  freq  type  err_floor  observed  modeled  residual
    1.0   1.0   5.0   0.0        4.8243    4.9036   -0.12094E-01
    ...

Columns
-------
0  site_index   (1-based float → int)
1  freq_index   (1-based float → int)
2  type_code    (int; same codes as OccamData)
3  error_floor  (float; often 0.0)
4  observed     (float; observed datum)
5  modeled      (float; forward-model prediction)
6  residual     (float; weighted residual)

Entry point
-----------
``OccamResponse.read(path)``
    Parse a ``.resp`` file and populate observed / predicted arrays.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import numpy as np

from .base import OccamBase

PathLike = Union[str, Path]

__all__ = ["OccamResponse"]

_N_COLS = 7   # expected number of columns per data row


# -----------------------------------------------------------------------
# Low-level parser
# -----------------------------------------------------------------------

def _parse_response(path: Path) -> np.ndarray:
    """Return a (n_data, 7) float array from a ``.resp`` file.

    Raises
    ------
    FileNotFoundError
    ValueError
        If no valid data rows can be parsed.
    """
    if not path.exists():
        raise FileNotFoundError(f"Occam response file not found: {path}")

    rows: list[list[float]] = []
    with path.open("r", errors="replace") as fh:
        for raw in fh:
            tokens = raw.strip().split()
            if len(tokens) == _N_COLS:
                try:
                    rows.append([float(t) for t in tokens])
                except ValueError:
                    pass

    if not rows:
        raise ValueError(f"No valid data rows found in response file: {path}")

    return np.array(rows, dtype=float)


# -----------------------------------------------------------------------
# OccamResponse
# -----------------------------------------------------------------------

class OccamResponse(OccamBase):
    """Occam2D response-file container.

    Attributes
    ----------
    data : np.ndarray, shape (n_data, 7)
        Full raw table from the ``.resp`` file.  Columns:
        site_index, freq_index, type_code, error_floor,
        observed, modeled, residual.
    observed : np.ndarray, shape (n_data,)
        Observed data values (column 4).
    modeled : np.ndarray, shape (n_data,)
        Forward-model predictions (column 5).
    residuals : np.ndarray, shape (n_data,)
        Weighted residuals (column 6).
    rms : float
        RMS of the weighted residuals,
        ``sqrt(mean(residuals ** 2))``.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.data:      np.ndarray = np.empty((0, _N_COLS))
        self.observed:  np.ndarray = np.array([])
        self.modeled:   np.ndarray = np.array([])
        self.residuals: np.ndarray = np.array([])
        self.rms:       float      = 0.0

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(
        cls,
        path: PathLike,
        data_fn: Optional[PathLike] = None,
        **kwargs,
    ) -> "OccamResponse":
        """Parse a ``.resp`` file.

        Parameters
        ----------
        path : path-like
            Path to the ``.resp`` file.
        data_fn : path-like, optional
            Unused; reserved for future cross-checking with the data file.

        Returns
        -------
        OccamResponse

        Raises
        ------
        FileNotFoundError
        ValueError
            If no valid data rows can be parsed.
        """
        p   = Path(path)
        arr = _parse_response(p)
        obj = cls(**kwargs)

        obj.data      = arr
        obj.observed  = arr[:, 4]
        obj.modeled   = arr[:, 5]
        obj.residuals = arr[:, 6]
        obj.rms       = float(np.sqrt(np.mean(obj.residuals ** 2)))

        if obj.verbose:
            obj.logger.info(
                "OccamResponse.read: %d data points, rms=%.4f from %s",
                obj.n_data, obj.rms, p,
            )
        return obj

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    @property
    def n_data(self) -> int:
        return len(self.data)

    @property
    def site_indices(self) -> np.ndarray:
        """1-based site indices (int)."""
        return self.data[:, 0].astype(int) if self.data.size else np.array([], dtype=int)

    @property
    def freq_indices(self) -> np.ndarray:
        """1-based frequency indices (int)."""
        return self.data[:, 1].astype(int) if self.data.size else np.array([], dtype=int)

    @property
    def type_codes(self) -> np.ndarray:
        """Unique data-type codes present in this response."""
        if not self.data.size:
            return np.array([], dtype=int)
        return np.unique(self.data[:, 2].astype(int))

    # ------------------------------------------------------------------
    # Derived per-site / per-frequency misfit
    # ------------------------------------------------------------------
    def misfit_per_site(self) -> dict[int, float]:
        """Return per-site RMS misfit keyed by 1-based site index."""
        if not self.data.size:
            return {}
        result: dict[int, float] = {}
        for s in np.unique(self.site_indices):
            mask = self.site_indices == s
            result[int(s)] = float(np.sqrt(np.mean(self.residuals[mask] ** 2)))
        return result

    def misfit_per_frequency(self) -> dict[int, float]:
        """Return per-frequency RMS misfit keyed by 1-based frequency index."""
        if not self.data.size:
            return {}
        result: dict[int, float] = {}
        for f in np.unique(self.freq_indices):
            mask = self.freq_indices == f
            result[int(f)] = float(np.sqrt(np.mean(self.residuals[mask] ** 2)))
        return result