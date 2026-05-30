# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamLog — parse Occam2D log files.

The log file (typically ``LogFile.logfile``) records per-iteration
convergence statistics written by the Fortran binary.

Entry point
-----------
``OccamLog.read(path)``
    Parse the log and populate per-iteration arrays.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np

from .base import OccamBase

PathLike = Union[str, Path]

__all__ = ["OccamLog"]


class OccamLog(OccamBase):
    """Occam2D log-file container.

    Attributes
    ----------
    iterations : np.ndarray, shape (n_iter,)
        Iteration numbers.
    rms : np.ndarray, shape (n_iter,)
        RMS misfit per iteration.
    roughness : np.ndarray, shape (n_iter,)
        Roughness (deGroot-Hedlin measure) per iteration.
    lagrange : np.ndarray, shape (n_iter,)
        Lagrange multiplier per iteration.
    stepsize : np.ndarray, shape (n_iter,)
        Stepsize cut count per iteration.
    fminocc : np.ndarray, shape (n_iter,)
        Minimum tolerance value.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.iterations: np.ndarray = np.array([], dtype=int)
        self.rms:        np.ndarray = np.array([])
        self.roughness:  np.ndarray = np.array([])
        self.lagrange:   np.ndarray = np.array([])
        self.stepsize:   np.ndarray = np.array([])
        self.fminocc:    np.ndarray = np.array([])

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamLog":
        """Parse an Occam log file at *path*.

        Populates all per-iteration arrays.
        """
        # TODO: implement — parse keyword lines:
        #   "OCCAM ITERATION", "RMS", "ROUGHNESS", "LAGRANGE", etc.
        raise NotImplementedError("OccamLog.read — not yet implemented")

    # ------------------------------------------------------------------
    # Derived
    # ------------------------------------------------------------------
    @property
    def converged(self) -> bool:
        """True if the final RMS misfit is ≤ target (≈ 1.0)."""
        return bool(self.rms.size and self.rms[-1] <= 1.0)

    @property
    def best_iteration(self) -> int:
        """Index of the iteration with lowest RMS misfit."""
        if not self.rms.size:
            return 0
        return int(np.argmin(self.rms))
