# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamResponse — read Occam2D response (.resp) files.

A response file mirrors the data-file format but contains forward-modelled
values alongside the observed data.  It is produced by the Fortran binary
after each iteration and is used for misfit diagnostics and plot overlays.

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


class OccamResponse(OccamBase):
    """Occam2D response-file container.

    Attributes
    ----------
    sites : list[str]
    offsets : np.ndarray
    frequencies : np.ndarray
    observed : np.ndarray, shape (n_data, 5)
        Same layout as ``OccamData.data_blocks``.
    predicted : np.ndarray, shape (n_data, 5)
        Forward-modelled response at the same data points.
    residuals : np.ndarray, shape (n_data,)
        Normalised residuals ``(obs - pred) / error``.
    rms : float
        Global RMS misfit over all data.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.sites:      list         = []
        self.offsets:    np.ndarray   = np.array([])
        self.frequencies: np.ndarray  = np.array([])
        self.observed:   np.ndarray   = np.empty((0, 5))
        self.predicted:  np.ndarray   = np.empty((0, 5))
        self.residuals:  np.ndarray   = np.array([])
        self.rms:        float        = 0.0

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, data_fn: Optional[PathLike] = None, **kwargs) -> "OccamResponse":
        """Parse a ``.resp`` file.

        Parameters
        ----------
        path : path-like
            Path to the ``.resp`` file.
        data_fn : path-like, optional
            Path to the matching data file; used to cross-check site names
            and frequency ordering.
        """
        # TODO: implement
        raise NotImplementedError("OccamResponse.read — not yet implemented")

    # ------------------------------------------------------------------
    # Derived quantities
    # ------------------------------------------------------------------
    def misfit_per_site(self) -> dict[str, float]:
        """Return per-station RMS misfit."""
        # TODO: implement after read() is available
        raise NotImplementedError

    def misfit_per_frequency(self) -> dict[float, float]:
        """Return per-frequency RMS misfit."""
        # TODO: implement after read() is available
        raise NotImplementedError
