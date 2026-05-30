# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamStartup and OccamIter — Startup / .iter file handling.

Both file types share the ``OCCAMITER_FLEX`` format header.  The
distinction is that a Startup file has ``Iteration: 0`` and a uniform
initial model vector, while an ``.iter`` file carries a real iteration
number and the model parameters from that iteration.

Classes
-------
OccamStartup
    Builds and writes the initial Startup file.
OccamIter
    Reads any ``.iter`` file (iteration N ≥ 1) produced by Occam.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import numpy as np

from .base   import OccamBase
from .config import OccamConfig

PathLike = Union[str, Path]

__all__ = ["OccamStartup", "OccamIter"]

_HEADER_FORMAT = "OCCAMITER_FLEX"


class OccamStartup(OccamBase):
    """Occam2D startup (iteration-0) file.

    Attributes
    ----------
    description : str
    model_file : str
    data_file : str
    max_iterations : int
    target_misfit : float
    roughness_type : int
    diagonal_penalties : int
    stepsize_cut_count : int
    debug_level : int
    lagrange_start : float
    initial_rho : float
        Uniform starting resistivity (Ω·m); log10 stored in param vector.
    n_params : int
        Set automatically from the model when writing.
    param_values : np.ndarray, shape (n_params,)
        Log10 resistivity values; uniform ``log10(initial_rho)`` at startup.
    """

    def __init__(
        self,
        config: Optional[OccamConfig] = None,
        description: str = "startup created by pycsamt",
        **kwargs,
    ):
        super().__init__(**kwargs)
        cfg = config or OccamConfig()
        self.config             = cfg
        self.description        = description
        self.model_file         = cfg.model_file
        self.data_file          = cfg.data_file
        self.max_iterations     = cfg.max_iterations
        self.target_misfit      = cfg.target_misfit
        self.roughness_type     = cfg.roughness_type
        self.diagonal_penalties = cfg.diagonal_penalties
        self.stepsize_cut_count = cfg.stepsize_cut_count
        self.debug_level        = cfg.debug_level
        self.lagrange_start     = cfg.lagrange_start
        self.initial_rho        = cfg.initial_rho
        self.n_params           = 0
        self.param_values:  np.ndarray = np.array([])

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_model(
        cls,
        model: "OccamModel",  # noqa: F821
        config: Optional[OccamConfig] = None,
        **kwargs,
    ) -> "OccamStartup":
        """Build a startup file from a populated model object.

        Sets the parameter vector to ``log10(config.initial_rho)``
        repeated ``model.n_params`` times.
        """
        # TODO: implement
        raise NotImplementedError("OccamStartup.from_model — not yet implemented")

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamStartup":
        """Parse a Startup file (Iteration == 0)."""
        # TODO: implement
        raise NotImplementedError("OccamStartup.read — not yet implemented")

    def write(self, path: PathLike) -> Path:
        """Write startup file to *path*."""
        # TODO: implement
        raise NotImplementedError("OccamStartup.write — not yet implemented")


class OccamIter(OccamBase):
    """Reader for Occam ``.iter`` files (iteration N ≥ 1).

    Attributes
    ----------
    iteration : int
    lagrange_value : float
    roughness_value : float
    misfit_value : float
    misfit_reached : bool
    n_params : int
    param_values : np.ndarray, shape (n_params,)
        Log10 resistivity values at this iteration.
    model_file : str
    data_file : str
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.iteration:      int         = 0
        self.lagrange_value: float       = 0.0
        self.roughness_value: float      = 0.0
        self.misfit_value:   float       = 0.0
        self.misfit_reached: bool        = False
        self.n_params:       int         = 0
        self.param_values:   np.ndarray  = np.array([])
        self.model_file:     str         = ""
        self.data_file:      str         = ""

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamIter":
        """Parse a ``.iter`` file at *path*."""
        # TODO: implement OCCAMITER_FLEX parser
        raise NotImplementedError("OccamIter.read — not yet implemented")

    def to_resistivity(self) -> "np.ndarray":
        """Return resistivity model (Ω·m) from log10 param_values."""
        return 10.0 ** self.param_values
