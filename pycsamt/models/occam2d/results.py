# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""InversionResult — unified post-inversion access layer.

Scans the working directory after an Occam2D run, loads all output files
(iter, response, log), provides 2-D resistivity model arrays, and
implements the ``iter2dat`` transformation (Bo Yang's format).

Usage
-----
>>> result = InversionResult(workdir="occam_run/")
>>> result.plot_model()
>>> result.plot_response()
>>> result.plot_misfit()
>>> result.iter2dat("final_model.dat")
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

import numpy as np

from .base     import OccamBase
from .log      import OccamLog
from .mesh     import OccamMesh
from .model    import OccamModel
from .response import OccamResponse
from .startup  import OccamIter

PathLike = Union[str, Path]

__all__ = ["InversionResult"]


class InversionResult(OccamBase):
    """Post-inversion result container.

    Parameters
    ----------
    workdir : path-like
        Directory produced by an ``OccamRunner`` run.
    iteration : int or None
        Which iteration to load.  ``None`` → load the last (highest)
        numbered ``.iter`` file.

    Attributes
    ----------
    workdir : Path
    log : OccamLog
    mesh : OccamMesh
    model : OccamModel
    best_iter : OccamIter
        The selected (or last) iteration.
    response : OccamResponse
        Response matching the selected iteration.
    iter_files : list[Path]
        All ``.iter`` files found, sorted by iteration number.
    resp_files : list[Path]
        All ``.resp`` files found, sorted by iteration number.
    rho_2d : np.ndarray, shape (n_zcells, n_xcells)
        Resistivity model (Ω·m) on the mesh grid.
    """

    def __init__(
        self,
        workdir: PathLike = ".",
        iteration: Optional[int] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.workdir    = Path(workdir)
        self._iteration = iteration

        self.log:       Optional[OccamLog]      = None
        self.mesh:      Optional[OccamMesh]     = None
        self.model:     Optional[OccamModel]    = None
        self.best_iter: Optional[OccamIter]     = None
        self.response:  Optional[OccamResponse] = None

        self.iter_files: List[Path] = []
        self.resp_files: List[Path] = []
        self.rho_2d:     Optional[np.ndarray] = None

        self._load()

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------
    def _load(self) -> None:
        """Scan *workdir* and load all available output files."""
        wd = self.workdir
        if not wd.is_dir():
            raise NotADirectoryError(f"workdir not found: {wd}")

        # TODO: implement
        # 1. Find and sort all .iter and .resp files
        # 2. Load OccamLog, OccamMesh, OccamModel
        # 3. Load the requested (or last) OccamIter + matching OccamResponse
        # 4. Build rho_2d from iter.param_values + model.param_map + mesh
        raise NotImplementedError("InversionResult._load — not yet implemented")

    # ------------------------------------------------------------------
    # iter2dat (Bo Yang's post-processing format)
    # ------------------------------------------------------------------
    def iter2dat(self, output_file: PathLike) -> Path:
        """Write a ``iter2dat``-style ASCII model file.

        The format stores (x, z, log10_rho) triplets suitable for
        contouring or import into third-party visualisation tools.

        Parameters
        ----------
        output_file : path-like
            Output path for the ``.dat`` file.

        Returns
        -------
        Path
            The written file path.
        """
        # TODO: implement
        # 1. Compute cell centres (x, z) from mesh node positions
        # 2. Map param_values → resistivity using model.param_map
        # 3. Write ASCII columns: x  z  log10(rho)
        raise NotImplementedError("InversionResult.iter2dat — not yet implemented")

    # ------------------------------------------------------------------
    # Plotting (delegates to plot.py)
    # ------------------------------------------------------------------
    def plot_model(self, **kwargs):
        """Plot the 2-D resistivity model (delegates to ``PlotModel``)."""
        from .plot import PlotModel
        return PlotModel(result=self, **kwargs).plot()

    def plot_response(self, **kwargs):
        """Plot observed vs predicted response (delegates to ``PlotResponse``)."""
        from .plot import PlotResponse
        return PlotResponse(result=self, **kwargs).plot()

    def plot_misfit(self, **kwargs):
        """Plot RMS misfit vs iteration (delegates to ``PlotMisfit``)."""
        from .plot import PlotMisfit
        return PlotMisfit(result=self, **kwargs).plot()

    def plot_pseudo(self, **kwargs):
        """Plot pseudosection of observed data (delegates to ``PlotPseudo``)."""
        from .plot import PlotPseudo
        return PlotPseudo(result=self, **kwargs).plot()

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    @property
    def final_rms(self) -> float:
        """RMS misfit of the selected iteration."""
        if self.best_iter is None:
            return float("nan")
        return self.best_iter.misfit_value

    @property
    def n_iterations(self) -> int:
        return len(self.iter_files)

    def summary(self) -> str:
        return (
            f"InversionResult\n"
            f"  workdir    : {self.workdir}\n"
            f"  iterations : {self.n_iterations}\n"
            f"  final RMS  : {self.final_rms:.4f}\n"
            f"  converged  : {self.log.converged if self.log else 'N/A'}\n"
        )
