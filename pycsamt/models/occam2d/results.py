# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""InversionResult — unified post-inversion access layer.

Scans the working directory after an Occam2D run, loads all output files
(iter, response, log), provides a 2-D log₁₀-resistivity model array, and
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

import re
from pathlib import Path
from typing import List, Optional, Union

import numpy as np

from .base       import OccamBase
from .log        import OccamLog
from .mesh       import OccamMesh
from .model      import OccamModel
from .response   import OccamResponse
from .startup    import OccamIter
from .validation import is_log_file, is_mesh_file, is_model_file

PathLike = Union[str, Path]

__all__ = ["InversionResult"]


def _iter_number(path: Path) -> int:
    """Extract the leading or embedded integer from a filename stem."""
    m = re.search(r"\d+", path.stem)
    return int(m.group()) if m else -1


def _scan_one(wd: Path, predicate, glob_pat: str) -> Optional[Path]:
    """Return the first file matching *predicate* under *wd*."""
    for p in sorted(wd.glob(glob_pat)):
        if predicate(p):
            return p
    return None


# -----------------------------------------------------------------------
# rho_2d builder
# -----------------------------------------------------------------------

def _build_rho_2d(
    param_values: np.ndarray,
    model: OccamModel,
    mesh: OccamMesh,
) -> np.ndarray:
    """Return a (n_zcells, n_xcells) array of log₁₀-resistivity.

    The mapping follows the PW2D convention: each column code in a layer
    gives the number of mesh x-cells that column spans; n_merge gives
    the number of mesh z-rows the layer covers.
    """
    n_z = mesh.n_zcells   # 31
    n_x = mesh.n_xcells   # 576
    grid = np.full((n_z, n_x), np.nan)

    z_row = 0
    occm  = 0
    for layer in model.layers:
        n_merge = layer["n_merge"]
        codes   = layer["params"]
        fp_sta  = z_row
        fp_end  = min(z_row + n_merge, n_z)
        z_row   = fp_end

        x_col = 0
        for code in codes:
            cp_sta = x_col
            cp_end = min(x_col + code, n_x)
            x_col  = cp_end
            if occm < len(param_values):
                grid[fp_sta:fp_end, cp_sta:cp_end] = param_values[occm]
            occm += 1

    return grid


# -----------------------------------------------------------------------
# InversionResult
# -----------------------------------------------------------------------

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
        Log₁₀-resistivity on the mesh grid (NaN for cells outside
        the model domain).
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

        # --- iter and resp files ---
        self.iter_files = sorted(wd.glob("*.iter"), key=_iter_number)
        self.resp_files = sorted(wd.glob("*.resp"), key=_iter_number)

        # --- log file ---
        log_p = _scan_one(wd, is_log_file, "*.logfile")
        if log_p is None:
            log_p = _scan_one(wd, is_log_file, "LogFile*")
        if log_p:
            try:
                self.log = OccamLog.read(log_p)
            except Exception:
                pass

        # --- mesh file ---
        mesh_p = _scan_one(wd, is_mesh_file, "*esh*")
        if mesh_p is None:
            mesh_p = _scan_one(wd, is_mesh_file, "*")
        if mesh_p:
            try:
                self.mesh = OccamMesh.read(mesh_p)
            except Exception:
                pass

        # --- model file ---
        model_p = _scan_one(wd, is_model_file, "*odel*")
        if model_p is None:
            model_p = _scan_one(wd, is_model_file, "*")
        if model_p:
            try:
                self.model = OccamModel.read(model_p)
            except Exception:
                pass

        # --- select iteration ---
        if not self.iter_files:
            return

        if self._iteration is None:
            iter_p = self.iter_files[-1]
        else:
            wanted = self._iteration
            matches = [p for p in self.iter_files if _iter_number(p) == wanted]
            iter_p  = matches[0] if matches else self.iter_files[-1]

        try:
            self.best_iter = OccamIter.read(iter_p)
        except Exception:
            pass

        # --- matching response file ---
        iter_n = _iter_number(iter_p)
        resp_matches = [p for p in self.resp_files if _iter_number(p) == iter_n]
        resp_p = resp_matches[0] if resp_matches else (
            self.resp_files[-1] if self.resp_files else None
        )
        if resp_p:
            try:
                self.response = OccamResponse.read(resp_p)
            except Exception:
                pass

        # --- build rho_2d ---
        if self.best_iter is not None and self.model is not None and self.mesh is not None:
            try:
                self.rho_2d = _build_rho_2d(
                    self.best_iter.param_values, self.model, self.mesh
                )
            except Exception:
                pass

        if self.verbose:
            self.logger.info(
                "InversionResult: %d iter files, loaded iter %d from %s",
                len(self.iter_files), iter_n, wd,
            )

    # ------------------------------------------------------------------
    # iter2dat (Bo Yang's post-processing format)
    # ------------------------------------------------------------------
    def iter2dat(self, output_file: PathLike) -> Path:
        """Write an ``iter2dat``-style ASCII model file.

        Columns: ``x_center  z_center  log10_rho``

        x is centred so the profile midpoint is at x = 0.
        z is positive downward (depth in metres).

        Parameters
        ----------
        output_file : path-like
            Output path for the ``.dat`` file.

        Returns
        -------
        Path
            The written file path.

        Raises
        ------
        RuntimeError
            If the result has not been fully loaded (missing mesh, model,
            or iter file).
        """
        if self.rho_2d is None or self.mesh is None:
            raise RuntimeError(
                "Cannot write iter2dat: rho_2d or mesh not available. "
                "Ensure workdir contains mesh, model, and iter files."
            )

        out = Path(output_file)
        out.parent.mkdir(parents=True, exist_ok=True)

        # Cell centres
        xn = self.mesh.x_nodes
        zn = self.mesh.z_nodes
        x_centers = (xn[:-1] + xn[1:]) / 2.0
        z_centers = (zn[:-1] + zn[1:]) / 2.0

        # Centre x around the profile midpoint
        x_centers = x_centers - x_centers.mean()

        n_z, n_x = self.rho_2d.shape
        lines: list[str] = []
        for iz in range(n_z):
            for ix in range(n_x):
                val = self.rho_2d[iz, ix]
                if not np.isnan(val):
                    lines.append(
                        f"{x_centers[ix]:>15.4f}  {z_centers[iz]:>15.4f}  "
                        f"{val:>12.5f}\n"
                    )

        with out.open("w") as fh:
            fh.writelines(lines)

        if self.verbose:
            self.logger.info("iter2dat: %d cells written to %s", len(lines), out)
        return out

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
