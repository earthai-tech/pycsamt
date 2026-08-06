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
from typing import Union

import numpy as np

from .base import OccamBase
from .data import OccamData
from .log import OccamLog
from .mesh import OccamMesh
from .model import OccamModel
from .response import OccamResponse
from .startup import OccamIter
from .validation import (
    is_data_file,
    is_log_file,
    is_mesh_file,
    is_model_file,
)

PathLike = Union[str, Path]

__all__ = ["InversionResult"]


def _iter_number(path: Path) -> int:
    """Extract the leading or embedded integer from a filename stem."""
    m = re.search(r"\d+", path.stem)
    return int(m.group()) if m else -1


def _scan_one(wd: Path, predicate, glob_pat: str) -> Path | None:
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

    The PW2D mapping uses each column code as the number of
    mesh x-cells spanned. ``n_merge`` gives the number of
    mesh z-rows covered by a model layer.
    """
    n_z = mesh.n_zcells  # 31
    n_x = mesh.n_xcells  # 576
    grid = np.full((n_z, n_x), np.nan)

    z_row = 0
    occm = 0
    for layer in model.layers:
        n_merge = layer["n_merge"]
        codes = layer["params"]
        fp_sta = z_row
        fp_end = min(z_row + n_merge, n_z)
        z_row = fp_end

        x_col = 0
        for code in codes:
            cp_sta = x_col
            cp_end = min(x_col + code, n_x)
            x_col = cp_end
            if occm < len(param_values):
                grid[fp_sta:fp_end, cp_sta:cp_end] = param_values[occm]
            occm += 1

    return grid


# -----------------------------------------------------------------------
# InversionResult
# -----------------------------------------------------------------------


class InversionResult(OccamBase):
    r"""Load and summarize a completed Occam2D inversion run.

    ``InversionResult`` is the post-processing access layer
    for an Occam2D working directory. It scans the directory
    for input and output files, loads the selected iteration,
    matches the response file, reads the log, and builds a
    two-dimensional log10-resistivity grid on the mesh.

    The reconstruction maps the iteration parameter vector to
    model layers and mesh cells. If :math:`m_j` is the
    log10-resistivity value assigned to model parameter
    :math:`j`, then the physical resistivity represented by a
    cell in that parameter group is

    .. math::

        \rho_j = 10^{m_j}.

    The stored :attr:`rho_2d` array keeps :math:`m_j`, not
    :math:`\rho_j`, because Occam iteration files store
    log10-resistivity values.

    Parameters
    ----------
    workdir : path-like, default "."
        Directory produced by an Occam2D run. It should
        contain an ``Occam2DMesh`` file, an ``Occam2DModel``
        file, a data file, a log file, ``.iter`` files, and
        matching ``.resp`` files.
    iteration : int or None, default None
        Iteration number to load. If ``None``, the highest
        numbered ``.iter`` file is selected. If the requested
        iteration is unavailable, the loader falls back to the
        last available iteration file.
    verbose : int or bool, default 0
        Verbosity level inherited from :class:`OccamBase`.
        Positive values enable progress messages through the
        instance logger.
    logger : logging.Logger, optional
        Logger used for progress and diagnostic messages. If
        omitted, a class-specific PyCSAMT logger is created.

    Attributes
    ----------
    workdir : pathlib.Path
        Directory scanned by the loader.
    log : OccamLog or None
        Parsed convergence log when a log file is found.
    mesh : OccamMesh or None
        Parsed finite-element mesh.
    model : OccamModel or None
        Parsed model-parameter definition.
    best_iter : OccamIter or None
        Selected iteration file. The name is kept for backward
        compatibility; it may be the requested iteration or
        last available iteration.
    response : OccamResponse or None
        Response file matching the selected iteration when
        available. If no exact match exists, the loader falls
        back to the last response file.
    data : OccamData or None
        Parsed observed-data file when available.
    iter_files : list of pathlib.Path
        All ``.iter`` files found in ``workdir``, sorted by
        embedded iteration number.
    resp_files : list of pathlib.Path
        All ``.resp`` files found in ``workdir``, sorted by
        embedded iteration number.
    rho_2d : numpy.ndarray of float or None
        Log10-resistivity grid with shape
        ``(mesh.n_zcells, mesh.n_xcells)``. Cells outside the
        model domain are stored as ``nan``.

    Notes
    -----
    The loader is deliberately tolerant. Missing optional
    files leave corresponding attributes as ``None`` instead
    of failing immediately. A missing working directory still
    raises :class:`NotADirectoryError` because there is no
    useful scan to perform.

    See Also
    --------
    OccamRunner
        Runs the executable that produces result files.
    OccamLog
        Parses convergence information loaded here.
    OccamResponse
        Parses modeled responses and weighted residuals.
    PlotModel
        Visualizes the reconstructed model grid.

    Examples
    --------
    Load the latest available iteration:

    >>> from pycsamt.models.occam2d import InversionResult
    >>> result = InversionResult(workdir="occam_run")
    >>> result.final_rms

    Load a specific iteration and export the model grid:

    >>> result = InversionResult("occam_run", iteration=17)
    >>> result.iter2dat("occam_run/final_model.dat")

    Access loaded components directly:

    >>> result.mesh.n_xcells
    >>> result.response.rms

    References
    ----------
    .. [InversionResult-1] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    .. [InversionResult-2] Constable, S. C., Parker, R. L., and Constable,
       C. G., "Occam's inversion: A practical algorithm for
       generating smooth models from electromagnetic sounding
       data", Geophysics, 52(3), 289-300, 1987.
    """

    def __init__(
        self,
        workdir: PathLike = ".",
        iteration: int | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.workdir = Path(workdir)
        self._iteration = iteration

        self.log: OccamLog | None = None
        self.mesh: OccamMesh | None = None
        self.model: OccamModel | None = None
        self.best_iter: OccamIter | None = None
        self.response: OccamResponse | None = None
        self.data: OccamData | None = None

        self.iter_files: list[Path] = []
        self.resp_files: list[Path] = []
        self.rho_2d: np.ndarray | None = None

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
            iter_p = matches[0] if matches else self.iter_files[-1]

        try:
            self.best_iter = OccamIter.read(iter_p)
        except Exception:
            pass

        # --- matching response file ---
        iter_n = _iter_number(iter_p)
        resp_matches = [
            p for p in self.resp_files if _iter_number(p) == iter_n
        ]
        resp_p = (
            resp_matches[0]
            if resp_matches
            else (self.resp_files[-1] if self.resp_files else None)
        )
        if resp_p:
            try:
                self.response = OccamResponse.read(resp_p)
            except Exception:
                pass

        # --- data file ---
        data_p = None
        if self.best_iter is not None and getattr(
            self.best_iter, "data_file", None
        ):
            candidate = wd / self.best_iter.data_file
            if candidate.exists():
                data_p = candidate
        if data_p is None:
            data_p = _scan_one(wd, is_data_file, "*.dat")
        if data_p:
            try:
                self.data = OccamData.read(data_p)
            except Exception:
                pass

        # --- build rho_2d ---
        if (
            self.best_iter is not None
            and self.model is not None
            and self.mesh is not None
        ):
            try:
                self.rho_2d = _build_rho_2d(
                    self.best_iter.param_values, self.model, self.mesh
                )
            except Exception:
                pass

        if self.verbose:
            self.logger.info(
                "InversionResult: %d iter files, loaded iter %d from %s",
                len(self.iter_files),
                iter_n,
                wd,
            )

    # ------------------------------------------------------------------
    # iter2dat (Bo Yang's post-processing format)
    # ------------------------------------------------------------------
    def iter2dat(self, output_file: PathLike) -> Path:
        r"""Write the selected model as a three-column ASCII file.

        The exported file contains one row for each finite
        cell in :attr:`rho_2d`. Columns are ``x_center``,
        ``z_center``, and ``log10_rho``. Horizontal
        coordinates are centered around the profile midpoint
        and depths are positive downward.

        This format is useful for external plotting tools and
        for workflows that expect Bo Yang-style ``iter2dat``
        output. The values remain in log10-resistivity units:

        .. math::

            m = \log_{10}(\rho).

        Parameters
        ----------
        output_file : path-like
            Destination path for the exported ASCII model.
            Parent directories are created when needed.

        Returns
        -------
        pathlib.Path
            Path to the file that was written.

        Raises
        ------
        RuntimeError
            Raised when the result is not fully loaded and the
            mesh or reconstructed grid is unavailable.

        See Also
        --------
        InversionResult.rho_2d
            Grid used to generate the exported values.
        PlotModel
            Visualizes the same reconstructed model grid.

        Examples
        --------
        >>> from pycsamt.models.occam2d import InversionResult
        >>> result = InversionResult("occam_run")
        >>> out = result.iter2dat("occam_run/final_model.dat")
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
            self.logger.info(
                "iter2dat: %d cells written to %s", len(lines), out
            )
        return out

    # ------------------------------------------------------------------
    # Plotting (delegates to plot.py)
    # ------------------------------------------------------------------
    def plot_model(self, **kwargs):
        """Plot the reconstructed 2-D resistivity model."""
        from .plot import PlotModel

        return PlotModel(result=self, **kwargs).plot()

    def plot_response(self, **kwargs):
        """Plot observed and modeled response curves."""
        from .plot import PlotResponse

        return PlotResponse(result=self, **kwargs).plot()

    def plot_misfit(self, **kwargs):
        """Plot RMS misfit as a function of iteration."""
        from .plot import PlotMisfit

        return PlotMisfit(result=self, **kwargs).plot()

    def plot_pseudo(self, **kwargs):
        """Plot an observed-data pseudosection."""
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
        """Number of iteration files discovered in ``workdir``."""
        return len(self.iter_files)

    def summary(self) -> str:
        """Return a short text summary of the loaded inversion."""
        return (
            f"InversionResult\n"
            f"  workdir    : {self.workdir}\n"
            f"  iterations : {self.n_iterations}\n"
            f"  final RMS  : {self.final_rms:.4f}\n"
            f"  converged  : {self.log.converged if self.log else 'N/A'}\n"
        )
