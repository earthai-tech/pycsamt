# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEmRunner — launch a ModEM inversion subprocess.

Supports serial and MPI execution for both 2-D and 3-D runs.

Command pattern
---------------
Serial 3D::

    Mod3DMT -I NLCG <model.ws> <data.dat> <control.inv> [covariance.cov]

MPI 3D::

    mpirun -np <n_procs> Mod3DMT -I NLCG <model.ws> <data.dat> <control.inv>

The binary is resolved in this order:

1. ``config.binary_3d`` / ``config.binary_2d`` (if it is an absolute path or on PATH)
2. ``workdir / config.binary_3d``
3. ``workdir / '_source' / '3D' / config.binary_3d``  (compiled in-tree)

Quick start
-----------
>>> from pycsamt.models.modem.runner import ModEmRunner
>>> runner = ModEmRunner(workdir="./run", config=cfg)
>>> result = runner.run(model="m0.ws", data="d0.dat", control="block2.inv")
"""

from __future__ import annotations

import shlex
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Union

from .base    import ModEmBase
from .config  import ModEmConfig
from .results import InversionResult

PathLike = Union[str, Path]

__all__ = ["ModEmRunner"]


def _resolve_binary(name: str, workdir: Path) -> Optional[Path]:
    """Return the resolved path to *name* or None if not found."""
    if shutil.which(name):
        return Path(name)
    candidates = [
        workdir / name,
        workdir / "_source" / "3D" / name,
        workdir / "_source" / "2D" / name,
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


class ModEmRunner(ModEmBase):
    """Launch and optionally monitor a ModEM inversion.

    Parameters
    ----------
    workdir : path-like
        Directory from which the binary is invoked (all relative paths
        in the command line are resolved here).
    config : ModEmConfig, optional
    """

    def __init__(
        self,
        workdir: PathLike,
        config: Optional[ModEmConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.workdir: Path = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.config: ModEmConfig = config or ModEmConfig()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        model:       PathLike,
        data:        PathLike,
        control:     PathLike,
        covariance:  Optional[PathLike] = None,
        *,
        mode:        Optional[str]  = None,
        use_mpi:     Optional[bool] = None,
        n_procs:     Optional[int]  = None,
        extra_args:  Optional[List[str]] = None,
        timeout:     Optional[int]  = None,
        load_result: bool = True,
    ) -> Optional[InversionResult]:
        """Run ModEM and optionally return an :class:`InversionResult`.

        Parameters
        ----------
        model : path-like
            Starting model file (relative to ``workdir``).
        data : path-like
            Observed-data file (relative to ``workdir``).
        control : path-like
            Control/inversion-parameters file (relative to ``workdir``).
        covariance : path-like, optional
            Covariance file (3-D only; relative to ``workdir``).
        mode : {'2d', '3d'}, optional
            Override ``config.mode``.
        use_mpi : bool, optional
            Override ``config.use_mpi``.
        n_procs : int, optional
            Override ``config.n_procs``.
        extra_args : list of str, optional
            Additional command-line arguments appended verbatim.
        timeout : int, optional
            Seconds to wait before killing the process (``None`` = unlimited).
        load_result : bool
            If ``True`` (default) load and return an :class:`InversionResult`
            from ``workdir`` after the run finishes.

        Returns
        -------
        InversionResult or None
        """
        cfg    = self.config
        _mode  = (mode or cfg.mode).lower()
        _mpi   = cfg.use_mpi if use_mpi is None else use_mpi
        _procs = cfg.n_procs if n_procs is None else n_procs

        bin_name = cfg.binary_3d if _mode == "3d" else cfg.binary_2d
        binary   = _resolve_binary(bin_name, self.workdir)
        if binary is None:
            raise FileNotFoundError(
                f"ModEM binary '{bin_name}' not found. "
                "Build it first (see pycsamt/models/modem/_source/README.txt)."
            )

        cmd: List[str] = []
        if _mpi:
            cmd += [cfg.mpi_command, "-np", str(_procs)]
        cmd += [str(binary), "-I", "NLCG",
                str(model), str(data), str(control)]
        if covariance is not None:
            cmd.append(str(covariance))
        if extra_args:
            cmd.extend(extra_args)

        if self.verbose:
            self.logger.info("ModEmRunner: %s", " ".join(shlex.quote(c) for c in cmd))

        proc = subprocess.run(
            cmd,
            cwd=str(self.workdir),
            timeout=timeout,
        )
        proc.check_returncode()

        if load_result:
            return InversionResult(self.workdir, config=cfg)
        return None

    # ------------------------------------------------------------------
    # Forward-only run
    # ------------------------------------------------------------------

    def run_forward(
        self,
        model:      PathLike,
        data:       PathLike,
        *,
        mode:       Optional[str]  = None,
        use_mpi:    Optional[bool] = None,
        n_procs:    Optional[int]  = None,
        timeout:    Optional[int]  = None,
        load_result: bool = True,
    ) -> Optional[InversionResult]:
        """Run a single forward calculation (``-F`` flag).

        Parameters and return value mirror :meth:`run`.
        """
        cfg    = self.config
        _mode  = (mode or cfg.mode).lower()
        _mpi   = cfg.use_mpi if use_mpi is None else use_mpi
        _procs = cfg.n_procs if n_procs is None else n_procs

        bin_name = cfg.binary_3d if _mode == "3d" else cfg.binary_2d
        binary   = _resolve_binary(bin_name, self.workdir)
        if binary is None:
            raise FileNotFoundError(
                f"ModEM binary '{bin_name}' not found."
            )

        cmd: List[str] = []
        if _mpi:
            cmd += [cfg.mpi_command, "-np", str(_procs)]
        cmd += [str(binary), "-F", str(model), str(data)]

        if self.verbose:
            self.logger.info("ModEmRunner.forward: %s", " ".join(shlex.quote(c) for c in cmd))

        proc = subprocess.run(cmd, cwd=str(self.workdir), timeout=timeout)
        proc.check_returncode()

        if load_result:
            return InversionResult(self.workdir, config=cfg)
        return None

    def command(
        self,
        model:      PathLike,
        data:       PathLike,
        control:    PathLike,
        covariance: Optional[PathLike] = None,
        *,
        mode:       Optional[str]  = None,
        use_mpi:    Optional[bool] = None,
        n_procs:    Optional[int]  = None,
    ) -> str:
        """Return the shell command string without executing it."""
        cfg    = self.config
        _mode  = (mode or cfg.mode).lower()
        _mpi   = cfg.use_mpi if use_mpi is None else use_mpi
        _procs = cfg.n_procs if n_procs is None else n_procs

        bin_name = cfg.binary_3d if _mode == "3d" else cfg.binary_2d
        cmd: List[str] = []
        if _mpi:
            cmd += [cfg.mpi_command, "-np", str(_procs)]
        cmd += [bin_name, "-I", "NLCG", str(model), str(data), str(control)]
        if covariance is not None:
            cmd.append(str(covariance))
        return " ".join(shlex.quote(c) for c in cmd)
