# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Subprocess launcher for ModEM inversion and forward runs."""

from __future__ import annotations

import shlex
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Union

from .base import ModEmBase
from .config import ModEmConfig
from .doc import _modem_param_docs as _params
from .results import InversionResult

PathLike = Union[str, Path]

__all__ = ["ModEmRunner"]


def _resolve_binary(name: str, workdir: Path) -> Path | None:
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
    def __init__(
        self,
        workdir: PathLike,
        config: ModEmConfig | None = None,
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
        model: PathLike,
        data: PathLike,
        control: PathLike,
        covariance: PathLike | None = None,
        *,
        mode: str | None = None,
        use_mpi: bool | None = None,
        n_procs: int | None = None,
        extra_args: Sequence[str] | None = None,
        timeout: int | None = None,
        load_result: bool = True,
    ) -> InversionResult | None:
        """Run a nonlinear ModEM inversion subprocess.

        This method builds the command line for an inversion run,
        launches it with :func:`subprocess.run`, checks the process
        return code, and optionally scans ``workdir`` into an
        :class:`~pycsamt.models.modem.results.InversionResult`.
        The executable is called with the inversion flag sequence
        ``-I NLCG`` so the run uses ModEM's nonlinear conjugate
        gradient inversion mode.

        Parameters
        ----------
        model : path-like
            Starting or current ModEM model file passed to the
            runner. The path may be absolute or relative to
            ``workdir``. In 2-D runs it usually points to a
            ``.rho`` file; in 3-D runs it usually points to a
            ``.ws`` model file.
        data : path-like
            Observed-data file passed to ModEM. The file should
            match the dimensionality, component selection, sign
            convention, period list, and units expected by the
            selected executable and control settings.
        control : path-like
            ModEM inversion-control file passed to inversion
            runs. It contains iteration limits, target misfit,
            lambda controls, output stem, line-search settings,
            and related nonlinear solver parameters.
        covariance : path-like, optional
            ModEM covariance file passed to 3-D inversion runs.
            The file describes smoothing strengths, smoothing
            iteration count, and active-cell masks. It is
            commonly omitted for 2-D workflows.
        mode : {"2d", "3d"}, optional
            Dimensionality override for this invocation. If omitted,
            ``config.mode`` is used. The value selects
            ``config.binary_2d`` or ``config.binary_3d`` and should
            match the model and data file formats.
        use_mpi : bool, optional
            MPI override for this invocation. If omitted,
            ``config.use_mpi`` is used. When true, the command is
            prefixed by ``config.mpi_command -np <n_procs>``.
        n_procs : int, optional
            Number of MPI processes requested for this invocation.
            If omitted, ``config.n_procs`` is used. The value is
            ignored when MPI execution is disabled.
        extra_args : sequence of str, optional
            Additional command-line arguments appended to the
            ModEM executable invocation. Use this for advanced
            executable flags while keeping standard file
            handling under the runner.
        timeout : float, optional
            Maximum run time in seconds for the ModEM process.
            If the process exceeds this duration, it is
            terminated by the caller. Leave as ``None`` for no
            Python-side timeout.
        load_result : bool, default True
            Whether to scan ``workdir`` and return an
            :class:`InversionResult` after the executable
            finishes. Set to ``False`` when the caller only
            needs process completion.

        Returns
        -------
        InversionResult or None
            Parsed result object when ``load_result`` is true.
            Otherwise ``None`` is returned after the subprocess
            completes successfully.

        Raises
        ------
        FileNotFoundError
            Raised when the selected ModEM executable cannot be
            resolved from ``PATH`` or the known local build
            locations under ``workdir``.
        subprocess.CalledProcessError
            Raised when ModEM exits with a non-zero return code.
        subprocess.TimeoutExpired
            Raised when ``timeout`` is set and the process does not
            finish before the timeout expires.

        Notes
        -----
        Relative paths are passed to ModEM exactly as supplied while
        the subprocess working directory is set to ``workdir``. This
        mirrors the way ModEM control files and examples are usually
        written: model, data, control, and covariance file names are
        interpreted relative to the run directory.

        Examples
        --------
        Run a serial 3-D inversion and load the result:

        >>> from pycsamt.models.modem.config import ModEmConfig
        >>> from pycsamt.models.modem.runner import ModEmRunner
        >>> cfg = ModEmConfig(mode="3d", use_mpi=False)
        >>> runner = ModEmRunner("modem_run", config=cfg)
        >>> result = runner.run(
        ...     model="m0.ws",
        ...     data="d0.dat",
        ...     control="control.inv",
        ...     covariance="covariance.cov",
        ... )

        Run with MPI and defer result loading:

        >>> cfg = ModEmConfig(mode="3d", use_mpi=True, n_procs=8)
        >>> runner = ModEmRunner("modem_run", config=cfg)
        >>> runner.run(
        ...     "m0.ws",
        ...     "d0.dat",
        ...     "control.inv",
        ...     covariance="covariance.cov",
        ...     load_result=False,
        ... )

        See Also
        --------
        command
            Build the same inversion command without executing it.
        run_forward
            Execute a forward-only ModEM response calculation.
        InversionResult
            Loader for logs, models, data, and covariance outputs.
        """
        cfg = self.config
        _mode = (mode or cfg.mode).lower()
        _mpi = cfg.use_mpi if use_mpi is None else use_mpi
        _procs = cfg.n_procs if n_procs is None else n_procs

        bin_name = cfg.binary_3d if _mode == "3d" else cfg.binary_2d
        binary = _resolve_binary(bin_name, self.workdir)
        if binary is None:
            msg = (
                f"ModEM binary '{bin_name}' not found. "
                "Build it first; see the ModEM source README."
            )
            raise FileNotFoundError(msg)

        cmd: list[str] = []
        if _mpi:
            cmd += [cfg.mpi_command, "-np", str(_procs)]
        cmd += [
            str(binary),
            "-I",
            "NLCG",
            str(model),
            str(data),
            str(control),
        ]
        if covariance is not None:
            cmd.append(str(covariance))
        if extra_args:
            cmd.extend(extra_args)

        if self.verbose:
            self.logger.info(
                "ModEmRunner: %s",
                " ".join(shlex.quote(c) for c in cmd),
            )

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
        model: PathLike,
        data: PathLike,
        *,
        mode: str | None = None,
        use_mpi: bool | None = None,
        n_procs: int | None = None,
        timeout: int | None = None,
        load_result: bool = True,
    ) -> InversionResult | None:
        """Run a forward-only ModEM response calculation.

        Forward mode evaluates predicted responses for an
        existing model and data file without updating the
        model. The executable is called with the ``-F`` flag.
        This is useful for checking a starting model, computing
        synthetic responses, or re-evaluating a final model
        after editing data weights.

        Parameters
        ----------
        model : path-like
            Starting, current, or final ModEM model file used
            for response calculation. The path may be absolute
            or relative to ``workdir``.
        data : path-like
            Data file that defines the stations, periods,
            components, and errors for which responses are
            calculated.
        mode : {"2d", "3d"}, optional
            Dimensionality override for this invocation. If omitted,
            ``config.mode`` selects the executable.
        use_mpi : bool, optional
            MPI override for this invocation. If omitted,
            ``config.use_mpi`` is used.
        n_procs : int, optional
            Number of MPI processes requested when MPI execution is
            enabled. If omitted, ``config.n_procs`` is used.
        timeout : float, optional
            Maximum run time in seconds for the ModEM process.
            Leave as ``None`` for no Python-side timeout.
        load_result : bool, default True
            Whether to scan ``workdir`` and return an
            :class:`InversionResult` after the executable
            finishes.

        Returns
        -------
        InversionResult or None
            Parsed result object when ``load_result`` is true.
            Otherwise ``None`` is returned after successful process
            completion.

        Raises
        ------
        FileNotFoundError
            Raised when the selected ModEM executable cannot be
            resolved.
        subprocess.CalledProcessError
            Raised when the forward executable exits with a non-zero
            return code.
        subprocess.TimeoutExpired
            Raised when the process exceeds ``timeout`` seconds.

        Examples
        --------
        >>> from pycsamt.models.modem.runner import ModEmRunner
        >>> runner = ModEmRunner("modem_run")
        >>> result = runner.run_forward("mi.ws", "d0.dat")
        """
        cfg = self.config
        _mode = (mode or cfg.mode).lower()
        _mpi = cfg.use_mpi if use_mpi is None else use_mpi
        _procs = cfg.n_procs if n_procs is None else n_procs

        bin_name = cfg.binary_3d if _mode == "3d" else cfg.binary_2d
        binary = _resolve_binary(bin_name, self.workdir)
        if binary is None:
            msg = f"ModEM binary '{bin_name}' not found."
            raise FileNotFoundError(msg)

        cmd: list[str] = []
        if _mpi:
            cmd += [cfg.mpi_command, "-np", str(_procs)]
        cmd += [str(binary), "-F", str(model), str(data)]

        if self.verbose:
            self.logger.info(
                "ModEmRunner.forward: %s",
                " ".join(shlex.quote(c) for c in cmd),
            )

        proc = subprocess.run(
            cmd,
            cwd=str(self.workdir),
            timeout=timeout,
        )
        proc.check_returncode()

        if load_result:
            return InversionResult(self.workdir, config=cfg)
        return None

    def command(
        self,
        model: PathLike,
        data: PathLike,
        control: PathLike,
        covariance: PathLike | None = None,
        *,
        mode: str | None = None,
        use_mpi: bool | None = None,
        n_procs: int | None = None,
    ) -> str:
        """Return the inversion command without executing ModEM.

        This helper is a dry-run view of :meth:`run`. It applies the
        same mode, MPI, process-count, and file-name choices, but it
        does not resolve the executable or start a subprocess.

        Parameters
        ----------
        model : path-like
            Starting or current ModEM model file to include in
            the command string.
        data : path-like
            Observed-data file to include in the command string.
        control : path-like
            Inversion-control file to include in the command
            string.
        covariance : path-like, optional
            Covariance file appended to the command when
            supplied.
        mode : {"2d", "3d"}, optional
            Dimensionality override for the command string.
        use_mpi : bool, optional
            MPI override for the command string.
        n_procs : int, optional
            MPI process-count override.

        Returns
        -------
        str
            Shell-quoted command string suitable for display,
            logging, or copying into a terminal.

        Examples
        --------
        >>> from pycsamt.models.modem.config import ModEmConfig
        >>> from pycsamt.models.modem.runner import ModEmRunner
        >>> cfg = ModEmConfig(mode="3d", use_mpi=True, n_procs=4)
        >>> runner = ModEmRunner("modem_run", config=cfg)
        >>> "mpirun" in runner.command("m0.ws", "d0.dat", "c.inv")
        True
        """
        cfg = self.config
        _mode = (mode or cfg.mode).lower()
        _mpi = cfg.use_mpi if use_mpi is None else use_mpi
        _procs = cfg.n_procs if n_procs is None else n_procs

        bin_name = cfg.binary_3d if _mode == "3d" else cfg.binary_2d
        cmd: list[str] = []
        if _mpi:
            cmd += [cfg.mpi_command, "-np", str(_procs)]
        cmd += [
            bin_name,
            "-I",
            "NLCG",
            str(model),
            str(data),
            str(control),
        ]
        if covariance is not None:
            cmd.append(str(covariance))
        return " ".join(shlex.quote(c) for c in cmd)


ModEmRunner.__doc__ = rf"""
Launch ModEM inversion and forward-modeling subprocesses.

``ModEmRunner`` is the execution layer of the ModEM wrapper. It
does not build input files itself; instead it receives model,
data, control, and optional covariance files created by
:class:`~pycsamt.models.modem.builder.InputBuilder` or by user
code, selects the configured ModEM executable, launches the
process from ``workdir``, and can load the finished run into an
:class:`~pycsamt.models.modem.results.InversionResult`.

For inversion runs the command has the logical form

.. code-block:: text

   Mod3DMT -I NLCG model.ws data.dat control.inv covariance.cov

or, when MPI execution is requested,

.. code-block:: text

   mpirun -np 8 Mod3DMT -I NLCG model.ws data.dat control.inv

The 2-D runner uses ``config.binary_2d`` and 2-D model files,
whereas the 3-D runner uses ``config.binary_3d`` and can pass a
covariance file. The forward-only path uses the ``-F`` flag and
requires only model and data files.

Parameters
----------
{_params.common.workdir}
{_params.common.config}
**kwargs : dict
    Additional keyword arguments passed to
    :class:`~pycsamt.models.modem.base.ModEmBase`, including
    logging and verbosity options.

Attributes
----------
workdir : pathlib.Path
    Directory from which the subprocess is launched. Relative
    file names in ModEM command lines are interpreted from this
    directory.
config : ModEmConfig
    Configuration used to choose mode, executable names, MPI
    command, process count, and result-loading behavior.

Notes
-----
Executable resolution is intentionally small and predictable.
The runner first accepts a name found on ``PATH``. It then
checks ``workdir / name`` and local source-build locations under
``workdir / "_source" / "3D"`` and ``workdir / "_source" /
"2D"``. Use an absolute executable path in the configuration
when running against a system installation outside the run
folder.

The runner delegates process execution to :func:`subprocess.run`.
Standard output and standard error are inherited from the
calling process unless the surrounding application redirects
them. A non-zero exit status is converted to
:class:`subprocess.CalledProcessError` by ``check_returncode``.

Examples
--------
Create a runner for a serial 3-D inversion:

>>> from pycsamt.models.modem.config import ModEmConfig
>>> from pycsamt.models.modem.runner import ModEmRunner
>>> cfg = ModEmConfig(mode="3d", use_mpi=False)
>>> runner = ModEmRunner("modem_run", config=cfg)
>>> cmd = runner.command("m0.ws", "d0.dat", "control.inv")
>>> "-I NLCG" in cmd
True

Run a forward response calculation from an existing model:

>>> runner = ModEmRunner("modem_run", config=cfg)
>>> result = runner.run_forward("mi.ws", "d0.dat")

See Also
--------
InputBuilder
    Build ModEM model, data, covariance, and control files.
ModEmConfig
    Store executable names, MPI settings, and inversion
    options shared by the runner.
InversionResult
    Load logs, models, data, and covariance after a run.

References
----------
.. [1] Kelbert, A., Meqbel, N., Egbert, G. D., and
   Tandon, K. (2014). ModEM: A modular system for inversion
   of electromagnetic geophysical data. *Computers &
   Geosciences*, 66, 40-53.
   doi:10.1016/j.cageo.2014.01.010.
"""
