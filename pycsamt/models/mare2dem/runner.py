# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Subprocess launcher for MARE2DEM inversion runs."""

from __future__ import annotations

import shlex
import shutil
import subprocess
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Union

from .._process import run_streamed
from .base import Mare2DEMBase
from .config import Mare2DEMConfig
from .doc import _mare2dem_param_docs as _params
from .source import SourceManager

PathLike = Union[str, Path]

__all__ = ["Mare2DEMRunner"]


def _resolve_binary(
    name: str,
    source_mgr: SourceManager,
) -> Path | None:
    """Return the MARE2DEM executable path or ``None``.

    Resolution order:

    1. Binary name on ``PATH``.
    2. Binary resolved by :meth:`SourceManager.resolve_binary`.
    """
    if shutil.which(name):
        return Path(name)
    return source_mgr.resolve_binary()


def _is_wsl(binary) -> bool:
    return isinstance(binary, str) and binary.startswith("wsl:")


def _run_stem(resistivity_stem) -> str:
    """MARE2DEM's command-line argument for a ``.resistivity`` file.

    Only a trailing ``.resistivity`` is removed: the iteration number is
    part of the stem (``mare2dem.0``), and ``Path.with_suffix("")`` would
    strip it as if ``.0`` were an extension -- MARE2DEM then stops with
    "no iteration number in resistivity file" (and exits 0).
    """
    name = Path(str(resistivity_stem)).name
    if name.lower().endswith(".resistivity"):
        name = name[: -len(".resistivity")]
    return name


def _wsl_command(cfg, workdir, resistivity_stem, use_mpi, n_procs,
                 extra_args) -> list[str]:
    """``wsl -e bash -lc ...`` for a binary built inside WSL2.

    MARE2DEM cannot be built as a native Windows program; the desktop
    Solver Builder builds it in WSL and registers it as
    ``wsl:/home/<user>/.local/share/pycsamt/mare2dem/build/MARE2DEM``.
    The Windows work directory is reached through ``/mnt/<drive>``, and the
    pycsamt-managed toolchain (for ``mpirun``) is put on ``PATH``.
    """
    from pycsamt.models.solver_build import to_wsl_path

    binary = cfg.binary[len("wsl:"):]
    stem = _run_stem(resistivity_stem)
    parts = []
    if use_mpi:
        parts += [cfg.mpi_command, "-np", str(n_procs)]
    parts += [binary, stem, *(extra_args or [])]
    inner = (
        'TC="${XDG_DATA_HOME:-$HOME/.local/share}/pycsamt/toolchain/'
        'mare2dem"; [ -d "$TC/bin" ] && export PATH="$TC/bin:$PATH"; '
        "export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 "
        "GFORTRAN_UNBUFFERED_ALL=y; "
        f"cd {shlex.quote(to_wsl_path(Path(workdir).resolve()))} && "
        + " ".join(shlex.quote(p) for p in parts)
    )
    return ["wsl", "-e", "bash", "-lc", inner]


class Mare2DEMRunner(Mare2DEMBase):
    """Launch MARE2DEM inversion subprocesses.

    Parameters
    ----------
    workdir : path-like
        Directory from which MARE2DEM is executed. The resistivity
        file stem passed to the binary is resolved relative to
        this directory.
    config : Mare2DEMConfig, optional
        Configuration supplying binary name, MPI settings, and
        source directory. A default :class:`Mare2DEMConfig` is used
        when omitted.
    **kwargs :
        Forwarded to :class:`Mare2DEMBase`.
    """

    def __init__(
        self,
        workdir: PathLike,
        config: Mare2DEMConfig | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.workdir: Path = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.config: Mare2DEMConfig = config or Mare2DEMConfig()
        self._source_mgr = SourceManager(config=self.config)

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(
        self,
        resistivity_stem: PathLike,
        *,
        use_mpi: bool | None = None,
        n_procs: int | None = None,
        extra_args: Sequence[str] | None = None,
        timeout: int | None = None,
        load_result: bool = True,
        on_output: Callable[[str], None] | None = None,
        cancel: Callable[[], bool] | None = None,
    ) -> InversionResult | None:
        """Run a MARE2DEM inversion subprocess.

        MARE2DEM receives one positional argument: the stem of the
        ``.resistivity`` file. It derives the data filename by
        replacing the extension with ``.emdata`` and the settings
        filename with ``.settings``.

        Parameters
        ----------
        resistivity_stem : path-like
            Stem or full path to the ``.resistivity`` file. MARE2DEM
            strips the extension itself; you may pass either
            ``"run"`` or ``"run.resistivity"``. Relative paths are
            interpreted from ``workdir``.
        use_mpi : bool, optional
            MPI override. Falls back to ``config.use_mpi``.
        n_procs : int, optional
            Number of MPI processes. Falls back to ``config.n_procs``.
        extra_args : sequence of str, optional
            Additional command-line arguments appended to the
            MARE2DEM invocation.
        timeout : float, optional
            Maximum run time in seconds. ``None`` means no timeout.
        load_result : bool, default True
            Whether to scan ``workdir`` and return an
            :class:`~pycsamt.models.mare2dem.results.InversionResult`
            after the run completes.
        on_output : callable, optional
            Called with every console line while MARE2DEM runs. Supplying
            it, or ``cancel``, runs the solver through
            :func:`pycsamt.models._process.run_streamed`.
        cancel : callable returning bool, optional
            Polled while the solver runs; returning ``True`` stops the
            process tree (``mpirun`` and its ranks) and raises
            :class:`~pycsamt.models._process.ProcessCancelled`.

        Returns
        -------
        InversionResult or None
            Parsed result object when ``load_result`` is ``True``.
            ``None`` otherwise.

        Raises
        ------
        FileNotFoundError
            When the MARE2DEM binary cannot be located. Build it
            first with :class:`SourceManager`.
        subprocess.CalledProcessError
            When MARE2DEM exits with a non-zero return code.
        subprocess.TimeoutExpired
            When ``timeout`` is set and the process exceeds it.

        Examples
        --------
        Serial run (special single-process build):

        >>> from pycsamt.models.mare2dem import Mare2DEMConfig, Mare2DEMRunner
        >>> cfg = Mare2DEMConfig(use_mpi=False)
        >>> runner = Mare2DEMRunner("./mare2dem_run", config=cfg)
        >>> result = runner.run("mare2dem")

        MPI run with 8 processes:

        >>> cfg = Mare2DEMConfig(use_mpi=True, n_procs=8)
        >>> runner = Mare2DEMRunner("./mare2dem_run", config=cfg)
        >>> result = runner.run("mare2dem")
        """
        cfg = self.config
        _mpi = cfg.use_mpi if use_mpi is None else use_mpi
        _procs = cfg.n_procs if n_procs is None else n_procs

        if _is_wsl(cfg.binary):
            cmd = _wsl_command(cfg, self.workdir, resistivity_stem, _mpi,
                               _procs, extra_args)
            if self.verbose:
                self.logger.info("Mare2DEMRunner (WSL): %s", cmd[-1])
            self._execute(cmd, None, timeout, on_output, cancel)
            if load_result:
                from .results import InversionResult

                return InversionResult(self.workdir, config=cfg)
            return None

        binary = _resolve_binary(cfg.binary, self._source_mgr)
        if binary is None:
            raise FileNotFoundError(
                f"MARE2DEM binary '{cfg.binary}' not found on PATH or in "
                "the source directory. "
                "Build it first:  SourceManager().download(); "
                "SourceManager().build()"
            )

        stem = _run_stem(resistivity_stem)

        cmd: list[str] = []
        if _mpi:
            cmd += [cfg.mpi_command, "-np", str(_procs)]
        cmd += [str(binary), stem]
        if extra_args:
            cmd.extend(extra_args)

        if self.verbose:
            self.logger.info(
                "Mare2DEMRunner: %s",
                " ".join(shlex.quote(c) for c in cmd),
            )

        self._execute(cmd, self.workdir, timeout, on_output, cancel)

        if load_result:
            from .results import InversionResult

            return InversionResult(self.workdir, config=cfg)
        return None

    @staticmethod
    def _execute(cmd, cwd, timeout, on_output, cancel) -> None:
        """Run *cmd*; stream it when a console/cancel hook is given."""
        if on_output is not None or cancel is not None:
            code = run_streamed(cmd, cwd=cwd, on_output=on_output,
                                cancel=cancel, timeout=timeout)
            if code:
                raise subprocess.CalledProcessError(code, cmd)
            return
        proc = subprocess.run(
            cmd,
            cwd=None if cwd is None else str(cwd),
            timeout=timeout,
        )
        proc.check_returncode()

    # ------------------------------------------------------------------
    # Dry-run helper
    # ------------------------------------------------------------------

    def command(
        self,
        resistivity_stem: PathLike,
        *,
        use_mpi: bool | None = None,
        n_procs: int | None = None,
    ) -> str:
        """Return the MARE2DEM command string without executing it.

        Parameters
        ----------
        resistivity_stem : path-like
            Resistivity file stem passed to MARE2DEM.
        use_mpi : bool, optional
            MPI override.
        n_procs : int, optional
            Process-count override.

        Returns
        -------
        str
            Shell-quoted command string for display or logging.

        Examples
        --------
        >>> from pycsamt.models.mare2dem import Mare2DEMConfig, Mare2DEMRunner
        >>> cfg = Mare2DEMConfig(use_mpi=True, n_procs=4)
        >>> runner = Mare2DEMRunner("./run", config=cfg)
        >>> "mpirun" in runner.command("mare2dem")
        True
        """
        cfg = self.config
        _mpi = cfg.use_mpi if use_mpi is None else use_mpi
        _procs = cfg.n_procs if n_procs is None else n_procs
        stem = _run_stem(resistivity_stem)

        cmd: list[str] = []
        if _mpi:
            cmd += [cfg.mpi_command, "-np", str(_procs)]
        cmd += [cfg.binary, stem]
        return " ".join(shlex.quote(c) for c in cmd)


Mare2DEMRunner.__doc__ = rf"""
Launch MARE2DEM inversion subprocesses.

``Mare2DEMRunner`` is the execution layer of the MARE2DEM wrapper.
It receives the stem of a ``.resistivity`` file prepared by
:class:`~pycsamt.models.mare2dem.builder.InputBuilder`, selects
the configured MARE2DEM executable, and launches the MPI process
from ``workdir``. After the run it optionally loads output into
an :class:`~pycsamt.models.mare2dem.results.InversionResult`.

The command has the logical form:

.. code-block:: text

    mpirun -np 8 MARE2DEM mare2dem

where ``mare2dem`` is the stem of ``mare2dem.resistivity``.

Parameters
----------
{_params.common.workdir}
{_params.common.config}
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
workdir : pathlib.Path
    Directory from which the subprocess is launched.
config : Mare2DEMConfig
    Configuration used for binary name, MPI, and source
    management.

Notes
-----
Binary resolution is performed by :func:`_resolve_binary`, which
first checks ``PATH``, then delegates to
:meth:`SourceManager.resolve_binary` for locally compiled
binaries.

See Also
--------
SourceManager
    Download and compile the MARE2DEM binary.
InputBuilder
    Write the resistivity model, data, and settings files.
InversionResult
    Load MARE2DEM output after the run.

References
----------
.. [Mare2DEMRunner-1] Key, K. (2016). MARE2DEM: A 2-D inversion code for
   controlled-source electromagnetic and magnetotelluric data.
   *Geophysical Journal International*, 207(1), 571-588.
   doi:10.1093/gji/ggw290.
"""
