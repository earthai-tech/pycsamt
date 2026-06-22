# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Subprocess launcher for MARE2DEM inversion runs."""

from __future__ import annotations

import shlex
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from .base import Mare2DEMBase
from .config import Mare2DEMConfig
from .doc import _mare2dem_param_docs as _params
from .source import SourceManager

PathLike = str | Path

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
    ) -> "InversionResult | None":
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

        binary = _resolve_binary(cfg.binary, self._source_mgr)
        if binary is None:
            raise FileNotFoundError(
                f"MARE2DEM binary '{cfg.binary}' not found on PATH or in "
                "the source directory. "
                "Build it first:  SourceManager().download(); "
                "SourceManager().build()"
            )

        stem = Path(str(resistivity_stem)).with_suffix("").name

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

        proc = subprocess.run(
            cmd,
            cwd=str(self.workdir),
            timeout=timeout,
        )
        proc.check_returncode()

        if load_result:
            from .results import InversionResult
            return InversionResult(self.workdir, config=cfg)
        return None

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
        stem = Path(str(resistivity_stem)).with_suffix("").name

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
.. [1] Key, K. (2016). MARE2DEM: A 2-D inversion code for
   controlled-source electromagnetic and magnetotelluric data.
   *Geophysical Journal International*, 207(1), 571-588.
   doi:10.1093/gji/ggw290.
"""
