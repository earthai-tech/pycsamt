# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamRunner — invoke the Occam2D Fortran binary from Python.

The runner discovers the compiled binary in this order:

1. ``binary_path`` constructor argument (explicit override).
2. ``workdir`` — a binary named ``Occam2D`` / ``Occam2D.exe`` in the run
   directory.
3. ``PATH`` environment variable (``shutil.which``).
4. ``_source/`` sibling directory — auto-compiles with ``make`` if
   ``gfortran`` is available (see ``compile()``).

Usage
-----
>>> runner = OccamRunner(workdir="occam_run/")
>>> runner.run(max_iter=100, target_misfit=1.0)

>>> # or async — returns immediately, poll runner.is_running
>>> runner.run_async()
>>> while runner.is_running:
...     time.sleep(5)
>>> print(runner.exit_code)
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Union

from .base import OccamBase

PathLike = Union[str, Path]

__all__ = ["OccamRunner"]

# Path to the bundled Fortran source relative to this file
_SOURCE_DIR = Path(__file__).parent / "_source"
_BINARY_NAME = "Occam2D" if sys.platform != "win32" else "Occam2D.exe"


class OccamRunner(OccamBase):
    """Run the Occam2D Fortran executable from Python.

    ``OccamRunner`` is the execution layer of the Occam2D
    workflow. It assumes that an input directory already
    contains the files written by :class:`InputBuilder`:
    ``OccamDataFile.dat``, ``Occam2DMesh``, ``Occam2DModel``,
    and ``Startup``. The runner resolves a compiled
    executable, can compile the bundled Fortran source,
    launches the solver in ``workdir``, and captures standard
    output and error streams.

    Binary discovery follows a deterministic order:

    1. explicit ``binary_path`` passed to the constructor;
    2. executable named ``Occam2D`` or ``Occam2D.exe`` in
       ``workdir``;
    3. executable found on the system ``PATH``;
    4. bundled ``_source`` directory, if automatic compilation
       is enabled.

    The synchronous :meth:`run` method blocks until Occam2D
    exits. The asynchronous :meth:`run_async` method returns a
    process handle and lets the caller poll :attr:`is_running`
    or call :meth:`wait`.

    Parameters
    ----------
    workdir : path-like, default "."
        Directory containing the Occam2D run files. The binary
        is executed with this directory as its current working
        directory, so relative names inside ``Startup`` are
        resolved there. Output logs are also written there.
    binary_path : path-like, optional
        Explicit path to a compiled Occam2D executable. Use
        this when the binary is stored outside ``workdir`` or
        is not available on ``PATH``. When omitted, discovery
        the order described above.
    startup_file : str, default "Startup"
        Name of the startup file passed to the executable.
        It is resolved relative to ``workdir``.
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
        Run directory where the executable is launched.
    binary : pathlib.Path or None
        Resolved path after :meth:`discover_binary`.
    process : subprocess.Popen or None
        Live background process created by :meth:`run_async`.
    exit_code : int or None
        Return code from the most recent completed run.
    stdout_log : pathlib.Path
        File where process standard output is captured.
    stderr_log : pathlib.Path
        File where process standard error is captured.

    Notes
    -----
    :meth:`run` and :meth:`run_async` do not build input
    files. Use :class:`InputBuilder` first when starting from
    EDI data. The optional ``max_iter`` and
    ``target_misfit`` arguments to :meth:`run` patch the
    startup file in place before launch.

    See Also
    --------
    InputBuilder
        Builds the data, mesh, model, and startup files.
    OccamStartup
        Represents startup and iteration parameter vectors.
    InversionResult
        Loads results produced by a completed run.

    Examples
    --------
    Run a prepared inversion directory synchronously:

    >>> from pycsamt.models.occam2d import OccamRunner
    >>> runner = OccamRunner(workdir="occam_run")
    >>> code = runner.run(max_iter=80, target_misfit=1.0)

    Use an explicit executable path:

    >>> runner = OccamRunner(
    ...     workdir="occam_run",
    ...     binary_path="/usr/local/bin/Occam2D",
    ... )
    >>> runner.discover_binary(auto_compile=False)

    Start a background run and wait for completion:

    >>> runner = OccamRunner(workdir="occam_run")
    >>> process = runner.run_async()
    >>> runner.wait()

    References
    ----------
    .. [1] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    .. [2] Constable, S. C., Parker, R. L., and Constable,
       C. G., "Occam's inversion: A practical algorithm for
       generating smooth models from electromagnetic sounding
       data", Geophysics, 52(3), 289-300, 1987.
    """

    def __init__(
        self,
        workdir: PathLike = ".",
        binary_path: PathLike | None = None,
        startup_file: str = "Startup",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.workdir = Path(workdir)
        self._binary_path = Path(binary_path) if binary_path else None
        self.startup_file = startup_file

        self.binary: Path | None = None
        self.process: subprocess.Popen | None = None
        self.exit_code: int | None = None
        self.stdout_log = self.workdir / "occam_stdout.log"
        self.stderr_log = self.workdir / "occam_stderr.log"

    # ------------------------------------------------------------------
    # Binary discovery
    # ------------------------------------------------------------------
    def discover_binary(self, auto_compile: bool = True) -> Path:
        """Locate or compile the Occam2D executable.

        The method resolves the executable path and stores it
        on :attr:`binary`. It first honors the explicit
        ``binary_path`` constructor argument, then checks
        ``workdir``, then the system ``PATH``. If those fail
        and ``auto_compile`` is ``True``, it calls
        :meth:`compile` for the bundled Fortran source.

        Parameters
        ----------
        auto_compile : bool, default True
            If ``True``, attempt to compile the bundled source
            when no executable is found. Compilation requires
            ``make`` and a Fortran compiler such as
            ``gfortran``.

        Returns
        -------
        pathlib.Path
            Resolved path to the executable.

        Raises
        ------
        FileNotFoundError
            Raised when no executable is found and compilation
            is disabled or does not produce a binary.
        RuntimeError
            Propagated from :meth:`compile` when the compiler
            is missing or ``make`` fails.

        See Also
        --------
        OccamRunner.compile
            Compiles the bundled Fortran source.
        OccamRunner.run
            Calls this method before launching the solver.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamRunner
        >>> runner = OccamRunner("occam_run")
        >>> binary = runner.discover_binary(False)
        """
        # 1. Explicit override
        if self._binary_path and self._binary_path.is_file():
            self.binary = self._binary_path
            return self.binary

        # 2. Workdir
        candidate = self.workdir / _BINARY_NAME
        if candidate.is_file():
            self.binary = candidate
            return self.binary

        # 3. PATH
        found = shutil.which(_BINARY_NAME)
        if found:
            self.binary = Path(found)
            return self.binary

        # 4. Auto-compile from _source/
        if auto_compile:
            self.logger.info(
                "Occam2D binary not found — attempting auto-compile from %s",
                _SOURCE_DIR,
            )
            compiled = self.compile()
            if compiled.is_file():
                self.binary = compiled
                return self.binary

        raise FileNotFoundError(
            f"Occam2D binary not found.  "
            f"Compile it manually:\n"
            f"  cd {_SOURCE_DIR}\n"
            f"  make FC90=gfortran\n"
            f"then copy the produced 'Occam2D' binary to your PATH or workdir."
        )

    # ------------------------------------------------------------------
    # Compilation
    # ------------------------------------------------------------------
    def compile(self, fc: str = "gfortran", flags: str = "-O2") -> Path:
        """Compile the bundled Occam2D Fortran source.

        Compilation is performed in the package ``_source``
        directory by invoking ``make`` with ``FC90`` and
        ``FCFLAGS`` variables. The resulting executable is
        expected to be named ``Occam2D`` there.
        This method does not copy the binary into ``workdir``;
        :meth:`discover_binary` uses that path directly.

        Parameters
        ----------
        fc : str, default "gfortran"
            Fortran compiler command passed to ``make`` as
            ``FC90``. Use this to select another compiler that
            understands the bundled source.
        flags : str, default "-O2"
            Compiler flags passed to ``make`` as ``FCFLAGS``.
            Optimization flags are usually sufficient; debug
            builds can pass flags such as ``"-g"``.

        Returns
        -------
        pathlib.Path
            Path to the compiled binary inside ``_source``.

        Raises
        ------
        FileNotFoundError
            Raised when the source directory is absent.
        RuntimeError
            Raised when the requested compiler is unavailable,
            ``make`` fails, or no executable is produced.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamRunner
        >>> runner = OccamRunner("occam_run")
        >>> binary = runner.compile("gfortran", "-O2")
        """
        if not _SOURCE_DIR.is_dir():
            raise FileNotFoundError(
                f"Fortran source directory not found: {_SOURCE_DIR}"
            )

        if shutil.which(fc) is None:
            raise RuntimeError(
                f"Fortran compiler '{fc}' not found on PATH.  "
                f"Install gfortran (e.g. 'sudo apt install gfortran') "
                f"or pass fc='<compiler>'."
            )

        self.logger.info("Compiling Occam2D with %s %s in %s", fc, flags, _SOURCE_DIR)
        result = subprocess.run(
            ["make", f"FC90={fc}", f"FCFLAGS={flags}"],
            cwd=_SOURCE_DIR,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Compilation failed:\n{result.stderr}")
        compiled = _SOURCE_DIR / "Occam2D"
        if not compiled.is_file():
            raise RuntimeError("make succeeded but 'Occam2D' binary was not produced.")
        return compiled

    # ------------------------------------------------------------------
    # Run (blocking)
    # ------------------------------------------------------------------
    def run(
        self,
        max_iter: int | None = None,
        target_misfit: float | None = None,
        auto_compile: bool = True,
    ) -> int:
        """Run Occam2D synchronously.

        This method blocks until the executable exits. It
        resolves the binary, optionally patches the startup
        file, launches ``Occam2D <startup_file>`` inside
        ``workdir``, and writes process streams to
        ``occam_stdout.log`` and ``occam_stderr.log``.

        Parameters
        ----------
        max_iter : int, optional
            Temporary override for the ``Iterations to run``
            field in the startup file. The override is written
            in place before launch, so the file must be
            writable.
        target_misfit : float, optional
            Temporary override for the ``Target Misfit`` field
            in the startup file. This changes the run-control
            file before launch.
        auto_compile : bool, default True
            Passed to :meth:`discover_binary`. If ``True``,
            missing binaries may trigger compilation.

        Returns
        -------
        int
            Process exit code. A value of ``0`` indicates that
            the executable returned successfully.

        Raises
        ------
        FileNotFoundError
            Raised when the binary or patched startup file
            cannot be found.
        RuntimeError
            Propagated from compilation if automatic
            compilation fails.

        See Also
        --------
        OccamRunner.run_async
            Starts the same executable without blocking.
        OccamRunner._patch_startup
            Applies ``max_iter`` and ``target_misfit``.
        InversionResult
            Loads output files after a successful run.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamRunner
        >>> runner = OccamRunner(workdir="occam_run")
        >>> code = runner.run(max_iter=100, target_misfit=1.0)
        """
        self.discover_binary(auto_compile=auto_compile)

        # Optionally patch startup file before running
        if max_iter is not None or target_misfit is not None:
            self._patch_startup(max_iter=max_iter, target_misfit=target_misfit)

        self.logger.info("Running Occam2D in %s", self.workdir)

        with (
            open(self.stdout_log, "w") as fout,
            open(self.stderr_log, "w") as ferr,
        ):
            proc = subprocess.run(
                [str(self.binary), self.startup_file],
                cwd=self.workdir,
                stdout=fout,
                stderr=ferr,
            )

        self.exit_code = proc.returncode
        if self.exit_code != 0:
            self.logger.warning(
                "Occam2D exited with code %d.  See %s",
                self.exit_code,
                self.stderr_log,
            )
        else:
            self.logger.info("Occam2D finished successfully.")

        return self.exit_code

    # ------------------------------------------------------------------
    # Run (async)
    # ------------------------------------------------------------------
    def run_async(
        self,
        auto_compile: bool = True,
    ) -> subprocess.Popen:
        """Start Occam2D in a background process.

        The method resolves the binary and launches the solver
        with :class:`subprocess.Popen`. It returns immediately
        with a process handle. Standard output and standard
        error are redirected to ``stdout_log`` and
        ``stderr_log``.

        Parameters
        ----------
        auto_compile : bool, default True
            Passed to :meth:`discover_binary`. If ``True``,
            missing binaries may trigger compilation.

        Returns
        -------
        subprocess.Popen
            Live process handle for the background run.

        Raises
        ------
        FileNotFoundError
            Raised when no executable can be found.
        RuntimeError
            Propagated from automatic compilation failures.

        See Also
        --------
        OccamRunner.wait
            Blocks until the background process completes.
        OccamRunner.is_running
            Reports whether the process is still active.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamRunner
        >>> runner = OccamRunner(workdir="occam_run")
        >>> process = runner.run_async()
        >>> runner.is_running
        """
        self.discover_binary(auto_compile=auto_compile)

        fout = open(self.stdout_log, "w")  # kept open by Popen
        ferr = open(self.stderr_log, "w")
        self.process = subprocess.Popen(
            [str(self.binary), self.startup_file],
            cwd=self.workdir,
            stdout=fout,
            stderr=ferr,
        )
        self.logger.info(
            "Occam2D started (pid=%d).  Logging to %s",
            self.process.pid,
            self.stdout_log,
        )
        return self.process

    def wait(self) -> int:
        """Block until the asynchronous run finishes.

        Returns
        -------
        int
            Exit code returned by the background process.

        Raises
        ------
        RuntimeError
            Raised when no process has been started with
            :meth:`run_async`.
        """
        if self.process is None:
            raise RuntimeError("No async run in progress.  Call run_async() first.")
        self.exit_code = self.process.wait()
        return self.exit_code

    @property
    def is_running(self) -> bool:
        """Whether an asynchronous Occam2D process is active."""
        return self.process is not None and self.process.poll() is None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _patch_startup(
        self,
        max_iter: int | None,
        target_misfit: float | None,
    ) -> None:
        """Patch startup iteration and target-misfit fields."""
        startup_path = self.workdir / self.startup_file
        if not startup_path.is_file():
            raise FileNotFoundError(f"Startup file not found: {startup_path}")
        lines = startup_path.read_text().splitlines()
        new_lines = []
        for line in lines:
            key = line.split(":")[0].strip().upper()
            if max_iter is not None and key == "ITERATIONS TO RUN":
                new_lines.append(f"Iterations to run:  {max_iter}")
            elif target_misfit is not None and key == "TARGET MISFIT":
                new_lines.append(f"Target Misfit:      {target_misfit:.4f}")
            else:
                new_lines.append(line)
        startup_path.write_text("\n".join(new_lines) + "\n")
