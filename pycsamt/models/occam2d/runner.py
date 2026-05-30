# -*- coding: utf-8 -*-
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

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional, Union

from .base import OccamBase

PathLike = Union[str, Path]

__all__ = ["OccamRunner"]

# Path to the bundled Fortran source relative to this file
_SOURCE_DIR = Path(__file__).parent / "_source"
_BINARY_NAME = "Occam2D" if sys.platform != "win32" else "Occam2D.exe"


class OccamRunner(OccamBase):
    """Subprocess wrapper for the Occam2D Fortran binary.

    Parameters
    ----------
    workdir : path-like
        Directory that contains the four Occam input files (data, mesh,
        model, startup).  The binary is executed here.
    binary_path : path-like, optional
        Explicit path to the compiled ``Occam2D`` executable.  When not
        supplied the runner uses the discovery order described above.
    startup_file : str
        Name of the startup file inside *workdir*.  Default ``"Startup"``.

    Attributes
    ----------
    binary : Path | None
        Resolved path to the binary (set after ``discover_binary()``).
    process : subprocess.Popen | None
        Live process handle (only set during ``run_async()``).
    exit_code : int | None
        Return code of the last completed run.
    stdout_log : Path
        Path where stdout is captured (``occam_stdout.log``).
    stderr_log : Path
        Path where stderr is captured (``occam_stderr.log``).
    """

    def __init__(
        self,
        workdir: PathLike = ".",
        binary_path: Optional[PathLike] = None,
        startup_file: str = "Startup",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.workdir      = Path(workdir)
        self._binary_path = Path(binary_path) if binary_path else None
        self.startup_file = startup_file

        self.binary:     Optional[Path]             = None
        self.process:    Optional[subprocess.Popen] = None
        self.exit_code:  Optional[int]              = None
        self.stdout_log  = self.workdir / "occam_stdout.log"
        self.stderr_log  = self.workdir / "occam_stderr.log"

    # ------------------------------------------------------------------
    # Binary discovery
    # ------------------------------------------------------------------
    def discover_binary(self, auto_compile: bool = True) -> Path:
        """Locate (or compile) the Occam2D binary.

        Sets and returns ``self.binary``.  Raises ``FileNotFoundError``
        if the binary cannot be found or compiled.

        Parameters
        ----------
        auto_compile : bool
            When ``True``, attempt ``compile()`` if the binary is not
            found anywhere else.
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
        """Compile the bundled Fortran source.

        Parameters
        ----------
        fc : str
            Fortran compiler command (default ``gfortran``).
        flags : str
            Compiler flags passed as ``FCFLAGS``.

        Returns
        -------
        Path
            Path to the compiled binary (inside ``_source/``).
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
            raise RuntimeError(
                f"Compilation failed:\n{result.stderr}"
            )
        compiled = _SOURCE_DIR / "Occam2D"
        if not compiled.is_file():
            raise RuntimeError(
                "make succeeded but 'Occam2D' binary was not produced."
            )
        return compiled

    # ------------------------------------------------------------------
    # Run (blocking)
    # ------------------------------------------------------------------
    def run(
        self,
        max_iter: Optional[int]    = None,
        target_misfit: Optional[float] = None,
        auto_compile: bool = True,
    ) -> int:
        """Run Occam2D synchronously (blocks until completion).

        Parameters
        ----------
        max_iter : int, optional
            Overrides the ``Iterations to run`` field in Startup before
            launching.  Requires the startup file to be writable.
        target_misfit : float, optional
            Overrides ``Target Misfit`` in Startup.
        auto_compile : bool
            Passed to ``discover_binary()``.

        Returns
        -------
        int
            Process exit code (0 = success).
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
                self.exit_code, self.stderr_log,
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
        """Start Occam2D in a background process and return immediately.

        Poll ``runner.is_running`` or call ``runner.wait()`` to block.
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
            self.process.pid, self.stdout_log,
        )
        return self.process

    def wait(self) -> int:
        """Block until the async run finishes.  Returns exit code."""
        if self.process is None:
            raise RuntimeError("No async run in progress.  Call run_async() first.")
        self.exit_code = self.process.wait()
        return self.exit_code

    @property
    def is_running(self) -> bool:
        """True if an async process is currently active."""
        return self.process is not None and self.process.poll() is None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _patch_startup(
        self,
        max_iter: Optional[int],
        target_misfit: Optional[float],
    ) -> None:
        """In-place edit of the Startup file to override iteration/misfit."""
        startup_path = self.workdir / self.startup_file
        if not startup_path.is_file():
            raise FileNotFoundError(
                f"Startup file not found: {startup_path}"
            )
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
