# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Validated, synchronous execution boundary for an Occam1D solver."""

from __future__ import annotations

import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Any

from ...compat.sklearn import validate_params
from .base import Occam1DBase
from .data import Occam1DData
from .model import Occam1DModel
from .schema import RUNNER_RUN_SCHEMA, RUNNER_SCHEMA
from .startup import Occam1DStartup

__all__ = [
    "Occam1DExecutionError",
    "Occam1DRunner",
    "Occam1DRunnerState",
    "Occam1DRunRecord",
]


class Occam1DRunnerState(str, Enum):
    """Lifecycle states for one runner's most recent execution."""

    IDLE = "idle"
    READY = "ready"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    TIMED_OUT = "timed_out"


@dataclass(frozen=True)
class Occam1DRunRecord:
    """Immutable facts captured for one completed execution attempt.

    Parameters
    ----------
    command : tuple of str
        Exact argument vector passed to :func:`subprocess.run`.
    workdir : pathlib.Path
        Resolved process working directory.
    exit_code : int
        Solver return code. ``-9`` represents the package timeout sentinel.
    runner_state : Occam1DRunnerState
        Terminal run state.
    elapsed_seconds : float
        Non-negative wall-clock duration measured with a monotonic clock.
    timeout_seconds : float or None
        Requested timeout, if any.
    stdout_log, stderr_log : pathlib.Path
        Captured process stream destinations.
    """

    command: tuple[str, ...]
    workdir: Path
    exit_code: int
    state: Occam1DRunnerState
    elapsed_seconds: float
    timeout_seconds: float | None
    stdout_log: Path
    stderr_log: Path

    @property
    def succeeded(self) -> bool:
        """Whether the executable returned exit code zero."""
        return self.state is Occam1DRunnerState.SUCCEEDED

    @property
    def timed_out(self) -> bool:
        """Whether Python terminated waiting after the requested timeout."""
        return self.state is Occam1DRunnerState.TIMED_OUT

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly execution mapping."""
        values = asdict(self)
        values["command"] = list(self.command)
        values["workdir"] = str(self.workdir)
        values["state"] = self.state.value
        values["stdout_log"] = str(self.stdout_log)
        values["stderr_log"] = str(self.stderr_log)
        values["succeeded"] = self.succeeded
        values["timed_out"] = self.timed_out
        return values


class Occam1DExecutionError(RuntimeError):
    """Raised when ``check=True`` and an Occam1D run does not succeed."""

    def __init__(self, record: Occam1DRunRecord, stderr_tail: str = ""):
        self.record = record
        self.stderr_tail = stderr_tail
        reason = (
            "timed out"
            if record.timed_out
            else f"exited with code {record.exit_code}"
        )
        message = (
            f"Occam1D {reason}; see stderr log {record.stderr_log}."
        )
        if stderr_tail:
            message += f" Last stderr output: {stderr_tail}"
        super().__init__(message)


class Occam1DRunner(Occam1DBase):
    """Discover, validate, and synchronously run an Occam1D executable.

    Parameters
    ----------
    workdir : path-like, default="."
        Directory containing startup, data, and layered-model inputs. It is
        resolved to an absolute path but need not exist until preflight/run.
    binary_path : path-like, optional
        Explicit executable. When supplied, failure to resolve that exact file
        raises; discovery does not silently fall back to another executable.
    startup_file : str, default="Startup"
        Startup filename or path passed as the executable's first argument.
        A relative value is resolved inside ``workdir``.
    binary_name : str, default="Occam1D"
        Name searched in ``workdir`` and then on ``PATH`` when no explicit
        executable is supplied.
    verbose, logger, path, metadata, stream
        Shared options inherited from :class:`Occam1DBase`.

    Attributes
    ----------
    state : Occam1DRunnerState
        Execution state independent of the base object's path lifecycle.
    binary : pathlib.Path or None
        Resolved executable after discovery.
    exit_code : int or None
        Most recent terminal return code.
    last_run : Occam1DRunRecord or None
        Immutable record of the most recent process attempt.
    stdout_log, stderr_log : pathlib.Path
        Process output files inside ``workdir``.

    Notes
    -----
    Logging and verbosity are independent. Process diagnostics use the logger;
    concise stage messages use the configured verbosity stream.

    ``max_iterations`` and ``target_misfit`` overrides are persisted to the
    startup file before launch. The binary is resolved first, so a missing
    executable does not modify run controls.

    The runner never invokes a shell. Extra arguments are appended to a list,
    avoiding command-string interpretation.

    Examples
    --------
    >>> runner = Occam1DRunner("occam1d-inversion")
    >>> runner.preflight()  # doctest: +SKIP
    >>> code = runner.run(timeout=3600, check=True)  # doctest: +SKIP
    >>> runner.last_run.succeeded  # doctest: +SKIP
    True
    """

    @validate_params(RUNNER_SCHEMA)
    def __init__(
        self,
        workdir=".",
        binary_path=None,
        startup_file="Startup",
        binary_name="Occam1D",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.workdir = self._normalize_path(workdir)
        self.binary_path = (
            self._normalize_path(binary_path)
            if binary_path is not None
            else None
        )
        self.startup_file = self._name("startup_file", startup_file)
        self.binary_name = self._name("binary_name", binary_name)
        self.binary: Path | None = None
        self.exit_code: int | None = None
        self.runner_state = Occam1DRunnerState.IDLE
        self.last_run: Occam1DRunRecord | None = None
        self.last_command: tuple[str, ...] = tuple()
        self.stdout_log = self.workdir / "occam1d_stdout.log"
        self.stderr_log = self.workdir / "occam1d_stderr.log"

    @staticmethod
    def _name(name: str, value: str) -> str:
        text = value.strip()
        if not text:
            raise ValueError(f"{name} cannot be empty.")
        if "\x00" in text or "\n" in text or "\r" in text:
            raise ValueError(f"{name} contains an invalid control character.")
        return text

    @property
    def startup_path(self) -> Path:
        """Resolved startup path used during preflight and execution."""
        candidate = Path(self.startup_file)
        if not candidate.is_absolute():
            candidate = self.workdir / candidate
        return candidate.resolve(strict=False)

    @property
    def is_ready(self) -> bool:
        """Whether the latest preflight completed without launching."""
        return self.runner_state is Occam1DRunnerState.READY

    @property
    def succeeded(self) -> bool:
        """Whether the most recent process attempt succeeded."""
        return self.runner_state is Occam1DRunnerState.SUCCEEDED

    def discover_binary(self) -> Path:
        """Resolve the explicit, local, or ``PATH`` Occam1D executable.

        Discovery order is deterministic: explicit ``binary_path``, platform
        names inside ``workdir``, then those names on ``PATH``.
        """
        if self.binary_path is not None:
            if not self.binary_path.exists():
                raise FileNotFoundError(
                    f"Explicit Occam1D executable not found: "
                    f"{self.binary_path}"
                )
            if not self.binary_path.is_file():
                raise IsADirectoryError(
                    f"Explicit Occam1D executable is not a file: "
                    f"{self.binary_path}"
                )
            self.binary = self.binary_path.resolve()
            return self.binary

        names = [self.binary_name]
        if sys.platform == "win32" and not self.binary_name.lower().endswith(
            ".exe"
        ):
            names.insert(0, self.binary_name + ".exe")
        for name in names:
            candidate = self.workdir / name
            if candidate.is_file():
                self.binary = candidate.resolve()
                return self.binary
        for name in names:
            found = shutil.which(name)
            if found:
                self.binary = Path(found).resolve()
                return self.binary
        raise FileNotFoundError(
            "Occam1D executable was not found. Pass binary_path, place "
            f"{self.binary_name!r} in the run directory, or add it to PATH."
        )

    def preflight(self) -> dict[str, Any]:
        """Validate the complete input set without launching the solver.

        Returns
        -------
        dict
            Startup, data, model objects and their resolved paths.

        Raises
        ------
        FileNotFoundError, IsADirectoryError
            If the working directory or a referenced input is unavailable.
        ValueError
            If a native input is corrupt or parameter/layer counts differ.
        """
        if not self.workdir.exists():
            raise FileNotFoundError(
                f"Occam1D run directory not found: {self.workdir}"
            )
        if not self.workdir.is_dir():
            raise NotADirectoryError(
                f"Occam1D workdir is not a directory: {self.workdir}"
            )
        startup = Occam1DStartup.read(self.startup_path)
        model_path = self._referenced_path(startup.model_file)
        data_path = self._referenced_path(startup.data_file)
        model = Occam1DModel.read(model_path)
        data = Occam1DData.read(data_path)
        startup.apply_to_model(model)
        self._bind_path(self.workdir)
        self.runner_state = Occam1DRunnerState.READY
        self.logger.debug(
            "Occam1D preflight passed for %s: %d data, %d layers.",
            self.workdir,
            data.n_data,
            model.n_layers,
        )
        return {
            "startup": startup,
            "startup_path": self.startup_path,
            "model": model,
            "model_path": model_path,
            "data": data,
            "data_path": data_path,
        }

    def _referenced_path(self, value: str) -> Path:
        """Resolve and require one startup-referenced regular file."""
        candidate = Path(value)
        if not candidate.is_absolute():
            candidate = self.workdir / candidate
        candidate = candidate.resolve(strict=False)
        if not candidate.exists():
            raise FileNotFoundError(
                f"Startup-referenced input file not found: {candidate}"
            )
        if not candidate.is_file():
            raise IsADirectoryError(
                f"Startup-referenced input is not a file: {candidate}"
            )
        return candidate

    def _patch_startup(self, max_iterations=None, target_misfit=None):
        """Persist validated overrides and return the startup object."""
        startup = Occam1DStartup.read(self.startup_path)
        if max_iterations is not None:
            startup.max_iterations = int(max_iterations)
        if target_misfit is not None:
            startup.target_misfit = float(target_misfit)
        startup.validate()
        startup.write(self.startup_path)
        return startup

    @staticmethod
    def _arguments(extra_args) -> tuple[str, ...]:
        """Normalize optional shell-free executable arguments."""
        if extra_args is None:
            return tuple()
        if isinstance(extra_args, (str, bytes)):
            raise TypeError("extra_args must be an iterable of arguments.")
        result = []
        for value in extra_args:
            text = str(value)
            if "\x00" in text:
                raise ValueError("extra_args cannot contain a null byte.")
            result.append(text)
        return tuple(result)

    @validate_params(RUNNER_RUN_SCHEMA)
    def run(
        self,
        *,
        max_iterations=None,
        target_misfit=None,
        timeout=None,
        check=False,
        extra_args=None,
    ) -> int:
        """Run the solver synchronously and return its exit code.

        Parameters
        ----------
        max_iterations : int, optional
            Positive persistent override for ``Iterations to run``.
        target_misfit : float, optional
            Positive persistent normalized RMS target override.
        timeout : float, optional
            Positive wall-clock limit in seconds. Timeout returns ``-9`` unless
            ``check=True``, in which case :class:`Occam1DExecutionError` is
            raised after the run record is stored.
        check : bool, default=False
            Raise :class:`Occam1DExecutionError` for any non-zero outcome.
        extra_args : iterable, optional
            Additional argument values appended without shell evaluation.

        Returns
        -------
        int
            Solver return code, or ``-9`` for a timeout.

        Raises
        ------
        Occam1DExecutionError
            If ``check=True`` and the solver fails or times out.
        OSError
            If the validated executable cannot be launched.
        """
        arguments = self._arguments(extra_args)
        self.preflight()
        self.discover_binary()
        if max_iterations is not None or target_misfit is not None:
            self._patch_startup(max_iterations, target_misfit)
            self.preflight()

        startup_argument = (
            str(self.startup_path)
            if Path(self.startup_file).is_absolute()
            else self.startup_file
        )
        command = (str(self.binary), startup_argument, *arguments)
        self.last_command = command
        self.exit_code = None
        self.last_run = None
        self.runner_state = Occam1DRunnerState.RUNNING
        self._emit(f"Running Occam1D in {self.workdir}")
        self.logger.info("Running command %s in %s", command, self.workdir)
        started = time.monotonic()
        terminal_state = Occam1DRunnerState.FAILED
        exit_code = -1
        self.workdir.mkdir(parents=True, exist_ok=True)
        with (
            self.stdout_log.open("w", encoding="utf8") as stdout,
            self.stderr_log.open("w", encoding="utf8") as stderr,
        ):
            try:
                completed = subprocess.run(
                    list(command),
                    cwd=self.workdir,
                    stdout=stdout,
                    stderr=stderr,
                    timeout=timeout,
                    check=False,
                )
                exit_code = int(completed.returncode)
                terminal_state = (
                    Occam1DRunnerState.SUCCEEDED
                    if exit_code == 0
                    else Occam1DRunnerState.FAILED
                )
            except subprocess.TimeoutExpired:
                exit_code = -9
                terminal_state = Occam1DRunnerState.TIMED_OUT
                self.logger.warning(
                    "Occam1D exceeded timeout=%s seconds in %s.",
                    timeout,
                    self.workdir,
                )
            except OSError:
                elapsed = max(0.0, time.monotonic() - started)
                self.exit_code = -1
                self.runner_state = Occam1DRunnerState.FAILED
                self.last_run = Occam1DRunRecord(
                    command=command,
                    workdir=self.workdir,
                    exit_code=-1,
                    state=Occam1DRunnerState.FAILED,
                    elapsed_seconds=elapsed,
                    timeout_seconds=(
                        float(timeout) if timeout is not None else None
                    ),
                    stdout_log=self.stdout_log,
                    stderr_log=self.stderr_log,
                )
                self.logger.exception(
                    "Could not launch Occam1D executable %s.", self.binary
                )
                raise
        elapsed = max(0.0, time.monotonic() - started)
        self.exit_code = exit_code
        self.runner_state = terminal_state
        self.last_run = Occam1DRunRecord(
            command=command,
            workdir=self.workdir,
            exit_code=exit_code,
            state=terminal_state,
            elapsed_seconds=elapsed,
            timeout_seconds=float(timeout) if timeout is not None else None,
            stdout_log=self.stdout_log,
            stderr_log=self.stderr_log,
        )
        if self.succeeded:
            self.logger.info(
                "Occam1D finished successfully in %.3f seconds.", elapsed
            )
            self._emit("Occam1D finished successfully")
        else:
            self.logger.warning(
                "Occam1D finished with state=%s and exit_code=%d; see %s.",
                self.runner_state.value,
                exit_code,
                self.stderr_log,
            )
        if check and not self.succeeded:
            raise Occam1DExecutionError(
                self.last_run, self.stderr_tail(lines=10)
            )
        return exit_code

    def stderr_tail(self, *, lines=20) -> str:
        """Return the final non-empty stderr lines from the latest log."""
        if isinstance(lines, bool) or not isinstance(lines, int):
            raise TypeError("lines must be a positive integer.")
        if lines < 1:
            raise ValueError("lines must be a positive integer.")
        if not self.stderr_log.is_file():
            return ""
        values = [
            line.strip()
            for line in self.stderr_log.read_text(
                encoding="utf8", errors="replace"
            ).splitlines()
            if line.strip()
        ]
        return "\n".join(values[-lines:])

    def summary(self) -> str:
        """Return a concise human-readable runner summary."""
        binary = str(self.binary) if self.binary is not None else "unresolved"
        return (
            f"Occam1DRunner(state={self.runner_state.value!r}, "
            f"workdir={str(self.workdir)!r}, binary={binary!r}, "
            f"exit_code={self.exit_code!r})"
        )

    def diagnostics(self) -> dict[str, Any]:
        """Extend lifecycle diagnostics with process execution state."""
        values = super().diagnostics()
        values.update(
            {
                "workdir": str(self.workdir),
                "startup_file": self.startup_file,
                "startup_path": str(self.startup_path),
                "binary_name": self.binary_name,
                "binary_path": (
                    str(self.binary_path)
                    if self.binary_path is not None
                    else None
                ),
                "binary": str(self.binary) if self.binary else None,
                "runner_state": self.runner_state.value,
                "exit_code": self.exit_code,
                "last_command": list(self.last_command),
                "stdout_log": str(self.stdout_log),
                "stderr_log": str(self.stderr_log),
                "last_run": (
                    self.last_run.to_dict()
                    if self.last_run is not None
                    else None
                ),
            }
        )
        return values
