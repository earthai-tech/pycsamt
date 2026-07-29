# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Adapter foundation for trusted external Maxwell solver executables.

Some verified EM forward/inversion codes worth wrapping (for example ModEM,
Occam2D, or MARE2DEM; see :mod:`pycsamt.models.modem`,
:mod:`pycsamt.models.occam2d`, :mod:`pycsamt.models.mare2dem`) are external
executables driven by input files and read back from output files, not
in-process Python callables. :class:`CallableMaxwellAdapter` in
:mod:`pycsamt.forward.maxwell.adapters` cannot wrap that shape of solver
directly. This module provides the shared, solver-independent mechanics for
that integration:

* resolving a configured executable from ``PATH`` or declared search
  directories (:func:`resolve_executable`);
* a reusable, zero-argument availability probe for backend registration
  (:func:`make_availability_probe`);
* a best-effort external solver version probe for provenance
  (:func:`probe_executable_version`);
* :class:`BaseExternalMaxwellAdapter`, which owns working-directory
  lifecycle, subprocess execution with timeout and retry, and diagnostic
  capture, so a concrete adapter implements only three solver-specific
  steps: writing input files, building the command line, and parsing
  output files back into a canonical
  :class:`~pycsamt.forward.maxwell.contracts.ForwardResult`.

No external solver is imported or executed by importing this module.
"""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
import tempfile
import time
from abc import abstractmethod
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Callable

from .adapters import (
    BackendExecutionError,
    BaseMaxwellAdapter,
    MaxwellAdapterError,
)
from .backends import BackendCapabilities
from .contracts import ForwardResult, MaxwellProblem

if TYPE_CHECKING:
    from .adapters import AdapterPolicy

__all__ = [
    "ExecutableNotFoundError",
    "ExternalProcessError",
    "ExternalRunPolicy",
    "ExternalRunResult",
    "BaseExternalMaxwellAdapter",
    "resolve_executable",
    "probe_executable_version",
    "make_availability_probe",
]


class ExecutableNotFoundError(MaxwellAdapterError):
    """Indicate that a configured external solver executable is missing.

    Examples
    --------
    >>> isinstance(ExecutableNotFoundError("missing"), MaxwellAdapterError)
    True
    """


class ExternalProcessError(BackendExecutionError):
    """Indicate that every attempt to run an external solver process failed.

    Parameters
    ----------
    message : str
        Human-readable summary, normally including the last attempt's exit
        code and a tail of its captured stderr.
    attempts : tuple of ExternalRunResult
        Every attempt made, in order, including the failing ones.

    Examples
    --------
    >>> error = ExternalProcessError("occam2d failed after 1 attempt", ())
    >>> error.attempts
    ()
    """

    def __init__(
        self, message: str, attempts: tuple[ExternalRunResult, ...]
    ) -> None:
        super().__init__(message)
        self.attempts = tuple(attempts)


def _positive_finite(
    value: float, name: str, *, allow_none: bool = False
) -> float | None:
    if value is None and allow_none:
        return None
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be finite and positive.")
    return result


def _non_negative_finite(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be finite and non-negative.")
    return result


@dataclass(frozen=True)
class ExternalRunPolicy:
    """Configure how an external solver executable is located and run.

    Parameters
    ----------
    executable : str
        Executable name resolved via ``PATH`` and ``search_paths``, or an
        absolute/relative path to the external solver binary.
    search_paths : sequence of str, optional
        Additional directories checked, in order, after ``PATH`` and
        before giving up.
    timeout_s : float or None, default=None
        Maximum wall-clock time allowed per attempt. ``None`` disables the
        per-attempt timeout.
    max_attempts : int, default=1
        Total attempts per solve, including the first. Values above one
        retry a failed or timed-out run.
    retry_backoff_s : float, default=1.0
        Base delay before each retry; attempt *n* (n > 1) waits
        ``retry_backoff_s * n`` seconds before relaunching.
    workdir : str or pathlib.Path or None, optional
        Fixed working directory reused across solves and owned by the
        caller; it is created if missing and never deleted by this
        adapter. ``None`` creates and deletes a private temporary
        directory per solve.
    keep_workdir_on_failure : bool, default=True
        Preserve a private temporary working directory (see ``workdir``)
        when every attempt fails, so its contents can be inspected. Has no
        effect when ``workdir`` is caller-supplied, since that directory
        is always preserved.
    extra_env : mapping of str to str, optional
        Extra environment variables merged over the current process
        environment for the subprocess only.
    capture_output : bool, default=True
        Capture stdout/stderr for diagnostics instead of inheriting the
        parent process streams.

    Examples
    --------
    >>> policy = ExternalRunPolicy("occam2d", max_attempts=2, timeout_s=60.0)
    >>> policy.max_attempts, policy.timeout_s
    (2, 60.0)
    """

    executable: str
    search_paths: tuple[str, ...] = ()
    timeout_s: float | None = None
    max_attempts: int = 1
    retry_backoff_s: float = 1.0
    workdir: str | None = None
    keep_workdir_on_failure: bool = True
    extra_env: Mapping[str, str] = field(default_factory=dict)
    capture_output: bool = True

    def __post_init__(self) -> None:
        executable = str(self.executable).strip()
        if not executable:
            raise ValueError("executable cannot be empty.")
        search_paths = tuple(str(value) for value in self.search_paths)
        timeout = _positive_finite(
            self.timeout_s, "timeout_s", allow_none=True
        )
        if (
            not isinstance(self.max_attempts, int)
            or isinstance(self.max_attempts, bool)
            or self.max_attempts < 1
        ):
            raise ValueError("max_attempts must be a positive integer.")
        backoff = _non_negative_finite(self.retry_backoff_s, "retry_backoff_s")
        workdir = None if self.workdir is None else str(self.workdir)
        extra_env = {
            str(key): str(value) for key, value in dict(self.extra_env).items()
        }
        object.__setattr__(self, "executable", executable)
        object.__setattr__(self, "search_paths", search_paths)
        object.__setattr__(self, "timeout_s", timeout)
        object.__setattr__(self, "retry_backoff_s", backoff)
        object.__setattr__(self, "workdir", workdir)
        object.__setattr__(self, "extra_env", MappingProxyType(extra_env))


@dataclass(frozen=True)
class ExternalRunResult:
    """Record one external-process execution attempt.

    Parameters
    ----------
    command : sequence of str
        Exact argv executed.
    returncode : int
        Process exit status; ``-1`` denotes a timeout.
    stdout, stderr : str
        Captured output; empty when ``capture_output`` was ``False`` or the
        process timed out before producing output.
    runtime_s : float
        Wall-clock duration of this attempt, in seconds.
    attempt : int
        1-based attempt number.
    workdir : str or pathlib.Path
        Directory the process was launched from.

    Examples
    --------
    >>> result = ExternalRunResult(("occam2d",), 0, "done", "", 0.5, 1, ".")
    >>> result.success
    True
    """

    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    runtime_s: float
    attempt: int
    workdir: str

    def __post_init__(self) -> None:
        command = tuple(str(value) for value in self.command)
        if not command:
            raise ValueError("command cannot be empty.")
        returncode = int(self.returncode)
        runtime = float(self.runtime_s)
        if not math.isfinite(runtime) or runtime < 0:
            raise ValueError("runtime_s must be finite and non-negative.")
        attempt = int(self.attempt)
        if attempt < 1:
            raise ValueError("attempt must be a positive integer.")
        object.__setattr__(self, "command", command)
        object.__setattr__(self, "returncode", returncode)
        object.__setattr__(self, "stdout", str(self.stdout))
        object.__setattr__(self, "stderr", str(self.stderr))
        object.__setattr__(self, "runtime_s", runtime)
        object.__setattr__(self, "attempt", attempt)
        object.__setattr__(self, "workdir", str(self.workdir))

    @property
    def success(self) -> bool:
        """Return whether the process exited with status zero.

        Returns
        -------
        bool
            ``True`` only for a normal, non-timed-out, zero exit status.

        Examples
        --------
        >>> ExternalRunResult(("a",), 1, "", "", 0.0, 1, ".").success
        False
        """
        return self.returncode == 0

    def tail(self, *, stream: str = "stderr", lines: int = 20) -> str:
        """Return the last lines of captured stdout or stderr.

        Parameters
        ----------
        stream : {"stderr", "stdout"}, default="stderr"
            Which captured stream to summarize.
        lines : int, default=20
            Maximum number of trailing lines returned.

        Returns
        -------
        str
            Newline-joined tail, or an empty string when nothing was
            captured.

        Raises
        ------
        ValueError
            If ``stream`` is not ``"stderr"`` or ``"stdout"``.

        Examples
        --------
        >>> result = ExternalRunResult(
        ...     ("a",), 1, "", "line1\\nline2\\nline3", 0.0, 1, "."
        ... )
        >>> result.tail(lines=2)
        'line2\\nline3'
        """
        if stream not in ("stderr", "stdout"):
            raise ValueError('stream must be "stderr" or "stdout".')
        text = self.stderr if stream == "stderr" else self.stdout
        content = [line for line in text.splitlines() if line.strip()]
        return "\n".join(content[-max(int(lines), 0) :])

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible representation of this attempt.

        Returns
        -------
        dict
            Command, exit status, captured output, timing, and attempt
            number. Useful for batch-run failure manifests.

        Examples
        --------
        >>> ExternalRunResult(("a",), 0, "ok", "", 1.0, 1, ".").to_dict()[
        ...     "success"
        ... ]
        True
        """
        return {
            "command": list(self.command),
            "returncode": self.returncode,
            "success": self.success,
            "stdout": self.stdout,
            "stderr": self.stderr,
            "runtime_s": self.runtime_s,
            "attempt": self.attempt,
            "workdir": self.workdir,
        }


def resolve_executable(
    name_or_path: str,
    *,
    search_paths: Sequence[str] = (),
) -> Path:
    """Resolve an external solver executable to a concrete file path.

    Parameters
    ----------
    name_or_path : str
        Executable name looked up on ``PATH``, or an absolute/relative path
        checked directly.
    search_paths : sequence of str, optional
        Extra directories checked, in order, after ``PATH`` and before
        giving up.

    Returns
    -------
    pathlib.Path
        Resolved, existing executable path.

    Raises
    ------
    ExecutableNotFoundError
        If the executable cannot be found as a direct path, on ``PATH``,
        or in any of ``search_paths``.

    Examples
    --------
    >>> resolve_executable("nope")  # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    pycsamt.forward.maxwell.external.ExecutableNotFoundError: ...
    """
    direct = Path(name_or_path)
    if direct.is_file():
        return direct
    found = shutil.which(str(name_or_path))
    if found:
        return Path(found)
    for directory in search_paths:
        candidate = Path(directory) / str(name_or_path)
        if candidate.is_file():
            return candidate
    raise ExecutableNotFoundError(
        f"executable {str(name_or_path)!r} was not found on PATH or in "
        f"{tuple(str(value) for value in search_paths)}."
    )


def probe_executable_version(
    executable: str,
    *,
    version_args: Sequence[str] = ("--version",),
    version_pattern: str | None = None,
    timeout_s: float = 5.0,
) -> str | None:
    """Best-effort external solver version string for provenance.

    Parameters
    ----------
    executable : str
        Resolved or resolvable executable path or name.
    version_args : sequence of str, default=("--version",)
        Arguments appended to the executable when probing.
    version_pattern : str, optional
        Regular expression whose first capture group extracts the version
        from the combined stdout/stderr. When omitted, the first non-empty
        output line is returned verbatim.
    timeout_s : float, default=5.0
        Maximum time allowed for the probe process.

    Returns
    -------
    str or None
        Detected version string, or ``None`` when the probe process
        cannot be started, times out, or produces no matching output.
        This function never raises for an ordinary probe failure; a
        concrete adapter must still supply a non-empty
        ``backend_version`` to
        :class:`~pycsamt.forward.maxwell.contracts.ForwardResult`
        when this returns ``None``.

    Examples
    --------
    >>> probe_executable_version("does-not-exist-xyz") is None
    True
    """
    try:
        completed = subprocess.run(
            [str(executable), *version_args],
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    combined = f"{completed.stdout}\n{completed.stderr}"
    if version_pattern is not None:
        match = re.search(version_pattern, combined)
        return match.group(1) if match else None
    for line in combined.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return None


def make_availability_probe(
    name_or_path: str,
    *,
    search_paths: Sequence[str] = (),
) -> Callable[[], tuple[bool, str | None]]:
    """Build a zero-argument availability probe for backend registration.

    Parameters
    ----------
    name_or_path : str
        Executable resolved the same way as :func:`resolve_executable`.
    search_paths : sequence of str, optional
        Extra directories checked before giving up.

    Returns
    -------
    callable
        Zero-argument function returning ``(available, reason)``, directly
        usable as
        :attr:`~pycsamt.forward.maxwell.backends.BackendRegistration.availability_probe`.

    Examples
    --------
    >>> probe = make_availability_probe("does-not-exist-xyz")
    >>> probe()[0]
    False
    """

    def _probe() -> tuple[bool, str | None]:
        try:
            resolve_executable(name_or_path, search_paths=search_paths)
        except ExecutableNotFoundError as exc:
            return False, str(exc)
        return True, None

    return _probe


class BaseExternalMaxwellAdapter(BaseMaxwellAdapter):
    """Base class for adapters that run a trusted external solver process.

    Concrete adapters wrap file-based external tools (for example ModEM,
    Occam2D, or MARE2DEM) launched as subprocesses rather than called as an
    in-process Python function. This class owns the shared,
    solver-independent mechanics: working-directory lifecycle, executable
    resolution, subprocess execution with timeout and retry, and captured
    diagnostics. A concrete subclass implements only three solver-specific
    extension points:

    ``_prepare_inputs(problem, workdir)``
        Write the external tool's input files for ``problem`` into
        ``workdir`` and return a context object (any value) carried
        through to the other two extension points.

    ``_build_command(problem, workdir, executable, context)``
        Return the argv sequence that runs the external tool against the
        files written by ``_prepare_inputs``.

    ``_parse_result(problem, workdir, run_result, context)``
        Read the external tool's output files from ``workdir`` and return
        a canonical :class:`~pycsamt.forward.maxwell.contracts.ForwardResult`.
        ``run_result`` is the successful :class:`ExternalRunResult`.

    Parameters
    ----------
    capabilities : BackendCapabilities
        Immutable declaration for the exact external solver version.
    run_policy : ExternalRunPolicy
        Executable resolution, timeout, retry, and working-directory
        rules.
    policy : AdapterPolicy or None, optional
        Solver-independent result acceptance policy (convergence,
        residual, and validity checks applied after ``_parse_result``
        returns).

    Examples
    --------
    A minimal concrete adapter (see the module docstring for context)
    would look like::

        class DemoExternalAdapter(BaseExternalMaxwellAdapter):
            def _prepare_inputs(self, problem, workdir):
                (workdir / "input.txt").write_text(str(problem.problem_hash))
                return None

            def _build_command(self, problem, workdir, executable, context):
                return [str(executable), "input.txt", "output.txt"]

            def _parse_result(
                self, problem, workdir, run_result, context
            ): ...  # read workdir / "output.txt" and build a ForwardResult
    """

    def __init__(
        self,
        capabilities: BackendCapabilities,
        run_policy: ExternalRunPolicy,
        policy: AdapterPolicy | None = None,
    ) -> None:
        super().__init__(capabilities, policy)
        if not isinstance(run_policy, ExternalRunPolicy):
            raise TypeError("run_policy must be an ExternalRunPolicy.")
        self._run_policy = run_policy

    @property
    def run_policy(self) -> ExternalRunPolicy:
        """Return the immutable external-process execution policy.

        Returns
        -------
        ExternalRunPolicy
            Executable resolution, timeout, retry, and
            working-directory rules used by every
            :meth:`~pycsamt.forward.maxwell.adapters.BaseMaxwellAdapter.solve`
            call.

        Examples
        --------
        See :class:`BaseExternalMaxwellAdapter` for a complete subclass
        example.
        """
        return self._run_policy

    def resolve_executable(self) -> Path:
        """Resolve this adapter's configured executable to a concrete path.

        Returns
        -------
        pathlib.Path
            Resolved executable path.

        Raises
        ------
        ExecutableNotFoundError
            If the executable cannot be found on ``PATH`` or in
            :attr:`ExternalRunPolicy.search_paths`.

        Examples
        --------
        See :class:`BaseExternalMaxwellAdapter` for a complete subclass
        example.
        """
        return resolve_executable(
            self._run_policy.executable,
            search_paths=self._run_policy.search_paths,
        )

    def _solve_backend(self, problem: MaxwellProblem) -> ForwardResult:
        executable = self.resolve_executable()
        workdir, owned = self._acquire_workdir()
        succeeded = False
        try:
            context = self._prepare_inputs(problem, workdir)
            command = self._build_command(
                problem, workdir, executable, context
            )
            run_result = self._run_with_retries(tuple(command), workdir)
            result = self._parse_result(problem, workdir, run_result, context)
            succeeded = True
            return result
        finally:
            self._release_workdir(workdir, owned, keep=not succeeded)

    def _acquire_workdir(self) -> tuple[Path, bool]:
        fixed = self._run_policy.workdir
        if fixed is not None:
            path = Path(fixed)
            path.mkdir(parents=True, exist_ok=True)
            return path, False
        return Path(tempfile.mkdtemp(prefix="pycsamt-maxwell-external-")), True

    def _release_workdir(
        self, workdir: Path, owned: bool, *, keep: bool
    ) -> None:
        if not owned:
            return
        if keep and self._run_policy.keep_workdir_on_failure:
            return
        shutil.rmtree(workdir, ignore_errors=True)

    def _run_with_retries(
        self, command: tuple[str, ...], workdir: Path
    ) -> ExternalRunResult:
        attempts: list[ExternalRunResult] = []
        for attempt in range(1, self._run_policy.max_attempts + 1):
            if attempt > 1:
                time.sleep(self._run_policy.retry_backoff_s * attempt)
            result = _run_subprocess(
                command, workdir, self._run_policy, attempt
            )
            attempts.append(result)
            if result.success:
                return result
        last = attempts[-1]
        raise ExternalProcessError(
            f"external solver {self.capabilities.name!r} failed after "
            f"{len(attempts)} attempt(s); last exit code {last.returncode}: "
            f"{last.tail()}",
            tuple(attempts),
        )

    @abstractmethod
    def _prepare_inputs(self, problem: MaxwellProblem, workdir: Path) -> Any:
        raise NotImplementedError

    @abstractmethod
    def _build_command(
        self,
        problem: MaxwellProblem,
        workdir: Path,
        executable: Path,
        context: Any,
    ) -> Sequence[str]:
        raise NotImplementedError

    @abstractmethod
    def _parse_result(
        self,
        problem: MaxwellProblem,
        workdir: Path,
        run_result: ExternalRunResult,
        context: Any,
    ) -> ForwardResult:
        raise NotImplementedError


def _run_subprocess(
    command: tuple[str, ...],
    workdir: Path,
    run_policy: ExternalRunPolicy,
    attempt: int,
) -> ExternalRunResult:
    env = None
    if run_policy.extra_env:
        env = dict(os.environ)
        env.update(run_policy.extra_env)
    start = time.monotonic()
    try:
        completed = subprocess.run(
            command,
            cwd=str(workdir),
            timeout=run_policy.timeout_s,
            capture_output=run_policy.capture_output,
            text=True,
            env=env,
        )
    except subprocess.TimeoutExpired as exc:
        runtime = time.monotonic() - start
        stdout = exc.stdout if isinstance(exc.stdout, str) else ""
        stderr = exc.stderr if isinstance(exc.stderr, str) else ""
        stderr = f"{stderr}\ntimed out after {run_policy.timeout_s}s".strip()
        return ExternalRunResult(
            command, -1, stdout, stderr, runtime, attempt, workdir
        )
    except OSError as exc:
        runtime = time.monotonic() - start
        return ExternalRunResult(
            command, -1, "", str(exc), runtime, attempt, workdir
        )
    runtime = time.monotonic() - start
    return ExternalRunResult(
        command,
        completed.returncode,
        completed.stdout or "",
        completed.stderr or "",
        runtime,
        attempt,
        workdir,
    )
