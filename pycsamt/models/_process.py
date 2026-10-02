# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Run a solver executable while streaming its console output.

The Occam2D, ModEM and MARE2DEM runners historically launched their
binaries with :func:`subprocess.run` and no captured output: a caller could
neither show progress nor stop a run.  :func:`run_streamed` is the shared
replacement they use when a caller passes ``on_output`` / ``cancel``:

* output lines (stdout and stderr merged) are delivered as they appear,
  decoded leniently from **binary** pipes;
* gfortran's pipe buffering is switched off (``GFORTRAN_UNBUFFERED_ALL``),
  otherwise a Fortran solver's lines only arrive in large bursts;
* ``cancel`` is polled every 0.2 s, independent of output, so a solver that
  prints nothing for minutes can still be stopped; the whole process tree
  (``mpirun`` and its ranks) is terminated.

Without ``on_output`` and ``cancel`` the runners keep their original
behaviour.
"""

from __future__ import annotations

import os
import queue
import signal
import subprocess
import sys
import threading
import time
from collections.abc import Callable, Sequence
from pathlib import Path

__all__ = ["ProcessCancelled", "run_streamed", "solver_env"]

_IS_WIN = sys.platform == "win32"


class ProcessCancelled(InterruptedError):
    """Raised by :func:`run_streamed` when ``cancel()`` returned True."""


def solver_env(env: dict | None = None) -> dict:
    """Return a copy of *env* (default: ``os.environ``) for a solver run.

    Unbuffered gfortran output makes iteration lines arrive when they are
    printed instead of when a 8 KiB pipe buffer fills.
    """
    out = dict(os.environ if env is None else env)
    out.setdefault("GFORTRAN_UNBUFFERED_ALL", "y")
    return out


def _kill_tree(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    try:
        if _IS_WIN:
            subprocess.run(
                ["taskkill", "/T", "/F", "/PID", str(proc.pid)],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                creationflags=subprocess.CREATE_NO_WINDOW,
            )
        else:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
    except Exception:
        proc.kill()
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        pass


def run_streamed(
    argv: Sequence[str],
    *,
    cwd: str | Path | None = None,
    env: dict | None = None,
    on_output: Callable[[str], None] | None = None,
    cancel: Callable[[], bool] | None = None,
    timeout: float | None = None,
    tee: str | Path | None = None,
) -> int:
    """Run *argv*, forwarding each output line to *on_output*.

    Parameters
    ----------
    argv : sequence of str
        Command and arguments.
    cwd : path-like, optional
        Working directory.
    env : dict, optional
        Environment; :func:`solver_env` is applied to it.
    on_output : callable, optional
        Called with every output line (trailing newline removed).
    cancel : callable returning bool, optional
        Polled every 0.2 s; when it returns True the process tree is
        terminated and :class:`ProcessCancelled` is raised.
    timeout : float, optional
        Wall-clock limit in seconds; the process tree is killed and
        :class:`subprocess.TimeoutExpired` is raised on expiry.
    tee : path-like, optional
        Also append every line to this file (e.g. the runner's stdout log).

    Returns
    -------
    int
        The process exit code.
    """
    kwargs: dict = {}
    if _IS_WIN:
        kwargs["creationflags"] = (subprocess.CREATE_NO_WINDOW
                                   | subprocess.CREATE_NEW_PROCESS_GROUP)
    else:
        kwargs["start_new_session"] = True
    proc = subprocess.Popen(
        [str(a) for a in argv],
        cwd=None if cwd is None else str(cwd),
        env=solver_env(env),
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        **kwargs,
    )
    lines: queue.Queue = queue.Queue()

    def _reader() -> None:
        for raw in iter(proc.stdout.readline, b""):
            lines.put(raw.decode("utf-8", "replace").rstrip("\r\n"))
        proc.stdout.close()
        lines.put(None)

    threading.Thread(target=_reader, daemon=True).start()
    sink = open(tee, "a", encoding="utf-8") if tee else None
    start = time.monotonic()
    try:
        while True:
            try:
                line = lines.get(timeout=0.2)
            except queue.Empty:
                line = ""
            else:
                if line is None:
                    break
                if sink is not None:
                    sink.write(line + "\n")
                    sink.flush()
                if on_output is not None:
                    on_output(line)
            if cancel is not None and cancel():
                _kill_tree(proc)
                raise ProcessCancelled("Run cancelled by the user.")
            if timeout is not None and time.monotonic() - start > timeout:
                _kill_tree(proc)
                raise subprocess.TimeoutExpired(list(argv), timeout)
        return proc.wait()
    finally:
        if sink is not None:
            sink.close()
        if proc.poll() is None:
            _kill_tree(proc)
