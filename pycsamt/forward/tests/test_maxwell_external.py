"""Tests for the external-solver adapter foundation."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    BackendCapabilities,
    BaseExternalMaxwellAdapter,
    ExecutableNotFoundError,
    ExternalProcessError,
    ExternalRunPolicy,
    ExternalRunResult,
    ForwardResult,
    MaxwellMesh,
    MaxwellProblem,
    ReceiverSet,
    SolverDiagnostics,
    make_availability_probe,
    probe_executable_version,
    resolve_executable,
)


def _problem(frequencies=(10, 1)):
    mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
    return MaxwellProblem(
        mesh,
        np.ones(mesh.shape),
        frequencies,
        ReceiverSet([[0.5, 0]], ["S"]),
        ("zxy",),
    )


def _capabilities(**changes):
    values = dict(
        name="demo-external",
        version="1",
        dimensions=(2,),
        components=("zxy",),
        verified_benchmarks=("half-space",),
    )
    values.update(changes)
    return BackendCapabilities(**values)


class _ScriptedAdapter(BaseExternalMaxwellAdapter):
    """Runs the current Python interpreter with an inline script.

    ``script`` receives the workdir as ``sys.argv[1]`` and must create a
    file named ``ok`` in that directory on success (the mechanism
    ``_parse_result`` uses to build a ForwardResult).
    """

    def __init__(self, script: str, *, workdirs: list | None = None, **kwargs):
        super().__init__(_capabilities(), **kwargs)
        self._script = script
        self._workdirs = workdirs

    def _prepare_inputs(self, problem, workdir):
        if self._workdirs is not None:
            self._workdirs.append(workdir)
        return None

    def _build_command(self, problem, workdir, executable, context):
        return [str(executable), "-c", self._script, str(workdir)]

    def _parse_result(self, problem, workdir, run_result, context):
        assert (workdir / "ok").exists()
        count = len(problem.frequencies_hz)
        diagnostics = SolverDiagnostics(
            np.full((count, 1), True),
            np.ones((count, 1), int),
            np.full((count, 1), 1e-9),
            run_result.runtime_s,
        )
        return ForwardResult(
            problem.problem_hash,
            problem.frequencies_hz,
            problem.receivers.names,
            problem.components,
            np.full((1, count, 1), 1 + 2j),
            None,
            self.capabilities.name,
            self.capabilities.version,
            diagnostics,
        )


_SUCCESS_SCRIPT = (
    "import sys, pathlib; pathlib.Path(sys.argv[1], 'ok').write_text('done')"
)
_FAILURE_SCRIPT = "import sys; sys.exit(3)"
_SLEEP_SCRIPT = "import time; time.sleep(5)"


# --------------------------------------------------------------------------- #
# resolve_executable / probe_executable_version / make_availability_probe
# --------------------------------------------------------------------------- #


def test_resolve_executable_finds_a_direct_path(tmp_path):
    target = tmp_path / "tool"
    target.write_text("#!/bin/sh\n")
    assert resolve_executable(str(target)) == target


def test_resolve_executable_searches_extra_directories(tmp_path):
    target = tmp_path / "tool"
    target.write_text("x")
    resolved = resolve_executable("tool", search_paths=(str(tmp_path),))
    assert resolved == target


def test_resolve_executable_raises_when_missing():
    with pytest.raises(ExecutableNotFoundError, match="was not found"):
        resolve_executable("does-not-exist-xyz")


def test_make_availability_probe_reports_missing_and_present(tmp_path):
    missing_probe = make_availability_probe("does-not-exist-xyz")
    available, reason = missing_probe()
    assert available is False
    assert "was not found" in reason

    target = tmp_path / "tool"
    target.write_text("x")
    present_probe = make_availability_probe("tool", search_paths=(str(tmp_path),))
    assert present_probe() == (True, None)


def test_probe_executable_version_returns_none_for_missing_executable():
    assert probe_executable_version("does-not-exist-xyz") is None


def test_probe_executable_version_extracts_pattern_from_real_process():
    version = probe_executable_version(
        sys.executable,
        version_args=("--version",),
        version_pattern=r"Python (\S+)",
    )
    assert version is not None
    assert version[0].isdigit()


# --------------------------------------------------------------------------- #
# ExternalRunPolicy / ExternalRunResult
# --------------------------------------------------------------------------- #


def test_run_policy_rejects_invalid_fields():
    with pytest.raises(ValueError, match="executable cannot be empty"):
        ExternalRunPolicy("")
    with pytest.raises(ValueError, match="max_attempts"):
        ExternalRunPolicy("tool", max_attempts=0)
    with pytest.raises(ValueError, match="timeout_s"):
        ExternalRunPolicy("tool", timeout_s=-1.0)


def test_run_policy_allows_none_timeout_and_defaults():
    policy = ExternalRunPolicy("tool")
    assert policy.timeout_s is None
    assert policy.max_attempts == 1
    assert policy.extra_env == {}


def test_run_result_success_and_tail_and_to_dict():
    result = ExternalRunResult(
        ("a", "b"), 1, "out", "line1\nline2\nline3", 0.25, 1, "."
    )
    assert result.success is False
    assert result.tail(lines=2) == "line2\nline3"
    payload = result.to_dict()
    assert payload["returncode"] == 1
    assert payload["success"] is False


def test_run_result_rejects_invalid_fields():
    with pytest.raises(ValueError, match="command cannot be empty"):
        ExternalRunResult((), 0, "", "", 0.0, 1, ".")
    with pytest.raises(ValueError, match="attempt must be"):
        ExternalRunResult(("a",), 0, "", "", 0.0, 0, ".")


# --------------------------------------------------------------------------- #
# BaseExternalMaxwellAdapter end-to-end (real subprocess, no mocking)
# --------------------------------------------------------------------------- #


def test_scripted_adapter_solves_successfully_and_cleans_up_temp_workdir():
    workdirs: list[Path] = []
    adapter = _ScriptedAdapter(
        _SUCCESS_SCRIPT,
        workdirs=workdirs,
        run_policy=ExternalRunPolicy(sys.executable),
    )
    result = adapter.solve(_problem())
    assert result.success
    assert len(workdirs) == 1
    assert not workdirs[0].exists()  # cleaned up after success


def test_scripted_adapter_preserves_workdir_on_exhausted_failure():
    workdirs: list[Path] = []
    adapter = _ScriptedAdapter(
        _FAILURE_SCRIPT,
        workdirs=workdirs,
        run_policy=ExternalRunPolicy(sys.executable, max_attempts=1),
    )
    with pytest.raises(ExternalProcessError) as excinfo:
        adapter.solve(_problem())
    assert len(excinfo.value.attempts) == 1
    assert excinfo.value.attempts[0].returncode == 3
    assert workdirs[0].exists()  # preserved for inspection


def test_scripted_adapter_deletes_workdir_on_failure_when_disabled():
    workdirs: list[Path] = []
    adapter = _ScriptedAdapter(
        _FAILURE_SCRIPT,
        workdirs=workdirs,
        run_policy=ExternalRunPolicy(
            sys.executable, max_attempts=1, keep_workdir_on_failure=False
        ),
    )
    with pytest.raises(ExternalProcessError):
        adapter.solve(_problem())
    assert not workdirs[0].exists()


def test_scripted_adapter_retries_and_records_every_attempt():
    workdirs: list[Path] = []
    adapter = _ScriptedAdapter(
        _FAILURE_SCRIPT,
        workdirs=workdirs,
        run_policy=ExternalRunPolicy(
            sys.executable, max_attempts=3, retry_backoff_s=0.0
        ),
    )
    with pytest.raises(ExternalProcessError) as excinfo:
        adapter.solve(_problem())
    assert len(excinfo.value.attempts) == 3
    assert all(attempt.returncode == 3 for attempt in excinfo.value.attempts)


def test_scripted_adapter_times_out_and_raises_external_process_error():
    adapter = _ScriptedAdapter(
        _SLEEP_SCRIPT,
        run_policy=ExternalRunPolicy(sys.executable, timeout_s=0.2, max_attempts=1),
    )
    with pytest.raises(ExternalProcessError) as excinfo:
        adapter.solve(_problem())
    assert excinfo.value.attempts[0].returncode == -1
    assert "timed out" in excinfo.value.attempts[0].stderr


def test_scripted_adapter_reports_executable_not_found_unwrapped():
    adapter = _ScriptedAdapter(
        _SUCCESS_SCRIPT, run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    with pytest.raises(ExecutableNotFoundError):
        adapter.solve(_problem())


def test_scripted_adapter_never_deletes_caller_supplied_workdir(tmp_path):
    fixed = tmp_path / "run"
    adapter = _ScriptedAdapter(
        _SUCCESS_SCRIPT, run_policy=ExternalRunPolicy(sys.executable, workdir=str(fixed))
    )
    result = adapter.solve(_problem())
    assert result.success
    assert fixed.exists()  # caller-owned; never cleaned up
    assert (fixed / "ok").exists()
