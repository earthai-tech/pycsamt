"""Preflight, discovery, execution, and state tests for Occam1DRunner."""

import io
import subprocess
from types import SimpleNamespace

import pytest

from pycsamt.models.occam1d import (
    Occam1DData,
    Occam1DExecutionError,
    Occam1DModel,
    Occam1DRunner,
    Occam1DRunnerState,
    Occam1DStartup,
)


def _prepared(tmp_path, *, parameter_count=3):
    workdir = tmp_path / "run"
    workdir.mkdir()
    Occam1DData(
        [100.0], [100.0], [45.0], [0.05], [2.0]
    ).write(workdir / "Data")
    model = Occam1DModel([0.0, 10.0, 30.0], [100.0] * 3)
    model.write(workdir / "Model")
    Occam1DStartup(
        "Model", "Data", [2.0] * parameter_count
    ).write(workdir / "Startup")
    binary = tmp_path / "occam-test"
    binary.write_text("test executable", encoding="utf8")
    return workdir, binary


def test_preflight_loads_and_cross_validates_all_inputs(tmp_path):
    workdir, binary = _prepared(tmp_path)
    runner = Occam1DRunner(workdir, binary_path=binary)

    inputs = runner.preflight()

    assert runner.runner_state is Occam1DRunnerState.READY
    assert runner.is_ready
    assert runner.is_bound
    assert inputs["data"].n_data == 2
    assert inputs["model"].n_layers == 3
    assert inputs["startup"].n_parameters == 3


def test_preflight_rejects_missing_directory_and_referenced_file(tmp_path):
    missing = Occam1DRunner(tmp_path / "missing")
    with pytest.raises(FileNotFoundError, match="run directory"):
        missing.preflight()

    workdir, binary = _prepared(tmp_path)
    (workdir / "Data").unlink()
    runner = Occam1DRunner(workdir, binary_path=binary)
    with pytest.raises(FileNotFoundError, match="referenced"):
        runner.preflight()


def test_preflight_rejects_parameter_layer_mismatch(tmp_path):
    workdir, binary = _prepared(tmp_path, parameter_count=2)
    runner = Occam1DRunner(workdir, binary_path=binary)

    with pytest.raises(ValueError, match="layer count"):
        runner.preflight()


def test_explicit_binary_is_authoritative(tmp_path, monkeypatch):
    workdir, _ = _prepared(tmp_path)
    monkeypatch.setattr("shutil.which", lambda name: "fallback")
    runner = Occam1DRunner(
        workdir, binary_path=tmp_path / "missing-binary"
    )

    with pytest.raises(FileNotFoundError, match="Explicit"):
        runner.discover_binary()


def test_binary_discovery_checks_workdir_then_path(tmp_path, monkeypatch):
    workdir, _ = _prepared(tmp_path)
    local = workdir / "custom-occam"
    local.write_text("binary", encoding="utf8")
    runner = Occam1DRunner(workdir, binary_name="custom-occam")
    assert runner.discover_binary() == local.resolve()

    local.unlink()
    external = tmp_path / "external-occam"
    external.write_text("binary", encoding="utf8")
    monkeypatch.setattr("shutil.which", lambda name: str(external))
    runner = Occam1DRunner(workdir, binary_name="custom-occam")
    assert runner.discover_binary() == external.resolve()


def test_successful_run_records_command_logs_and_diagnostics(
    tmp_path, monkeypatch
):
    workdir, binary = _prepared(tmp_path)
    stream = io.StringIO()

    def fake_run(command, **kwargs):
        kwargs["stdout"].write("solver output\n")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("subprocess.run", fake_run)
    runner = Occam1DRunner(
        workdir, binary_path=binary, verbose=1, stream=stream
    )

    code = runner.run(extra_args=["--flag", 3])

    assert code == 0
    assert runner.succeeded
    assert runner.runner_state is Occam1DRunnerState.SUCCEEDED
    assert runner.last_command[-2:] == ("--flag", "3")
    assert runner.last_run.succeeded
    assert runner.last_run.elapsed_seconds >= 0
    assert runner.stdout_log.read_text(encoding="utf8") == "solver output\n"
    assert "Running Occam1D" in stream.getvalue()
    assert runner.diagnostics()["last_run"]["succeeded"] is True
    assert "succeeded" in runner.summary()


def test_nonzero_exit_and_check_exception_include_stderr(
    tmp_path, monkeypatch
):
    workdir, binary = _prepared(tmp_path)

    def fake_run(command, **kwargs):
        kwargs["stderr"].write("line one\nfatal solver failure\n")
        return SimpleNamespace(returncode=7)

    monkeypatch.setattr("subprocess.run", fake_run)
    runner = Occam1DRunner(workdir, binary_path=binary)

    assert runner.run() == 7
    assert runner.runner_state is Occam1DRunnerState.FAILED
    assert runner.stderr_tail(lines=1) == "fatal solver failure"
    with pytest.raises(Occam1DExecutionError, match="code 7") as error:
        runner.run(check=True)
    assert error.value.record.exit_code == 7
    assert "fatal solver failure" in str(error.value)


def test_timeout_records_sentinel_and_check_raises(tmp_path, monkeypatch):
    workdir, binary = _prepared(tmp_path)

    def timeout(command, **kwargs):
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])

    monkeypatch.setattr("subprocess.run", timeout)
    runner = Occam1DRunner(workdir, binary_path=binary)

    assert runner.run(timeout=0.1) == -9
    assert runner.runner_state is Occam1DRunnerState.TIMED_OUT
    assert runner.last_run.timed_out
    with pytest.raises(Occam1DExecutionError, match="timed out"):
        runner.run(timeout=0.1, check=True)


def test_launch_oserror_is_recorded_and_propagated(tmp_path, monkeypatch):
    workdir, binary = _prepared(tmp_path)

    def fail(command, **kwargs):
        raise PermissionError("not executable")

    monkeypatch.setattr("subprocess.run", fail)
    runner = Occam1DRunner(workdir, binary_path=binary)

    with pytest.raises(PermissionError, match="not executable"):
        runner.run()
    assert runner.runner_state is Occam1DRunnerState.FAILED
    assert runner.last_run.exit_code == -1


def test_runtime_overrides_are_persisted_after_binary_resolution(
    tmp_path, monkeypatch
):
    workdir, binary = _prepared(tmp_path)
    monkeypatch.setattr(
        "subprocess.run", lambda command, **kwargs: SimpleNamespace(returncode=0)
    )
    runner = Occam1DRunner(workdir, binary_path=binary)

    runner.run(max_iterations=9, target_misfit=1.25)
    startup = Occam1DStartup.read(workdir / "Startup")

    assert startup.max_iterations == 9
    assert startup.target_misfit == pytest.approx(1.25)


def test_extra_arguments_and_tail_count_are_validated(tmp_path):
    workdir, binary = _prepared(tmp_path)
    runner = Occam1DRunner(workdir, binary_path=binary)

    with pytest.raises(TypeError, match="iterable"):
        runner.run(extra_args="--bad")
    with pytest.raises(ValueError, match="null"):
        runner.run(extra_args=["bad\x00argument"])
    with pytest.raises(ValueError, match="positive"):
        runner.stderr_tail(lines=0)


@pytest.mark.parametrize("field", ["startup_file", "binary_name"])
def test_constructor_rejects_unsafe_names(tmp_path, field):
    values = {field: "bad\nname"}

    with pytest.raises(ValueError, match="control character"):
        Occam1DRunner(tmp_path, **values)
