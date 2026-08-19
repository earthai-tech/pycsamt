"""Tests for Occam1D execution outputs and text exports."""

from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DConfig,
    Occam1DData,
    Occam1DLog,
    Occam1DModel,
    Occam1DResponse,
    Occam1DResult,
    Occam1DRunner,
    Occam1DStartup,
)


def _make_run(path):
    path.mkdir()
    data = Occam1DData(
        [100.0, 10.0],
        [100.0, 200.0],
        [45.0, 50.0],
        [0.05, 0.05],
        [2.0, 2.0],
        station="S01",
    )
    data.write(path / "Occam1DData")
    model = Occam1DModel.build(4, 5.0, 500.0)
    model.write(path / "Occam1DModel")
    startup = Occam1DStartup.from_model(
        model, Occam1DConfig(n_layers=4, depth_max=500.0)
    )
    startup.write(path / "Startup")
    for iteration, rho in ((1, 50.0), (2, 80.0)):
        item = Occam1DStartup.from_model(model)
        item.iteration = iteration
        item.parameters[:] = np.log10(rho)
        item.write(path / f"ITER_{iteration:02d}")
    (path / "RESP_02.resp").write_text(
        "1 1 103 0.1 2.0 2.1 -1.0\n"
        "1 1 104 2.0 45.0 44.0 0.5\n",
        encoding="utf8",
    )
    (path / "occam1d.log").write_text(
        "** ITERATION 1 **\n"
        "AND IS = 2.0\n"
        "ROUGHNESS IS = 10\n"
        "** ITERATION 2 **\n"
        "AND IS = 1.0\n"
        "ROUGHNESS IS = 12\n",
        encoding="utf8",
    )


def test_response_converts_log_resistivity(tmp_path):
    path = tmp_path / "response.resp"
    path.write_text("1 1 103 0.1 2.0 2.1 -1.0\n", encoding="utf8")
    response = Occam1DResponse.read(path)
    observed, modeled = response.physical_values()
    assert observed[0] == pytest.approx(100.0)
    assert modeled[0] == pytest.approx(10**2.1)
    assert response.rms == pytest.approx(1.0)


def test_log_selects_lowest_rms(tmp_path):
    path = tmp_path / "run.log"
    path.write_text(
        "** ITERATION 1 **\nAND IS = 2.5\n"
        "** ITERATION 2 **\nAND IS = 0.9\n",
        encoding="utf8",
    )
    log = Occam1DLog.read(path)
    assert log.best_iteration == 2
    assert log.converged()


def test_result_loads_best_iteration_and_exports_text(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    result = Occam1DResult(run)
    assert result.iteration == 2
    np.testing.assert_allclose(result.resistivity, 80.0)
    outputs = result.export_text(tmp_path / "model-text")
    assert set(outputs) == {"model", "summary", "response", "iterations"}
    assert all(path.is_file() for path in outputs.values())


def test_runner_patches_startup_and_captures_exit_code(tmp_path, monkeypatch):
    run = tmp_path / "run"
    _make_run(run)
    binary = tmp_path / "Occam1D-test"
    binary.write_text("test executable marker", encoding="utf8")
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("subprocess.run", fake_run)
    runner = Occam1DRunner(run, binary_path=binary)
    code = runner.run(max_iterations=7, target_misfit=1.2)
    assert code == 0
    assert calls[0][0][1] == "Startup"
    startup = Occam1DStartup.read(run / "Startup")
    assert startup.max_iterations == 7
    assert startup.target_misfit == pytest.approx(1.2)
