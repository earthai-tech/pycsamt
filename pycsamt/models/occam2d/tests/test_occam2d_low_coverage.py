"""Focused branch coverage for Occam2D runner and validators."""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.models.occam2d import runner as runner_module
from pycsamt.models.occam2d.runner import OccamRunner
from pycsamt.models.occam2d.plot import (
    PlotStation1DFit,
    _extract_1d_column,
    _extract_fit_data,
    _station_idx,
    plot_station_1d_fit,
)
from pycsamt.models.occam2d.validation import (
    OccamFileType,
    detect_file_type,
    is_data_file,
    is_iter_file,
    is_log_file,
    is_mesh_file,
    is_model_file,
    is_response_file,
    is_startup_file,
)


@pytest.mark.parametrize(
    ("text", "predicate", "kind"),
    [
        ("Format: OCCAM2MTDATA_1.0\n", is_data_file, OccamFileType.DATA),
        ("Format: OCCAM2MTMOD_1.0\n", is_model_file, OccamFileType.MODEL),
        ("Format: OCCAMITER_FLEX\nIteration: 0\n", is_startup_file, OccamFileType.STARTUP),
        ("Format: OCCAMITER_FLEX\nIteration: 2\n", is_iter_file, OccamFileType.ITER),
        ("mesh file\n", is_mesh_file, OccamFileType.MESH),
        ("Format: OCCAM2MTRESP_1.0\n", is_response_file, OccamFileType.RESPONSE),
        ("** iteration 3\n", is_log_file, OccamFileType.LOG),
    ],
)
def test_validation_recognizes_every_format(tmp_path, text, predicate, kind):
    path = tmp_path / kind
    path.write_text("\n" + text)
    assert predicate(path)
    assert detect_file_type(path) == kind


def test_validation_negative_and_malformed_cases(tmp_path):
    numeric_mesh = tmp_path / "numeric.mesh"
    numeric_mesh.write_text("4 5 anything\n")
    assert is_mesh_file(numeric_mesh)
    unknown = tmp_path / "unknown"
    unknown.write_text("hello\n")
    assert detect_file_type(unknown) == OccamFileType.UNKNOWN
    malformed = tmp_path / "bad.iter"
    malformed.write_text("Format: OCCAMITER_FLEX\nIteration: nope\n")
    assert not is_iter_file(malformed)
    assert not is_startup_file(tmp_path / "missing")
    assert not is_iter_file(tmp_path / "missing")
    assert not is_mesh_file(tmp_path / "missing")
    assert not is_model_file(tmp_path / "missing")
    assert not is_response_file(tmp_path / "missing")
    assert not is_log_file(tmp_path / "missing")


def test_runner_binary_discovery_paths(tmp_path, monkeypatch):
    explicit = tmp_path / "explicit"
    explicit.write_text("binary")
    assert OccamRunner(tmp_path, binary_path=explicit).discover_binary(False) == explicit

    local = tmp_path / runner_module._BINARY_NAME
    local.write_text("binary")
    assert OccamRunner(tmp_path).discover_binary(False) == local
    local.unlink()
    monkeypatch.setattr(runner_module.shutil, "which", lambda name: "/bin/occam")
    assert OccamRunner(tmp_path).discover_binary(False).name == "occam"
    monkeypatch.setattr(runner_module.shutil, "which", lambda name: None)
    compiled = tmp_path / "compiled"
    compiled.write_text("binary")
    runner = OccamRunner(tmp_path)
    monkeypatch.setattr(runner, "compile", lambda: compiled)
    assert runner.discover_binary(True) == compiled


def test_runner_compile_failures_and_success(tmp_path, monkeypatch):
    runner = OccamRunner(tmp_path)
    source = tmp_path / "source"
    monkeypatch.setattr(runner_module, "_SOURCE_DIR", source)
    with pytest.raises(FileNotFoundError, match="source directory"):
        runner.compile()
    source.mkdir()
    monkeypatch.setattr(runner_module.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match="compiler"):
        runner.compile("missing")
    monkeypatch.setattr(runner_module.shutil, "which", lambda name: "compiler")
    monkeypatch.setattr(runner_module.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=1, stderr="bad"))
    with pytest.raises(RuntimeError, match="Compilation failed"):
        runner.compile()
    monkeypatch.setattr(runner_module.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=0, stderr=""))
    with pytest.raises(RuntimeError, match="not produced"):
        runner.compile()
    binary = source / runner_module._BINARY_NAME
    binary.write_text("ok")
    assert runner.compile() == binary


def test_runner_sync_timeout_failure_and_patch(tmp_path, monkeypatch):
    binary = tmp_path / "solver"
    binary.write_text("ok")
    startup = tmp_path / "Startup"
    startup.write_text("Iterations to run: 1\nTarget Misfit: 2\nOther: keep\n")
    runner = OccamRunner(tmp_path, binary_path=binary, verbose=True)
    monkeypatch.setattr(runner_module.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=3))
    assert runner.run(max_iter=7, target_misfit=1.25, auto_compile=False) == 3
    text = startup.read_text()
    assert "Iterations to run:  7" in text and "Target Misfit:      1.2500" in text

    def timeout(*args, **kwargs):
        raise runner_module.subprocess.TimeoutExpired(args[0], 1)
    monkeypatch.setattr(runner_module.subprocess, "run", timeout)
    assert runner.run(timeout=1, auto_compile=False) == -9
    with pytest.raises(FileNotFoundError, match="Startup"):
        OccamRunner(tmp_path, startup_file="Absent")._patch_startup(1, None)


def test_runner_forward_errors_async_and_wait(tmp_path, monkeypatch):
    binary = tmp_path / "solver"
    binary.write_text("ok")
    runner = OccamRunner(tmp_path, binary_path=binary)
    with pytest.raises(FileNotFoundError, match="Startup"):
        runner.run_forward("Forward", auto_compile=False)
    (tmp_path / "Startup").write_text("startup")
    monkeypatch.setattr(runner_module.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=2))
    with pytest.raises(RuntimeError, match="exited with code"):
        runner.run_forward("Forward", auto_compile=False)
    monkeypatch.setattr(runner_module.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=0))
    with pytest.raises(RuntimeError, match="did not create"):
        runner.run_forward("Forward", auto_compile=False)

    class FakeProcess:
        pid = 123
        def poll(self): return None
        def wait(self): return 0
    monkeypatch.setattr(runner_module.subprocess, "Popen", lambda *a, **k: FakeProcess())
    process = runner.run_async(auto_compile=False)
    assert process is runner.process and runner.is_running
    assert runner.wait() == 0
    with pytest.raises(RuntimeError, match="No async"):
        OccamRunner(tmp_path).wait()


def _plot_result():
    response_rows = []
    data_rows = []
    for code, obs, mod, err in [(1, 2.0, 2.1, 0.1), (2, 45, 46, 2), (5, 2.3, 2.2, 0), (6, 40, 41, 0)]:
        for freq in (1, 2):
            response_rows.append([1, freq, code, 0, obs, mod, 0])
            data_rows.append([1, freq, code, obs, err])
    data = SimpleNamespace(
        sites=["SITE01"],
        frequencies=np.array([1.0, 10.0]),
        offsets=np.array([500.0]),
        data_blocks=np.asarray(data_rows, dtype=float),
    )
    return SimpleNamespace(
        response=SimpleNamespace(data=np.asarray(response_rows), n_data=len(response_rows)),
        data=data,
        mesh=SimpleNamespace(x_nodes=np.array([0, 1000]), z_nodes=np.array([0, 100, 500, 1000])),
        rho_2d=np.array([[1.0], [np.nan], [2.0]]),
        best_iter=SimpleNamespace(iteration=4),
    )


def test_station_1d_plot_helpers_and_full_plot():
    result = _plot_result()
    assert _station_idx(result, "SITE01") == 1
    assert _station_idx(result, 1) == 1
    with pytest.raises(ValueError, match="not found"):
        _station_idx(result, "absent")
    with pytest.raises(RuntimeError, match="No site names"):
        _station_idx(SimpleNamespace(data=None), "SITE01")
    fit = _extract_fit_data(result, 1)
    assert set(fit) == {"TE", "TM"}
    ztop, zbot, rho = _extract_1d_column(result, 1)
    assert len(ztop) == len(zbot) == len(rho) == 3
    assert _extract_fit_data(SimpleNamespace(response=None), 1) == {}
    assert _extract_1d_column(SimpleNamespace(), 1) == (None, None, None)
    no_offsets = SimpleNamespace(mesh=result.mesh, rho_2d=result.rho_2d, data=SimpleNamespace(offsets=np.array([])))
    assert _extract_1d_column(no_offsets, 1) == (None, None, None)

    fig = plot_station_1d_fit(
        result, station="SITE01", depth_max=0.6, rho_lim=(1, 1000),
        phase_lim=(30, 50), rho_depth_lim=(1, 200), title="Fit",
    )
    assert len(fig.axes) == 3
    plt.close(fig)
    with pytest.raises(RuntimeError, match="no InversionResult"):
        PlotStation1DFit().plot()


def test_station_1d_plot_without_model_or_matching_modes():
    result = _plot_result()
    result.mesh = None
    fig = PlotStation1DFit(result, modes=["XX"], phase_lim=(40, 50)).plot()
    assert "no model" in fig.axes[2].texts[0].get_text()
    plt.close(fig)
