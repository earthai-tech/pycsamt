"""Tests for OccamRunner."""

from types import SimpleNamespace

import pytest


def test_runner_imports():
    from pycsamt.models.occam2d.runner import OccamRunner

    assert OccamRunner is not None


def test_runner_no_binary_raises(tmp_path):
    """discover_binary should raise when binary is absent and compile fails."""

    from pycsamt.models.occam2d.runner import OccamRunner

    runner = OccamRunner(workdir=tmp_path)
    # If gfortran is not on PATH, compile() will raise RuntimeError
    # If it is, it should still fail because _source/ is bundled — just
    # verify that discover_binary() raises something meaningful.
    with pytest.raises((FileNotFoundError, RuntimeError)):
        runner.discover_binary(auto_compile=False)


def test_runner_is_not_running_by_default(tmp_path):
    from pycsamt.models.occam2d.runner import OccamRunner

    runner = OccamRunner(workdir=tmp_path)
    assert not runner.is_running


def test_runner_forward_mode_uses_native_flag(tmp_path, monkeypatch):
    from pycsamt.models.occam2d.runner import OccamRunner

    binary = tmp_path / "Occam2D.exe"
    binary.write_bytes(b"solver")
    (tmp_path / "TruthStartup").write_text("startup", encoding="ascii")
    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["cwd"] = kwargs["cwd"]
        (tmp_path / "TruthForward.fwd").write_text(
            "forward",
            encoding="ascii",
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("subprocess.run", fake_run)
    runner = OccamRunner(
        workdir=tmp_path,
        binary_path=binary,
        startup_file="TruthStartup",
    )

    output = runner.run_forward("TruthForward", auto_compile=False)

    assert output == tmp_path / "TruthForward.fwd"
    assert captured["command"] == [
        str(binary),
        "-F",
        "TruthStartup",
        "TruthForward",
    ]
    assert captured["cwd"] == tmp_path
    assert runner.exit_code == 0


@pytest.mark.parametrize("output_root", ["", "nested/name", "name.fwd"])
def test_runner_forward_rejects_invalid_output_root(tmp_path, output_root):
    from pycsamt.models.occam2d.runner import OccamRunner

    runner = OccamRunner(workdir=tmp_path)
    with pytest.raises(ValueError, match="local filename root"):
        runner.run_forward(output_root)
