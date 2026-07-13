"""Phase 5 tests — ModEmRunner (command assembly, no subprocess execution)."""

from pathlib import Path

import pytest

_EX3D = (
    Path(__file__).parents[4]
    / "ModEMv626"
    / "ModEM"
    / "examples"
    / "3D_MT"
    / "BLOCK2"
)


# ---------------------------------------------------------------------------
# command() — pure string assembly, no binary required
# ---------------------------------------------------------------------------


def test_runner_command_serial_3d(tmp_path):
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.runner import ModEmRunner

    cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
    r = ModEmRunner(workdir=tmp_path, config=cfg)
    cmd = r.command("m0.ws", "d0.dat", "ctrl.inv")
    assert "Mod3DMT" in cmd
    assert "-I" in cmd
    assert "NLCG" in cmd
    assert "m0.ws" in cmd
    assert "d0.dat" in cmd
    assert "ctrl.inv" in cmd
    assert "mpirun" not in cmd


def test_runner_command_mpi_3d(tmp_path):
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.runner import ModEmRunner

    cfg = ModEmConfig(
        mode="3d",
        use_mpi=True,
        n_procs=8,
        mpi_command="mpirun",
        binary_3d="Mod3DMT",
    )
    r = ModEmRunner(workdir=tmp_path, config=cfg)
    cmd = r.command("m0.ws", "d0.dat", "ctrl.inv")
    assert "mpirun" in cmd
    assert "-np" in cmd
    assert "8" in cmd
    assert "Mod3DMT" in cmd


def test_runner_command_serial_2d(tmp_path):
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.runner import ModEmRunner

    cfg = ModEmConfig(mode="2d", use_mpi=False, binary_2d="Mod2DMT")
    r = ModEmRunner(workdir=tmp_path, config=cfg)
    cmd = r.command("m0.rho", "d0.dat", "ctrl.inv", mode="2d")
    assert "Mod2DMT" in cmd
    assert "m0.rho" in cmd


def test_runner_command_with_covariance(tmp_path):
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.runner import ModEmRunner

    cfg = ModEmConfig(mode="3d", binary_3d="Mod3DMT")
    r = ModEmRunner(workdir=tmp_path, config=cfg)
    cmd = r.command("m0.ws", "d0.dat", "ctrl.inv", covariance="model.cov")
    assert "model.cov" in cmd


def test_runner_command_mode_override(tmp_path):
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.runner import ModEmRunner

    cfg = ModEmConfig(mode="3d", binary_2d="Mod2DMT")
    r = ModEmRunner(workdir=tmp_path, config=cfg)
    cmd = r.command("m0.rho", "d0.dat", "ctrl.inv", mode="2d")
    assert "Mod2DMT" in cmd


# ---------------------------------------------------------------------------
# run() — binary not available, expect FileNotFoundError
# ---------------------------------------------------------------------------


def test_runner_run_missing_binary_raises(tmp_path):
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.runner import ModEmRunner

    cfg = ModEmConfig(mode="3d", binary_3d="__nonexistent_modem_binary__")
    r = ModEmRunner(workdir=tmp_path, config=cfg)
    with pytest.raises(FileNotFoundError):
        r.run("m0.ws", "d0.dat", "ctrl.inv")


def test_runner_forward_missing_binary_raises(tmp_path):
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.runner import ModEmRunner

    cfg = ModEmConfig(mode="3d", binary_3d="__nonexistent_modem_binary__")
    r = ModEmRunner(workdir=tmp_path, config=cfg)
    with pytest.raises(FileNotFoundError):
        r.run_forward("m0.ws", "d0.dat")


# ---------------------------------------------------------------------------
# workdir creation
# ---------------------------------------------------------------------------


def test_runner_creates_workdir(tmp_path):
    from pycsamt.models.modem.runner import ModEmRunner

    new_dir = tmp_path / "new_run"
    ModEmRunner(workdir=new_dir)
    assert new_dir.is_dir()
