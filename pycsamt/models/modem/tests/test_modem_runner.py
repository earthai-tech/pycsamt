"""Self-contained tests for pycsamt.models.modem.runner.ModEmRunner.

Complements test_modem_phase5.py (command assembly + missing-binary
errors) by mocking subprocess.run/shutil.which to exercise the
success paths of run(), run_forward(), and _resolve_binary(), without
depending on a real ModEM executable or ModEMv626 example data.
"""

from __future__ import annotations

import subprocess
from unittest.mock import patch

import pytest

from pycsamt.models.modem.config import ModEmConfig
from pycsamt.models.modem.results import InversionResult
from pycsamt.models.modem.runner import ModEmRunner, _resolve_binary


class _FakeCompletedProcess:
    def __init__(self, returncode=0):
        self.returncode = returncode

    def check_returncode(self):
        if self.returncode != 0:
            raise subprocess.CalledProcessError(self.returncode, "cmd")


# ---------------------------------------------------------------------------
# _resolve_binary
# ---------------------------------------------------------------------------


class TestResolveBinary:
    def test_found_on_path(self, tmp_path):
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="/usr/bin/Mod3DMT"
        ):
            result = _resolve_binary("Mod3DMT", tmp_path)
        assert result.name == "Mod3DMT"

    def test_found_in_workdir(self, tmp_path):
        exe = tmp_path / "Mod3DMT"
        exe.write_text("binary")
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value=None
        ):
            result = _resolve_binary("Mod3DMT", tmp_path)
        assert result == exe

    def test_found_in_source_3d(self, tmp_path):
        d = tmp_path / "_source" / "3D"
        d.mkdir(parents=True)
        exe = d / "Mod3DMT"
        exe.write_text("binary")
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value=None
        ):
            result = _resolve_binary("Mod3DMT", tmp_path)
        assert result == exe

    def test_found_in_source_2d(self, tmp_path):
        d = tmp_path / "_source" / "2D"
        d.mkdir(parents=True)
        exe = d / "Mod2DMT"
        exe.write_text("binary")
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value=None
        ):
            result = _resolve_binary("Mod2DMT", tmp_path)
        assert result == exe

    def test_not_found_returns_none(self, tmp_path):
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value=None
        ):
            result = _resolve_binary("__nope__", tmp_path)
        assert result is None


# ---------------------------------------------------------------------------
# run() success paths
# ---------------------------------------------------------------------------


class TestRunSuccess:
    def test_run_serial_3d_loads_result(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            result = r.run("m0.ws", "d0.dat", "ctrl.inv")
        assert isinstance(result, InversionResult)
        cmd = mock_run.call_args.args[0]
        assert "Mod3DMT" in cmd
        assert mock_run.call_args.kwargs["cwd"] == str(tmp_path)

    def test_run_load_result_false_returns_none(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ):
            result = r.run("m0.ws", "d0.dat", "ctrl.inv", load_result=False)
        assert result is None

    def test_run_mpi(self, tmp_path):
        cfg = ModEmConfig(
            mode="3d", use_mpi=True, n_procs=4, mpi_command="mpirun",
            binary_3d="Mod3DMT",
        )
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            r.run("m0.ws", "d0.dat", "ctrl.inv", load_result=False)
        cmd = mock_run.call_args.args[0]
        assert cmd[:3] == ["mpirun", "-np", "4"]

    def test_run_use_mpi_override(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            r.run(
                "m0.ws", "d0.dat", "ctrl.inv",
                use_mpi=True, n_procs=2, load_result=False,
            )
        cmd = mock_run.call_args.args[0]
        assert cmd[:3] == [cfg.mpi_command, "-np", "2"]

    def test_run_with_covariance_writes_fwd_control(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            r.run(
                "m0.ws", "d0.dat", "ctrl.inv",
                covariance="cov.cov", load_result=False,
            )
        cmd = mock_run.call_args.args[0]
        assert cmd[-2:] == [cfg.fwd_control_file, "cov.cov"]
        assert (tmp_path / cfg.fwd_control_file).exists()

    def test_run_with_explicit_fwd_control_no_write(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            r.run(
                "m0.ws", "d0.dat", "ctrl.inv",
                fwd_control="my_fwd.ctrl", covariance="cov.cov",
                load_result=False,
            )
        cmd = mock_run.call_args.args[0]
        assert cmd[-2:] == ["my_fwd.ctrl", "cov.cov"]
        assert not (tmp_path / cfg.fwd_control_file).exists()

    def test_run_extra_args_appended(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            r.run(
                "m0.ws", "d0.dat", "ctrl.inv",
                extra_args=["--foo", "bar"], load_result=False,
            )
        cmd = mock_run.call_args.args[0]
        assert cmd[-2:] == ["--foo", "bar"]

    def test_run_nonzero_returncode_raises(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(1),
        ):
            with pytest.raises(subprocess.CalledProcessError):
                r.run("m0.ws", "d0.dat", "ctrl.inv")

    def test_run_verbose_logs_command(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg, verbose=1)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ):
            r.run("m0.ws", "d0.dat", "ctrl.inv", load_result=False)

    def test_run_2d_mode(self, tmp_path):
        cfg = ModEmConfig(mode="2d", use_mpi=False, binary_2d="Mod2DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod2DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            r.run("m0.rho", "d0.dat", "ctrl.inv", load_result=False)
        cmd = mock_run.call_args.args[0]
        assert "Mod2DMT" in cmd

    def test_run_timeout_passed_through(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            r.run("m0.ws", "d0.dat", "ctrl.inv", timeout=30, load_result=False)
        assert mock_run.call_args.kwargs["timeout"] == 30

    def test_run_timeout_expired_propagates(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            side_effect=subprocess.TimeoutExpired("cmd", 5),
        ):
            with pytest.raises(subprocess.TimeoutExpired):
                r.run("m0.ws", "d0.dat", "ctrl.inv", timeout=5)


# ---------------------------------------------------------------------------
# run_forward() success paths
# ---------------------------------------------------------------------------


class TestRunForwardSuccess:
    def test_run_forward_loads_result(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            result = r.run_forward("mi.ws", "d0.dat")
        assert isinstance(result, InversionResult)
        cmd = mock_run.call_args.args[0]
        assert "-F" in cmd

    def test_run_forward_mpi(self, tmp_path):
        cfg = ModEmConfig(
            mode="3d", use_mpi=True, n_procs=6, mpi_command="mpirun",
            binary_3d="Mod3DMT",
        )
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            r.run_forward("mi.ws", "d0.dat", load_result=False)
        cmd = mock_run.call_args.args[0]
        assert cmd[:3] == ["mpirun", "-np", "6"]

    def test_run_forward_load_result_false(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ):
            result = r.run_forward("mi.ws", "d0.dat", load_result=False)
        assert result is None

    def test_run_forward_nonzero_returncode_raises(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(1),
        ):
            with pytest.raises(subprocess.CalledProcessError):
                r.run_forward("mi.ws", "d0.dat")

    def test_run_forward_verbose_logs(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg, verbose=1)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod3DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ):
            r.run_forward("mi.ws", "d0.dat", load_result=False)

    def test_run_forward_2d_mode_override(self, tmp_path):
        cfg = ModEmConfig(mode="3d", use_mpi=False, binary_2d="Mod2DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value="Mod2DMT"
        ), patch(
            "pycsamt.models.modem.runner.subprocess.run",
            return_value=_FakeCompletedProcess(0),
        ) as mock_run:
            r.run_forward("m0.rho", "d0.dat", mode="2d", load_result=False)
        cmd = mock_run.call_args.args[0]
        assert "Mod2DMT" in cmd

    def test_run_forward_missing_binary_raises(self, tmp_path):
        cfg = ModEmConfig(mode="3d", binary_3d="__nonexistent_modem_binary__")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        with patch(
            "pycsamt.models.modem.runner.shutil.which", return_value=None
        ):
            with pytest.raises(FileNotFoundError):
                r.run_forward("m0.ws", "d0.dat")


# ---------------------------------------------------------------------------
# command() -- fwd_control default-name branch on covariance-only
# ---------------------------------------------------------------------------


class TestCommandFwdControlBranch:
    def test_command_with_explicit_fwd_control(self, tmp_path):
        cfg = ModEmConfig(mode="3d", binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        cmd = r.command(
            "m0.ws", "d0.dat", "c.inv",
            covariance="cov.cov", fwd_control="explicit.ctrl",
        )
        assert "explicit.ctrl" in cmd
        assert cmd.endswith("c.inv explicit.ctrl cov.cov")
        assert not (tmp_path / cfg.fwd_control_file).exists()

    def test_command_covariance_without_fwd_control_uses_default_name(
        self, tmp_path
    ):
        cfg = ModEmConfig(mode="3d", binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        cmd = r.command("m0.ws", "d0.dat", "c.inv", covariance="cov.cov")
        assert cfg.fwd_control_file in cmd
        assert not (tmp_path / cfg.fwd_control_file).exists()

    def test_command_no_covariance_no_fwd_control(self, tmp_path):
        cfg = ModEmConfig(mode="3d", binary_3d="Mod3DMT")
        r = ModEmRunner(workdir=tmp_path, config=cfg)
        cmd = r.command("m0.ws", "d0.dat", "c.inv")
        assert cfg.fwd_control_file not in cmd
