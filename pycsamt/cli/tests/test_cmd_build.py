# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt build`` command group.

These scripts compile external Fortran solvers, which is not something a
unit test suite should attempt. Coverage here is limited to what is safe
and deterministic everywhere: help wiring, script discovery, and the
locate/dispatch logic in ``_base.py`` -- exercised via ``--help``, which
every underlying script handles without touching a toolchain or the
filesystem.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from pycsamt.cli import main
from pycsamt.cli.commands.build import _base as build_base
from pycsamt.cli.commands.build._base import (
    find_bash,
    run_solver_script,
    solver_build_dir,
)

_SOLVERS = ("modem2d", "modem3d", "occam2d", "mare2dem")


class TestBuildGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["build", "--help"])
        assert result.exit_code == 0
        for sub in _SOLVERS:
            assert sub in result.output

    @pytest.mark.parametrize("sub", _SOLVERS)
    def test_each_subcommand_help_forwards_to_script(
        self, runner: CliRunner, sub: str, capfd: pytest.CaptureFixture[str]
    ) -> None:
        if find_bash() is None:
            pytest.skip("no bash-compatible shell on PATH")
        result = runner.invoke(main, ["build", sub, "--help"])
        assert result.exit_code == 0
        # The script's help text is written by a real subprocess directly
        # to the inherited file descriptor, so it shows up at the fd level
        # (capfd) rather than in CliRunner's own sys.stdout capture.
        out = capfd.readouterr().out
        assert sub in out
        assert "Usage" in out or "usage" in out


class TestSolverBuildDir:
    def test_directory_exists(self) -> None:
        assert solver_build_dir().is_dir()

    @pytest.mark.parametrize("sub", _SOLVERS)
    def test_each_script_exists(self, sub: str) -> None:
        assert (solver_build_dir() / f"{sub}.sh").is_file()

    def test_dispatcher_exists(self) -> None:
        assert (solver_build_dir() / "build.sh").is_file()


class TestRunSolverScript:
    def test_unknown_solver_reports_missing_script(self) -> None:
        code = run_solver_script("does-not-exist", ())
        assert code == 1

    def test_missing_bash_is_reported_not_raised(self, monkeypatch) -> None:
        monkeypatch.setattr(
            "pycsamt.cli.commands.build._base.find_bash", lambda: None
        )
        code = run_solver_script("occam2d", ())
        assert code == 1

    def test_keyboard_interrupt_returns_130(self, monkeypatch) -> None:
        def _raise(cmd):
            raise KeyboardInterrupt

        monkeypatch.setattr(build_base.subprocess, "call", _raise)
        code = run_solver_script("occam2d", ())
        assert code == 130


class TestFindBash:
    def test_returns_first_candidate_on_path(self, monkeypatch) -> None:
        monkeypatch.setattr(build_base.shutil, "which", lambda c: f"/usr/bin/{c}")
        assert find_bash() == "/usr/bin/bash"

    def test_falls_through_to_second_candidate(self, monkeypatch) -> None:
        def _which(candidate):
            return None if candidate == "bash" else f"/usr/bin/{candidate}"

        monkeypatch.setattr(build_base.shutil, "which", _which)
        assert find_bash() == "/usr/bin/sh"

    def test_windows_fallback_found(self, monkeypatch) -> None:
        monkeypatch.setattr(build_base.shutil, "which", lambda c: None)
        monkeypatch.setattr(build_base.sys, "platform", "win32")
        monkeypatch.setattr(build_base.Path, "exists", lambda self: True)
        assert find_bash() == build_base._WINDOWS_BASH_CANDIDATES[0]

    def test_windows_fallback_not_found(self, monkeypatch) -> None:
        monkeypatch.setattr(build_base.shutil, "which", lambda c: None)
        monkeypatch.setattr(build_base.sys, "platform", "win32")
        monkeypatch.setattr(build_base.Path, "exists", lambda self: False)
        assert find_bash() is None

    def test_non_windows_no_fallback(self, monkeypatch) -> None:
        monkeypatch.setattr(build_base.shutil, "which", lambda c: None)
        monkeypatch.setattr(build_base.sys, "platform", "linux")
        assert find_bash() is None
