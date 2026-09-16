# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ``pycsamt.cli._base`` — the root Click group.

Covers the pieces not exercised incidentally by other ``test_cmd_*``
modules: the UTF-8 stdio shim, the rich banner/help panel shown for a
bare ``pycsamt`` invocation, ``~/.pycsamt.toml`` loading at startup, and
the ``-v``/``--no-color`` root flags.
"""

from __future__ import annotations

import pathlib
from pathlib import Path

import pytest
from click.testing import CliRunner

import pycsamt.cli._base as _base
from pycsamt.cli import main


# ---------------------------------------------------------------------------
# _ensure_utf8_stdio
# ---------------------------------------------------------------------------


class _FakeStream:
    def __init__(self, raise_on_reconfigure: bool = False) -> None:
        self.reconfigured = False
        self._raise = raise_on_reconfigure

    def reconfigure(self, **kwargs) -> None:
        if self._raise:
            raise RuntimeError("cannot reconfigure")
        self.reconfigured = True


class TestEnsureUtf8Stdio:
    def test_reconfigures_streams_with_the_attribute(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        out, err = _FakeStream(), _FakeStream()
        monkeypatch.setattr(_base.sys, "stdout", out)
        monkeypatch.setattr(_base.sys, "stderr", err)
        _base._ensure_utf8_stdio()
        assert out.reconfigured and err.reconfigured

    def test_swallows_reconfigure_errors(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(
            _base.sys, "stdout", _FakeStream(raise_on_reconfigure=True)
        )
        monkeypatch.setattr(
            _base.sys, "stderr", _FakeStream(raise_on_reconfigure=True)
        )
        _base._ensure_utf8_stdio()  # must not raise

    def test_skips_streams_without_reconfigure_or_none(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        class Plain:
            pass

        monkeypatch.setattr(_base.sys, "stdout", Plain())
        monkeypatch.setattr(_base.sys, "stderr", None)
        _base._ensure_utf8_stdio()  # must not raise


# ---------------------------------------------------------------------------
# _load_toml_config
# ---------------------------------------------------------------------------


class TestLoadTomlConfig:
    def test_missing_file_returns_empty(self, tmp_path: Path) -> None:
        result = _base._load_toml_config(tmp_path / "nope.toml")
        assert result == {}

    def test_existing_file_is_parsed(self, tmp_path: Path) -> None:
        toml_path = tmp_path / "cfg.toml"
        toml_path.write_text("[plot]\ndpi = 111\n", encoding="utf-8")
        result = _base._load_toml_config(toml_path)
        assert result == {"plot": {"dpi": 111}}


# ---------------------------------------------------------------------------
# Bare invocation → rich banner + help panel
# ---------------------------------------------------------------------------


class TestBareInvocation:
    def test_no_args_shows_help_and_exits_zero(self, runner: CliRunner) -> None:
        result = runner.invoke(main, [])
        assert result.exit_code == 0
        assert "pyCSAMT" in result.output

    def test_no_args_lists_subcommands(self, runner: CliRunner) -> None:
        result = runner.invoke(main, [])
        assert "avg" in result.output
        assert "config" in result.output

    def test_narrow_terminal_skips_banner(
        self, monkeypatch: pytest.MonkeyPatch, runner: CliRunner
    ) -> None:
        class _NarrowConsole:
            class _Size:
                width = 10

            size = _Size()

            def print(self, *a, **k):
                pass

        monkeypatch.setattr(_base, "Console", lambda: _NarrowConsole())
        result = runner.invoke(main, [])
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# ~/.pycsamt.toml loading at CLI startup
# ---------------------------------------------------------------------------


class TestTomlConfigAtStartup:
    def test_toml_config_applied_at_startup(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, runner: CliRunner
    ) -> None:
        from pycsamt.api.plot import PLOT_CONFIG

        fake_home = tmp_path
        (fake_home / ".pycsamt.toml").write_text(
            "[plot]\ndpi = 234\n", encoding="utf-8"
        )
        monkeypatch.setattr(
            pathlib.Path, "home", classmethod(lambda cls: fake_home)
        )
        original = PLOT_CONFIG.dpi
        try:
            result = runner.invoke(main, ["config", "get", "plot.dpi"])
            assert result.exit_code == 0
            assert "234" in result.output
        finally:
            PLOT_CONFIG.dpi = original

    def test_bad_toml_does_not_crash_cli(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, runner: CliRunner
    ) -> None:
        fake_home = tmp_path
        (fake_home / ".pycsamt.toml").write_text(
            "[plot]\ndpi = 234\n", encoding="utf-8"
        )
        monkeypatch.setattr(
            pathlib.Path, "home", classmethod(lambda cls: fake_home)
        )

        def _boom(*a, **k):
            raise RuntimeError("bad config section")

        monkeypatch.setattr(
            "pycsamt.cli.commands.config.load_all_config", _boom
        )
        result = runner.invoke(main, ["avg", "--help"])
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# Root -v / -vv / --no-color flags
# ---------------------------------------------------------------------------


class TestRootFlags:
    def test_verbose_info_level(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["-v", "avg", "--help"])
        assert result.exit_code == 0

    def test_verbose_debug_level(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["-vv", "avg", "--help"])
        assert result.exit_code == 0

    def test_no_color_flag(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["--no-color", "avg", "--help"])
        assert result.exit_code == 0

    def test_version_flag(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "pycsamt" in result.output.lower()
