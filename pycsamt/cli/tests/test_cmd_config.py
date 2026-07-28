# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt config`` command group.

Strategy
--------
All tests use an isolated TOML path so they never touch ~/.pycsamt.toml.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pycsamt.cli import main
from pycsamt.cli.commands.config._base import (
    _read_toml,
    _write_toml,
    apply_section,
    coerce,
    load_all_config,
    parse_key,
)

# ---------------------------------------------------------------------------
# Fixture: redirect TOML_PATH to a tmp file
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def isolated_toml(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect TOML_PATH to a temporary file for every test.

    config.py accesses TOML_PATH through _b.TOML_PATH (the _base module),
    so patching _base is sufficient.
    """
    fake_toml = tmp_path / ".pycsamt.toml"
    import pycsamt.cli.commands.config._base as _b

    monkeypatch.setattr(_b, "TOML_PATH", fake_toml)
    return fake_toml


# ---------------------------------------------------------------------------
# Unit tests for helpers
# ---------------------------------------------------------------------------


class TestParseKey:
    def test_two_parts(self) -> None:
        section, key = parse_key("plot.dpi")
        assert section == "plot"
        assert key == "dpi"

    def test_three_parts(self) -> None:
        section, key = parse_key("control.rho.view")
        assert section == "control"
        assert key == "rho__view"

    def test_four_parts(self) -> None:
        section, key = parse_key("style.mt.xy.color")
        assert section == "style"
        assert key == "mt__xy__color"

    def test_unknown_section_raises(self) -> None:
        import click

        with pytest.raises(click.BadParameter):
            parse_key("unknown.field")

    def test_missing_dot_raises(self) -> None:
        import click

        with pytest.raises(click.BadParameter):
            parse_key("plotdpi")


class TestCoerce:
    def test_bool_true(self) -> None:
        assert coerce("true") is True
        assert coerce("True") is True
        assert coerce("yes") is True

    def test_bool_false(self) -> None:
        assert coerce("false") is False
        assert coerce("no") is False

    def test_int(self) -> None:
        assert coerce("300") == 300
        assert isinstance(coerce("300"), int)

    def test_float(self) -> None:
        assert coerce("3.14") == pytest.approx(3.14)

    def test_string(self) -> None:
        assert coerce("pdf") == "pdf"
        assert coerce("log10") == "log10"


class TestTomlReadWrite:
    def test_round_trip(self, isolated_toml: Path) -> None:
        data = {"plot": {"dpi": 300, "fmt": "pdf"}}
        _write_toml(data)
        result = _read_toml()
        assert result == data

    def test_read_missing_returns_empty(self) -> None:
        result = _read_toml()
        assert result == {}


class TestApplySection:
    def test_apply_plot_dpi(self) -> None:
        from pycsamt.api.plot import PLOT_CONFIG

        original = PLOT_CONFIG.dpi
        apply_section("plot", {"dpi": 99})
        assert PLOT_CONFIG.dpi == 99
        # restore
        PLOT_CONFIG.dpi = original

    def test_apply_control_rho_view(self) -> None:
        from pycsamt.api.control import PYCSAMT_CONTROL

        original = PYCSAMT_CONTROL.rho.view
        apply_section("control", {"rho__view": "linear"})
        assert PYCSAMT_CONTROL.rho.view == "linear"
        # restore
        PYCSAMT_CONTROL.rho.view = original

    def test_apply_view_backend(self) -> None:
        from pycsamt.api.view.config import PYCSAMT_API_VIEW

        original = PYCSAMT_API_VIEW.backend
        apply_section("view", {"backend": "pandas"})
        assert PYCSAMT_API_VIEW.backend == "pandas"
        # restore
        apply_section("view", {"backend": original})

    def test_apply_pipe_error_mode(self) -> None:
        from pycsamt.api.pipe.config import PYCSAMT_PIPE

        original = PYCSAMT_PIPE.on_step_error
        apply_section("pipe", {"on_step_error": "skip"})
        assert PYCSAMT_PIPE.on_step_error == "skip"
        # restore
        PYCSAMT_PIPE.on_step_error = original

    def test_load_all_config(self) -> None:
        from pycsamt.api.plot import PLOT_CONFIG

        original = PLOT_CONFIG.dpi
        load_all_config({"plot": {"dpi": 77}})
        assert PLOT_CONFIG.dpi == 77
        PLOT_CONFIG.dpi = original


# ---------------------------------------------------------------------------
# CLI help wiring
# ---------------------------------------------------------------------------


class TestConfigGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "--help"])
        assert result.exit_code == 0
        for sub in (
            "list",
            "get",
            "set",
            "unset",
            "reset",
            "show",
            "env",
            "style",
            "interp",
            "agent",
        ):
            assert sub in result.output

    @pytest.mark.parametrize(
        "sub",
        [
            "list",
            "get",
            "set",
            "unset",
            "reset",
            "show",
            "env",
            "style",
            "interp",
        ],
    )
    def test_each_subcommand_help(self, runner: CliRunner, sub: str) -> None:
        result = runner.invoke(main, ["config", sub, "--help"])
        assert result.exit_code == 0

    def test_agent_subcommand_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "agent", "--help"])
        assert result.exit_code == 0
        assert "status" in result.output
        assert "set-key" in result.output


# ---------------------------------------------------------------------------
# config set / get / list / unset / reset
# ---------------------------------------------------------------------------


class TestConfigSetGet:
    def test_set_plot_dpi(self, runner: CliRunner, isolated_toml: Path) -> None:
        result = runner.invoke(main, ["config", "set", "plot.dpi", "250"])
        assert result.exit_code == 0
        assert "250" in result.output

    def test_set_persists_to_toml(self, runner: CliRunner, isolated_toml: Path) -> None:
        runner.invoke(main, ["config", "set", "plot.dpi", "250"])
        data = _read_toml()
        assert data.get("plot", {}).get("dpi") == 250

    def test_set_string_value(self, runner: CliRunner, isolated_toml: Path) -> None:
        result = runner.invoke(main, ["config", "set", "plot.fmt", "pdf"])
        assert result.exit_code == 0
        data = _read_toml()
        assert data["plot"]["fmt"] == "pdf"

    def test_set_bool_value(self, runner: CliRunner, isolated_toml: Path) -> None:
        result = runner.invoke(main, ["config", "set", "control.phase.wrap", "true"])
        assert result.exit_code == 0
        data = _read_toml()
        assert data["control"]["phase__wrap"] is True

    def test_set_dry_run_no_write(self, runner: CliRunner, isolated_toml: Path) -> None:
        result = runner.invoke(main, ["config", "set", "plot.dpi", "999", "--dry-run"])
        assert result.exit_code == 0
        assert "dry-run" in result.output.lower()
        assert not isolated_toml.exists()

    def test_get_existing_key(self, runner: CliRunner, isolated_toml: Path) -> None:
        from pycsamt.api.plot import PLOT_CONFIG

        result = runner.invoke(main, ["config", "get", "plot.dpi"])
        assert result.exit_code == 0
        assert str(PLOT_CONFIG.dpi) in result.output

    def test_get_json(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "get", "plot.dpi", "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "plot.dpi" in data

    def test_get_unknown_key_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "get", "plot.nosuchkey"])
        assert result.exit_code != 0

    def test_set_unknown_section_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "set", "nosection.key", "val"])
        assert result.exit_code != 0


class TestConfigList:
    def test_list_empty_config(self, runner: CliRunner, isolated_toml: Path) -> None:
        result = runner.invoke(main, ["config", "list"])
        assert result.exit_code == 0

    def test_list_after_set(self, runner: CliRunner, isolated_toml: Path) -> None:
        runner.invoke(main, ["config", "set", "plot.dpi", "300"])
        result = runner.invoke(main, ["config", "list"])
        assert result.exit_code == 0
        assert "dpi" in result.output or "300" in result.output

    def test_list_section_filter(self, runner: CliRunner, isolated_toml: Path) -> None:
        runner.invoke(main, ["config", "set", "plot.dpi", "300"])
        result = runner.invoke(main, ["config", "list", "plot"])
        assert result.exit_code == 0

    def test_list_json(self, runner: CliRunner, isolated_toml: Path) -> None:
        runner.invoke(main, ["config", "set", "plot.dpi", "300"])
        result = runner.invoke(main, ["config", "list", "plot", "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "plot" in data


class TestConfigUnsetReset:
    def test_unset_existing_key(self, runner: CliRunner, isolated_toml: Path) -> None:
        runner.invoke(main, ["config", "set", "plot.dpi", "300"])
        result = runner.invoke(main, ["config", "unset", "plot.dpi"])
        assert result.exit_code == 0
        data = _read_toml()
        assert "dpi" not in data.get("plot", {})

    def test_unset_nonexistent_warns(
        self, runner: CliRunner, isolated_toml: Path
    ) -> None:
        result = runner.invoke(main, ["config", "unset", "plot.dpi"])
        # should exit 0 with a warning (not a hard failure)
        assert result.exit_code == 0

    def test_reset_section(self, runner: CliRunner, isolated_toml: Path) -> None:
        runner.invoke(main, ["config", "set", "plot.dpi", "300"])
        result = runner.invoke(main, ["config", "reset", "plot", "--yes"])
        assert result.exit_code == 0
        data = _read_toml()
        assert "plot" not in data

    def test_reset_all(self, runner: CliRunner, isolated_toml: Path) -> None:
        runner.invoke(main, ["config", "set", "plot.dpi", "300"])
        result = runner.invoke(main, ["config", "reset", "--yes"])
        assert result.exit_code == 0
        assert not isolated_toml.exists()


class TestConfigShow:
    def test_show_no_error(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "show"])
        assert result.exit_code == 0

    def test_show_section(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "show", "plot"])
        assert result.exit_code == 0

    def test_show_json(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "show", "plot", "--format", "json"])
        assert result.exit_code == 0
        assert json.loads(result.output)  # valid JSON


class TestConfigEnv:
    def test_env_no_error(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "env"])
        assert result.exit_code == 0
        assert "PYCSAMT" in result.output or "API_KEY" in result.output

    def test_env_section_filter(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "env", "--section", "plot"])
        assert result.exit_code == 0
        assert "PYCSAMT_DPI" in result.output or "PYCSAMT_FMT" in result.output


class TestConfigPresets:
    def test_style_preset(self, runner: CliRunner, isolated_toml: Path) -> None:
        result = runner.invoke(main, ["config", "style", "publication"])
        assert result.exit_code == 0
        data = _read_toml()
        assert data.get("style", {}).get("preset") == "publication"

    def test_style_preset_no_persist(
        self, runner: CliRunner, isolated_toml: Path
    ) -> None:
        result = runner.invoke(main, ["config", "style", "dark", "--no-persist"])
        assert result.exit_code == 0
        assert not isolated_toml.exists()

    def test_interp_preset(self, runner: CliRunner, isolated_toml: Path) -> None:
        result = runner.invoke(main, ["config", "interp", "accessible"])
        assert result.exit_code == 0
        data = _read_toml()
        assert data.get("interp", {}).get("preset") == "accessible"


class TestConfigAgent:
    def test_agent_status(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "agent", "status"])
        assert result.exit_code == 0

    def test_agent_status_json(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["config", "agent", "status", "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "provider" in data

    @pytest.mark.parametrize("provider", ["claude", "openai", "gemini"])
    def test_set_key_instructions(self, runner: CliRunner, provider: str) -> None:
        result = runner.invoke(main, ["config", "agent", "set-key", provider])
        assert result.exit_code == 0
        # Should show the env var name
        expected = {
            "claude": "ANTHROPIC_API_KEY",
            "openai": "OPENAI_API_KEY",
            "gemini": "GOOGLE_API_KEY",
        }[provider]
        assert expected in result.output
