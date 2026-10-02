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


@pytest.fixture(autouse=True)
def _isolate_config_singletons() -> None:
    """Restore every ``_SINGLETON_MAP`` singleton pyCSAMT config sections
    resolve to, after each test.

    ``config set`` / ``config style`` / ``config interp`` etc. run through
    Click's ``CliRunner``, i.e. in-process -- so they mutate pyCSAMT's real
    global config singletons (``PLOT_CONFIG``, ``PYCSAMT_STYLE``, ...), not
    a subprocess-local copy. In real CLI usage each invocation is its own
    short-lived process and this never matters; here, left unrestored, a
    value set by one test (e.g. ``config set plot.fmt pdf``, or
    ``config style dark``) leaks into every unrelated test that runs later
    in the same pytest session -- ``isolated_toml`` above only isolates the
    TOML file, not this in-memory state. Snapshot and restore every
    singleton ``apply_section`` can reach, regardless of which subcommand a
    given test invokes.
    """
    import copy

    from pycsamt.cli.commands.config._base import (
        _SINGLETON_MAP,
        get_singleton,
    )

    snapshots = {
        section: copy.deepcopy(vars(get_singleton(section)))
        for section in _SINGLETON_MAP
    }
    yield
    for section, snapshot in snapshots.items():
        obj = get_singleton(section)
        vars(obj).clear()
        vars(obj).update(snapshot)


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

    def test_load_all_config_skips_non_dict_section(self) -> None:
        # "plot" mapped to a non-dict value must be skipped, not raise.
        load_all_config({"plot": "not-a-dict", "control": {}})

    def test_apply_section_empty_kwargs_is_noop(self) -> None:
        apply_section("plot", {})

    def test_apply_section_unknown_section_is_noop(self) -> None:
        # No branch matches; must fall through silently.
        apply_section("totally_unknown_section", {"foo": "bar"})

    def test_apply_section_cli_flat(self) -> None:
        from pycsamt.api.cli.config import PYCSAMT_CLI

        original = PYCSAMT_CLI.log.level
        try:
            apply_section("cli", {"log__level": 2})
            assert PYCSAMT_CLI.log.level == 2
        finally:
            PYCSAMT_CLI.log.level = original

    def test_apply_section_cli_nested_legacy_form(self) -> None:
        from pycsamt.api.cli.config import PYCSAMT_CLI

        original = PYCSAMT_CLI.log.level
        try:
            apply_section("cli", {"log": {"level": 1}})
            assert PYCSAMT_CLI.log.level == 1
        finally:
            PYCSAMT_CLI.log.level = original

    def test_apply_section_style_preset_only(self) -> None:
        apply_section("style", {"preset": "publication"})

    def test_apply_section_style_extra_kwargs(self) -> None:
        # No preset key: goes straight to configure_style(**kw).
        apply_section("style", {"multiline__mode": "gradient"})

    def test_apply_section_section_view(self) -> None:
        apply_section("section_view", {"figsize": "10,8"})

    def test_apply_section_station(self) -> None:
        apply_section("station", {"density": 5})

    def test_apply_section_interp_preset_only(self) -> None:
        apply_section("interp", {"preset": "accessible"})

    def test_apply_section_interp_extra_kwargs(self) -> None:
        apply_section("interp", {"cmap": "viridis"})

    def test_apply_section_agent_provider(self) -> None:
        # provider is a read-only property; restoration is handled by the
        # autouse _isolate_config_singletons fixture (it snapshots/restores
        # __dict__ directly, bypassing the property).
        apply_section("agent", {"provider": "claude"})

    def test_apply_section_agent_model_only(self) -> None:
        from pycsamt.api.agents import AGENT_CONFIG

        apply_section("agent", {"model": "claude-x"})
        assert AGENT_CONFIG.model == "claude-x"

    def test_apply_section_agent_budget(self) -> None:
        from pycsamt.api.agents import AGENT_CONFIG

        apply_section("agent", {"budget_usd": 3.5})
        assert AGENT_CONFIG.remaining_usd is not None

    def test_apply_section_error_is_caught_and_warned(
        self, capsys: pytest.CaptureFixture
    ) -> None:
        apply_section("view", {"backend": "bogus_backend_xyz"})
        captured = capsys.readouterr()
        assert "Warning" in captured.err


class TestSingletonGetters:
    def test_get_singleton_all_known_sections(self) -> None:
        from pycsamt.cli.commands.config._base import (
            _SINGLETON_MAP,
            get_singleton,
        )

        for section in _SINGLETON_MAP:
            assert get_singleton(section) is not None

    def test_section_summary_unknown_section_is_unavailable(self) -> None:
        from pycsamt.cli.commands.config._base import section_summary

        assert section_summary("no_such_section") == "(unavailable)"

    def test_section_summary_falls_back_to_repr_when_summary_raises(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import pycsamt.cli.commands.config._base as _b
        from pycsamt.cli.commands.config._base import section_summary

        class _Weird:
            def summary(self):
                raise RuntimeError("boom")

            def __repr__(self):
                return "weird-repr"

        monkeypatch.setattr(_b, "get_singleton", lambda s: _Weird())
        assert section_summary("plot") == "weird-repr"

    def test_section_summary_falls_back_to_str_when_repr_raises(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import pycsamt.cli.commands.config._base as _b
        from pycsamt.cli.commands.config._base import section_summary

        class _Weird:
            def __repr__(self):
                raise RuntimeError("no repr")

            def __str__(self):
                return "weird-str"

        monkeypatch.setattr(_b, "get_singleton", lambda s: _Weird())
        assert section_summary("plot") == "weird-str"


class TestTomlIoEdgeCases:
    def test_read_toml_corrupted_file_returns_empty(
        self, isolated_toml: Path
    ) -> None:
        isolated_toml.write_text("not [ valid toml =", encoding="utf-8")
        assert _read_toml() == {}

    def test_read_toml_missing_both_backends_returns_empty(
        self, isolated_toml: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import builtins

        isolated_toml.write_text("[plot]\ndpi = 1\n", encoding="utf-8")
        real_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name in ("tomllib", "tomli"):
                raise ImportError("simulated: no toml backend")
            return real_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        assert _read_toml() == {}

    def test_write_toml_uses_tomli_w_when_available(
        self, isolated_toml: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import sys
        import types

        calls: dict = {}
        fake_mod = types.ModuleType("tomli_w")

        def fake_dump(data, fh):
            calls["data"] = data
            fh.write(b"[fake]\n")

        fake_mod.dump = fake_dump
        monkeypatch.setitem(sys.modules, "tomli_w", fake_mod)

        _write_toml({"plot": {"dpi": 9}})
        assert calls["data"] == {"plot": {"dpi": 9}}
        assert isolated_toml.read_bytes() == b"[fake]\n"

    def test_write_toml_fallback_skips_non_dict_and_empty_sections(
        self, isolated_toml: Path
    ) -> None:
        _write_toml(
            {"bad": "not-a-dict", "empty": {}, "good": {"k": "v"}}
        )
        text = isolated_toml.read_text(encoding="utf-8")
        assert "[good]" in text
        assert "[bad]" not in text
        assert "[empty]" not in text

    def test_write_toml_fallback_else_branch_for_other_types(
        self, isolated_toml: Path
    ) -> None:
        _write_toml({"sec": {"k": [1, 2, 3]}})
        text = isolated_toml.read_text(encoding="utf-8")
        assert "k =" in text


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
