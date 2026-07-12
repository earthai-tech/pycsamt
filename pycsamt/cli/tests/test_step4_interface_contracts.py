from __future__ import annotations

import json

from click.testing import CliRunner

from pycsamt.cli import main


def test_root_cli_registers_public_command_groups():
    expected = {
        "avg",
        "config",
        "convert",
        "edi",
        "forward",
        "info",
        "interp",
        "invert",
        "jones",
        "map",
        "pipe",
        "site",
        "survey",
        "tdem",
        "transform",
    }

    assert expected <= set(main.commands)


def test_root_help_and_version_exit_zero():
    runner = CliRunner()

    help_result = runner.invoke(main, ["--help"])
    version_result = runner.invoke(main, ["--version"])

    assert help_result.exit_code == 0
    assert "pyCSAMT" in help_result.output or "pycsamt" in help_result.output
    assert version_result.exit_code == 0
    assert "pycsamt" in version_result.output.lower()


def test_unknown_command_fails_with_click_usage_message():
    result = CliRunner().invoke(main, ["does-not-exist"])

    assert result.exit_code != 0
    assert "No such command" in result.output
    assert result.exception is not None


def test_config_get_json_returns_machine_readable_payload():
    result = CliRunner().invoke(
        main,
        ["config", "get", "plot.dpi", "--format", "json"],
    )

    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert set(payload) == {"plot.dpi"}
    assert isinstance(payload["plot.dpi"], int)


def test_forward_scenarios_unknown_name_fails_cleanly():
    result = CliRunner().invoke(
        main,
        ["forward", "model", "geology", "--name", "missing_scenario"],
    )

    assert result.exit_code != 0
    assert "missing_scenario" in result.output
