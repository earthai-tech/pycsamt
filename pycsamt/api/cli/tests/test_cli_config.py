from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.api.cli.config import (
    PYCSAMT_CLI,
    BuildConfig,
    LogConfig,
    OutputConfig,
    PyCSAMTCLI,
    configure_cli,
    reset_cli,
)


# ─────────────────────────────────────────────────────────────────────────
# Sub-config dataclasses
# ─────────────────────────────────────────────────────────────────────────


def test_log_config_accepts_valid_levels():
    for level in (0, 1, 2):
        assert LogConfig(level=level).level == level


def test_log_config_rejects_invalid_level():
    with pytest.raises(ValueError):
        LogConfig(level=99)


def test_output_config_accepts_valid_formats_and_converts_dir():
    cfg = OutputConfig(format="json", dir="some/dir")
    assert cfg.format == "json"
    assert cfg.dir == Path("some/dir")


def test_output_config_rejects_invalid_format():
    with pytest.raises(ValueError):
        OutputConfig(format="xml")


def test_build_config_rejects_non_positive_n_jobs():
    with pytest.raises(ValueError):
        BuildConfig(n_jobs=0)


def test_build_config_converts_cache_dir_to_path():
    cfg = BuildConfig(cache_dir="cache")
    assert cfg.cache_dir == Path("cache")


def test_build_config_cache_dir_none_by_default():
    assert BuildConfig().cache_dir is None


# ─────────────────────────────────────────────────────────────────────────
# PyCSAMTCLI
# ─────────────────────────────────────────────────────────────────────────


@pytest.fixture
def cli() -> PyCSAMTCLI:
    return PyCSAMTCLI()


def test_configure_sets_nested_attributes(cli):
    cli.configure(log__level=1, output__format="json", build__n_jobs=4)
    assert cli.log.level == 1
    assert cli.output.format == "json"
    assert cli.build.n_jobs == 4


def test_context_restores_settings_after_block(cli):
    original_level = cli.log.level
    with cli.context(log__level=2) as ctx:
        assert ctx is cli
        assert cli.log.level == 2
    assert cli.log.level == original_level


def test_context_restores_settings_on_exception(cli):
    with pytest.raises(RuntimeError):
        with cli.context(log__level=2):
            raise RuntimeError("boom")
    assert cli.log.level == 0


def test_context_with_no_overrides_changes_nothing(cli):
    with cli.context():
        assert cli.log.level == 0


def test_reset_restores_defaults(cli):
    cli.configure(log__level=2, build__n_jobs=8)
    cli.reset()
    assert cli.log.level == 0
    assert cli.build.n_jobs == 1


def test_summary_and_repr_contain_all_settings(cli):
    text = cli.summary()
    assert "log.level" in text
    assert "output.format" in text
    assert "build.n_jobs" in text
    assert repr(cli) == text


# ─────────────────────────────────────────────────────────────────────────
# load_env
# ─────────────────────────────────────────────────────────────────────────


def test_load_env_sets_verbose_level(cli, monkeypatch):
    monkeypatch.setenv("PYCSAMT_VERBOSE", "2")
    cli.load_env()
    assert cli.log.level == 2


def test_load_env_disables_color_on_any_nonempty_value(cli, monkeypatch):
    monkeypatch.setenv("PYCSAMT_NO_COLOR", "1")
    cli.load_env()
    assert cli.log.color is False


def test_load_env_sets_output_format(cli, monkeypatch):
    monkeypatch.setenv("PYCSAMT_OUTPUT", "csv")
    cli.load_env()
    assert cli.output.format == "csv"


def test_load_env_sets_output_dir(cli, monkeypatch, tmp_path):
    monkeypatch.setenv("PYCSAMT_OUTPUT_DIR", str(tmp_path))
    cli.load_env()
    assert cli.output.dir == tmp_path


def test_load_env_sets_build_n_jobs(cli, monkeypatch):
    monkeypatch.setenv("PYCSAMT_JOBS", "6")
    cli.load_env()
    assert cli.build.n_jobs == 6


def test_load_env_leaves_defaults_when_unset(cli, monkeypatch):
    for var in (
        "PYCSAMT_VERBOSE",
        "PYCSAMT_NO_COLOR",
        "PYCSAMT_OUTPUT",
        "PYCSAMT_OUTPUT_DIR",
        "PYCSAMT_JOBS",
    ):
        monkeypatch.delenv(var, raising=False)
    cli.load_env()
    assert cli.log.level == 0
    assert cli.log.color is True
    assert cli.output.format == "text"
    assert cli.build.n_jobs == 1


# ─────────────────────────────────────────────────────────────────────────
# Module-level singleton + convenience wrappers
# ─────────────────────────────────────────────────────────────────────────


def test_configure_cli_and_reset_cli_use_global_singleton():
    try:
        configure_cli(log__level=1)
        assert PYCSAMT_CLI.log.level == 1
        reset_cli()
        assert PYCSAMT_CLI.log.level == 0
    finally:
        reset_cli()
