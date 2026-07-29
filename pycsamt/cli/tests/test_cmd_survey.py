# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ``pycsamt survey`` command group."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from pycsamt.cli import main
from pycsamt.cli.survey import SurveyContext


class _FakeSite:
    def __init__(self, name: str) -> None:
        self.name = name


class _FakeSites:
    def __init__(self, names: list[str]) -> None:
        self._items = [_FakeSite(n) for n in names]

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self):
        return iter(self._items)


def _make_fake_sites(n: int = 5, names: list[str] | None = None) -> _FakeSites:
    names = names or [f"S{i:02d}" for i in range(1, n + 1)]
    return _FakeSites(names)


class TestSurveySet:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["survey", "set", "--help"])
        assert result.exit_code == 0
        assert "EDI_DIR" in result.output

    def test_set_writes_context(
        self, runner: CliRunner, isolated_home: Path, edi_dir: Path
    ) -> None:
        fake_sites = _make_fake_sites(3)
        with patch("pycsamt.cli.survey._build_sites", return_value=fake_sites):
            result = runner.invoke(main, ["survey", "set", str(edi_dir)])
        assert result.exit_code == 0
        ctx = SurveyContext.load()
        assert ctx is not None
        assert ctx.n_stations == 3

    def test_set_nonexistent_path_fails(
        self, runner: CliRunner, isolated_home: Path, tmp_path: Path
    ) -> None:
        missing = tmp_path / "no_such_dir"
        result = runner.invoke(main, ["survey", "set", str(missing)])
        assert result.exit_code != 0

    def test_set_force_flag_rebuilds(
        self, runner: CliRunner, isolated_home: Path, edi_dir: Path
    ) -> None:
        call_count = {"n": 0}
        fake_sites = _make_fake_sites(2)

        def counting_build(p, verbose=0):
            call_count["n"] += 1
            return fake_sites

        with patch("pycsamt.cli.survey._build_sites", side_effect=counting_build):
            runner.invoke(main, ["survey", "set", str(edi_dir)])
            runner.invoke(main, ["survey", "set", str(edi_dir)])  # should use cache
            runner.invoke(
                main, ["survey", "set", str(edi_dir), "--force"]
            )  # force rebuild

        assert call_count["n"] == 2


class TestSurveyShow:
    def test_show_no_context(self, runner: CliRunner, isolated_home: Path) -> None:
        result = runner.invoke(main, ["survey", "show"])
        assert result.exit_code == 0
        assert "No active survey" in result.output

    def test_show_with_context(
        self, runner: CliRunner, isolated_home: Path, edi_dir: Path
    ) -> None:
        fake_sites = _make_fake_sites(4, ["A1", "A2", "A3", "A4"])
        with patch("pycsamt.cli.survey._build_sites", return_value=fake_sites):
            runner.invoke(main, ["survey", "set", str(edi_dir)])

        result = runner.invoke(main, ["survey", "show"])
        assert result.exit_code == 0

    def test_show_json_format(
        self, runner: CliRunner, isolated_home: Path, edi_dir: Path
    ) -> None:
        fake_sites = _make_fake_sites(2)
        with patch("pycsamt.cli.survey._build_sites", return_value=fake_sites):
            runner.invoke(main, ["survey", "set", str(edi_dir)])

        result = runner.invoke(main, ["survey", "show", "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "survey_path" in data
        assert "n_stations" in data


class TestSurveyClear:
    def test_clear_no_context(self, runner: CliRunner, isolated_home: Path) -> None:
        result = runner.invoke(main, ["survey", "clear", "--yes"])
        assert result.exit_code == 0
        assert "nothing to clear" in result.output.lower()

    def test_clear_with_context(
        self, runner: CliRunner, isolated_home: Path, edi_dir: Path
    ) -> None:
        fake_sites = _make_fake_sites(3)
        with patch("pycsamt.cli.survey._build_sites", return_value=fake_sites):
            runner.invoke(main, ["survey", "set", str(edi_dir)])

        assert SurveyContext.load() is not None
        result = runner.invoke(main, ["survey", "clear", "--yes"])
        assert result.exit_code == 0
        assert SurveyContext.load() is None


class TestSurveyRebuild:
    def test_rebuild_no_context_fails(
        self, runner: CliRunner, isolated_home: Path
    ) -> None:
        result = runner.invoke(main, ["survey", "rebuild"])
        assert result.exit_code != 0

    def test_rebuild_triggers_fresh_build(
        self, runner: CliRunner, isolated_home: Path, edi_dir: Path
    ) -> None:
        call_count = {"n": 0}
        fake_sites = _make_fake_sites(2)

        def counting_build(p, verbose=0):
            call_count["n"] += 1
            return fake_sites

        with patch("pycsamt.cli.survey._build_sites", side_effect=counting_build):
            runner.invoke(main, ["survey", "set", str(edi_dir)])
            runner.invoke(main, ["survey", "rebuild", "--force"])

        # set + rebuild = 2 builds
        assert call_count["n"] == 2


class TestSurveyCache:
    def test_cache_list_empty(self, runner: CliRunner, isolated_home: Path) -> None:
        result = runner.invoke(main, ["survey", "cache", "list"])
        assert result.exit_code == 0
        assert "empty" in result.output.lower() or "No cache" in result.output

    def test_cache_list_after_set(
        self, runner: CliRunner, isolated_home: Path, edi_dir: Path
    ) -> None:
        fake_sites = _make_fake_sites(2)
        with patch("pycsamt.cli.survey._build_sites", return_value=fake_sites):
            runner.invoke(main, ["survey", "set", str(edi_dir)])

        result = runner.invoke(main, ["survey", "cache", "list"])
        assert result.exit_code == 0

    def test_cache_purge_active(
        self, runner: CliRunner, isolated_home: Path, edi_dir: Path
    ) -> None:
        fake_sites = _make_fake_sites(2)
        with patch("pycsamt.cli.survey._build_sites", return_value=fake_sites):
            runner.invoke(main, ["survey", "set", str(edi_dir)])

        result = runner.invoke(main, ["survey", "cache", "purge", "--yes"])
        assert result.exit_code == 0
        assert "purged" in result.output.lower()
