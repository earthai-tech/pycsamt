# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ``pycsamt info`` command."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from pycsamt.cli import main


class TestInfoCommand:
    # ------------------------------------------------------------------
    # Help
    # ------------------------------------------------------------------

    def test_help_exits_zero(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["info", "--help"])
        assert result.exit_code == 0
        assert "EDI_FILE_OR_DIR" in result.output

    # ------------------------------------------------------------------
    # No argument + no active survey → UsageError
    # ------------------------------------------------------------------

    def test_no_args_no_survey_raises(
        self, runner: CliRunner, isolated_home: Path
    ) -> None:
        result = runner.invoke(main, ["info"])
        assert result.exit_code != 0
        assert "No active survey" in result.output or "No active survey" in (
            result.exception or ""
        )

    # ------------------------------------------------------------------
    # Explicit path — live EDI files
    # ------------------------------------------------------------------

    def test_explicit_single_edi(
        self, runner: CliRunner, single_edi: Path
    ) -> None:
        result = runner.invoke(main, ["info", str(single_edi)])
        assert result.exit_code == 0
        assert "Station" in result.output or "Frequencies" in result.output

    def test_explicit_edi_dir_text(
        self, runner: CliRunner, edi_dir: Path
    ) -> None:
        result = runner.invoke(main, ["info", str(edi_dir)])
        assert result.exit_code == 0
        len(list(edi_dir.glob("*.edi")))
        # At least one station block in the output
        assert (
            result.output.count("Station") >= 1
            or result.output.count("File") >= 1
        )

    def test_explicit_edi_dir_json(
        self, runner: CliRunner, edi_dir: Path
    ) -> None:
        import json

        result = runner.invoke(
            main, ["info", str(edi_dir), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) == len(list(edi_dir.glob("*.edi")))

    def test_explicit_edi_dir_csv(
        self, runner: CliRunner, edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["info", str(edi_dir), "--format", "csv"]
        )
        assert result.exit_code == 0
        lines = [l for l in result.output.splitlines() if l.strip()]
        assert len(lines) >= 2  # header + at least one data row

    # ------------------------------------------------------------------
    # Station filter
    # ------------------------------------------------------------------

    def test_station_filter_matching(
        self, runner: CliRunner, edi_dir: Path
    ) -> None:
        edi_files = sorted(edi_dir.glob("*.edi"))
        first_stem = edi_files[0].stem.upper()
        result = runner.invoke(
            main, ["info", str(edi_dir), "--stations", first_stem]
        )
        assert result.exit_code == 0

    def test_station_filter_no_match(
        self, runner: CliRunner, edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["info", str(edi_dir), "--stations", "DOESNOTEXIST_XYZ"],
        )
        assert result.exit_code != 0

    # ------------------------------------------------------------------
    # Verbosity
    # ------------------------------------------------------------------

    def test_verbose_flag(self, runner: CliRunner, edi_dir: Path) -> None:
        result = runner.invoke(main, ["info", str(edi_dir), "-v"])
        assert result.exit_code == 0

    # ------------------------------------------------------------------
    # --survey flag resolves correctly
    # ------------------------------------------------------------------

    def test_survey_flag_uses_given_dir(
        self, runner: CliRunner, edi_dir: Path, isolated_home: Path
    ) -> None:
        from unittest.mock import MagicMock

        fake_sites = MagicMock()
        fake_sites.__len__ = lambda _: 0
        fake_sites.__iter__ = lambda _: iter([])
        # With an explicit --survey, info should at least not crash on path resolution
        result = runner.invoke(main, ["info", "--survey", str(edi_dir)])
        # May exit 1 if no paths resolved, but should not raise Python exception
        assert result.exception is None or isinstance(
            result.exception, SystemExit
        )
