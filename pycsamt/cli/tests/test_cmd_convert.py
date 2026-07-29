# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ``pycsamt convert`` command."""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from pycsamt.cli import main


class TestConvertCommand:
    # ------------------------------------------------------------------
    # Help
    # ------------------------------------------------------------------

    def test_help_exits_zero(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["convert", "--help"])
        assert result.exit_code == 0
        assert ".j" in result.output
        assert ".avg" in result.output

    # ------------------------------------------------------------------
    # Missing source → non-zero exit
    # ------------------------------------------------------------------

    def test_missing_source_fails(self, runner: CliRunner, tmp_path: Path) -> None:
        missing = tmp_path / "no_such"
        result = runner.invoke(main, ["convert", str(missing)])
        assert result.exit_code != 0

    # ------------------------------------------------------------------
    # EDI pass-through (copy)
    # ------------------------------------------------------------------

    def test_edi_passthrough_single_file(
        self, runner: CliRunner, single_edi: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            ["convert", str(single_edi), "--output-dir", str(out)],
        )
        assert result.exit_code == 0
        out_edi = out / (single_edi.stem + ".edi")
        assert out_edi.exists()

    def test_edi_passthrough_directory(
        self, runner: CliRunner, edi_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            ["convert", str(edi_dir), "--output-dir", str(out)],
        )
        assert result.exit_code == 0
        n_in = len(list(edi_dir.glob("*.edi")))
        n_out = len(list(out.glob("*.edi")))
        assert n_out == n_in

    # ------------------------------------------------------------------
    # Overwrite protection
    # ------------------------------------------------------------------

    def test_existing_output_skipped_without_overwrite(
        self, runner: CliRunner, single_edi: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        out.mkdir()
        existing = out / (single_edi.stem + ".edi")
        existing.write_text("existing content")

        runner.invoke(
            main,
            ["convert", str(single_edi), "--output-dir", str(out)],
        )
        # file content should NOT be replaced
        assert existing.read_text() == "existing content"

    def test_existing_output_replaced_with_overwrite(
        self, runner: CliRunner, single_edi: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        out.mkdir()
        existing = out / (single_edi.stem + ".edi")
        existing.write_text("OLD")

        runner.invoke(
            main,
            [
                "convert",
                str(single_edi),
                "--output-dir",
                str(out),
                "--overwrite",
            ],
        )
        assert existing.read_text() != "OLD"

    # ------------------------------------------------------------------
    # Dry run
    # ------------------------------------------------------------------

    def test_dry_run_writes_nothing(
        self, runner: CliRunner, edi_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            ["convert", str(edi_dir), "--output-dir", str(out), "--dry-run"],
        )
        assert result.exit_code == 0
        assert not out.exists() or not list(out.glob("*.edi"))

    def test_dry_run_lists_would_convert(
        self, runner: CliRunner, edi_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            ["convert", str(edi_dir), "--output-dir", str(out), "--dry-run"],
        )
        assert result.exit_code == 0
        n_edi = len(list(edi_dir.glob("*.edi")))
        # output should mention each file
        assert result.output.count("→") == n_edi

    # ------------------------------------------------------------------
    # Output formats
    # ------------------------------------------------------------------

    def test_json_format(
        self, runner: CliRunner, edi_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            [
                "convert",
                str(edi_dir),
                "--output-dir",
                str(out),
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert all("status" in r for r in data)

    def test_csv_format(self, runner: CliRunner, edi_dir: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            [
                "convert",
                str(edi_dir),
                "--output-dir",
                str(out),
                "--format",
                "csv",
            ],
        )
        assert result.exit_code == 0
        lines = [l for l in result.output.splitlines() if l.strip()]
        assert len(lines) >= 2  # header + at least 1 data row

    # ------------------------------------------------------------------
    # No supported files
    # ------------------------------------------------------------------

    def test_empty_dir_fails(self, runner: CliRunner, tmp_path: Path) -> None:
        src = tmp_path / "empty"
        src.mkdir()
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            ["convert", str(src), "--output-dir", str(out)],
        )
        assert result.exit_code != 0
