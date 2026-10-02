# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ``pycsamt convert`` command."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
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

    def test_edi_passthrough_directory_with_explicit_to_edi(
        self, runner: CliRunner, edi_dir: Path, tmp_path: Path
    ) -> None:
        """An explicit ``--to edi`` on a directory must still run the
        legacy batch path rather than being rejected as a single-file
        transfer-function conversion (regression guard)."""
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            [
                "convert",
                str(edi_dir),
                "--to",
                "edi",
                "--output-dir",
                str(out),
            ],
        )
        assert result.exit_code == 0
        n_in = len(list(edi_dir.glob("*.edi")))
        n_out = len(list(out.glob("*.edi")))
        assert n_out == n_in

    def test_directory_with_to_emtf_xml_is_rejected_clearly(
        self, runner: CliRunner, edi_dir: Path, tmp_path: Path
    ) -> None:
        """A target format other than EDI has no directory-batch support;
        it must fail with a clear message, not silently misbehave."""
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            [
                "convert",
                str(edi_dir),
                "--to",
                "emtf-xml",
                "--output-dir",
                str(out),
            ],
        )
        assert result.exit_code != 0
        assert "SOURCE must be a file" in result.output

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

    # ------------------------------------------------------------------
    # Jones J-file conversion (real bundled sample: data/j/kb0-s001.txt)
    # ------------------------------------------------------------------

    def test_j_file_success(
        self, runner: CliRunner, j_single_file: Path, tmp_path: Path
    ) -> None:
        src = tmp_path / "kb0-s001.j"
        src.write_text(j_single_file.read_text(encoding="utf-8"), encoding="utf-8")
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            ["convert", str(src), "--output-dir", str(out), "--format", "json"],
        )
        assert result.exit_code == 0, result.output
        data = json.loads(result.output)
        assert data[0]["status"] == "ok"
        assert data[0]["station"] == "KB0001"
        edi = out / "kb0-s001.edi"
        assert edi.exists()
        assert 'DATAID="KB0001"' in edi.read_text(encoding="utf-8")

    def test_j_file_verbose_echoes_progress(
        self, runner: CliRunner, j_single_file: Path, tmp_path: Path
    ) -> None:
        src = tmp_path / "kb0-s001.j"
        src.write_text(j_single_file.read_text(encoding="utf-8"), encoding="utf-8")
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            ["convert", str(src), "--output-dir", str(out), "-v"],
        )
        assert result.exit_code == 0, result.output
        assert "→" in result.output

    def test_j_file_malformed_reports_error(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        src = tmp_path / "bad.j"
        src.write_text("not a real j-file\n1 2 3\n", encoding="utf-8")
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            ["convert", str(src), "--output-dir", str(out), "--format", "json"],
        )
        assert result.exit_code != 0
        data = json.loads(result.output)
        assert data[0]["status"] == "error"

    def test_j_file_error_text_format_lists_message(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        src = tmp_path / "bad.j"
        src.write_text("garbage\n", encoding="utf-8")
        out = tmp_path / "out"
        result = runner.invoke(
            main, ["convert", str(src), "--output-dir", str(out)]
        )
        assert result.exit_code != 0
        assert "Errors" in result.output

    # ------------------------------------------------------------------
    # Zonge AVG conversion (real bundled sample data/avg/K1.AVG, K2.AVG)
    # ------------------------------------------------------------------

    def test_avg_file_success_multi_station(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        avg_src = Path("data/avg/K2.AVG")
        if not avg_src.exists():
            pytest.skip("data/avg/K2.AVG not found")
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            [
                "convert",
                str(avg_src),
                "--output-dir",
                str(out),
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0, result.output
        data = json.loads(result.output)
        assert data[0]["status"] == "ok"
        assert len(data[0]["stations"]) == 28
        assert len(list(out.glob("*.edi"))) == 28

    def test_avg_file_malformed_reports_error(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        src = tmp_path / "bad.avg"
        src.write_text("not a real avg file\ngarbage\n", encoding="utf-8")
        out = tmp_path / "out"
        result = runner.invoke(
            main, ["convert", str(src), "--output-dir", str(out), "--format", "json"]
        )
        assert result.exit_code != 0
        data = json.loads(result.output)
        assert data[0]["status"] == "error"

    def test_avg_csv_format(self, runner: CliRunner, tmp_path: Path) -> None:
        avg_src = Path("data/avg/K1.AVG")
        if not avg_src.exists():
            pytest.skip("data/avg/K1.AVG not found")
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            [
                "convert",
                str(avg_src),
                "--output-dir",
                str(out),
                "--format",
                "csv",
            ],
        )
        assert result.exit_code == 0, result.output
        lines = [l for l in result.output.splitlines() if l.strip()]
        assert len(lines) >= 2

    def test_avg_verbose_echoes_station_count(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        avg_src = Path("data/avg/K1.AVG")
        if not avg_src.exists():
            pytest.skip("data/avg/K1.AVG not found")
        out = tmp_path / "out"
        result = runner.invoke(
            main, ["convert", str(avg_src), "--output-dir", str(out), "-v"]
        )
        assert result.exit_code == 0, result.output
        assert "→" in result.output

    # ------------------------------------------------------------------
    # EDI pass-through verbose flag
    # ------------------------------------------------------------------

    def test_edi_passthrough_verbose_echoes_copy(
        self, runner: CliRunner, single_edi: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        result = runner.invoke(
            main, ["convert", str(single_edi), "--output-dir", str(out), "-v"]
        )
        assert result.exit_code == 0
        assert "copy" in result.output.lower()

    # ------------------------------------------------------------------
    # Single unsupported-extension file (not a directory)
    # ------------------------------------------------------------------

    def test_single_file_unsupported_extension_fails(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        src = tmp_path / "notes.txt"
        src.write_text("irrelevant", encoding="utf-8")
        out = tmp_path / "out"
        result = runner.invoke(
            main, ["convert", str(src), "--output-dir", str(out)]
        )
        assert result.exit_code != 0
        assert "No convertible files" in result.output

    # ------------------------------------------------------------------
    # Defensive "no converter for extension" branch: _SUPPORTED_EXTS and
    # _CONVERTERS are always kept in sync in real usage, so this is only
    # reachable by deliberately desyncing them.
    # ------------------------------------------------------------------

    def test_no_converter_for_extension_branch(
        self,
        runner: CliRunner,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        import pycsamt.cli.commands.convert as _convert_mod

        monkeypatch.setattr(
            _convert_mod,
            "_SUPPORTED_EXTS",
            _convert_mod._SUPPORTED_EXTS | {".zzz"},
        )
        src = tmp_path / "weird.zzz"
        src.write_text("irrelevant", encoding="utf-8")
        out = tmp_path / "out"
        result = runner.invoke(
            main,
            ["convert", str(src), "--output-dir", str(out), "--format", "json"],
        )
        assert result.exit_code != 0
        data = json.loads(result.output)
        assert data[0]["status"] == "error"
        assert "no converter" in data[0]["message"].lower()

    # ------------------------------------------------------------------
    # _fmt_csv helper — empty input
    # ------------------------------------------------------------------

    def test_fmt_csv_empty_results(self) -> None:
        from pycsamt.cli.commands.convert import _fmt_csv

        assert _fmt_csv([]) == ""
