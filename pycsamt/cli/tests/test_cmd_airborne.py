# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt airborne`` command group.

Test strategy
-------------
Uses the real, bundled EMTF-XML sample surveys under ``data/ZTEM``,
``data/mobileMT``, and ``data/AFMAG`` (see the ``*_xml_dir`` fixtures in
``conftest.py``) rather than mocks -- these are small, fast to load, and
exercise the real ``ensure_asites`` / ``emtools`` code paths end to end.
Tests skip automatically when a sample directory is absent.
"""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from pycsamt.cli import main

# ---------------------------------------------------------------------------
# pycsamt airborne  (group help)
# ---------------------------------------------------------------------------


class TestAirborneGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["airborne", "--help"])
        assert result.exit_code == 0
        for sub in ("info", "diagnose"):
            assert sub in result.output

    def test_info_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["airborne", "info", "--help"])
        assert result.exit_code == 0
        for opt in ("SOURCE", "--site", "--format"):
            assert opt in result.output

    def test_diagnose_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["airborne", "diagnose", "--help"])
        assert result.exit_code == 0
        for opt in ("SOURCE", "--technology", "--line", "--component"):
            assert opt in result.output


# ---------------------------------------------------------------------------
# pycsamt airborne info
# ---------------------------------------------------------------------------


class TestAirborneInfo:
    def test_nonexistent_source_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["airborne", "info", "/nonexistent/path"])
        assert result.exit_code != 0

    def test_ztem_collection_text(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        result = runner.invoke(main, ["airborne", "info", str(ztem_xml_dir)])
        assert result.exit_code == 0
        assert "ztem" in result.output
        assert "Sites:" in result.output

    def test_ztem_collection_json(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["airborne", "info", str(ztem_xml_dir), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["n_sites"] > 0
        assert "ztem" in data["technologies"]
        assert len(data["sites"]) == data["n_sites"]

    def test_ztem_collection_csv(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["airborne", "info", str(ztem_xml_dir), "--format", "csv"]
        )
        assert result.exit_code == 0
        first_line = result.output.strip().splitlines()[0]
        assert "technology" in first_line

    def test_single_site(self, runner: CliRunner, ztem_xml_dir: Path) -> None:
        first_xml = sorted(ztem_xml_dir.glob("*.xml"))[0]
        site_name = first_xml.stem
        result = runner.invoke(
            main,
            ["airborne", "info", str(ztem_xml_dir), "--site", site_name],
        )
        assert result.exit_code == 0
        assert site_name in result.output

    def test_single_site_not_found(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["airborne", "info", str(ztem_xml_dir), "--site", "NOPE"],
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower()

    def test_single_file_source(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        first_xml = sorted(ztem_xml_dir.glob("*.xml"))[0]
        result = runner.invoke(main, ["airborne", "info", str(first_xml)])
        assert result.exit_code == 0
        assert "Sites:        1" in result.output

    def test_mobilemt_technology_detected(
        self, runner: CliRunner, mobilemt_xml_dir: Path
    ) -> None:
        result = runner.invoke(main, ["airborne", "info", str(mobilemt_xml_dir)])
        assert result.exit_code == 0
        assert "mobilemt" in result.output

    def test_afmag_original_technology_detected(
        self, runner: CliRunner, afmag_original_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["airborne", "info", str(afmag_original_xml_dir)]
        )
        assert result.exit_code == 0
        assert "afmag_original" in result.output

    def test_afmag_airmt_technology_detected(
        self, runner: CliRunner, afmag_airmt_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["airborne", "info", str(afmag_airmt_xml_dir)]
        )
        assert result.exit_code == 0
        assert "afmag_airmt" in result.output


# ---------------------------------------------------------------------------
# pycsamt airborne diagnose
# ---------------------------------------------------------------------------


class TestAirborneDiagnose:
    def test_nonexistent_source_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["airborne", "diagnose", "/nonexistent/path"]
        )
        assert result.exit_code != 0

    def test_ztem_multiline_requires_line(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        """A multi-line ZTEM survey without --line must fail helpfully."""
        result = runner.invoke(main, ["airborne", "diagnose", str(ztem_xml_dir)])
        assert result.exit_code != 0
        assert "flight lines" in result.output
        assert "--line" in result.output

    def test_ztem_with_line_index(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["airborne", "diagnose", str(ztem_xml_dir), "--line", "0"]
        )
        assert result.exit_code == 0
        first_line = result.output.strip().splitlines()[0]
        assert "divergence" in first_line

    def test_ztem_bad_line_index_fails(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["airborne", "diagnose", str(ztem_xml_dir), "--line", "999"]
        )
        assert result.exit_code != 0

    def test_ztem_component_tzy(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            [
                "airborne", "diagnose", str(ztem_xml_dir),
                "--line", "1", "--component", "tzy",
            ],
        )
        assert result.exit_code == 0

    def test_ztem_json_output(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            [
                "airborne", "diagnose", str(ztem_xml_dir),
                "--line", "0", "--format", "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) > 0
        assert "divergence_real" in data[0]

    def test_mobilemt_admittance_table(
        self, runner: CliRunner, mobilemt_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["airborne", "diagnose", str(mobilemt_xml_dir)]
        )
        assert result.exit_code == 0
        assert "Yxx_real" in result.output

    def test_mobilemt_save_csv(
        self, runner: CliRunner, mobilemt_xml_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "admittance.csv"
        result = runner.invoke(
            main,
            [
                "airborne", "diagnose", str(mobilemt_xml_dir),
                "--output", str(out),
            ],
        )
        assert result.exit_code == 0
        assert out.exists()
        assert out.stat().st_size > 0

    def test_afmag_original_tilt_table(
        self, runner: CliRunner, afmag_original_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["airborne", "diagnose", str(afmag_original_xml_dir)]
        )
        assert result.exit_code == 0
        assert "tilt_deg" in result.output

    def test_afmag_airmt_tilt_table(
        self, runner: CliRunner, afmag_airmt_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["airborne", "diagnose", str(afmag_airmt_xml_dir)]
        )
        assert result.exit_code == 0
        assert "tilt_resultant_deg" in result.output

    def test_explicit_technology_override(
        self, runner: CliRunner, mobilemt_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            [
                "airborne", "diagnose", str(mobilemt_xml_dir),
                "--technology", "mobilemt",
            ],
        )
        assert result.exit_code == 0

    def test_invalid_technology_choice_fails(
        self, runner: CliRunner, mobilemt_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            [
                "airborne", "diagnose", str(mobilemt_xml_dir),
                "--technology", "not-a-technology",
            ],
        )
        assert result.exit_code != 0

    def test_verbose_flag(
        self, runner: CliRunner, mobilemt_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["airborne", "diagnose", str(mobilemt_xml_dir), "-v"]
        )
        assert result.exit_code == 0
