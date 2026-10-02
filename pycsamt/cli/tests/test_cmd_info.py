# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ``pycsamt info`` command."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
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

    def test_explicit_single_edi(self, runner: CliRunner, single_edi: Path) -> None:
        result = runner.invoke(main, ["info", str(single_edi)])
        assert result.exit_code == 0
        assert "Station" in result.output or "Frequencies" in result.output

    def test_explicit_edi_dir_text(self, runner: CliRunner, edi_dir: Path) -> None:
        result = runner.invoke(main, ["info", str(edi_dir)])
        assert result.exit_code == 0
        len(list(edi_dir.glob("*.edi")))
        # At least one station block in the output
        assert result.output.count("Station") >= 1 or result.output.count("File") >= 1

    def test_explicit_edi_dir_json(self, runner: CliRunner, edi_dir: Path) -> None:
        import json

        result = runner.invoke(main, ["info", str(edi_dir), "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) == len(list(edi_dir.glob("*.edi")))

    def test_explicit_edi_dir_csv(self, runner: CliRunner, edi_dir: Path) -> None:
        result = runner.invoke(main, ["info", str(edi_dir), "--format", "csv"])
        assert result.exit_code == 0
        lines = [l for l in result.output.splitlines() if l.strip()]
        assert len(lines) >= 2  # header + at least one data row

    # ------------------------------------------------------------------
    # Station filter
    # ------------------------------------------------------------------

    def test_station_filter_matching(self, runner: CliRunner, edi_dir: Path) -> None:
        edi_files = sorted(edi_dir.glob("*.edi"))
        first_stem = edi_files[0].stem.upper()
        result = runner.invoke(main, ["info", str(edi_dir), "--stations", first_stem])
        assert result.exit_code == 0

    def test_station_filter_no_match(self, runner: CliRunner, edi_dir: Path) -> None:
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
        assert result.exception is None or isinstance(result.exception, SystemExit)


# ---------------------------------------------------------------------------
# EMTF-XML transfer functions (real bundled airborne sample data)
# ---------------------------------------------------------------------------


class TestInfoEmtfXml:
    def test_single_xml_text_shows_format(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        xml_file = sorted(ztem_xml_dir.glob("*.xml"))[0]
        result = runner.invoke(main, ["info", str(xml_file)])
        assert result.exit_code == 0, result.output
        assert "Format" in result.output
        assert "xml" in result.output.lower()

    def test_single_xml_verbose_shows_covariance_and_orientation(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        xml_file = sorted(ztem_xml_dir.glob("*.xml"))[0]
        result = runner.invoke(main, ["info", str(xml_file), "-v"])
        assert result.exit_code == 0, result.output
        assert "Orientation" in result.output
        assert "Estimates" in result.output

    def test_xml_dir_json_has_orientation_and_tags(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["info", str(ztem_xml_dir), "--format", "json"]
        )
        assert result.exit_code == 0, result.output
        data = json.loads(result.output)
        assert len(data) == len(list(ztem_xml_dir.glob("*.xml")))
        first = data[0]
        assert first["format"] == "emtf_xml"
        assert first["orientation"] is not None
        assert "ztem" in first["tags"] or "tipper" in first["tags"]

    def test_xml_dir_csv_output(
        self, runner: CliRunner, ztem_xml_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["info", str(ztem_xml_dir), "--format", "csv"]
        )
        assert result.exit_code == 0, result.output
        lines = [l for l in result.output.splitlines() if l.strip()]
        assert "station" in lines[0]
        assert len(lines) >= 2


# ---------------------------------------------------------------------------
# _edi_info / _parse_edi_header — malformed-input branches
# ---------------------------------------------------------------------------


class TestEdiHeaderParsing:
    def test_bad_elev_and_freq_tokens_are_skipped(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        edi_text = (
            ">HEAD\n"
            '  DATAID="BADVALS"\n'
            "  LAT=30.0\n"
            "  LONG=100.0\n"
            "  ELEV=not_a_number\n"
            ">FREQ\n"
            "  not_a_freq also_not_a_freq\n"
        )
        p = tmp_path / "bad.edi"
        p.write_text(edi_text, encoding="utf-8")
        result = runner.invoke(main, ["info", str(p), "--format", "json"])
        assert result.exit_code == 0, result.output
        data = json.loads(result.output)[0]
        assert data["station"] == "BADVALS"
        assert data["elevation"] is None
        assert data["n_frequencies"] is None

        text_result = runner.invoke(main, ["info", str(p)])
        assert text_result.exit_code == 0, text_result.output
        assert "Frequencies: 0" in text_result.output

    def test_edi_info_records_error_when_header_parse_always_fails(
        self, runner: CliRunner, single_edi: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import pycsamt.cli.commands.info as info_mod

        def _boom(path, record):
            raise RuntimeError("simulated parse failure")

        monkeypatch.setattr(info_mod, "_parse_edi_header", _boom)
        result = runner.invoke(main, ["info", str(single_edi)])
        assert result.exit_code == 0, result.output
        assert "ERROR" in result.output
        assert "simulated parse failure" in result.output


# ---------------------------------------------------------------------------
# _collect_edi_paths — sites-based resolution branches
# ---------------------------------------------------------------------------


class TestCollectEdiPaths:
    def test_sites_with_missing_edi_attribute_are_skipped(self) -> None:
        from pycsamt.cli.commands.info import _collect_edi_paths

        class _NoEdiSite:
            pass

        paths = _collect_edi_paths(None, sites=[_NoEdiSite(), _NoEdiSite()])
        assert paths == []

    def test_sites_with_real_edi_path_are_collected(
        self, edi_dir: Path, single_edi: Path
    ) -> None:
        from pycsamt.cli.commands.info import _collect_edi_paths

        class _FakeEdi:
            def __init__(self, path):
                self.path = path

        class _FakeSite:
            def __init__(self, path):
                self.edi = _FakeEdi(path)

        other = sorted(edi_dir.glob("*.edi"))[-1]
        paths = _collect_edi_paths(
            None, sites=[_FakeSite(str(single_edi)), _FakeSite(str(other))]
        )
        assert set(paths) == {single_edi, other}

    def test_no_source_no_sites_returns_empty(self) -> None:
        from pycsamt.cli.commands.info import _collect_edi_paths

        assert _collect_edi_paths(None, sites=None) == []

    def test_single_edi_file_as_source(self, single_edi: Path) -> None:
        from pycsamt.cli.commands.info import _collect_edi_paths

        assert _collect_edi_paths(single_edi) == [single_edi]

    def test_non_edi_file_as_source_returns_empty(self, tmp_path: Path) -> None:
        from pycsamt.cli.commands.info import _collect_edi_paths

        p = tmp_path / "notes.txt"
        p.write_text("hi")
        assert _collect_edi_paths(p) == []


class TestSitesBasedResolution:
    def test_fresh_flag_falls_back_to_tf_paths_when_sites_yield_nothing(
        self, runner: CliRunner, edi_dir: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import pycsamt.cli.commands.info as info_mod

        class _NoEdiSite:
            pass

        fake_sites = [_NoEdiSite(), _NoEdiSite()]
        monkeypatch.setattr(
            info_mod, "resolve_survey", lambda *a, **k: fake_sites
        )
        result = runner.invoke(main, ["info", str(edi_dir), "--fresh"])
        assert result.exit_code == 0, result.output
        assert result.output.count("Station") >= 1

    def test_no_explicit_falls_back_to_context_survey_path(
        self, runner: CliRunner, edi_dir: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import pycsamt.cli.commands.info as info_mod

        class _NoEdiSite:
            pass

        class _FakeContext:
            survey_path = edi_dir

        monkeypatch.setattr(
            info_mod, "resolve_survey", lambda *a, **k: [_NoEdiSite()]
        )
        monkeypatch.setattr(
            info_mod.SurveyContext, "load", staticmethod(lambda: _FakeContext())
        )
        result = runner.invoke(main, ["info", "--fresh"])
        assert result.exit_code == 0, result.output
        assert result.output.count("Station") >= 1

    def test_no_explicit_no_active_context_reports_none_found(
        self,
        runner: CliRunner,
        edi_dir: Path,
        isolated_home: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        import pycsamt.cli.commands.info as info_mod

        class _NoEdiSite:
            pass

        monkeypatch.setattr(
            info_mod, "resolve_survey", lambda *a, **k: [_NoEdiSite()]
        )
        # isolated_home guarantees no ~/.pycsamt/context.json exists, so
        # SurveyContext.load() genuinely returns None here.
        result = runner.invoke(main, ["info", "--fresh"])
        assert result.exit_code != 0
        assert "No supported EDI/XML" in result.output

    def test_fresh_flag_real_sites_resolution_end_to_end(
        self, runner: CliRunner, edi_dir: Path
    ) -> None:
        result = runner.invoke(main, ["info", str(edi_dir), "--fresh"])
        assert result.exit_code == 0, result.output
        assert result.output.count("Station") >= 1


# ---------------------------------------------------------------------------
# Formatter helpers — direct unit coverage
# ---------------------------------------------------------------------------


class TestFormatters:
    def test_fmt_csv_empty_records_returns_empty_string(self) -> None:
        from pycsamt.cli.commands.info import _fmt_csv

        assert _fmt_csv([]) == ""
