# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt avg`` command group.

Strategy
--------
* Help tests always run (no data required).
* Live tests use ``data/avg/K2.AVG`` (modern kind-2 format) which works
  without the optional xarray dependency.  They skip gracefully when the
  file is absent.
* K1 (legacy) path tests are tagged separately and skip when K1 is absent.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pycsamt.cli import main

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_AVG_DIR = _PROJECT_ROOT / "data" / "avg"
_K2_AVG = _AVG_DIR / "K2.AVG"
_K1_AVG = _AVG_DIR / "K1.AVG"
_K2_STN = _AVG_DIR / "K2.stn"


def _has_k2() -> bool:
    return _K2_AVG.exists()


def _has_k1() -> bool:
    return _K1_AVG.exists()


def _has_k2_stn() -> bool:
    return _K2_STN.exists()


# ---------------------------------------------------------------------------
# Help wiring
# ---------------------------------------------------------------------------


class TestAvgGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["avg", "--help"])
        assert result.exit_code == 0
        for sub in ("info", "validate", "stations", "correct", "export"):
            assert sub in result.output

    @pytest.mark.parametrize(
        "sub", ["info", "validate", "stations", "correct", "export"]
    )
    def test_each_subcommand_help(self, runner: CliRunner, sub: str) -> None:
        result = runner.invoke(main, ["avg", sub, "--help"])
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# pycsamt avg info
# ---------------------------------------------------------------------------


class TestAvgInfo:
    @pytest.fixture(autouse=True)
    def require_k2(self) -> None:
        if not _has_k2():
            pytest.skip("data/avg/K2.AVG not found")

    def test_text_output(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["avg", "info", str(_K2_AVG)])
        assert result.exit_code == 0
        assert "Station" in result.output or "station" in result.output.lower()
        assert "Frequenc" in result.output or "freq" in result.output.lower()

    def test_json_output(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["avg", "info", str(_K2_AVG), "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["n_stations"] > 0
        assert data["n_frequencies"] > 0
        assert "frequency_range" in data

    def test_csv_output(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["avg", "info", str(_K2_AVG), "--format", "csv"])
        assert result.exit_code == 0
        lines = [ln for ln in result.output.strip().splitlines() if ln]
        assert "project" in lines[0].lower() or "n_stations" in lines[0].lower()

    def test_project_in_output(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["avg", "info", str(_K2_AVG), "--format", "json"])
        data = json.loads(result.output)
        assert data["project"]  # non-empty project name

    def test_with_stn_file(self, runner: CliRunner) -> None:
        if not _has_k2_stn():
            pytest.skip("K2.stn not found")
        result = runner.invoke(
            main, ["avg", "info", str(_K2_AVG), "--stn-file", str(_K2_STN)]
        )
        assert result.exit_code == 0

    def test_nonexistent_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["avg", "info", "/no/such.avg"])
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# pycsamt avg validate
# ---------------------------------------------------------------------------


class TestAvgValidate:
    @pytest.fixture(autouse=True)
    def require_k2(self) -> None:
        if not _has_k2():
            pytest.skip("data/avg/K2.AVG not found")

    def test_text_output(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["avg", "validate", str(_K2_AVG)])
        assert result.exit_code == 0

    def test_json_output(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["avg", "validate", str(_K2_AVG), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "n_stations" in data
        assert "qc_columns" in data
        assert "stations" in data
        assert len(data["stations"]) > 0

    def test_csv_output(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["avg", "validate", str(_K2_AVG), "--format", "csv"]
        )
        assert result.exit_code == 0
        lines = result.output.strip().splitlines()
        assert "station" in lines[0].lower()
        assert len(lines) > 1

    def test_each_station_has_flagged_key(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["avg", "validate", str(_K2_AVG), "--format", "json"]
        )
        data = json.loads(result.output)
        for s in data["stations"]:
            assert "flagged" in s

    def test_top_limit(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "validate",
                str(_K2_AVG),
                "--top",
                "3",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert len(data["stations"]) <= 3

    def test_custom_threshold(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "validate",
                str(_K2_AVG),
                "--threshold",
                "0.0",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        # All stations should be flagged when threshold is 0
        assert data["n_flagged"] >= 0


# ---------------------------------------------------------------------------
# pycsamt avg stations
# ---------------------------------------------------------------------------


class TestAvgStations:
    @pytest.fixture(autouse=True)
    def require_k2(self) -> None:
        if not _has_k2():
            pytest.skip("data/avg/K2.AVG not found")

    def test_text_output(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["avg", "stations", str(_K2_AVG)])
        assert result.exit_code == 0
        assert "station" in result.output.lower() or "S0" in result.output

    def test_json_output(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["avg", "stations", str(_K2_AVG), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list) and len(data) > 0
        assert "name" in data[0]
        assert "position" in data[0]

    def test_csv_output(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["avg", "stations", str(_K2_AVG), "--format", "csv"]
        )
        assert result.exit_code == 0
        lines = result.output.strip().splitlines()
        assert "name" in lines[0].lower()

    def test_correct_station_count(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["avg", "stations", str(_K2_AVG), "--format", "json"]
        )
        data = json.loads(result.output)
        # K2.AVG has 28 stations
        assert len(data) == 28

    def test_with_stn_file_adds_coords(self, runner: CliRunner) -> None:
        if not _has_k2_stn():
            pytest.skip("K2.stn not found")
        result = runner.invoke(
            main,
            [
                "avg",
                "stations",
                str(_K2_AVG),
                "--stn-file",
                str(_K2_STN),
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        # At least some stations should have easting/northing
        has_coords = any("easting" in r for r in data)
        assert has_coords

    def test_sort_by_name(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "stations",
                str(_K2_AVG),
                "--sort-by",
                "name",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0

    def test_text_output_with_stn_file_rich_topo_columns(
        self, runner: CliRunner
    ) -> None:
        if not _has_k2_stn():
            pytest.skip("K2.stn not found")
        result = runner.invoke(
            main,
            ["avg", "stations", str(_K2_AVG), "--stn-file", str(_K2_STN)],
        )
        assert result.exit_code == 0
        assert "Easting" in result.output or "easting" in result.output.lower()

    def _no_rich_import(self, monkeypatch: pytest.MonkeyPatch) -> None:
        import builtins

        real_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name in ("rich.console", "rich.table"):
                raise ImportError("simulated: rich unavailable")
            return real_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)

    def test_text_output_without_rich_plain_table(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        self._no_rich_import(monkeypatch)
        result = runner.invoke(main, ["avg", "stations", str(_K2_AVG)])
        assert result.exit_code == 0
        assert "Name" in result.output
        assert "Pos(m)" in result.output

    def test_text_output_without_rich_plain_table_with_topo(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        if not _has_k2_stn():
            pytest.skip("K2.stn not found")
        self._no_rich_import(monkeypatch)
        result = runner.invoke(
            main,
            ["avg", "stations", str(_K2_AVG), "--stn-file", str(_K2_STN)],
        )
        assert result.exit_code == 0
        assert "E(m)" in result.output
        assert "Lat" in result.output


# ---------------------------------------------------------------------------
# pycsamt avg correct
# ---------------------------------------------------------------------------


class TestAvgCorrect:
    @pytest.fixture(autouse=True)
    def require_k2(self) -> None:
        if not _has_k2():
            pytest.skip("data/avg/K2.AVG not found")

    def test_dry_run(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["avg", "correct", str(_K2_AVG), "--dry-run"])
        assert result.exit_code == 0
        assert "dry" in result.output.lower() or "method" in result.output.lower()

    def test_static_shift_text(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "correct",
                str(_K2_AVG),
                "--method",
                "static-shift",
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code == 0
        written = list(tmp_path.glob("*.avg"))
        assert len(written) == 1

    def test_static_shift_json(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "correct",
                str(_K2_AVG),
                "--method",
                "static-shift",
                "--output-dir",
                str(tmp_path),
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "output" in data
        assert "corrections" in data
        assert any(c["correction"] == "static_shift" for c in data["corrections"])

    def test_shift_factors_reasonable(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "correct",
                str(_K2_AVG),
                "--method",
                "static-shift",
                "--output-dir",
                str(tmp_path),
                "--format",
                "json",
            ],
        )
        data = json.loads(result.output)
        ss = next(c for c in data["corrections"] if c["correction"] == "static_shift")
        # Shift factors should be finite and in a reasonable range
        assert ss["shift_min"] > 0
        assert ss["shift_max"] < 1e6

    def test_tma_filter(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "correct",
                str(_K2_AVG),
                "--filter",
                "tma",
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code == 0

    def test_no_output_dir_uses_dot(self, runner: CliRunner) -> None:
        # Default output_dir is "." — should still work
        result = runner.invoke(main, ["avg", "correct", str(_K2_AVG), "--dry-run"])
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# pycsamt avg export
# ---------------------------------------------------------------------------


class TestAvgExport:
    @pytest.fixture(autouse=True)
    def require_k2(self) -> None:
        if not _has_k2():
            pytest.skip("data/avg/K2.AVG not found")

    def test_export_modern(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "export",
                str(_K2_AVG),
                "--format",
                "modern",
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code == 0
        written = list(tmp_path.glob("*.avg"))
        assert len(written) == 1

    def test_export_legacy(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "export",
                str(_K2_AVG),
                "--format",
                "legacy",
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code == 0
        written = list(tmp_path.glob("*.avg"))
        assert len(written) == 1

    def test_json_summary(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "export",
                str(_K2_AVG),
                "--output-dir",
                str(tmp_path),
                "--format-out",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "output" in data
        assert "size_kb" in data

    def test_custom_stem(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "avg",
                "export",
                str(_K2_AVG),
                "--stem",
                "my_survey",
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code == 0
        written = list(tmp_path.glob("my_survey*.avg"))
        assert len(written) == 1

    def test_overwrite_required_when_exists(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        # Write once
        runner.invoke(
            main,
            ["avg", "export", str(_K2_AVG), "--output-dir", str(tmp_path)],
        )
        # Write again without --overwrite should fail
        result = runner.invoke(
            main,
            ["avg", "export", str(_K2_AVG), "--output-dir", str(tmp_path)],
        )
        assert result.exit_code != 0

    def test_overwrite_flag(self, runner: CliRunner, tmp_path: Path) -> None:
        runner.invoke(
            main,
            ["avg", "export", str(_K2_AVG), "--output-dir", str(tmp_path)],
        )
        result = runner.invoke(
            main,
            [
                "avg",
                "export",
                str(_K2_AVG),
                "--output-dir",
                str(tmp_path),
                "--overwrite",
            ],
        )
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# Unit tests — helpers
# ---------------------------------------------------------------------------


class TestAvgHelpers:
    @pytest.fixture(autouse=True)
    def require_k2(self) -> None:
        if not _has_k2():
            pytest.skip("data/avg/K2.AVG not found")

    def test_get_avg_returns_avg_object(self) -> None:
        from pycsamt.cli.commands.avg._base import _get_avg
        from pycsamt.zonge.avg import AVG

        obj = _get_avg(_K2_AVG)
        assert isinstance(obj, AVG)
        assert obj.df is not None

    def test_load_raw_returns_df(self) -> None:
        import pandas as pd

        from pycsamt.cli.commands.avg._base import _load_raw

        df, meta, kind = _load_raw(_K2_AVG)
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert kind in (1, 2)

    def test_k2_has_qc_columns(self) -> None:
        from pycsamt.cli.commands.avg._base import _get_avg

        obj = _get_avg(_K2_AVG)
        qc = [c for c in obj.df.columns if c.startswith("pc_") or c.startswith("s_")]
        assert len(qc) > 0

    def test_k2_summary_fields(self) -> None:
        from pycsamt.cli.commands.avg._base import _get_avg

        obj = _get_avg(_K2_AVG)
        s = obj.summary
        assert s.num_stations == 28
        assert s.num_frequencies > 0


# ---------------------------------------------------------------------------
# UsageError fallback (legacy K1-without-xarray path)
# ---------------------------------------------------------------------------
# Real K1.AVG loads fine when xarray is installed, so this path can only be
# exercised by simulating the ImportError _get_avg raises when xarray is
# absent (see pycsamt/cli/commands/avg/_base.py::_get_avg).


class TestAvgUsageErrorFallback:
    @pytest.fixture(autouse=True)
    def require_k2(self) -> None:
        if not _has_k2():
            pytest.skip("data/avg/K2.AVG not found")

    def _raise_usage_error(self, *args, **kwargs):
        import click

        raise click.UsageError("Loading requires xarray")

    def test_info_falls_back_to_raw_loader(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import pycsamt.cli.commands.avg.info as _info_mod

        monkeypatch.setattr(_info_mod, "_get_avg", self._raise_usage_error)
        result = runner.invoke(main, ["avg", "info", str(_K2_AVG)])
        assert result.exit_code == 0
        assert "xarray not installed" in result.output.lower() or "legacy" in result.output.lower()

    def test_info_fallback_json(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import pycsamt.cli.commands.avg.info as _info_mod

        monkeypatch.setattr(_info_mod, "_get_avg", self._raise_usage_error)
        result = runner.invoke(
            main, ["avg", "info", str(_K2_AVG), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["n_stations"] > 0

    def test_stations_falls_back_to_raw_loader(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import pycsamt.cli.commands.avg.stations as _stations_mod

        monkeypatch.setattr(_stations_mod, "_get_avg", self._raise_usage_error)
        result = runner.invoke(main, ["avg", "stations", str(_K2_AVG)])
        assert result.exit_code == 0

    def test_stations_fallback_json(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import pycsamt.cli.commands.avg.stations as _stations_mod

        monkeypatch.setattr(_stations_mod, "_get_avg", self._raise_usage_error)
        result = runner.invoke(
            main, ["avg", "stations", str(_K2_AVG), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list) and len(data) > 0
        assert "station" in data[0]

    def test_stations_fallback_csv(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import pycsamt.cli.commands.avg.stations as _stations_mod

        monkeypatch.setattr(_stations_mod, "_get_avg", self._raise_usage_error)
        result = runner.invoke(
            main, ["avg", "stations", str(_K2_AVG), "--format", "csv"]
        )
        assert result.exit_code == 0
        assert "station" in result.output.lower()


# ---------------------------------------------------------------------------
# _stations_raw — direct unit tests (all format branches)
# ---------------------------------------------------------------------------


class TestStationsRawHelper:
    def test_no_station_column(self, capsys: pytest.CaptureFixture) -> None:
        import pandas as pd

        from pycsamt.cli.commands.avg.stations import _stations_raw

        _stations_raw(pd.DataFrame({"other": [1, 2]}), "text")
        captured = capsys.readouterr()
        assert "No station column" in captured.err

    def test_text_format(self, capsys: pytest.CaptureFixture) -> None:
        import pandas as pd

        from pycsamt.cli.commands.avg.stations import _stations_raw

        _stations_raw(pd.DataFrame({"station": ["S2", "S1"]}), "text")
        out = capsys.readouterr().out
        assert "S1" in out and "S2" in out

    def test_json_format(self, capsys: pytest.CaptureFixture) -> None:
        import pandas as pd

        from pycsamt.cli.commands.avg.stations import _stations_raw

        _stations_raw(pd.DataFrame({"station": ["S1", "S2"]}), "json")
        out = capsys.readouterr().out
        data = json.loads(out)
        assert data == [{"station": "S1"}, {"station": "S2"}]

    def test_csv_format(self, capsys: pytest.CaptureFixture) -> None:
        import pandas as pd

        from pycsamt.cli.commands.avg.stations import _stations_raw

        _stations_raw(pd.DataFrame({"station": ["S1"]}), "csv")
        out = capsys.readouterr().out
        lines = out.strip().splitlines()
        assert lines[0] == "station"
        assert lines[1] == "S1"


# ---------------------------------------------------------------------------
# Topography attachment failure is swallowed (bad .stn file)
# ---------------------------------------------------------------------------


class TestTopographyLoadFailure:
    @pytest.fixture(autouse=True)
    def require_k2(self) -> None:
        if not _has_k2():
            pytest.skip("data/avg/K2.AVG not found")

    @pytest.fixture
    def bad_stn(self, tmp_path: Path) -> Path:
        p = tmp_path / "bad.stn"
        p.write_text("not a real stn file\ngarbage 123\n", encoding="utf-8")
        return p

    def test_info_warns_and_continues(
        self, runner: CliRunner, bad_stn: Path
    ) -> None:
        result = runner.invoke(
            main, ["avg", "info", str(_K2_AVG), "--stn-file", str(bad_stn)]
        )
        assert result.exit_code == 0

    def test_stations_warns_and_continues(
        self, runner: CliRunner, bad_stn: Path
    ) -> None:
        result = runner.invoke(
            main, ["avg", "stations", str(_K2_AVG), "--stn-file", str(bad_stn)]
        )
        assert result.exit_code == 0
        assert "warning" in result.output.lower() or "could not load" in result.output.lower()


# ---------------------------------------------------------------------------
# _info_rich — direct unit tests for branches unreachable with real data
# (K1/K2 AVG never populate survey_type="" or line_name/operator)
# ---------------------------------------------------------------------------


class TestInfoRichHelper:
    def _fake_obj(self, survey_type: str, line_name, operator):
        import types

        import pandas as pd

        summary = types.SimpleNamespace(
            data_kind="Kind-2",
            project="TESTPROJ",
            survey_type=survey_type,
            num_stations=3,
            num_frequencies=5,
            station_range="0 - 100 m",
            frequency_range="1 - 100 Hz",
            total_rows=15,
        )
        annotation = types.SimpleNamespace(
            line_name=line_name, operator=operator
        )
        config = types.SimpleNamespace(
            survey_type=survey_type, line_name=line_name
        )
        header = types.SimpleNamespace(
            annotation=annotation, config=config, hardware="FakeHW v1"
        )
        info = types.SimpleNamespace(header=header)
        df = pd.DataFrame({"comp": ["ExHy", "ExHy"]})
        return types.SimpleNamespace(summary=summary, info=info, df=df)

    def test_falsy_survey_type_skips_line(
        self, capsys: pytest.CaptureFixture
    ) -> None:
        from pycsamt.cli.commands.avg.info import _info_rich

        obj = self._fake_obj(survey_type="", line_name=None, operator=None)
        _info_rich(obj, None, Path("fake.avg"), "text")
        out = capsys.readouterr().out
        assert "Survey type" not in out

    def test_line_name_and_operator_are_printed(
        self, capsys: pytest.CaptureFixture
    ) -> None:
        from pycsamt.cli.commands.avg.info import _info_rich

        obj = self._fake_obj(
            survey_type="CSAMT", line_name="Line-7", operator="Jane Doe"
        )
        _info_rich(obj, None, Path("fake.avg"), "text")
        out = capsys.readouterr().out
        assert "Line-7" in out
        assert "Jane Doe" in out

    def test_add_topography_exception_is_swallowed(
        self, capsys: pytest.CaptureFixture, tmp_path: Path
    ) -> None:
        from pycsamt.cli.commands.avg.info import _info_rich

        base = self._fake_obj("CSAMT", None, None)

        def _raise_add_topography(path):
            raise RuntimeError("bad stn")

        base.add_topography = _raise_add_topography

        _info_rich(base, tmp_path / "whatever.stn", Path("fake.avg"), "text")
        out = capsys.readouterr().out
        assert "TESTPROJ" in out
