# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt jones`` command group.

Strategy
--------
* Help tests always run (no data required).
* All live tests use synthetic J-files created in tmp_path so they
  run without repository data.  Tests that require the optional
  ``data/j/`` or ``data/jc/`` directories skip gracefully when absent.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pycsamt.cli import main

# ---------------------------------------------------------------------------
# Minimal valid J-file content
# ---------------------------------------------------------------------------

_J_RXY_CONTENT = """\
#WRITTEN BY PYCSAMT-TEST: S01 2024-01-01
>LATITUDE=30.0
>LONGITUDE=100.0
>ELEVATION=500.0
>AZIMUTH=0.0
S01
RXY
3
1.0 100.0 45.0 110.0 90.0 50.0 40.0 1.0 1.0
2.0 200.0 50.0 220.0 180.0 60.0 40.0 1.0 1.0
5.0 300.0 55.0 330.0 270.0 65.0 45.0 1.0 1.0
S01
RYX
3
1.0 120.0 35.0 132.0 108.0 45.0 25.0 1.0 1.0
2.0 220.0 40.0 242.0 198.0 50.0 30.0 1.0 1.0
5.0 320.0 45.0 352.0 288.0 55.0 35.0 1.0 1.0
"""

_J_S02_CONTENT = """\
#WRITTEN BY PYCSAMT-TEST: S02 2024-01-01
>LATITUDE=31.0
>LONGITUDE=101.0
>ELEVATION=600.0
>AZIMUTH=10.0
S02
RXY
2
1.0 80.0 40.0 88.0 72.0 48.0 32.0 1.0 1.0
2.0 160.0 45.0 176.0 144.0 55.0 35.0 1.0 1.0
S02
RYX
2
1.0 90.0 30.0 99.0 81.0 40.0 20.0 1.0 1.0
2.0 180.0 35.0 198.0 162.0 45.0 25.0 1.0 1.0
"""


@pytest.fixture(scope="module")
def j_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A temporary directory with two valid J-files."""
    d = tmp_path_factory.mktemp("jones_data")
    (d / "S01.j").write_text(_J_RXY_CONTENT)
    (d / "S02.j").write_text(_J_S02_CONTENT)
    return d


@pytest.fixture(scope="module")
def j_single(j_dir: Path) -> Path:
    """Single J-file path."""
    return j_dir / "S01.j"


# ---------------------------------------------------------------------------
# Data paths (optional real data)
# ---------------------------------------------------------------------------

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_J_DATA       = _PROJECT_ROOT / "data" / "j"
_JC_DATA      = _PROJECT_ROOT / "data" / "jc"


def _has_j_data() -> bool:
    return _J_DATA.exists() and bool(list(_J_DATA.glob("*")))


def _has_jc_data() -> bool:
    return _JC_DATA.exists() and bool(list(_JC_DATA.glob("*")))


# ---------------------------------------------------------------------------
# Help wiring
# ---------------------------------------------------------------------------

class TestJonesGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["jones", "--help"])
        assert result.exit_code == 0
        for sub in ("info", "validate", "stations", "blocks", "select"):
            assert sub in result.output

    @pytest.mark.parametrize("sub", [
        "info", "validate", "stations", "blocks", "select"
    ])
    def test_each_subcommand_help(self, runner: CliRunner, sub: str) -> None:
        result = runner.invoke(main, ["jones", sub, "--help"])
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# pycsamt jones info
# ---------------------------------------------------------------------------

class TestJonesInfo:
    def test_single_file_text(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(main, ["jones", "info", str(j_single)])
        assert result.exit_code == 0
        assert "S01" in result.output
        assert "n_freq" in result.output

    def test_single_file_json(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(
            main, ["jones", "info", str(j_single), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["station"] == "S01"
        assert data["lat"] == pytest.approx(30.0)
        assert data["n_freq"] == 3

    def test_single_file_csv(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(
            main, ["jones", "info", str(j_single), "--format", "csv"]
        )
        assert result.exit_code == 0
        lines = [l for l in result.output.strip().splitlines() if l]
        assert "station" in lines[0]
        assert len(lines) >= 2

    def test_directory_text(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(main, ["jones", "info", str(j_dir)])
        assert result.exit_code == 0

    def test_directory_json(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(
            main, ["jones", "info", str(j_dir), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) == 2
        stations = {r["station"] for r in data}
        assert "S01" in stations and "S02" in stations

    def test_top_limit(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(
            main, ["jones", "info", str(j_dir), "--top", "1", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert len(data) == 1


# ---------------------------------------------------------------------------
# pycsamt jones validate
# ---------------------------------------------------------------------------

class TestJonesValidate:
    def test_single_valid(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(
            main, ["jones", "validate", str(j_single)]
        )
        assert result.exit_code == 0
        assert "ok" in result.output.lower() or "valid" in result.output.lower()

    def test_directory_valid(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(main, ["jones", "validate", str(j_dir)])
        assert result.exit_code == 0

    def test_json_output(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(
            main, ["jones", "validate", str(j_dir), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "n_ok" in data and data["n_ok"] == 2
        assert data["n_fail"] == 0

    def test_no_deep(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(
            main, ["jones", "validate", str(j_single), "--no-deep"]
        )
        assert result.exit_code == 0

    def test_invalid_file_exits_1(self, runner: CliRunner, tmp_path: Path) -> None:
        bad = tmp_path / "bad.j"
        bad.write_text("this is not a j file\n")
        result = runner.invoke(main, ["jones", "validate", str(bad)])
        assert result.exit_code != 0

    def test_nonexistent_exits(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["jones", "validate", "/no/such/path"])
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# pycsamt jones stations
# ---------------------------------------------------------------------------

class TestJonesStations:
    def test_text_output(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(main, ["jones", "stations", str(j_dir)])
        assert result.exit_code == 0
        assert "2 station" in result.output

    def test_json_output(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(
            main, ["jones", "stations", str(j_dir), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert len(data) == 2
        assert all("station" in r and "lat" in r for r in data)

    def test_csv_output(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(
            main, ["jones", "stations", str(j_dir), "--format", "csv"]
        )
        assert result.exit_code == 0
        lines = result.output.strip().splitlines()
        assert "station" in lines[0].lower()

    def test_sort_by_lat(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(
            main, ["jones", "stations", str(j_dir),
                   "--sort-by", "lat", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        lats = [r["lat"] for r in data if r.get("lat") is not None]
        if len(lats) > 1:
            assert lats == sorted(lats)

    def test_top_limit(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(
            main, ["jones", "stations", str(j_dir), "--top", "1", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert len(data) == 1


# ---------------------------------------------------------------------------
# pycsamt jones blocks
# ---------------------------------------------------------------------------

class TestJonesBlocks:
    def test_text_output(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(main, ["jones", "blocks", str(j_single)])
        assert result.exit_code == 0
        assert "S01" in result.output
        assert "RXY" in result.output or "RYX" in result.output

    def test_json_output(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(
            main, ["jones", "blocks", str(j_single), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["n_blocks"] == 2
        assert data["station"] == "S01"
        assert all("token" in b and "nrows" in b for b in data["blocks"])

    def test_csv_output(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(
            main, ["jones", "blocks", str(j_single), "--format", "csv"]
        )
        assert result.exit_code == 0
        lines = result.output.strip().splitlines()
        assert "token" in lines[0].lower()

    def test_kind_filter(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(
            main, ["jones", "blocks", str(j_single), "--kind", "R", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert all(b["kind"] == "R" for b in data["blocks"])

    def test_comp_filter(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(
            main, ["jones", "blocks", str(j_single),
                   "--comp", "XY", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert all(b["comp"] == "XY" for b in data["blocks"])

    def test_qa_flag(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(
            main, ["jones", "blocks", str(j_single), "--qa"]
        )
        assert result.exit_code == 0

    def test_no_match_shows_error(self, runner: CliRunner, j_single: Path) -> None:
        result = runner.invoke(
            main, ["jones", "blocks", str(j_single), "--kind", "T"]
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# pycsamt jones select
# ---------------------------------------------------------------------------

class TestJonesSelect:
    def test_dry_run_all(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(
            main, ["jones", "select", str(j_dir), "--dry-run"]
        )
        assert result.exit_code == 0
        assert "2" in result.output  # 2 stations match

    def test_dry_run_json(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(
            main, ["jones", "select", str(j_dir), "--dry-run", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["n_selected"] == 2

    def test_select_by_station(
        self, runner: CliRunner, j_dir: Path, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            main, ["jones", "select", str(j_dir),
                   "--stations", "S01", "--output-dir", str(tmp_path)]
        )
        assert result.exit_code == 0
        written = list(tmp_path.glob("*.j"))
        assert len(written) >= 1

    def test_select_all_exports(
        self, runner: CliRunner, j_dir: Path, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            main, ["jones", "select", str(j_dir), "--output-dir", str(tmp_path)]
        )
        assert result.exit_code == 0
        written = list(tmp_path.glob("*.j"))
        assert len(written) == 2

    def test_select_json_output(
        self, runner: CliRunner, j_dir: Path, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            main, ["jones", "select", str(j_dir),
                   "--output-dir", str(tmp_path), "--format", "json"]
        )
        assert result.exit_code == 0
        # tqdm may prefix the JSON — find the first '{'
        raw = result.output
        json_start = raw.find("{")
        assert json_start >= 0, f"No JSON in output: {raw[:200]!r}"
        data = json.loads(raw[json_start:])
        assert "n_written" in data
        assert data["n_written"] == 2

    def test_has_r_filter(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(
            main, ["jones", "select", str(j_dir), "--has-r", "--dry-run",
                   "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["n_selected"] == 2  # both have R data

    def test_has_z_filter_dry_run(self, runner: CliRunner, j_dir: Path) -> None:
        # Our synthetic files have Res (converted from R) but may or may not
        # have explicit Z. Just verify the filter runs without crashing.
        result = runner.invoke(
            main, ["jones", "select", str(j_dir), "--has-z", "--dry-run",
                   "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "n_selected" in data

    def test_no_output_dir_fails(self, runner: CliRunner, j_dir: Path) -> None:
        result = runner.invoke(main, ["jones", "select", str(j_dir)])
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# Unit tests — low-level helpers
# ---------------------------------------------------------------------------

class TestJonesHelpers:
    def test_get_jfile(self, j_single: Path) -> None:
        from pycsamt.cli.commands.jones._base import _get_jfile
        from pycsamt.jones.j import JFile
        jf = _get_jfile(j_single)
        assert isinstance(jf, JFile)
        assert jf.station == "S01"
        assert jf.n_freq == 3

    def test_get_collection(self, j_dir: Path) -> None:
        from pycsamt.cli.commands.jones._base import _get_collection
        from pycsamt.jones.collection import JCollection
        coll = _get_collection(j_dir)
        assert isinstance(coll, JCollection)
        assert len(coll) == 2

    def test_jfile_has_res(self, j_single: Path) -> None:
        from pycsamt.cli.commands.jones._base import _get_jfile
        jf = _get_jfile(j_single)
        assert jf.Res is not None

    def test_jfile_coordinates(self, j_single: Path) -> None:
        from pycsamt.cli.commands.jones._base import _get_jfile
        jf = _get_jfile(j_single)
        assert jf.lat == pytest.approx(30.0)
        assert jf.lon == pytest.approx(100.0)
        assert jf.elev == pytest.approx(500.0)

    def test_blocks_qa_summary(self, j_single: Path) -> None:
        from pycsamt.cli.commands.jones._base import _get_jfile
        jf = _get_jfile(j_single)
        blks = jf.blocks
        assert blks is not None and blks.n == 2
        for b in blks.blocks:
            qa = b.qa_summary()
            assert isinstance(qa, dict)
