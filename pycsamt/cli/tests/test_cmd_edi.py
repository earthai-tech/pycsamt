# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt edi`` command group.

Strategy
--------
* Help tests always run (no data required).
* Live tests use data/ directories that already exist in the repo;
  each test skips gracefully when the directory is absent.
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
_EDI_3EDIS    = _PROJECT_ROOT / "data" / "3edis"
_EDI_TIPPER   = _PROJECT_ROOT / "data" / "AMT" / "TIPPER"
_EDI_WILLY    = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_EDI_AMT      = _PROJECT_ROOT / "data" / "AMT"


def _has(p: Path) -> bool:
    return p.exists() and bool(list(p.glob("*.edi")))


# ---------------------------------------------------------------------------
# Help wiring
# ---------------------------------------------------------------------------

class TestEdiGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["edi", "--help"])
        assert result.exit_code == 0
        for sub in ("info", "validate", "stations", "profile", "rotate", "select"):
            assert sub in result.output

    @pytest.mark.parametrize("sub", [
        "info", "validate", "stations", "profile", "rotate", "select"
    ])
    def test_each_subcommand_help(self, runner: CliRunner, sub: str) -> None:
        result = runner.invoke(main, ["edi", sub, "--help"])
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# pycsamt edi info
# ---------------------------------------------------------------------------

class TestEdiInfo:
    @pytest.fixture(autouse=True)
    def require_data(self) -> None:
        if not _has(_EDI_3EDIS):
            pytest.skip("data/3edis/ not found")

    def test_single_file_text(self, runner: CliRunner) -> None:
        edi_path = sorted(_EDI_3EDIS.glob("*.edi"))[0]
        result = runner.invoke(main, ["edi", "info", str(edi_path)])
        assert result.exit_code == 0
        assert "Station" in result.output or "station" in result.output.lower()
        assert "n_freq" in result.output

    def test_single_file_json(self, runner: CliRunner) -> None:
        edi_path = sorted(_EDI_3EDIS.glob("*.edi"))[0]
        result = runner.invoke(
            main, ["edi", "info", str(edi_path), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "n_freq" in data
        assert "has_tipper" in data

    def test_directory_text(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["edi", "info", str(_EDI_3EDIS)])
        assert result.exit_code == 0

    def test_directory_json(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "info", str(_EDI_3EDIS), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) > 0
        assert "n_freq" in data[0]

    def test_directory_csv(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "info", str(_EDI_3EDIS), "--format", "csv"]
        )
        assert result.exit_code == 0
        lines = [ln for ln in result.output.strip().splitlines() if ln]
        assert "station" in lines[0]
        assert len(lines) >= 2

    def test_top_limit(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "info", str(_EDI_3EDIS), "--top", "1", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert len(data) == 1


# ---------------------------------------------------------------------------
# pycsamt edi validate
# ---------------------------------------------------------------------------

class TestEdiValidate:
    @pytest.fixture(autouse=True)
    def require_data(self) -> None:
        if not _has(_EDI_3EDIS):
            pytest.skip("data/3edis/ not found")

    def test_directory_all_valid(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["edi", "validate", str(_EDI_3EDIS)])
        assert result.exit_code == 0
        assert "ok" in result.output or "valid" in result.output.lower()

    def test_json_output(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "validate", str(_EDI_3EDIS), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "n_ok" in data and "n_fail" in data
        assert data["n_ok"] > 0

    def test_single_file_valid(self, runner: CliRunner) -> None:
        edi_path = sorted(_EDI_3EDIS.glob("*.edi"))[0]
        result = runner.invoke(main, ["edi", "validate", str(edi_path)])
        assert result.exit_code == 0

    def test_no_deep_mode(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "validate", str(_EDI_3EDIS), "--no-deep"]
        )
        assert result.exit_code == 0

    def test_nonexistent_exits_nonzero(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["edi", "validate", "/no/such/path"])
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# pycsamt edi stations
# ---------------------------------------------------------------------------

class TestEdiStations:
    @pytest.fixture(autouse=True)
    def require_data(self) -> None:
        if not (_has(_EDI_WILLY) or _has(_EDI_3EDIS)):
            pytest.skip("no EDI data directory found")
        self.edi_dir = _EDI_WILLY if _has(_EDI_WILLY) else _EDI_3EDIS

    def test_text_output(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["edi", "stations", str(self.edi_dir)])
        assert result.exit_code == 0

    def test_json_output(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "stations", str(self.edi_dir), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list) and len(data) > 0
        assert "station" in data[0]

    def test_csv_output(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "stations", str(self.edi_dir), "--format", "csv"]
        )
        assert result.exit_code == 0
        lines = result.output.strip().splitlines()
        assert "station" in lines[0].lower()

    def test_top_limit(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "stations", str(self.edi_dir),
                   "--top", "2", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert len(data) <= 2

    def test_sort_by_lat(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "stations", str(self.edi_dir),
                   "--sort-by", "lat", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        lats = [r.get("lat") for r in data if r.get("lat") is not None]
        if len(lats) > 1:
            assert lats == sorted(lats)


# ---------------------------------------------------------------------------
# pycsamt edi profile
# ---------------------------------------------------------------------------

class TestEdiProfile:
    @pytest.fixture(autouse=True)
    def require_data(self) -> None:
        if not _has(_EDI_WILLY):
            pytest.skip("data/AMT/WILLY_DATA/L18PLT/ not found")

    def test_text_output(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["edi", "profile", str(_EDI_WILLY)])
        assert result.exit_code == 0
        assert "Bearing" in result.output or "bearing" in result.output.lower()

    def test_json_output(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "profile", str(_EDI_WILLY), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "bearing_deg" in data and "step_m" in data
        assert data["n_stations"] > 0

    def test_bearing_method_linear(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "profile", str(_EDI_WILLY),
                   "--bearing-method", "linear", "--format", "json"]
        )
        assert result.exit_code == 0
        assert json.loads(result.output).get("bearing_method") == "linear"

    def test_distances_flag(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "profile", str(_EDI_WILLY),
                   "--distances", "--format", "json"]
        )
        assert result.exit_code == 0
        assert "station_table" in json.loads(result.output)


# ---------------------------------------------------------------------------
# pycsamt edi rotate
# ---------------------------------------------------------------------------

class TestEdiRotate:
    @pytest.fixture(autouse=True)
    def require_data(self) -> None:
        if not _has(_EDI_3EDIS):
            pytest.skip("data/3edis/ not found")

    def test_dry_run(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "rotate", str(_EDI_3EDIS), "--angle", "30", "--dry-run"]
        )
        assert result.exit_code == 0
        assert "dry" in result.output.lower()

    def test_rotate_writes_edis(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main, ["edi", "rotate", str(_EDI_3EDIS),
                   "--angle", "30", "--output-dir", str(tmp_path)]
        )
        assert result.exit_code == 0
        assert len(list(tmp_path.glob("*.edi"))) > 0

    def test_json_output(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main, ["edi", "rotate", str(_EDI_3EDIS),
                   "--angle", "45", "--output-dir", str(tmp_path),
                   "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "n_written" in data and data["angle_deg"] == 45.0

    def test_angle_required(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main, ["edi", "rotate", str(_EDI_3EDIS), "--output-dir", str(tmp_path)]
        )
        assert result.exit_code != 0

    def test_no_output_dir_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "rotate", str(_EDI_3EDIS), "--angle", "30"]
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# pycsamt edi select
# ---------------------------------------------------------------------------

class TestEdiSelect:
    @pytest.fixture(autouse=True)
    def require_data(self) -> None:
        if not _has(_EDI_3EDIS):
            pytest.skip("data/3edis/ not found")

    def test_dry_run_all(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "select", str(_EDI_3EDIS), "--dry-run"]
        )
        assert result.exit_code == 0
        assert "match" in result.output.lower() or "station" in result.output.lower()

    def test_dry_run_json(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "select", str(_EDI_3EDIS), "--dry-run", "--format", "json"]
        )
        assert result.exit_code == 0
        assert "n_selected" in json.loads(result.output)

    def test_select_by_station(self, runner: CliRunner, tmp_path: Path) -> None:
        from pycsamt.seg.collection import EDICollection
        coll = EDICollection.from_sources(_EDI_3EDIS)
        names = [e.station for e in coll if e.station]
        if not names:
            pytest.skip("No named stations in 3edis/")
        result = runner.invoke(
            main, ["edi", "select", str(_EDI_3EDIS),
                   "--stations", names[0], "--output-dir", str(tmp_path)]
        )
        assert result.exit_code == 0
        assert len(list(tmp_path.glob("*.edi"))) >= 1

    def test_select_json_output(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main, ["edi", "select", str(_EDI_3EDIS),
                   "--output-dir", str(tmp_path), "--format", "json"]
        )
        assert result.exit_code == 0
        raw = result.output
        json_start = raw.find("{")
        assert json_start >= 0, f"No JSON in output: {raw[:200]!r}"
        assert "n_written" in json.loads(raw[json_start:])

    def test_no_output_dir_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["edi", "select", str(_EDI_3EDIS)])
        assert result.exit_code != 0

    def test_has_tipper_filter_dry_run(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["edi", "select", str(_EDI_3EDIS),
                   "--has-tipper", "--dry-run", "--format", "json"]
        )
        assert result.exit_code == 0
        assert json.loads(result.output)["n_selected"] >= 0


# ---------------------------------------------------------------------------
# Unit tests — low-level helpers (import from edi, not seg)
# ---------------------------------------------------------------------------

class TestEdiHelpers:
    @pytest.fixture(autouse=True)
    def require_data(self) -> None:
        if not _has(_EDI_3EDIS):
            pytest.skip("data/3edis/ not found")

    def test_get_collection_returns_nonempty(self) -> None:
        from pycsamt.cli.commands.edi._base import (
            _get_collection,
        )
        coll = _get_collection(_EDI_3EDIS)
        assert len(coll) > 0

    def test_get_edi_returns_edifile(self) -> None:
        from pycsamt.cli.commands.edi._base import _get_edi
        from pycsamt.seg.edi import EDIFile
        edi_path = sorted(_EDI_3EDIS.glob("*.edi"))[0]
        edi = _get_edi(edi_path)
        assert isinstance(edi, EDIFile)
        assert edi.Z.n_freq > 0

    def test_rotate_impedance_changes_z(self) -> None:
        import numpy as np

        from pycsamt.cli.commands.edi._base import _get_edi
        from pycsamt.seg.ops import rotate_impedance
        edi_path = sorted(_EDI_3EDIS.glob("*.edi"))[0]
        edi = _get_edi(edi_path)
        z_orig = edi.Z.z.copy()
        z_rot  = rotate_impedance(z_orig, 45.0)
        assert z_rot.shape == z_orig.shape
        assert not np.allclose(z_rot, z_orig)

    def test_rotate_back_gives_original(self) -> None:
        import numpy as np

        from pycsamt.cli.commands.edi._base import _get_edi
        from pycsamt.seg.ops import rotate_impedance
        edi_path = sorted(_EDI_3EDIS.glob("*.edi"))[0]
        edi = _get_edi(edi_path)
        z_orig = edi.Z.z.copy()
        z_rot  = rotate_impedance(z_orig, 45.0)
        z_back = rotate_impedance(z_rot, -45.0)
        assert np.allclose(z_back, z_orig, atol=1e-10)
